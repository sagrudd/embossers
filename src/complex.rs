//! Linguistic sequence complexity (EMBOSS `complex`).
//!
//! This module computes *linguistic complexity* by sliding a window across
//! input sequences and, for each k in `[jmin, jmax]`, measuring the fraction
//! of distinct k‑mers observed relative to the theoretical maximum. The final
//! complexity is the product of the per‑k fractions averaged over windows.
//!
//! ### Definition
//! For a window `W` of length `L` and a k‑mer size `k`,
//! ```text
//! U_k(W) = |distinct k-mers in W restricted to A,C,G,T| / min(4^k, L-k+1)
//! ```
//! The **linguistic complexity** is `∏_{k=jmin..jmax} mean_W U_k(W)`.
//!
//! ### Simulations
//! If `sim > 0`, random DNA windows are generated either uniformly or using
//! empirical base frequencies, and the same U_k statistics are computed to
//! give a baseline (mean ± sd) for comparison.
//!
//! ### Complexity & Performance
//! Time is `O(N * (jmax-jmin+1))` per window; k‑mer counting uses a hashed
//! slice view and is fast for typical `j ≤ 12` and `L ≤ 1e5`.
//!
//! ### Examples
//! ```rust
//! use embossers::{ComplexOptions, compute_complexity};
//! let opts = ComplexOptions { lwin: 8, step: 4, jmin: 2, jmax: 4, sim: 0, freq_weighted_sim: false };
//! let (c, rows) = compute_complexity("ACGTACGT", &opts).unwrap();
//! assert!(c > 0.0);
//! assert!(!rows.is_empty());
//! ```
//!
use std::collections::HashSet;
use crate::common::EmbossersError;

/// Options controlling the linguistic complexity computation.

#[derive(Clone, Debug)]
/// Options that control sliding window size, step and the k‑mer range.
///
/// * `lwin` — window length (L). If an input sequence is shorter than `lwin`,
///   a single window covering the whole sequence is used.
/// * `step` — slide increment between consecutive windows.
/// * `jmin`, `jmax` — inclusive bounds for k‑mer sizes considered.
/// * `sim` — number of Monte‑Carlo simulations for a baseline (0 disables).
/// * `freq_weighted_sim` — if `true`, simulations match empirical A/C/G/T
///   frequencies; otherwise uniform 0.25 for each base.
pub struct ComplexOptions {
    /// Sliding window length (L). If the sequence is shorter than `lwin`, one
    /// window covering the entire sequence is used.
    pub lwin: usize,
    /// Step in bases to slide the window.
    pub step: usize,
    /// Minimum word length (inclusive).
    pub jmin: usize,
    /// Maximum word length (inclusive).
    pub jmax: usize,
    /// Number of simulations to run (0 = none).
    pub sim: usize,
    /// If `true`, simulations match the per-base frequencies of each input
    /// sequence; otherwise simulated bases are sampled uniformly over A,C,G,T.
    pub freq_weighted_sim: bool,
}

impl Default for ComplexOptions {
    fn default() -> Self {
        Self { lwin: 100, step: 5, jmin: 4, jmax: 6, sim: 0, freq_weighted_sim: false }
    }
}

/// Result rows for a Uj table.
#[derive(Clone, Debug)]
/// A single result line for a given `j` (k‑mer length).
pub struct UjRow {
    /// Word length (k).
    pub j: usize,
    /// Mean `U_j` across windows of the real sequence(s).
    pub uj_real_mean: f64,
    /// (If simulations requested) mean `U_j` across windows and simulations.
    pub uj_sim_mean: Option<f64>,
    /// (If simulations requested) standard deviation of `U_j` across simulations.
    pub uj_sim_std: Option<f64>,
}

/// Compute linguistic complexity for a single sequence with options.

/// Compute linguistic complexity for **one** sequence and return the scalar
/// complexity along with per‑k rows.
///
/// See [`compute_complexity_multi`] for a multi‑sequence aggregate.
pub fn compute_complexity(seq: &str, opts: &ComplexOptions) -> Result<(f64, Vec<UjRow>), EmbossersError> {
    compute_complexity_multi([seq].into_iter(), opts)
}

/// Compute linguistic complexity across **multiple sequences** collectively.

/// Compute linguistic complexity across multiple sequences collectively by
/// concatenating their windows for the per‑k means used in the product.
pub fn compute_complexity_multi<'a>(
    seqs: impl IntoIterator<Item = &'a str>,
    opts: &ComplexOptions
) -> Result<(f64, Vec<UjRow>), EmbossersError> {
    if opts.jmin == 0 || opts.jmax == 0 || opts.jmin > opts.jmax {
        return Err(EmbossersError::InvalidJRange { jmin: opts.jmin, jmax: opts.jmax });
    }
    if opts.lwin == 0 || opts.step == 0 {
        return Err(EmbossersError::InvalidWindow { lwin: opts.lwin, step: opts.step });
    }

    // Collect windows
    let mut windows: Vec<String> = Vec::new();
    let mut global_jmax = 0usize;
    let mut global_jmin = usize::MAX;
    for s in seqs {
        let seq_u = s.to_ascii_uppercase();
        let n = seq_u.len();
        if n == 0 { continue; }
        let lwin = opts.lwin.min(n);
        let jmax = opts.jmax.min(lwin);
        let jmin = opts.jmin.min(jmax);
        global_jmax = global_jmax.max(jmax);
        global_jmin = global_jmin.min(jmin);
        if n <= lwin {
            windows.push(seq_u);
        } else {
            let mut i = 0usize;
            while i + lwin <= n {
                windows.push(seq_u[i..i + lwin].to_string());
                i += opts.step;
            }
            if i < n && windows.is_empty() {
                windows.push(seq_u[n - lwin..].to_string());
            }
        }
    }
    if windows.is_empty() {
        return Ok((0.0, Vec::new()));
    }
    let jmin = global_jmin;
    let jmax = global_jmax;

    // Real Uj per window.
    let mut uj_sums = vec![0.0; jmax - jmin + 1];
    for w in &windows {
        let l = w.len();
        for (idx, j) in (jmin..=jmax).enumerate() {
            let denom = max_vocab(l, j);
            let uj = if denom == 0 { 0.0 } else { distinct_kmers(w, j) as f64 / denom as f64 };
            uj_sums[idx] += uj;
        }
    }
    let win_count = windows.len() as f64;
    let uj_means: Vec<f64> = uj_sums.iter().map(|s| s / win_count).collect();
    let real_complexity = uj_means.iter().product::<f64>();

    // Simulations
    let mut rows: Vec<UjRow> = vec![];
    if opts.sim == 0 {
        for (i, j) in (jmin..=jmax).enumerate() {
            rows.push(UjRow { j, uj_real_mean: uj_means[i], uj_sim_mean: None, uj_sim_std: None });
        }
        return Ok((real_complexity, rows));
    }

    // Empirical base frequencies if requested
    let concat: String = windows.join("");
    let (pa, pc, pg, pt) = if opts.freq_weighted_sim { base_freqs(&concat) } else { (0.25,0.25,0.25,0.25) };

    let mut rng = rand::thread_rng();
    let mut sim_uj_acc = vec![vec![]; jmax - jmin + 1];

    for _ in 0..opts.sim {
        let mut uj_sums_sim = vec![0.0; jmax - jmin + 1];
        for w in &windows {
            let l = w.len();
            let simw = random_dna(l, pa, pc, pg, pt, &mut rng);
            for (idx, j) in (jmin..=jmax).enumerate() {
                let denom = max_vocab(l, j);
                let uj = if denom == 0 { 0.0 } else { distinct_kmers(&simw, j) as f64 / denom as f64 };
                uj_sums_sim[idx] += uj;
            }
        }
        for (idx, sum) in uj_sums_sim.into_iter().enumerate() {
            sim_uj_acc[idx].push(sum / win_count);
        }
    }

    for (i, j) in (jmin..=jmax).enumerate() {
        let vec = &sim_uj_acc[i];
        let mean = vec.iter().copied().sum::<f64>() / vec.len() as f64;
        let var = vec.iter().map(|x| (x-mean)*(x-mean)).sum::<f64>() / (vec.len().saturating_sub(1).max(1)) as f64;
        rows.push(UjRow { j, uj_real_mean: uj_means[i], uj_sim_mean: Some(mean), uj_sim_std: Some(var.sqrt()) });
    }

    Ok((real_complexity, rows))
}

/// Compute the maximum possible vocabulary size for `j`-mers in a window of
/// length `l`, as `min(4^j, l - j + 1)`.

/// The theoretical maximum number of distinct k‑mers `min(4^j, l-j+1)` for a
/// window of length `l`. Returns 0 when `l < j`.
pub fn max_vocab(l: usize, j: usize) -> usize {
    if l < j { return 0; }
    let pow = 4usize.saturating_pow(j as u32);
    pow.min(l - j + 1)
}

/// Count the number of distinct `j`-mers over `ACGT` in `window`.
fn distinct_kmers(window: &str, j: usize) -> usize {
    if window.len() < j { return 0; }
    let bytes = window.as_bytes();
    let mut set: HashSet<&[u8]> = HashSet::new();
    let mut i = 0;
    while i + j <= bytes.len() {
        let k = &bytes[i..i+j];
        if k.iter().all(|&b| matches!(b, b'A'|b'C'|b'G'|b'T')) {
            set.insert(k);
        }
        i += 1;
    }
    set.len()
}

/// Compute base frequencies (A,C,G,T) in a sequence, ignoring other symbols.
fn base_freqs(seq: &str) -> (f64,f64,f64,f64) {
    let mut a=0usize; let mut c=0usize; let mut g=0usize; let mut t=0usize;
    for b in seq.bytes() {
        match b { b'A' => a+=1, b'C'=>c+=1, b'G'=>g+=1, b'T'=>t+=1, _=>{} }
    }
    let tot = (a+c+g+t).max(1) as f64;
    (a as f64/tot, c as f64/tot, g as f64/tot, t as f64/tot)
}

/// Generate a random DNA string of length `l` with given base probabilities.
fn random_dna<R: rand::Rng>(l: usize, pa: f64, pc: f64, pg: f64, _pt: f64, rng: &mut R) -> String {
    // No explicit `use rand::Rng;` needed; bound provides method resolution.
    let mut s = String::with_capacity(l);
    for _ in 0..l {
        let x: f64 = rng.gen();
        let mut acc = pa;
        let ch = if x < acc { 'A' } else { acc+=pc; if x<acc {'C'} else { acc+=pg; if x<acc {'G'} else {'T'} } };
        s.push(ch);
    }
    s
}
