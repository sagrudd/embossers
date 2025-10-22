//! # embossers
//!
//! EMBOSS-inspired utilities implemented in Rust.
//!
//! This crate provides a Rust implementation of the **EMBOSS `complex`** algorithm,
//! computing *linguistic sequence complexity* for nucleotide sequences using
//! a sliding-window approach and configurable k-mer sizes.
//!
//! ## What is linguistic sequence complexity?
//! For a word size *j*, the *vocabulary usage* `U_j` in a window of length *L* is:
//!
//! ```text
//! U_j = observed_unique_jmers / min(4^j, L - j + 1)
//! ```
//!
//! The window complexity is the product of `U_j` over a range `j = jmin..=jmax`.
//! The sequence complexity is the mean of window complexities across all windows.
//!
//! This mirrors the behaviour described in the EMBOSS `complex` manual while
//! being explicit about assumptions and edge cases.
//!
//! ### Example
//! ```
//! use embossers::{ComplexOptions, compute_complexity};
//! let seq = "ACGTACGTACGT";
//! let opts = ComplexOptions{ lwin: 12, step: 6, jmin: 2, jmax: 3, sim: 0, freq_weighted_sim: false };
//! let (c, uj_rows) = compute_complexity(seq, &opts).unwrap();
//! assert!(c >= 0.0 && c <= 1.0);
//! assert!(!uj_rows.is_empty());
//! ```
#![cfg_attr(docsrs, feature(doc_cfg, doc_auto_cfg))]

use std::collections::HashSet;

/// Errors that can be returned by the `complex` implementation.
#[derive(thiserror::Error, Debug)]
pub enum ComplexError {
    /// Returned if `jmin` > `jmax` or either is less than 1.
    #[error("invalid j-range: jmin={jmin}, jmax={jmax}")]
    InvalidJRange { jmin: usize, jmax: usize },
    /// Returned if `lwin` or `step` is zero.
    #[error("window length and step must be > 0 (lwin={lwin}, step={step})")]
    InvalidWindow { lwin: usize, step: usize },
}

/// A single FASTA sequence (identifier and uppercase sequence letters).
#[derive(Clone, Debug)]
pub struct FastaRecord {
    /// Identifier from the FASTA header (text after '>').
    pub id: String,
    /// Raw sequence (uppercase). Non-ACGT symbols are kept as-is; k-mers that
    /// include non-ACGT are **ignored** when computing `U_j`.
    pub seq: String,
}

/// Options controlling the linguistic complexity computation.
#[derive(Clone, Debug)]
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
///
/// Returns `(complexity, uj_rows)` where `complexity` is the product of
/// mean `U_j` values across windows, and `uj_rows` summarizes `U_j` across
/// windows (and simulations if requested).
///
/// # Notes
/// - K-mers containing characters outside `A|C|G|T` are skipped.
/// - If `jmax` exceeds the window length `L`, it is clamped to `L`.
pub fn compute_complexity(seq: &str, opts: &ComplexOptions) -> Result<(f64, Vec<UjRow>), ComplexError> {
    compute_complexity_multi([seq].into_iter(), opts)
}

/// Compute linguistic complexity across **multiple sequences** collectively.
///
/// This corresponds to EMBOSS `-omnia`, where windows from all input sequences
/// are considered together for reporting.
pub fn compute_complexity_multi<'a>(
    seqs: impl IntoIterator<Item = &'a str>,
    opts: &ComplexOptions
) -> Result<(f64, Vec<UjRow>), ComplexError> {
    if opts.jmin == 0 || opts.jmax == 0 || opts.jmin > opts.jmax {
        return Err(ComplexError::InvalidJRange { jmin: opts.jmin, jmax: opts.jmax });
    }
    if opts.lwin == 0 || opts.step == 0 {
        return Err(ComplexError::InvalidWindow { lwin: opts.lwin, step: opts.step });
    }

    // First pass to collect all windows (as owned Strings to keep lifetimes simple)
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
            // ensure at least one trailing window if step doesn't land on the end
            if i < n && windows.is_empty() {
                windows.push(seq_u[n - lwin..].to_string());
            }
        }
    }
    if windows.is_empty() {
        // No windows -> define complexity as 0 with empty table.
        return Ok((0.0, Vec::new()));
    }
    let jmin = global_jmin;
    let jmax = global_jmax;

    // Precompute real Uj per window.
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

    // Empirical base frequencies if requested — computed over concatenated sequence
    let concat: String = windows.join("");
    let (pa, pc, pg, pt) = if opts.freq_weighted_sim { base_freqs(&concat) } else { (0.25,0.25,0.25,0.25) };

    let mut rng = rand::thread_rng();
    let mut sim_uj_acc = vec![vec![]; jmax - jmin + 1];

    for _ in 0..opts.sim {
        let mut uj_sums_sim = vec![0.0; jmax - jmin + 1];
        for w in &windows {
            let l = w.len();
            // Generate a synthetic window with same length
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
    let mut s = String::with_capacity(l);
    for _ in 0..l {
        let x: f64 = rand::Rng::gen(rng);
        let mut acc = pa;
        let ch = if x < acc { 'A' } else { acc+=pc; if x<acc {'C'} else { acc+=pg; if x<acc {'G'} else {'T'} } };
        s.push(ch);
    }
    s
}

/// Parse a simple FASTA string into records.
/// This is deliberately minimal and tolerant of whitespace.
pub fn parse_fasta(text: &str) -> Vec<FastaRecord> {
    let mut out: Vec<FastaRecord> = vec![];
    let mut id = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !id.is_empty() { out.push(FastaRecord{ id: id.clone(), seq: seq.to_ascii_uppercase() }); seq.clear(); }
            id = rest.trim().split_whitespace().next().unwrap_or("").to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !id.is_empty() { out.push(FastaRecord{ id, seq: seq.to_ascii_uppercase() }); }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn max_vocab_basic() {
        assert_eq!(max_vocab(10, 1), 4.min(10));
        assert_eq!(max_vocab(10, 3), 4usize.pow(3).min(8));
    }

    #[test]
    fn compute_simple() {
        let seq = "ACGTACGTACGT";
        let opts = ComplexOptions{ lwin: 12, step: 5, jmin: 2, jmax: 3, sim:0, freq_weighted_sim:false };
        let (c, rows) = compute_complexity(seq, &opts).unwrap();
        assert!(c>0.0 && c<=1.0);
        assert_eq!(rows[0].j, 2);
    }
}
