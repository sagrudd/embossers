//! EST-to-genome spliced alignment (EMBOSS `est2genome`).
//!
//! This implementation performs **semi‑global** alignment of an EST/cDNA
//! sequence against a genomic sequence with standard affine‑gap penalties and
//! post‑hoc *intron* detection. Any gap in the EST (a run of `D` operations)
//! whose length ≥ `intron_min` is classified as an intron. Canonical splice
//! sites (`GT-AG` on the genome) are detected and annotated; an optional
//! `splice_bonus` can be applied to the total reported score.
//!
//! ### Notes
//! - Semi‑global means the EST is aligned end‑to‑end while genome ends are
//!   free (unpenalized). This mirrors typical use of `est2genome`.
//! - We expose conservative scoring knobs and keep DNA scoring by default.
//! - For exact parity with EMBOSS output formatting and nuanced heuristics,
//!   adjust parameters accordingly or extend this module (e.g., GC-AG/AT-AC).
//!
//! ### Example
//! ```rust,no_run
//! use embossers::{est2genome, Est2GenomeParams, WaterMatrix};
//! let g = "ACGT"*100; // pretend genome
//! let e = "ACGTACGT";
//! let params = Est2GenomeParams{ matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 }, ..Default::default() };
//! let aln = est2genome(&g, &e, &params).unwrap();
//! assert!(aln.align_a.len() == aln.align_b.len());
//! ```
//!
use crate::common::{EmbossersError, WaterMatrix, blosum62_score, equals_case_insensitive};

/// Parameters for `est2genome` spliced alignment.
#[derive(Clone, Debug)]
pub struct Est2GenomeParams {
    /// Scoring matrix (DNA by default).
    pub matrix: WaterMatrix,
    /// Gap open penalty for short indels.
    pub gap_open: f32,
    /// Gap extension penalty for short indels.
    pub gap_extend: f32,
    /// Integer scale to preserve fractional penalties.
    pub scale: f32,
    /// Minimum length (in genomic bases) to classify a deletion in the EST as an intron.
    pub intron_min: usize,
    /// Bonus (added to final score) for each canonical `GT-AG` intron detected.
    pub splice_bonus: i32,
}

impl Default for Est2GenomeParams {
    fn default() -> Self {
        Self {
            matrix: WaterMatrix::Dna { match_score: 2, mismatch: -1 },
            gap_open: 10.0,
            gap_extend: 0.5,
            scale: 2.0,
            intron_min: 20,
            splice_bonus: 5,
        }
    }
}

/// A detected intron on the genomic sequence.
#[derive(Clone, Debug)]
pub struct Intron {
    /// Genomic coordinates **inclusive** `[start, end]` (0‑based) of the intron span.
    pub start: usize,
    pub end: usize,
    /// Whether the intron shows canonical `GT-AG` splice dinucleotides.
    pub canonical_gt_ag: bool,
}

/// Output of an EST→genome alignment.
#[derive(Clone, Debug)]
pub struct Est2GenomeAlignment {
    /// Total score (integer, scaled) including any `splice_bonus` applied.
    pub score: i32,
    /// Aligned genome string (with `-` for insertions in genome).
    pub align_a: String,
    /// Aligned EST string (with `-` for deletions/introns).
    pub align_b: String,
    /// CIGAR‑like string. Introns rendered as `N` blocks (e.g., `50M100N50M`).
    pub cigar: String,
    /// Introns detected with coordinates and canonical flags.
    pub introns: Vec<Intron>,
    /// Percent identity and percent gaps across aligned columns.
    pub pct_identity: f64,
    pub pct_gaps: f64,
}

/// Align an EST (`est`) to a genomic sequence (`genome`) with semi‑global DP,
/// then detect and annotate introns from long deletions in the EST.
pub fn est2genome(genome: &str, est: &str, params: &Est2GenomeParams) -> Result<Est2GenomeAlignment, EmbossersError> {
    if genome.is_empty() || est.is_empty() {
        return Err(EmbossersError::InvalidSequence("empty genome or est"));
    }
    let a: Vec<char> = genome.chars().collect();
    let b: Vec<char> = est.chars().collect();
    let n = a.len();
    let m = b.len();

    // scaling
    let scale = params.scale.max(1.0);
    let go = (params.gap_open * scale).round() as i32;
    let ge = (params.gap_extend * scale).round() as i32;

    // scorer
    let score_pair = |x: char, y: char| -> i32 {
        match &params.matrix {
            WaterMatrix::Dna{match_score, mismatch} => if x == y { *match_score } else { *mismatch },
            WaterMatrix::Blosum62 => blosum62_score(x, y),
        }
    };

    let neg_inf = i32::MIN / 4;
    // Semi‑global DP matrices (global on EST; genome ends free).
    // H = best, E = gap in genome (insertion in EST), F = gap in EST (deletion/intron).
    let mut h = vec![vec![neg_inf; m+1]; n+1];
    let mut e = vec![vec![neg_inf; m+1]; n+1];
    let mut f = vec![vec![neg_inf; m+1]; n+1];
    let mut tb = vec![vec![0u8; m+1]; n+1]; // 1=diag, 2=up(F), 3=left(E)

    h[0][0] = 0;
    // top row (j>0): standard global penalties (consumes EST, genome empty)
    for j in 1..=m {
        e[0][j] = if j==1 { h[0][j-1] - go } else { e[0][j-1] - ge };
        h[0][j] = e[0][j];
        tb[0][j] = 3;
    }
    // first column (i>0): free gaps in genome prefix (semi‑global)
    for i in 1..=n {
        h[i][0] = 0;
        f[i][0] = 0;
        tb[i][0] = 2;
    }

    for i in 1..=n {
        for j in 1..=m {
            e[i][j] = (h[i][j-1] - go).max(e[i][j-1] - ge);
            f[i][j] = (h[i-1][j] - go).max(f[i-1][j] - ge);
            let diag = h[i-1][j-1] + score_pair(a[i-1], b[j-1]);
            let (best, dir) = [
                (diag, 1u8),
                (f[i][j], 2u8),
                (e[i][j], 3u8),
            ].into_iter().max_by_key(|(v,_)| *v).unwrap();
            h[i][j] = best;
            tb[i][j] = dir;
        }
    }

    // traceback from best cell in last column (free genome suffix)
    let mut best_i = 0usize;
    let mut best_score = neg_inf;
    for i in 0..=n {
        if h[i][m] > best_score { best_score = h[i][m]; best_i = i; }
    }

    let mut i = best_i;
    let mut j = m;
    let mut a_aln: Vec<char> = Vec::new();
    let mut b_aln: Vec<char> = Vec::new();
    let mut cigar_ops: Vec<(char, usize)> = Vec::new(); // M/I/D (later D>=intron_min -> N)
    while j>0 { // must consume full EST
        if i>0 && j>0 && tb[i][j] == 1 {
            a_aln.push(a[i-1]); b_aln.push(b[j-1]); push_cigar(&mut cigar_ops, 'M', 1); i-=1; j-=1;
        } else if i>0 && (j==0 || tb[i][j] == 2) {
            a_aln.push(a[i-1]); b_aln.push('-'); push_cigar(&mut cigar_ops, 'D', 1); i-=1;
        } else if j>0 && (i==0 || tb[i][j] == 3) {
            a_aln.push('-'); b_aln.push(b[j-1]); push_cigar(&mut cigar_ops, 'I', 1); j-=1;
        } else {
            // Fallback to diagonal to avoid infinite loop
            if i>0 && j>0 { a_aln.push(a[i-1]); b_aln.push(b[j-1]); push_cigar(&mut cigar_ops, 'M', 1); i-=1; j-=1; }
            else { break; }
        }
    }
    // If genome prefix remains (i>0) it's free; prepend as deletions in EST (not necessary to print)
    while i>0 {
        a_aln.push(a[i-1]); b_aln.push('-'); push_cigar(&mut cigar_ops, 'D', 1); i-=1;
    }

    a_aln.reverse(); b_aln.reverse();
    cigar_ops.reverse();

    // Detect introns (convert long D runs to N; annotate canonical GT-AG)
    let (cigar, introns, splice_bonus_total) = build_spliced_cigar_and_introns(&a, &a_aln, &b_aln, &cigar_ops, params.intron_min, params.splice_bonus);

    // Identity/gaps
    let (mut ident, mut gaps) = (0usize, 0usize);
    for (&x,&y) in a_aln.iter().zip(b_aln.iter()) {
        if x == '-' || y == '-' { gaps += 1; } else if equals_case_insensitive(x,y) { ident += 1; }
    }
    let cols = a_aln.len().max(1);
    let pct_identity = (ident as f64)*100.0/(cols as f64);
    let pct_gaps = (gaps as f64)*100.0/(cols as f64);

    Ok(Est2GenomeAlignment{
        score: best_score + splice_bonus_total,
        align_a: a_aln.into_iter().collect(),
        align_b: b_aln.into_iter().collect(),
        cigar,
        introns,
        pct_identity,
        pct_gaps,
    })
}

fn push_cigar(ops: &mut Vec<(char,usize)>, op: char, k: usize) {
    if let Some(last) = ops.last_mut() {
        if last.0 == op { last.1 += k; return; }
    }
    ops.push((op, k));
}

fn build_spliced_cigar_and_introns(
    genome: &[char],
    aln_a: &[char],
    aln_b: &[char],
    cigar_ops: &[(char,usize)],
    intron_min: usize,
    splice_bonus: i32,
) -> (String, Vec<Intron>, i32) {
    // Build CIGAR while converting long D to N.
    let mut cigar = String::new();
    let mut introns: Vec<Intron> = Vec::new();
    // track mapping to genomic indices
    let mut g_idx: isize = -1; // will increment when aln_a char != '-'
    let mut bonus_total = 0i32;
    let mut ai = 0usize;
    for &(op, len) in cigar_ops {
        match op {
            'M' => { cigar.push_str(&format!("{}M", len)); for _ in 0..len { ai+=1; if aln_a[ai-1] != '-' { g_idx += 1; } } },
            'I' => { cigar.push_str(&format!("{}I", len)); ai += len; /* aln_a has '-' in these positions; g_idx unchanged */ },
            'D' => {
                // Count real genomic span covered by these D's in alignment
                let mut span = 0usize;
                let start_g = (g_idx+1) as isize;
                for _ in 0..len { // among these D's, aln_a has bases, aln_b has '-' so genome advances
                    ai += 1;
                    if aln_a[ai-1] != '-' { g_idx += 1; span += 1; }
                }
                if span >= intron_min {
                    // canonical check on genome coordinates [start_g, start_g+span-1]
                    let (start, end) = (start_g as usize, (start_g as usize)+span-1);
                    let canonical = is_canonical_gt_ag(genome, start, end);
                    if canonical { bonus_total += splice_bonus; }
                    introns.push(Intron{ start, end, canonical_gt_ag: canonical });
                    cigar.push_str(&format!("{}N", span));
                } else {
                    cigar.push_str(&format!("{}D", span));
                }
            }
            _ => {}
        }
    }
    (cigar, introns, bonus_total)
}

fn is_canonical_gt_ag(genome: &[char], start: usize, end: usize) -> bool {
    if start+1 < genome.len() && end >= 1 && end < genome.len() {
        let donor = (genome[start], genome[start+1]);
        let acceptor = (genome[end-1], genome[end]);
        donor.0.to_ascii_uppercase()=='G' && donor.1.to_ascii_uppercase()=='T' &&
        acceptor.0.to_ascii_uppercase()=='A' && acceptor.1.to_ascii_uppercase()=='G'
    } else { false }
}
