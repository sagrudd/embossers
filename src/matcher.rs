//! Waterman–Eggert-style local alignments (EMBOSS `matcher`).
#![allow(clippy::needless_return)]
//!
//! This module enumerates the **top-k non-overlapping local alignments** between
//! two sequences using a declumped Smith–Waterman (Waterman–Eggert) approach.
//! Each iteration runs a local DP with affine gaps, performs a traceback to get
//! the best alignment, and then **forbids the visited DP cells** so subsequent
//! passes yield alternative, non-overlapping alignments.
//!
//! This mirrors EMBOSS **matcher** behaviour: multiple best local alignments
//! under affine gaps and a configurable substitution matrix. See EMBOSS docs
//! for reference behaviour and parameters.
//!
//! ### Example
//! ```rust,no_run
//! use embossers::{matcher, MatcherParams, WaterMatrix};
//! let params = MatcherParams{ alternatives: 3, matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 }, ..Default::default() };
//! let outs = matcher("GATTACA", "CAGTTA", &params).unwrap();
//! assert!(!outs.is_empty());
//! ```
//!
use std::collections::HashSet;
use crate::common::{EmbossersError, WaterMatrix, blosum62_score, equals_case_insensitive};

/// Parameters for `matcher` (Waterman–Eggert local alignments).
#[derive(Clone, Debug)]
pub struct MatcherParams {
    /// Scoring matrix (DNA or BLOSUM62).
    pub matrix: WaterMatrix,
    /// Gap open penalty.
    pub gap_open: f32,
    /// Gap extension penalty.
    pub gap_extend: f32,
    /// Integer scale factor to preserve fractional penalties (e.g., 2.0 maps 0.5→1).
    pub scale: f32,
    /// Number of non-overlapping alternatives to return (like `-alternatives` / `-alt`).
    pub alternatives: usize,
    /// Minimum score threshold (alignments scoring below are ignored).
    pub min_score: i32,
}

impl Default for MatcherParams {
    fn default() -> Self {
        Self {
            matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 },
            gap_open: 10.0,
            gap_extend: 0.5,
            scale: 2.0,
            alternatives: 1,
            min_score: 1,
        }
    }
}

/// One local alignment hit.
#[derive(Clone, Debug)]
pub struct MatcherHit {
    /// Alignment score (scaled as per parameters).
    pub score: i32,
    /// 0-based inclusive coordinates in A and B for the aligned subranges.
    pub a_start: usize, pub a_end: usize,
    pub b_start: usize, pub b_end: usize,
    /// Aligned strings with `-` for gaps.
    pub align_a: String,
    pub align_b: String,
    /// CIGAR-like string (e.g., `40M1I10M2D`).
    pub cigar: String,
    /// Percent identity and percent gaps across aligned columns.
    pub pct_identity: f64,
    pub pct_gaps: f64,
}

/// Compute up to `params.alternatives` **non-overlapping** local alignments using
/// a Waterman–Eggert-style iterative DP with affine gaps.
pub fn matcher(a: &str, b: &str, params: &MatcherParams) -> Result<Vec<MatcherHit>, EmbossersError> {
    if a.is_empty() || b.is_empty() {
        return Err(EmbossersError::InvalidSequence("empty sequence"));
    }
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();
    let n = a_chars.len();
    let m = b_chars.len();

    let scale = params.scale.max(1.0);
    let go = (params.gap_open * scale).round() as i32;
    let ge = (params.gap_extend * scale).round() as i32;
    let score_pair = |x: char, y: char| -> i32 {
        match &params.matrix {
            WaterMatrix::Dna{match_score, mismatch} => if x == y { *match_score } else { *mismatch },
            WaterMatrix::Blosum62 => blosum62_score(x, y),
        }
    };

    let mut out: Vec<MatcherHit> = Vec::new();
    let mut forbidden: HashSet<(usize,usize)> = HashSet::new();

    for _k in 0..params.alternatives {
        // Smith–Waterman with affine gaps and per-cell forbidding.
        // Matrices are local: reset to 0 if negative.
        let neg_inf = i32::MIN/4;
        let mut h = vec![vec![0i32; m+1]; n+1];
        let mut e = vec![vec![neg_inf; m+1]; n+1];
        let mut f = vec![vec![neg_inf; m+1]; n+1];
        // Track best cell
        let mut best = 0i32;
        let mut bi = 0usize; let mut bj = 0usize;

        for i in 1..=n {
            for j in 1..=m {
                // If this cell is forbidden from previous alignments, zero it.
                if forbidden.contains(&(i,j)) {
                    h[i][j] = 0; e[i][j] = neg_inf; f[i][j] = neg_inf;
                    continue;
                }
                e[i][j] = (h[i][j-1] - go).max(e[i][j-1] - ge);
                f[i][j] = (h[i-1][j] - go).max(f[i-1][j] - ge);
                let diag = h[i-1][j-1] + score_pair(a_chars[i-1], b_chars[j-1]);
                let val = 0.max(diag).max(e[i][j]).max(f[i][j]);
                h[i][j] = val;
                if val > best {
                    best = val; bi = i; bj = j;
                }
            }
        }
        if best < params.min_score || best == 0 {
            break;
        }

        // Traceback until score falls to 0
        let mut i = bi; let mut j = bj;
        let mut aln_a: Vec<char> = Vec::new();
        let mut aln_b: Vec<char> = Vec::new();
        let mut cigar_ops: Vec<(char,usize)> = Vec::new();
        let mut used_cells: Vec<(usize,usize)> = Vec::new();

        while i>0 && j>0 && h[i][j] > 0 {
            used_cells.push((i,j));
            // Prefer diagonal if it produces the current cell
            let score_here = h[i][j];
            let diag_score = h[i-1][j-1] + score_pair(a_chars[i-1], b_chars[j-1]);
            if score_here == diag_score {
                aln_a.push(a_chars[i-1]); aln_b.push(b_chars[j-1]); push_cigar(&mut cigar_ops, 'M', 1);
                i-=1; j-=1;
            } else if score_here == e[i][j] {
                // came from left -> insertion in A (gap in A => A has '-', B has base)
                aln_a.push('-'); aln_b.push(b_chars[j-1]); push_cigar(&mut cigar_ops, 'I', 1);
                j-=1;
            } else if score_here == f[i][j] {
                // came from up -> deletion in A (gap in B => A has base, B has '-')
                aln_a.push(a_chars[i-1]); aln_b.push('-'); push_cigar(&mut cigar_ops, 'D', 1);
                i-=1;
            } else {
                // fallback (shouldn't happen often): break
                break;
            }
        }
        let a_end = i; // actually these are index after consumption; will recompute below
        aln_a.reverse(); aln_b.reverse();
        cigar_ops.reverse();

        // Coordinates: find start/end in original strings by counting non-gaps in aln_a/aln_b
        let mut a_pos = (bi, 0usize); // we'll reconstruct below properly
        let (a_start, a_end, b_start, b_end) = compute_local_coords(&aln_a, &aln_b, bi, bj, &a_chars, &b_chars);

        // Compute stats
        let (mut ident, mut gaps) = (0usize, 0usize);
        for (&x,&y) in aln_a.iter().zip(aln_b.iter()) {
            if x=='-' || y=='-' { gaps += 1; } else if equals_case_insensitive(x,y) { ident += 1; }
        }
        let cols = aln_a.len().max(1);
        let pct_identity = (ident as f64)*100.0/(cols as f64);
        let pct_gaps = (gaps as f64)*100.0/(cols as f64);
        let cigar = cigar_ops.into_iter().map(|(op,len)| format!("{len}{op}")).collect::<String>();

        // Save
        out.push(MatcherHit{
            score: best,
            a_start, a_end, b_start, b_end,
            align_a: aln_a.iter().collect(),
            align_b: aln_b.iter().collect(),
            cigar,
            pct_identity, pct_gaps,
        });

        // Forbid all cells on the used traceback path so we don't reuse them.
        for cell in used_cells { forbidden.insert(cell); }
    }

    Ok(out)
}

fn push_cigar(ops: &mut Vec<(char,usize)>, op: char, k: usize) {
    if let Some(last) = ops.last_mut() { if last.0 == op { last.1 += k; return; } }
    ops.push((op, k));
}

fn compute_local_coords(
    aln_a: &[char], aln_b: &[char],
    end_i: usize, end_j: usize,
    a: &[char], b: &[char],
) -> (usize, usize, usize, usize) {
    // Walk backwards from (end_i, end_j) by counting consumed residues during traceback length.
    let mut i = end_i;
    let mut j = end_j;
    // Count residues consumed in the traceback
    let consumed_a = aln_a.iter().filter(|&&c| c!='-').count();
    let consumed_b = aln_b.iter().filter(|&&c| c!='-').count();
    let a_end = end_i; let b_end = end_j;
    let a_start = a_end - consumed_a;
    let b_start = b_end - consumed_b;
    (a_start, a_end.saturating_sub(1), b_start, b_end.saturating_sub(1))
}
