//! Sequence merging by local overlap (EMBOSS `merger`/`megamerger`‑like).
//!
//! A pragmatic reimplementation that finds the best **local** alignment between
//! two sequences (using our `matcher` engine) and merges them into a consensus
//! over the overlap, concatenating the non-overlapping flanks. DNA uses IUPAC
//! ambiguity on mismatches; proteins emit `X` for conflicting residues.
//!
use crate::common::{EmbossersError, WaterMatrix, blosum62_score};
use crate::matcher::{matcher, MatcherParams, MatcherHit};

/// Parameters shared by merger/megamerger.
#[derive(Clone, Debug)]
pub struct MergeParams {
    /// Scoring matrix for local alignment (DNA or BLOSUM62).
    pub matrix: WaterMatrix,
    /// Gap open penalty.
    pub gap_open: f32,
    /// Gap extension penalty.
    pub gap_extend: f32,
    /// Integer scaling.
    pub scale: f32,
    /// Minimum local score to accept merging.
    pub min_score: i32,
}

impl Default for MergeParams {
    fn default() -> Self {
        Self {
            matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 },
            gap_open: 10.0,
            gap_extend: 0.5,
            scale: 2.0,
            min_score: 1,
        }
    }
}

/// Result of a merge operation.
#[derive(Clone, Debug)]
pub struct MergeResult {
    /// Consensus merged sequence.
    pub merged: String,
    /// Score of the overlap alignment used.
    pub score: i32,
    /// A and B overlap ranges used (0‑based inclusive).
    pub a_start: usize, pub a_end: usize,
    pub b_start: usize, pub b_end: usize,
}

/// Merge two **DNA** sequences via local overlap (IUPAC in conflicts).
pub fn merger(a: &str, b: &str, params: &MergeParams) -> Result<MergeResult, EmbossersError> {
    merge_impl(a, b, params, true)
}

/// Merge two **protein** sequences via local overlap (`X` in conflicts).
pub fn megamerger(a: &str, b: &str, params: &MergeParams) -> Result<MergeResult, EmbossersError> {
    merge_impl(a, b, params, false)
}

fn merge_impl(a: &str, b: &str, params: &MergeParams, dna: bool) -> Result<MergeResult, EmbossersError> {
    let mparams = MatcherParams{
        matrix: params.matrix.clone(),
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        scale: params.scale,
        alternatives: 1,
        min_score: params.min_score,
    };
    let hits = matcher(a, b, &mparams)?;
    if hits.is_empty() { return Err(EmbossersError::InvalidSequence("no sufficient overlap")); }
    let h = &hits[0];

    // Build consensus over the aligned strings
    let aa: Vec<char> = h.align_a.chars().collect();
    let bb: Vec<char> = h.align_b.chars().collect();
    let mut ovl = String::with_capacity(aa.len());
    for (&x,&y) in aa.iter().zip(bb.iter()) {
        if x=='-' { ovl.push(y); }
        else if y=='-' { ovl.push(x); }
        else if x.eq_ignore_ascii_case(&y) { ovl.push(x.to_ascii_uppercase()); }
        else {
            if dna { ovl.push(iupac_pair(x,y)); } else { ovl.push('X'); }
        }
    }

    // Stitch: prefix + ovl + suffix
    let mut merged = String::new();
    // Prefix: up to a_start and b_start choose the one that starts earlier
    if h.a_start <= h.b_start {
        merged.push_str(&a[..h.a_start]);
    } else {
        merged.push_str(&b[..h.b_start]);
    }
    merged.push_str(&ovl);
    // Suffix: take the tail from the sequence that extends further
    let a_tail = &a[h.a_end+1..];
    let b_tail = &b[h.b_end+1..];
    if a_tail.len() >= b_tail.len() {
        merged.push_str(a_tail);
    } else {
        merged.push_str(b_tail);
    }

    Ok(MergeResult{
        merged,
        score: h.score,
        a_start: h.a_start, a_end: h.a_end,
        b_start: h.b_start, b_end: h.b_end,
    })
}

fn iupac_pair(a: char, b: char) -> char {
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    let (a, b) = (if a=='U' {'T'} else {a}, if b=='U' {'T'} else {b});
    match (a,b) {
        ('A','C') | ('C','A') => 'M',
        ('A','G') | ('G','A') => 'R',
        ('A','T') | ('T','A') => 'W',
        ('C','G') | ('G','C') => 'S',
        ('C','T') | ('T','C') => 'Y',
        ('G','T') | ('T','G') => 'K',
        _ => 'N'
    }
}
