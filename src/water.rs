//! Smith–Waterman local alignment with affine gaps (EMBOSS `water`).
use crate::common::{EmbossersError, WaterMatrix, blosum62_score, equals_case_insensitive};

/// Parameters for `water` Smith–Waterman.
#[derive(Clone, Debug)]
pub struct WaterParams {
    /// Scoring matrix to use.
    pub matrix: WaterMatrix,
    /// Gap open penalty (can be fractional; default 10.0 to mimic EMBOSS).
    pub gap_open: f32,
    /// Gap extension penalty (can be fractional; default 0.5 to mimic EMBOSS).
    pub gap_extend: f32,
    /// Integer scale factor to preserve fractional penalties (e.g. 2.0 so 0.5→1).
    pub scale: f32,
}

impl Default for WaterParams {
    fn default() -> Self {
        Self {
            matrix: WaterMatrix::Dna { match_score: 2, mismatch: -1 },
            gap_open: 10.0,
            gap_extend: 0.5,
            scale: 2.0,
        }
    }
}

/// A local alignment computed by Smith–Waterman with affine gaps.
#[derive(Clone, Debug)]
pub struct WaterAlignment {
    /// Total alignment **score** (integer, scaled as provided).
    pub score: i32,
    /// Start (inclusive) and end (exclusive) in sequence A for the aligned region.
    pub range_a: (usize, usize),
    /// Start (inclusive) and end (exclusive) in sequence B for the aligned region.
    pub range_b: (usize, usize),
    /// Aligned string for A (including gaps `-`).
    pub align_a: String,
    /// Aligned string for B (including gaps `-`).
    pub align_b: String,
    /// CIGAR-like operations (e.g. `10M1I5M2D`).
    pub cigar: String,
    /// Percent identity over aligned columns (0..=100).
    pub pct_identity: f64,
    /// Percent gaps over aligned columns (0..=100).
    pub pct_gaps: f64,
}

/// Run Smith–Waterman local alignment with affine gaps.
pub fn water(a: &str, b: &str, params: &WaterParams) -> Result<WaterAlignment, EmbossersError> {
    if a.is_empty() || b.is_empty() {
        return Err(EmbossersError::InvalidSequence("empty sequence"));
    }
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let n = a.len();
    let m = b.len();
    // scaling
    let scale = params.scale.max(1.0);
    let go = (params.gap_open * scale).round() as i32;
    let ge = (params.gap_extend * scale).round() as i32;

    // Scorer
    let score_pair = |x: char, y: char| -> i32 {
        match &params.matrix {
            WaterMatrix::Dna{match_score, mismatch} => if x == y { *match_score } else { *mismatch },
            WaterMatrix::Blosum62 => blosum62_score(x, y),
        }
    };

    let neg_inf = i32::MIN / 4;
    // Matrices: H=best, E=gap in A, F=gap in B
    let mut h = vec![vec![0i32; m+1]; n+1];
    let mut e = vec![vec![neg_inf; m+1]; n+1];
    let mut f = vec![vec![neg_inf; m+1]; n+1];
    // Traceback: 0=stop, 1=diag, 2=up (gap in B), 3=left (gap in A)
    let mut tb = vec![vec![0u8; m+1]; n+1];

    let mut max_i = 0usize;
    let mut max_j = 0usize;
    let mut max_score = 0i32;

    for i in 1..=n {
        for j in 1..=m {
            // extend or open gaps
            e[i][j] = (e[i][j-1] - ge).max(h[i][j-1] - go);
            f[i][j] = (f[i-1][j] - ge).max(h[i-1][j] - go);
            let diag = h[i-1][j-1] + score_pair(a[i-1], b[j-1]);
            let (best, dir) = [
                (0, 0u8),         // local alignment allows restart at 0
                (diag, 1u8),      // match/mismatch
                (f[i][j], 2u8),   // gap in B (move up)
                (e[i][j], 3u8),   // gap in A (move left)
            ].into_iter().max_by_key(|(v,_)| *v).unwrap();
            h[i][j] = best;
            tb[i][j] = dir;
            if best > max_score {
                max_score = best;
                max_i = i; max_j = j;
            }
        }
    }

    // Traceback from max cell
    let mut i = max_i;
    let mut j = max_j;
    let mut a_aln = Vec::new();
    let mut b_aln = Vec::new();
    let mut cigar_ops: Vec<(char, usize)> = Vec::new(); // M/I/D
    while i>0 && j>0 {
        match tb[i][j] {
            0 => break,
            1 => { // diag
                a_aln.push(a[i-1]); b_aln.push(b[j-1]); push_cigar(&mut cigar_ops, 'M', 1); i-=1; j-=1;
            }
            2 => { // up (gap in B)
                a_aln.push(a[i-1]); b_aln.push('-'); push_cigar(&mut cigar_ops, 'D', 1); i-=1;
            }
            3 => { // left (gap in A)
                a_aln.push('-'); b_aln.push(b[j-1]); push_cigar(&mut cigar_ops, 'I', 1); j-=1;
            }
            _ => break,
        }
    }
    a_aln.reverse(); b_aln.reverse();
    cigar_ops.reverse();
    let cigar = cigar_ops.into_iter().map(|(op,len)| format!("{len}{op}")).collect::<String>();

    // Compute ranges
    let start_a = i;
    let start_b = j;
    let end_a = max_i;
    let end_b = max_j;

    // Compute identity/gaps
    let (mut ident, mut gaps) = (0usize, 0usize);
    for (&x,&y) in a_aln.iter().zip(b_aln.iter()) {
        if x == '-' || y == '-' { gaps += 1; }
        else if equals_case_insensitive(x,y) { ident += 1; }
    }
    let cols = a_aln.len().max(1);
    let pct_identity = (ident as f64) * 100.0 / (cols as f64);
    let pct_gaps = (gaps as f64) * 100.0 / (cols as f64);

    Ok(WaterAlignment {
        score: max_score,
        range_a: (start_a, end_a),
        range_b: (start_b, end_b),
        align_a: a_aln.into_iter().collect(),
        align_b: b_aln.into_iter().collect(),
        cigar,
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
