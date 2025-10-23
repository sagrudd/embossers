//! Global alignment using Myers–Miller (linear space) a.k.a. EMBOSS `stretcher`.
//!
//! This module implements a divide‑and‑conquer Needleman–Wunsch with **affine gaps**
//! and **linear memory** inspired by Myers–Miller. For small problems we fall back
//! to the full‑matrix implementation in [`crate::needle`].
//!
//! ### When to use
//! Use `stretcher` when sequences are long and a full DP matrix would be too
//! memory‑hungry; time is still `O(n*m)` but memory stays near `O(n+m)`.
//!
//! ### Example
//! ```rust,no_run
//! use embossers::{stretcher, StretcherParams, WaterMatrix};
//! let p = StretcherParams{ matrix: WaterMatrix::Dna{ match_score: 1, mismatch: -1 }, ..Default::default() };
//! let aln = stretcher("GATTACA", "GCATGCU", &p).unwrap();
//! assert_eq!(aln.align_a.len(), aln.align_b.len());
//! ```
//!
use crate::common::{EmbossersError, WaterMatrix, blosum62_score, equals_case_insensitive};
use crate::needle::{needle, NeedleParams};

/// Parameters for `stretcher` global alignment.
#[derive(Clone, Debug)]
pub struct StretcherParams {
    /// Scoring matrix to use.
    pub matrix: WaterMatrix,
    /// Gap open penalty (can be fractional; default 10.0 to mimic EMBOSS).
    pub gap_open: f32,
    /// Gap extension penalty (can be fractional; default 0.5 to mimic EMBOSS).
    pub gap_extend: f32,
    /// Integer scale factor to preserve fractional penalties (e.g. 2.0 so 0.5→1).
    pub scale: f32,
    /// Size threshold to switch to full DP (smaller is faster recursion).
    pub fallback_threshold: usize,
}

impl Default for StretcherParams {
    fn default() -> Self {
        Self {
            matrix: WaterMatrix::Dna { match_score: 1, mismatch: -1 },
            gap_open: 10.0,
            gap_extend: 0.5,
            scale: 2.0,
            fallback_threshold: 256, // switch to full DP when n*m <= this
        }
    }
}

/// Output of a `stretcher` global alignment.
#[derive(Clone, Debug)]
pub struct StretcherAlignment {
    /// Total alignment **score** (integer, scaled as provided).
    pub score: i32,
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

/// Run a Myers–Miller/Hirschberg global alignment with affine gaps in linear space.
pub fn stretcher(a: &str, b: &str, params: &StretcherParams) -> Result<StretcherAlignment, EmbossersError> {
    if a.is_empty() || b.is_empty() {
        return Err(EmbossersError::InvalidSequence("empty sequence"));
    }
    // Prepare scoring
    let scale = params.scale.max(1.0);
    let go = (params.gap_open * scale).round() as i32;
    let ge = (params.gap_extend * scale).round() as i32;
    let score_pair = |x: char, y: char| -> i32 {
        match &params.matrix {
            WaterMatrix::Dna{match_score, mismatch} => if x == y { *match_score } else { *mismatch },
            WaterMatrix::Blosum62 => blosum62_score(x, y),
        }
    };

    // Fallback to full DP for small problems: reuse `needle` implementation.
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();
    if a_chars.len() * b_chars.len() <= params.fallback_threshold {
        let np = NeedleParams{ matrix: params.matrix.clone(), gap_open: params.gap_open, gap_extend: params.gap_extend, scale: params.scale };
        let n = needle(a, b, &np)?;
        return Ok(StretcherAlignment{
            score: n.score,
            align_a: n.align_a,
            align_b: n.align_b,
            cigar: n.cigar,
            pct_identity: n.pct_identity,
            pct_gaps: n.pct_gaps,
        });
    }

    // Divide-and-conquer alignment strings
    let (aln_a, aln_b) = hirschberg_affine(&a_chars, &b_chars, go, ge, &score_pair, params.fallback_threshold);

    // Score and stats
    let mut cigar_ops: Vec<(char, usize)> = Vec::new();
    let mut ident: usize = 0;
    let mut gaps: usize = 0;
    for (&x,&y) in aln_a.iter().zip(aln_b.iter()) {
        if x=='-' && y=='-' { continue; } // shouldn't happen
        if x=='-' || y=='-' {
            gaps += 1;
            // Score gaps using affine accounting
            // We need to add open+extend runs; build CIGAR to do both tasks once.
            push_cigar(&mut cigar_ops, if x=='-' { 'I' } else { 'D' }, 1);
        } else {
            if equals_case_insensitive(x,y) { ident += 1; }
            push_cigar(&mut cigar_ops, 'M', 1);
        }
    }
    // Compute score from CIGAR (to avoid double-charging opens during traversal)
    let score = score_from_cigar(&aln_a, &aln_b, go, ge, &score_pair);
    let cigar = cigar_ops.into_iter().map(|(op,len)| format!("{len}{op}")).collect::<String>();
    let cols = aln_a.len().max(1);
    let pct_identity = (ident as f64)*100.0/(cols as f64);
    let pct_gaps = (gaps as f64)*100.0/(cols as f64);

    Ok(StretcherAlignment{
        score,
        align_a: aln_a.into_iter().collect(),
        align_b: aln_b.into_iter().collect(),
        cigar,
        pct_identity,
        pct_gaps,
    })
}

fn score_from_cigar<F>(aln_a: &[char], aln_b: &[char], go: i32, ge: i32, s: &F) -> i32
where F: Fn(char,char)->i32 {
    let mut sc = 0i32;
    let mut i = 0usize;
    let n = aln_a.len();
    while i < n {
        if aln_a[i] == '-' && aln_b[i] == '-' { i+=1; continue; }
        if aln_a[i] == '-' || aln_b[i] == '-' {
            // run of gaps
            let mut len = 0usize;
            let is_ins = aln_a[i] == '-';
            while i < n && (aln_a[i] == '-' || aln_b[i] == '-') && (aln_a[i] == '-' ) == is_ins {
                len+=1; i+=1;
            }
            sc -= go + (len as i32 -1).max(0)*ge;
        } else {
            sc += s(aln_a[i], aln_b[i]);
            i+=1;
        }
    }
    sc
}

fn push_cigar(ops: &mut Vec<(char,usize)>, op: char, k: usize) {
    if let Some(last) = ops.last_mut() { if last.0 == op { last.1 += k; return; } }
    ops.push((op, k));
}

fn hirschberg_affine<F>(a: &[char], b: &[char], go: i32, ge: i32, s: &F, thresh: usize) -> (Vec<char>, Vec<char>)
where F: Fn(char,char)->i32 + Copy {
    let n = a.len();
    let m = b.len();
    if n==0 {
        return (vec!['-'; m], b.to_vec());
    }
    if m==0 {
        return (a.to_vec(), vec!['-'; n]);
    }
    if n * m <= thresh {
        // Small: full DP via needle to retrieve alignment strings.
        // We cannot call needle directly with arbitrary s; so we inline a small DP here.
        return full_dp_affine(a, b, go, ge, s);
    }
    let mid = n/2;
    let (f_m, f_x, f_y) = linear_forward(&a[..mid], b, go, ge, s);
    let (b_m, b_x, b_y) = linear_backward(&a[mid..], b, go, ge, s);
    // Choose split j maximizing combined score; subtract `go` if continuing X or Y past the split
    let mut best = i32::MIN/4;
    let mut split = 0usize;
    for j in 0..=m {
        let v = max3(f_m[j] + b_m[m-j], f_x[j] + b_x[m-j] - go, f_y[j] + b_y[m-j] - go);
        if v > best { best = v; split = j; }
    }
    let (a1,b1) = hirschberg_affine(&a[..mid], &b[..split], go, ge, s, thresh);
    let (a2,b2) = hirschberg_affine(&a[mid..], &b[split..], go, ge, s, thresh);
    let mut aa = a1; aa.extend(a2);
    let mut bb = b1; bb.extend(b2);
    (aa, bb)
}

fn max3(a:i32,b:i32,c:i32)->i32 { a.max(b).max(c) }

fn linear_forward<F>(a: &[char], b: &[char], go: i32, ge: i32, s: &F) -> (Vec<i32>, Vec<i32>, Vec<i32>)
where F: Fn(char,char)->i32 + Copy {
    let n = a.len();
    let m = b.len();
    let neg_inf = i32::MIN/4;
    let mut m_prev = vec![neg_inf; m+1];
    let mut x_prev = vec![neg_inf; m+1];
    let mut y_prev = vec![neg_inf; m+1];
    m_prev[0] = 0;
    // row 0 initialization for X
    for j in 1..=m {
        x_prev[j] = -go - (j as i32 -1)*ge;
        // m_prev[j] and y_prev[j] already -inf
    }
    for i in 1..=n {
        let mut m_cur = vec![neg_inf; m+1];
        let mut x_cur = vec![neg_inf; m+1];
        let mut y_cur = vec![neg_inf; m+1];
        // column 0
        y_cur[0] = -go - (i as i32 -1)*ge;
        // m_cur[0], x_cur[0] stay -inf
        for j in 1..=m {
            // X: gap in A (left)
            x_cur[j] = (m_cur[j-1] - go).max(x_cur[j-1] - ge);
            // Y: gap in B (up)
            y_cur[j] = (m_prev[j] - go).max(y_prev[j] - ge);
            // M: diag
            let diag_base = m_prev[j-1].max(x_prev[j-1]).max(y_prev[j-1]);
            m_cur[j] = diag_base + s(a[i-1], b[j-1]);
        }
        m_prev = m_cur; x_prev = x_cur; y_prev = y_cur;
    }
    (m_prev, x_prev, y_prev)
}

fn linear_backward<F>(a: &[char], b: &[char], go: i32, ge: i32, s: &F) -> (Vec<i32>, Vec<i32>, Vec<i32>)
where F: Fn(char,char)->i32 + Copy {
    // compute on reversed suffix to get costs from (mid, j) to end
    let ar: Vec<char> = a.iter().rev().copied().collect();
    let br: Vec<char> = b.iter().rev().copied().collect();
    let (m, x, y) = linear_forward(&ar, &br, go, ge, s);
    // m,x,y correspond to arrays at start of reversed -> which equals end of original
    (m.into_iter().rev().collect(), x.into_iter().rev().collect(), y.into_iter().rev().collect())
}

fn full_dp_affine<F>(a: &[char], b: &[char], go: i32, ge: i32, s: &F) -> (Vec<char>, Vec<char>)
where F: Fn(char,char)->i32 + Copy {
    // Small full DP with affine gaps; return alignment strings (global).
    let n = a.len();
    let m = b.len();
    let neg_inf = i32::MIN/4;
    // mat: best; gap_a: gap in A (left); gap_b: gap in B (up)
    let mut mat = vec![vec![neg_inf; m+1]; n+1];
    let mut gap_a = vec![vec![neg_inf; m+1]; n+1];
    let mut gap_b = vec![vec![neg_inf; m+1]; n+1];
    let mut tb  = vec![vec![0u8;    m+1]; n+1]; // 1=diag, 2=up(gap_b), 3=left(gap_a)

    mat[0][0] = 0;
    for j in 1..=m {
        gap_a[0][j] = -go - (j as i32 - 1)*ge;
        mat[0][j]   = gap_a[0][j];
        tb[0][j]    = 3;
    }
    for i in 1..=n {
        gap_b[i][0] = -go - (i as i32 - 1)*ge;
        mat[i][0]   = gap_b[i][0];
        tb[i][0]    = 2;
    }
    for i in 1..=n {
        for j in 1..=m {
            gap_a[i][j] = (mat[i][j-1] - go).max(gap_a[i][j-1] - ge);
            gap_b[i][j] = (mat[i-1][j] - go).max(gap_b[i-1][j] - ge);
            let diag_base = mat[i-1][j-1].max(gap_a[i-1][j-1]).max(gap_b[i-1][j-1]);
            let diag_score = diag_base + s(a[i-1], b[j-1]);
            let (best, dir) = [(diag_score,1u8),(gap_b[i][j],2u8),(gap_a[i][j],3u8)]
                .into_iter().max_by_key(|(v,_)| *v).unwrap();
            mat[i][j] = best;
            tb[i][j] = dir;
        }
    }
    // Traceback
    let mut i = n; let mut j = m;
    let mut aa: Vec<char> = Vec::new();
    let mut bb: Vec<char> = Vec::new();
    while i>0 || j>0 {
        if i>0 && j>0 && tb[i][j]==1 {
            aa.push(a[i-1]); bb.push(b[j-1]); i-=1; j-=1;
        } else if i>0 && (j==0 || tb[i][j]==2) {
            aa.push(a[i-1]); bb.push('-'); i-=1;
        } else {
            aa.push('-'); bb.push(b[j-1]); j-=1;
        }
    }
    aa.reverse(); bb.reverse();
    (aa, bb)
}
