//! Exact word matches (EMBOSS `wordmatch`).
//!
//! Finds all **exact** matches of at least `wordlen` between two sequences by
//! k‑mer seeding and right‑extension to maximal runs. Outputs regions in both
//! sequences and optional GFF3 lines can be produced by the CLI layer.
//!
use std::collections::HashMap;
use crate::common::EmbossersError;

/// Parameters for exact matching.
#[derive(Clone, Debug)]
pub struct WordMatchParams {
    /// Minimum exact word length.
    pub wordlen: usize,
}

impl Default for WordMatchParams {
    fn default() -> Self { Self { wordlen: 10 } }
}

/// One exact match region.
#[derive(Clone, Debug)]
pub struct WordMatch {
    /// 0‑based inclusive coordinates in A.
    pub a_start: usize, pub a_end: usize,
    /// 0‑based inclusive coordinates in B.
    pub b_start: usize, pub b_end: usize,
    /// Length of the exact match.
    pub len: usize,
    /// The exact matching substring.
    pub seq: String,
}

/// Find all exact match regions of length ≥ `wordlen`.
pub fn wordmatch(a: &str, b: &str, params: &WordMatchParams) -> Result<Vec<WordMatch>, EmbossersError> {
    if a.len() < params.wordlen || b.len() < params.wordlen {
        return Ok(Vec::new());
    }
    // Index all k‑mers in A
    let mut index: HashMap<&str, Vec<usize>> = HashMap::new();
    for i in 0..=a.len()-params.wordlen {
        let k = &a[i..i+params.wordlen];
        index.entry(k).or_default().push(i);
    }
    // Scan B and extend matches
    let mut out: Vec<WordMatch> = Vec::new();
    for j in 0..=b.len()-params.wordlen {
        let k = &b[j..j+params.wordlen];
        if let Some(pos) = index.get(k) {
            for &i in pos {
                // Extend right to maximal run
                let mut len = params.wordlen;
                while i+len < a.len() && j+len < b.len() && &a[i+len..i+len+1] == &b[j+len..j+len+1] {
                    len += 1;
                }
                let seq = a[i..i+len].to_string();
                out.push(WordMatch{ a_start: i, a_end: i+len-1, b_start: j, b_end: j+len-1, len, seq });
            }
        }
    }
    // Optionally deduplicate identical segments
    out.sort_by_key(|w| (w.a_start, w.b_start, w.len));
    out.dedup_by(|x,y| x.a_start==y.a_start && x.b_start==y.b_start && x.len==y.len);
    Ok(out)
}
