//! Within-sequence repeat finder (EMBASSY `domainrep`‑like).
//!
//! Detects repeated blocks in a single sequence by **k‑mer seeding** followed by
//! **maximal right extension**. This reports coordinates of repeated segments and
//! can serve as a quick domain repeat discovery tool.
//!
//! ### Example
//! ```rust,no_run
//! use embossers::{domainrep, DomainRepParams};
//! let s = "ABCDABCDXYZ";
//! let p = DomainRepParams{ wordlen: 2, min_len: 4 };
//! let reps = domainrep(s, &p).unwrap();
//! assert!(!reps.is_empty());
//! ```
//!
use std::collections::HashMap;
use crate::common::EmbossersError;

/// Parameters for [`domainrep`].
#[derive(Clone, Debug)]
pub struct DomainRepParams {
    /// Seed size (k-mer length).
    pub wordlen: usize,
    /// Minimum repeated block length to report.
    pub min_len: usize,
}

impl Default for DomainRepParams {
    fn default() -> Self { Self { wordlen: 3, min_len: 12 } }
}

/// A repeated block in the sequence.
#[derive(Clone, Debug)]
pub struct RepeatBlock {
    /// First copy coordinates [start,end] inclusive.
    pub a_start: usize, pub a_end: usize,
    /// Second copy coordinates [start,end] inclusive.
    pub b_start: usize, pub b_end: usize,
    /// Length of the repeat.
    pub len: usize,
}

/// Find repeated blocks using k‑mer seeding and maximal right extension.
pub fn domainrep(seq: &str, params: &DomainRepParams) -> Result<Vec<RepeatBlock>, EmbossersError> {
    if seq.len() < params.wordlen { return Ok(Vec::new()); }
    let mut index: HashMap<&str, Vec<usize>> = HashMap::new();
    for i in 0..=seq.len() - params.wordlen {
        index.entry(&seq[i..i+params.wordlen]).or_default().push(i);
    }
    let mut out: Vec<RepeatBlock> = Vec::new();
    for j in 0..=seq.len() - params.wordlen {
        if let Some(list) = index.get(&seq[j..j+params.wordlen]) {
            for &i in list {
                if i >= j { continue; } // avoid duplicates and self
                // extend right
                let mut len = params.wordlen;
                while i+len < seq.len() && j+len < seq.len() && &seq[i+len..i+len+1] == &seq[j+len..j+len+1] {
                    len += 1;
                }
                if len >= params.min_len {
                    out.push(RepeatBlock{ a_start: i, a_end: i+len-1, b_start: j, b_end: j+len-1, len });
                }
            }
        }
    }
    // Merge duplicates
    out.sort_by_key(|r| (r.a_start, r.b_start, r.len));
    out.dedup_by(|x,y| x.a_start==y.a_start && x.b_start==y.b_start && x.len==y.len);
    Ok(out)
}
