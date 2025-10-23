//! Many-to-many global alignments (EMBOSS `needleall`).
//!
//! This module runs [`needle`](crate::needle::needle) across **all pairs** of
//! two FASTA sets and returns a flat list of results. It mirrors the spirit
//! of EMBOSS **needleall**: rigorous Needleman–Wunsch global alignment with
//! affine gaps, configurable matrix, and per-pair reporting.
//!
//! ### Examples
//! ```rust,no_run
//! use embossers::{parse_fasta, NeedleParams, WaterMatrix, needleall_pairs};
//! let a = parse_fasta(r#">a1\nACGT\n>a2\nACGA\n"#);
//! let b = parse_fasta(r#">b1\nACGG\n"#);
//! let p = NeedleParams{ matrix: WaterMatrix::Dna{ match_score: 1, mismatch: -1 }, ..Default::default() };
//! let results = needleall_pairs(&a, &b, &p).unwrap();
//! assert_eq!(results.len(), a.len()*b.len());
//! ```
//!
use crate::{NeedleParams, needle, EmbossersError, FastaRecord};

/// A single pairwise result from `needleall_pairs`.
#[derive(Clone, Debug)]
pub struct NeedleAllResult {
    /// Query (left) record ID.
    pub a_id: String,
    /// Subject (right) record ID.
    pub b_id: String,
    /// Score from the global alignment.
    pub score: i32,
    /// Percent identity over aligned columns.
    pub pct_identity: f64,
    /// Percent gaps over aligned columns.
    pub pct_gaps: f64,
    /// Aligned strings and CIGAR (optional to serialize/print).
    pub align_a: String,
    pub align_b: String,
    pub cigar: String,
}

/// Run Needleman–Wunsch across **all pairs** of `left × right` and return results
/// in row-major order (for each `a` in `left`, iterate all `b` in `right`).
pub fn needleall_pairs(left: &[FastaRecord], right: &[FastaRecord], params: &NeedleParams)
    -> Result<Vec<NeedleAllResult>, EmbossersError>
{
    if left.is_empty() || right.is_empty() {
        return Err(EmbossersError::InvalidSequence("empty input set"));
    }
    let mut out = Vec::with_capacity(left.len()*right.len());
    for a in left {
        for b in right {
            let aln = needle(&a.seq, &b.seq, params)?;
            out.push(NeedleAllResult{
                a_id: a.id.clone(),
                b_id: b.id.clone(),
                score: aln.score,
                pct_identity: aln.pct_identity,
                pct_gaps: aln.pct_gaps,
                align_a: aln.align_a,
                align_b: aln.align_b,
                cigar: aln.cigar,
            });
        }
    }
    Ok(out)
}
