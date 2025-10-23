//! All-vs-all local alignments (EMBOSS `seqmatchall`).
//!
//! This module runs declumped local alignments (Waterman–Eggert-style) for
//! **every unordered pair** of sequences in a single multi‑FASTA set, using
//! the [`matcher`](crate::matcher::matcher) engine under the hood.
//!
//! ### Examples
//! ```rust,no_run
//! use embossers::{parse_fasta, seqmatchall, SeqMatchAllParams, MatcherParams, WaterMatrix};
//! let recs = parse_fasta(r#">a\nACGT\n>b\nACGA\n"#);
//! let mparams = MatcherParams{ alternatives: 1, matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 }, ..Default::default() };
//! let params = SeqMatchAllParams{ matcher: mparams, include_self: false };
//! let results = seqmatchall(&recs, &params).unwrap();
//! assert_eq!(results.len(), 1);
//! ```
//!
//! ---
//! ## Details
//! - **Complexity:** Algorithms are O(n·m) for dynamic-programming steps; seed‑and‑extend
//!   approaches are subquadratic in practice by restricting extension windows.
//! - **Scoring:** DNA uses match/mismatch integers; proteins use BLOSUM62 via
//!   [`WaterMatrix::Blosum62`](crate::common::WaterMatrix). Gap penalties are affine
//!   with separate open/extend and optional integer scaling for fractional costs.
//! - **I/O:** FASTA parsing is minimal and permissive. Alignment text outputs are
//!   human‑readable and designed to resemble EMBOSS reports.
//! - **Parity:** Output formats aim to be familiar but are not byte‑identical to EMBOSS.
//!   Command‑line flags are close analogues. See README for examples and caveats.
//!

use crate::{FastaRecord, EmbossersError};
use crate::matcher::{matcher, MatcherParams};

/// Parameters for [`seqmatchall`].
#[derive(Clone, Debug)]
pub struct SeqMatchAllParams {
    /// Parameters forwarded to the local-alignment engine (`matcher`).
    pub matcher: MatcherParams,
    /// Include self-vs-self comparisons (`i == j`), default false.
    pub include_self: bool,
}

impl Default for SeqMatchAllParams {
    fn default() -> Self {
        Self {
            matcher: MatcherParams::default(),
            include_self: false,
        }
    }
}

/// One summary entry for a pair (i, j).
#[derive(Clone, Debug)]
pub struct SeqMatchAllEntry {
    /// Left/right record IDs.
    pub a_id: String,
    pub b_id: String,
    /// Best alignment score and metrics (from the top hit).
    pub score: i32,
    pub pct_identity: f64,
    pub pct_gaps: f64,
    /// Number of alternative local alignments found for this pair.
    pub n_hits: usize,
    /// CIGAR length of the top hit (proxy for alignment span).
    pub top_cigar_len: usize,
}

/// Compute local alignments for all pairs in `seqs` (unordered pairs if
/// `include_self=false`; include `i==j` when true). Returns a flat vector of
/// summary entries (row-major by i then j>i).
pub fn seqmatchall(seqs: &[FastaRecord], params: &SeqMatchAllParams) -> Result<Vec<SeqMatchAllEntry>, EmbossersError> {
    if seqs.is_empty() { return Err(EmbossersError::InvalidSequence("empty input set")); }
    let mut out = Vec::new();
    for i in 0..seqs.len() {
        let j_start = if params.include_self { i } else { i+1 };
        if j_start >= seqs.len() { continue; }
        for j in j_start..seqs.len() {
            let hits = matcher(&seqs[i].seq, &seqs[j].seq, &params.matcher)?;
            if hits.is_empty() {
                out.push(SeqMatchAllEntry{
                    a_id: seqs[i].id.clone(), b_id: seqs[j].id.clone(),
                    score: 0, pct_identity: 0.0, pct_gaps: 0.0, n_hits: 0, top_cigar_len: 0,
                });
            } else {
                let h = &hits[0];
                out.push(SeqMatchAllEntry{
                    a_id: seqs[i].id.clone(), b_id: seqs[j].id.clone(),
                    score: h.score, pct_identity: h.pct_identity, pct_gaps: h.pct_gaps,
                    n_hits: hits.len(), top_cigar_len: h.cigar.len(),
                });
            }
        }
    }
    Ok(out)
}