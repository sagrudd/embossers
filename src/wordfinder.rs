//! Query‑vs‑set local alignments (EMBOSS `wordfinder`‑like).
//!
//! Runs [`supermatcher`](crate::supermatcher::supermatcher) for a **single query**
//! against each sequence in a set, producing top local alignments per target.
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
use crate::supermatcher::{supermatcher, SuperMatcherParams};
use crate::matcher::MatcherHit;

/// Run `supermatcher` for a query against a set of target sequences.
pub fn wordfinder(query: &FastaRecord, targets: &[FastaRecord], params: &SuperMatcherParams)
    -> Result<Vec<(String, Vec<MatcherHit>)>, EmbossersError>
{
    let mut out = Vec::new();
    for t in targets {
        let hits = supermatcher(&query.seq, &t.seq, params)?;
        out.push((t.id.clone(), hits));
    }
    Ok(out)
}