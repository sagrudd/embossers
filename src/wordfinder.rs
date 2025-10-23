//! Query‑vs‑set local alignments (EMBOSS `wordfinder`‑like).
//!
//! Runs [`supermatcher`](crate::supermatcher::supermatcher) for a **single query**
//! against each sequence in a set, producing top local alignments per target.
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
