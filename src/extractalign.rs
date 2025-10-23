//! Extract alignment columns (EMBOSS `extractalign`‑like).
//!
//! Given a gapped alignment in **FASTA** or **simple** (`id<space>row`) format,
//! extract one or more **1‑based inclusive** column ranges and emit a new simple
//! alignment limited to those columns (in order).
//!
//! ### Example
//! ```rust
//! use embossers::extractalign_ranges;
//! let aln = "a A-CGT\nb ATCG-\n";
//! let out = extractalign_ranges(aln, &[(2,4)]).unwrap();
//! assert!(out.contains("a CGT"));
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

use crate::common::EmbossersError;

/// Parse an alignment (gapped FASTA or simple) into `(ids, rows)`.
fn parse_alignment(s: &str) -> Result<(Vec<String>, Vec<String>), EmbossersError> {
    let t = s.trim_start();
    if t.starts_with('>') {
        let recs = crate::parse_fasta(s);
        let ids = recs.iter().map(|r| r.id.clone()).collect::<Vec<_>>();
        let rows = recs.iter().map(|r| r.seq.clone()).collect::<Vec<_>>();
        return Ok((ids, rows));
    }
    let mut ids = Vec::new();
    let mut rows = Vec::new();
    for line in s.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        if let Some((id, row)) = line.split_once(' ') {
            ids.push(id.to_string());
            rows.push(row.trim().to_string());
        }
    }
    if ids.is_empty() { return Err(EmbossersError::InvalidSequence("could not parse alignment")); }
    Ok((ids, rows))
}

/// Extract the specified **1‑based inclusive** column ranges from `aln_text` and
/// return a new alignment in **simple** format (`id<space>row\n`).
pub fn extractalign_ranges(aln_text: &str, ranges: &[(usize,usize)]) -> Result<String, EmbossersError> {
    let (ids, rows) = parse_alignment(aln_text)?;
    let cols = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    // Normalize and clamp ranges
    let mut rgs: Vec<(usize,usize)> = Vec::new();
    for &(s,e) in ranges {
        if s==0 { continue; }
        let a = s.min(e).max(1);
        let b = s.max(e).min(cols);
        if a<=b { rgs.push((a,b)); }
    }
    if rgs.is_empty() { return Err(EmbossersError::InvalidSequence("no valid ranges")); }
    // Build output rows
    let mut out = String::new();
    for (i,row) in rows.iter().enumerate() {
        let chars: Vec<char> = row.chars().collect();
        let mut newrow: Vec<char> = Vec::new();
        for (a,b) in &rgs {
            let a0 = a-1;
            let b0 = *b;
            for c in a0..b0 {
                newrow.push(*chars.get(c).unwrap_or(&'-'));
            }
        }
        out.push_str(&ids[i]); out.push(' '); out.push_str(&newrow.iter().collect::<String>()); out.push('\n');
    }
    Ok(out)
}