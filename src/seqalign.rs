//! Extend an existing alignment with new sequences (EMBOSS `seqalign`‑like).
//!
//! This provides a pragmatic substitute for EMBASSY **seqalign**: take a base
//! multiple alignment (EMBOSS "simple" format or FASTA with gaps), compute a
//! **consensus**, then align each new sequence globally to the consensus and
//! insert columns accordingly to produce an **extended alignment**.
//!
//! This is not a structural alignment tool; it’s a fast sequence‑only stand‑in
//! suitable for everyday workflows or CI. It preserves all original columns and
//! adds gaps to fit the new sequence while trying to maintain overall column
//! homology.
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

use crate::common::{EmbossersError, WaterMatrix};
use crate::needle::{needle, NeedleParams};

/// Extend an alignment (string with gapped sequences) by aligning a new sequence
/// to the alignment consensus and inserting columns as needed.
pub fn seqalign_extend(base_alignment: &str, new_seq_id: &str, new_seq: &str, matrix: &WaterMatrix) -> Result<String, EmbossersError> {
    let (ids, rows) = parse_simple_alignment(base_alignment)?;
    let consensus = consensus_of_rows(&rows);
    let np = NeedleParams{ matrix: matrix.clone(), gap_open: 10.0, gap_extend: 0.5, scale: 2.0 };
    let aln = needle(&consensus, new_seq, &np)?;
    // Integrate new sequence according to alignment of consensus vs new_seq
    let extended = integrate_row(&ids, &rows, &aln.align_a, &aln.align_b, new_seq_id, new_seq);
    Ok(extended)
}

fn parse_simple_alignment(s: &str) -> Result<(Vec<String>, Vec<String>), EmbossersError> {
    let mut ids = Vec::new();
    let mut rows = Vec::new();
    for line in s.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') { continue; }
        if let Some((id, seq)) = t.split_once(' ') {
            ids.push(id.to_string());
            rows.push(seq.trim().to_string());
        } else if t.starts_with('>') {
            // FASTA with gaps
            // very small parser: capture id and next line as gapped row
            let id = t.trim_start_matches('>').split_whitespace().next().unwrap_or("seq").to_string();
            ids.push(id);
        } else {
            if let Some(last) = rows.last_mut() { last.push_str(t) } else { rows.push(t.to_string()); }
        }
    }
    if ids.len() != rows.len() { return Err(EmbossersError::InvalidSequence("could not parse simple alignment")); }
    Ok((ids, rows))
}

fn consensus_of_rows(rows: &[String]) -> String {
    if rows.is_empty() { return String::new(); }
    let cols = rows[0].len();
    let mut out = String::with_capacity(cols);
    for c in 0..cols {
        let mut counts = [0usize; 26];
        let mut best = ('-', 0usize);
        for r in rows {
            let ch = r.chars().nth(c).unwrap_or('-');
            if ch == '-' { continue; }
            let idx = (ch.to_ascii_uppercase() as u8).saturating_sub(b'A') as usize;
            if idx < 26 { counts[idx]+=1; if counts[idx] > best.1 { best=(ch, counts[idx]); } }
        }
        out.push(if best.1==0 { '-' } else { best.0 });
    }
    out
}

fn integrate_row(ids: &[String], rows: &[String], cons_aln: &str, new_aln: &str, new_id: &str, _new_seq: &str) -> String {
    // cons_aln/new_aln are aligned strings (with gaps). We need to inject gaps
    // into all existing rows where cons_aln has '-' so that new_aln fits column‑wise.
    let ca: Vec<char> = cons_aln.chars().collect();
    let nb: Vec<char> = new_aln.chars().collect();
    // Build insertion columns
    let mut insert_cols: Vec<usize> = Vec::new();
    for (i, &ch) in ca.iter().enumerate() {
        if ch == '-' { insert_cols.push(i); }
    }
    // Build extended rows by inserting '-' at each insert column
    let mut extended_rows: Vec<String> = Vec::new();
    for r in rows {
        let mut out: Vec<char> = Vec::new();
        let src: Vec<char> = r.chars().collect();
        let mut si = 0usize;
        for i in 0..ca.len() {
            if insert_cols.binary_search(&i).is_ok() {
                out.push('-');
            } else {
                out.push(*src.get(si).unwrap_or(&'-'));
                si+=1;
            }
        }
        extended_rows.push(out.iter().collect());
    }
    // Build the new sequence row by removing gaps from new_aln (should be new_seq with '-'s)
    let new_row: String = nb.iter().collect();
    // Compose "simple" format (id + space + row)
    let mut s = String::new();
    for (id, row) in ids.iter().zip(extended_rows.iter()) {
        s.push_str(id); s.push(' '); s.push_str(row); s.push('\n');
    }
    s.push_str(new_id); s.push(' '); s.push_str(&new_row); s.push('\n');
    s
}