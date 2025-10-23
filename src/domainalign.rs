//! Progressive, consensus-guided multiple alignment (EMBASSY `domainalign`‑like).
//!
//! This is a pragmatic stand‑in for the `domalign` tools: it builds an alignment
//! by iteratively aligning each input sequence **globally** to the **current
//! consensus** and inserting columns as needed, preserving existing columns.
//!
//! For robustness it reuses the library's `needle` (Needleman–Wunsch) engine.
//!
//! ### Example
//! ```rust,no_run
//! use embossers::{parse_fasta, domainalign, DomainAlignParams, WaterMatrix};
//! let recs = parse_fasta(">a\nACGT\n>b\nACGA\n>c\nACG-\n");
//! let p = DomainAlignParams{ matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 }, ..Default::default() };
//! let daf = domainalign(&recs, &p).unwrap();
//! assert!(daf.contains("a "));
//! ```
//!
use crate::{FastaRecord};
use crate::common::{EmbossersError, WaterMatrix};
use crate::needle::{needle, NeedleParams};

/// Parameters for [`domainalign`].
#[derive(Clone, Debug)]
pub struct DomainAlignParams {
    /// Pairwise scoring matrix.
    pub matrix: WaterMatrix,
    /// Gap open penalty.
    pub gap_open: f32,
    /// Gap extension penalty.
    pub gap_extend: f32,
    /// Integer scaling for fractional penalties.
    pub scale: f32,
}

impl Default for DomainAlignParams {
    fn default() -> Self {
        Self {
            matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 },
            gap_open: 10.0,
            gap_extend: 0.5,
            scale: 2.0,
        }
    }
}

/// Build a "simple" alignment ("id<space>row\n") from sequences by iteratively
/// aligning each to the running consensus using global alignment.
pub fn domainalign(seqs: &[FastaRecord], params: &DomainAlignParams) -> Result<String, EmbossersError> {
    if seqs.is_empty() { return Err(EmbossersError::InvalidSequence("empty input")); }
    // Start with first seq as the initial alignment (single row)
    let mut ids: Vec<String> = vec![seqs[0].id.clone()];
    let mut rows: Vec<String> = vec![seqs[0].seq.clone()];
    for rec in &seqs[1..] {
        // Build consensus from current rows
        let cons = consensus_of_rows(&rows);
        // Align consensus vs new sequence
        let np = NeedleParams{ matrix: params.matrix.clone(), gap_open: params.gap_open, gap_extend: params.gap_extend, scale: params.scale };
        let aln = needle(&cons, &rec.seq, &np)?;
        // Insert columns according to gaps in consensus alignment, then append the new row
        let (new_rows, new_row) = integrate_rows(&rows, &aln.align_a, &aln.align_b, &rec.seq);
        rows = new_rows;
        rows.push(new_row);
        ids.push(rec.id.clone());
    }
    // Compose simple alignment
    let mut out = String::new();
    for (id,row) in ids.iter().zip(rows.iter()) {
        out.push_str(id); out.push(' '); out.push_str(row); out.push('\n');
    }
    Ok(out)
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

fn integrate_rows(rows: &[String], cons_aln: &str, new_aln: &str, new_seq: &str) -> (Vec<String>, String) {
    // Insert gap columns for all rows where cons_aln has '-'.
    let ca: Vec<char> = cons_aln.chars().collect();
    let nb: Vec<char> = new_aln.chars().collect();
    let mut insert_cols: Vec<usize> = Vec::new();
    for (i, &ch) in ca.iter().enumerate() { if ch=='-' { insert_cols.push(i); } }

    let mut new_rows: Vec<String> = Vec::new();
    for r in rows {
        let src: Vec<char> = r.chars().collect();
        let mut out: Vec<char> = Vec::with_capacity(ca.len());
        let mut si = 0usize;
        for i in 0..ca.len() {
            if insert_cols.binary_search(&i).is_ok() { out.push('-'); } else { out.push(*src.get(si).unwrap_or(&'-')); si+=1; }
        }
        new_rows.push(out.iter().collect());
    }
    let new_row: String = nb.iter().collect();
    (new_rows, new_row)
}
