//! Alignment statistics (EMBASSY `oalistat`‑like).
//!
//! Reads a gapped FASTA or "simple" alignment and computes summary metrics:
//! number of sequences, columns, gap fraction, mean pairwise identity
//! (ignoring gap–gap), and count of conserved columns.
//!
//! ### Example
//! ```rust,no_run
//! use embossers::oalistat;
//! let s = ">a\nA-CG\n>b\nATCG\n";
//! let stats = oalistat(&s).unwrap();
//! assert_eq!(stats.nseqs, 2);
//! ```
//!
use crate::common::EmbossersError;

/// Summary statistics for a multiple sequence alignment.
#[derive(Clone, Debug)]
pub struct OaStats {
    /// Number of sequences in the alignment.
    pub nseqs: usize,
    /// Number of alignment columns (max row length).
    pub cols: usize,
    /// Fraction of columns that are all gaps.
    pub gap_columns_fraction: f64,
    /// Mean pairwise identity across all pairs (excluding gap-gap).
    pub mean_pairwise_identity: f64,
    /// Number of strictly conserved columns (all non-gaps are identical).
    pub conserved_columns: usize,
}

/// Compute statistics for a gapped alignment given as a string (FASTA-with-gaps
/// or simple "id<space>row").
pub fn oalistat(s: &str) -> Result<OaStats, EmbossersError> {
    let (ids, rows) = parse_alignment(s)?;
    if rows.is_empty() { return Err(EmbossersError::InvalidSequence("empty alignment")); }
    let n = rows.len();
    let cols = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut gap_cols = 0usize;
    let mut conserved = 0usize;

    // pairwise identity
    let mut pairs = 0usize;
    let mut ident_sum = 0usize;

    for c in 0..cols {
        let mut all_gap = true;
        let mut symbol: Option<char> = None;
        let mut same = true;
        for r in &rows {
            let ch = r.chars().nth(c).unwrap_or('-');
            if ch != '-' {
                all_gap = false;
                if let Some(s) = symbol {
                    if !eq_case_insensitive(s, ch) { same = false; }
                } else {
                    symbol = Some(ch);
                }
            }
        }
        if all_gap { gap_cols += 1; }
        if same && !all_gap { conserved += 1; }
    }

    // Mean pairwise identity
    for i in 0..n {
        for j in i+1..n {
            let (mut id, mut denom) = (0usize, 0usize);
            let a: Vec<char> = rows[i].chars().collect();
            let b: Vec<char> = rows[j].chars().collect();
            let m = cols;
            for c in 0..m {
                let x = *a.get(c).unwrap_or(&'-');
                let y = *b.get(c).unwrap_or(&'-');
                if x=='-' && y=='-' { continue; }
                denom += 1;
                if x!='-' && y!='-' && eq_case_insensitive(x,y) { id += 1; }
            }
            if denom>0 { ident_sum += (id*1000)/denom; pairs += 1; }
        }
    }
    let mpi = if pairs>0 { (ident_sum as f64)/(pairs as f64)/10.0 } else { 0.0 };
    let gap_frac = if cols>0 { (gap_cols as f64)/(cols as f64) } else { 0.0 };

    Ok(OaStats{
        nseqs: n, cols,
        gap_columns_fraction: gap_frac,
        mean_pairwise_identity: mpi,
        conserved_columns: conserved,
    })
}

fn eq_case_insensitive(a: char, b: char) -> bool {
    a.to_ascii_uppercase() == b.to_ascii_uppercase()
}

fn parse_alignment(s: &str) -> Result<(Vec<String>, Vec<String>), EmbossersError> {
    let t = s.trim_start();
    if t.starts_with('>') {
        // FASTA with gaps
        let recs = crate::parse_fasta(s);
        let ids = recs.iter().map(|r| r.id.clone()).collect::<Vec<_>>();
        let rows = recs.iter().map(|r| r.seq.clone()).collect::<Vec<_>>();
        return Ok((ids, rows));
    }
    // Simple: lines "id<space>ROW"
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
