//! Consensus sequence utilities (EMBOSS `cons`/`consambig`‑like).
//!
//! - `consensus_majority`: majority rule consensus with a threshold (DNA/protein).
//! - `consensus_ambig_dna`: DNA consensus using IUPAC ambiguity codes.
//!
use crate::common::EmbossersError;

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

/// Majority rule consensus (like `cons`). If the top symbol's fraction is below
/// `threshold` (0.0..=1.0), emits `'N'` for DNA or `'X'` for protein. Gap-only
/// columns yield `'-'`.
pub fn consensus_majority(aln_text: &str, is_dna: bool, threshold: f64) -> Result<String, EmbossersError> {
    let (_ids, rows) = parse_alignment(aln_text)?;
    if rows.is_empty() { return Err(EmbossersError::InvalidSequence("empty alignment")); }
    let cols = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut out = String::with_capacity(cols);
    for c in 0..cols {
        let mut counts = [0usize; 26];
        let mut n_non_gap = 0usize;
        for r in &rows {
            let ch = r.chars().nth(c).unwrap_or('-');
            if ch=='-' { continue; }
            n_non_gap += 1;
            let u = ch.to_ascii_uppercase();
            let u = if is_dna && u=='U' { 'T' } else { u };
            let idx = (u as u8).saturating_sub(b'A') as usize;
            if idx < 26 { counts[idx]+=1; }
        }
        if n_non_gap==0 { out.push('-'); continue; }
        // pick best
        let mut best = ('N', 0usize);
        for i in 0..26 {
            if counts[i] > best.1 {
                best = ((b'A'+(i as u8)) as char, counts[i]);
            }
        }
        let frac = (best.1 as f64)/(n_non_gap as f64);
        if frac >= threshold {
            out.push(best.0);
        } else {
            out.push(if is_dna { 'N' } else { 'X' });
        }
    }
    Ok(out)
}

/// DNA consensus using IUPAC ambiguity symbols from all **observed** bases
/// (A,C,G,T/U). Gap-only columns emit `'-'`. If no DNA symbols present, `N`.
pub fn consensus_ambig_dna(aln_text: &str) -> Result<String, EmbossersError> {
    let (_ids, rows) = parse_alignment(aln_text)?;
    if rows.is_empty() { return Err(EmbossersError::InvalidSequence("empty alignment")); }
    let cols = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut out = String::with_capacity(cols);
    for c in 0..cols {
        let mut a=c==0; // dummy init to avoid warnings
        let mut seen = [false; 4]; // A C G T
        for r in &rows {
            let ch = r.chars().nth(c).unwrap_or('-');
            if ch=='-' { continue; }
            let u = ch.to_ascii_uppercase();
            match u {
                'A' => seen[0]=true,
                'C' => seen[1]=true,
                'G' => seen[2]=true,
                'T' | 'U' => seen[3]=true,
                _ => {}
            }
        }
        let code = iupac_from_seen(seen);
        out.push(code);
    }
    Ok(out)
}

fn iupac_from_seen(seen: [bool;4]) -> char {
    let a=seen[0]; let c=seen[1]; let g=seen[2]; let t=seen[3];
    match (a,c,g,t) {
        (false,false,false,false) => '-',
        (true,false,false,false) => 'A',
        (false,true,false,false) => 'C',
        (false,false,true,false) => 'G',
        (false,false,false,true) => 'T',
        (true,false,true,false) => 'R',  // A/G
        (false,true,false,true) => 'Y',  // C/T
        (false,true,true,false) => 'S',  // C/G
        (true,false,false,true) => 'W',  // A/T
        (false,false,true,true) => 'K',  // G/T
        (true,true,false,false) => 'M',  // A/C
        (false,true,true,true) => 'B',   // C/G/T
        (true,false,true,true) => 'D',   // A/G/T
        (true,true,false,true) => 'H',   // A/C/T
        (true,true,true,false) => 'V',   // A/C/G
        (true,true,true,true) => 'N',
    }
}
