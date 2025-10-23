//! SIM4-like spliced alignment (EMBOSS `esim4`).
#![allow(clippy::needless_return)]
//!
//! This module wraps [`crate::est2genome::est2genome`] && adds:
//!
//! - **Strand autodetection**: try EST forward && reverse-complement; keep best.
//! - **Splice classification**: canonical `GT-AG`, plus optional `GC-AG`, `AT-AC`.
////! - **Exon table**: derived from intron breaks; genomic coordinates reported.
//!
//! The goal is practical compatibility with EMBOSS `esim4` behaviour, not an
//! exact port of the historic SIM4 codebase. Parameters are intentionally
//! conservative && mirror common defaults.
//!
//! ### Example
//! ```rust,no_run
//! use embossers::{esim4, Esim4Params, WaterMatrix};
//! let genome = "ACGTACGTACGT";
//! let est = "ACGTACGT";
//! let params = Esim4Params{ matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 }, ..Default::default() };
//! let aln = esim4(genome, est, &params).unwrap();
//! assert!(aln.align.align_a.len() == aln.align.align_b.len());
//! ```
//!
use crate::common::{EmbossersError, WaterMatrix};
use crate::est2genome::{est2genome, Est2GenomeParams, Est2GenomeAlignment};

/// Which strand to align the EST on.
#[derive(Clone, Debug)]
pub enum StrandMode {
    /// Try both forward && reverse-complement && choose the higher score.
    Auto,
    /// Align EST as provided.
    Forward,
    /// Align reverse-complement of EST.
    Reverse,
}

/// Parameters for `esim4`.
#[derive(Clone, Debug)]
pub struct Esim4Params {
    /// Scoring matrix (DNA recommended).
    pub matrix: WaterMatrix,
    /// Gap open penalty (short indels).
    pub gap_open: f32,
    /// Gap extension penalty (short indels).
    pub gap_extend: f32,
    /// Integer scale to preserve fractional penalties.
    pub scale: f32,
    /// Minimum genomic gap length to classify as intron.
    pub intron_min: usize,
    /// Bonus for each canonical splice (`GT-AG` by default).
    pub splice_bonus: i32,
    /// Accept `GC-AG` as canonical.
    pub accept_gc_ag: bool,
    /// Accept `AT-AC` as canonical.
    pub accept_at_ac: bool,
    /// Whether to try reverse complement && pick the best.
    pub strand: StrandMode,
}

impl Default for Esim4Params {
    fn default() -> Self {
        Self {
            matrix: WaterMatrix::Dna{ match_score: 2, mismatch: -1 },
            gap_open: 10.0,
            gap_extend: 0.5,
            scale: 2.0,
            intron_min: 20,
            splice_bonus: 5,
            accept_gc_ag: true,
            accept_at_ac: false,
            strand: StrandMode::Auto,
        }
    }
}

/// Classification of an intron splice pair.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpliceClass { GtAg, GcAg, AtAc, NonCanonical }

/// A single exon mapped on the genomic sequence.
#[derive(Clone, Debug)]
pub struct Exon {
    /// Genomic coordinates **inclusive** `[start, end]` (0‑based) of the exon span.
    pub start: usize,
    pub end: usize,
}

/// Result of `esim4` alignment.
#[derive(Clone, Debug)]
pub struct Esim4Alignment {
    /// Strand used: `+` for forward, `-` for reverse-complement.
    pub strand: char,
    /// Underlying EST→genome alignment.
    pub align: Est2GenomeAlignment,
    /// Exon table derived from intron breaks.
    pub exons: Vec<Exon>,
    /// Splice class counts.
    pub splice_counts: (usize, usize, usize, usize), // GT-AG, GC-AG, AT-AC, NonCanonical
}

/// Align EST to genome with splice awareness && strand autodetection.
pub fn esim4(genome: &str, est: &str, params: &Esim4Params) -> Result<Esim4Alignment, EmbossersError> {
    if genome.is_empty() || est.is_empty() {
        return Err(EmbossersError::InvalidSequence("empty genome or est"));
    }
    // Prepare est2genome params
    let base = Est2GenomeParams{
        matrix: params.matrix.clone(),
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        scale: params.scale,
        intron_min: params.intron_min,
        splice_bonus: 0, // we apply bonuses after classifying (supports GC-AG/AT-AC)
    };

    let try_forward = || est2genome(genome, est, &base);
    let try_reverse = || {
        let rc = revcomp(est);
        est2genome(genome, &rc, &base)
    };

    let (strand, aln) = match params.strand {
        StrandMode::Forward => ('+', try_forward()?),
        StrandMode::Reverse => ('-', try_reverse()?),
        StrandMode::Auto => {
            let f = try_forward()?;
            let r = try_reverse()?;
            if r.score > f.score { ('-', r) } else { ('+', f) }
        }
    };

    // Reclassify introns && apply splice bonuses per class
    let genome_chars: Vec<char> = genome.chars().collect();
    let mut gt_ag = 0usize; let mut gc_ag = 0usize; let mut at_ac = 0usize; let mut non = 0usize;
    let mut bonus_total = 0i32;

    for intr in &aln.introns {
        let class = classify_splice(&genome_chars, intr.start, intr.end);
        match class {
            SpliceClass::GtAg => { gt_ag += 1; bonus_total += params.splice_bonus; }
            SpliceClass::GcAg => { gc_ag += 1; if params.accept_gc_ag { bonus_total += params.splice_bonus; } }
            SpliceClass::AtAc => { at_ac += 1; if params.accept_at_ac { bonus_total += params.splice_bonus; } }
            SpliceClass::NonCanonical => { non += 1; }
        }
    }
    // Apply bonus to score
    let mut align = aln;
    align.score += bonus_total;

    // Build exon table: contiguous aligned regions between introns
    let exons = build_exons_from_alignment(&align.align_a, &align.align_b, &genome_chars, &align.introns);

    Ok(Esim4Alignment{
        strand,
        align,
        exons,
        splice_counts: (gt_ag, gc_ag, at_ac, non),
    })
}

fn revcomp(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars().rev() {
        let c = match ch.to_ascii_uppercase() {
            'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'U' => 'A',
            'R' => 'Y', 'Y' => 'R', 'S' => 'S', 'W' => 'W',
            'K' => 'M', 'M' => 'K', 'B' => 'V', 'D' => 'H', 'H' => 'D', 'V' => 'B',
            'N' => 'N', _ => 'N',
        };
        out.push(c);
    }
    out
}

fn classify_splice(genome: &[char], start: usize, end: usize) -> SpliceClass {
    // donor: start..start+1 ; acceptor: end-1..end
    if start+1 < genome.len() && end >= 1 && end < genome.len() {
        let d0 = genome[start].to_ascii_uppercase();
        let d1 = genome[start+1].to_ascii_uppercase();
        let a0 = genome[end-1].to_ascii_uppercase();
        let a1 = genome[end].to_ascii_uppercase();
        if d0=='G' && d1=='T' && a0=='A' && a1=='G' { return SpliceClass::GtAg; }
        if d0=='G' && d1=='C' && a0=='A' && a1=='G' { return SpliceClass::GcAg; }
        if d0=='A' && d1=='T' && a0=='A' && a1=='C' { return SpliceClass::AtAc; }
    }
    SpliceClass::NonCanonical
}

use crate::est2genome::Intron;

fn build_exons_from_alignment(a_aln: &str, b_aln: &str, genome: &[char], introns: &[Intron]) -> Vec<Exon> {
    // Build mapping from aligned genome positions to genomic indices
    let a: Vec<char> = a_aln.chars().collect();
    let b: Vec<char> = b_aln.chars().collect();
    let mut g_idx: isize = -1;
    let mut exon_ranges: Vec<(usize,usize)> = Vec::new();
    let mut cur_start: Option<usize> = None;

    // We'll walk columns && start an exon when we see a match/mismatch (both not '-')
    // && end it when an intron boundary is crossed (long deletion already condensed into introns[]).
    let mut col_to_genome: Vec<Option<usize>> = Vec::with_capacity(a.len());
    for ch in &a {
        if *ch != '-' { g_idx += 1; col_to_genome.push(Some(g_idx as usize)); } else { col_to_genome.push(None); }
    }

    // Mark intron spans in genomic coordinates for quick lookup
    let mut in_intron = vec![false; genome.len()];
    for intr in introns {
        for i in intr.start..=intr.end {
            if i < in_intron.len() { in_intron[i] = true; }
        }
    }

    for (col, (&ca, &cb)) in a.iter().zip(b.iter()).enumerate() {
        let g = col_to_genome[col];
        let aligned = ca!='-' && cb!='-';
        if aligned {
            if let Some(gp) = g {
                if !in_intron[gp] {
                    if cur_start.is_none() { cur_start = Some(gp); }
                } else {
                    if let Some(s) = cur_start { exon_ranges.push((s, gp.saturating_sub(1))); cur_start=None; }
                }
            }
        } else {
            // gap breaks exon if currently tracking
            if let Some(s) = cur_start {
                let last_g = g.unwrap_or_else(|| g_idx as usize); // best-effort
                exon_ranges.push((s, last_g));
                cur_start=None;
            }
        }
    }
    if let Some(s) = cur_start {
        let end = g_idx as usize;
        exon_ranges.push((s, end));
    }

    // Merge small gaps && normalize
    let mut exons: Vec<Exon> = Vec::new();
    for (s,e) in exon_ranges {
        if e >= s { exons.push(Exon{ start: s, end: e }); }
    }
    exons
}
