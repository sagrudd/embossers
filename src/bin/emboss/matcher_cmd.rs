//! CLI subcommand for `emboss matcher` (Waterman–Eggert local alignments).
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use embossers::*;

/// Options for the `matcher` subcommand.
#[derive(Debug, Args)]
pub struct MatcherCmd {
    /// A-sequence FASTA file (first record used).
    #[arg(long, value_name="FILE")]
    pub asequence: PathBuf,
    /// B-sequence FASTA file (first record used).
    #[arg(long, value_name="FILE")]
    pub bsequence: PathBuf,
    /// Scoring matrix to use.
    #[arg(long, value_enum, default_value_t=MatrixChoice::Dna)]
    pub matrix: MatrixChoice,
    /// DNA match score (when --matrix dna).
    #[arg(long, default_value_t=2)]
    pub match_score: i32,
    /// DNA mismatch penalty (negative, when --matrix dna).
    #[arg(long, default_value_t=-1)]
    pub mismatch: i32,
    /// Gap open penalty.
    #[arg(long, default_value_t=10.0)]
    pub gapopen: f32,
    /// Gap extension penalty.
    #[arg(long, default_value_t=0.5)]
    pub gapextend: f32,
    /// Number of alternative local alignments to report.
    #[arg(long = "alternatives", alias = "alt", default_value_t=1)]
    pub alternatives: usize,
    /// Minimum score threshold; stop when best < threshold.
    #[arg(long, default_value_t=1)]
    pub minscore: i32,
    /// Output file for the report.
    #[arg(long, default_value="matcher.txt")]
    pub outfile: PathBuf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum MatrixChoice { Dna, Blosum62 }

pub fn run(cmd: MatcherCmd) -> Result<()> {
    let read_first = |path: &PathBuf| -> Result<FastaRecord> {
        let mut s = String::new();
        File::open(path).with_context(|| format!("open FASTA: {}", path.display()))?.read_to_string(&mut s)?;
        let recs = parse_fasta(&s);
        recs.into_iter().next().ok_or_else(|| anyhow::anyhow!("no FASTA records in {}", path.display()))
    };
    let a = read_first(&cmd.asequence)?;
    let b = read_first(&cmd.bsequence)?;

    let matrix = match cmd.matrix {
        MatrixChoice::Dna => WaterMatrix::Dna{ match_score: cmd.match_score, mismatch: cmd.mismatch },
        MatrixChoice::Blosum62 => WaterMatrix::Blosum62,
    };
    let params = MatcherParams{
        matrix,
        gap_open: cmd.gapopen,
        gap_extend: cmd.gapextend,
        scale: 2.0,
        alternatives: cmd.alternatives,
        min_score: cmd.minscore,
    };

    let hits = matcher(&a.seq, &b.seq, &params)?;

    // Write report
    let mut f = std::fs::File::create(&cmd.outfile).with_context(|| format!("create {}", cmd.outfile.display()))?;
    use std::io::Write;
    writeln!(f, "# EMBOSS-like MATCHER (Rust) report")?;
    writeln!(f, "# A: {} (len {})", a.id, a.seq.len())?;
    writeln!(f, "# B: {} (len {})", b.id, b.seq.len())?;
    writeln!(f, "Alternatives: {}", hits.len())?;
    for (idx, h) in hits.iter().enumerate() {
        writeln!(f, "
== Alignment {} ==", idx+1)?;
        writeln!(f, "Score: {}", h.score)?;
        writeln!(f, "A-range: {}..{}   B-range: {}..{}", h.a_start, h.a_end, h.b_start, h.b_end)?;
        writeln!(f, "Identity: {:.2}%   Gaps: {:.2}%", h.pct_identity, h.pct_gaps)?;
        writeln!(f, "CIGAR: {}", h.cigar)?;
        // blocks
        let a_chars: Vec<char> = h.align_a.chars().collect();
        let b_chars: Vec<char> = h.align_b.chars().collect();
        let mut i = 0usize;
        while i < a_chars.len() {
            let end = (i+60).min(a_chars.len());
            let a_block: String = a_chars[i..end].iter().collect();
            let b_block: String = b_chars[i..end].iter().collect();
            let mid: String = a_chars[i..end].iter().zip(b_chars[i..end].iter()).map(|(x,y)| {
                if *x=='-' || *y=='-' { ' ' } else if x.eq_ignore_ascii_case(y) { '|' } else { '.' }
            }).collect();
            writeln!(f, "A {}", a_block)?;
            writeln!(f, "  {}", mid)?;
            writeln!(f, "B {}", b_block)?;
            writeln!(f, "")?;
            i = end;
        }
    }
    Ok(())
}
