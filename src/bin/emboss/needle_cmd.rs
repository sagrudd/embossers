use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use embossers::*;

/// Options for the `needle` subcommand.
#[derive(Debug, Args)]
pub struct NeedleCmd {
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
    #[arg(long, default_value_t=1)]
    pub match_score: i32,
    /// DNA mismatch penalty (negative, when --matrix dna).
    #[arg(long, default_value_t=-1)]
    pub mismatch: i32,
    /// Gap open penalty (EMBOSS default ~10.0).
    #[arg(long, default_value_t=10.0)]
    pub gapopen: f32,
    /// Gap extension penalty (EMBOSS default ~0.5).
    #[arg(long, default_value_t=0.5)]
    pub gapextend: f32,
    /// Output file for a human-readable alignment.
    #[arg(long, default_value="needle.txt")]
    pub outfile: PathBuf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum MatrixChoice { Dna, Blosum62 }

pub fn run(cmd: NeedleCmd) -> Result<()> {
    // Read first record from each FASTA
    let read_first = |path: &PathBuf| -> Result<FastaRecord> {
        let mut s = String::new();
        File::open(path).with_context(|| format!("open FASTA: {}", path.display()))?.read_to_string(&mut s)?;
        let recs = parse_fasta(&s);
        recs.into_iter().next().ok_or_else(|| anyhow::anyhow!("no FASTA records in {}", path.display()))
    };
    let a = read_first(&cmd.asequence)?;
    let b = read_first(&cmd.bsequence)?;

    // Build params
    let matrix = match cmd.matrix {
        MatrixChoice::Dna => WaterMatrix::Dna{ match_score: cmd.match_score, mismatch: cmd.mismatch },
        MatrixChoice::Blosum62 => WaterMatrix::Blosum62,
    };
    let params = NeedleParams{ matrix, gap_open: cmd.gapopen, gap_extend: cmd.gapextend, ..Default::default() };

    let aln = needle(&a.seq, &b.seq, &params)?;

    // Pretty output
    let mut f = std::fs::File::create(&cmd.outfile).with_context(|| format!("create {}", cmd.outfile.display()))?;
    use std::io::Write;
    writeln!(f, "# EMBOSS-like NEEDLE (Rust) result")?;
    writeln!(f, "# A: {} (len {})", a.id, a.seq.len())?;
    writeln!(f, "# B: {} (len {})", b.id, b.seq.len())?;
    writeln!(f, "Score: {}", aln.score)?;
    writeln!(f, "Identity: {:.2}%   Gaps: {:.2}%", aln.pct_identity, aln.pct_gaps)?;
    writeln!(f, "CIGAR: {}", aln.cigar)?;
    writeln!(f, "")?;
    // Blocked alignment printing (60 cols)
    let a_chars: Vec<char> = aln.align_a.chars().collect();
    let b_chars: Vec<char> = aln.align_b.chars().collect();
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
    Ok(())
}
