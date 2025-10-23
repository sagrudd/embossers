//! CLI for `emboss merger` (DNA) and `emboss megamerger` (protein).
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum, Subcommand};
use embossers::*;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum MatrixChoice { Dna, Blosum62 }

#[derive(Debug, Args)]
pub struct MergeCommon {
    /// A-sequence FASTA (first record used).
    #[arg(long, value_name="FILE")]
    pub asequence: PathBuf,
    /// B-sequence FASTA (first record used).
    #[arg(long, value_name="FILE")]
    pub bsequence: PathBuf,
    /// Gap open.
    #[arg(long, default_value_t=10.0)]
    pub gapopen: f32,
    /// Gap extend.
    #[arg(long, default_value_t=0.5)]
    pub gapextend: f32,
    /// Minimum local score to accept merging.
    #[arg(long, default_value_t=1)]
    pub minscore: i32,
    /// Output FASTA file.
    #[arg(long, default_value="merged.fa")]
    pub outfile: PathBuf,
    /// FASTA ID for merged sequence.
    #[arg(long, default_value="merged")]
    pub id: String,
    /// Matrix to use (DNA: match/mismatch only; protein: BLOSUM62).
    #[arg(long, value_enum, default_value_t=MatrixChoice::Dna)]
    pub matrix: MatrixChoice,
    /// DNA match score (when --matrix dna).
    #[arg(long, default_value_t=2)]
    pub match_score: i32,
    /// DNA mismatch penalty (negative, when --matrix dna).
    #[arg(long, default_value_t=-1)]
    pub mismatch: i32,
}

fn read_first(path: &PathBuf) -> Result<FastaRecord> {
    let mut s = String::new();
    File::open(path).with_context(|| format!("open FASTA: {}", path.display()))?.read_to_string(&mut s)?;
    let recs = parse_fasta(&s);
    recs.into_iter().next().ok_or_else(|| anyhow::anyhow!("no FASTA records in {}", path.display()))
}

pub fn run_merger(cmd: MergeCommon, protein: bool) -> Result<()> {
    let a = read_first(&cmd.asequence)?;
    let b = read_first(&cmd.bsequence)?;
    let matrix = match cmd.matrix {
        MatrixChoice::Dna => WaterMatrix::Dna{ match_score: cmd.match_score, mismatch: cmd.mismatch },
        MatrixChoice::Blosum62 => WaterMatrix::Blosum62,
    };
    let p = MergeParams{ matrix, gap_open: cmd.gapopen, gap_extend: cmd.gapextend, scale: 2.0, min_score: cmd.minscore };
    let res = if protein { megamerger(&a.seq, &b.seq, &p)? } else { merger(&a.seq, &b.seq, &p)? };
    std::fs::write(&cmd.outfile, format!(">{}\n{}\n", cmd.id, res.merged)).with_context(|| format!("write {}", cmd.outfile.display()))?;
    Ok(())
}
