//! CLI for `emboss domainalign` (progressive, consensus-guided MSA).
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use embossers::*;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum MatrixChoice { Dna, Blosum62 }

#[derive(Debug, Args)]
pub struct DomainAlignCmd {
    /// Input multi-FASTA with domain sequences.
    #[arg(long, value_name="FILE")]
    pub seqs: PathBuf,
    /// Output alignment (simple format: id<space>row).
    #[arg(long, default_value="domainalign.daf")]
    pub outdaf: PathBuf,
    /// Matrix to use for pairwise alignment.
    #[arg(long, value_enum, default_value_t=MatrixChoice::Dna)]
    pub matrix: MatrixChoice,
    /// DNA match score (when --matrix dna).
    #[arg(long, default_value_t=2)]
    pub match_score: i32,
    /// DNA mismatch penalty (negative, when --matrix dna).
    #[arg(long, default_value_t=-1)]
    pub mismatch: i32,
    /// Gap open.
    #[arg(long, default_value_t=10.0)]
    pub gapopen: f32,
    /// Gap extend.
    #[arg(long, default_value_t=0.5)]
    pub gapextend: f32,
}

pub fn run(cmd: DomainAlignCmd) -> Result<()> {
    let mut s = String::new();
    File::open(&cmd.seqs).with_context(|| format!("open FASTA: {}", cmd.seqs.display()))?.read_to_string(&mut s)?;
    let recs = parse_fasta(&s);
    if recs.is_empty() { anyhow::bail!("no FASTA records in {}", cmd.seqs.display()); }
    let matrix = match cmd.matrix {
        MatrixChoice::Dna => WaterMatrix::Dna{ match_score: cmd.match_score, mismatch: cmd.mismatch },
        MatrixChoice::Blosum62 => WaterMatrix::Blosum62,
    };
    let p = DomainAlignParams{ matrix, gap_open: cmd.gapopen, gap_extend: cmd.gapextend, scale: 2.0 };
    let daf = domainalign(&recs, &p)?;
    std::fs::write(&cmd.outdaf, daf).with_context(|| format!("write {}", cmd.outdaf.display()))?;
    Ok(())
}
