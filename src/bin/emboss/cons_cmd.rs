//! CLI for `emboss cons` (majority consensus).
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use embossers::*;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum ConsType { Dna, Protein }

#[derive(Debug, Args)]
pub struct ConsCmd {
    /// Alignment file (gapped FASTA or simple).
    #[arg(long, value_name="FILE")]
    pub align: PathBuf,
    /// Output FASTA.
    #[arg(long, default_value="consensus.fa")]
    pub outfile: PathBuf,
    /// Sequence type.
    #[arg(long="type", value_enum, default_value_t=ConsType::Dna)]
    pub seqtype: ConsType,
    /// Majority threshold (0.0..=1.0).
    #[arg(long, default_value_t=0.5)]
    pub threshold: f64,
    /// Consensus ID for FASTA header.
    #[arg(long, default_value="consensus")]
    pub id: String,
}

pub fn run(cmd: ConsCmd) -> Result<()> {
    let mut s = String::new();
    File::open(&cmd.align).with_context(|| format!("open alignment: {}", cmd.align.display()))?.read_to_string(&mut s)?;
    let is_dna = matches!(cmd.seqtype, ConsType::Dna);
    let cons = consensus_majority(&s, is_dna, cmd.threshold)?;
    std::fs::write(&cmd.outfile, format!(">{}\n{}\n", cmd.id, cons)).with_context(|| format!("write {}", cmd.outfile.display()))?;
    Ok(())
}
