//! CLI for `emboss consambig` (DNA IUPAC consensus).
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::Args;
use embossers::*;

#[derive(Debug, Args)]
pub struct ConsAmbigCmd {
    /// Alignment file (gapped FASTA or simple).
    #[arg(long, value_name="FILE")]
    pub align: PathBuf,
    /// Output FASTA.
    #[arg(long, default_value="consambig.fa")]
    pub outfile: PathBuf,
    /// Consensus ID for FASTA header.
    #[arg(long, default_value="consambig")]
    pub id: String,
}

pub fn run(cmd: ConsAmbigCmd) -> Result<()> {
    let mut s = String::new();
    File::open(&cmd.align).with_context(|| format!("open alignment: {}", cmd.align.display()))?.read_to_string(&mut s)?;
    let cons = consensus_ambig_dna(&s)?;
    std::fs::write(&cmd.outfile, format!(">{}\n{}\n", cmd.id, cons)).with_context(|| format!("write {}", cmd.outfile.display()))?;
    Ok(())
}
