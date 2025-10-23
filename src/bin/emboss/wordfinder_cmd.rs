//! CLI for `emboss wordfinder` (query vs set using supermatcher).
//! Reads a query FASTA and a multi-FASTA set; writes per-target hits.
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::Args;
use embossers::*;

#[derive(Debug, Args)]
pub struct WordFinderCmd {
    /// Query FASTA (first record used).
    #[arg(long, value_name="FILE")]
    pub query: PathBuf,
    /// Target sequences (multi-FASTA).
    #[arg(long, value_name="FILE")]
    pub targets: PathBuf,
    /// Word length for seeding.
    #[arg(long, default_value_t=11)]
    pub wordlen: usize,
    /// Flank (half-window) around seed to extend.
    #[arg(long, default_value_t=200)]
    pub flank: usize,
    /// Minimum score to report.
    #[arg(long, default_value_t=1)]
    pub minscore: i32,
    /// Output file.
    #[arg(long, default_value="wordfinder.txt")]
    pub outfile: PathBuf,
}

pub fn run(cmd: WordFinderCmd) -> Result<()> {
    let mut s = String::new();
    File::open(&cmd.targets).with_context(|| format!("open FASTA: {}", cmd.targets.display()))?.read_to_string(&mut s)?;
    let targets = parse_fasta(&s);
    if targets.is_empty() { anyhow::bail!("no FASTA records in {}", cmd.targets.display()); }

    let mut q = String::new();
    File::open(&cmd.query).with_context(|| format!("open FASTA: {}", cmd.query.display()))?.read_to_string(&mut q)?;
    let qrec = parse_fasta(&q).into_iter().next().ok_or_else(|| anyhow::anyhow!("no FASTA records in {}", cmd.query.display()))?;

    let p = SuperMatcherParams{ wordlen: cmd.wordlen, flank: cmd.flank, min_score: cmd.minscore, ..Default::default() };
    let results = wordfinder(&qrec, &targets, &p)?;

    let mut f = std::fs::File::create(&cmd.outfile).with_context(|| format!("create {}", cmd.outfile.display()))?;
    use std::io::Write;
    writeln!(f, "# WORDFINDER query={}  targets={}", qrec.id, cmd.targets.display())?;
    for (tid, hits) in results {
        writeln!(f, "\n> Target {}", tid)?;
        writeln!(f, "Hits: {}", hits.len())?;
        for (i,h) in hits.iter().enumerate() {
            writeln!(f, "  [{}] score={}  A {}..{}  B {}..{}", i+1, h.score, h.a_start, h.a_end, h.b_start, h.b_end)?;
        }
    }
    Ok(())
}
