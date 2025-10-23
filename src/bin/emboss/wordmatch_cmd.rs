//! CLI for `emboss wordmatch` (exact matches >= wordlen).
//! Writes an alignment-like report and optional GFF3 feature files.
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::Args;
use embossers::*;

#[derive(Debug, Args)]
pub struct WordMatchCmd {
    /// A-sequence FASTA (first record used).
    #[arg(long, value_name="FILE")]
    pub asequence: PathBuf,
    /// B-sequence FASTA (first record used).
    #[arg(long, value_name="FILE")]
    pub bsequence: PathBuf,
    /// Word size threshold.
    #[arg(long, default_value_t=10)]
    pub word: usize,
    /// Output report.
    #[arg(long, default_value="wordmatch.txt")]
    pub outfile: PathBuf,
    /// Optional GFF3 features for A.
    #[arg(long)]
    pub gff_a: Option<PathBuf>,
    /// Optional GFF3 features for B.
    #[arg(long)]
    pub gff_b: Option<PathBuf>,
}

pub fn run(cmd: WordMatchCmd) -> Result<()> {
    let read_first = |p: &PathBuf| -> Result<FastaRecord> {
        let mut s = String::new();
        File::open(p).with_context(|| format!("open FASTA: {}", p.display()))?.read_to_string(&mut s)?;
        let recs = parse_fasta(&s);
        recs.into_iter().next().ok_or_else(|| anyhow::anyhow!("no FASTA records in {}", p.display()))
    };
    let a = read_first(&cmd.asequence)?;
    let b = read_first(&cmd.bsequence)?;

    let params = WordMatchParams{ wordlen: cmd.word };
    let matches = wordmatch(&a.seq, &b.seq, &params)?;

    // Report
    let mut f = std::fs::File::create(&cmd.outfile).with_context(|| format!("create {}", cmd.outfile.display()))?;
    use std::io::Write;
    writeln!(f, "# WORDMATCH A={}  B={}  word>={}", a.id, b.id, cmd.word)?;
    writeln!(f, "Matches: {}", matches.len())?;
    for (i, m) in matches.iter().enumerate() {
        writeln!(f, "[{}] A {}..{}  B {}..{}  len={}  {}", i+1, m.a_start, m.a_end, m.b_start, m.b_end, m.len, m.seq)?;
    }

    // Optional GFF3
    if let Some(pa) = &cmd.gff_a {
        let mut fa = std::fs::File::create(pa).with_context(|| format!("create {}", pa.display()))?;
        writeln!(fa, "##gff-version 3")?;
        for (i, m) in matches.iter().enumerate() {
            writeln!(fa, "{}\twordmatch\texact_match\t{}\t{}\t.\t+\t.\tID=wm{}_A;Note=len {}", a.id, m.a_start+1, m.a_end+1, i+1, m.len)?;
        }
    }
    if let Some(pb) = &cmd.gff_b {
        let mut fb = std::fs::File::create(pb).with_context(|| format!("create {}", pb.display()))?;
        writeln!(fb, "##gff-version 3")?;
        for (i, m) in matches.iter().enumerate() {
            writeln!(fb, "{}\twordmatch\texact_match\t{}\t{}\t.\t+\t.\tID=wm{}_B;Note=len {}", b.id, m.b_start+1, m.b_end+1, i+1, m.len)?;
        }
    }

    Ok(())
}
