//! CLI for `emboss supermatcher`.
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::Args;
use embossers::*;

/// Options for `supermatcher`.
#[derive(Debug, Args)]
pub struct SuperMatcherCmd {
    /// A-sequence FASTA (first record used).
    #[arg(long, value_name="FILE")]
    pub asequence: PathBuf,
    /// B-sequence FASTA (first record used).
    #[arg(long, value_name="FILE")]
    pub bsequence: PathBuf,
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
    #[arg(long, default_value="supermatcher.txt")]
    pub outfile: PathBuf,
}

pub fn run(cmd: SuperMatcherCmd) -> Result<()> {
    let read_first = |p: &PathBuf| -> Result<FastaRecord> {
        let mut s = String::new();
        File::open(p).with_context(|| format!("open FASTA: {}", p.display()))?.read_to_string(&mut s)?;
        let recs = parse_fasta(&s);
        recs.into_iter().next().ok_or_else(|| anyhow::anyhow!("no FASTA records in {}", p.display()))
    };
    let a = read_first(&cmd.asequence)?;
    let b = read_first(&cmd.bsequence)?;

    let params = SuperMatcherParams{ wordlen: cmd.wordlen, flank: cmd.flank, min_score: cmd.minscore, ..Default::default() };
    let hits = supermatcher(&a.seq, &b.seq, &params)?;

    let mut f = std::fs::File::create(&cmd.outfile).with_context(|| format!("create {}", cmd.outfile.display()))?;
    use std::io::Write;
    writeln!(f, "# SUPERMATCHER results  A={}  B={}", a.id, b.id)?;
    writeln!(f, "Hits: {}", hits.len())?;
    for (i,h) in hits.iter().enumerate() {
        writeln!(f, "\n== Hit {} ==", i+1)?;
        writeln!(f, "Score: {}", h.score)?;
        writeln!(f, "A {}..{}   B {}..{}", h.a_start, h.a_end, h.b_start, h.b_end)?;
        writeln!(f, "Identity: {:.2}%   Gaps: {:.2}%", h.pct_identity, h.pct_gaps)?;
        writeln!(f, "CIGAR: {}", h.cigar)?;
        // 60-col block
        let aa: Vec<char> = h.align_a.chars().collect();
        let bb: Vec<char> = h.align_b.chars().collect();
        let mut k=0usize;
        while k<aa.len() {
            let end=(k+60).min(aa.len());
            let a_block: String = aa[k..end].iter().collect();
            let b_block: String = bb[k..end].iter().collect();
            let mid: String = aa[k..end].iter().zip(bb[k..end].iter()).map(|(x,y)| {
                if *x=='-' || *y=='-' { ' ' } else if x.eq_ignore_ascii_case(y) { '|' } else { '.' }
            }).collect();
            writeln!(f, "A {}", a_block)?;
            writeln!(f, "  {}", mid)?;
            writeln!(f, "B {}", b_block)?;
            writeln!(f, "")?;
            k=end;
        }
    }
    Ok(())
}
