//! CLI for `emboss extractalign`.
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::Args;
use embossers::*;

#[derive(Debug, Args)]
pub struct ExtractAlignCmd {
    /// Alignment file (gapped FASTA or simple id<space>row).
    #[arg(long, value_name="FILE")]
    pub align: PathBuf,
    /// 1-based inclusive ranges: e.g. "5-80,120-150".
    #[arg(long, value_name="RANGES")]
    pub range: String,
    /// Output alignment (simple format).
    #[arg(long, default_value="extractalign.daf")]
    pub outfile: PathBuf,
}

pub fn run(cmd: ExtractAlignCmd) -> Result<()> {
    let mut s = String::new();
    File::open(&cmd.align).with_context(|| format!("open alignment: {}", cmd.align.display()))?.read_to_string(&mut s)?;
    let ranges = parse_ranges(&cmd.range)?;
    let out = extractalign_ranges(&s, &ranges)?;
    std::fs::write(&cmd.outfile, out).with_context(|| format!("write {}", cmd.outfile.display()))?;
    Ok(())
}

fn parse_ranges(s: &str) -> Result<Vec<(usize,usize)>> {
    let mut v = Vec::new();
    for part in s.split(',') {
        let t = part.trim();
        if t.is_empty() { continue; }
        let (a,b) = if let Some((l,r)) = t.split_once('-') {
            (l.trim().parse::<usize>()?, r.trim().parse::<usize>()?)
        } else {
            let x = t.parse::<usize>()?; (x,x)
        };
        v.push((a,b));
    }
    if v.is_empty() { anyhow::bail!("no valid ranges: {}", s); }
    Ok(v)
}
