//! CLI for `emboss domainrep` (within-sequence domain repeat finder).
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::Args;
use embossers::*;

#[derive(Debug, Args)]
pub struct DomainRepCmd {
    /// Input sequence FASTA (first record used).
    #[arg(long, value_name="FILE")]
    pub sequence: PathBuf,
    /// Seed size (k-mer length).
    #[arg(long, default_value_t=3)]
    pub wordlen: usize,
    /// Minimum repeat length.
    #[arg(long, default_value_t=12)]
    pub minlen: usize,
    /// Output report.
    #[arg(long, default_value="domainrep.tsv")]
    pub outfile: PathBuf,
    /// Optional GFF3 output.
    #[arg(long)]
    pub gff: Option<PathBuf>,
}

pub fn run(cmd: DomainRepCmd) -> Result<()> {
    let mut s = String::new();
    File::open(&cmd.sequence).with_context(|| format!("open FASTA: {}", cmd.sequence.display()))?.read_to_string(&mut s)?;
    let rec = parse_fasta(&s).into_iter().next().ok_or_else(|| anyhow::anyhow!("no FASTA records in {}", cmd.sequence.display()))?;
    let p = DomainRepParams{ wordlen: cmd.wordlen, min_len: cmd.minlen };
    let reps = domainrep(&rec.seq, &p)?;

    // TSV
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(&cmd.outfile)?;
    w.write_record(["a_start","a_end","b_start","b_end","len"])?;
    for r in &reps {
        w.write_record([r.a_start.to_string(), r.a_end.to_string(), r.b_start.to_string(), r.b_end.to_string(), r.len.to_string()])?;
    }
    w.flush()?;

    if let Some(gff) = &cmd.gff {
        use std::io::Write;
        let mut f = std::fs::File::create(gff).with_context(|| format!("create {}", gff.display()))?;
        writeln!(f, "##gff-version 3")?;
        for (i, r) in reps.iter().enumerate() {
            writeln!(f, "{}\tdomainrep\trepeat\t{}\t{}\t.\t+\t.\tID=domrep{};Note=len {}", rec.id, r.a_start+1, r.a_end+1, i+1, r.len)?;
            writeln!(f, "{}\tdomainrep\trepeat\t{}\t{}\t.\t+\t.\tID=domrep{}b;Note=len {}", rec.id, r.b_start+1, r.b_end+1, i+1, r.len)?;
        }
    }
    Ok(())
}
