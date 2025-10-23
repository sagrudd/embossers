use std::fs::File;
use std::io::{self, Read, Write};
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::Args;
use embossers::*;

/// Options for the `complex` subcommand.
#[derive(Debug, Args)]
pub struct ComplexCmd {
    /// Input FASTA files. If omitted, reads FASTA from stdin.
    #[arg(long, value_name="FILES", num_args=1.., value_delimiter=' ')]
    pub sequence: Option<Vec<PathBuf>>,
    /// Sliding window length (L). Default: 100
    #[arg(long, default_value_t=100)]
    pub lwin: usize,
    /// Step size for sliding window. Default: 5
    #[arg(long, default_value_t=5)]
    pub step: usize,
    /// Minimum k-mer length (inclusive). Default: 4
    #[arg(long, default_value_t=4)]
    pub jmin: usize,
    /// Maximum k-mer length (inclusive). Default: 6
    #[arg(long, default_value_t=6)]
    pub jmax: usize,
    /// Calculate complexity aggregating windows across *all* sequences.
    #[arg(long, default_value_t=false)]
    pub omnia: bool,
    /// Number of random simulations to run (0 = none).
    #[arg(long, default_value_t=0)]
    pub sim: usize,
    /// If set, simulations use empirical base frequencies instead of uniform.
    #[arg(long, default_value_t=false)]
    pub freq: bool,
    /// If set, write a Uj table file. (Placeholder for future expansion)
    #[arg(long, default_value_t=false)]
    pub print: bool,
    /// Path to the Uj table file. Used when --print is set.
    #[arg(long, default_value="complex.ujtable")]
    pub ujtablefile: PathBuf,
    /// Path to the summary output file (TSV).
    #[arg(long, default_value="complex.tsv")]
    pub outfile: PathBuf,
    /// Optional path to echo the input sequences in FASTA format.
    #[arg(long)]
    pub outseq: Option<PathBuf>,
}

pub fn run(cmd: ComplexCmd) -> Result<()> {
    // Read sequences from files or stdin
    let fasta_text = if let Some(files) = &cmd.sequence {
        let mut buf = String::new();
        for p in files {
            let mut s = String::new();
            File::open(p).with_context(|| format!("open FASTA: {}", p.display()))?.read_to_string(&mut s)?;
            buf.push_str(&s);
            if !s.ends_with('\n') { buf.push('\n'); }
        }
        buf
    } else {
        let mut s = String::new();
        io::stdin().read_to_string(&mut s)?;
        s
    };

    // Parse FASTA (very tolerant)
    let records = parse_fasta(&fasta_text);
    if records.is_empty() {
        anyhow::bail!("no FASTA records found");
    }

    // Optionally echo sequences
    if let Some(path) = &cmd.outseq {
        let mut f = File::create(path).with_context(|| format!("create outseq {}", path.display()))?;
        for r in &records {
            writeln!(f, ">{}", r.id)?;
            writeln!(f, "{}", r.seq)?;
        }
    }

    let opts = ComplexOptions {
        lwin: cmd.lwin,
        step: cmd.step,
        jmin: cmd.jmin,
        jmax: cmd.jmax,
        sim: cmd.sim,
        freq_weighted_sim: cmd.freq,
    };

    // Compute
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(&cmd.outfile)
        .with_context(|| format!("create outfile {}", cmd.outfile.display()))?;
    // Header
    w.write_record(["id","windows","lwin","step","jmin","jmax","complexity"])?;

    // If omnia, compute aggregate across all sequences
    if cmd.omnia {
        for r in &records {
            let (c, _rows) = compute_complexity(&r.seq, &opts)?;
            w.write_record([&r.id, "-", &opts.lwin.to_string(), &opts.step.to_string(), &opts.jmin.to_string(), &opts.jmax.to_string(), &format!("{:.6}", c)])?;
        }
        let (c_all, _rows) = compute_complexity_multi(records.iter().map(|r| r.seq.as_str()), &opts)?;
        w.write_record(["ALL", "-", &opts.lwin.to_string(), &opts.step.to_string(), &opts.jmin.to_string(), &opts.jmax.to_string(), &format!("{:.6}", c_all)])?;
        w.flush()?;
        return Ok(());
    }

    // Not omnia: per-sequence lines
    for r in &records {
        let (c, _rows) = compute_complexity(&r.seq, &opts)?;
        w.write_record([&r.id, "-", &opts.lwin.to_string(), &opts.step.to_string(), &opts.jmin.to_string(), &opts.jmax.to_string(), &format!("{:.6}", c)])?;
    }
    w.flush()?;
    Ok(())
}
