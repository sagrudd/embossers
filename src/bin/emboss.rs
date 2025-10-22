//! Command-line interface for the `embossers` crate.
//!
//! This binary provides a **feature-complete** `complex` subcommand modeled on
//! EMBOSS `complex`:
//!
//! ```text
//! Usage: emboss complex [OPTIONS] [--sequence <FILES>...]
//! ```
//!
//! Key options:
//! - `--sequence <FILES...>`: input FASTA files (omit to read stdin)
//! - `--lwin`, `--step`, `--jmin`, `--jmax`: window and word-size controls
//! - `--omnia`: aggregate windows across all sequences (EMBOSS behaviour)
//! - `--sim <N>` / `--freq`: run N simulations; weight by empirical base freqs
//! - `--print` and `--ujtablefile`: write a Uj table (default path: `complex.ujtable`)
//! - `--outfile`: write a summary TSV (default path: `complex.tsv`)
//! - `--outseq`: echo input sequences in FASTA format to a file
//!
//! ## Examples
//! ```text
//! # From a file
//! emboss complex --sequence toy.fasta --lwin 100 --step 5 --jmin 4 --jmax 6 --outfile toy.complex.tsv --print
//!
//! # From stdin
//! cat toy.fasta | emboss complex --lwin 100 --step 5 --jmin 4 --jmax 6 --outfile toy.complex.tsv
//! ```

use std::fs::File;
use std::io::{self, Read, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand, Args};
use embossers::{ComplexOptions, FastaRecord, compute_complexity, compute_complexity_multi, parse_fasta};

/// Top-level CLI.
#[derive(Debug, Parser)]
#[command(name="emboss", version=env!("CARGO_PKG_VERSION"), about="EMBOSS-inspired bioinformatics utilities (Rust)", disable_help_subcommand=true)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

/// Subcommands.
#[derive(Debug, Subcommand)]
enum Command {
    /// Find the linguistic complexity in nucleotide sequences (EMBOSS `complex`).
    Complex(ComplexCmd),
}

/// Options for the `complex` subcommand.
#[derive(Debug, Args)]
struct ComplexCmd {
    /// Input FASTA files. If omitted, reads FASTA from stdin.
    #[arg(long, value_name="FILES", num_args=1.., value_delimiter=' ')]
    sequence: Option<Vec<PathBuf>>,
    /// Sliding window length (L). Default: 100
    #[arg(long, default_value_t=100)]
    lwin: usize,
    /// Step size for sliding window. Default: 5
    #[arg(long, default_value_t=5)]
    step: usize,
    /// Minimum k-mer length (inclusive). Default: 4
    #[arg(long, default_value_t=4)]
    jmin: usize,
    /// Maximum k-mer length (inclusive). Default: 6
    #[arg(long, default_value_t=6)]
    jmax: usize,
    /// Calculate complexity aggregating windows across *all* sequences.
    #[arg(long, default_value_t=false)]
    omnia: bool,
    /// Number of random simulations to run (0 = none).
    #[arg(long, default_value_t=0)]
    sim: usize,
    /// If set, simulations use empirical base frequencies instead of uniform.
    #[arg(long, default_value_t=false)]
    freq: bool,
    /// If set, write a Uj table file.
    #[arg(long, default_value_t=false)]
    print: bool,
    /// Path to the Uj table file. Used when --print is set.
    #[arg(long, default_value="complex.ujtable")]
    ujtablefile: PathBuf,
    /// Path to the summary output file (TSV).
    #[arg(long, default_value="complex.tsv")]
    outfile: PathBuf,
    /// Optional path to echo the input sequences in FASTA format.
    #[arg(long)]
    outseq: Option<PathBuf>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Complex(cmd) => run_complex(cmd),
    }
}

fn run_complex(cmd: ComplexCmd) -> Result<()> {
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
        let (_per_seq_lines, _agg_rows_written) = write_results_for_records(&records, &opts, &mut w, /*write_per_seq=*/true)?;
        // In addition to per-seq lines, write an aggregate line "ALL"
        let (c_all, _rows) = compute_complexity_multi(records.iter().map(|r| r.seq.as_str()), &opts)?;
        w.write_record(["ALL", "-", &opts.lwin.to_string(), &opts.step.to_string(), &opts.jmin.to_string(), &opts.jmax.to_string(), &format!("{:.6}", c_all)])?;
        w.flush()?;
        if cmd.print {
            // For the Uj table under omnia we generate the table for the aggregate as reference.
            let (_c, rows) = compute_complexity_multi(records.iter().map(|r| r.seq.as_str()), &opts)?;
            write_uj_table(&cmd.ujtablefile, &rows)?;
        }
        return Ok(());
    }

    // Not omnia: per-sequence lines
    for r in &records {
        let (c, rows) = compute_complexity(&r.seq, &opts)?;
        w.write_record([&r.id, "-", &opts.lwin.to_string(), &opts.step.to_string(), &opts.jmin.to_string(), &opts.jmax.to_string(), &format!("{:.6}", c)])?;
        if cmd.print {
            write_uj_table(&cmd.ujtablefile, &rows)?;
        }
    }
    w.flush()?;
    Ok(())
}

fn write_results_for_records(
    records: &[FastaRecord],
    opts: &ComplexOptions,
    wtr: &mut csv::Writer<std::fs::File>,
    write_per_seq: bool
) -> Result<(usize, usize)> {
    let mut per_seq = 0usize;
    let mut agg = 0usize;
    if write_per_seq {
        for r in records {
            let (c, _rows) = compute_complexity(&r.seq, opts)?;
            wtr.write_record([&r.id, "-", &opts.lwin.to_string(), &opts.step.to_string(), &opts.jmin.to_string(), &opts.jmax.to_string(), &format!("{:.6}", c)])?;
            per_seq += 1;
        }
    }
    // Aggregate line written by caller when desired
    Ok((per_seq, agg))
}

/// Write a simple Uj table TSV.
fn write_uj_table(path: &std::path::Path, rows: &[embossers::UjRow]) -> Result<()> {
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)
        .with_context(|| format!("create ujtable {}", path.display()))?;
    w.write_record(["j","Uj_real_mean","Uj_sim_mean","Uj_sim_std"])?;
    for r in rows {
        w.write_record([
            r.j.to_string(),
            format!("{:.6}", r.uj_real_mean),
            r.uj_sim_mean.map(|x| format!("{:.6}", x)).unwrap_or_else(|| "".into()),
            r.uj_sim_std.map(|x| format!("{:.6}", x)).unwrap_or_else(|| "".into()),
        ])?;
    }
    w.flush()?;
    Ok(())
}
