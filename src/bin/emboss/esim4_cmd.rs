//! CLI subcommand implementation. Use via `emboss esim4`.
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use embossers::*;

/// Options for the `esim4` subcommand.
#[derive(Debug, Args)]
pub struct Esim4Cmd {
    /// Genomic DNA FASTA file (first record used).
    #[arg(long, value_name="FILE")]
    pub genome: PathBuf,
    /// EST/cDNA FASTA file (first record used).
    #[arg(long, value_name="FILE")]
    pub est: PathBuf,
    /// Scoring matrix to use.
    #[arg(long, value_enum, default_value_t=MatrixChoice::Dna)]
    pub matrix: MatrixChoice,
    /// DNA match score (when --matrix dna).
    #[arg(long, default_value_t=2)]
    pub match_score: i32,
    /// DNA mismatch penalty (negative, when --matrix dna).
    #[arg(long, default_value_t=-1)]
    pub mismatch: i32,
    /// Gap open penalty (short indels).
    #[arg(long, default_value_t=10.0)]
    pub gapopen: f32,
    /// Gap extension penalty (short indels).
    #[arg(long, default_value_t=0.5)]
    pub gapextend: f32,
    /// Minimum genomic gap length to classify as intron.
    #[arg(long, default_value_t=20)]
    pub intron_min: usize,
    /// Bonus applied to the score for each accepted canonical intron.
    #[arg(long, default_value_t=5)]
    pub splice_bonus: i32,
    /// Accept GC-AG as canonical.
    #[arg(long, default_value_t=true)]
    pub accept_gc_ag: bool,
    /// Accept AT-AC as canonical.
    #[arg(long, default_value_t=false)]
    pub accept_at_ac: bool,
    /// Strand handling.
    #[arg(long, value_enum, default_value_t=StrandChoice::Auto)]
    pub strand: StrandChoice,
    /// Output file for a human-readable alignment report.
    #[arg(long, default_value="esim4.txt")]
    pub outfile: PathBuf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum MatrixChoice { Dna, Blosum62 }

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum StrandChoice { Auto, Forward, Reverse }

pub fn run(cmd: Esim4Cmd) -> Result<()> {
    let read_first = |path: &PathBuf| -> Result<FastaRecord> {
        let mut s = String::new();
        File::open(path).with_context(|| format!("open FASTA: {}", path.display()))?.read_to_string(&mut s)?;
        let recs = parse_fasta(&s);
        recs.into_iter().next().ok_or_else(|| anyhow::anyhow!("no FASTA records in {}", path.display()))
    };
    let g = read_first(&cmd.genome)?;
    let e = read_first(&cmd.est)?;

    let matrix = match cmd.matrix {
        MatrixChoice::Dna => WaterMatrix::Dna{ match_score: cmd.match_score, mismatch: cmd.mismatch },
        MatrixChoice::Blosum62 => WaterMatrix::Blosum62,
    };
    let strand = match cmd.strand {
        StrandChoice::Auto => StrandMode::Auto,
        StrandChoice::Forward => StrandMode::Forward,
        StrandChoice::Reverse => StrandMode::Reverse,
    };
    let params = Esim4Params{
        matrix,
        gap_open: cmd.gapopen,
        gap_extend: cmd.gapextend,
        scale: 2.0,
        intron_min: cmd.intron_min,
        splice_bonus: cmd.splice_bonus,
        accept_gc_ag: cmd.accept_gc_ag,
        accept_at_ac: cmd.accept_at_ac,
        strand,
    };

    let aln = esim4(&g.seq, &e.seq, &params)?;

    // Write pretty output
    let mut f = std::fs::File::create(&cmd.outfile).with_context(|| format!("create {}", cmd.outfile.display()))?;
    use std::io::Write;
    writeln!(f, "# EMBOSS-like ESIM4 (Rust) result")?;
    writeln!(f, "# Genome: {} (len {})", g.id, g.seq.len())?;
    writeln!(f, "# EST:    {} (len {})", e.id, e.seq.len())?;
    writeln!(f, "Strand: {}", aln.strand)?;
    writeln!(f, "Score: {}", aln.align.score)?;
    writeln!(f, "Identity: {:.2}%   Gaps: {:.2}%", aln.align.pct_identity, aln.align.pct_gaps)?;
    writeln!(f, "CIGAR: {}", aln.align.cigar)?;
    let (gt,gc,at,non) = aln.splice_counts;
    writeln!(f, "Splice counts: GT-AG={} GC-AG={} AT-AC={} NonCanonical={}", gt,gc,at,non)?;
    if aln.exons.is_empty() { writeln!(f, "Exons: 0")?; } else {
        writeln!(f, "Exons: {}", aln.exons.len())?;
        for (i,x) in aln.exons.iter().enumerate() {
            writeln!(f, "  {}: {}..{}", i+1, x.start, x.end)?;
        }
    }
    writeln!(f, "")?;
    // Alignment blocks (60-col)
    let a: Vec<char> = aln.align.align_a.chars().collect();
    let b: Vec<char> = aln.align.align_b.chars().collect();
    let mut i = 0usize;
    while i < a.len() {
        let end = (i+60).min(a.len());
        let a_block: String = a[i..end].iter().collect();
        let b_block: String = b[i..end].iter().collect();
        let mid: String = a[i..end].iter().zip(b[i..end].iter()).map(|(x,y)| {
            if *x=='-' || *y=='-' { ' ' } else if x.eq_ignore_ascii_case(y) { '|' } else { '.' }
        }).collect();
        writeln!(f, "G {}", a_block)?;
        writeln!(f, "  {}", mid)?;
        writeln!(f, "E {}", b_block)?;
        writeln!(f, "")?;
        i = end;
    }
    Ok(())
}
