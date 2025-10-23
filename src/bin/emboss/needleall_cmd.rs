//! CLI subcommand implementation. Use via `emboss needleall`.
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use embossers::*;

/// Options for the `needleall` subcommand.
#[derive(Debug, Args)]
pub struct NeedleAllCmd {
    /// Multi-FASTA file of left/query sequences.
    #[arg(long, value_name="FILE")]
    pub aseqs: PathBuf,
    /// Multi-FASTA file of right/subject sequences.
    #[arg(long, value_name="FILE")]
    pub bseqs: PathBuf,
    /// Scoring matrix to use.
    #[arg(long, value_enum, default_value_t=MatrixChoice::Dna)]
    pub matrix: MatrixChoice,
    /// DNA match score (when --matrix dna).
    #[arg(long, default_value_t=1)]
    pub match_score: i32,
    /// DNA mismatch penalty (negative, when --matrix dna).
    #[arg(long, default_value_t=-1)]
    pub mismatch: i32,
    /// Gap open penalty.
    #[arg(long, default_value_t=10.0)]
    pub gapopen: f32,
    /// Gap extension penalty.
    #[arg(long, default_value_t=0.5)]
    pub gapextend: f32,
    /// Summary TSV file with pairwise scores and metrics.
    #[arg(long, default_value="needleall.tsv")]
    pub summary: PathBuf,
    /// Optional directory to write per-pair alignments (one file per pair).
    #[arg(long)]
    pub outdir: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum MatrixChoice { Dna, Blosum62 }

pub fn run(cmd: NeedleAllCmd) -> Result<()> {
    // Read multi-FASTA sets
    let read_all = |path: &PathBuf| -> Result<Vec<FastaRecord>> {
        let mut s = String::new();
        File::open(path).with_context(|| format!("open FASTA: {}", path.display()))?.read_to_string(&mut s)?;
        let recs = parse_fasta(&s);
        if recs.is_empty() { anyhow::bail!("no FASTA records in {}", path.display()); }
        Ok(recs)
    };
    let left = read_all(&cmd.aseqs)?;
    let right = read_all(&cmd.bseqs)?;

    // Params
    let matrix = match cmd.matrix {
        MatrixChoice::Dna => WaterMatrix::Dna{ match_score: cmd.match_score, mismatch: cmd.mismatch },
        MatrixChoice::Blosum62 => WaterMatrix::Blosum62,
    };
    let params = NeedleParams{ matrix, gap_open: cmd.gapopen, gap_extend: cmd.gapextend, ..Default::default() };

    // Compute
    let results = needleall_pairs(&left, &right, &params)?;

    // Summary TSV
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(&cmd.summary)
        .with_context(|| format!("create {}", cmd.summary.display()))?;
    w.write_record(["a_id","b_id","score","pct_identity","pct_gaps","cigar_len"])?;
    for r in &results {
        w.write_record([
            &r.a_id, &r.b_id,
            &r.score.to_string(),
            &format!("{:.2}", r.pct_identity),
            &format!("{:.2}", r.pct_gaps),
            &r.cigar.len().to_string()
        ])?;
    }
    w.flush()?;

    // Optional per-pair alignments
    if let Some(dir) = &cmd.outdir {
        std::fs::create_dir_all(dir).with_context(|| format!("mkdir {}", dir.display()))?;
        for r in &results {
            let fname = format!("{}_vs_{}.txt", sanitize(&r.a_id), sanitize(&r.b_id));
            let path = dir.join(fname);
            let mut f = std::fs::File::create(&path).with_context(|| format!("create {}", path.display()))?;
            use std::io::Write;
            writeln!(f, "# NEEDLEALL pair")?;
            writeln!(f, "A: {}", r.a_id)?;
            writeln!(f, "B: {}", r.b_id)?;
            writeln!(f, "Score: {}", r.score)?;
            writeln!(f, "Identity: {:.2}%   Gaps: {:.2}%", r.pct_identity, r.pct_gaps)?;
            writeln!(f, "CIGAR: {}", r.cigar)?;
            writeln!(f, "")?;
            // 60-col blocks
            let a: Vec<char> = r.align_a.chars().collect();
            let b: Vec<char> = r.align_b.chars().collect();
            let mut i = 0usize;
            while i < a.len() {
                let end = (i+60).min(a.len());
                let a_block: String = a[i..end].iter().collect();
                let b_block: String = b[i..end].iter().collect();
                let mid: String = a[i..end].iter().zip(b[i..end].iter()).map(|(x,y)| {
                    if *x=='-' || *y=='-' { ' ' } else if x.eq_ignore_ascii_case(y) { '|' } else { '.' }
                }).collect();
                writeln!(f, "A {}", a_block)?;
                writeln!(f, "  {}", mid)?;
                writeln!(f, "B {}", b_block)?;
                writeln!(f, "")?;
                i = end;
            }
        }
    }

    Ok(())
}

fn sanitize(s: &str) -> String {
    let mut t = String::with_capacity(s.len());
    for ch in s.chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '.' | '-' | '_') { t.push(ch); } else { t.push('_'); }
    }
    t
}
