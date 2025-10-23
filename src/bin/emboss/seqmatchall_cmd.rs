//! CLI subcommand for `emboss seqmatchall` (all-vs-all Waterman–Eggert local alignments).
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use embossers::*;

/// Options for the `seqmatchall` subcommand.
#[derive(Debug, Args)]
pub struct SeqMatchAllCmd {
    /// Multi-FASTA file containing all sequences.
    #[arg(long, value_name="FILE")]
    pub seqs: PathBuf,
    /// Scoring matrix to use.
    #[arg(long, value_enum, default_value_t=MatrixChoice::Dna)]
    pub matrix: MatrixChoice,
    /// DNA match score (when --matrix dna).
    #[arg(long, default_value_t=2)]
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
    /// Number of alternative local alignments per pair.
    #[arg(long = "alternatives", alias = "alt", default_value_t=1)]
    pub alternatives: usize,
    /// Minimum score threshold for local alignments.
    #[arg(long, default_value_t=1)]
    pub minscore: i32,
    /// Include self-vs-self comparisons (i==j).
    #[arg(long)]
    pub include_self: bool,
    /// Summary TSV file.
    #[arg(long, default_value="seqmatchall.tsv")]
    pub summary: PathBuf,
    /// Optional directory to write per-pair top alignment files.
    #[arg(long)]
    pub outdir: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum MatrixChoice { Dna, Blosum62 }

pub fn run(cmd: SeqMatchAllCmd) -> Result<()> {
    // Read multi-FASTA
    let mut s = String::new();
    File::open(&cmd.seqs).with_context(|| format!("open FASTA: {}", cmd.seqs.display()))?.read_to_string(&mut s)?;
    let recs = parse_fasta(&s);
    if recs.is_empty() { anyhow::bail!("no FASTA records in {}", cmd.seqs.display()); }

    // Build matcher params
    let matrix = match cmd.matrix {
        MatrixChoice::Dna => WaterMatrix::Dna{ match_score: cmd.match_score, mismatch: cmd.mismatch },
        MatrixChoice::Blosum62 => WaterMatrix::Blosum62,
    };
    let mparams = MatcherParams{
        matrix,
        gap_open: cmd.gapopen,
        gap_extend: cmd.gapextend,
        scale: 2.0,
        alternatives: cmd.alternatives,
        min_score: cmd.minscore,
    };
    let params = SeqMatchAllParams{ matcher: mparams.clone(), include_self: cmd.include_self };

    // Compute
    let results = seqmatchall(&recs, &params)?;

    // Summary TSV
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(&cmd.summary)
        .with_context(|| format!("create {}", cmd.summary.display()))?;
    w.write_record(["a_id","b_id","score","pct_identity","pct_gaps","n_hits","top_cigar_len"])?;
    for r in &results {
        w.write_record([
            &r.a_id, &r.b_id,
            &r.score.to_string(),
            &format!("{:.2}", r.pct_identity),
            &format!("{:.2}", r.pct_gaps),
            &r.n_hits.to_string(),
            &r.top_cigar_len.to_string(),
        ])?;
    }
    w.flush()?;

    // Optional per-pair top alignment files
    if let Some(dir) = &cmd.outdir {
        std::fs::create_dir_all(dir).with_context(|| format!("mkdir {}", dir.display()))?;
        for i in 0..recs.len() {
            let j_start = if cmd.include_self { i } else { i+1 };
            if j_start >= recs.len() { continue; }
            for j in j_start..recs.len() {
                let hits = matcher(&recs[i].seq, &recs[j].seq, &mparams)?;
                if hits.is_empty() { continue; }
                let h = &hits[0];
                let fname = format!("{}_vs_{}.txt", sanitize(&recs[i].id), sanitize(&recs[j].id));
                let path = dir.join(fname);
                let mut f = std::fs::File::create(&path).with_context(|| format!("create {}", path.display()))?;
                use std::io::Write;
                writeln!(f, "# SEQMATCHALL pair")?;
                writeln!(f, "A: {}", recs[i].id)?;
                writeln!(f, "B: {}", recs[j].id)?;
                writeln!(f, "Score: {}", h.score)?;
                writeln!(f, "Identity: {:.2}%   Gaps: {:.2}%", h.pct_identity, h.pct_gaps)?;
                writeln!(f, "CIGAR: {}", h.cigar)?;
                writeln!(f, "")?;
                // 60-col alignment block
                let a: Vec<char> = h.align_a.chars().collect();
                let b: Vec<char> = h.align_b.chars().collect();
                let mut k = 0usize;
                while k < a.len() {
                    let end = (k+60).min(a.len());
                    let a_block: String = a[k..end].iter().collect();
                    let b_block: String = b[k..end].iter().collect();
                    let mid: String = a[k..end].iter().zip(b[k..end].iter()).map(|(x,y)| {
                        if *x=='-' || *y=='-' { ' ' } else if x.eq_ignore_ascii_case(y) { '|' } else { '.' }
                    }).collect();
                    writeln!(f, "A {}", a_block)?;
                    writeln!(f, "  {}", mid)?;
                    writeln!(f, "B {}", b_block)?;
                    writeln!(f, "")?;
                    k = end;
                }
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
