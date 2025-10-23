//! CLI for `emboss seqalign` (extend alignment with sequences).
//!
//! This simplified stand-in reads a base alignment (EMBOSS "simple" or gapped FASTA)
//! and a directory of FASTA sequences, then appends each sequence to the alignment
//! by consensus-guided global alignment.
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use embossers::*;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum MatrixChoice { Dna, Blosum62 }

#[derive(Debug, Args)]
pub struct SeqAlignCmd {
    /// Base alignment file (EMBOSS 'simple' or gapped FASTA).
    #[arg(long, value_name="FILE")]
    pub daf: PathBuf,
    /// Directory containing FASTA sequences to add (DHF-like).
    #[arg(long, value_name="DIR")]
    pub dhfdir: PathBuf,
    /// Output alignment file.
    #[arg(long, default_value="seqalign.daf")]
    pub outdaf: PathBuf,
    /// Matrix to use for consensus-vs-sequence alignment.
    #[arg(long, value_enum, default_value_t=MatrixChoice::Dna)]
    pub matrix: MatrixChoice,
    /// DNA match score (when --matrix dna).
    #[arg(long, default_value_t=2)]
    pub match_score: i32,
    /// DNA mismatch penalty (negative, when --matrix dna).
    #[arg(long, default_value_t=-1)]
    pub mismatch: i32,
}

pub fn run(cmd: SeqAlignCmd) -> Result<()> {
    // Read base alignment
    let mut s = String::new();
    File::open(&cmd.daf).with_context(|| format!("open alignment: {}", cmd.daf.display()))?.read_to_string(&mut s)?;

    // Iterate FASTA files in dhfdir (non-recursive)
    let mut out = s.clone();
    for entry in std::fs::read_dir(&cmd.dhfdir)? {
        let e = entry?;
        if !e.file_type()?.is_file() { continue; }
        let path = e.path();
        if let Some(ext) = path.extension() {
            let ext = ext.to_string_lossy().to_ascii_lowercase();
            if ext!="fa" && ext!="fasta" && ext!="fna" { continue; }
        } else { continue; }
        let mut fs = String::new();
        File::open(&path).with_context(|| format!("open FASTA: {}", path.display()))?.read_to_string(&mut fs)?;
        let recs = parse_fasta(&fs);
        if let Some(rec) = recs.into_iter().next() {
            let matrix = match cmd.matrix {
                MatrixChoice::Dna => WaterMatrix::Dna{ match_score: cmd.match_score, mismatch: cmd.mismatch },
                MatrixChoice::Blosum62 => WaterMatrix::Blosum62,
            };
            let extended = seqalign_extend(&out, &rec.id, &rec.seq, &matrix)?;
            out = extended;
        }
    }
    std::fs::write(&cmd.outdaf, out).with_context(|| format!("write {}", cmd.outdaf.display()))?;
    Ok(())
}
