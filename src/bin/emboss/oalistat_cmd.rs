//! CLI for `emboss oalistat` (alignment statistics).
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::Args;
use embossers::*;

#[derive(Debug, Args)]
pub struct OaListatCmd {
    /// Alignment file (gapped FASTA or simple id<space>row format).
    #[arg(long, value_name="FILE")]
    pub align: PathBuf,
    /// Output report (TSV).
    #[arg(long, default_value="oalistat.tsv")]
    pub outfile: PathBuf,
}

pub fn run(cmd: OaListatCmd) -> Result<()> {
    let mut s = String::new();
    File::open(&cmd.align).with_context(|| format!("open alignment: {}", cmd.align.display()))?.read_to_string(&mut s)?;
    let st = oalistat(&s)?;

    // TSV
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').from_path(&cmd.outfile)?;
    w.write_record(["nseqs","cols","gap_columns_fraction","mean_pairwise_identity","conserved_columns"])?;
    w.write_record([
        st.nseqs.to_string(),
        st.cols.to_string(),
        format!("{:.4}", st.gap_columns_fraction),
        format!("{:.2}", st.mean_pairwise_identity),
        st.conserved_columns.to_string(),
    ])?;
    w.flush()?;
    Ok(())
}
