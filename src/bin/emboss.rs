//! Command-line interface for the `embossers` crate.
//!
//! Subcommands are implemented in separate files (modules) under `src/bin/emboss/`:
//! - `complex_cmd.rs`
//! - `water_cmd.rs` (alias: `waters`)
//! - `needle_cmd.rs`
//!
use clap::{Parser, Subcommand};
use anyhow::Result;

#[derive(Debug, Parser)]
#[command(name="emboss", version=env!("CARGO_PKG_VERSION"), about="EMBOSS-inspired bioinformatics utilities (Rust)", disable_help_subcommand=true)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Find the linguistic complexity in nucleotide sequences (EMBOSS `complex`).
    Complex(complex_cmd::ComplexCmd),
    /// Smith–Waterman local alignment (EMBOSS `water`). Alias: `waters`.
    #[command(visible_alias = "waters")]
    Water(water_cmd::WaterCmd),
    /// Needleman–Wunsch global alignment (EMBOSS `needle`).
    Needle(needle_cmd::NeedleCmd),
    /// EST→genome spliced alignment (EMBOSS `est2genome`).
    Est2genome(est2genome_cmd::Est2GenomeCmd),
    /// Many-to-many global alignments (EMBOSS `needleall`).
    Needleall(needleall_cmd::NeedleAllCmd),
    /// Global alignment using Myers–Miller (linear space). EMBOSS `stretcher`.
    Stretcher(stretcher_cmd::StretcherCmd),
    /// SIM4-like spliced EST→genome alignment. EMBOSS `esim4`.
    Esim4(esim4_cmd::Esim4Cmd),
    /// Waterman–Eggert local alignments (EMBOSS `matcher`).
    Matcher(matcher_cmd::MatcherCmd),
    /// All-vs-all Waterman–Eggert local alignments (EMBOSS `seqmatchall`).
    Seqmatchall(seqmatchall_cmd::SeqMatchAllCmd),
    /// Seed-and-extend local alignments (EMBOSS `supermatcher`).
    Supermatcher(supermatcher_cmd::SuperMatcherCmd),
    /// Query vs set seed-and-extend (EMBOSS `wordfinder`).
    Wordfinder(wordfinder_cmd::WordFinderCmd),
    /// Exact word matches (EMBOSS `wordmatch`).
    Wordmatch(wordmatch_cmd::WordMatchCmd),
    /// Extend alignment with sequences (EMBOSS `seqalign`).
    Seqalign(seqalign_cmd::SeqAlignCmd),
    /// Progressive, consensus-guided MSA (EMBASSY domalign).
    Domainalign(domainalign_cmd::DomainAlignCmd),
    /// Within-sequence repeat finder (EMBASSY domainrep).
    Domainrep(domainrep_cmd::DomainRepCmd),
    /// Alignment statistics (EMBASSY oalistat).
    Oalistat(oalistat_cmd::OaListatCmd),
}

#[path = "emboss/complex_cmd.rs"] mod complex_cmd;
#[path = "emboss/water_cmd.rs"] mod water_cmd;
#[path = "emboss/needle_cmd.rs"] mod needle_cmd;
#[path = "emboss/est2genome_cmd.rs"] mod est2genome_cmd;
#[path = "emboss/needleall_cmd.rs"] mod needleall_cmd;
#[path = "emboss/stretcher_cmd.rs"] mod stretcher_cmd;
#[path = "emboss/esim4_cmd.rs"] mod esim4_cmd;
#[path = "emboss/matcher_cmd.rs"] mod matcher_cmd;
#[path = "emboss/seqmatchall_cmd.rs"] mod seqmatchall_cmd;
#[path = "emboss/seqalign_cmd.rs"] mod seqalign_cmd;
#[path = "emboss/oalistat_cmd.rs"] mod oalistat_cmd;
#[path = "emboss/domainrep_cmd.rs"] mod domainrep_cmd;
#[path = "emboss/domainalign_cmd.rs"] mod domainalign_cmd;
#[path = "emboss/wordmatch_cmd.rs"] mod wordmatch_cmd;
#[path = "emboss/wordfinder_cmd.rs"] mod wordfinder_cmd;
#[path = "emboss/supermatcher_cmd.rs"] mod supermatcher_cmd;

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Complex(cmd) => complex_cmd::run(cmd),
        Command::Water(cmd) => water_cmd::run(cmd),
        Command::Needle(cmd) => needle_cmd::run(cmd),
        Command::Est2genome(cmd) => est2genome_cmd::run(cmd),
        Command::Needleall(cmd) => needleall_cmd::run(cmd),
        Command::Stretcher(cmd) => stretcher_cmd::run(cmd),
        Command::Esim4(cmd) => esim4_cmd::run(cmd),
        Command::Matcher(cmd) => matcher_cmd::run(cmd),
        Command::Seqmatchall(cmd) => seqmatchall_cmd::run(cmd),
        Command::Supermatcher(cmd) => supermatcher_cmd::run(cmd),
        Command::Wordfinder(cmd) => wordfinder_cmd::run(cmd),
        Command::Wordmatch(cmd) => wordmatch_cmd::run(cmd),
        Command::Seqalign(cmd) => seqalign_cmd::run(cmd),
        Command::Domainalign(cmd) => domainalign_cmd::run(cmd),
        Command::Domainrep(cmd) => domainrep_cmd::run(cmd),
        Command::Oalistat(cmd) => oalistat_cmd::run(cmd),
    }
}
