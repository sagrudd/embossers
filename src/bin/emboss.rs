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
}

#[path = "emboss/complex_cmd.rs"] mod complex_cmd;
#[path = "emboss/water_cmd.rs"] mod water_cmd;
#[path = "emboss/needle_cmd.rs"] mod needle_cmd;

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Complex(cmd) => complex_cmd::run(cmd),
        Command::Water(cmd) => water_cmd::run(cmd),
        Command::Needle(cmd) => needle_cmd::run(cmd),
    }
}
