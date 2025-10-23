//! # embossers
//!
//! A small, batteries‑included crate that reimplements several **EMBOSS** utilities
//! in pure Rust with a simple CLI. It’s designed for teaching, scripting, and
//! reproducible pipelines where a lightweight dependency helps.
//!
//! ## Modules
//! - [`complex`] — linguistic sequence complexity (EMBOSS **`complex`**).
//! - [`water`]   — Smith–Waterman local alignment with affine gaps (EMBOSS **`water`**).
//! - [`needle`]  — Needleman–Wunsch global alignment with affine gaps (EMBOSS **`needle`**).
//!
//! ## CLI
//! The `emboss` binary provides subcommands that mirror EMBOSS. See
//! `emboss --help` and subcommand help, or the README for full examples.
//!
//! ## Example (library)
//! ```rust
//! use embossers::{ComplexOptions, compute_complexity, water, WaterParams, WaterMatrix};
//! let (c, _rows) = compute_complexity("ACGTACGT", &ComplexOptions::default()).unwrap();
//! let w = water("PAWHEAE", "HEAGAWGHEE", &WaterParams{ matrix: WaterMatrix::Blosum62, ..Default::default() }).unwrap();
//! assert!(c > 0.0 && !w.align_a.is_empty());
//! ```
//!
#![cfg_attr(docsrs, feature(doc_cfg, doc_auto_cfg))]

pub mod common;
pub mod complex;
pub mod water;
pub mod needle;
pub mod stretcher;
pub mod est2genome;
pub mod needleall;
pub mod esim4;
pub mod matcher;

// Re-export the public API for ergonomic access.
pub use common::{EmbossersError, FastaRecord, parse_fasta, WaterMatrix};
pub use complex::{ComplexOptions, UjRow, compute_complexity, compute_complexity_multi, max_vocab};
pub use water::{WaterParams, WaterAlignment, water};
pub use needle::{NeedleParams, NeedleAlignment, needle};
pub use est2genome::{Est2GenomeParams, Est2GenomeAlignment, est2genome};

// needleall exports
pub use needleall::{needleall_pairs, NeedleAllResult};

// stretcher exports
pub use stretcher::{StretcherParams, StretcherAlignment, stretcher};

// esim4 exports
pub use esim4::{Esim4Params, Esim4Alignment, Exon, StrandMode, SpliceClass, esim4};

// matcher exports
pub use matcher::{MatcherParams, MatcherHit, matcher};
