//! # embossers
//!
//! EMBOSS-inspired utilities implemented in Rust.
//!
//! Algorithms (library + CLI):
//! - **`complex`**: Linguistic sequence complexity (EMBOSS `complex`)
//! - **`water`**: Smith–Waterman local alignment with affine gaps (EMBOSS `water`)
//! - **`needle`**: Needleman–Wunsch global alignment with affine gaps (EMBOSS `needle`)
//!
//! Each module is fully `rustdoc`-documented and tested with small examples.
//!
//! See the `emboss` binary for CLI usage; each subcommand is maintained in
//! its own Rust file under `src/bin/emboss/`.
//!
#![cfg_attr(docsrs, feature(doc_cfg, doc_auto_cfg))]

pub mod common;
pub mod complex;
pub mod water;
pub mod needle;

// Re-export the public API for ergonomic access.
pub use common::{EmbossersError, FastaRecord, parse_fasta, WaterMatrix};
pub use complex::{ComplexOptions, UjRow, compute_complexity, compute_complexity_multi, max_vocab};
pub use water::{WaterParams, WaterAlignment, water};
pub use needle::{NeedleParams, NeedleAlignment, needle};
