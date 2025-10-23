# Changelog

## [0.1.9] - 2025-10-23
### Fixed
- Corrected CLI layout: subcommand files moved under `src/bin/emboss/` as modules of the `emboss` binary to avoid extra bins and missing `main` errors.
- Removed an unused import warning in `complex.rs`.

## [0.1.8] - 2025-10-23
### Added
- **`needle`** (Needleman–Wunsch global alignment, affine gaps) implemented in library and CLI.
- Split all algorithms and CLI subcommands into their own Rust files for maintainability.

## [0.1.7] - 2025-10-23
### Changed
- Refreshed implementation; removed all stubs. Added CLI alias `waters` for `water`.

## [0.1.10] - 2025-10-23
### Fixed
- Adjusted `src/bin/emboss.rs` to use `#[path = "emboss/..." ]` module attributes so subcommands in `src/bin/emboss/` compile correctly.

## [0.1.11] - 2025-10-23
### Documentation
- Expanded rustdoc throughout (`common`, `complex`, `water`, `needle`) with algorithm notes, examples, and complexity.
- README updated with build instructions and usage examples for all subcommands.

## [0.1.12] - 2025-10-23
### Fixed
- Escaped/adjusted doc examples in `src/common.rs` to avoid unescaped quotes/newline escapes that some compilers interpret during attribute expansion.

## [0.1.13] - 2025-10-23
### Added
- **est2genome**: splice-aware EST↔genome alignment (semi-global on genome) with intron detection and canonical `GT-AG` annotation. Available in library and CLI.

## [0.1.14] - 2025-10-23
### Added
- **needleall**: many-to-many Needleman–Wunsch (global) alignments across two FASTA sets, summary TSV plus optional per-pair alignment files. CLI subcommand `emboss needleall`.

## [0.1.15] - 2025-10-23
### Fixed
- `needleall` CLI: corrected `sanitize()` (replaced Pythonic `in` with `matches!`).
- `est2genome`: silenced an unused parameter warning by prefixing with `_`.

## [0.1.16] - 2025-10-23
### Added
- **stretcher**: Myers–Miller/Hirschberg-style global alignment (linear space) with affine gaps.
  Library function `stretcher()` and CLI `emboss stretcher`.

## [0.1.17] - 2025-10-23
### Fixed
- `stretcher`: removed unused variables/assignments, corrected minor style warnings, and renamed internal variables to snake_case.

## [0.1.18] - 2025-10-23
### Added
- **esim4**: SIM4-like spliced EST→genome alignment built on top of `est2genome`, with strand auto-detection, expanded splice site classification (GT-AG, GC-AG, AT-AC), exon table, and CLI subcommand `emboss esim4`.

## [0.1.19] - 2025-10-23
### Fixed
- `stretcher`: corrected identifier shadowing introduced by a prior refactor; reimplemented small-DP branch with distinct matrix names to avoid clashing with `m` (length). Also removed stray unused variable.
- `esim4`: removed an unused import and tiny style cleanups.

## [0.1.20] - 2025-10-23
### Fixed
- `stretcher`: restored missing fallback `NeedleParams` variable (`np`) used when problem size is small.
- `esim4`: renamed splice enum variants to UpperCamelCase (`GtAg`, `GcAg`, `AtAc`) and updated matches to silence warnings.

## [0.1.21] - 2025-10-23
### Added
- **matcher**: Waterman–Eggert-style local alignments (declumped), returning the top `--alternatives` non-overlapping local alignments. Library API and CLI subcommand `emboss matcher`.

## [0.1.22] - 2025-10-23
### Added
- **seqmatchall**: all-vs-all Waterman–Eggert local alignments across a multi-FASTA set, built on the `matcher` engine. CLI subcommand `emboss seqmatchall` writes a summary TSV and optional per-pair alignment files.

## [0.1.23] - 2025-10-23
### Fixed
- `seqmatchall` CLI: avoid moving `MatcherParams` by cloning into `SeqMatchAllParams` so it can be reused for per-pair files.
- `matcher`: removed stray unused variables and silenced minor warnings in local coordinate computation.

## [0.1.24] - 2025-10-23
### Added
- **supermatcher**: seed-and-extend local alignments (word hits + SW) across two sets.
- **wordfinder**: query vs set, seed-and-extend local alignments (word hits + SW).
- **wordmatch**: exact-match regions of minimum word size between two sequences; optional GFF outputs.
- **seqalign**: extend an existing multiple alignment with sequences by consensus-guided profile alignment.

## [0.1.25] - 2025-10-23
### Fixed
- `seqalign`: replaced nonexistent `EmbossersError::ParseError` with `InvalidSequence` and tidied warnings (unused param, needless `mut`).

## [0.1.26] - 2025-10-23
### Added
- **domainalign**: progressive, consensus-guided multiple alignment (EMBASSY domalign-like), via iterative global alignment to consensus.
- **domainrep**: within-sequence repeat finder using k-mer seeding and maximal right extension; reports repeat blocks and optional GFF.
- **oalistat**: alignment statistics for gapped FASTA / simple alignment files (sequences, columns, %gap, pairwise identity, conserved columns).
