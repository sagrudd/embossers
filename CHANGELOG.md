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
