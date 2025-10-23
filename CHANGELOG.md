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
