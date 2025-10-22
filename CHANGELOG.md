# Changelog

All notable changes to this project will be documented in this file.

## [0.1.4] - 2025-10-22
### Fixed
- Removed any stray characters at the start of `src/lib.rs`.
- Corrected doc comments to use outer `///` only within items.
- Removed unused `rand::prelude` import to avoid warning.
- Silenced `pt` unused variable by naming it `_pt`.

## [0.1.3] - 2025-10-22
### Fixed
- Removed stray character at start of `src/lib.rs`.
- Switched accidental inner doc comments to outer `///` lines.
- Added `rand` dependency; fixed imports.

## [0.1.2] - 2025-10-22
### Added
- Feature-complete `emboss complex` CLI implemented with Clap v4.
- Over-documented library (`rustdoc`) and examples.
- `compute_complexity_multi` supporting EMBOSS `-omnia` behaviour.
