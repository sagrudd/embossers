# embossers

EMBOSS-inspired bioinformatics utilities in Rust. This crate includes a CLI binary **`emboss`** with:

- **`complex`** — linguistic sequence complexity (EMBOSS `complex`)
- **`water`** — Smith–Waterman local alignment with affine gaps (EMBOSS `water`) — alias **`waters`**
- **`needle`** — Needleman–Wunsch global alignment with affine gaps (EMBOSS `needle`)

Each subcommand is maintained in its own file under `src/bin/emboss/` and each algorithm is in its own module under `src/`.

## Install

```toml
# Cargo.toml
embossers = "0.1"
```

## CLI usage

```text
# complex
emboss complex --sequence input.fasta --lwin 100 --step 5 --jmin 4 --jmax 6 --outfile complex.tsv

# water (local alignment)
emboss water --asequence a.fasta --bsequence b.fasta --matrix blosum62 --gapopen 10 --gapextend 0.5 --outfile water.txt

# needle (global alignment)
emboss needle --asequence a.fasta --bsequence b.fasta --matrix dna --match-score 1 --mismatch -1 --gapopen 10 --gapextend 0.5 --outfile needle.txt
```

## Library quick start

```rust
use embossers::{ComplexOptions, compute_complexity, water, WaterParams, WaterMatrix, needle, NeedleParams};

let (c, _rows) = compute_complexity("ACGTACGT", &ComplexOptions::default()).unwrap();

let wparams = WaterParams{ matrix: WaterMatrix::Blosum62, ..Default::default() };
let w = water("PAWHEAE", "HEAGAWGHEE", &wparams).unwrap();

let nparams = NeedleParams{ matrix: WaterMatrix::Dna{ match_score: 1, mismatch: -1 }, ..Default::default() };
let n = needle("GATTACA", "GCATGCU", &nparams).unwrap();
```

## SemVer

Follows [Semantic Versioning](https://semver.org/). Current version: **v0.1.10** (2025-10-23).
