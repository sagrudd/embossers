# embossers

EMBOSS-inspired bioinformatics utilities in Rust. This crate includes a CLI binary **`emboss`** with a feature-complete **`complex`** subcommand that computes *linguistic sequence complexity* for nucleotide sequences.

## Install

```toml
# Cargo.toml
embossers = "0.1"
```

## CLI usage

```text
emboss complex --sequence input.fasta --lwin 100 --step 5 --jmin 4 --jmax 6 --outfile results.tsv --print --ujtablefile complex.ujtable
# or, from stdin
cat input.fasta | emboss complex --lwin 100 --step 5 --jmin 4 --jmax 6 --outfile results.tsv
```

Key flags: `--sequence`, `--lwin`, `--step`, `--jmin`, `--jmax`, `--omnia`, `--sim`, `--freq`, `--print`, `--outfile`, `--ujtablefile`, `--outseq`.

## Library example

```rust
use embossers::{ComplexOptions, compute_complexity};
let seq = "ACGTACGTACGT";
let opts = ComplexOptions{ lwin: 12, step: 6, jmin: 2, jmax: 3, sim: 0, freq_weighted_sim: false };
let (c, _rows) = compute_complexity(seq, &opts).unwrap();
println!("complexity: {c}");
```

## SemVer

Follows [Semantic Versioning](https://semver.org/). Current version: **v0.1.4** (2025-10-22).
