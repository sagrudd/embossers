# embossers

EMBOSS-inspired bioinformatics utilities in Rust. This crate provides a library **and** a CLI binary **`emboss`** implementing:

- **`complex`** — linguistic sequence complexity (EMBOSS `complex`)
- **`water`** — Smith–Waterman local alignment with affine gaps (EMBOSS `water`) — alias **`waters`**
- **`needle`** — Needleman–Wunsch global alignment with affine gaps (EMBOSS `needle`)

Each subcommand lives in `src/bin/emboss/`, and each algorithm lives in `src/`.

---

## Build & Install

### Prerequisites
- Rust toolchain (stable). If needed: <https://rustup.rs/>

### Build in place
```bash
git clone <your-repo> embossers
cd embossers
cargo build --release
```

### Run without installing
```bash
cargo run --bin emboss -- --help
```

### Install the CLI
```bash
cargo install --path .
# now `emboss` is on your PATH
```

---

## CLI Usage

### `complex` — linguistic complexity
Compute complexity per FASTA record (windows aggregated per record).

```bash
# One or more FASTA files
emboss complex   --sequence genome1.fa genome2.fa   --lwin 100 --step 5 --jmin 4 --jmax 6   --outfile complex.tsv

# Read FASTA from stdin
cat seqs.fa | emboss complex --lwin 200 --step 10 --jmin 3 --jmax 8 --outfile out.tsv

# Include random simulations with empirical base frequencies
emboss complex --sequence seqs.fa --sim 50 --freq --outfile complex.tsv
```

**Output:** a TSV summary with columns: `id, windows, lwin, step, jmin, jmax, complexity`.

---

### `water` — Smith–Waterman local alignment
```bash
emboss water   --asequence a.fa   --bsequence b.fa   --matrix blosum62   --gapopen 10.0 --gapextend 0.5   --outfile water.txt

# Alias also works
emboss waters --asequence a.fa --bsequence b.fa --outfile water.txt

# DNA scoring variant
emboss water --asequence a.fa --bsequence b.fa --matrix dna --match-score 2 --mismatch -1
```

**Output:** human-readable alignment file with score, identity/gaps, CIGAR and 60‑column blocks.

---

### `needle` — Needleman–Wunsch global alignment
---

### `needleall`
---

### `stretcher` — global alignment (linear-space, Myers–Miller)
Use a divide‑and‑conquer Needleman–Wunsch with affine gaps. Memory stays near `O(n+m)`.

```bash
emboss stretcher   --asequence a.fa   --bsequence b.fa   --matrix dna --match-score 1 --mismatch -1   --gapopen 10.0 --gapextend 0.5   --fallback-threshold 256   --outfile stretcher.txt
```

**Output:** human-readable alignment file with score, identity/gaps, and CIGAR.
 — many-to-many global alignments
Run Needleman–Wunsch alignments for **every pair** across two multi‑FASTA sets.
Outputs a summary TSV and (optionally) per‑pair alignment files.

```bash
emboss needleall   --aseqs setA.fa   --bseqs setB.fa   --matrix dna --match-score 1 --mismatch -1   --gapopen 10.0 --gapextend 0.5   --summary needleall.tsv   --outdir pairwise_alignments
```

**Summary TSV columns:** `a_id, b_id, score, pct_identity, pct_gaps, cigar_len`

---

### `est2genome`
---

### `esim4` — SIM4-like spliced EST ↔ genome alignment
A strand-aware wrapper around `est2genome` that tries both forward and reverse‑complement
(by default), classifies splice sites (GT‑AG, GC‑AG, AT‑AC), applies splice bonuses, and
reports an exon table.

```bash
emboss esim4   --genome genome.fa   --est read.fa   --matrix dna --match-score 2 --mismatch -1   --gapopen 10.0 --gapextend 0.5   --intron-min 20 --splice-bonus 5   --accept-gc-ag true --accept-at-ac false   --strand auto   --outfile esim4.txt
```

**Output:** Similar to `est2genome` but adds **strand**, **splice class counts**, and an **exon table**.
 — splice-aware EST ↔ genome alignment (semi‑global)
```bash
emboss est2genome   --genome genome.fa   --est read.fa   --matrix dna --match-score 2 --mismatch -1   --gapopen 10.0 --gapextend 0.5   --intron-min 20 --splice-bonus 5   --outfile est2genome.txt
```

**Output:** human-readable alignment file with score, identity/gaps, CIGAR (using `N` for introns), and an intron table with genomic coordinates and canonical `GT-AG` flags.

```bash
emboss needle   --asequence a.fa   --bsequence b.fa   --matrix dna --match-score 1 --mismatch -1   --gapopen 10.0 --gapextend 0.5   --outfile needle.txt

# Protein scoring with BLOSUM62
emboss needle --asequence a.fa --bsequence b.fa --matrix blosum62 --gapopen 10 --gapextend 0.5
```

**Output:** human-readable alignment file similar to `water`.

---

## Library Quick Start

```rust
use embossers::{ComplexOptions, compute_complexity, water, WaterParams, WaterMatrix, needle, NeedleParams};

let (c, rows) = compute_complexity("ACGTACGT", &ComplexOptions::default()).unwrap();
let w = water("PAWHEAE", "HEAGAWGHEE", &WaterParams{ matrix: WaterMatrix::Blosum62, ..Default::default() }).unwrap();
let n = needle("GATTACA", "GCATGCU", &NeedleParams{ matrix: WaterMatrix::Dna{ match_score: 1, mismatch: -1 }, ..Default::default() }).unwrap();
```

---

## Compatibility Notes
- Scoring and defaults are chosen to be close to EMBOSS. Output formatting is readable and compact; for byte‑for‑byte parity with EMBOSS reports, open an issue.

## License
MIT OR Apache‑2.0

## Version
**v0.1.20** (2025‑10‑23)
