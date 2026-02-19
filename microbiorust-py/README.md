# microbiorust 🦀

**Python bindings for [microBioRust](https://github.com/microBioRust/microBioRust) — a high-performance, modular bioinformatics toolkit written in Rust.**

`microbiorust` provides fast and memory-efficient bioinformatics functionality to Python users by leveraging the power of Rust, exposed through [PyO3](https://github.com/PyO3/pyo3). This package aims to offer an alternative to libraries like Biopython, with a focus on speed, correctness, and extensibility.

---

## Installation

```bash
pip install microbiorust
```

Wheels are available for Linux, macOS and Windows (Python 3.10+). No Rust toolchain required.
(no requirement to install Rust)

### Build from source

If you prefer to build from source using [maturin](https://maturin.rs):

```bash
pip install maturin
git clone https://github.com/microBioRust/microBioRust
cd microbiorust-py
maturin develop --features extension-module
```

To verify the Python module functions are correctly exposed from Rust:

```bash
cargo test
```

---

## Features

- **Fast parsers** for GenBank and EMBL formats
- **Fast parsers** for BLAST XML and tabular formats
- **Fast parser** for MSA alignments — subset, get_consensus
- Output to GFF3, FAA and FFN formats
- Accurate feature extraction and translation
- Sequence metrics: hydrophobicity, amino acid counts and percentages
- Seamless Python API for easy integration into existing pipelines
- Built with Rust for memory safety and performance

---

## Modules

### `microbiorust gbk` — GenBank format

```python
from microbiorust import gbk

# Extract protein sequences to FASTA
gbk.gbk_to_faa("input.gbk", "output.faa")

# Extract nucleotide sequences to FASTA
gbk.gbk_to_fna("input.gbk", "output.fna")

# Count protein sequences
count = gbk.gbk_to_faa_count("input.gbk")

# Convert annotations to GFF3
gbk.gbk_to_gff("input.gbk", "output.gff")
```

---

### `microbiorust embl` — EMBL format

```python
from microbiorust import embl

# Extract protein sequences to FASTA
embl.embl_to_faa("input.embl", "output.faa")

# Extract nucleotide sequences to FASTA
embl.embl_to_fna("input.embl", "output.fna")

# Convert annotations to GFF3
embl.embl_to_gff("input.embl", "output.gff")
```

---

### `microbiorust seqmetrics` — Sequence metrics

```python
from microbiorust import seqmetrics

sequence = "MKTLLLTLVVVTIVCLDLGAVGNGSSLSEDKDNVHK"

# Hydrophobicity score
window_size = 5
score = seqmetrics.hydrophobicity(sequence, window_size)

# Amino acid counts
counts = seqmetrics.amino_counts(sequence)

# Amino acid percentages
percentages = seqmetrics.amino_percentage(sequence)
```

---

### `microbiorust align` — Multiple sequence alignment

```python
from microbiorust import align

# Subset a fasta format MSA by row and column e.g.
align.subset_msa_alignment("input.fasta", "ids.txt", "output.fasta")
where the first tuple (0,10) is a row-wise subset and
the second tuple (0,100) is a column-wise subset
```

---

## Why Rust?

Rust gives microbiorust **C-level performance** with memory safety — no segfaults, no GIL limitations, and no need for NumPy or Pandas for core parsing operations. Large GenBank or EMBL files are parsed significantly faster than equivalent pure-Python implementations.

---

## Documentation

Full documentation: [https://microbiorust.github.io/docs/](https://microbiorust.github.io/docs/)

Source: [https://github.com/microBioRust/microBioRust](https://github.com/microBioRust/microBioRust)

---

## License

MIT
