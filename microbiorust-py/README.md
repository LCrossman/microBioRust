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
- **Fast parser** for MSA alignments — subset, purge_gaps, get_consensus
- Write directly to GFF3, FAA, FNA and FFN formats
- Typed collections, return type enforces what data is returned
- Accurate feature extraction, gene, product, strand, start, stop, codon_start
- Native JSON serialization to instantly export extracted data structures to standard JSON strings
- Sequence metrics: hydrophobicity, amino acid counts and percentages
- Python API for easy integration into existing pipelines
- Built with Rust for memory safety and performance

---

## Modules

### `microbiorust gbk` — GenBank format

```python
import microbiorust as mb

# Write directly to file — most efficient for large files
# all functions are also available for embl format parse_embl, embl_to_faa, embl_to_ffn etc.)
collection = mb.parse_gbk("genome.gbk")
collection.write_faa("output.faa")
collection.write_ffn("output.ffn")
collection.write_fna("output.fna")

# Flat access across whole genome file — returns FaaCollection
faa = mb.gbk_to_faa("genome.gbk")
# print valid protein fasta
for info in faa.values():
    print(f">{info.locus_tag}\n{info.faa}")

# Per-contig record access
for record in collection.values():
    # prints the contig id and sequence
    print(record.id(), record.sequence())

    # protein sequences
    for info in record.faa().values():
        # prints the protein fasta for each predicted protein in the record
        print(f">{info.locus_tag}\n{info.faa}")

    # nucleotide sequences
    for info in record.ffn().values():
        # prints the nucleotide fasta sequence of each predicted gene
        print(f">{info.locus_tag}\n{info.ffn}")

    # features
    features = record.features()
    # prints the features of each predicted gene by locus tag key
    if "b3304" in features:
        feat = features["b3304"]
        print(f"Gene: {feat.gene}, Product: {feat.product}")
        print(f"Location: {feat.start}..{feat.stop}, Strand: {feat.strand}")

# Convert collection to JSON string
json_str = collection.to_json()
print(json_str)
# Parse JSON string into Python dictionary
data = json.loads(json_str)

# Count proteins without loading sequences
count = mb.gbk_to_faa_count("genome.gbk")

# Convert annotations from gbk or embl to GFF3
mb.gbk_to_gff("genome.gbk", dna=True)

---

### EMBL format: illustrates use by calling on the submodule, can also be called directly as mb.embl_to_faa etc.

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


### `microbiorust.blast` — BLAST results

```python
import microbiorust

results = microbiorust.parse_tabular("blast_results.tab")
for hit in results:
    print(hit["qseqid"], hit["pident"], hit["bitscore"])
```

## Choice of the usage pattern

| Goal | Use |
|---|---|
| Write everything to file | `collection.write_faa()` / `write_ffn()` / `write_fna()` |
| Get all proteins across a whole genome file | `gbk_to_faa()` |
| Work per genome contig record | `parse_gbk()` then `record.faa()` or `record.ffn()` |
| Features and sequences together | `parse_gbk()` then `record.sequences()` + `record.features()` |
| Count proteins without loading | `gbk_to_faa_count()` |
| Convert collection to JSON string | `collection.to_json()` |
| Parse JSON string into Python dictionary | `json.loads()` |
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
