# microbiorust-py
    
**Python bindings for [microBioRust](https://github.com/LCrossman/microBioRust) — a high-performance, modular bioinformatics toolkit written in Rust.**
 
microBioRust is the core Rust Crate
 
`microbiorust-py` provides fast and memory-efficient bioinformatics functionality to Python users by leveraging the power of Rust, exposed through [PyO3](https://github.com/PyO3/pyo3). This package aims to offer an alternative to libraries like Biopython, with a focus on speed, correctness, and extensibility.
 
---

## Features

- **Fast parsers** for GenBank and EMBL formats  
- Output to GFF, FAA, and FFN  
- Accurate feature extraction and translation  
- Seamless Python API for easy integration into existing pipelines  
- Built with Rust for safety and performance
- **Fast parsers** for BLAST in XML and tabular formats
- **Fast parser** for MSA alignments, subset, get_consensus
- Sequence metrics: hydrophobicity,amino_counts,amino_percentage
---

## Installation

Build the PyModule from source using `maturin` (recommended) - You will need to use the --features flag below:

```bash
pip install maturin
maturin develop --features extension-module
```
You can run
```bash
cargo test
```
which is the Rust test to see if the pyfunctions have been successfully added to the PyModule

Example usage in Python:

```python
import microbiorust
result = microbiorust.gbk_to_faa("test_input.gbk")
for r in result:
   print(r)
```

