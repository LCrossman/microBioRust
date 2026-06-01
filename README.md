<img src="assets/MICROBIO B.svg#gh-light-mode-only" alt="microBioRust logo light mode" width="200">
<img src="assets/BIO W.png#gh-dark-mode-only" alt="microBioRust logo dark mode" width="200">

[![Docs](https://img.shields.io/badge/docs-mkdocs-blue.svg)](https://lcrossman.github.io/microBioRust-docs)

[![Crates.io Version](https://img.shields.io/crates/v/microBioRust?style=flat)](https://crates.io/crates/microBioRust)

See the microbiorust-py benchmarks with airspeed velocity:  
[![Benchmarks](https://img.shields.io/badge/benchmarks-passing-brightgreen)](https://lcrossman.github.io/microBioRust/)

## A Rust bioinformatics crate aimed at Microbial genomics<br>

The aim of this crate is to provide Microbiology friendly Rust functions for bioinformatics.<br>

Some concepts with many thanks to Rust-bio<br>
Please see the Roadmap for futher details [here](ROADMAP.md)

Check out the [docs here](https://microBioRust.github.io/microBioRust)

To install Rust - please see here [Rust install](https://www.rust-lang.org/tools/install) or with Conda<br>
If you would like to contribute please follow the [Rust code of conduct](https://www.rust-lang.org/policies/code-of-conduct)

Questions and comments - please join the Discord server :) [here](https://discord.gg/xP2ngwTttz)


Currently there is functionality for:<br>
````
 1. A Genbank parser and conversions to GFF, protein fasta (faa), DNA sequence (.fna), Gene sequence (.ffn) parser, JSON

 2. An Embl parser and conversions to GFF or GBK , .faa, .fna, .ffn, JSON

 3. Parsing of multiple sequence alignments

 4. Parsing BLAST or Diamond/MMSeq2 output formats in XML/Tabular 

 5. Calculate many sequence metrics e.g. hydrophobicity, distance measures

 6. A demo Heatmap plot with wasm and d3.js

````

To see more on how to use the project please have a look at usage: [here](docs/usage.md)

To use a specific workspace (such as microBioRust, seqmetrics, microbiorust-py or heatmap) clone the project, cd into the specific directory required and build the project from there.

For more background please see <https://LCrossman.github.io/microBioRust_details>
