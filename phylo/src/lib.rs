/*!
 The aim of this crate is to provide Microbiology friendly Rust functions for bioinformatics.

 With phylo module you should be able to

 1. Parse newick trees into Rust equivalent data structure
 2. Conversion of newick tree format to other file formats like phyloXML
 3. Neighbor joining, bootstrapping from many files, visualizing, rooting etc.
*/
#![allow(non_snake_case)]
pub mod parser;
