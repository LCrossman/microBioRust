//!  The aim of this crate is to provide Microbiology friendly Rust functions for bioinformatics.
//!  
//!  
//!  With the genbank parser, you are able to parse a genbank format file, then write into gff3 format
//!
//!  It is also possible to print the DNA sequences extracted from the coding sequences (genes, ffn format),
//!  plus the protein fasta sequences (faa format).
//!
//!  Additionally, you can create new features and records and save them either in genbank or gff3 format
//!
#![allow(non_snake_case)]
pub mod record;
pub mod embl;
pub mod gbk;
