#![feature(prelude_import)]
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
#[prelude_import]
use std::prelude::rust_2021::*;
#[macro_use]
extern crate std;
pub mod embl {
    //! # An EMBL format to GFF parser
    //!
    //!
    //! You are able to parse genbank and save as a GFF (gff3) format as well as extracting DNA sequences, gene DNA sequences (ffn) and protein fasta sequences (faa)
    //!
    //! You can also create new records and save as a embl (gbk) format
    //!
    //! ## Detailed Explanation
    //!
    //!
    //! The Embl parser contains:
    //!
    //! Records - a top level structure which consists of either one record (single embl) or multiple instances of record (multi-embl).
    //!
    //! Each Record contains:
    //!
    //! 1. A source, ```SourceAttributes```, construct(enum) of counter (source name), start, stop [of source or contig], organism, mol_type, strain, type_material, db_xref
    //! 2. Features, ```FeatureAttributes```, construct(enum) of counter (locus tag), gene (if present), product, codon start, strand, start, stop [of cds/gene]
    //! 3. Sequence features, ```SequenceAttributes```, construct(enum) of counter (locus tag), sequence_ffn (DNA gene sequence) sequence_faa (protein translation), strand, codon start, start, stop [cds/gene]
    //! 4. The DNA sequence of the whole record (or contig)
    //!
    //!  Example to extract and print all the protein sequence fasta, example using getters (or get_ functionality), simplified embl! macro
    //!
    //!```rust
    //! use clap::Parser;
    //! use std::fs::File;
    //! use microBioRust::embl::Reader;
    //! use std::io;
    //! use microBioRust::embl;
    //!
    //!
    //! #[derive(Parser, Debug)]
    //! #[clap(author, version, about)]
    //! struct Arguments {
    //! #[clap(short, long)]
    //! filename: String,
    //! }
    //!
    //! pub fn genbank_to_faa() -> Result<(), anyhow::Error> {
    //!            let args = Arguments::parse();
    //!            let records = embl!(&args.filename);
    //!            for record in records {
    //!               for (k, v) in &record.cds.attributes {
    //!                  if let Some(seq) = record.seq_features.get_sequence_faa(k) {
    //!                     println!(">{}|{}\n{}", &record.id, &k, seq);
    //!                     }
    //!                  }
    //!            }
    //!            return Ok(());
    //!  }
    //!```
    //!
    //! Example to extract protein sequence from embl file, debugging use
    //!```rust
    //! use clap::Parser;
    //! use std::fs::File;
    //! use microBioRust::embl::Reader;
    //! use std::io;
    //!
    //! #[derive(Parser, Debug)]
    //! #[clap(author, version, about)]
    //! struct Arguments {
    //! #[clap(short, long)]
    //! filename: String,
    //! }
    //!
    //! pub fn embl_to_faa() -> Result<(), anyhow::Error> {
    //!            let args = Arguments::parse();
    //!            let file_embl = File::open(args.filename)?;
    //!            let mut reader = Reader::new(file_embl);
    //!            let mut records = reader.records();
    //!            loop {
    //!                //collect from each record advancing on a next record basis, count cds records
    //!                match records.next() {	
    //!                    Some(Ok(mut record)) => {
    //!		                     for (k, v) in &record.cds.attributes {
    //!		                         match record.seq_features.get_sequence_faa(&k) {
    //!		                             Some(value) => { let seq_faa = value.to_string();
    //!				                              println!(">{}|{}\n{}", &record.id, &k, seq_faa);
    //!						              },
    //!				             _ => (),
    //!				             };
    //!		                         }
    //!                                      },
    //!	               Some(Err(e)) => { println!("Error encountered - an err {:?}", e); },
    //!	               None => break,
    //!	               }
    //!                 }
    //!            return Ok(());
    //!  }
    //!```
    //!
    //!
    //!  Example to save a provided multi- or single genbank file as a GFF file (by joining any multi-genbank)
    //!
    //!
    //! ```rust
    //!    use microBioRust::embl::gff_write;
    //!    use microBioRust::embl::Reader;
    //!    use microBioRust::embl::Record;
    //!    use std::collections::BTreeMap;
    //!    use std::fs::File;
    //!    use clap::Parser;
    //!    use std::io;
    //!
    //!   #[derive(Parser, Debug)]
    //!   #[clap(author, version, about)]
    //!   struct Arguments {
    //!   #[clap(short, long)]
    //!   filename: String,
    //!   }
    //!
    //!    pub fn embl_to_gff() -> io::Result<()> {
    //!        let args = Arguments::parse();
    //!        let file_embl = File::open(&args.filename)?;
    //!        let prev_start: u32 = 0;
    //!        let mut prev_end: u32 = 0;
    //!        let mut reader = Reader::new(file_embl);
    //!        let mut records = reader.records();
    //!        let mut read_counter: u32 = 0;
    //!        let mut seq_region: BTreeMap<String, (u32,u32)> = BTreeMap::new();
    //!        let mut record_vec: Vec<Record> = Vec::new();
    //!        loop {
    //!            match records.next() {	
    //!                Some(Ok(mut record)) => {
    //!	               //println!("next record");
    //!                    //println!("Record id: {:?}", record.id);
    //!		       let source = record.source_map.source_name.clone().expect("issue collecting source name");
    //!		       let beginning = match record.source_map.get_start(&source) {
    //!		                        Some(value) => value.get_value(),
    //!				        _ => 0,
    //!					};
    //!		       let ending = match record.source_map.get_stop(&source) {
    //!		                        Some(value) => value.get_value(),
    //!					_ => 0,
    //!					};
    //!		       if ending + prev_end < beginning + prev_end {
    //!		          println!("debug: end value smaller is than the start {:?}", beginning);
    //!		          }
    //!		       seq_region.insert(source, (beginning + prev_end, ending + prev_end));
    //!		       record_vec.push(record);
    //!                    // Add additional fields to print if needed
    //!		       read_counter+=1;
    //!		       prev_end+=ending; // create the joined record if there are multiple
    //!                    },
    //!	            Some(Err(e)) => { println!("theres an err {:?}", e); },
    //!	            None => {
    //!	               println!("finished iteration");
    //!	                     break; },
    //!	            }
    //!            }
    //!        let output_file = format!("{}.gff", &args.filename);
    //!        gff_write(seq_region.clone(), record_vec, &output_file, true);
    //!        println!("Total records processed: {}", read_counter);
    //!        return Ok(());
    //!    }
    //!```
    //! Example to create a completely new record, use of setters or set_ functionality
    //!
    //! To write into GFF format requires gff_write(seq_region, record_vec, filename, true or false)
    //!
    //! The seq_region is the region of interest to save with name and DNA coordinates such as ``` seqregion.entry("source_1".to_string(), (1,897))```
    //! This makes it possible to save the whole file or to subset it
    //!
    //! record_vec is a list of the records.  If there is only one record, include this as a vec using ``` vec![record] ```
    //!
    //! The boolean true/false describes whether the DNA sequence should be included in the GFF3 file
    //!
    //! To write into embl format requires embl_write(seq_region, record_vec, filename), no true or false since embl format will include the DNA sequence
    //!
    //!
    //! ```rust
    //!    use microBioRust::embl::gff_write;
    //!    use microBioRust::embl::RangeValue;
    //!    use microBioRust::embl::Record;
    //!    use std::collections::BTreeMap;
    //!
    //!     pub fn create_new_record() -> Result<(), anyhow::Error> {
    //!         let filename = format!("new_record.gff");
    //!	    let mut record = Record::new();
    //!	    let mut seq_region: BTreeMap<String, (u32,u32)> = BTreeMap::new();
    //!         //example from E.coli K12
    //!	    seq_region.insert("source_1".to_string(), (1,897));
    //!         //Add the source into SourceAttributes
    //!         record.source_map
    //!	         .set_counter("source_1".to_string())
    //!	         .set_start(RangeValue::Exact(1))
    //!	         .set_stop(RangeValue::Exact(897))
    //!	         .set_organism("Escherichia coli".to_string())
    //!	         .set_mol_type("DNA".to_string())
    //!	         .set_strain("K-12 substr. MG1655".to_string())
    //!		 .set_type_material("type strain of Escherichia coli K12".to_string())
    //!	         .set_db_xref("PRJNA57779".to_string());
    //!         //Add the features into FeatureAttributes, here we are setting two features, i.e. coding sequences or genes
    //!	    record.cds
    //!                  .set_counter("b3304".to_string())
    //!                  .set_start(RangeValue::Exact(1))
    //!                  .set_stop(RangeValue::Exact(354))
    //!                  .set_gene("rplR".to_string())
    //!                  .set_product("50S ribosomal subunit protein L18".to_string())
    //!                  .set_codon_start(1)
    //!                  .set_strand(-1);
    //!	    record.cds
    //!                  .set_counter("b3305".to_string())
    //!                  .set_start(RangeValue::Exact(364))
    //!                  .set_stop(RangeValue::Exact(897))
    //!                  .set_gene("rplF".to_string())
    //!                  .set_product("50S ribosomal subunit protein L6".to_string())
    //!                  .set_codon_start(1)
    //!                  .set_strand(-1);
    //!         //Add the sequences for the coding sequence (CDS) into SequenceAttributes
    //!	    record.seq_features
    //!	         .set_counter("b3304".to_string())
    //!		 .set_start(RangeValue::Exact(1))
    //!                 .set_stop(RangeValue::Exact(354))
    //!                 .set_sequence_ffn("ATGGATAAGAAATCTGCTCGTATCCGTCGTGCGACCCGCGCACGCCGCAAGCTCCAGGAG
    //!CTGGGCGCAACTCGCCTGGTGGTACATCGTACCCCGCGTCACATTTACGCACAGGTAATT
    //!GCACCGAACGGTTCTGAAGTTCTGGTAGCTGCTTCTACTGTAGAAAAAGCTATCGCTGAA
    //!CAACTGAAGTACACCGGTAACAAAGACGCGGCTGCAGCTGTGGGTAAAGCTGTCGCTGAA
    //!CGCGCTCTGGAAAAAGGCATCAAAGATGTATCCTTTGACCGTTCCGGGTTCCAATATCAT
    //!GGTCGTGTCCAGGCACTGGCAGATGCTGCCCGTGAAGCTGGCCTTCAGTTCTAA".to_string())
    //!                 .set_sequence_faa("MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAE
    //!QLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF".to_string())
    //!                 .set_codon_start(1)
    //!                 .set_strand(-1);
    //!	    record.seq_features
    //!	         .set_counter("bb3305".to_string())
    //!		 .set_start(RangeValue::Exact(364))
    //!                 .set_stop(RangeValue::Exact(897))
    //!                 .set_sequence_ffn("ATGTCTCGTGTTGCTAAAGCACCGGTCGTTGTTCCTGCCGGCGTTGACGTAAAAATCAAC
    //!GGTCAGGTTATTACGATCAAAGGTAAAAACGGCGAGCTGACTCGTACTCTCAACGATGCT
    //!GTTGAAGTTAAACATGCAGATAATACCCTGACCTTCGGTCCGCGTGATGGTTACGCAGAC
    //!GGTTGGGCACAGGCTGGTACCGCGCGTGCCCTGCTGAACTCAATGGTTATCGGTGTTACC
    //!GAAGGCTTCACTAAGAAGCTGCAGCTGGTTGGTGTAGGTTACCGTGCAGCGGTTAAAGGC
    //!AATGTGATTAACCTGTCTCTGGGTTTCTCTCATCCTGTTGACCATCAGCTGCCTGCGGGT
    //!ATCACTGCTGAATGTCCGACTCAGACTGAAATCGTGCTGAAAGGCGCTGATAAGCAGGTG
    //!ATCGGCCAGGTTGCAGCGGATCTGCGCGCCTACCGTCGTCCTGAGCCTTATAAAGGCAAG
    //!GGTGTTCGTTACGCCGACGAAGTCGTGCGTACCAAAGAGGCTAAGAAGAAGTAA".to_string())
    //!                 .set_sequence_faa("MSRVAKAPVVVPAGVDVKINGQVITIKGKNGELTRTLNDAVEVKHADNTLTFGPRDGYAD
    //!GWAQAGTARALLNSMVIGVTEGFTKKLQLVGVGYRAAVKGNVINLSLGFSHPVDHQLPAG
    //!ITAECPTQTEIVLKGADKQVIGQVAADLRAYRRPEPYKGKGVRYADEVVRTKEAKKK".to_string())
    //!                 .set_codon_start(1)
    //!                 .set_strand(-1);
    //!         //Add the full sequence of the entire record into the record.sequence
    //!	    record.sequence = "TTAGAACTGAAGGCCAGCTTCACGGGCAGCATCTGCCAGTGCCTGGACACGACCATGATA
    //!TTGGAACCCGGAACGGTCAAAGGATACATCTTTGATGCCTTTTTCCAGAGCGCGTTCAGC
    //!GACAGCTTTACCCACAGCTGCAGCCGCGTCTTTGTTACCGGTGTACTTCAGTTGTTCAGC
    //!GATAGCTTTTTCTACAGTAGAAGCAGCTACCAGAACTTCAGAACCGTTCGGTGCAATTAC
    //!CTGTGCGTAAATGTGACGCGGGGTACGATGTACCACCAGGCGAGTTGCGCCCAGCTCCTG
    //!GAGCTTGCGGCGTGCGCGGGTCGCACGACGGATACGAGCAGATTTCTTATCCATAGTGTT
    //!ACCTTACTTCTTCTTAGCCTCTTTGGTACGCACGACTTCGTCGGCGTAACGAACACCCTT
    //!GCCTTTATAAGGCTCAGGACGACGGTAGGCGCGCAGATCCGCTGCAACCTGGCCGATCAC
    //!CTGCTTATCAGCGCCTTTCAGCACGATTTCAGTCTGAGTCGGACATTCAGCAGTGATACC
    //!CGCAGGCAGCTGATGGTCAACAGGATGAGAGAAACCCAGAGACAGGTTAATCACATTGCC
    //!TTTAACCGCTGCACGGTAACCTACACCAACCAGCTGCAGCTTCTTAGTGAAGCCTTCGGT
    //!AACACCGATAACCATTGAGTTCAGCAGGGCACGCGCGGTACCAGCCTGTGCCCAACCGTC
    //!TGCGTAACCATCACGCGGACCGAAGGTCAGGGTATTATCTGCATGTTTAACTTCAACAGC
    //!ATCGTTGAGAGTACGAGTCAGCTCGCCGTTTTTACCTTTGATCGTAATAACCTGACCGTT
    //!GATTTTTACGTCAACGCCGGCAGGAACAACGACCGGTGCTTTAGCAACACGAGACAT".to_string();
    //!           gff_write(seq_region, vec![record], &filename, true);
    //!	   return Ok(());
    //!      }
    //!```
    //!
    use std::io::{self, Write};
    use std::fs;
    use regex::Regex;
    use std::vec::Vec;
    use std::str;
    use std::convert::AsRef;
    use protein_translate::translate;
    use std::path::Path;
    use bio::alphabets::dna::revcomp;
    use anyhow::anyhow;
    use std::collections::BTreeMap;
    use std::fs::{OpenOptions, File};
    use anyhow::Context;
    use std::collections::HashSet;
    use paste::paste;
    use std::convert::TryInto;
    use chrono::prelude::*;
    /// import macro to create get_ functions for the values
    use crate::create_getters;
    /// import macro to create the set_ functions for the values in a Builder format
    use crate::create_builder;
    /// An EMBL reader.
    pub struct Records<B>
    where
        B: io::BufRead,
    {
        reader: Reader<B>,
        error_has_occurred: bool,
    }
    #[automatically_derived]
    impl<B: ::core::fmt::Debug> ::core::fmt::Debug for Records<B>
    where
        B: io::BufRead,
    {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "Records",
                "reader",
                &self.reader,
                "error_has_occurred",
                &&self.error_has_occurred,
            )
        }
    }
    impl<B> Records<B>
    where
        B: io::BufRead,
    {
        #[allow(unused_mut)]
        pub fn new(mut reader: Reader<B>) -> Self {
            Records {
                reader: reader,
                error_has_occurred: false,
            }
        }
    }
    impl<B> Iterator for Records<B>
    where
        B: io::BufRead,
    {
        type Item = Result<Record, anyhow::Error>;
        fn next(&mut self) -> Option<Result<Record, anyhow::Error>> {
            if self.error_has_occurred {
                {
                    ::std::io::_print(
                        format_args!("error was encountered in iteration\n"),
                    );
                };
                None
            } else {
                let mut record = Record::new();
                match self.reader.read(&mut record) {
                    Ok(_) => if record.is_empty() { None } else { Some(Ok(record)) }
                    Err(err) => {
                        self.error_has_occurred = true;
                        Some(
                            Err(
                                ::anyhow::Error::msg(
                                    ::alloc::__export::must_use({
                                        let res = ::alloc::fmt::format(
                                            format_args!("next record read error {0:?}", err),
                                        );
                                        res
                                    }),
                                ),
                            ),
                        )
                    }
                }
            }
        }
    }
    pub trait EmblRead {
        fn read(&mut self, record: &mut Record) -> Result<Record, anyhow::Error>;
    }
    ///per line reader for the file
    pub struct Reader<B> {
        reader: B,
        line_buffer: String,
    }
    #[automatically_derived]
    impl<B: ::core::fmt::Debug> ::core::fmt::Debug for Reader<B> {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "Reader",
                "reader",
                &self.reader,
                "line_buffer",
                &&self.line_buffer,
            )
        }
    }
    #[automatically_derived]
    impl<B: ::core::default::Default> ::core::default::Default for Reader<B> {
        #[inline]
        fn default() -> Reader<B> {
            Reader {
                reader: ::core::default::Default::default(),
                line_buffer: ::core::default::Default::default(),
            }
        }
    }
    impl Reader<io::BufReader<fs::File>> {
        /// Read Embl from given file path in given format.
        pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(
            path: P,
        ) -> anyhow::Result<Self> {
            fs::File::open(&path)
                .map(Reader::new)
                .with_context(|| ::alloc::__export::must_use({
                    let res = ::alloc::fmt::format(
                        format_args!("Failed to read Embl from {0:#?}", path),
                    );
                    res
                }))
        }
    }
    impl<R> Reader<io::BufReader<R>>
    where
        R: io::Read,
    {
        pub fn new(reader: R) -> Self {
            Reader {
                reader: io::BufReader::new(reader),
                line_buffer: String::new(),
            }
        }
    }
    impl<B> Reader<B>
    where
        B: io::BufRead,
    {
        pub fn from_bufread(bufreader: B) -> Self {
            Reader {
                reader: bufreader,
                line_buffer: String::new(),
            }
        }
        pub fn records(self) -> Records<B> {
            Records {
                reader: self,
                error_has_occurred: false,
            }
        }
    }
    ///main embl parser
    impl<'a, B> EmblRead for Reader<B>
    where
        B: io::BufRead,
    {
        #[allow(unused_mut)]
        #[allow(unused_variables)]
        #[allow(unused_assignments)]
        fn read(&mut self, record: &mut Record) -> Result<Record, anyhow::Error> {
            record.rec_clear();
            let mut sequences = String::new();
            let mut source_map = SourceAttributeBuilder::new();
            let mut cds = FeatureAttributeBuilder::new();
            let mut seq_features = SequenceAttributeBuilder::new();
            let mut cds_counter: i32 = 0;
            let mut source_counter: i32 = 0;
            let mut prev_end: u32 = 0;
            let mut organism = String::new();
            let mut mol_type = String::new();
            let mut strain = String::new();
            let mut source_name = String::new();
            let mut type_material = String::new();
            let mut theend: u32 = 0;
            let mut thestart: u32 = 0;
            let mut db_xref = String::new();
            if self.line_buffer.is_empty() {
                self.reader.read_line(&mut self.line_buffer)?;
                if self.line_buffer.is_empty() {
                    return Ok(record.to_owned());
                }
            }
            'outer: while !self.line_buffer.is_empty() {
                if self.line_buffer.starts_with("ID") {
                    record.rec_clear();
                    let mut header_fields: Vec<&str> = self
                        .line_buffer
                        .split_whitespace()
                        .collect();
                    let header_len = header_fields.len();
                    let mut header_iter = header_fields.iter();
                    header_iter.next();
                    record.id = header_iter
                        .next()
                        .ok_or_else(|| ::anyhow::__private::must_use({
                            let error = ::anyhow::__private::format_err(
                                format_args!("missing record id"),
                            );
                            error
                        }))?
                        .to_string();
                    if record.id.ends_with(";") {
                        record.id.pop();
                    }
                    header_iter.next();
                    header_iter.next();
                    header_iter.next();
                    header_iter.next();
                    header_iter.next();
                    header_iter.next();
                    header_iter.next();
                    let lens = header_iter
                        .next()
                        .ok_or_else(|| ::anyhow::__private::must_use({
                            let error = ::anyhow::__private::format_err(
                                format_args!("missing record length"),
                            );
                            error
                        }))?
                        .to_string();
                    record.length = lens.trim().parse::<u32>()?;
                    self.line_buffer.clear();
                }
                if self.line_buffer.starts_with("FT   source") {
                    let re = Regex::new(r"([0-9]+)[[:punct:]]+([0-9]+)")?;
                    let location = re
                        .captures(&self.line_buffer)
                        .ok_or_else(|| ::anyhow::__private::must_use({
                            let error = ::anyhow::__private::format_err(
                                format_args!("missing location"),
                            );
                            error
                        }))?;
                    let start = &location[1];
                    let end = &location[2];
                    thestart = start.trim().parse::<u32>()?;
                    source_counter += 1;
                    source_name = ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("source_{0}_{1}", record.id, source_counter),
                            );
                            res
                        })
                        .to_string();
                    thestart += prev_end;
                    theend = end.trim().parse::<u32>()? + prev_end;
                    loop {
                        self.line_buffer.clear();
                        self.reader.read_line(&mut self.line_buffer)?;
                        if self.line_buffer.starts_with("FT   CDS") {
                            record
                                .source_map
                                .set_counter(source_name.to_string())
                                .set_start(RangeValue::Exact(thestart))
                                .set_stop(RangeValue::Exact(theend))
                                .set_organism(organism.clone())
                                .set_mol_type(mol_type.clone())
                                .set_strain(strain.clone())
                                .set_type_material(type_material.clone())
                                .set_db_xref(db_xref.clone());
                            continue 'outer;
                        }
                        if self.line_buffer.contains("/organism") {
                            let org: Vec<&str> = self.line_buffer.split('\"').collect();
                            organism = org[1].to_string();
                        }
                        if self.line_buffer.contains("/mol_type") {
                            let mol: Vec<&str> = self.line_buffer.split('\"').collect();
                            mol_type = mol[1].to_string();
                        }
                        if self.line_buffer.contains("/strain") {
                            let stra: Vec<&str> = self.line_buffer.split('\"').collect();
                            strain = stra[1].to_string();
                        }
                        if self.line_buffer.contains("/type_material") {
                            let mat: Vec<&str> = self.line_buffer.split('\"').collect();
                            type_material = mat[1].to_string();
                        }
                        if self.line_buffer.contains("/db_xref") {
                            let db: Vec<&str> = self.line_buffer.split('\"').collect();
                            db_xref = db[1].to_string();
                        }
                    }
                }
                if self.line_buffer.starts_with("FT   CDS") {
                    let mut startiter: Vec<_> = Vec::new();
                    let mut enditer: Vec<_> = Vec::new();
                    let mut thestart: u32 = 0;
                    let mut thend: u32 = 0;
                    let mut joined: bool = false;
                    let joined = if self.line_buffer.contains("join") {
                        true
                    } else {
                        false
                    };
                    let re = Regex::new(r"([0-9]+)[[:punct:]]+([0-9]+)")?;
                    for cap in re.captures_iter(&self.line_buffer) {
                        cds_counter += 1;
                        thestart = cap[1]
                            .parse()
                            .expect("failed to match and parse numerical start");
                        theend = cap[2]
                            .parse()
                            .expect("failed to match and parse numerical end");
                        startiter.push(thestart);
                        enditer.push(theend);
                    }
                    let mut gene = String::new();
                    let mut product = String::new();
                    let strand: i8 = if self.line_buffer.contains("complement") {
                        -1
                    } else {
                        1
                    };
                    let mut locus_tag = String::new();
                    let mut codon_start: u8 = 1;
                    loop {
                        self.line_buffer.clear();
                        self.reader.read_line(&mut self.line_buffer)?;
                        if self.line_buffer.contains("/locus_tag=") {
                            let loctag: Vec<&str> = self
                                .line_buffer
                                .split('\"')
                                .collect();
                            locus_tag = loctag[1].to_string();
                        }
                        if self.line_buffer.contains("/codon_start") {
                            let codstart: Vec<&str> = self
                                .line_buffer
                                .split('=')
                                .collect();
                            let valstart = codstart[1].trim().parse::<u8>()?;
                            codon_start = valstart;
                        }
                        if self.line_buffer.contains("/gene=") {
                            let gen: Vec<&str> = self.line_buffer.split('\"').collect();
                            gene = gen[1].to_string();
                        }
                        if self.line_buffer.contains("/product") {
                            let prod: Vec<&str> = self.line_buffer.split('\"').collect();
                            product = substitute_odd_punctuation(prod[1].to_string())?;
                        }
                        if self.line_buffer.starts_with("FT   CDS")
                            || self.line_buffer.starts_with("SQ   Sequence")
                            || self.line_buffer.starts_with("FT   intron")
                            || self.line_buffer.starts_with("FT   exon")
                            || self.line_buffer.starts_with("     misc_feature")
                        {
                            if locus_tag.is_empty() {
                                locus_tag = ::alloc::__export::must_use({
                                        let res = ::alloc::fmt::format(
                                            format_args!("CDS_{0}", cds_counter),
                                        );
                                        res
                                    })
                                    .to_string();
                            }
                            if joined {
                                for (i, m) in startiter.iter().enumerate() {
                                    let loc_tag = ::alloc::__export::must_use({
                                        let res = ::alloc::fmt::format(
                                            format_args!("{0}_{1}", locus_tag.clone(), i),
                                        );
                                        res
                                    });
                                    record
                                        .cds
                                        .set_counter(loc_tag)
                                        .set_start(RangeValue::Exact(*m))
                                        .set_stop(RangeValue::Exact(enditer[i]))
                                        .set_gene(gene.to_string())
                                        .set_product(product.to_string())
                                        .set_codon_start(codon_start)
                                        .set_strand(strand);
                                }
                                continue 'outer;
                            } else {
                                record
                                    .cds
                                    .set_counter(locus_tag.clone())
                                    .set_start(RangeValue::Exact(thestart))
                                    .set_stop(RangeValue::Exact(theend))
                                    .set_gene(gene.to_string())
                                    .set_product(product.to_string())
                                    .set_codon_start(codon_start)
                                    .set_strand(strand);
                                continue 'outer;
                            }
                        }
                    }
                }
                if self.line_buffer.starts_with("SQ   Sequence") {
                    let mut sequences = String::new();
                    let result_seq = loop {
                        self.line_buffer.clear();
                        self.reader.read_line(&mut self.line_buffer)?;
                        if self.line_buffer.starts_with("//") {
                            break sequences;
                        } else {
                            let s: Vec<&str> = self
                                .line_buffer
                                .split_whitespace()
                                .collect();
                            let sequence = if s.len() > 1 {
                                s[0..s.len() - 1].join("")
                            } else {
                                String::new()
                            };
                            sequences.push_str(&sequence);
                        }
                    };
                    record.sequence = result_seq.to_string();
                    let mut iterablecount: u32 = 0;
                    for (key, val) in record.cds.iter_sorted() {
                        let (
                            mut a,
                            mut b,
                            mut c,
                            mut d,
                        ): (Option<u32>, Option<u32>, Option<i8>, Option<u8>) = (
                            None,
                            None,
                            None,
                            None,
                        );
                        for value in val {
                            match value {
                                FeatureAttributes::Start { value } => {
                                    a = match value {
                                        RangeValue::Exact(v) => Some(*v),
                                        RangeValue::LessThan(v) => Some(*v),
                                        RangeValue::GreaterThan(v) => Some(*v),
                                    };
                                }
                                FeatureAttributes::Stop { value } => {
                                    b = match value {
                                        RangeValue::Exact(v) => Some(*v),
                                        RangeValue::LessThan(v) => Some(*v),
                                        RangeValue::GreaterThan(v) => Some(*v),
                                    };
                                }
                                FeatureAttributes::Strand { value } => {
                                    c = match value {
                                        value => Some(*value),
                                    };
                                }
                                FeatureAttributes::CodonStart { value } => {
                                    d = match value {
                                        value => Some(value.clone()),
                                    };
                                }
                                _ => {}
                            }
                        }
                        let sta = a
                            .map(|o| o as usize)
                            .ok_or(
                                ::anyhow::__private::must_use({
                                    let error = ::anyhow::__private::format_err(
                                        format_args!("No value for start"),
                                    );
                                    error
                                }),
                            )?;
                        let sto = b
                            .map(|t| t as usize)
                            .ok_or(
                                ::anyhow::__private::must_use({
                                    let error = ::anyhow::__private::format_err(
                                        format_args!("No value for stop"),
                                    );
                                    error
                                }),
                            )? - 1;
                        let stra = c
                            .map(|u| u as i8)
                            .ok_or(
                                ::anyhow::__private::must_use({
                                    let error = ::anyhow::__private::format_err(
                                        format_args!("No value for strand"),
                                    );
                                    error
                                }),
                            )?;
                        let cod = d
                            .map(|v| v as usize - 1)
                            .ok_or(
                                ::anyhow::__private::must_use({
                                    let error = ::anyhow::__private::format_err(
                                        format_args!("No value for strand"),
                                    );
                                    error
                                }),
                            )?;
                        let star = sta.try_into()?;
                        let stow = sto.try_into()?;
                        let codd = cod.try_into()?;
                        let mut sliced_sequence: &str = "";
                        if stra == -1 {
                            if cod > 1 {
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "reverse strand coding start more than one {0:?}\n",
                                            &iterablecount,
                                        ),
                                    );
                                };
                                if sto + 1 <= record.sequence.len() {
                                    sliced_sequence = &record.sequence[sta + cod..sto + 1];
                                } else {
                                    sliced_sequence = &record.sequence[sta + cod..sto];
                                }
                            } else {
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "record sta {0:?} sto {1:?} cod {2:?} stra {3:?} record.seq length {4:?}\n",
                                            &sta,
                                            &sto,
                                            &cod,
                                            &stra,
                                            &record.sequence.len(),
                                        ),
                                    );
                                };
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "sliced sta {0:?} sliced sto {1:?} record.id {2:?}\n",
                                            sta,
                                            sto,
                                            &record.id,
                                        ),
                                    );
                                };
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "iterable count is {0:?} reverse strand codon start one\n",
                                            &iterablecount,
                                        ),
                                    );
                                };
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "this is the sequence len {0:?}\n",
                                            &record.sequence.len(),
                                        ),
                                    );
                                };
                                if sto + 1 <= record.sequence.len() {
                                    sliced_sequence = &record.sequence[sta..sto + 1];
                                } else {
                                    sliced_sequence = &record.sequence[sta..sto];
                                }
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "iterable count after is {0:?}\n",
                                            &iterablecount,
                                        ),
                                    );
                                };
                            }
                            let cds_char = sliced_sequence;
                            let prot_seq = translate(&revcomp(cds_char.as_bytes()));
                            let parts: Vec<&str> = prot_seq.split('*').collect();
                            {
                                ::std::io::_print(
                                    format_args!("this is the prot_seq {0:?}\n", &prot_seq),
                                );
                            };
                            record
                                .seq_features
                                .set_counter(key.to_string())
                                .set_start(RangeValue::Exact(star))
                                .set_stop(RangeValue::Exact(stow))
                                .set_sequence_ffn(cds_char.to_string())
                                .set_sequence_faa(parts[0].to_string())
                                .set_codon_start(codd)
                                .set_strand(stra);
                        } else {
                            if cod > 1 {
                                sliced_sequence = &record.sequence[sta + cod - 1..sto];
                            } else {
                                sliced_sequence = &record.sequence[sta - 1..sto];
                            }
                            let cds_char = sliced_sequence;
                            let prot_seq = translate(cds_char.as_bytes());
                            let parts: Vec<&str> = prot_seq.split('*').collect();
                            record
                                .seq_features
                                .set_counter(key.to_string())
                                .set_start(RangeValue::Exact(star))
                                .set_stop(RangeValue::Exact(stow))
                                .set_sequence_ffn(cds_char.to_string())
                                .set_sequence_faa(parts[0].to_string())
                                .set_codon_start(codd)
                                .set_strand(stra);
                        }
                    }
                    return Ok(record.to_owned());
                }
                self.line_buffer.clear();
                self.reader.read_line(&mut self.line_buffer)?;
            }
            Ok(record.to_owned())
        }
    }
    ///stores a value for start or stop (end) which can be denoted as a < value or > value.
    pub enum RangeValue {
        Exact(u32),
        LessThan(u32),
        GreaterThan(u32),
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for RangeValue {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match self {
                RangeValue::Exact(__self_0) => {
                    ::core::fmt::Formatter::debug_tuple_field1_finish(
                        f,
                        "Exact",
                        &__self_0,
                    )
                }
                RangeValue::LessThan(__self_0) => {
                    ::core::fmt::Formatter::debug_tuple_field1_finish(
                        f,
                        "LessThan",
                        &__self_0,
                    )
                }
                RangeValue::GreaterThan(__self_0) => {
                    ::core::fmt::Formatter::debug_tuple_field1_finish(
                        f,
                        "GreaterThan",
                        &__self_0,
                    )
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::hash::Hash for RangeValue {
        #[inline]
        fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            ::core::hash::Hash::hash(&__self_discr, state);
            match self {
                RangeValue::Exact(__self_0) => ::core::hash::Hash::hash(__self_0, state),
                RangeValue::LessThan(__self_0) => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                RangeValue::GreaterThan(__self_0) => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::marker::StructuralPartialEq for RangeValue {}
    #[automatically_derived]
    impl ::core::cmp::PartialEq for RangeValue {
        #[inline]
        fn eq(&self, other: &RangeValue) -> bool {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            let __arg1_discr = ::core::intrinsics::discriminant_value(other);
            __self_discr == __arg1_discr
                && match (self, other) {
                    (RangeValue::Exact(__self_0), RangeValue::Exact(__arg1_0)) => {
                        __self_0 == __arg1_0
                    }
                    (RangeValue::LessThan(__self_0), RangeValue::LessThan(__arg1_0)) => {
                        __self_0 == __arg1_0
                    }
                    (
                        RangeValue::GreaterThan(__self_0),
                        RangeValue::GreaterThan(__arg1_0),
                    ) => __self_0 == __arg1_0,
                    _ => unsafe { ::core::intrinsics::unreachable() }
                }
        }
    }
    #[automatically_derived]
    impl ::core::cmp::Eq for RangeValue {
        #[inline]
        #[doc(hidden)]
        #[coverage(off)]
        fn assert_receiver_is_total_eq(&self) -> () {
            let _: ::core::cmp::AssertParamIsEq<u32>;
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for RangeValue {
        #[inline]
        fn clone(&self) -> RangeValue {
            match self {
                RangeValue::Exact(__self_0) => {
                    RangeValue::Exact(::core::clone::Clone::clone(__self_0))
                }
                RangeValue::LessThan(__self_0) => {
                    RangeValue::LessThan(::core::clone::Clone::clone(__self_0))
                }
                RangeValue::GreaterThan(__self_0) => {
                    RangeValue::GreaterThan(::core::clone::Clone::clone(__self_0))
                }
            }
        }
    }
    impl RangeValue {
        pub fn get_value(&self) -> u32 {
            match self {
                RangeValue::Exact(value) => *value,
                RangeValue::LessThan(value) => *value,
                RangeValue::GreaterThan(value) => *value,
            }
        }
    }
    ///stores the details of the source features in genbank (contigs)
    pub enum SourceAttributes {
        Start { value: RangeValue },
        Stop { value: RangeValue },
        Organism { value: String },
        MolType { value: String },
        Strain { value: String },
        CultureCollection { value: String },
        TypeMaterial { value: String },
        DbXref { value: String },
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for SourceAttributes {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match self {
                SourceAttributes::Start { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Start",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::Stop { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Stop",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::Organism { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Organism",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::MolType { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "MolType",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::Strain { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Strain",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::CultureCollection { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "CultureCollection",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::TypeMaterial { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "TypeMaterial",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::DbXref { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "DbXref",
                        "value",
                        &__self_0,
                    )
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::cmp::Eq for SourceAttributes {
        #[inline]
        #[doc(hidden)]
        #[coverage(off)]
        fn assert_receiver_is_total_eq(&self) -> () {
            let _: ::core::cmp::AssertParamIsEq<RangeValue>;
            let _: ::core::cmp::AssertParamIsEq<String>;
        }
    }
    #[automatically_derived]
    impl ::core::marker::StructuralPartialEq for SourceAttributes {}
    #[automatically_derived]
    impl ::core::cmp::PartialEq for SourceAttributes {
        #[inline]
        fn eq(&self, other: &SourceAttributes) -> bool {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            let __arg1_discr = ::core::intrinsics::discriminant_value(other);
            __self_discr == __arg1_discr
                && match (self, other) {
                    (
                        SourceAttributes::Start { value: __self_0 },
                        SourceAttributes::Start { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::Stop { value: __self_0 },
                        SourceAttributes::Stop { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::Organism { value: __self_0 },
                        SourceAttributes::Organism { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::MolType { value: __self_0 },
                        SourceAttributes::MolType { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::Strain { value: __self_0 },
                        SourceAttributes::Strain { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::CultureCollection { value: __self_0 },
                        SourceAttributes::CultureCollection { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::TypeMaterial { value: __self_0 },
                        SourceAttributes::TypeMaterial { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::DbXref { value: __self_0 },
                        SourceAttributes::DbXref { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    _ => unsafe { ::core::intrinsics::unreachable() }
                }
        }
    }
    #[automatically_derived]
    impl ::core::hash::Hash for SourceAttributes {
        #[inline]
        fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            ::core::hash::Hash::hash(&__self_discr, state);
            match self {
                SourceAttributes::Start { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::Stop { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::Organism { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::MolType { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::Strain { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::CultureCollection { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::TypeMaterial { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::DbXref { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for SourceAttributes {
        #[inline]
        fn clone(&self) -> SourceAttributes {
            match self {
                SourceAttributes::Start { value: __self_0 } => {
                    SourceAttributes::Start {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::Stop { value: __self_0 } => {
                    SourceAttributes::Stop {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::Organism { value: __self_0 } => {
                    SourceAttributes::Organism {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::MolType { value: __self_0 } => {
                    SourceAttributes::MolType {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::Strain { value: __self_0 } => {
                    SourceAttributes::Strain {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::CultureCollection { value: __self_0 } => {
                    SourceAttributes::CultureCollection {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::TypeMaterial { value: __self_0 } => {
                    SourceAttributes::TypeMaterial {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::DbXref { value: __self_0 } => {
                    SourceAttributes::DbXref {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
            }
        }
    }
    impl SourceAttributeBuilder {
        pub fn get_start(&self, key: &str) -> Option<&RangeValue> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::Start { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_stop(&self, key: &str) -> Option<&RangeValue> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::Stop { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_organism(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::Organism { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_mol_type(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::MolType { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_strain(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::Strain { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_type_material(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::TypeMaterial { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_db_xref(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::DbXref { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
    }
    ///builder for the source information on a per record basis
    pub struct SourceAttributeBuilder {
        pub source_attributes: BTreeMap<String, HashSet<SourceAttributes>>,
        pub source_name: Option<String>,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for SourceAttributeBuilder {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "SourceAttributeBuilder",
                "source_attributes",
                &self.source_attributes,
                "source_name",
                &&self.source_name,
            )
        }
    }
    #[automatically_derived]
    impl ::core::default::Default for SourceAttributeBuilder {
        #[inline]
        fn default() -> SourceAttributeBuilder {
            SourceAttributeBuilder {
                source_attributes: ::core::default::Default::default(),
                source_name: ::core::default::Default::default(),
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for SourceAttributeBuilder {
        #[inline]
        fn clone(&self) -> SourceAttributeBuilder {
            SourceAttributeBuilder {
                source_attributes: ::core::clone::Clone::clone(&self.source_attributes),
                source_name: ::core::clone::Clone::clone(&self.source_name),
            }
        }
    }
    impl SourceAttributeBuilder {
        pub fn set_source_name(&mut self, name: String) {
            self.source_name = Some(name);
        }
        pub fn get_source_name(&self) -> Option<&String> {
            self.source_name.as_ref()
        }
        pub fn add_source_attribute(
            &mut self,
            key: String,
            attribute: SourceAttributes,
        ) {
            self.source_attributes
                .entry(key)
                .or_insert_with(HashSet::new)
                .insert(attribute);
        }
        pub fn get_source_attributes(
            &self,
            key: &str,
        ) -> Option<&HashSet<SourceAttributes>> {
            self.source_attributes.get(key)
        }
    }
    impl SourceAttributeBuilder {
        pub fn new() -> Self {
            SourceAttributeBuilder {
                source_attributes: BTreeMap::new(),
                source_name: None,
            }
        }
        pub fn set_counter(&mut self, counter: String) -> &mut Self {
            self.source_name = Some(counter);
            self
        }
        pub fn insert_to(&mut self, value: SourceAttributes) {
            if let Some(counter) = &self.source_name {
                self.source_attributes
                    .entry(counter.to_string())
                    .or_insert_with(HashSet::new)
                    .insert(value);
            } else {
                {
                    ::core::panicking::panic_fmt(format_args!("Counter key not set"));
                };
            }
        }
        pub fn set_start(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(SourceAttributes::Start { value });
            self
        }
        pub fn set_stop(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(SourceAttributes::Stop { value });
            self
        }
        pub fn set_organism(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::Organism {
                value,
            });
            self
        }
        pub fn set_mol_type(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::MolType { value });
            self
        }
        pub fn set_strain(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::Strain { value });
            self
        }
        pub fn set_type_material(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::TypeMaterial {
                value,
            });
            self
        }
        pub fn set_db_xref(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::DbXref { value });
            self
        }
        pub fn build(self) -> BTreeMap<String, HashSet<SourceAttributes>> {
            self.source_attributes
        }
        pub fn iter_sorted(
            &self,
        ) -> std::collections::btree_map::Iter<String, HashSet<SourceAttributes>> {
            self.source_attributes.iter()
        }
        pub fn default() -> Self {
            SourceAttributeBuilder {
                source_attributes: BTreeMap::new(),
                source_name: None,
            }
        }
    }
    ///attributes for each feature, cds or gene
    pub enum FeatureAttributes {
        Start { value: RangeValue },
        Stop { value: RangeValue },
        Gene { value: String },
        Product { value: String },
        CodonStart { value: u8 },
        Strand { value: i8 },
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for FeatureAttributes {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match self {
                FeatureAttributes::Start { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Start",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::Stop { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Stop",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::Gene { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Gene",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::Product { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Product",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::CodonStart { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "CodonStart",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::Strand { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Strand",
                        "value",
                        &__self_0,
                    )
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::cmp::Eq for FeatureAttributes {
        #[inline]
        #[doc(hidden)]
        #[coverage(off)]
        fn assert_receiver_is_total_eq(&self) -> () {
            let _: ::core::cmp::AssertParamIsEq<RangeValue>;
            let _: ::core::cmp::AssertParamIsEq<String>;
            let _: ::core::cmp::AssertParamIsEq<u8>;
            let _: ::core::cmp::AssertParamIsEq<i8>;
        }
    }
    #[automatically_derived]
    impl ::core::hash::Hash for FeatureAttributes {
        #[inline]
        fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            ::core::hash::Hash::hash(&__self_discr, state);
            match self {
                FeatureAttributes::Start { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::Stop { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::Gene { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::Product { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::CodonStart { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::Strand { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::marker::StructuralPartialEq for FeatureAttributes {}
    #[automatically_derived]
    impl ::core::cmp::PartialEq for FeatureAttributes {
        #[inline]
        fn eq(&self, other: &FeatureAttributes) -> bool {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            let __arg1_discr = ::core::intrinsics::discriminant_value(other);
            __self_discr == __arg1_discr
                && match (self, other) {
                    (
                        FeatureAttributes::Start { value: __self_0 },
                        FeatureAttributes::Start { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::Stop { value: __self_0 },
                        FeatureAttributes::Stop { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::Gene { value: __self_0 },
                        FeatureAttributes::Gene { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::Product { value: __self_0 },
                        FeatureAttributes::Product { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::CodonStart { value: __self_0 },
                        FeatureAttributes::CodonStart { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::Strand { value: __self_0 },
                        FeatureAttributes::Strand { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    _ => unsafe { ::core::intrinsics::unreachable() }
                }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for FeatureAttributes {
        #[inline]
        fn clone(&self) -> FeatureAttributes {
            match self {
                FeatureAttributes::Start { value: __self_0 } => {
                    FeatureAttributes::Start {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::Stop { value: __self_0 } => {
                    FeatureAttributes::Stop {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::Gene { value: __self_0 } => {
                    FeatureAttributes::Gene {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::Product { value: __self_0 } => {
                    FeatureAttributes::Product {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::CodonStart { value: __self_0 } => {
                    FeatureAttributes::CodonStart {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::Strand { value: __self_0 } => {
                    FeatureAttributes::Strand {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
            }
        }
    }
    impl FeatureAttributeBuilder {
        pub fn get_start(&self, key: &str) -> Option<&RangeValue> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Start { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_stop(&self, key: &str) -> Option<&RangeValue> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Stop { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_gene(&self, key: &str) -> Option<&String> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Gene { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_product(&self, key: &str) -> Option<&String> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Product { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_codon_start(&self, key: &str) -> Option<&u8> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::CodonStart { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_strand(&self, key: &str) -> Option<&i8> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Strand { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
    }
    ///builder for the feature information on a per coding sequence (CDS) basis
    pub struct FeatureAttributeBuilder {
        pub attributes: BTreeMap<String, HashSet<FeatureAttributes>>,
        locus_tag: Option<String>,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for FeatureAttributeBuilder {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "FeatureAttributeBuilder",
                "attributes",
                &self.attributes,
                "locus_tag",
                &&self.locus_tag,
            )
        }
    }
    #[automatically_derived]
    impl ::core::default::Default for FeatureAttributeBuilder {
        #[inline]
        fn default() -> FeatureAttributeBuilder {
            FeatureAttributeBuilder {
                attributes: ::core::default::Default::default(),
                locus_tag: ::core::default::Default::default(),
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for FeatureAttributeBuilder {
        #[inline]
        fn clone(&self) -> FeatureAttributeBuilder {
            FeatureAttributeBuilder {
                attributes: ::core::clone::Clone::clone(&self.attributes),
                locus_tag: ::core::clone::Clone::clone(&self.locus_tag),
            }
        }
    }
    impl FeatureAttributeBuilder {
        pub fn new() -> Self {
            FeatureAttributeBuilder {
                attributes: BTreeMap::new(),
                locus_tag: None,
            }
        }
        pub fn set_counter(&mut self, counter: String) -> &mut Self {
            self.locus_tag = Some(counter);
            self
        }
        pub fn insert_to(&mut self, value: FeatureAttributes) {
            if let Some(counter) = &self.locus_tag {
                self.attributes
                    .entry(counter.to_string())
                    .or_insert_with(HashSet::new)
                    .insert(value);
            } else {
                {
                    ::core::panicking::panic_fmt(format_args!("Counter key not set"));
                };
            }
        }
        pub fn set_start(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(FeatureAttributes::Start { value });
            self
        }
        pub fn set_stop(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(FeatureAttributes::Stop { value });
            self
        }
        pub fn set_gene(&mut self, value: String) -> &mut Self {
            self.insert_to(FeatureAttributes::Gene { value });
            self
        }
        pub fn set_product(&mut self, value: String) -> &mut Self {
            self.insert_to(FeatureAttributes::Product {
                value,
            });
            self
        }
        pub fn set_codon_start(&mut self, value: u8) -> &mut Self {
            self.insert_to(FeatureAttributes::CodonStart {
                value,
            });
            self
        }
        pub fn set_strand(&mut self, value: i8) -> &mut Self {
            self.insert_to(FeatureAttributes::Strand { value });
            self
        }
        pub fn build(self) -> BTreeMap<String, HashSet<FeatureAttributes>> {
            self.attributes
        }
        pub fn iter_sorted(
            &self,
        ) -> std::collections::btree_map::Iter<String, HashSet<FeatureAttributes>> {
            self.attributes.iter()
        }
        pub fn default() -> Self {
            FeatureAttributeBuilder {
                attributes: BTreeMap::new(),
                locus_tag: None,
            }
        }
    }
    ///stores the sequences of the coding sequences (genes) and proteins. Also stores start, stop, codon_start and strand information
    pub enum SequenceAttributes {
        Start { value: RangeValue },
        Stop { value: RangeValue },
        SequenceFfn { value: String },
        SequenceFaa { value: String },
        CodonStart { value: u8 },
        Strand { value: i8 },
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for SequenceAttributes {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match self {
                SequenceAttributes::Start { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Start",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::Stop { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Stop",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::SequenceFfn { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "SequenceFfn",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::SequenceFaa { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "SequenceFaa",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::CodonStart { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "CodonStart",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::Strand { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Strand",
                        "value",
                        &__self_0,
                    )
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::cmp::Eq for SequenceAttributes {
        #[inline]
        #[doc(hidden)]
        #[coverage(off)]
        fn assert_receiver_is_total_eq(&self) -> () {
            let _: ::core::cmp::AssertParamIsEq<RangeValue>;
            let _: ::core::cmp::AssertParamIsEq<String>;
            let _: ::core::cmp::AssertParamIsEq<u8>;
            let _: ::core::cmp::AssertParamIsEq<i8>;
        }
    }
    #[automatically_derived]
    impl ::core::marker::StructuralPartialEq for SequenceAttributes {}
    #[automatically_derived]
    impl ::core::cmp::PartialEq for SequenceAttributes {
        #[inline]
        fn eq(&self, other: &SequenceAttributes) -> bool {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            let __arg1_discr = ::core::intrinsics::discriminant_value(other);
            __self_discr == __arg1_discr
                && match (self, other) {
                    (
                        SequenceAttributes::Start { value: __self_0 },
                        SequenceAttributes::Start { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::Stop { value: __self_0 },
                        SequenceAttributes::Stop { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::SequenceFfn { value: __self_0 },
                        SequenceAttributes::SequenceFfn { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::SequenceFaa { value: __self_0 },
                        SequenceAttributes::SequenceFaa { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::CodonStart { value: __self_0 },
                        SequenceAttributes::CodonStart { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::Strand { value: __self_0 },
                        SequenceAttributes::Strand { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    _ => unsafe { ::core::intrinsics::unreachable() }
                }
        }
    }
    #[automatically_derived]
    impl ::core::hash::Hash for SequenceAttributes {
        #[inline]
        fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            ::core::hash::Hash::hash(&__self_discr, state);
            match self {
                SequenceAttributes::Start { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::Stop { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::SequenceFfn { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::SequenceFaa { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::CodonStart { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::Strand { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for SequenceAttributes {
        #[inline]
        fn clone(&self) -> SequenceAttributes {
            match self {
                SequenceAttributes::Start { value: __self_0 } => {
                    SequenceAttributes::Start {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::Stop { value: __self_0 } => {
                    SequenceAttributes::Stop {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::SequenceFfn { value: __self_0 } => {
                    SequenceAttributes::SequenceFfn {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::SequenceFaa { value: __self_0 } => {
                    SequenceAttributes::SequenceFaa {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::CodonStart { value: __self_0 } => {
                    SequenceAttributes::CodonStart {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::Strand { value: __self_0 } => {
                    SequenceAttributes::Strand {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
            }
        }
    }
    impl SequenceAttributeBuilder {
        pub fn get_start(&self, key: &str) -> Option<&RangeValue> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::Start { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_stop(&self, key: &str) -> Option<&RangeValue> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::Stop { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_sequence_ffn(&self, key: &str) -> Option<&String> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::SequenceFfn { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_sequence_faa(&self, key: &str) -> Option<&String> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::SequenceFaa { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_codon_start(&self, key: &str) -> Option<&u8> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::CodonStart { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_strand(&self, key: &str) -> Option<&i8> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::Strand { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
    }
    ///builder for the sequence information on a per coding sequence (CDS) basis
    pub struct SequenceAttributeBuilder {
        pub seq_attributes: BTreeMap<String, HashSet<SequenceAttributes>>,
        locus_tag: Option<String>,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for SequenceAttributeBuilder {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "SequenceAttributeBuilder",
                "seq_attributes",
                &self.seq_attributes,
                "locus_tag",
                &&self.locus_tag,
            )
        }
    }
    #[automatically_derived]
    impl ::core::default::Default for SequenceAttributeBuilder {
        #[inline]
        fn default() -> SequenceAttributeBuilder {
            SequenceAttributeBuilder {
                seq_attributes: ::core::default::Default::default(),
                locus_tag: ::core::default::Default::default(),
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for SequenceAttributeBuilder {
        #[inline]
        fn clone(&self) -> SequenceAttributeBuilder {
            SequenceAttributeBuilder {
                seq_attributes: ::core::clone::Clone::clone(&self.seq_attributes),
                locus_tag: ::core::clone::Clone::clone(&self.locus_tag),
            }
        }
    }
    impl SequenceAttributeBuilder {
        pub fn new() -> Self {
            SequenceAttributeBuilder {
                seq_attributes: BTreeMap::new(),
                locus_tag: None,
            }
        }
        pub fn set_counter(&mut self, counter: String) -> &mut Self {
            self.locus_tag = Some(counter);
            self
        }
        pub fn insert_to(&mut self, value: SequenceAttributes) {
            if let Some(counter) = &self.locus_tag {
                self.seq_attributes
                    .entry(counter.to_string())
                    .or_insert_with(HashSet::new)
                    .insert(value);
            } else {
                {
                    ::core::panicking::panic_fmt(format_args!("Counter key not set"));
                };
            }
        }
        pub fn set_start(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(SequenceAttributes::Start { value });
            self
        }
        pub fn set_stop(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(SequenceAttributes::Stop { value });
            self
        }
        pub fn set_sequence_ffn(&mut self, value: String) -> &mut Self {
            self.insert_to(SequenceAttributes::SequenceFfn {
                value,
            });
            self
        }
        pub fn set_sequence_faa(&mut self, value: String) -> &mut Self {
            self.insert_to(SequenceAttributes::SequenceFaa {
                value,
            });
            self
        }
        pub fn set_codon_start(&mut self, value: u8) -> &mut Self {
            self.insert_to(SequenceAttributes::CodonStart {
                value,
            });
            self
        }
        pub fn set_strand(&mut self, value: i8) -> &mut Self {
            self.insert_to(SequenceAttributes::Strand {
                value,
            });
            self
        }
        pub fn build(self) -> BTreeMap<String, HashSet<SequenceAttributes>> {
            self.seq_attributes
        }
        pub fn iter_sorted(
            &self,
        ) -> std::collections::btree_map::Iter<String, HashSet<SequenceAttributes>> {
            self.seq_attributes.iter()
        }
        pub fn default() -> Self {
            SequenceAttributeBuilder {
                seq_attributes: BTreeMap::new(),
                locus_tag: None,
            }
        }
    }
    ///product lines can contain difficult to parse punctuation such as biochemical symbols like unclosed single quotes, superscripts, single and double brackets etc.
    ///here we substitute these for an underscore
    pub fn substitute_odd_punctuation(input: String) -> Result<String, anyhow::Error> {
        let re = Regex::new(r"[/?()',`]|[α-ωΑ-Ω]")?;
        let cleaned = input.trim_end_matches(&['\r', '\n'][..]);
        Ok(re.replace_all(cleaned, "_").to_string())
    }
    ///GFF3 field9 construct
    pub struct GFFInner {
        id: String,
        name: String,
        locus_tag: String,
        gene: String,
        product: String,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for GFFInner {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field5_finish(
                f,
                "GFFInner",
                "id",
                &self.id,
                "name",
                &self.name,
                "locus_tag",
                &self.locus_tag,
                "gene",
                &self.gene,
                "product",
                &&self.product,
            )
        }
    }
    impl GFFInner {
        pub fn new(
            id: String,
            name: String,
            locus_tag: String,
            gene: String,
            product: String,
        ) -> Self {
            GFFInner {
                id,
                name,
                locus_tag,
                gene,
                product,
            }
        }
    }
    ///The main GFF3 construct
    pub struct GFFOuter<'a> {
        seqid: String,
        source: String,
        type_val: String,
        start: u32,
        end: u32,
        score: f64,
        strand: String,
        phase: u8,
        attributes: &'a GFFInner,
    }
    #[automatically_derived]
    impl<'a> ::core::fmt::Debug for GFFOuter<'a> {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            let names: &'static _ = &[
                "seqid",
                "source",
                "type_val",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ];
            let values: &[&dyn ::core::fmt::Debug] = &[
                &self.seqid,
                &self.source,
                &self.type_val,
                &self.start,
                &self.end,
                &self.score,
                &self.strand,
                &self.phase,
                &&self.attributes,
            ];
            ::core::fmt::Formatter::debug_struct_fields_finish(
                f,
                "GFFOuter",
                names,
                values,
            )
        }
    }
    impl<'a> GFFOuter<'a> {
        pub fn new(
            seqid: String,
            source: String,
            type_val: String,
            start: u32,
            end: u32,
            score: f64,
            strand: String,
            phase: u8,
            attributes: &'a GFFInner,
        ) -> Self {
            GFFOuter {
                seqid,
                source,
                type_val,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            }
        }
        pub fn field9_attributes_build(&self) -> String {
            let mut full_field9 = Vec::new();
            if !self.attributes.id.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("id={0}", self.attributes.id),
                            );
                            res
                        }),
                    );
            }
            if !self.attributes.name.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("name={0}", self.attributes.name),
                            );
                            res
                        }),
                    );
            }
            if !self.attributes.gene.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("gene={0}", self.attributes.gene),
                            );
                            res
                        }),
                    );
            }
            if !self.attributes.locus_tag.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("locus_tag={0}", self.attributes.locus_tag),
                            );
                            res
                        }),
                    );
            }
            if !self.attributes.product.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("product={0}", self.attributes.product),
                            );
                            res
                        }),
                    );
            }
            full_field9.join(";")
        }
    }
    ///formats the translation string which can be mulitple lines, for embl
    pub fn format_translation(translation: &str) -> String {
        let mut formatted = String::new();
        let cleaned_translation = translation.replace("\n", "");
        formatted.push_str("                     /translation=\"");
        let line_length: usize = 60;
        let final_num = line_length - 15;
        formatted
            .push_str(
                &::alloc::__export::must_use({
                    let res = ::alloc::fmt::format(
                        format_args!("{0}\n", &cleaned_translation[0..final_num]),
                    );
                    res
                }),
            );
        for i in (47..translation.len()).step_by(60) {
            let end = i + 60 - 1;
            let valid_end = if end >= translation.len() {
                &cleaned_translation.len() - 1
            } else {
                end
            };
            formatted
                .push_str(
                    &::alloc::__export::must_use({
                        let res = ::alloc::fmt::format(
                            format_args!(
                                "                     {0}",
                                &cleaned_translation[i..valid_end],
                            ),
                        );
                        res
                    }),
                );
            {
                ::std::io::_print(
                    format_args!(
                        "cleaned translation leng is {0:?}\n",
                        &cleaned_translation[i..valid_end].len(),
                    ),
                );
            };
            if *&cleaned_translation[i..valid_end].len() < 59 {
                formatted.push('\"');
            } else {
                formatted.push('\n');
            }
        }
        formatted
    }
    ///writes the DNA sequence in gbk format with numbering
    pub fn write_gbk_format_sequence(sequence: &str, file: &mut File) -> io::Result<()> {
        file.write_fmt(format_args!("ORIGIN\n"))?;
        let mut formatted = String::new();
        let cleaned_input = sequence.replace("\n", "");
        let mut index = 1;
        for (_i, chunk) in cleaned_input.as_bytes().chunks(60).enumerate() {
            formatted
                .push_str(
                    &::alloc::__export::must_use({
                        let res = ::alloc::fmt::format(format_args!("{0:>5} ", index));
                        res
                    }),
                );
            for (j, sub_chunk) in chunk.chunks(10).enumerate() {
                if j > 0 {
                    formatted.push(' ');
                }
                formatted.push_str(&String::from_utf8_lossy(sub_chunk));
            }
            formatted.push('\n');
            index += 60;
        }
        file.write_fmt(format_args!("{0:>6}\n", &formatted))?;
        file.write_fmt(format_args!("//\n"))?;
        Ok(())
    }
    ///saves the parsed data in genbank format
    pub fn gbk_write(
        seq_region: BTreeMap<String, (u32, u32)>,
        record_vec: Vec<Record>,
        filename: &str,
    ) -> io::Result<()> {
        let now = Local::now();
        let formatted_date = now.format("%d-%b-%Y").to_string().to_uppercase();
        let mut file = OpenOptions::new()
            .write(true)
            .append(true)
            .create(true)
            .open(filename)?;
        for (i, (key, _val)) in seq_region.iter().enumerate() {
            let strain = match &record_vec[i].source_map.get_strain(key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let organism = match &record_vec[i].source_map.get_organism(key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let mol_type = match &record_vec[i].source_map.get_mol_type(key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let type_material = match &record_vec[i].source_map.get_type_material(&key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let db_xref = match &record_vec[i].source_map.get_db_xref(key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let source_stop = match &record_vec[i].source_map.get_stop(key) {
                Some(value) => value.get_value(),
                None => {
                    {
                        {
                            ::std::io::_print(format_args!("stop value not found\n"));
                        };
                        None
                    }
                        .expect("stop value not received")
                }
            };
            file.write_fmt(
                format_args!(
                    "LOCUS       {0}             {1} bp    DNA     linear CON {2}\n",
                    &key,
                    &record_vec[i].sequence.len(),
                    &formatted_date,
                ),
            )?;
            file.write_fmt(format_args!("DEFINITION  {0} {1}.\n", &organism, &strain))?;
            file.write_fmt(format_args!("ACCESSION   {0}\n", &key))?;
            file.write_fmt(format_args!("KEYWORDS    .\n"))?;
            file.write_fmt(format_args!("SOURCE      {0} {1}\n", &organism, &strain))?;
            file.write_fmt(format_args!("  ORGANISM  {0} {1}\n", &organism, &strain))?;
            file.write_fmt(format_args!("FEATURES             Location/Qualifiers\n"))?;
            file.write_fmt(format_args!("     source          1..{0}\n", &source_stop))?;
            file.write_fmt(
                format_args!("                     /organism=\"{0}\"\n", &strain),
            )?;
            file.write_fmt(
                format_args!("                     /mol_type=\"{0}\"\n", &mol_type),
            )?;
            file.write_fmt(
                format_args!("                     /strain=\"{0}\"\n", &strain),
            )?;
            if type_material != *"Unknown".to_string() {
                file.write_fmt(
                    format_args!(
                        "                     /type_material=\"{0}\"\n",
                        &type_material,
                    ),
                )?;
            }
            file.write_fmt(
                format_args!("                     /db_xref=\"{0}\"\n", &db_xref),
            )?;
            for (locus_tag, _value) in &record_vec[i].cds.attributes {
                let start = match &record_vec[i].cds.get_start(locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("start value not found\n"));
                            };
                            None
                        }
                            .expect("start value not received")
                    }
                };
                let stop = match &record_vec[i].cds.get_stop(locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("stop value not found\n"));
                            };
                            None
                        }
                            .expect("stop value not received")
                    }
                };
                let product = match &record_vec[i].cds.get_product(locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown product".to_string(),
                };
                let strand = match &record_vec[i].cds.get_strand(locus_tag) {
                    Some(value) => **value,
                    None => 0,
                };
                let codon_start = match &record_vec[i].cds.get_codon_start(locus_tag) {
                    Some(value) => **value,
                    None => 0,
                };
                let gene = match &record_vec[i].cds.get_gene(locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown".to_string(),
                };
                let translation = match &record_vec[i]
                    .seq_features
                    .get_sequence_faa(locus_tag)
                {
                    Some(value) => value.to_string(),
                    None => "unknown".to_string(),
                };
                if strand == 1 {
                    file.write_fmt(
                        format_args!("     gene            {0}..{1}\n", &start, &stop),
                    )?;
                } else {
                    file.write_fmt(
                        format_args!(
                            "     gene            complement({0}..{1})\n",
                            &start,
                            &stop,
                        ),
                    )?;
                }
                file.write_fmt(
                    format_args!("                     /locus_tag=\"{0}\"\n", &locus_tag),
                )?;
                if strand == 1 {
                    file.write_fmt(
                        format_args!("     CDS             {0}..{1}\n", &start, &stop),
                    )?;
                } else {
                    file.write_fmt(
                        format_args!(
                            "     CDS             complement({0}..{1})\n",
                            &start,
                            &stop,
                        ),
                    )?;
                }
                file.write_fmt(
                    format_args!("                     /locus_tag=\"{0}\"\n", &locus_tag),
                )?;
                file.write_fmt(
                    format_args!(
                        "                     /codon_start=\"{0}\"\n",
                        &codon_start,
                    ),
                )?;
                if gene != "unknown" {
                    file.write_fmt(
                        format_args!("                     /gene=\"{0}\"\n", &gene),
                    )?;
                }
                if translation != "unknown" {
                    let formatted_translation = format_translation(&translation);
                    file.write_fmt(format_args!("{0}\n", &formatted_translation))?;
                }
                file.write_fmt(
                    format_args!("                     /product=\"{0}\"\n", &product),
                )?;
            }
            write_gbk_format_sequence(&record_vec[i].sequence, &mut file)?;
        }
        Ok(())
    }
    ///saves the parsed data in gff3 format
    #[allow(unused_assignments)]
    #[allow(unused_variables)]
    pub fn gff_write(
        seq_region: BTreeMap<String, (u32, u32)>,
        mut record_vec: Vec<Record>,
        filename: &str,
        dna: bool,
    ) -> io::Result<()> {
        let mut file = OpenOptions::new().append(true).create(true).open(filename)?;
        if file.metadata()?.len() == 0 {
            file.write_fmt(format_args!("##gff-version 3\n"))?;
        }
        let mut full_seq = String::new();
        let mut prev_end: u32 = 0;
        for (k, v) in seq_region.iter() {
            file.write_fmt(
                format_args!("##sequence-region\t{0}\t{1}\t{2}\n", &k, v.0, v.1),
            )?;
        }
        for ((source_name, (seq_start, seq_end)), record) in seq_region
            .iter()
            .zip(record_vec.drain(..))
        {
            if dna == true {
                full_seq.push_str(&record.sequence);
            }
            for (locus_tag, _valu) in &record.cds.attributes {
                let start = match record.cds.get_start(&locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("start value not found\n"));
                            };
                            None
                        }
                            .expect("start value not received")
                    }
                };
                let stop = match record.cds.get_stop(&locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("stop value not found\n"));
                            };
                            None
                        }
                            .expect("stop value not received")
                    }
                };
                let gene = match record.cds.get_gene(&locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown".to_string(),
                };
                let product = match record.cds.get_product(&locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown product".to_string(),
                };
                let strand = match record.cds.get_strand(&locus_tag) {
                    Some(valu) => {
                        match valu {
                            1 => "+".to_string(),
                            -1 => "-".to_string(),
                            _ => {
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "unexpected strand value {0} for locus_tag {1}\n",
                                            valu,
                                            &locus_tag,
                                        ),
                                    );
                                };
                                "unknownstrand".to_string()
                            }
                        }
                    }
                    None => "unknownvalue".to_string(),
                };
                let phase = match record.cds.get_codon_start(&locus_tag) {
                    Some(valuer) => {
                        match valuer {
                            1 => 0,
                            2 => 1,
                            3 => 2,
                            _ => {
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "unexpected phase value {0} in the bagging area for locus_tag {1}\n",
                                            valuer,
                                            &locus_tag,
                                        ),
                                    );
                                };
                                1
                            }
                        }
                    }
                    None => 1,
                };
                let gff_inner = GFFInner::new(
                    locus_tag.to_string(),
                    source_name.clone(),
                    locus_tag.to_string(),
                    gene,
                    product,
                );
                let gff_outer = GFFOuter::new(
                    source_name.clone(),
                    ".".to_string(),
                    "CDS".to_string(),
                    start + prev_end,
                    stop + prev_end,
                    0.0,
                    strand,
                    phase,
                    &gff_inner,
                );
                let field9_attributes = gff_outer.field9_attributes_build();
                file.write_fmt(
                    format_args!(
                        "{0}\t{1}\t{2}\t{3:?}\t{4:?}\t{5}\t{6}\t{7}\t{8}\n",
                        gff_outer.seqid,
                        gff_outer.source,
                        gff_outer.type_val,
                        gff_outer.start,
                        gff_outer.end,
                        gff_outer.score,
                        gff_outer.strand,
                        gff_outer.phase,
                        field9_attributes,
                    ),
                )?;
            }
            prev_end = *seq_end;
        }
        if dna {
            file.write_fmt(format_args!("##FASTA\n"))?;
            file.write_fmt(format_args!("{0}\n", full_seq))?;
        }
        Ok(())
    }
    ///internal record containing data from a single source or contig.  Has multiple features.
    pub struct Record {
        pub id: String,
        pub length: u32,
        pub sequence: String,
        pub start: usize,
        pub end: usize,
        pub strand: i32,
        pub cds: FeatureAttributeBuilder,
        pub source_map: SourceAttributeBuilder,
        pub seq_features: SequenceAttributeBuilder,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for Record {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            let names: &'static _ = &[
                "id",
                "length",
                "sequence",
                "start",
                "end",
                "strand",
                "cds",
                "source_map",
                "seq_features",
            ];
            let values: &[&dyn ::core::fmt::Debug] = &[
                &self.id,
                &self.length,
                &self.sequence,
                &self.start,
                &self.end,
                &self.strand,
                &self.cds,
                &self.source_map,
                &&self.seq_features,
            ];
            ::core::fmt::Formatter::debug_struct_fields_finish(
                f,
                "Record",
                names,
                values,
            )
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for Record {
        #[inline]
        fn clone(&self) -> Record {
            Record {
                id: ::core::clone::Clone::clone(&self.id),
                length: ::core::clone::Clone::clone(&self.length),
                sequence: ::core::clone::Clone::clone(&self.sequence),
                start: ::core::clone::Clone::clone(&self.start),
                end: ::core::clone::Clone::clone(&self.end),
                strand: ::core::clone::Clone::clone(&self.strand),
                cds: ::core::clone::Clone::clone(&self.cds),
                source_map: ::core::clone::Clone::clone(&self.source_map),
                seq_features: ::core::clone::Clone::clone(&self.seq_features),
            }
        }
    }
    impl Record {
        /// Create a new instance.
        pub fn new() -> Self {
            Record {
                id: "".to_owned(),
                length: 0,
                sequence: "".to_owned(),
                start: 0,
                end: 0,
                strand: 0,
                source_map: SourceAttributeBuilder::new(),
                cds: FeatureAttributeBuilder::new(),
                seq_features: SequenceAttributeBuilder::new(),
            }
        }
        pub fn is_empty(&mut self) -> bool {
            self.id.is_empty() && self.length == 0
        }
        pub fn check(&mut self) -> Result<(), &str> {
            if self.id().is_empty() {
                return Err("Expecting id for Embl record.");
            }
            Ok(())
        }
        pub fn id(&mut self) -> &str {
            &self.id
        }
        pub fn length(&mut self) -> u32 {
            self.length
        }
        pub fn sequence(&mut self) -> &str {
            &self.sequence
        }
        pub fn start(&mut self) -> u32 {
            self.start.try_into().unwrap()
        }
        pub fn end(&mut self) -> u32 {
            self.end.try_into().unwrap()
        }
        pub fn strand(&mut self) -> i32 {
            self.strand
        }
        pub fn cds(&mut self) -> FeatureAttributeBuilder {
            self.cds.clone()
        }
        pub fn source_map(&mut self) -> SourceAttributeBuilder {
            self.source_map.clone()
        }
        pub fn seq_features(&mut self) -> SequenceAttributeBuilder {
            self.seq_features.clone()
        }
        fn rec_clear(&mut self) {
            self.id.clear();
            self.length = 0;
            self.sequence.clear();
            self.start = 0;
            self.end = 0;
            self.strand = 0;
            self.source_map = SourceAttributeBuilder::new();
            self.cds = FeatureAttributeBuilder::new();
            self.seq_features = SequenceAttributeBuilder::new();
        }
    }
    impl Default for Record {
        fn default() -> Self {
            Self::new()
        }
    }
    #[allow(dead_code)]
    pub struct Config {
        filename: String,
    }
    impl Config {
        pub fn new(args: &[String]) -> Result<Config, &str> {
            if args.len() < 2 {
                {
                    ::core::panicking::panic_fmt(
                        format_args!("not enough arguments, please provide filename"),
                    );
                };
            }
            let filename = args[1].clone();
            Ok(Config { filename })
        }
    }
}
pub mod gbk {
    //! # A Genbank to GFF parser
    //!
    //!
    //! You are able to parse genbank and save as a GFF (gff3) format as well as extracting DNA sequences, gene DNA sequences (ffn) and protein fasta sequences (faa)
    //!
    //! You can also create new records and save as a genbank (gbk) format
    //!
    //! ## Detailed Explanation
    //!
    //!
    //! The Genbank parser contains:
    //!
    //! Records - a top level structure which consists of either one record (single genbank) or multiple instances of record (multi-genbank).
    //!
    //! Each Record contains:
    //!
    //! 1. A source, ```SourceAttributes```, construct(enum) of counter (source name), start, stop [of source or contig], organism, mol_type, strain, type_material, db_xref
    //! 2. Features, ```FeatureAttributes```, construct(enum) of counter (locus tag), gene (if present), product, codon start, strand, start, stop [of cds/gene]
    //! 3. Sequence features, ```SequenceAttributes```, construct(enum) of counter (locus tag), sequence_ffn (DNA gene sequence) sequence_faa (protein translation), strand, codon start, start, stop [cds/gene]
    //! 4. The DNA sequence of the whole record (or contig)
    //!
    //!  Example to extract and print all the protein sequence fasta, example using getters or get_ functionality
    //!
    //!
    //!```rust
    //! use clap::Parser;
    //! use std::fs::File;
    //! use microBioRust::gbk::Reader;
    //! use std::io;
    //!
    //! #[derive(Parser, Debug)]
    //! #[clap(author, version, about)]
    //! struct Arguments {
    //! #[clap(short, long)]
    //! filename: String,
    //! }
    //!
    //! pub fn genbank_to_faa() -> Result<(), anyhow::Error> {
    //!            let args = Arguments::parse();
    //!            let file_gbk = File::open(args.filename)?;
    //!            let mut reader = Reader::new(file_gbk);
    //!            let mut records = reader.records();
    //!            loop {
    //!                //collect from each record advancing on a next record basis, count cds records
    //!                match records.next() {	
    //!                    Some(Ok(mut record)) => {
    //!		                     for (k, v) in &record.cds.attributes {
    //!		                         match record.seq_features.get_sequence_faa(&k) {
    //!		                             Some(value) => { let seq_faa = value.to_string();
    //!				                              println!(">{}|{}\n{}", &record.id, &k, seq_faa);
    //!						              },
    //!				             _ => (),
    //!				             };
    //!		                         }
    //!                                      },
    //!	               Some(Err(e)) => { println!("Error encountered - an err {:?}", e); },
    //!	               None => break,
    //!	               }
    //!                 }
    //!            return Ok(());
    //!  }
    //!```
    //!
    //!  Example to extract the protein sequences with simplified genbank! macro use
    //!
    //!```rust
    //! use clap::Parser;
    //! use std::fs::File;
    //! use microBioRust::gbk::Reader;
    //! use std::io;
    //! use microBioRust::genbank;
    //!
    //!
    //! #[derive(Parser, Debug)]
    //! #[clap(author, version, about)]
    //! struct Arguments {
    //! #[clap(short, long)]
    //! filename: String,
    //! }
    //!
    //! pub fn genbank_to_faa() -> Result<(), anyhow::Error> {
    //!            let args = Arguments::parse();
    //!            let records = genbank!(&args.filename);
    //!            for record in records {
    //!	          for (k, v) in &record.cds.attributes {
    //!                  if let Some(seq) = record.seq_features.get_sequence_faa(k) {
    //!		        println!(">{}|{}\n{}", &record.id, &k, seq);
    //!                     }
    //!                  }
    //!            }
    //!            return Ok(());
    //!  }
    //!
    //!```
    //!  Example to save a provided multi- or single genbank file as a GFF file (by joining any multi-genbank)
    //!
    //! ```rust
    //!    use microBioRust::gbk::gff_write;
    //!    use microBioRust::gbk::Reader;
    //!    use microBioRust::gbk::Record;
    //!    use std::collections::BTreeMap;
    //!    use std::fs::File;
    //!    use clap::Parser;
    //!    use std::io;
    //!
    //!   #[derive(Parser, Debug)]
    //!   #[clap(author, version, about)]
    //!   struct Arguments {
    //!   #[clap(short, long)]
    //!   filename: String,
    //!   }
    //!
    //!    pub fn genbank_to_gff() -> io::Result<()> {
    //!        let args = Arguments::parse();
    //!        let file_gbk = File::open(&args.filename)?;
    //!        let prev_start: u32 = 0;
    //!        let mut prev_end: u32 = 0;
    //!        let mut reader = Reader::new(file_gbk);
    //!        let mut records = reader.records();
    //!        let mut read_counter: u32 = 0;
    //!        let mut seq_region: BTreeMap<String, (u32,u32)> = BTreeMap::new();
    //!        let mut record_vec: Vec<Record> = Vec::new();
    //!        loop {
    //!            match records.next() {	
    //!                Some(Ok(mut record)) => {
    //!	               println!("next record");
    //!                    println!("Record id: {:?}", record.id);
    //!		       let source = record.source_map.source_name.clone().expect("issue collecting source name");
    //!		       let beginning = match record.source_map.get_start(&source) {
    //!		                        Some(value) => value.get_value(),
    //!				        _ => 0,
    //!					};
    //!		       let ending = match record.source_map.get_stop(&source) {
    //!		                        Some(value) => value.get_value(),
    //!					_ => 0,
    //!					};
    //!		       if ending + prev_end < beginning + prev_end {
    //!		          println!("debug: end value smaller is than the start {:?}", beginning);
    //!		          }
    //!		       seq_region.insert(source, (beginning + prev_end, ending + prev_end));
    //!		       record_vec.push(record);
    //!                    // Add additional fields to print if needed
    //!		       read_counter+=1;
    //!		       prev_end+=ending; // create the joined record if there are multiple
    //!                    },
    //!	            Some(Err(e)) => { println!("theres an err {:?}", e); },
    //!	            None => {
    //!	               println!("finished iteration");
    //!	                     break; },
    //!	            }
    //!            }
    //!        let output_file = format!("{}.gff", &args.filename);
    //!        if std::path::Path::new(&output_file).exists() {
    //!           println!("Deleting existing file: {}", &output_file);
    //!           std::fs::remove_file(&output_file).expect("NOOO");
    //!           }
    //!        gff_write(seq_region.clone(), record_vec, &output_file, true);
    //!        println!("Total records processed: {}", read_counter);
    //!        return Ok(());
    //!    }
    //!```
    //! Example to create a completely new record, use of setters or set_ functionality
    //!
    //! To write into GFF format requires gff_write(seq_region, record_vec, filename, true or false)
    //!
    //! The seq_region is the region of interest to save with name and DNA coordinates such as ``` seqregion.entry("source_1".to_string(), (1,897))```
    //! This makes it possible to save the whole file or to subset it
    //!
    //! record_vec is a list of the records.  If there is only one record, include this as a vec using ``` vec![record] ```
    //!
    //! The boolean true/false describes whether the DNA sequence should be included in the GFF3 file
    //!
    //! To write into genbank format requires gbk_write(seq_region, record_vec, filename), no true or false since genbank format will include the DNA sequence
    //!
    //!
    //! ```rust
    //!    use microBioRust::gbk::gff_write;
    //!    use microBioRust::gbk::RangeValue;
    //!    use microBioRust::gbk::Record;
    //!    use std::fs::File;
    //!    use std::collections::BTreeMap;
    //!
    //!     pub fn create_new_record() -> Result<(), anyhow::Error> {
    //!         let filename = format!("new_record.gff");
    //!         if std::path::Path::new(&filename).exists() {
    //!           std::fs::remove_file(&filename)?;
    //!           }
    //!	    let mut record = Record::new();
    //!	    let mut seq_region: BTreeMap<String, (u32,u32)> = BTreeMap::new();
    //!         //example from E.coli K12
    //!	    seq_region.insert("source_1".to_string(), (1,897));
    //!         //Add the source into SourceAttributes
    //!         record.source_map
    //!	         .set_counter("source_1".to_string())
    //!	         .set_start(RangeValue::Exact(1))
    //!	         .set_stop(RangeValue::Exact(897))
    //!	         .set_organism("Escherichia coli".to_string())
    //!	         .set_mol_type("DNA".to_string())
    //!	         .set_strain("K-12 substr. MG1655".to_string())
    //!		 .set_type_material("type strain of Escherichia coli K12".to_string())
    //!	         .set_db_xref("PRJNA57779".to_string());
    //!         //Add the features into FeatureAttributes, here we are setting two features, i.e. coding sequences or genes
    //!	    record.cds
    //!                  .set_counter("b3304".to_string())
    //!                  .set_start(RangeValue::Exact(1))
    //!                  .set_stop(RangeValue::Exact(354))
    //!                  .set_gene("rplR".to_string())
    //!                  .set_product("50S ribosomal subunit protein L18".to_string())
    //!                  .set_codon_start(1)
    //!                  .set_strand(-1);
    //!	    record.cds
    //!                  .set_counter("b3305".to_string())
    //!                  .set_start(RangeValue::Exact(364))
    //!                  .set_stop(RangeValue::Exact(897))
    //!                  .set_gene("rplF".to_string())
    //!                  .set_product("50S ribosomal subunit protein L6".to_string())
    //!                  .set_codon_start(1)
    //!                  .set_strand(-1);
    //!         //Add the sequences for the coding sequence (CDS) into SequenceAttributes
    //!	    record.seq_features
    //!	         .set_counter("b3304".to_string())
    //!		 .set_start(RangeValue::Exact(1))
    //!                 .set_stop(RangeValue::Exact(354))
    //!                 .set_sequence_ffn("ATGGATAAGAAATCTGCTCGTATCCGTCGTGCGACCCGCGCACGCCGCAAGCTCCAGGAG
    //!CTGGGCGCAACTCGCCTGGTGGTACATCGTACCCCGCGTCACATTTACGCACAGGTAATT
    //!GCACCGAACGGTTCTGAAGTTCTGGTAGCTGCTTCTACTGTAGAAAAAGCTATCGCTGAA
    //!CAACTGAAGTACACCGGTAACAAAGACGCGGCTGCAGCTGTGGGTAAAGCTGTCGCTGAA
    //!CGCGCTCTGGAAAAAGGCATCAAAGATGTATCCTTTGACCGTTCCGGGTTCCAATATCAT
    //!GGTCGTGTCCAGGCACTGGCAGATGCTGCCCGTGAAGCTGGCCTTCAGTTCTAA".to_string())
    //!                 .set_sequence_faa("MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAE
    //!QLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF".to_string())
    //!                 .set_codon_start(1)
    //!                 .set_strand(-1);
    //!	    record.seq_features
    //!	         .set_counter("bb3305".to_string())
    //!		 .set_start(RangeValue::Exact(364))
    //!                 .set_stop(RangeValue::Exact(897))
    //!                 .set_sequence_ffn("ATGTCTCGTGTTGCTAAAGCACCGGTCGTTGTTCCTGCCGGCGTTGACGTAAAAATCAAC
    //!GGTCAGGTTATTACGATCAAAGGTAAAAACGGCGAGCTGACTCGTACTCTCAACGATGCT
    //!GTTGAAGTTAAACATGCAGATAATACCCTGACCTTCGGTCCGCGTGATGGTTACGCAGAC
    //!GGTTGGGCACAGGCTGGTACCGCGCGTGCCCTGCTGAACTCAATGGTTATCGGTGTTACC
    //!GAAGGCTTCACTAAGAAGCTGCAGCTGGTTGGTGTAGGTTACCGTGCAGCGGTTAAAGGC
    //!AATGTGATTAACCTGTCTCTGGGTTTCTCTCATCCTGTTGACCATCAGCTGCCTGCGGGT
    //!ATCACTGCTGAATGTCCGACTCAGACTGAAATCGTGCTGAAAGGCGCTGATAAGCAGGTG
    //!ATCGGCCAGGTTGCAGCGGATCTGCGCGCCTACCGTCGTCCTGAGCCTTATAAAGGCAAG
    //!GGTGTTCGTTACGCCGACGAAGTCGTGCGTACCAAAGAGGCTAAGAAGAAGTAA".to_string())
    //!                 .set_sequence_faa("MSRVAKAPVVVPAGVDVKINGQVITIKGKNGELTRTLNDAVEVKHADNTLTFGPRDGYAD
    //!GWAQAGTARALLNSMVIGVTEGFTKKLQLVGVGYRAAVKGNVINLSLGFSHPVDHQLPAG
    //!ITAECPTQTEIVLKGADKQVIGQVAADLRAYRRPEPYKGKGVRYADEVVRTKEAKKK".to_string())
    //!                 .set_codon_start(1)
    //!                 .set_strand(-1);
    //!         //Add the full sequence of the entire record into the record.sequence
    //!	    record.sequence = "TTAGAACTGAAGGCCAGCTTCACGGGCAGCATCTGCCAGTGCCTGGACACGACCATGATA
    //!TTGGAACCCGGAACGGTCAAAGGATACATCTTTGATGCCTTTTTCCAGAGCGCGTTCAGC
    //!GACAGCTTTACCCACAGCTGCAGCCGCGTCTTTGTTACCGGTGTACTTCAGTTGTTCAGC
    //!GATAGCTTTTTCTACAGTAGAAGCAGCTACCAGAACTTCAGAACCGTTCGGTGCAATTAC
    //!CTGTGCGTAAATGTGACGCGGGGTACGATGTACCACCAGGCGAGTTGCGCCCAGCTCCTG
    //!GAGCTTGCGGCGTGCGCGGGTCGCACGACGGATACGAGCAGATTTCTTATCCATAGTGTT
    //!ACCTTACTTCTTCTTAGCCTCTTTGGTACGCACGACTTCGTCGGCGTAACGAACACCCTT
    //!GCCTTTATAAGGCTCAGGACGACGGTAGGCGCGCAGATCCGCTGCAACCTGGCCGATCAC
    //!CTGCTTATCAGCGCCTTTCAGCACGATTTCAGTCTGAGTCGGACATTCAGCAGTGATACC
    //!CGCAGGCAGCTGATGGTCAACAGGATGAGAGAAACCCAGAGACAGGTTAATCACATTGCC
    //!TTTAACCGCTGCACGGTAACCTACACCAACCAGCTGCAGCTTCTTAGTGAAGCCTTCGGT
    //!AACACCGATAACCATTGAGTTCAGCAGGGCACGCGCGGTACCAGCCTGTGCCCAACCGTC
    //!TGCGTAACCATCACGCGGACCGAAGGTCAGGGTATTATCTGCATGTTTAACTTCAACAGC
    //!ATCGTTGAGAGTACGAGTCAGCTCGCCGTTTTTACCTTTGATCGTAATAACCTGACCGTT
    //!GATTTTTACGTCAACGCCGGCAGGAACAACGACCGGTGCTTTAGCAACACGAGACAT".to_string();
    //!           gff_write(seq_region, vec![record], &filename, true);
    //!	   return Ok(());
    //!      }
    //!```
    //!
    use std::io::{self, Write};
    use std::fs;
    use regex::Regex;
    use itertools::Itertools;
    use std::vec::Vec;
    use std::str;
    use std::convert::AsRef;
    use protein_translate::translate;
    use std::path::Path;
    use bio::alphabets::dna::revcomp;
    use anyhow::anyhow;
    use std::collections::BTreeMap;
    use std::fs::{OpenOptions, File};
    use anyhow::Context;
    use std::collections::HashSet;
    use paste::paste;
    use std::convert::TryInto;
    use chrono::prelude::*;
    /// A Gbk reader.
    #[allow(unused_mut)]
    pub struct Records<B>
    where
        B: io::BufRead,
    {
        reader: Reader<B>,
        error_has_occurred: bool,
    }
    #[automatically_derived]
    #[allow(unused_mut)]
    impl<B: ::core::fmt::Debug> ::core::fmt::Debug for Records<B>
    where
        B: io::BufRead,
    {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "Records",
                "reader",
                &self.reader,
                "error_has_occurred",
                &&self.error_has_occurred,
            )
        }
    }
    impl<B> Records<B>
    where
        B: io::BufRead,
    {
        #[allow(unused_mut)]
        pub fn new(mut reader: Reader<B>) -> Self {
            Records {
                reader: reader,
                error_has_occurred: false,
            }
        }
    }
    impl<B> Iterator for Records<B>
    where
        B: io::BufRead,
    {
        type Item = Result<Record, anyhow::Error>;
        fn next(&mut self) -> Option<Self::Item> {
            if self.error_has_occurred {
                {
                    ::std::io::_print(
                        format_args!("error was encountered in iteration\n"),
                    );
                };
                None
            } else {
                let mut record = Record::new();
                match self.reader.read(&mut record) {
                    Ok(_) => if record.is_empty() { None } else { Some(Ok(record)) }
                    Err(err) => {
                        self.error_has_occurred = true;
                        Some(
                            Err(
                                ::anyhow::Error::msg(
                                    ::alloc::__export::must_use({
                                        let res = ::alloc::fmt::format(
                                            format_args!("next record read error {0:?}", err),
                                        );
                                        res
                                    }),
                                ),
                            ),
                        )
                    }
                }
            }
        }
    }
    pub trait GbkRead {
        fn read(&mut self, record: &mut Record) -> Result<Record, anyhow::Error>;
    }
    ///per line reader for the file
    pub struct Reader<B> {
        reader: B,
        line_buffer: String,
    }
    #[automatically_derived]
    impl<B: ::core::fmt::Debug> ::core::fmt::Debug for Reader<B> {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "Reader",
                "reader",
                &self.reader,
                "line_buffer",
                &&self.line_buffer,
            )
        }
    }
    #[automatically_derived]
    impl<B: ::core::default::Default> ::core::default::Default for Reader<B> {
        #[inline]
        fn default() -> Reader<B> {
            Reader {
                reader: ::core::default::Default::default(),
                line_buffer: ::core::default::Default::default(),
            }
        }
    }
    impl Reader<io::BufReader<fs::File>> {
        /// Read Gbk from given file path in given format.
        pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(
            path: P,
        ) -> anyhow::Result<Self> {
            fs::File::open(&path)
                .map(Reader::new)
                .with_context(|| ::alloc::__export::must_use({
                    let res = ::alloc::fmt::format(
                        format_args!("Failed to read Gbk from {0:#?}", path),
                    );
                    res
                }))
        }
    }
    impl<R> Reader<io::BufReader<R>>
    where
        R: io::Read,
    {
        pub fn new(reader: R) -> Self {
            Reader {
                reader: io::BufReader::new(reader),
                line_buffer: String::new(),
            }
        }
    }
    impl<B> Reader<B>
    where
        B: io::BufRead,
    {
        pub fn from_bufread(bufreader: B) -> Self {
            Reader {
                reader: bufreader,
                line_buffer: String::new(),
            }
        }
        pub fn records(self) -> Records<B> {
            Records {
                reader: self,
                error_has_occurred: false,
            }
        }
    }
    ///main gbk parser
    impl<'a, B> GbkRead for Reader<B>
    where
        B: io::BufRead,
    {
        #[allow(unused_mut)]
        #[allow(unused_variables)]
        #[allow(unused_assignments)]
        fn read(&mut self, record: &mut Record) -> Result<Record, anyhow::Error> {
            record.rec_clear();
            let mut sequences = String::new();
            let mut source_map = SourceAttributeBuilder::new();
            let mut cds = FeatureAttributeBuilder::new();
            let mut seq_features = SequenceAttributeBuilder::new();
            let mut cds_counter: i32 = 0;
            let mut source_counter: i32 = 0;
            let mut prev_end: u32 = 0;
            let mut organism = String::new();
            let mut mol_type = String::new();
            let mut strain = String::new();
            let mut source_name = String::new();
            let mut type_material = String::new();
            let mut theend: u32 = 0;
            let mut thestart: u32 = 0;
            let mut db_xref = String::new();
            if self.line_buffer.is_empty() {
                self.reader.read_line(&mut self.line_buffer)?;
                if self.line_buffer.is_empty() {
                    return Ok(record.to_owned());
                }
            }
            'outer: while !self.line_buffer.is_empty() {
                if self.line_buffer.starts_with("LOCUS") {
                    record.rec_clear();
                    let mut header_fields: Vec<&str> = self
                        .line_buffer
                        .split_whitespace()
                        .collect();
                    let mut header_iter = header_fields.iter();
                    header_iter.next();
                    record.id = header_iter
                        .next()
                        .ok_or_else(|| ::anyhow::__private::must_use({
                            let error = ::anyhow::__private::format_err(
                                format_args!("missing record id"),
                            );
                            error
                        }))?
                        .to_string();
                    let lens = header_iter
                        .next()
                        .ok_or_else(|| ::anyhow::__private::must_use({
                            let error = ::anyhow::__private::format_err(
                                format_args!("missing record length"),
                            );
                            error
                        }))?
                        .to_string();
                    record.length = lens.trim().parse::<u32>()?;
                    self.line_buffer.clear();
                }
                if self.line_buffer.starts_with("     source") {
                    let re = Regex::new(r"([0-9]+)[[:punct:]]+([0-9]+)")?;
                    let location = re
                        .captures(&self.line_buffer)
                        .ok_or_else(|| ::anyhow::__private::must_use({
                            let error = ::anyhow::__private::format_err(
                                format_args!("missing location"),
                            );
                            error
                        }))?;
                    let start = &location[1];
                    let end = &location[2];
                    thestart = start.trim().parse::<u32>()?;
                    source_counter += 1;
                    source_name = ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("source_{0}_{1}", record.id, source_counter),
                            );
                            res
                        })
                        .to_string();
                    thestart += prev_end;
                    theend = end.trim().parse::<u32>()? + prev_end;
                    loop {
                        self.line_buffer.clear();
                        self.reader.read_line(&mut self.line_buffer)?;
                        if self.line_buffer.starts_with("     CDS") {
                            record
                                .source_map
                                .set_counter(source_name.to_string())
                                .set_start(RangeValue::Exact(thestart))
                                .set_stop(RangeValue::Exact(theend))
                                .set_organism(organism.clone())
                                .set_mol_type(mol_type.clone())
                                .set_strain(strain.clone())
                                .set_type_material(type_material.clone())
                                .set_db_xref(db_xref.clone());
                            continue 'outer;
                        }
                        if self.line_buffer.contains("/organism") {
                            let org: Vec<&str> = self.line_buffer.split('\"').collect();
                            organism = org[1].to_string();
                        }
                        if self.line_buffer.contains("/mol_type") {
                            let mol: Vec<&str> = self.line_buffer.split('\"').collect();
                            mol_type = mol[1].to_string();
                        }
                        if self.line_buffer.contains("/strain") {
                            let stra: Vec<&str> = self.line_buffer.split('\"').collect();
                            strain = stra[1].to_string();
                        }
                        if self.line_buffer.contains("/type_material") {
                            let mat: Vec<&str> = self.line_buffer.split('\"').collect();
                            type_material = mat[1].to_string();
                        }
                        if self.line_buffer.contains("/db_xref") {
                            let db: Vec<&str> = self.line_buffer.split('\"').collect();
                            db_xref = db[1].to_string();
                        }
                    }
                }
                if self.line_buffer.starts_with("     CDS") {
                    let mut startiter: Vec<_> = Vec::new();
                    let mut enditer: Vec<_> = Vec::new();
                    let mut thestart: u32 = 0;
                    let mut thend: u32 = 0;
                    let mut joined: bool = false;
                    let joined = if self.line_buffer.contains("join") {
                        true
                    } else {
                        false
                    };
                    let re = Regex::new(r"([0-9]+)[[:punct:]]+([0-9]+)")?;
                    for cap in re.captures_iter(&self.line_buffer) {
                        cds_counter += 1;
                        thestart = cap[1]
                            .parse()
                            .expect("failed to match and parse numerical start");
                        theend = cap[2]
                            .parse()
                            .expect("failed to match and parse numerical end");
                        startiter.push(thestart);
                        enditer.push(theend);
                    }
                    let mut gene = String::new();
                    let mut product = String::new();
                    let strand: i8 = if self.line_buffer.contains("complement") {
                        -1
                    } else {
                        1
                    };
                    let mut locus_tag = String::new();
                    let mut codon_start: u8 = 1;
                    loop {
                        self.line_buffer.clear();
                        self.reader.read_line(&mut self.line_buffer)?;
                        if self.line_buffer.contains("/locus_tag=") {
                            let loctag: Vec<&str> = self
                                .line_buffer
                                .split('\"')
                                .collect();
                            locus_tag = loctag[1].to_string();
                        }
                        if self.line_buffer.contains("/codon_start") {
                            let codstart: Vec<&str> = self
                                .line_buffer
                                .split('=')
                                .collect();
                            let valstart = codstart[1].trim().parse::<u8>()?;
                            codon_start = valstart;
                        }
                        if self.line_buffer.contains("/gene=") {
                            let gen: Vec<&str> = self.line_buffer.split('\"').collect();
                            gene = gen[1].to_string();
                        }
                        if self.line_buffer.contains("/product") {
                            let prod: Vec<&str> = self.line_buffer.split('\"').collect();
                            product = substitute_odd_punctuation(prod[1].to_string())?;
                        }
                        if self.line_buffer.starts_with("     CDS")
                            || self.line_buffer.starts_with("ORIGIN")
                            || self.line_buffer.starts_with("     gene")
                            || self.line_buffer.starts_with("     misc_feature")
                        {
                            if locus_tag.is_empty() {
                                locus_tag = ::alloc::__export::must_use({
                                        let res = ::alloc::fmt::format(
                                            format_args!("CDS_{0}", cds_counter),
                                        );
                                        res
                                    })
                                    .to_string();
                            }
                            if joined {
                                for (i, m) in startiter.iter().enumerate() {
                                    let loc_tag = ::alloc::__export::must_use({
                                        let res = ::alloc::fmt::format(
                                            format_args!("{0}_{1}", locus_tag.clone(), i),
                                        );
                                        res
                                    });
                                    record
                                        .cds
                                        .set_counter(loc_tag)
                                        .set_start(RangeValue::Exact(*m))
                                        .set_stop(RangeValue::Exact(enditer[i]))
                                        .set_gene(gene.to_string())
                                        .set_product(product.to_string())
                                        .set_codon_start(codon_start)
                                        .set_strand(strand);
                                }
                                continue 'outer;
                            } else {
                                record
                                    .cds
                                    .set_counter(locus_tag.clone())
                                    .set_start(RangeValue::Exact(thestart))
                                    .set_stop(RangeValue::Exact(theend))
                                    .set_gene(gene.to_string())
                                    .set_product(product.to_string())
                                    .set_codon_start(codon_start)
                                    .set_strand(strand);
                                continue 'outer;
                            }
                        }
                    }
                }
                if self.line_buffer.starts_with("ORIGIN") {
                    let mut sequences = String::new();
                    let result_seq = loop {
                        self.line_buffer.clear();
                        self.reader.read_line(&mut self.line_buffer)?;
                        if self.line_buffer.starts_with("//") {
                            break sequences;
                        } else {
                            let s: Vec<&str> = self
                                .line_buffer
                                .split_whitespace()
                                .collect();
                            let s = &s[1..];
                            let sequence = s.iter().join("");
                            sequences.push_str(&sequence);
                        }
                    };
                    record.sequence = result_seq.to_string();
                    let mut iterablecount: u32 = 0;
                    for (key, val) in record.cds.iter_sorted() {
                        let (
                            mut a,
                            mut b,
                            mut c,
                            mut d,
                        ): (Option<u32>, Option<u32>, Option<i8>, Option<u8>) = (
                            None,
                            None,
                            None,
                            None,
                        );
                        for value in val {
                            match value {
                                FeatureAttributes::Start { value } => {
                                    a = match value {
                                        RangeValue::Exact(v) => Some(*v),
                                        RangeValue::LessThan(v) => Some(*v),
                                        RangeValue::GreaterThan(v) => Some(*v),
                                    };
                                }
                                FeatureAttributes::Stop { value } => {
                                    b = match value {
                                        RangeValue::Exact(v) => Some(*v),
                                        RangeValue::LessThan(v) => Some(*v),
                                        RangeValue::GreaterThan(v) => Some(*v),
                                    };
                                }
                                FeatureAttributes::Strand { value } => {
                                    c = match value {
                                        value => Some(*value),
                                    };
                                }
                                FeatureAttributes::CodonStart { value } => {
                                    d = match value {
                                        value => Some(value.clone()),
                                    };
                                }
                                _ => {}
                            }
                        }
                        let sta = a
                            .map(|o| o as usize)
                            .ok_or(
                                ::anyhow::__private::must_use({
                                    let error = ::anyhow::__private::format_err(
                                        format_args!("No value for start"),
                                    );
                                    error
                                }),
                            )?;
                        let sto = b
                            .map(|t| t as usize)
                            .ok_or(
                                ::anyhow::__private::must_use({
                                    let error = ::anyhow::__private::format_err(
                                        format_args!("No value for stop"),
                                    );
                                    error
                                }),
                            )? - 1;
                        let stra = c
                            .map(|u| u as i8)
                            .ok_or(
                                ::anyhow::__private::must_use({
                                    let error = ::anyhow::__private::format_err(
                                        format_args!("No value for strand"),
                                    );
                                    error
                                }),
                            )?;
                        let cod = d
                            .map(|v| v as usize - 1)
                            .ok_or(
                                ::anyhow::__private::must_use({
                                    let error = ::anyhow::__private::format_err(
                                        format_args!("No value for strand"),
                                    );
                                    error
                                }),
                            )?;
                        let star = sta.try_into()?;
                        let stow = sto.try_into()?;
                        let codd = cod.try_into()?;
                        let mut sliced_sequence: &str = "";
                        if stra == -1 {
                            if cod > 1 {
                                if sto + 1 <= record.sequence.len() {
                                    sliced_sequence = &record.sequence[sta + cod..sto + 1];
                                } else {
                                    sliced_sequence = &record.sequence[sta + cod..sto];
                                }
                            } else {
                                if sto + 1 <= record.sequence.len() {
                                    sliced_sequence = &record.sequence[sta..sto + 1];
                                } else {
                                    sliced_sequence = &record.sequence[sta..sto];
                                }
                            }
                            let cds_char = sliced_sequence;
                            let prot_seq = translate(&revcomp(cds_char.as_bytes()));
                            let parts: Vec<&str> = prot_seq.split('*').collect();
                            record
                                .seq_features
                                .set_counter(key.to_string())
                                .set_start(RangeValue::Exact(star))
                                .set_stop(RangeValue::Exact(stow))
                                .set_sequence_ffn(cds_char.to_string())
                                .set_sequence_faa(parts[0].to_string())
                                .set_codon_start(codd)
                                .set_strand(stra);
                        } else {
                            if cod > 1 {
                                sliced_sequence = &record.sequence[sta + cod - 1..sto];
                            } else {
                                sliced_sequence = &record.sequence[sta - 1..sto];
                            }
                            let cds_char = sliced_sequence;
                            let prot_seq = translate(cds_char.as_bytes());
                            let parts: Vec<&str> = prot_seq.split('*').collect();
                            record
                                .seq_features
                                .set_counter(key.to_string())
                                .set_start(RangeValue::Exact(star))
                                .set_stop(RangeValue::Exact(stow))
                                .set_sequence_ffn(cds_char.to_string())
                                .set_sequence_faa(parts[0].to_string())
                                .set_codon_start(codd)
                                .set_strand(stra);
                        }
                    }
                    return Ok(record.to_owned());
                }
                self.line_buffer.clear();
                self.reader.read_line(&mut self.line_buffer)?;
            }
            Ok(record.to_owned())
        }
    }
    ///stores a value for start or stop (end) which can be denoted as a < value or > value.
    pub enum RangeValue {
        Exact(u32),
        LessThan(u32),
        GreaterThan(u32),
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for RangeValue {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match self {
                RangeValue::Exact(__self_0) => {
                    ::core::fmt::Formatter::debug_tuple_field1_finish(
                        f,
                        "Exact",
                        &__self_0,
                    )
                }
                RangeValue::LessThan(__self_0) => {
                    ::core::fmt::Formatter::debug_tuple_field1_finish(
                        f,
                        "LessThan",
                        &__self_0,
                    )
                }
                RangeValue::GreaterThan(__self_0) => {
                    ::core::fmt::Formatter::debug_tuple_field1_finish(
                        f,
                        "GreaterThan",
                        &__self_0,
                    )
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::hash::Hash for RangeValue {
        #[inline]
        fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            ::core::hash::Hash::hash(&__self_discr, state);
            match self {
                RangeValue::Exact(__self_0) => ::core::hash::Hash::hash(__self_0, state),
                RangeValue::LessThan(__self_0) => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                RangeValue::GreaterThan(__self_0) => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::marker::StructuralPartialEq for RangeValue {}
    #[automatically_derived]
    impl ::core::cmp::PartialEq for RangeValue {
        #[inline]
        fn eq(&self, other: &RangeValue) -> bool {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            let __arg1_discr = ::core::intrinsics::discriminant_value(other);
            __self_discr == __arg1_discr
                && match (self, other) {
                    (RangeValue::Exact(__self_0), RangeValue::Exact(__arg1_0)) => {
                        __self_0 == __arg1_0
                    }
                    (RangeValue::LessThan(__self_0), RangeValue::LessThan(__arg1_0)) => {
                        __self_0 == __arg1_0
                    }
                    (
                        RangeValue::GreaterThan(__self_0),
                        RangeValue::GreaterThan(__arg1_0),
                    ) => __self_0 == __arg1_0,
                    _ => unsafe { ::core::intrinsics::unreachable() }
                }
        }
    }
    #[automatically_derived]
    impl ::core::cmp::Eq for RangeValue {
        #[inline]
        #[doc(hidden)]
        #[coverage(off)]
        fn assert_receiver_is_total_eq(&self) -> () {
            let _: ::core::cmp::AssertParamIsEq<u32>;
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for RangeValue {
        #[inline]
        fn clone(&self) -> RangeValue {
            match self {
                RangeValue::Exact(__self_0) => {
                    RangeValue::Exact(::core::clone::Clone::clone(__self_0))
                }
                RangeValue::LessThan(__self_0) => {
                    RangeValue::LessThan(::core::clone::Clone::clone(__self_0))
                }
                RangeValue::GreaterThan(__self_0) => {
                    RangeValue::GreaterThan(::core::clone::Clone::clone(__self_0))
                }
            }
        }
    }
    impl RangeValue {
        pub fn get_value(&self) -> u32 {
            match self {
                RangeValue::Exact(value) => *value,
                RangeValue::LessThan(value) => *value,
                RangeValue::GreaterThan(value) => *value,
            }
        }
    }
    pub enum SourceAttributes {
        Start { value: RangeValue },
        Stop { value: RangeValue },
        Organism { value: String },
        MolType { value: String },
        Strain { value: String },
        CultureCollection { value: String },
        TypeMaterial { value: String },
        DbXref { value: String },
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for SourceAttributes {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match self {
                SourceAttributes::Start { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Start",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::Stop { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Stop",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::Organism { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Organism",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::MolType { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "MolType",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::Strain { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Strain",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::CultureCollection { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "CultureCollection",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::TypeMaterial { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "TypeMaterial",
                        "value",
                        &__self_0,
                    )
                }
                SourceAttributes::DbXref { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "DbXref",
                        "value",
                        &__self_0,
                    )
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::cmp::Eq for SourceAttributes {
        #[inline]
        #[doc(hidden)]
        #[coverage(off)]
        fn assert_receiver_is_total_eq(&self) -> () {
            let _: ::core::cmp::AssertParamIsEq<RangeValue>;
            let _: ::core::cmp::AssertParamIsEq<String>;
        }
    }
    #[automatically_derived]
    impl ::core::marker::StructuralPartialEq for SourceAttributes {}
    #[automatically_derived]
    impl ::core::cmp::PartialEq for SourceAttributes {
        #[inline]
        fn eq(&self, other: &SourceAttributes) -> bool {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            let __arg1_discr = ::core::intrinsics::discriminant_value(other);
            __self_discr == __arg1_discr
                && match (self, other) {
                    (
                        SourceAttributes::Start { value: __self_0 },
                        SourceAttributes::Start { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::Stop { value: __self_0 },
                        SourceAttributes::Stop { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::Organism { value: __self_0 },
                        SourceAttributes::Organism { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::MolType { value: __self_0 },
                        SourceAttributes::MolType { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::Strain { value: __self_0 },
                        SourceAttributes::Strain { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::CultureCollection { value: __self_0 },
                        SourceAttributes::CultureCollection { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::TypeMaterial { value: __self_0 },
                        SourceAttributes::TypeMaterial { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SourceAttributes::DbXref { value: __self_0 },
                        SourceAttributes::DbXref { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    _ => unsafe { ::core::intrinsics::unreachable() }
                }
        }
    }
    #[automatically_derived]
    impl ::core::hash::Hash for SourceAttributes {
        #[inline]
        fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            ::core::hash::Hash::hash(&__self_discr, state);
            match self {
                SourceAttributes::Start { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::Stop { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::Organism { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::MolType { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::Strain { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::CultureCollection { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::TypeMaterial { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SourceAttributes::DbXref { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for SourceAttributes {
        #[inline]
        fn clone(&self) -> SourceAttributes {
            match self {
                SourceAttributes::Start { value: __self_0 } => {
                    SourceAttributes::Start {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::Stop { value: __self_0 } => {
                    SourceAttributes::Stop {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::Organism { value: __self_0 } => {
                    SourceAttributes::Organism {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::MolType { value: __self_0 } => {
                    SourceAttributes::MolType {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::Strain { value: __self_0 } => {
                    SourceAttributes::Strain {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::CultureCollection { value: __self_0 } => {
                    SourceAttributes::CultureCollection {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::TypeMaterial { value: __self_0 } => {
                    SourceAttributes::TypeMaterial {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SourceAttributes::DbXref { value: __self_0 } => {
                    SourceAttributes::DbXref {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
            }
        }
    }
    impl SourceAttributeBuilder {
        pub fn get_start(&self, key: &str) -> Option<&RangeValue> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::Start { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_stop(&self, key: &str) -> Option<&RangeValue> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::Stop { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_organism(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::Organism { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_mol_type(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::MolType { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_strain(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::Strain { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_type_material(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::TypeMaterial { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_db_xref(&self, key: &str) -> Option<&String> {
            self.source_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SourceAttributes::DbXref { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
    }
    ///builder for the source information on a per record basis
    pub struct SourceAttributeBuilder {
        pub source_attributes: BTreeMap<String, HashSet<SourceAttributes>>,
        pub source_name: Option<String>,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for SourceAttributeBuilder {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "SourceAttributeBuilder",
                "source_attributes",
                &self.source_attributes,
                "source_name",
                &&self.source_name,
            )
        }
    }
    #[automatically_derived]
    impl ::core::default::Default for SourceAttributeBuilder {
        #[inline]
        fn default() -> SourceAttributeBuilder {
            SourceAttributeBuilder {
                source_attributes: ::core::default::Default::default(),
                source_name: ::core::default::Default::default(),
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for SourceAttributeBuilder {
        #[inline]
        fn clone(&self) -> SourceAttributeBuilder {
            SourceAttributeBuilder {
                source_attributes: ::core::clone::Clone::clone(&self.source_attributes),
                source_name: ::core::clone::Clone::clone(&self.source_name),
            }
        }
    }
    impl SourceAttributeBuilder {
        pub fn set_source_name(&mut self, name: String) {
            self.source_name = Some(name);
        }
        pub fn get_source_name(&self) -> Option<&String> {
            self.source_name.as_ref()
        }
        pub fn add_source_attribute(
            &mut self,
            key: String,
            attribute: SourceAttributes,
        ) {
            self.source_attributes
                .entry(key)
                .or_insert_with(HashSet::new)
                .insert(attribute);
        }
        pub fn get_source_attributes(
            &self,
            key: &str,
        ) -> Option<&HashSet<SourceAttributes>> {
            self.source_attributes.get(key)
        }
    }
    impl SourceAttributeBuilder {
        pub fn new() -> Self {
            SourceAttributeBuilder {
                source_attributes: BTreeMap::new(),
                source_name: None,
            }
        }
        pub fn set_counter(&mut self, counter: String) -> &mut Self {
            self.source_name = Some(counter);
            self
        }
        pub fn insert_to(&mut self, value: SourceAttributes) {
            if let Some(counter) = &self.source_name {
                self.source_attributes
                    .entry(counter.to_string())
                    .or_insert_with(HashSet::new)
                    .insert(value);
            } else {
                {
                    ::core::panicking::panic_fmt(format_args!("Counter key not set"));
                };
            }
        }
        pub fn set_start(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(SourceAttributes::Start { value });
            self
        }
        pub fn set_stop(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(SourceAttributes::Stop { value });
            self
        }
        pub fn set_organism(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::Organism {
                value,
            });
            self
        }
        pub fn set_mol_type(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::MolType { value });
            self
        }
        pub fn set_strain(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::Strain { value });
            self
        }
        pub fn set_type_material(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::TypeMaterial {
                value,
            });
            self
        }
        pub fn set_db_xref(&mut self, value: String) -> &mut Self {
            self.insert_to(SourceAttributes::DbXref { value });
            self
        }
        pub fn build(self) -> BTreeMap<String, HashSet<SourceAttributes>> {
            self.source_attributes
        }
        pub fn iter_sorted(
            &self,
        ) -> std::collections::btree_map::Iter<String, HashSet<SourceAttributes>> {
            self.source_attributes.iter()
        }
        pub fn default() -> Self {
            SourceAttributeBuilder {
                source_attributes: BTreeMap::new(),
                source_name: None,
            }
        }
    }
    ///attributes for each feature, cds or gene
    pub enum FeatureAttributes {
        Start { value: RangeValue },
        Stop { value: RangeValue },
        Gene { value: String },
        Product { value: String },
        CodonStart { value: u8 },
        Strand { value: i8 },
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for FeatureAttributes {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match self {
                FeatureAttributes::Start { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Start",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::Stop { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Stop",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::Gene { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Gene",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::Product { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Product",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::CodonStart { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "CodonStart",
                        "value",
                        &__self_0,
                    )
                }
                FeatureAttributes::Strand { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Strand",
                        "value",
                        &__self_0,
                    )
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::cmp::Eq for FeatureAttributes {
        #[inline]
        #[doc(hidden)]
        #[coverage(off)]
        fn assert_receiver_is_total_eq(&self) -> () {
            let _: ::core::cmp::AssertParamIsEq<RangeValue>;
            let _: ::core::cmp::AssertParamIsEq<String>;
            let _: ::core::cmp::AssertParamIsEq<u8>;
            let _: ::core::cmp::AssertParamIsEq<i8>;
        }
    }
    #[automatically_derived]
    impl ::core::hash::Hash for FeatureAttributes {
        #[inline]
        fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            ::core::hash::Hash::hash(&__self_discr, state);
            match self {
                FeatureAttributes::Start { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::Stop { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::Gene { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::Product { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::CodonStart { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                FeatureAttributes::Strand { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::marker::StructuralPartialEq for FeatureAttributes {}
    #[automatically_derived]
    impl ::core::cmp::PartialEq for FeatureAttributes {
        #[inline]
        fn eq(&self, other: &FeatureAttributes) -> bool {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            let __arg1_discr = ::core::intrinsics::discriminant_value(other);
            __self_discr == __arg1_discr
                && match (self, other) {
                    (
                        FeatureAttributes::Start { value: __self_0 },
                        FeatureAttributes::Start { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::Stop { value: __self_0 },
                        FeatureAttributes::Stop { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::Gene { value: __self_0 },
                        FeatureAttributes::Gene { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::Product { value: __self_0 },
                        FeatureAttributes::Product { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::CodonStart { value: __self_0 },
                        FeatureAttributes::CodonStart { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        FeatureAttributes::Strand { value: __self_0 },
                        FeatureAttributes::Strand { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    _ => unsafe { ::core::intrinsics::unreachable() }
                }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for FeatureAttributes {
        #[inline]
        fn clone(&self) -> FeatureAttributes {
            match self {
                FeatureAttributes::Start { value: __self_0 } => {
                    FeatureAttributes::Start {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::Stop { value: __self_0 } => {
                    FeatureAttributes::Stop {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::Gene { value: __self_0 } => {
                    FeatureAttributes::Gene {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::Product { value: __self_0 } => {
                    FeatureAttributes::Product {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::CodonStart { value: __self_0 } => {
                    FeatureAttributes::CodonStart {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                FeatureAttributes::Strand { value: __self_0 } => {
                    FeatureAttributes::Strand {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
            }
        }
    }
    impl FeatureAttributeBuilder {
        pub fn get_start(&self, key: &str) -> Option<&RangeValue> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Start { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_stop(&self, key: &str) -> Option<&RangeValue> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Stop { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_gene(&self, key: &str) -> Option<&String> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Gene { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_product(&self, key: &str) -> Option<&String> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Product { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_codon_start(&self, key: &str) -> Option<&u8> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::CodonStart { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_strand(&self, key: &str) -> Option<&i8> {
            self.attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let FeatureAttributes::Strand { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
    }
    ///builder for the feature information on a per coding sequence (CDS) basis
    pub struct FeatureAttributeBuilder {
        pub attributes: BTreeMap<String, HashSet<FeatureAttributes>>,
        locus_tag: Option<String>,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for FeatureAttributeBuilder {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "FeatureAttributeBuilder",
                "attributes",
                &self.attributes,
                "locus_tag",
                &&self.locus_tag,
            )
        }
    }
    #[automatically_derived]
    impl ::core::default::Default for FeatureAttributeBuilder {
        #[inline]
        fn default() -> FeatureAttributeBuilder {
            FeatureAttributeBuilder {
                attributes: ::core::default::Default::default(),
                locus_tag: ::core::default::Default::default(),
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for FeatureAttributeBuilder {
        #[inline]
        fn clone(&self) -> FeatureAttributeBuilder {
            FeatureAttributeBuilder {
                attributes: ::core::clone::Clone::clone(&self.attributes),
                locus_tag: ::core::clone::Clone::clone(&self.locus_tag),
            }
        }
    }
    impl FeatureAttributeBuilder {
        pub fn new() -> Self {
            FeatureAttributeBuilder {
                attributes: BTreeMap::new(),
                locus_tag: None,
            }
        }
        pub fn set_counter(&mut self, counter: String) -> &mut Self {
            self.locus_tag = Some(counter);
            self
        }
        pub fn insert_to(&mut self, value: FeatureAttributes) {
            if let Some(counter) = &self.locus_tag {
                self.attributes
                    .entry(counter.to_string())
                    .or_insert_with(HashSet::new)
                    .insert(value);
            } else {
                {
                    ::core::panicking::panic_fmt(format_args!("Counter key not set"));
                };
            }
        }
        pub fn set_start(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(FeatureAttributes::Start { value });
            self
        }
        pub fn set_stop(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(FeatureAttributes::Stop { value });
            self
        }
        pub fn set_gene(&mut self, value: String) -> &mut Self {
            self.insert_to(FeatureAttributes::Gene { value });
            self
        }
        pub fn set_product(&mut self, value: String) -> &mut Self {
            self.insert_to(FeatureAttributes::Product {
                value,
            });
            self
        }
        pub fn set_codon_start(&mut self, value: u8) -> &mut Self {
            self.insert_to(FeatureAttributes::CodonStart {
                value,
            });
            self
        }
        pub fn set_strand(&mut self, value: i8) -> &mut Self {
            self.insert_to(FeatureAttributes::Strand { value });
            self
        }
        pub fn build(self) -> BTreeMap<String, HashSet<FeatureAttributes>> {
            self.attributes
        }
        pub fn iter_sorted(
            &self,
        ) -> std::collections::btree_map::Iter<String, HashSet<FeatureAttributes>> {
            self.attributes.iter()
        }
        pub fn default() -> Self {
            FeatureAttributeBuilder {
                attributes: BTreeMap::new(),
                locus_tag: None,
            }
        }
    }
    ///stores the sequences of the coding sequences (genes) and proteins. Also stores start, stop, codon_start and strand information
    pub enum SequenceAttributes {
        Start { value: RangeValue },
        Stop { value: RangeValue },
        SequenceFfn { value: String },
        SequenceFaa { value: String },
        CodonStart { value: u8 },
        Strand { value: i8 },
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for SequenceAttributes {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            match self {
                SequenceAttributes::Start { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Start",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::Stop { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Stop",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::SequenceFfn { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "SequenceFfn",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::SequenceFaa { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "SequenceFaa",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::CodonStart { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "CodonStart",
                        "value",
                        &__self_0,
                    )
                }
                SequenceAttributes::Strand { value: __self_0 } => {
                    ::core::fmt::Formatter::debug_struct_field1_finish(
                        f,
                        "Strand",
                        "value",
                        &__self_0,
                    )
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::cmp::Eq for SequenceAttributes {
        #[inline]
        #[doc(hidden)]
        #[coverage(off)]
        fn assert_receiver_is_total_eq(&self) -> () {
            let _: ::core::cmp::AssertParamIsEq<RangeValue>;
            let _: ::core::cmp::AssertParamIsEq<String>;
            let _: ::core::cmp::AssertParamIsEq<u8>;
            let _: ::core::cmp::AssertParamIsEq<i8>;
        }
    }
    #[automatically_derived]
    impl ::core::marker::StructuralPartialEq for SequenceAttributes {}
    #[automatically_derived]
    impl ::core::cmp::PartialEq for SequenceAttributes {
        #[inline]
        fn eq(&self, other: &SequenceAttributes) -> bool {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            let __arg1_discr = ::core::intrinsics::discriminant_value(other);
            __self_discr == __arg1_discr
                && match (self, other) {
                    (
                        SequenceAttributes::Start { value: __self_0 },
                        SequenceAttributes::Start { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::Stop { value: __self_0 },
                        SequenceAttributes::Stop { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::SequenceFfn { value: __self_0 },
                        SequenceAttributes::SequenceFfn { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::SequenceFaa { value: __self_0 },
                        SequenceAttributes::SequenceFaa { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::CodonStart { value: __self_0 },
                        SequenceAttributes::CodonStart { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    (
                        SequenceAttributes::Strand { value: __self_0 },
                        SequenceAttributes::Strand { value: __arg1_0 },
                    ) => __self_0 == __arg1_0,
                    _ => unsafe { ::core::intrinsics::unreachable() }
                }
        }
    }
    #[automatically_derived]
    impl ::core::hash::Hash for SequenceAttributes {
        #[inline]
        fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
            let __self_discr = ::core::intrinsics::discriminant_value(self);
            ::core::hash::Hash::hash(&__self_discr, state);
            match self {
                SequenceAttributes::Start { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::Stop { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::SequenceFfn { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::SequenceFaa { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::CodonStart { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
                SequenceAttributes::Strand { value: __self_0 } => {
                    ::core::hash::Hash::hash(__self_0, state)
                }
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for SequenceAttributes {
        #[inline]
        fn clone(&self) -> SequenceAttributes {
            match self {
                SequenceAttributes::Start { value: __self_0 } => {
                    SequenceAttributes::Start {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::Stop { value: __self_0 } => {
                    SequenceAttributes::Stop {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::SequenceFfn { value: __self_0 } => {
                    SequenceAttributes::SequenceFfn {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::SequenceFaa { value: __self_0 } => {
                    SequenceAttributes::SequenceFaa {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::CodonStart { value: __self_0 } => {
                    SequenceAttributes::CodonStart {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
                SequenceAttributes::Strand { value: __self_0 } => {
                    SequenceAttributes::Strand {
                        value: ::core::clone::Clone::clone(__self_0),
                    }
                }
            }
        }
    }
    impl SequenceAttributeBuilder {
        pub fn get_start(&self, key: &str) -> Option<&RangeValue> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::Start { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_stop(&self, key: &str) -> Option<&RangeValue> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::Stop { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_sequence_ffn(&self, key: &str) -> Option<&String> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::SequenceFfn { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_sequence_faa(&self, key: &str) -> Option<&String> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::SequenceFaa { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_codon_start(&self, key: &str) -> Option<&u8> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::CodonStart { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
        pub fn get_strand(&self, key: &str) -> Option<&i8> {
            self.seq_attributes
                .get(key)
                .and_then(|set| {
                    set.iter()
                        .find_map(|attr| {
                            if let SequenceAttributes::Strand { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                })
        }
    }
    ///builder for the sequence information on a per coding sequence (CDS) basis
    pub struct SequenceAttributeBuilder {
        pub seq_attributes: BTreeMap<String, HashSet<SequenceAttributes>>,
        pub locus_tag: Option<String>,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for SequenceAttributeBuilder {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field2_finish(
                f,
                "SequenceAttributeBuilder",
                "seq_attributes",
                &self.seq_attributes,
                "locus_tag",
                &&self.locus_tag,
            )
        }
    }
    #[automatically_derived]
    impl ::core::default::Default for SequenceAttributeBuilder {
        #[inline]
        fn default() -> SequenceAttributeBuilder {
            SequenceAttributeBuilder {
                seq_attributes: ::core::default::Default::default(),
                locus_tag: ::core::default::Default::default(),
            }
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for SequenceAttributeBuilder {
        #[inline]
        fn clone(&self) -> SequenceAttributeBuilder {
            SequenceAttributeBuilder {
                seq_attributes: ::core::clone::Clone::clone(&self.seq_attributes),
                locus_tag: ::core::clone::Clone::clone(&self.locus_tag),
            }
        }
    }
    impl SequenceAttributeBuilder {
        pub fn new() -> Self {
            SequenceAttributeBuilder {
                seq_attributes: BTreeMap::new(),
                locus_tag: None,
            }
        }
        pub fn set_counter(&mut self, counter: String) -> &mut Self {
            self.locus_tag = Some(counter);
            self
        }
        pub fn insert_to(&mut self, value: SequenceAttributes) {
            if let Some(counter) = &self.locus_tag {
                self.seq_attributes
                    .entry(counter.to_string())
                    .or_insert_with(HashSet::new)
                    .insert(value);
            } else {
                {
                    ::core::panicking::panic_fmt(format_args!("Counter key not set"));
                };
            }
        }
        pub fn set_start(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(SequenceAttributes::Start { value });
            self
        }
        pub fn set_stop(&mut self, value: RangeValue) -> &mut Self {
            self.insert_to(SequenceAttributes::Stop { value });
            self
        }
        pub fn set_sequence_ffn(&mut self, value: String) -> &mut Self {
            self.insert_to(SequenceAttributes::SequenceFfn {
                value,
            });
            self
        }
        pub fn set_sequence_faa(&mut self, value: String) -> &mut Self {
            self.insert_to(SequenceAttributes::SequenceFaa {
                value,
            });
            self
        }
        pub fn set_codon_start(&mut self, value: u8) -> &mut Self {
            self.insert_to(SequenceAttributes::CodonStart {
                value,
            });
            self
        }
        pub fn set_strand(&mut self, value: i8) -> &mut Self {
            self.insert_to(SequenceAttributes::Strand {
                value,
            });
            self
        }
        pub fn build(self) -> BTreeMap<String, HashSet<SequenceAttributes>> {
            self.seq_attributes
        }
        pub fn iter_sorted(
            &self,
        ) -> std::collections::btree_map::Iter<String, HashSet<SequenceAttributes>> {
            self.seq_attributes.iter()
        }
        pub fn default() -> Self {
            SequenceAttributeBuilder {
                seq_attributes: BTreeMap::new(),
                locus_tag: None,
            }
        }
    }
    ///product lines can contain difficult to parse punctuation such as biochemical symbols like unclosed single quotes, superscripts, single and double brackets etc.
    ///here we substitute these for an underscore
    pub fn substitute_odd_punctuation(input: String) -> Result<String, anyhow::Error> {
        let re = Regex::new(r"[/?()',`]|[α-ωΑ-Ω]")?;
        let cleaned = input.trim_end_matches(&['\r', '\n'][..]);
        Ok(re.replace_all(cleaned, "_").to_string())
    }
    ///GFF3 field9 construct
    pub struct GFFInner {
        pub id: String,
        pub name: String,
        pub locus_tag: String,
        pub gene: String,
        pub product: String,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for GFFInner {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            ::core::fmt::Formatter::debug_struct_field5_finish(
                f,
                "GFFInner",
                "id",
                &self.id,
                "name",
                &self.name,
                "locus_tag",
                &self.locus_tag,
                "gene",
                &self.gene,
                "product",
                &&self.product,
            )
        }
    }
    impl GFFInner {
        pub fn new(
            id: String,
            name: String,
            locus_tag: String,
            gene: String,
            product: String,
        ) -> Self {
            GFFInner {
                id,
                name,
                locus_tag,
                gene,
                product,
            }
        }
    }
    ///The main GFF3 construct
    pub struct GFFOuter<'a> {
        pub seqid: String,
        pub source: String,
        pub type_val: String,
        pub start: u32,
        pub end: u32,
        pub score: f64,
        pub strand: String,
        pub phase: u8,
        pub attributes: &'a GFFInner,
    }
    #[automatically_derived]
    impl<'a> ::core::fmt::Debug for GFFOuter<'a> {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            let names: &'static _ = &[
                "seqid",
                "source",
                "type_val",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ];
            let values: &[&dyn ::core::fmt::Debug] = &[
                &self.seqid,
                &self.source,
                &self.type_val,
                &self.start,
                &self.end,
                &self.score,
                &self.strand,
                &self.phase,
                &&self.attributes,
            ];
            ::core::fmt::Formatter::debug_struct_fields_finish(
                f,
                "GFFOuter",
                names,
                values,
            )
        }
    }
    impl<'a> GFFOuter<'a> {
        pub fn new(
            seqid: String,
            source: String,
            type_val: String,
            start: u32,
            end: u32,
            score: f64,
            strand: String,
            phase: u8,
            attributes: &'a GFFInner,
        ) -> Self {
            GFFOuter {
                seqid,
                source,
                type_val,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            }
        }
        pub fn field9_attributes_build(&self) -> String {
            let mut full_field9 = Vec::new();
            if !self.attributes.id.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("id={0}", self.attributes.id),
                            );
                            res
                        }),
                    );
            }
            if !self.attributes.name.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("name={0}", self.attributes.name),
                            );
                            res
                        }),
                    );
            }
            if !self.attributes.gene.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("gene={0}", self.attributes.gene),
                            );
                            res
                        }),
                    );
            }
            if !self.attributes.locus_tag.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("locus_tag={0}", self.attributes.locus_tag),
                            );
                            res
                        }),
                    );
            }
            if !self.attributes.product.is_empty() {
                full_field9
                    .push(
                        ::alloc::__export::must_use({
                            let res = ::alloc::fmt::format(
                                format_args!("product={0}", self.attributes.product),
                            );
                            res
                        }),
                    );
            }
            full_field9.join(";")
        }
    }
    ///formats the translation string which can be multiple lines, for gbk
    pub fn format_translation(translation: &str) -> String {
        let mut formatted = String::new();
        let cleaned_translation = translation.replace("\n", "");
        formatted.push_str("                     /translation=\"");
        let line_length: usize = 60;
        let final_num = line_length - 15;
        formatted
            .push_str(
                &::alloc::__export::must_use({
                    let res = ::alloc::fmt::format(
                        format_args!("{0}\n", &cleaned_translation[0..final_num]),
                    );
                    res
                }),
            );
        for i in (47..translation.len()).step_by(60) {
            let end = i + 60 - 1;
            let valid_end = if end >= translation.len() {
                &cleaned_translation.len() - 1
            } else {
                end
            };
            formatted
                .push_str(
                    &::alloc::__export::must_use({
                        let res = ::alloc::fmt::format(
                            format_args!(
                                "                     {0}",
                                &cleaned_translation[i..valid_end],
                            ),
                        );
                        res
                    }),
                );
            {
                ::std::io::_print(
                    format_args!(
                        "cleaned translation leng is {0:?}\n",
                        &cleaned_translation[i..valid_end].len(),
                    ),
                );
            };
            if *&cleaned_translation[i..valid_end].len() < 59 {
                formatted.push('\"');
            } else {
                formatted.push('\n');
            }
        }
        formatted
    }
    ///writes the DNA sequence in gbk format with numbering
    pub fn write_gbk_format_sequence(sequence: &str, file: &mut File) -> io::Result<()> {
        file.write_fmt(format_args!("ORIGIN\n"))?;
        let mut formatted = String::new();
        let cleaned_input = sequence.replace("\n", "");
        let mut index = 1;
        for (_i, chunk) in cleaned_input.as_bytes().chunks(60).enumerate() {
            formatted
                .push_str(
                    &::alloc::__export::must_use({
                        let res = ::alloc::fmt::format(format_args!("{0:>5} ", index));
                        res
                    }),
                );
            for (j, sub_chunk) in chunk.chunks(10).enumerate() {
                if j > 0 {
                    formatted.push(' ');
                }
                formatted.push_str(&String::from_utf8_lossy(sub_chunk));
            }
            formatted.push('\n');
            index += 60;
        }
        file.write_fmt(format_args!("{0:>6}\n", &formatted))?;
        file.write_fmt(format_args!("//\n"))?;
        Ok(())
    }
    ///saves the parsed data in genbank format
    pub fn gbk_write(
        seq_region: BTreeMap<String, (u32, u32)>,
        record_vec: Vec<Record>,
        filename: &str,
    ) -> io::Result<()> {
        let now = Local::now();
        let formatted_date = now.format("%d-%b-%Y").to_string().to_uppercase();
        let mut file = OpenOptions::new()
            .write(true)
            .append(true)
            .create(true)
            .open(filename)?;
        for (i, (key, _val)) in seq_region.iter().enumerate() {
            let strain = match &record_vec[i].source_map.get_strain(key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let organism = match &record_vec[i].source_map.get_organism(key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let mol_type = match &record_vec[i].source_map.get_mol_type(key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let type_material = match &record_vec[i].source_map.get_type_material(&key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let db_xref = match &record_vec[i].source_map.get_db_xref(key) {
                Some(value) => value.to_string(),
                None => "Unknown".to_string(),
            };
            let source_stop = match &record_vec[i].source_map.get_stop(key) {
                Some(value) => value.get_value(),
                None => {
                    {
                        {
                            ::std::io::_print(format_args!("stop value not found\n"));
                        };
                        None
                    }
                        .expect("stop value not received")
                }
            };
            file.write_fmt(
                format_args!(
                    "LOCUS       {0}             {1} bp    DNA     linear CON {2}\n",
                    &key,
                    &record_vec[i].sequence.len(),
                    &formatted_date,
                ),
            )?;
            file.write_fmt(format_args!("DEFINITION  {0} {1}.\n", &organism, &strain))?;
            file.write_fmt(format_args!("ACCESSION   {0}\n", &key))?;
            file.write_fmt(format_args!("KEYWORDS    .\n"))?;
            file.write_fmt(format_args!("SOURCE      {0} {1}\n", &organism, &strain))?;
            file.write_fmt(format_args!("  ORGANISM  {0} {1}\n", &organism, &strain))?;
            file.write_fmt(format_args!("FEATURES             Location/Qualifiers\n"))?;
            file.write_fmt(format_args!("     source          1..{0}\n", &source_stop))?;
            file.write_fmt(
                format_args!("                     /organism=\"{0}\"\n", &strain),
            )?;
            file.write_fmt(
                format_args!("                     /mol_type=\"{0}\"\n", &mol_type),
            )?;
            file.write_fmt(
                format_args!("                     /strain=\"{0}\"\n", &strain),
            )?;
            if type_material != *"Unknown".to_string() {
                file.write_fmt(
                    format_args!(
                        "                     /type_material=\"{0}\"\n",
                        &type_material,
                    ),
                )?;
            }
            file.write_fmt(
                format_args!("                     /db_xref=\"{0}\"\n", &db_xref),
            )?;
            for (locus_tag, _value) in &record_vec[i].cds.attributes {
                let start = match &record_vec[i].cds.get_start(locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("start value not found\n"));
                            };
                            None
                        }
                            .expect("start value not received")
                    }
                };
                let stop = match &record_vec[i].cds.get_stop(locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("stop value not found\n"));
                            };
                            None
                        }
                            .expect("stop value not received")
                    }
                };
                let product = match &record_vec[i].cds.get_product(locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown product".to_string(),
                };
                let strand = match &record_vec[i].cds.get_strand(locus_tag) {
                    Some(value) => **value,
                    None => 0,
                };
                let codon_start = match &record_vec[i].cds.get_codon_start(locus_tag) {
                    Some(value) => **value,
                    None => 0,
                };
                let gene = match &record_vec[i].cds.get_gene(locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown".to_string(),
                };
                let translation = match &record_vec[i]
                    .seq_features
                    .get_sequence_faa(locus_tag)
                {
                    Some(value) => value.to_string(),
                    None => "unknown".to_string(),
                };
                if strand == 1 {
                    file.write_fmt(
                        format_args!("     gene            {0}..{1}\n", &start, &stop),
                    )?;
                } else {
                    file.write_fmt(
                        format_args!(
                            "     gene            complement({0}..{1})\n",
                            &start,
                            &stop,
                        ),
                    )?;
                }
                file.write_fmt(
                    format_args!("                     /locus_tag=\"{0}\"\n", &locus_tag),
                )?;
                if strand == 1 {
                    file.write_fmt(
                        format_args!("     CDS             {0}..{1}\n", &start, &stop),
                    )?;
                } else {
                    file.write_fmt(
                        format_args!(
                            "     CDS             complement({0}..{1})\n",
                            &start,
                            &stop,
                        ),
                    )?;
                }
                file.write_fmt(
                    format_args!("                     /locus_tag=\"{0}\"\n", &locus_tag),
                )?;
                file.write_fmt(
                    format_args!(
                        "                     /codon_start=\"{0}\"\n",
                        &codon_start,
                    ),
                )?;
                if gene != "unknown" {
                    file.write_fmt(
                        format_args!("                     /gene=\"{0}\"\n", &gene),
                    )?;
                }
                if translation != "unknown" {
                    let formatted_translation = format_translation(&translation);
                    file.write_fmt(format_args!("{0}\n", &formatted_translation))?;
                }
                file.write_fmt(
                    format_args!("                     /product=\"{0}\"\n", &product),
                )?;
            }
            write_gbk_format_sequence(&record_vec[i].sequence, &mut file)?;
        }
        Ok(())
    }
    ///saves the parsed data in gff3 format
    #[allow(unused_assignments)]
    #[allow(unused_variables)]
    pub fn gff_write(
        seq_region: BTreeMap<String, (u32, u32)>,
        mut record_vec: Vec<Record>,
        filename: &str,
        dna: bool,
    ) -> io::Result<()> {
        let mut file = OpenOptions::new().append(true).create(true).open(filename)?;
        if file.metadata()?.len() == 0 {
            file.write_fmt(format_args!("##gff-version 3\n"))?;
        }
        let mut full_seq = String::new();
        let mut prev_end: u32 = 0;
        for (k, v) in seq_region.iter() {
            file.write_fmt(
                format_args!("##sequence-region\t{0}\t{1}\t{2}\n", &k, v.0, v.1),
            )?;
        }
        for ((source_name, (seq_start, seq_end)), record) in seq_region
            .iter()
            .zip(record_vec.drain(..))
        {
            if dna == true {
                full_seq.push_str(&record.sequence);
            }
            for (locus_tag, _valu) in &record.cds.attributes {
                let start = match record.cds.get_start(locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("start value not found\n"));
                            };
                            None
                        }
                            .expect("start value not received")
                    }
                };
                let stop = match record.cds.get_stop(locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("stop value not found\n"));
                            };
                            None
                        }
                            .expect("stop value not received")
                    }
                };
                let gene = match record.cds.get_gene(locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown".to_string(),
                };
                let product = match record.cds.get_product(locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown product".to_string(),
                };
                let strand = match record.cds.get_strand(locus_tag) {
                    Some(valu) => {
                        match valu {
                            1 => "+".to_string(),
                            -1 => "-".to_string(),
                            _ => {
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "unexpected strand value {0} for locus_tag {1}\n",
                                            valu,
                                            locus_tag,
                                        ),
                                    );
                                };
                                "unknownstrand".to_string()
                            }
                        }
                    }
                    None => "unknownvalue".to_string(),
                };
                let phase = match record.cds.get_codon_start(locus_tag) {
                    Some(valuer) => {
                        match valuer {
                            1 => 0,
                            2 => 1,
                            3 => 2,
                            _ => {
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "unexpected phase value {0} in the bagging area for locus_tag {1}\n",
                                            valuer,
                                            locus_tag,
                                        ),
                                    );
                                };
                                1
                            }
                        }
                    }
                    None => 1,
                };
                let gff_inner = GFFInner::new(
                    locus_tag.to_string(),
                    source_name.clone(),
                    locus_tag.to_string(),
                    gene,
                    product,
                );
                let gff_outer = GFFOuter::new(
                    source_name.clone(),
                    ".".to_string(),
                    "CDS".to_string(),
                    start + prev_end,
                    stop + prev_end,
                    0.0,
                    strand,
                    phase,
                    &gff_inner,
                );
                let field9_attributes = gff_outer.field9_attributes_build();
                file.write_fmt(
                    format_args!(
                        "{0}\t{1}\t{2}\t{3:?}\t{4:?}\t{5}\t{6}\t{7}\t{8}\n",
                        gff_outer.seqid,
                        gff_outer.source,
                        gff_outer.type_val,
                        gff_outer.start,
                        gff_outer.end,
                        gff_outer.score,
                        gff_outer.strand,
                        gff_outer.phase,
                        field9_attributes,
                    ),
                )?;
            }
            prev_end = *seq_end;
        }
        if dna {
            file.write_fmt(format_args!("##FASTA\n"))?;
            file.write_fmt(format_args!("{0}\n", full_seq))?;
        }
        Ok(())
    }
    ///saves the parsed data in gff3 format
    #[allow(unused_assignments)]
    pub fn orig_gff_write(
        seq_region: BTreeMap<String, (u32, u32)>,
        record_vec: Vec<Record>,
        filename: &str,
        dna: bool,
    ) -> io::Result<()> {
        let mut file = OpenOptions::new().append(true).create(true).open(filename)?;
        if file.metadata()?.len() == 0 {
            file.write_fmt(format_args!("##gff-version 3\n"))?;
        }
        let mut source_name = String::new();
        let mut full_seq = String::new();
        let mut prev_end: u32 = 0;
        for (k, v) in seq_region.iter() {
            file.write_fmt(
                format_args!("##sequence-region\t{0}\t{1}\t{2}\n", &k, v.0, v.1),
            )?;
        }
        for (i, (key, val)) in seq_region.iter().enumerate() {
            source_name = key.to_string();
            if dna == true {
                full_seq.push_str(&record_vec[i].sequence);
            }
            for (locus_tag, _valu) in &record_vec[i].cds.attributes {
                let start = match record_vec[i].cds.get_start(locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("start value not found\n"));
                            };
                            None
                        }
                            .expect("start value not received")
                    }
                };
                let stop = match record_vec[i].cds.get_stop(locus_tag) {
                    Some(value) => value.get_value(),
                    None => {
                        {
                            {
                                ::std::io::_print(format_args!("stop value not found\n"));
                            };
                            None
                        }
                            .expect("stop value not received")
                    }
                };
                let gene = match record_vec[i].cds.get_gene(locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown".to_string(),
                };
                let product = match record_vec[i].cds.get_product(locus_tag) {
                    Some(value) => value.to_string(),
                    None => "unknown product".to_string(),
                };
                let strand = match record_vec[i].cds.get_strand(locus_tag) {
                    Some(valu) => {
                        match valu {
                            1 => "+".to_string(),
                            -1 => "-".to_string(),
                            _ => {
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "unexpected strand value {0} for locus_tag {1}\n",
                                            valu,
                                            locus_tag,
                                        ),
                                    );
                                };
                                "unknownstrand".to_string()
                            }
                        }
                    }
                    None => "unknownvalue".to_string(),
                };
                let phase = match record_vec[i].cds.get_codon_start(locus_tag) {
                    Some(valuer) => {
                        match valuer {
                            1 => 0,
                            2 => 1,
                            3 => 2,
                            _ => {
                                {
                                    ::std::io::_print(
                                        format_args!(
                                            "unexpected phase value {0} in the bagging area for locus_tag {1}\n",
                                            valuer,
                                            locus_tag,
                                        ),
                                    );
                                };
                                1
                            }
                        }
                    }
                    None => 1,
                };
                let gff_inner = GFFInner::new(
                    locus_tag.to_string(),
                    source_name.clone(),
                    locus_tag.to_string(),
                    gene,
                    product,
                );
                let gff_outer = GFFOuter::new(
                    source_name.clone(),
                    ".".to_string(),
                    "CDS".to_string(),
                    start + prev_end,
                    stop + prev_end,
                    0.0,
                    strand,
                    phase,
                    &gff_inner,
                );
                let field9_attributes = gff_outer.field9_attributes_build();
                file.write_fmt(
                    format_args!(
                        "{0}\t{1}\t{2}\t{3:?}\t{4:?}\t{5}\t{6}\t{7}\t{8}\n",
                        gff_outer.seqid,
                        gff_outer.source,
                        gff_outer.type_val,
                        gff_outer.start,
                        gff_outer.end,
                        gff_outer.score,
                        gff_outer.strand,
                        gff_outer.phase,
                        field9_attributes,
                    ),
                )?;
            }
            prev_end = val.1;
        }
        if dna {
            file.write_fmt(format_args!("##FASTA\n"))?;
            file.write_fmt(format_args!("{0}\n", full_seq))?;
        }
        Ok(())
    }
    ///internal record containing data from a single source or contig.  Has multiple features.
    pub struct Record {
        pub id: String,
        pub length: u32,
        pub sequence: String,
        pub start: usize,
        pub end: usize,
        pub strand: i32,
        pub cds: FeatureAttributeBuilder,
        pub source_map: SourceAttributeBuilder,
        pub seq_features: SequenceAttributeBuilder,
    }
    #[automatically_derived]
    impl ::core::fmt::Debug for Record {
        #[inline]
        fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
            let names: &'static _ = &[
                "id",
                "length",
                "sequence",
                "start",
                "end",
                "strand",
                "cds",
                "source_map",
                "seq_features",
            ];
            let values: &[&dyn ::core::fmt::Debug] = &[
                &self.id,
                &self.length,
                &self.sequence,
                &self.start,
                &self.end,
                &self.strand,
                &self.cds,
                &self.source_map,
                &&self.seq_features,
            ];
            ::core::fmt::Formatter::debug_struct_fields_finish(
                f,
                "Record",
                names,
                values,
            )
        }
    }
    #[automatically_derived]
    impl ::core::clone::Clone for Record {
        #[inline]
        fn clone(&self) -> Record {
            Record {
                id: ::core::clone::Clone::clone(&self.id),
                length: ::core::clone::Clone::clone(&self.length),
                sequence: ::core::clone::Clone::clone(&self.sequence),
                start: ::core::clone::Clone::clone(&self.start),
                end: ::core::clone::Clone::clone(&self.end),
                strand: ::core::clone::Clone::clone(&self.strand),
                cds: ::core::clone::Clone::clone(&self.cds),
                source_map: ::core::clone::Clone::clone(&self.source_map),
                seq_features: ::core::clone::Clone::clone(&self.seq_features),
            }
        }
    }
    impl Record {
        /// Create a new instance.
        pub fn new() -> Self {
            Record {
                id: "".to_owned(),
                length: 0,
                sequence: "".to_owned(),
                start: 0,
                end: 0,
                strand: 0,
                source_map: SourceAttributeBuilder::new(),
                cds: FeatureAttributeBuilder::new(),
                seq_features: SequenceAttributeBuilder::new(),
            }
        }
        pub fn is_empty(&mut self) -> bool {
            self.id.is_empty() && self.length == 0
        }
        pub fn check(&mut self) -> Result<(), &str> {
            if self.id().is_empty() {
                return Err("Expecting id for Gbk record.");
            }
            Ok(())
        }
        pub fn id(&mut self) -> &str {
            &self.id
        }
        pub fn length(&mut self) -> u32 {
            self.length
        }
        pub fn sequence(&mut self) -> &str {
            &self.sequence
        }
        pub fn start(&mut self) -> u32 {
            self.start.try_into().unwrap()
        }
        pub fn end(&mut self) -> u32 {
            self.end.try_into().unwrap()
        }
        pub fn strand(&mut self) -> i32 {
            self.strand
        }
        pub fn cds(&mut self) -> FeatureAttributeBuilder {
            self.cds.clone()
        }
        pub fn source_map(&mut self) -> SourceAttributeBuilder {
            self.source_map.clone()
        }
        pub fn seq_features(&mut self) -> SequenceAttributeBuilder {
            self.seq_features.clone()
        }
        fn rec_clear(&mut self) {
            self.id.clear();
            self.length = 0;
            self.sequence.clear();
            self.start = 0;
            self.end = 0;
            self.strand = 0;
            self.source_map = SourceAttributeBuilder::new();
            self.cds = FeatureAttributeBuilder::new();
            self.seq_features = SequenceAttributeBuilder::new();
        }
    }
    impl Default for Record {
        fn default() -> Self {
            Self::new()
        }
    }
    #[allow(dead_code)]
    pub struct Config {
        filename: String,
    }
    impl Config {
        pub fn new(args: &[String]) -> Result<Config, &str> {
            if args.len() < 2 {
                {
                    ::core::panicking::panic_fmt(
                        format_args!("not enough arguments, please provide filename"),
                    );
                };
            }
            let filename = args[1].clone();
            Ok(Config { filename })
        }
    }
}
