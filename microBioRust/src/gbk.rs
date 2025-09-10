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
//! use std::{
//!     fs::File,
//!     io,
//! };
//! use microBioRust::gbk::Reader;
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
//! use std::{
//!     fs::File,
//!     io,
//! };
//! use microBioRust::{
//!     gbk::Reader,
//!     genbank,
//! };
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
//!    use microBioRust::gbk::{gff_write, Reader, Record};
//!    use std::collections::BTreeMap;
//!    use std::{
//!        fs::File,
//!        io,
//!    };
//!    use clap::Parser;
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
//!    use microBioRust::gbk::{gff_write, RangeValue, Record};
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

use std::{
	io::{self, Write},
	fs::{self, OpenOptions, File},
	vec::Vec,
	str,
	convert::{AsRef, TryInto},
	path::Path,
	collections::{BTreeMap, HashSet},
};
use regex::Regex;
use itertools::Itertools;
use protein_translate::translate;
use bio::alphabets::dna::revcomp;
use anyhow::{anyhow, Context};
use paste::paste;
use chrono::prelude::*;


/// macro to create get_ functions for the values
#[macro_export]
macro_rules! create_getters {
    // macro for creating get methods
    ($struct_name:ident, $attributes:ident, $enum_name:ident, $( $field:ident { value: $type:ty } ),* ) => {
		impl $struct_name {
            $(
	        // creates a get method for each of the fields in the SourceAttributes, FeatureAttributes and SequenceAttributes
	        paste! {
                  pub fn [<get_$field:snake>](&self, key: &str) -> Option<&$type> {
                    // Get the HashSet for the key (e.g., "source_1")
                    self.$attributes.get(key).and_then(|set| {
                        // Iterate over the HashSet to find the correct SourceAttributes value
                        set.iter().find_map(|attr| {
                            if let $enum_name::$field { value } = attr {
                                Some(value)
                            } else {
                                None
                            }
                        })
                    })
                }
	      }
            )*
        }
    };
}

/// macro to create the set_ functions for the values in a Builder format
#[macro_export] 
macro_rules! create_builder {
    // Macro for creating attribute builders for SourceAttributes, FeatureAttributes and SequenceAttributes
    ($builder_name:ident, $attributes:ident, $enum_name:ident, $counter_name:ident, $( $field:ident { value: $type:ty } ),* ) => {
        impl $builder_name {
            pub fn new() -> Self {
                $builder_name {
                    $attributes: BTreeMap::new(),
                    $counter_name: None,
                }
            }
            //sets the key for the BTreeMap 
            pub fn set_counter(&mut self, counter: String) -> &mut Self {
                self.$counter_name = Some(counter);
		self
            }    
            //function to insert the fields from the enum into the attributes
            pub fn insert_to(&mut self, value: $enum_name) {
	        if let Some(counter) = &self.$counter_name {
		    self.$attributes
		        .entry(counter.to_string())
                        .or_insert_with(HashSet::new)
                        .insert(value);
		    }
		else {
		    panic!("Counter key not set"); // Needs better error handling
		    }
            }
            // function to set each of the alternative fields in the builder
            $(
	      paste! { 
	        pub fn [<set_$field:snake>](&mut self, value: $type) -> &mut Self {
	           self.insert_to($enum_name::$field { value });
		   self
	           }
		}
	    )*
	    // build function to the attributes
	    pub fn build(self) -> BTreeMap<String, HashSet<$enum_name>> {
	        self.$attributes
            }
	    // function to iterate immutably through the BTreeMap as required
	    pub fn iter_sorted(&'_ self) -> std::collections::btree_map::Iter<String, HashSet<$enum_name>> {
	        self.$attributes.iter()
	    }
	    //default function
	    pub fn default() -> Self {
	        $builder_name {
		    $attributes: BTreeMap::new(),
		    $counter_name: None,
		    }
		}
            }
     };
}

#[macro_export]
macro_rules! genbank {
    ($filename:expr) => {{
        use std::fs::File;
        use std::io::BufReader;
        let file = File::open($filename)
            .unwrap_or_else(|e| panic!("Could not open file {}: {}", $filename, e));
        let mut reader = $crate::gbk::Reader::new(file);
        let mut vec = Vec::new();
        for rec in reader.records() {
            match rec {
                Ok(r) => { //println!("this is r {:?}", &r);
		           vec.push(r);
			   }
                Err(e) => panic!("Error reading record: {:?}", e),
            }
        }
        vec
    }};
}


//const MAX_GBK_BUFFER_SIZE: usize = 512;
/// A Gbk reader.

#[derive(Debug)]
#[allow(unused_mut)]
pub struct Records<B>
where
    B: io::BufRead,
{
    reader: Reader<B>,
    error_has_occurred: bool,
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
	        println!("error was encountered in iteration");
		None
	        } else {
                let mut record = Record::new();
                match self.reader.read(&mut record) {
	            Ok(_) => { if record.is_empty() {
		       None }
		       else {
		           Some(Ok(record))
			   }
		        }
		    Err(err) => {
		        //println!("we encountered an error {:?}", &err);
		        self.error_has_occurred = true;
		        Some(Err(anyhow!("next record read error {:?}",err)))
		        }
                    }
                }
     }
}

pub trait GbkRead {
    fn read(&mut self, record: &mut Record) -> Result<Record, anyhow::Error>;
}

///per line reader for the file 
#[derive(Debug, Default)]
pub struct Reader<B> {
    reader: B,
    line_buffer: String,
}

impl Reader<io::BufReader<fs::File>> {
    /// Read Gbk from given file path in given format.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .with_context(|| format!("Failed to read Gbk from {:#?}", path))
    }
}

impl<R> Reader<io::BufReader<R>>
where
     R: io::Read,
{
    //// Create a new Gbk reader given an instance of `io::Read` in given format
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
    //return an iterator over the records of the genbank file
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
	//println!("reading new record");
	//initialise variables
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
	//check if there are any more lines, if not return the record as is
    	if self.line_buffer.is_empty() {
	    self.reader.read_line(&mut self.line_buffer)?;
	    if self.line_buffer.is_empty() {
	        return Ok(record.to_owned());
	        }
            }
	//main loop to populate the attributes and iterate through the file
	'outer: while !self.line_buffer.is_empty() {
	    //println!("is line buffer {:?}", &self.line_buffer);
	    //collect the header fields
	    if self.line_buffer.starts_with("LOCUS") {
			record.rec_clear();
	            	let mut header_fields: Vec<&str> = self.line_buffer.split_whitespace().collect();
	                let mut header_iter = header_fields.iter();
	                header_iter.next();
	                record.id = header_iter.next()
						.ok_or_else(|| anyhow::anyhow!("missing record id"))? // Get &str or error
						.to_string();
	                let lens = header_iter.next()
						.ok_or_else(|| anyhow::anyhow!("missing record length"))? // Get &str or error
						.to_string();
	                record.length = lens.trim().parse::<u32>()?;
			self.line_buffer.clear();
			}
	    //collect the source fields and populate the source_map and source_attributes
	    if self.line_buffer.starts_with("     source") {
	        let re = Regex::new(r"([0-9]+)[[:punct:]]+([0-9]+)")?;
		let location = re.captures(&self.line_buffer).ok_or_else(|| anyhow::anyhow!("missing location"))?;
		let start = &location[1];
		let end = &location[2];
		thestart = start.trim().parse::<u32>()?;
		source_counter+=1;
		source_name = format!("source_{}_{}",record.id,source_counter).to_string();
		thestart += prev_end;
		theend = end.trim().parse::<u32>()? + prev_end;
		//println!("so the start and end are {:?} {:?}", &thestart, &theend);
		loop {
		    self.line_buffer.clear();
		    self.reader.read_line(&mut self.line_buffer)?;
		    if self.line_buffer.starts_with("     CDS") {
		            //println!("this source name {:?} start {:?} end {:?} organism {:?} mol_type {:?} strain {:?} type_material {:?} db_xref {:?}", &source_name,&thestart, &theend, &organism, &mol_type, &strain, &type_material, &db_xref);
		            record.source_map
			       .set_counter(source_name.to_string())
			       .set_start(RangeValue::Exact(thestart))
			       .set_stop(RangeValue::Exact(theend))
			       .set_organism(organism.clone())
			       .set_mol_type(mol_type.clone())
			       .set_strain(strain.clone())
			      // culture_collection.clone()
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
		//    if self.line_buffer.contains("/culture_collection") {
		//        let cc: Vec<&str> = self.line_buffer.split('\"').collect();
	//		culture_collection = cc[1].to_string();
	//		}
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
	    //populate the FeatureAttributes and the coding sequence annotation
	    if self.line_buffer.starts_with("     CDS") {
	        let mut startiter: Vec<_> = Vec::new();
		let mut enditer: Vec<_> = Vec::new();
		let mut thestart: u32 = 0;
		let mut thend: u32 = 0;
		let mut joined: bool = false;
	        //gather the feature coordinates
		let joined = if self.line_buffer.contains("join") { true } else { false };
	        let re = Regex::new(r"([0-9]+)[[:punct:]]+([0-9]+)")?;
		//let matches: Vec<&regex::Captures> = re.captures_iter(&self.line_buffer).collect();
		for cap in re.captures_iter(&self.line_buffer) {
		   cds_counter+=1;
		   thestart = cap[1].parse().expect("failed to match and parse numerical start");
		   theend = cap[2].parse().expect("failed to match and parse numerical end");
		   startiter.push(thestart);
		   enditer.push(theend);
		   }
		let mut gene = String::new();
		let mut product = String::new();
		let strand: i8 = if self.line_buffer.contains("complement") {-1} else {1};
                let mut locus_tag = String::new();
                let mut codon_start: u8 = 1;
		//loop to populate the feature attributes, when complete it calls to the outer loop directly to prevent reading a new line into self.line_buffer
		loop {
		        self.line_buffer.clear();
			self.reader.read_line(&mut self.line_buffer)?;
                        if self.line_buffer.contains("/locus_tag=") {
                            let loctag: Vec<&str> = self.line_buffer.split('\"').collect();
                            locus_tag = loctag[1].to_string();
			    //println!("designated locus tag {:?}", &locus_tag);
                            }
                        if self.line_buffer.contains("/codon_start") {
                            let codstart: Vec<&str> = self.line_buffer.split('=').collect();
                            let valstart = codstart[1].trim().parse::<u8>()?;
                            codon_start = valstart;
			    //println!("designated codon start {:?} {:?}", &codon_start, &locus_tag);
                            }
                        if self.line_buffer.contains("/gene=") {
		            let gen: Vec<&str> = self.line_buffer.split('\"').collect();
			    gene = gen[1].to_string();
			    //println!("gene designated {:?} {:?}", &gene, &locus_tag);
			    }
			if self.line_buffer.contains("/product") {
		            let prod: Vec<&str> = self.line_buffer.split('\"').collect();
			    product = substitute_odd_punctuation(prod[1].to_string())?;
			    //println!("designated product {:?} {:?}", &product, &locus_tag);
			    }
			if self.line_buffer.starts_with("     CDS") || self.line_buffer.starts_with("ORIGIN") || self.line_buffer.starts_with("     gene") || self.line_buffer.starts_with("     misc_feature") {
			    if locus_tag.is_empty() {
			         locus_tag = format!("CDS_{}",cds_counter).to_string();
				 }
			    if joined {
				 //println!("currently the start is {:?} and the stop is {:?}", &startiter, &enditer);
				 for (i, m) in startiter.iter().enumerate() {
				      let loc_tag = format!("{}_{}",locus_tag.clone(),i);
                                      //check we may need to add or subtract one to m
				      record.cds
				          .set_counter(loc_tag)
					  .set_start(RangeValue::Exact(*m))
					  .set_stop(RangeValue::Exact(enditer[i]))
					  .set_gene(gene.to_string())
					  .set_product(product.to_string())
					  .set_codon_start(codon_start)
					  .set_strand(strand);
				      }
				 continue 'outer;
			         }
			    else {
			         record.cds
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
	     }     }
	    //check if we have reached the DNA sequence section and populate the record sequences field if so.  Returns the record on finding end of record mark
	    if self.line_buffer.starts_with("ORIGIN") {
	        let mut sequences = String::new();
	        let result_seq = loop {  
		     self.line_buffer.clear();
		     self.reader.read_line(&mut self.line_buffer)?;
                     if self.line_buffer.starts_with("//") {
		         break sequences;
                     } else {
	                 let s: Vec<&str> = self.line_buffer.split_whitespace().collect();
		         let s = &s[1..];
		         let sequence = s.iter().join("");
		         sequences.push_str(&sequence);
                         }     
	             };
		record.sequence = result_seq.to_string();
		let mut iterablecount: u32 = 0;
		//Fields are completed and populated for the FeatureAttributes, collect and populate the SequenceAttributes fields
	        for (key,val) in record.cds.iter_sorted() {
	              let (mut a, mut b, mut c, mut d): (Option<u32>, Option<u32>, Option<i8>, Option<u8>) = (None, None, None, None);
	              for value in val {
		          //println!("this is key {:?} value {:?}", &key, &value);
	                  match value {
		               FeatureAttributes::Start { value } => a = match value {
		                   RangeValue::Exact(v) => Some(*v),
                                   RangeValue::LessThan(v) => Some(*v), // Assign the value even if it's <value
                                   RangeValue::GreaterThan(v) => Some(*v), //Assign the value even it's > value
                                   },
		               FeatureAttributes::Stop { value } => b = match value {
		                   RangeValue::Exact(v) => Some(*v),
                                   RangeValue::LessThan(v) => Some(*v), // Assign the value even if it's <value
                                   RangeValue::GreaterThan(v) => Some(*v), //Assign the value even if it's > value
                                   },
		               FeatureAttributes::Strand { value } => c = match value {
		                  value => Some(*value),
			          },
		               FeatureAttributes::CodonStart { value } => d = match value {
		                  value => Some(value.clone()),
			          },
		               _ => (),
		               }
	                 }
	                 let sta = a.map(|o| o as usize)
						 .ok_or(anyhow!("No value for start"))?;
	                 let sto = b.map(|t| t as usize)
						 .ok_or(anyhow!("No value for stop"))? - 1;
	                 let stra = c.map(|u| u as i8)
						 .ok_or(anyhow!("No value for strand"))?;
	                 let cod = d.map(|v| v as usize - 1)
						 .ok_or(anyhow!("No value for strand"))?;

	                 let star = sta.try_into()?;
	                 let stow = sto.try_into()?;
	                 let codd = cod.try_into()?;
	                 let mut sliced_sequence: &str = "";
	                 //collects the DNA sequence and translations on the correct strand
	                 if stra == -1 {
	                    if cod > 1 {
			       //println!("reverse strand coding start more than one {:?}", &iterablecount);
			       if sto + 1 <= record.sequence.len() {
		                   sliced_sequence = &record.sequence[sta+cod..sto+1];
				   }
				   else {
				   sliced_sequence = &record.sequence[sta+cod..sto];
				   }
		               }
		            else {
			       //println!("record sta {:?} sto {:?} cod {:?} stra {:?} record.seq length {:?}", &sta, &sto, &cod, &stra, &record.sequence.len());
			       //println!("sliced sta {:?} sliced sto {:?} record.id {:?}", sta, sto, &record.id);
			       //println!("iterable count is {:?} reverse strand codon start one", &iterablecount);
			       if sto + 1 <= record.sequence.len() {
	                                        sliced_sequence = &record.sequence[sta..sto+1];
				                  }
			       else {
			           sliced_sequence = &record.sequence[sta..sto];
				   }
			       println!("iterable count after is {:?}", &iterablecount);
		               }
	                 let cds_char = sliced_sequence;
		         let prot_seq =  translate(&revcomp(cds_char.as_bytes()));
		         let parts: Vec<&str> = prot_seq.split('*').collect();
		         record.seq_features
			             .set_counter(key.to_string())
				     .set_start(RangeValue::Exact(star))
                                     .set_stop(RangeValue::Exact(stow))
                                     .set_sequence_ffn(cds_char.to_string())
                                     .set_sequence_faa(parts[0].to_string())
                                     .set_codon_start(codd)
                                     .set_strand(stra);
	                 } else {
	                      if cod > 1 {
			          //println!("forward strand codon value more than one cnt {:?}", &iterablecount);
		                  sliced_sequence = &record.sequence[sta+cod-1..sto];
		                  }
		              else {
			          println!("forward strand codon value one cnt {:?}", &iterablecount);
		                  sliced_sequence = &record.sequence[sta-1..sto];
		                  }
		         let cds_char = sliced_sequence;
		         let prot_seq = translate(cds_char.as_bytes());
		         let parts: Vec<&str> = prot_seq.split('*').collect();
		         record.seq_features
			             .set_counter(key.to_string())
				     .set_start(RangeValue::Exact(star))
                                     .set_stop(RangeValue::Exact(stow))
                                     .set_sequence_ffn(cds_char.to_string())
                                     .set_sequence_faa(parts[0].to_string())
                                     .set_codon_start(codd)
                                     .set_strand(stra);
		                      }
                              }
                 	//return the record when completed
		        return Ok(record.to_owned());
                        }
	 //clear the line buffer and read the next to continue back to the outer loop
	 self.line_buffer.clear();
	 self.reader.read_line(&mut self.line_buffer)?;
        }
       Ok(record.to_owned())
     }
}

pub use crate::record::RangeValue;

//stores the details of the source features in genbank (contigs)
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub enum SourceAttributes {
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Organism { value: String },
    MolType { value: String},
    Strain { value: String},
    CultureCollection { value: String},
    TypeMaterial { value: String},
    DbXref { value:String}
}

//macro for creating the getters
create_getters!(
    SourceAttributeBuilder,
    source_attributes,
    SourceAttributes,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Organism { value: String },
    MolType { value: String},
    Strain { value: String},
    // CultureCollection { value: String},
    TypeMaterial { value: String},
    DbXref { value:String}
);

///builder for the source information on a per record basis
#[derive(Debug, Default, Clone)]
pub struct SourceAttributeBuilder {
    pub source_attributes: BTreeMap<String, HashSet<SourceAttributes>>,
    pub source_name: Option<String>,
}

impl SourceAttributeBuilder {
    // Method to set source name
    pub fn set_source_name(&mut self, name: String) {
        self.source_name = Some(name);
    }

    // Method to get source name
    pub fn get_source_name(&self) -> Option<&String> {
        self.source_name.as_ref()
    }

    // Method to add source attributes
    pub fn add_source_attribute(&mut self, key: String, attribute: SourceAttributes) {
        self.source_attributes
            .entry(key)
            .or_insert_with(HashSet::new)
            .insert(attribute);
    }

    // Method to retrieve source attributes for a given key
    pub fn get_source_attributes(&self, key: &str) -> Option<&HashSet<SourceAttributes>> {
        self.source_attributes.get(key)
    }
}


create_builder!(
    SourceAttributeBuilder,
    source_attributes,
    SourceAttributes,
    source_name,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Organism { value: String },
    MolType { value: String},
    Strain { value: String},
    // CultureCollection { value: String},
    TypeMaterial { value: String},
    DbXref { value:String}
);

///attributes for each feature, cds or gene
#[derive(Debug, Eq, Hash, PartialEq, Clone)]
pub enum FeatureAttributes {
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Gene { value: String },
    Product { value: String },
    CodonStart { value: u8 },
    Strand { value: i8 },
 //   ec_number { value: String }
}


create_getters!(
    FeatureAttributeBuilder,
    attributes,
    FeatureAttributes,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Gene { value: String },
    Product { value: String },
    CodonStart { value: u8 },
    Strand { value: i8 }
);

///builder for the feature information on a per coding sequence (CDS) basis
#[derive(Debug, Default, Clone)]
pub struct FeatureAttributeBuilder {
    pub attributes: BTreeMap<String, HashSet<FeatureAttributes>>,
    locus_tag: Option<String>,
}

create_builder!(
    FeatureAttributeBuilder,
    attributes,
    FeatureAttributes,
    locus_tag,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Gene { value: String },
    Product { value: String },
    CodonStart { value: u8 },
    Strand { value: i8 }
);

///stores the sequences of the coding sequences (genes) and proteins. Also stores start, stop, codon_start and strand information
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub enum SequenceAttributes {
    Start { value: RangeValue },
    Stop { value: RangeValue },
    SequenceFfn { value: String },
    SequenceFaa { value: String },
    CodonStart { value: u8 },
    Strand { value: i8 },
}

create_getters!(
    SequenceAttributeBuilder,
    seq_attributes,
    SequenceAttributes,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    SequenceFfn { value: String},
    SequenceFaa { value: String},
    CodonStart { value: u8},
    Strand { value: i8}
);

///builder for the sequence information on a per coding sequence (CDS) basis
#[derive(Debug, Default, Clone)]
pub struct SequenceAttributeBuilder {
    pub seq_attributes: BTreeMap<String, HashSet<SequenceAttributes>>,
    pub locus_tag: Option<String>,
}

create_builder!(
    SequenceAttributeBuilder,
    seq_attributes,
    SequenceAttributes,
    locus_tag,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    SequenceFfn { value: String},
    SequenceFaa { value: String},
    CodonStart { value: u8 },
    Strand { value: i8 }
);

///product lines can contain difficult to parse punctuation such as biochemical symbols like unclosed single quotes, superscripts, single and double brackets etc.
///here we substitute these for an underscore
pub fn substitute_odd_punctuation(input: String) -> Result<String, anyhow::Error> {
    let re = Regex::new(r"[/?()',`]|[α-ωΑ-Ω]")?;

	// Strip either \r\n or \n more elegantly
	let cleaned = input.trim_end_matches(&['\r', '\n'][..]);

	Ok(re.replace_all(cleaned, "_").to_string())
}

///GFF3 field9 construct
#[derive(Debug)]
pub struct GFFInner {
    pub id: String,
    pub name: String,
    pub locus_tag: String,
    pub gene: String,
   // Inference: String,
   // Parent: String,
 //   db_xref: String,
    pub product: String,
   // is_circular: bool,
}

impl GFFInner {
    pub fn new(
      id: String,
      name: String,
      locus_tag: String,
      gene: String,
   //   Inference: String,
  //    Parent: String,
  //    db_xref: String,
      product: String,
    ) -> Self {
       GFFInner {
          id, name, locus_tag, gene, product,
          }
    }
}

///The main GFF3 construct
#[derive(Debug)]
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
       attributes: &'a GFFInner
    ) -> Self {
       GFFOuter {
          seqid, source, type_val, start, end, score, strand, phase, attributes,
          }
    }
    pub fn field9_attributes_build(&self) -> String {
       let mut full_field9 = Vec::new();
       if !self.attributes.id.is_empty() {
          full_field9.push(format!("id={}",self.attributes.id));
	  }
       if !self.attributes.name.is_empty() {
          full_field9.push(format!("name={}", self.attributes.name));
	  }
       if !self.attributes.gene.is_empty() {
          full_field9.push(format!("gene={}",self.attributes.gene));
	  }
   //    if !self.attributes.Inference.is_empty() {
   //       full_field9.push(format!("inference={}",self.attributes.Inference));
//	  }
       if !self.attributes.locus_tag.is_empty() {
          full_field9.push(format!("locus_tag={}",self.attributes.locus_tag));
	  }
       if !self.attributes.product.is_empty() {
          full_field9.push(format!("product={}",self.attributes.product));
	  }
   //    if !self.attributes.Parent.is_empty() {
   //       full_field9.push(format!("Parent={}",self.attributes.Parent));
//	  }
//       if !self.attributes.db_xref.is_empty() {
//          full_field9.push(format!("db_xref={}",self.attributes.db_xref));
//	  }
       full_field9.join(";")
       }
}

///formats the translation string which can be multiple lines, for gbk
pub fn format_translation(translation: &str) -> String {
	       //create method to add the protein sequence into the translation qualifier with correct line lengths
	       let mut formatted = String::new();
	       let cleaned_translation = translation.replace("\n", "");
	       formatted.push_str("                     /translation=\"");
	       let line_length: usize = 60;
	       let final_num = line_length - 15;
	       formatted.push_str(&format!("{}\n",&cleaned_translation[0..final_num]));
	       for i in (47..translation.len()).step_by(60) {
	             let end = i+60 -1;
		     let valid_end = if end >= translation.len() { &cleaned_translation.len() -1 } else { end };
                     formatted.push_str(&format!("                     {}",&cleaned_translation[i..valid_end]));
		     println!("cleaned translation leng is {:?}", &cleaned_translation[i..valid_end].len());
		     if *&cleaned_translation[i..valid_end].len() < 59 {
		        formatted.push('\"');
			}
		     else {
		        formatted.push('\n');
			}
		}
	       formatted
}

///writes the DNA sequence in gbk format with numbering
pub fn write_gbk_format_sequence(sequence: &str,file: &mut File) -> io::Result<()> {
       //function to write gbk format sequence
       writeln!(file, "ORIGIN")?;
       let mut formatted = String::new();
       let cleaned_input = sequence.replace("\n", "");
       let mut index = 1;
       for (_i, chunk) in cleaned_input.as_bytes().chunks(60).enumerate() {
           formatted.push_str(&format!("{:>5} ", index));
	   for (j, sub_chunk) in chunk.chunks(10).enumerate() {
	      if j > 0 {
	          formatted.push(' ');
		  }
	      formatted.push_str(&String::from_utf8_lossy(sub_chunk));
	    }
	    formatted.push('\n');
	    index+=60;
       }
       writeln!(file, "{:>6}", &formatted)?;
       writeln!(file, "//")?;
   Ok(())
}

///saves the parsed data in genbank format
//writes a genbank or multi-genbank file
pub fn gbk_write(seq_region: BTreeMap<String, (u32,u32)>, record_vec: Vec<Record>, filename: &str) -> io::Result<()> {
       let now = Local::now();
       let formatted_date = now.format("%d-%b-%Y").to_string().to_uppercase();
       let mut file = OpenOptions::new()
           .write(true)     // Allow writing to the file
           .append(true)    // Enable appending to the file
           .create(true)    // Create the file if it doesn't exist
           .open(filename)?;
       for (i, (key, _val)) in seq_region.iter().enumerate() {
	   let strain  = match &record_vec[i].source_map.get_strain(key) {
	          Some(value) => value.to_string(),
	          None => "Unknown".to_string(),
	          };
	   //write lines for the header
	   let organism  = match &record_vec[i].source_map.get_organism(key) {
	          Some(value) => value.to_string(),
	          None => "Unknown".to_string(),
	          };
           let mol_type  = match &record_vec[i].source_map.get_mol_type(key) {
	          Some(value) => value.to_string(),
	          None => "Unknown".to_string(),
	          };
	   let type_material  = match &record_vec[i].source_map.get_type_material(&key) {
	          Some(value) => value.to_string(),
	          None => "Unknown".to_string(),
                  };
	   let db_xref = match &record_vec[i].source_map.get_db_xref(key) {
	          Some(value) => value.to_string(),
	          None => "Unknown".to_string(),
	          };
	   let source_stop  = match &record_vec[i].source_map.get_stop(key) {
	          Some(value) => value.get_value(),
	          None => { println!("stop value not found");
	                   None }.expect("stop value not received")
	          };
           writeln!(file, "LOCUS       {}             {} bp    DNA     linear CON {}", &key,&record_vec[i].sequence.len(),&formatted_date)?;
	   writeln!(file, "DEFINITION  {} {}.", &organism, &strain)?;
	   writeln!(file, "ACCESSION   {}", &key)?;
	   writeln!(file, "KEYWORDS    .")?;
	   writeln!(file, "SOURCE      {} {}", &organism,&strain)?;
	   writeln!(file, "  ORGANISM  {} {}", &organism,&strain)?;
	   //write lines for the source
	   writeln!(file, "FEATURES             Location/Qualifiers")?;
	   writeln!(file, "     source          1..{}", &source_stop)?;
	   writeln!(file, "                     /organism=\"{}\"",&strain)?;
	   writeln!(file, "                     /mol_type=\"{}\"",&mol_type)?;
	   writeln!(file, "                     /strain=\"{}\"",&strain)?;
	   if type_material != *"Unknown".to_string() {
	       writeln!(file, "                     /type_material=\"{}\"",&type_material)?;
	       }
	   writeln!(file, "                     /db_xref=\"{}\"",&db_xref)?;
	   //write lines for each CDS
	   for (locus_tag, _value) in &record_vec[i].cds.attributes {
	      let start  = match &record_vec[i].cds.get_start(locus_tag) {
	          Some(value) => value.get_value(),
	          None => { println!("start value not found");
	                   None }.expect("start value not received")
	          };
	      let stop  = match &record_vec[i].cds.get_stop(locus_tag) {
	          Some(value) => value.get_value(),
	          None => { println!("stop value not found");
	                   None }.expect("stop value not received")
	          };
	      let product  = match &record_vec[i].cds.get_product(locus_tag) {
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
	       let gene  = match &record_vec[i].cds.get_gene(locus_tag) {
	                         Some(value) => value.to_string(),
	                         None => "unknown".to_string(),
	                         };
	       let translation = match &record_vec[i].seq_features.get_sequence_faa(locus_tag) {
	                         Some(value) => value.to_string(),
				 None => "unknown".to_string(),
				 };
	       if strand == 1 {
	           writeln!(file, "     gene            {}..{}",&start,&stop)?;
                   } else {
	           writeln!(file, "     gene            complement({}..{})",&start,&stop)?;
		   }
	       writeln!(file, "                     /locus_tag=\"{}\"",&locus_tag)?;
	       if strand == 1 {
	           writeln!(file, "     CDS             {}..{}",&start,&stop)?;
		   }
	       else {
	           writeln!(file, "     CDS             complement({}..{})",&start,&stop)?;
		   }
	       writeln!(file, "                     /locus_tag=\"{}\"",&locus_tag)?;
	       writeln!(file, "                     /codon_start=\"{}\"", &codon_start)?;
	       if gene != "unknown" {
	           writeln!(file, "                     /gene=\"{}\"", &gene)?;
		   }
	       if translation != "unknown" {
  	           let formatted_translation = format_translation(&translation);
		   writeln!(file, "{}", &formatted_translation)?;
		   }
	       writeln!(file, "                     /product=\"{}\"",&product)?;
	       }
	  write_gbk_format_sequence(&record_vec[i].sequence, &mut file)?;
	  }
	  Ok(())
}

///saves the parsed data in gff3 format
//writes a gff3 file from a genbank
#[allow(unused_assignments)]
#[allow(unused_variables)]
pub fn gff_write(seq_region: BTreeMap<String, (u32, u32)>, mut record_vec: Vec<Record>, filename: &str, dna: bool) -> io::Result<()> {
       let mut file = OpenOptions::new()
           //.write(true)     // Allow writing to the file
           .append(true)    // Enable appending to the file
           .create(true)    // Create the file if it doesn't exist
           .open(filename)?;
       if file.metadata()?.len() == 0 {
           writeln!(file, "##gff-version 3")?;
	   }
       let mut full_seq = String::new();
       let mut prev_end: u32 = 0;
       //println!("this is the full seq_region {:?}", &seq_region);
       for (k, v) in seq_region.iter() {
          writeln!(file, "##sequence-region\t{}\t{}\t{}", &k, v.0, v.1)?;
	  }
       for ((source_name, (seq_start, seq_end)), record) in seq_region.iter().zip(record_vec.drain(..)) {
	  if dna == true {
	     full_seq.push_str(&record.sequence);
             }
           for (locus_tag, _valu) in &record.cds.attributes {
               let start  = match record.cds.get_start(locus_tag) {
	          Some(value) => value.get_value(),
	          None => { println!("start value not found");
	                   None }.expect("start value not received")
	          };
	      let stop  = match record.cds.get_stop(locus_tag) {
	          Some(value) => value.get_value(),
	          None => { println!("stop value not found");
	                   None }.expect("stop value not received")
	          };
	      let gene  = match record.cds.get_gene(locus_tag) {
	          Some(value) => value.to_string(),
	          None => "unknown".to_string(),
	          };
	      let product  = match record.cds.get_product(locus_tag) {
	          Some(value) => value.to_string(),
	          None => "unknown product".to_string(),
	          };
	      let strand  = match record.cds.get_strand(locus_tag) {
	          Some(valu) => {
	             match valu {
		        1 => "+".to_string(),
		        -1 => "-".to_string(),
		        _ => { println!("unexpected strand value {} for locus_tag {}", valu, locus_tag);
		            "unknownstrand".to_string() }
		     }
	          },
	          None => "unknownvalue".to_string(),
	       };
	       let phase = match record.cds.get_codon_start(locus_tag) {
	          Some(valuer) => {
	             match valuer {
		        1 => 0,
		        2 => 1,
		        3 => 2,
		        _ => { println!("unexpected phase value {} in the bagging area for locus_tag {}", valuer, locus_tag);
		            1 }
		     }
	          },
	          None => 1,
	       };
              let gff_inner = GFFInner::new(
                 locus_tag.to_string(),
	         source_name.clone(),
	         locus_tag.to_string(),
	         gene,
	    //  &record.cds.get_Inference(&locus_tag),
	    //  &record.cds.get_Parent(&locus_tag),
	   //   db_xref,
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
	      //println!("{}\t{}\t{}\t{:?}\t{:?}\t{}\t{}\t{}\t{}", gff_outer.seqid, gff_outer.source, gff_outer.type_val, gff_outer.start, gff_outer.end, gff_outer.score, gff_outer.strand, gff_outer.phase, field9_attributes);
              writeln!(file, "{}\t{}\t{}\t{:?}\t{:?}\t{}\t{}\t{}\t{}", gff_outer.seqid, gff_outer.source, gff_outer.type_val, gff_outer.start, gff_outer.end, gff_outer.score, gff_outer.strand, gff_outer.phase, field9_attributes)?;
          
	  }
	  prev_end = *seq_end;
	  }
          if dna {
             writeln!(file, "##FASTA")?;
	     //writeln!(file, ">{}\n",&filename.to_string())?;
             writeln!(file, "{}", full_seq)?;
	     }
          Ok(())
}
	       	       
///saves the parsed data in gff3 format
//writes a gff3 file from a genbank
#[allow(unused_assignments)]
pub fn orig_gff_write(seq_region: BTreeMap<String, (u32, u32)>, record_vec: Vec<Record>, filename: &str, dna: bool) -> io::Result<()> {
       let mut file = OpenOptions::new()
           //.write(true)     // Allow writing to the file
           .append(true)    // Enable appending to the file
           .create(true)    // Create the file if it doesn't exist
           .open(filename)?;
       if file.metadata()?.len() == 0 {
           writeln!(file, "##gff-version 3")?;
	   }
       let mut source_name = String::new();
       let mut full_seq = String::new();
       let mut prev_end: u32 = 0;
       //println!("this is the full seq_region {:?}", &seq_region);
       for (k, v) in seq_region.iter() {
          writeln!(file, "##sequence-region\t{}\t{}\t{}", &k, v.0, v.1)?;
	  }
       for (i, (key, val)) in seq_region.iter().enumerate() {
	  source_name = key.to_string();
	  if dna == true {
	     full_seq.push_str(&record_vec[i].sequence);
             }
           for (locus_tag, _valu) in &record_vec[i].cds.attributes {
               let start  = match record_vec[i].cds.get_start(locus_tag) {
	          Some(value) => value.get_value(),
	          None => { println!("start value not found");
	                   None }.expect("start value not received")
	          };
	      let stop  = match record_vec[i].cds.get_stop(locus_tag) {
	          Some(value) => value.get_value(),
	          None => { println!("stop value not found");
	                   None }.expect("stop value not received")
	          };
	      let gene  = match record_vec[i].cds.get_gene(locus_tag) {
	          Some(value) => value.to_string(),
	          None => "unknown".to_string(),
	          };
	      let product  = match record_vec[i].cds.get_product(locus_tag) {
	          Some(value) => value.to_string(),
	          None => "unknown product".to_string(),
	          };
	      let strand  = match record_vec[i].cds.get_strand(locus_tag) {
	          Some(valu) => {
	             match valu {
		        1 => "+".to_string(),
		        -1 => "-".to_string(),
		        _ => { println!("unexpected strand value {} for locus_tag {}", valu, locus_tag);
		            "unknownstrand".to_string() }
		     }
	          },
	          None => "unknownvalue".to_string(),
	       };
	       let phase = match record_vec[i].cds.get_codon_start(locus_tag) {
	          Some(valuer) => {
	             match valuer {
		        1 => 0,
		        2 => 1,
		        3 => 2,
		        _ => { println!("unexpected phase value {} in the bagging area for locus_tag {}", valuer, locus_tag);
		            1 }
		     }
	          },
	          None => 1,
	       };
              let gff_inner = GFFInner::new(
                 locus_tag.to_string(),
	         source_name.clone(),
	         locus_tag.to_string(),
	         gene,
	    //  &record.cds.get_Inference(&locus_tag),
	    //  &record.cds.get_Parent(&locus_tag),
	   //   db_xref,
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
	      //println!("{}\t{}\t{}\t{:?}\t{:?}\t{}\t{}\t{}\t{}", gff_outer.seqid, gff_outer.source, gff_outer.type_val, gff_outer.start, gff_outer.end, gff_outer.score, gff_outer.strand, gff_outer.phase, field9_attributes);
              writeln!(file, "{}\t{}\t{}\t{:?}\t{:?}\t{}\t{}\t{}\t{}", gff_outer.seqid, gff_outer.source, gff_outer.type_val, gff_outer.start, gff_outer.end, gff_outer.score, gff_outer.strand, gff_outer.phase, field9_attributes)?;
          
	  }
	  prev_end = val.1;
	  }
          if dna {
             writeln!(file, "##FASTA")?;
	     //writeln!(file, ">{}\n",&filename.to_string())?;
             writeln!(file, "{}", full_seq)?;
	     }
          Ok(())
}

///internal record containing data from a single source or contig.  Has multiple features.
//sets up a record
#[derive(Debug, Clone)]
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

// Provide a type alias and conversion to a generic record to aid interoperability
pub type GenericRecordGbk = crate::record::GenericRecord<SourceAttributeBuilder, FeatureAttributeBuilder, SequenceAttributeBuilder>;

impl From<&Record> for GenericRecordGbk {
    fn from(r: &Record) -> Self {
        Self {
            id: r.id.clone(),
            seq: r.sequence.clone(),
            seqid: r.id.clone(),
            start: r.start as u32,
            end: r.end as u32,
            strand: r.strand,
            source: r.source_map.clone(),
            cds: r.cds.clone(),
            seq_features: r.seq_features.clone(),
        }
    }
}

#[allow(dead_code)]
pub struct Config {
    filename: String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &str> {
    if args.len() < 2 {
        panic!("not enough arguments, please provide filename");
    }
    let filename = args[1].clone();

    Ok(Config { filename })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    #[allow(unused_mut)]
    #[allow(unused_variables)]
    #[allow(dead_code)]
    #[allow(unused_assignments)]
    #[allow(unused_imports)]
    fn test_read_file() {
       let content = std::fs::read_to_string("K12_ribo.gbk").expect("error reading file");
       assert!(content.contains("LOCUS"));
       assert!(content.len() > 0);
       }
    #[test]
    #[allow(unused_mut)]
    #[allow(unused_variables)]
    #[allow(dead_code)]
    #[allow(unused_assignments)]
    #[allow(unused_imports)]
    fn test_parse_gbk() {
       let file_gbk = "K12_ribo.gbk";
       let records = genbank!(&file_gbk);
       assert!(records.len() > 0);
       }
    #[test]
    #[allow(unused_mut)]
    #[allow(unused_variables)]
    #[allow(dead_code)]
    #[allow(unused_assignments)]
    #[allow(unused_imports)]
    fn test_parse_source_attributes() {
       let file_gbk = "K12_ribo.gbk";
       let records = genbank!(&file_gbk);
       if let Some(record) = records.first() {
          if let Some((key, val)) = record.source_map.source_attributes.first_key_value() {
	      assert_eq!(key, &"source_NC_000913_1".to_string());
	      }
          }
       }
    #[test]
    #[allow(unused_mut)]
    #[allow(unused_variables)]
    #[allow(dead_code)]
    #[allow(unused_assignments)]
    #[allow(unused_imports)]
    fn test_parse_cds_attributes() {
       let file_gbk = "K12_ribo.gbk";
       let records = genbank!(&file_gbk);
       if let Some(record) = records.first() {
          if let Some((locus_tag, vals)) = record.cds.attributes.first_key_value() {
                assert_eq!(locus_tag, &"b3304".to_string());
                assert_eq!(record.cds.get_gene(&locus_tag).as_deref(), Some(&"rplR".to_string()));
		}
          }
       }
    #[test]
    #[allow(unused_mut)]
    #[allow(unused_variables)]
    #[allow(dead_code)]
    #[allow(unused_assignments)]
    #[allow(unused_imports)]
    fn test_parse_sequence_attributes() {
        let file_gbk = "K12_ribo.gbk";
	let records = genbank!(&file_gbk);
	if let Some(record) = records.first() {
	   if let Some((key, vals)) = record.cds.attributes.first_key_value() {
	       assert_eq!(key, &"b3304".to_string());
	       assert_eq!(record.seq_features.get_sequence_faa(&key), Some(&"MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAEQLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF".to_string()));
	       }
	   }
	}
}
    
