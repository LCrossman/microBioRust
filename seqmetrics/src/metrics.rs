//!  The purpose of seqmetrics is to allow calculations of protein sequence parameters
//!  We define a set of amino acid getters to allow getting of the protein sequence
//!  as either the three letter code such as Trp, Phe, the full name such as alanine, glycine
//!  or the single letter code such as Y, A
//!
//!   Example function to calculate transmembrane statistics
//!
//! ```rust
//!
//! use clap::Parser;
//! use std::fs::File;
//! use microBioRust::gbk::Reader;
//! use std::io;
//! use std::collections::HashMap;
//! use microBioRust_seqmetrics::metrics::hydrophobicity;
//!
//! pub fn suggest_transmembrane_domains() -> Result<(), anyhow::Error> {
//!            let file_gbk = File::open("test_output.gbk")?;
//!            let mut reader = Reader::new(file_gbk);
//!            let mut records = reader.records();
//!            loop {  
//!                match records.next() {  
//!                    Some(Ok(mut record)) => {
//!                       //println!("next record");
//!                       //println!("Record id: {:?}", record.id);
//!                       for (k, v) in &record.cds.attributes {
//!                           match record.seq_features.get_sequence_faa(&k) {
//!                                     Some(value) => { let seq_faa = value.to_string();
//!				                      println!("{:?}", &seq_faa);
//!				                      let hydro_values = hydrophobicity(&seq_faa, 21);
//!						      let mut result = String::new();
//!						      for hydro in hydro_values {
//!						           if hydro > 1.6 {
//!						               println!("possible transmembrane region - score {}",&hydro);  
//!							       }
//!						           else {
//!						               ()
//!							   }
//!						      }
//!                                                 },
//!                                     _ => (),
//!                                     };
//!                       
//!                           }
//!                    },
//!                    Some(Err(e)) => { println!("theres an err {:?}", e); },
//!                    None => {
//!                             println!("finished iteration");
//!                             break;
//!                             },
//!                    };
//!               }
//!            return Ok(());
//!   }
//!```   
//!
//!   Example function to calculate the molecular weight of a protein sequence
//!
//!```rust
//!
//! use clap::Parser;
//! use std::fs::File;
//! use microBioRust::gbk::Reader;
//! use std::io;
//! use std::collections::HashMap;
//! use microBioRust_seqmetrics::metrics::molecular_weight;
//!
//! pub fn collect_molecular_weight() -> Result<(), anyhow::Error> {
//!            let file_gbk = File::open("test_output.gbk")?;
//!            let mut reader = Reader::new(file_gbk);
//!            let mut records = reader.records();
//!	    let mut molecular_weight_total: f64 = 0.0;
//!            loop {  
//!                match records.next() {  
//!                    Some(Ok(mut record)) => {
//!                       //println!("next record");
//!                       //println!("Record id: {:?}", record.id);
//!                      for (k, v) in &record.cds.attributes {
//!                           match record.seq_features.get_sequence_faa(&k) {
//!                                     Some(value) => {
//!                                                     let seq_faa = value.to_string();
//!				                        println!("id: {:?}", &k);
//!				                        molecular_weight_total = molecular_weight(&seq_faa);
//!                                                     println!(">{}|{}\n{}", &record.id, &k, molecular_weight_total);
//!                                                    },
//!                                     _ => (),
//!                                     };
//!                       
//!                           }
//!                    },
//!                    Some(Err(e)) => { println!("theres an err {:?}", e); },
//!                    None => {
//!                              println!("finished iteration");
//!                              break;
//!                            },
//!                    }
//!               }
//!            return Ok(());
//!   }
//!```
//!
//!   Example function to count amino acids of each protein as raw counts, see below to generate percentages per protein
//!
//!```rust
//!
//! use clap::Parser;
//! use std::fs::File;
//! use microBioRust::gbk::Reader;
//! use std::io;
//! use std::collections::HashMap;
//! use microBioRust_seqmetrics::metrics::amino_counts;
//!
//! pub fn count_aminos() -> Result<(), anyhow::Error> {
//!            let file_gbk = File::open("test_output.gbk")?;
//!            let mut reader = Reader::new(file_gbk);
//!            let mut records = reader.records();
//!	    let mut results: HashMap<char, u64> = HashMap::new();
//!            loop {  
//!                match records.next() {  
//!                   Some(Ok(mut record)) => {
//!                       //println!("next record");
//!                       //println!("Record id: {:?}", record.id);
//!                       for (k, v) in &record.cds.attributes {
//!                           match record.seq_features.get_sequence_faa(&k) {
//!                                     Some(value) => { let seq_faa = value.to_string();
//!				                      println!("id: {:?}", &k);
//!				                      results = amino_counts(&seq_faa);
//!                                                      println!(">{}|{}\n{:?}", &record.id, &k, results);
//!                                                      },
//!                                     _ => (),
//!                                     };
//!                       
//!                           }
//!                    },
//!                    Some(Err(e)) => { println!("theres an err {:?}", e); },
//!                    None => {
//!                             println!("finished iteration");
//!                             break;
//!                             },
//!                    }
//!               }
//!            return Ok(());
//!   }
//!```
//!  Example function to calculate and print out the aromaticity of each protein
//!
//!```rust
//! use clap::Parser;
//! use std::fs::File;
//! use microBioRust::gbk::Reader;
//! use std::io;
//! use std::collections::HashMap;
//! use microBioRust_seqmetrics::metrics::amino_percentage;
//!
//! pub fn aromaticity() -> Result<(), anyhow::Error> {
//!        // calculated as in biopython with aromaticity according to Lobry, 1994 as the relative freq of Phe+Trp+Tyr
//!        let file_gbk = File::open("test_output.gbk")?;
//!	let mut reader = Reader::new(file_gbk);
//!	let mut records = reader.records();
//!	let mut results: HashMap<char, f64> = HashMap::new();
//!	loop {
//!	   match records.next() {
//!	      Some(Ok(mut record)) => {
//!	          for (k, v) in &record.cds.attributes {
//!		     match record.seq_features.get_sequence_faa(&k) {
//!		         Some(value) => {  let seq_faa = value.to_string();
//!			                   results = amino_percentage(&seq_faa);
//!					   let aromatic_aas = vec!['Y','W','F'];
//!					   let aromaticity: f64 = aromatic_aas.iter()
//!					       .filter_map(|&amino| results.get(&amino))
//!					       .map(|&perc| perc / 100.0)
//!					       .sum();
//!					   println!("aromaticity for {} {} is {}",&record.id, &k, &aromaticity);
//!					  },
//!			_ => (),
//!			};
//!		   }
//!	         },
//!	    Some(Err(e)) => { println!("theres an error {:?}", e); },
//!	    None => {
//!                   println!("finished iteration");
//!	              break;
//!		    },
//!	    }
//!       }
//!      return Ok(());
//!   }
//!```
//!  The purpose of hamming.rs is to allow calculations of the hamming distances between sequences
//!  The Hamming distance is the minimum number of substitutions required to change one string to another
//!  It is one of several string metrics for measuring the edit distance between two sequences.
//!  It does not encompass or take into account any biology
//!  It is named after the American Mathematician Richard Hamming (wikipedia)
//!  This is aimed essentially at protein fasta sequences
//!  
//!
//!  ```
//!  use microBioRust_seqmetrics::hamming::hamming_matrix;
//!  use microBioRust_seqmetrics::write_dst_csv::write_distances_csv;
//!  use tokio::fs::File;
//!  use std::collections::HashMap;
//!  use bio::io::fasta;
//!  use tokio::io;
//!  use tokio::io::{AsyncWriteExt, BufWriter};
//!
//!
//!  #[tokio::main]
//!  async fn main() -> Result<(), anyhow::Error> {
//!            let reader = fasta::Reader::new(std::io::stdin());
//!            let records: Vec<_> = reader.records().collect::<Result<_, _>>()?;
//!	    println!("gathering records");
//!            let sequences: Vec<String> = records
//!	                          .iter()
//!				  .map(|rec| String::from_utf8_lossy(rec.seq()).to_string())
//!				  .collect();
//!            let ids: Vec<String> = records
//!	                          .iter()
//!				  .map(|rec| rec.id().to_string())
//!				  .collect();
//!	    println!("gathered ids");
//!	    let distances = hamming_matrix(&sequences).await?;
//!	    write_distances_csv(ids, distances, "hamming_dists.csv").await?;
//!
//!         Ok(())
//!  }
//!  ```

#![allow(unused_imports)]
use crate::hamming::hamming_matrix;
use crate::write_dst_csv::write_distances_csv;
use bio::io::fasta;
use microBioRust::gbk::Reader;
use std::collections::HashMap;
use std::fs::File;

// Define a macro to generate the getters for each amino acid
#[macro_export]
macro_rules! amino_acid_getters {
    ($struct_name:ident, $( ($field:ident, $full_name:ident, $three_letter:ident, $single_letter:ident) ),* ) => {
        #[allow(non_snake_case)]
	#[allow(dead_code)]
        impl $struct_name {
            $(
	        // Capital full name getter
		fn $field(&self) -> f64 {
		   self.$field
		}
                // Full name getter
                fn $full_name(&self) -> f64 {
                    self.$field
                }
                // Three-letter code getter
                fn $three_letter(&self) -> f64 {
                    self.$field
                }
                // Single-letter code getter
                fn $single_letter(&self) -> f64 {
                    self.$field
                }
            )*
        }
    };
}

#[allow(non_snake_case)]
#[allow(dead_code)]
pub struct MolWeights {
    Alanine: f64,
    Arginine: f64,
    Asparagine: f64,
    Aspartate: f64,
    Cysteine: f64,
    Glutamate: f64,
    Glutamine: f64,
    Glycine: f64,
    Histidine: f64,
    Isoleucine: f64,
    Leucine: f64,
    Lysine: f64,
    Methionine: f64,
    Phenylalanine: f64,
    Proline: f64,
    Serine: f64,
    Threonine: f64,
    Tryptophan: f64,
    Tyrosine: f64,
    Valine: f64,
}

#[allow(non_snake_case)]
#[allow(dead_code)]
impl MolWeights {
    fn new() -> Self {
        Self {
            //masses from NIST chemistry webbook US Dept of commerce
            Alanine: 89.0932,        //C3H7NO2
            Arginine: 174.2010,      //C6H14N4O2
            Asparagine: 132.1179,    //C4H8N2O3
            Aspartate: 133.1027,     //C4H7NO4
            Cysteine: 121.158,       //C3H7NO2S
            Glutamate: 147.1293,     //C5H9NO4
            Glutamine: 146.1445,     //C5H10N2O3
            Glycine: 75.0666,        //C2H5NO2
            Histidine: 155.1546,     //C6H9N3O2
            Isoleucine: 131.1729,    //C6H13NO2
            Leucine: 131.1729,       //C6H13NO2
            Lysine: 146.1876,        //C6H14N2O2
            Methionine: 149.211,     //C5H11NO2S
            Phenylalanine: 165.1891, //C9H11NO2
            Proline: 115.1305,       //C5H9NO2
            Serine: 105.0926,        //C3H7NO2
            Threonine: 119.1192,     //C4H9NO3
            Tryptophan: 204.2252,    //C11H12N2O2
            Tyrosine: 181.1885,      //C9H11NO3
            Valine: 117.1463,        //C5H11NO2
        }
    }
}

amino_acid_getters!(
    MolWeights,
    (Alanine, alanine, Ala, A),
    (Arginine, arginine, Arg, R),
    (Asparagine, asparagine, Asn, N),
    (Aspartate, aspartate, Asp, D),
    (Cysteine, cysteine, Cys, C),
    (Glutamine, glutamine, Gln, Q),
    (Glutamate, glutamate, Glu, E),
    (Glycine, glycine, Gly, G),
    (Histidine, histidine, His, H),
    (Isoleucine, isoleucine, Ile, I),
    (Leucine, leucine, Leu, L),
    (Lysine, lysine, Lys, K),
    (Methionine, methionine, Met, M),
    (Phenylalanine, phenylalanine, Phe, F),
    (Proline, proline, Pro, P),
    (Serine, serine, Ser, S),
    (Threonine, threonine, Thr, T),
    (Tryptophan, tryptophan, Trp, W),
    (Tyrosine, tyrosine, Tyr, Y),
    (Valine, valine, Val, V)
);

#[allow(non_snake_case)]
#[allow(dead_code)]
#[allow(unused_variables)]
pub fn molecular_weight(protein_seq: &str) -> f64 {
    let amino_weights: MolWeights = MolWeights::new();
    let mut total_weight = 0.0;
    for ch in protein_seq.chars() {
        match ch {
            'A' => total_weight += amino_weights.A(),
            'R' => total_weight += amino_weights.R(),
            'N' => total_weight += amino_weights.N(),
            'D' => total_weight += amino_weights.D(),
            'C' => total_weight += amino_weights.C(),
            'Q' => total_weight += amino_weights.Q(),
            'E' => total_weight += amino_weights.E(),
            'G' => total_weight += amino_weights.G(),
            'H' => total_weight += amino_weights.H(),
            'I' => total_weight += amino_weights.I(),
            'L' => total_weight += amino_weights.L(),
            'K' => total_weight += amino_weights.K(),
            'M' => total_weight += amino_weights.M(),
            'F' => total_weight += amino_weights.F(),
            'P' => total_weight += amino_weights.P(),
            'S' => total_weight += amino_weights.S(),
            'T' => total_weight += amino_weights.T(),
            'W' => total_weight += amino_weights.W(),
            'Y' => total_weight += amino_weights.Y(),
            'V' => total_weight += amino_weights.V(),
            _ => continue,
        }
    }
    let result_weight = total_weight - ((protein_seq.len() - 1) as f64 * 18.02);
    result_weight
}

#[allow(non_snake_case)]
#[allow(dead_code)]
pub struct Hydrophobicity {
    Alanine: f64,
    Arginine: f64,
    Asparagine: f64,
    Aspartate: f64,
    Cysteine: f64,
    Glutamate: f64,
    Glutamine: f64,
    Glycine: f64,
    Histidine: f64,
    Isoleucine: f64,
    Leucine: f64,
    Lysine: f64,
    Methionine: f64,
    Phenylalanine: f64,
    Proline: f64,
    Serine: f64,
    Threonine: f64,
    Tryptophan: f64,
    Tyrosine: f64,
    Valine: f64,
}

impl Hydrophobicity {
    #[allow(non_snake_case)]
    fn new_KD() -> Self {
        Self {
            //Kyte-Doolittle values from the Qiagen resources website
            Alanine: 1.80,
            Arginine: -4.50,
            Asparagine: -3.50,
            Aspartate: -3.50,
            Cysteine: 2.50,
            Glutamate: -3.50,
            Glutamine: -3.50,
            Glycine: -0.40,
            Histidine: -3.20,
            Isoleucine: 4.50,
            Leucine: 3.80,
            Lysine: -3.90,
            Methionine: 1.90,
            Phenylalanine: 2.80,
            Proline: -1.60,
            Serine: -0.80,
            Threonine: -0.70,
            Tryptophan: -0.90,
            Tyrosine: -1.30,
            Valine: 4.20,
        }
    }
}

amino_acid_getters!(
    Hydrophobicity,
    (Alanine, alanine, Ala, A),
    (Arginine, arginine, Arg, R),
    (Asparagine, asparagine, Asn, N),
    (Aspartate, aspartate, Asp, D),
    (Cysteine, cysteine, Cys, C),
    (Glutamine, glutamine, Gln, Q),
    (Glutamate, glutamate, Glu, E),
    (Glycine, glycine, Gly, G),
    (Histidine, histidine, His, H),
    (Isoleucine, isoleucine, Ile, I),
    (Leucine, leucine, Leu, L),
    (Lysine, lysine, Lys, K),
    (Methionine, methionine, Met, M),
    (Phenylalanine, phenylalanine, Phe, F),
    (Proline, proline, Pro, P),
    (Serine, serine, Ser, S),
    (Threonine, threonine, Thr, T),
    (Tryptophan, trytophan, Trp, W),
    (Tyrosine, tyrosine, Tyr, Y),
    (Valine, valine, Val, V)
);

#[allow(non_snake_case)]
#[allow(dead_code)]
#[allow(unused_mut)]
#[allow(unused_variables)]
pub fn hydrophobicity(protein_seq: &str, window_size: usize) -> Vec<f64> {
    let mut hydrophobicity: Hydrophobicity = Hydrophobicity::new_KD();
    let mut total_hydrophobicity: Vec<f64> = Vec::new();
    let mut window_values: f64 = 0.0;
    let mut windows: Vec<String> = protein_seq
        .chars()
        .collect::<Vec<_>>()
        .windows(window_size)
        .map(|window| window.iter().collect())
        .collect();
    for (index, window) in windows.iter().enumerate() {
        for ch in window.chars() {
            match ch {
                'A' => window_values += hydrophobicity.A(),
                'R' => window_values += hydrophobicity.R(),
                'N' => window_values += hydrophobicity.N(),
                'D' => window_values += hydrophobicity.D(),
                'C' => window_values += hydrophobicity.C(),
                'Q' => window_values += hydrophobicity.Q(),
                'E' => window_values += hydrophobicity.E(),
                'G' => window_values += hydrophobicity.G(),
                'H' => window_values += hydrophobicity.H(),
                'I' => window_values += hydrophobicity.I(),
                'L' => window_values += hydrophobicity.L(),
                'K' => window_values += hydrophobicity.K(),
                'M' => window_values += hydrophobicity.M(),
                'F' => window_values += hydrophobicity.F(),
                'P' => window_values += hydrophobicity.P(),
                'S' => window_values += hydrophobicity.S(),
                'T' => window_values += hydrophobicity.T(),
                'W' => window_values += hydrophobicity.W(),
                'Y' => window_values += hydrophobicity.Y(),
                'V' => window_values += hydrophobicity.V(),
                _ => continue,
            }
        }
        total_hydrophobicity.push(window_values);
    }
    total_hydrophobicity
}

pub fn amino_counts(protein_seq: &str) -> HashMap<char, u64> {
    let mut counts: HashMap<char, u64> = HashMap::new();
    for c in protein_seq.chars() {
        *counts.entry(c).or_insert(0) += 1;
    }
    counts
}

#[allow(unused_mut)]
#[allow(unused_variables)]
#[allow(unused_assignments)]
#[allow(dead_code)]
pub fn amino_percentage(protein_seq: &str) -> HashMap<char, f64> {
    let mut percentages: HashMap<char, f64> = HashMap::new();
    let counts = amino_counts(protein_seq);
    let seq_len: f64 = (protein_seq.len() as f64) as f64;
    percentages = counts
        .iter()
        .map(|(k, &value)| {
            let percentage = (value as f64 / seq_len) * 100.0;
            (k.clone(), percentage)
        })
        .collect();
    percentages
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(unused_mut)]
    #[allow(dead_code)]
    #[allow(unused_variables)]
    pub fn suggest_transmembrane_domains() -> Result<(), anyhow::Error> {
        let file_gbk = File::open("test_output.gbk")?;
        let mut reader = Reader::new(file_gbk);
        let mut records = reader.records();
        loop {
            match records.next() {
                Some(Ok(mut record)) => {
                    //println!("next record");
                    //println!("Record id: {:?}", record.id);
                    for (k, v) in &record.cds.attributes {
                        match record.seq_features.get_sequence_faa(&k) {
                            Some(value) => {
                                let seq_faa = value.to_string();
                                println!("{:?}", &seq_faa);
                                let hydro_values = hydrophobicity(&seq_faa, 21);
                                let mut result = String::new();
                                for hydro in hydro_values {
                                    if hydro > 1.6 {
                                        println!(
                                            "possible transmembrane region - score {}",
                                            &hydro
                                        );
                                    } else {
                                        ()
                                    }
                                }
                            }
                            _ => (),
                        };
                    }
                }
                Some(Err(e)) => {
                    println!("theres an err {:?}", e);
                }
                None => {
                    println!("finished iteration");
                    break;
                }
            }
        }
        return Ok(());
    }

    #[test]
    #[allow(unused_mut)]
    #[allow(dead_code)]
    #[allow(unused_variables)]
    #[allow(unused_assignments)]
    pub fn collect_molecular_weight() -> Result<(), anyhow::Error> {
        let file_gbk = File::open("test_output.gbk")?;
        let mut reader = Reader::new(file_gbk);
        let mut records = reader.records();
        let mut molecular_weight_total: f64 = 0.0;
        loop {
            match records.next() {
                Some(Ok(mut record)) => {
                    //println!("next record");
                    //println!("Record id: {:?}", record.id);
                    for (k, v) in &record.cds.attributes {
                        match record.seq_features.get_sequence_faa(&k) {
                            Some(value) => {
                                let seq_faa = value.to_string();
                                println!("id: {:?}", &k);
                                molecular_weight_total = molecular_weight(&seq_faa);
                                println!(">{}|{}\n{}", &record.id, &k, molecular_weight_total);
                            }
                            _ => (),
                        };
                    }
                }
                Some(Err(e)) => {
                    println!("theres an err {:?}", e);
                }
                None => {
                    println!("finished iteration");
                    break;
                }
            }
        }
        return Ok(());
    }

    #[test]
    #[allow(unused_mut)]
    #[allow(unused_variables)]
    #[allow(unused_assignments)]
    pub fn count_aminos() -> Result<(), anyhow::Error> {
        let file_gbk = File::open("test_output.gbk")?;
        let mut reader = Reader::new(file_gbk);
        let mut records = reader.records();
        let mut results: HashMap<char, u64> = HashMap::new();
        loop {
            match records.next() {
                Some(Ok(record)) => {
                    for (k, _v) in &record.cds.attributes {
                        match record.seq_features.get_sequence_faa(&k) {
                            Some(value) => {
                                let seq_faa = value.to_string();
                                println!("id: {:?}", &k);
                                results = amino_counts(&seq_faa);
                                println!(">{}|{}\n{:?}", &record.id, &k, results);
                            }
                            _ => (),
                        };
                    }
                }
                Some(Err(e)) => {
                    println!("theres an err {:?}", e);
                }
                None => {
                    println!("finished iteration");
                    break;
                }
            }
        }
        return Ok(());
    }

    #[test]
    #[allow(dead_code)]
    #[allow(unused_mut)]
    #[allow(unused_variables)]
    #[allow(unused_assignments)]
    pub fn aromaticity() -> Result<(), anyhow::Error> {
        // calculated as in biopython with aromaticity according to Lobry, 1994 as the relative freq of Phe+Trp+Tyr
        let file_gbk = File::open("test_output.gbk")?;
        let reader = Reader::new(file_gbk);
        let mut records = reader.records();
        let mut results: HashMap<char, f64> = HashMap::new();
        loop {
            match records.next() {
                Some(Ok(record)) => {
                    for (k, _v) in &record.cds.attributes {
                        match record.seq_features.get_sequence_faa(&k) {
                            Some(value) => {
                                let seq_faa = value.to_string();
                                results = amino_percentage(&seq_faa);
                                let aromatic_aas = vec!['Y', 'W', 'F'];
                                let aromaticity: f64 = aromatic_aas
                                    .iter()
                                    .filter_map(|&amino| results.get(&amino))
                                    .map(|&perc| perc / 100.0)
                                    .sum();
                                println!(
                                    "aromaticity for {} {} is {}",
                                    &record.id, &k, &aromaticity
                                );
                            }
                            _ => (),
                        };
                    }
                }
                Some(Err(e)) => {
                    println!("theres an error {:?}", e);
                }
                None => {
                    println!("finished iteration");
                    break;
                }
            }
        }
        return Ok(());
    }
    #[tokio::test]
    pub async fn main() -> Result<(), anyhow::Error> {
        let reader = fasta::Reader::new(File::open("test_hamming.aln")?);
        let records: Vec<_> = reader.records().collect::<Result<_, _>>()?;
        let sequences: Vec<String> = records
            .iter()
            .map(|rec| String::from_utf8_lossy(rec.seq()).to_string())
            .collect();
        let ids: Vec<String> = records.iter().map(|rec| rec.id().to_string()).collect();
        let distances = hamming_matrix(&sequences).await?;
        let _ = write_distances_csv(ids, distances, "hamming_dists.csv");
        Ok(())
    }
}
