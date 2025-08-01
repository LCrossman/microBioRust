//!  The purpose of microbiorust-py is to allow access to microBioRust
//!  from Python - so Python InterOperability via pyo3
//!  This is a collection of pyfunctions to build a PyModule
//!  From Rust you can check the pyfunctions are registered to the PyModule using:
//!  
//!  cargo test
//!  
//!  To use the PyModule in Python you will need to run
//!  
//!  maturin develop
//!  
//!  Once developed, the PyModule can be loaded into Python and used:
//!  
//!  import microbiorust
//!  result = gbk_to_faa("test_input.gbk")
//!  for r in result:
//!      print(r)
//!  gbk_to_gff("test_input.gbk")
//!  
//!  Other pyfunctions that can be run include gbk_to_faa, embl_to_faa, gbk_to_gff and embl_to_gff


use pyo3::{
   prelude::*,
   types::PyModule,
};
use microBioRust::genbank;
use std::{
   collections::BTreeMap,
   io::{self, Write},
   fs::OpenOptions,
};
use microBioRust::gbk::{Record, Reader, RangeValue, gff_write};
use microBioRust::embl;
use microBioRust::embl::gff_write as embl_gff_write;

#[pyfunction]
fn gbk_to_faa(filename: &str) -> PyResult<Vec<String>> {
    let records = genbank!(&filename);
    let mut result = Vec::new();
    for record in records {
        for (k, _v) in &record.cds.attributes {
            if let Some(seq) = record.seq_features.get_sequence_faa(k) {
                result.push(format!(">{}|{}\n{}", &record.id, &k, seq));
            }
        }
    }
    Ok(result)
}

#[pyfunction]
fn embl_to_faa(filename: &str) -> PyResult<Vec<String>> {
    let records = genbank!(&filename);
    let mut result = Vec::new();
    for record in records {
        for (k, _v) in &record.cds.attributes {
            if let Some(seq) = record.seq_features.get_sequence_faa(k) {
                result.push(format!(">{}|{}\n{}", &record.id, &k, seq));
            }
        }
    }
    Ok(result)
}

#[allow(unused_variables)]
#[pyfunction]
fn gbk_to_gff(filename: &str, dna: bool) -> PyResult<()> {
    let records = genbank!(&filename);
    let prev_start: u32 = 0;
    let mut prev_end: u32 = 0;
    let mut seq_region: BTreeMap<String, (u32, u32)> = BTreeMap::new();
    let mut record_vec = Vec::new();
    let mut read_counter = 0;
    for record in records {
        if let Some(ref source) = record.source_map.source_name {
            let beginning = record
                .source_map
                .get_start(&source)
                .map_or(0, |v| v.get_value());
            let ending = record
                .source_map
                .get_stop(&source)
                .map_or(0, |v| v.get_value());
            if (ending + prev_end) < (beginning + prev_end) {
                println!(
                    "debug: end value is smaller than the start value at {:?}",
                    beginning
                );
            }
            seq_region.insert(
                source.to_string(),
                (beginning + prev_end, ending + prev_end),
            );
            record_vec.push(record);
            read_counter += 1;
            prev_end += ending; // this is to create the joined record if there are multiple
        } else {
            println!("missing record source name, skipping");
        }
    }
    let output_file = format!("{}.gff", &filename);
    if std::path::Path::new(&output_file).exists() {
        println!("deleting existing file {:?}", &output_file);
        std::fs::remove_file(&output_file).expect("Issue deleting output filename");
    }
    let _ = gff_write(seq_region.clone(), record_vec, &output_file, dna);
    println!("total records processed: {}", read_counter);
    return Ok(());
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[pyfunction]
fn embl_to_gff(filename: &str, dna: bool) -> PyResult<()> {
    let records = embl!(&filename);
    let prev_start: u32 = 0;
    let mut prev_end: u32 = 0;
    let mut seq_region: BTreeMap<String, (u32, u32)> = BTreeMap::new();
    let mut record_vec = Vec::new();
    let mut read_counter = 0;
    for record in records {
        if let Some(ref source) = record.source_map.source_name {
            let beginning = record
                .source_map
                .get_start(&source)
                .map_or(0, |v| v.get_value());
            let ending = record
                .source_map
                .get_stop(&source)
                .map_or(0, |v| v.get_value());
            if (ending + prev_end) < (beginning + prev_end) {
                println!(
                    "debug: end value is smaller than the start value at {:?}",
                    beginning
                );
            }
            seq_region.insert(
                source.to_string(),
                (beginning + prev_end, ending + prev_end),
            );
            record_vec.push(record);
            read_counter += 1;
            prev_end += ending; // this is to create the joined record if there are multiple
        } else {
            println!("missing record source name, skipping");
        }
    }
    let output_file = format!("{}.gff", &filename);
    if std::path::Path::new(&output_file).exists() {
        println!("deleting existing file {:?}", &output_file);
        std::fs::remove_file(&output_file).expect("Issue deleting output filename");
    }
    let _ = embl_gff_write(seq_region.clone(), record_vec, &output_file, dna);
    println!("total records processed: {}", read_counter);
    return Ok(());
}

#[pymodule]
fn microbiorust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gbk_to_faa, m)?)?;
    m.add_function(wrap_pyfunction!(embl_to_faa, m)?)?;
    m.add_function(wrap_pyfunction!(gbk_to_gff, m)?)?;
    m.add_function(wrap_pyfunction!(embl_to_gff, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::microbiorust;
    use pyo3::prelude::*;
    use pyo3::types::PyModule;

    #[test]
    fn test_functions_are_registered() {
        Python::with_gil(|py| {
            let m = PyModule::new(py, "microbiorust").unwrap();
            microbiorust(&m).unwrap();
            for func in &["gbk_to_faa", "embl_to_faa", "gbk_to_gff", "embl_to_gff"] {
                assert!(m.getattr(func).is_ok(), "Function {} not found", func);
            }
        });
    }
}
