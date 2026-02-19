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
//!  from microbiorust import gbk_to_faa
//!  result = gbk_to_faa("test_input.gbk")
//!  for r in result:
//!      print(r)
//!  gbk_to_gff("test_input.gbk")
//!  
//!  Other pyfunctions that can be run include gbk_to_faa, embl_to_faa, gbk_to_gff, embl_to_gff, amino_counts, amino_percentage, hydrophobicity
//!
//!  from microbiorust import amino_percentage
//!  result = amino_percentage("MSNTQKKNVPELRFPGFEGEWEEKKLGDLTTKIGSGKTPKGGSENYTNKGIPFLRSQNIRNGKLNLNDLVYISKDIDDEMKNSRTY")
//!
//!  print(result)
//!
//!  Example for use with a Multiple Sequence Alignment (load in fasta format)
//!
//!  from microbiorust import load_msa_auto
//!
//!  result = load_msa_auto("test_alignment.aln")
//!  //to subset the alignment pass the row and column indices
//!  sub_alignment = result.subset(10,30)
//!  //to save to file in usual clustal type format
//!  sub_alignment.display_interleaved()
//!
//!
#![allow(unused_imports)]
#[macro_use]
mod macros;

use pyo3::{
   prelude::*,
   types::PyModule,
};
use pythonize::pythonize;
use tokio::io::AsyncBufReadExt;
use std::collections::HashMap;
use std::{
   collections::BTreeMap,
   io::{self, Write},
   fs::OpenOptions,
};
use microBioRust::genbank;
use microBioRust::gbk::{Record, Reader, RangeValue, gff_write};
use microBioRust::embl;
use std::ops::Range;
use microBioRust::embl::gff_write as embl_gff_write;
use microBioRust::blast::*;
use microBioRust::align::{Alignment, AlignmentKind, load_msa_auto};
use microBioRust_seqmetrics::metrics::hydrophobicity as rust_hydrophobicity;
use microBioRust_seqmetrics::metrics::amino_counts as rust_amino_counts;
use microBioRust_seqmetrics::metrics::amino_percentage as rust_amino_percentage;


#[pyfunction]
pub fn gbk_to_faa(filename: &str) -> PyResult<Vec<String>> {  
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
pub fn gbk_to_faa_count(filename: &str) -> PyResult<usize> {
    let records = genbank!(filename);
    let mut count = 0;
    for record in records {
        for (_k, _) in &record.cds.attributes {
            count += 1;
        }
    }
    Ok(count)
}
#[pyfunction]
pub fn gbk_to_fna(filename: &str) -> PyResult<Vec<String>> {  
    let records = genbank!(&filename);
    let mut result = Vec::new();
    for record in records {
        result.push(format!(">{}\n{}", &record.id, &record.sequence));
    }
    Ok(result)
}
#[pyfunction]
pub fn gbk_to_ffn(filename: &str) -> PyResult<Vec<String>> {  
    let records = genbank!(&filename);
    let mut result = Vec::new();
    for record in records {
        for (k, _v) in &record.cds.attributes {
            if let Some(seq) = record.seq_features.get_sequence_ffn(k) {
                result.push(format!(">{}|{}\n{}", &record.id, &k, seq));
            }
        }
    }
    Ok(result)
}

#[pyfunction]
pub fn embl_to_faa(filename: &str) -> PyResult<Vec<String>> {
    let records = embl!(&filename);
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
pub fn embl_to_fna(filename: &str) -> PyResult<Vec<String>> {
    let records = embl!(&filename);
    let mut result = Vec::new();
    for record in records {
        result.push(format!(">{}\n{}", &record.id, &record.sequence));
        }
    Ok(result)
}

pub fn process_alignment<T>(aln: Alignment<T>, rows: Range<usize>, cols: Range<usize>) -> Vec<String> {
    let sub = aln.subset(rows.clone(), cols.clone());
    let sequences: Vec<String> = sub.sequence_data.chunks_exact(sub.cols).map(|chunk| String::from_utf8_lossy(chunk).into_owned()).collect();
    
    sub.ids
        .iter()
        .zip(sequences.iter())
        .map(|(id, seq)| format!(">{}\n{}", id, seq))
        .collect()
}

#[pyfunction]
pub fn subset_msa_alignment(filename: &str, row_index: (usize,usize), col_index: (usize,usize)) -> PyResult<Vec<String>> {
    let records = load_msa_auto(&filename).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError,_>(format!("{}", e)))?;
    let subrows = row_index.0..row_index.1;
    let subcols = col_index.0..col_index.1;

    let result = match records {
          AlignmentKind::DNA(aln) => process_alignment(aln, subrows, subcols),
	  AlignmentKind::Protein(aln) => process_alignment(aln, subrows, subcols),
	  };
	         
    Ok(result)
}

#[allow(unused_variables)]
#[pyfunction]
pub fn gbk_to_gff(filename: &str, dna: bool) -> PyResult<()> {
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
pub fn embl_to_gff(filename: &str, dna: bool) -> PyResult<()> {
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

#[allow(unused_imports)]
#[allow(unused_variables)]
#[pyfunction]
fn hydrophobicity(seq: &str, window_size: usize) -> Vec<f64> {
    rust_hydrophobicity(seq, window_size)
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[pyfunction]
fn amino_percentage(seq: &str) -> HashMap<char,f64> {
    rust_amino_percentage(seq)
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[pyfunction]
fn amino_counts(seq: &str) -> HashMap<char, u64> {
    rust_amino_counts(seq)
}

// Bridge function for GAP PURGING
#[pyfunction]
pub fn purge_gaps(input: &str, output: &str, threshold: f32) -> PyResult<()> {
    let kind = load_msa_auto(input).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;
    
    match kind {
        AlignmentKind::DNA(mut aln) => {
            aln.purge_gappy_columns(threshold);
            aln.write_fasta(output, Some(60)).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;
        },
        AlignmentKind::Protein(mut aln) => {
            aln.purge_gappy_columns(threshold);
            aln.write_fasta(output, Some(60)).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;
        },
    }
    Ok(())
}

#[pyfunction]
pub fn get_consensus(filename: &str) -> PyResult<String> {
    let kind = load_msa_auto(filename).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;
    
    if let AlignmentKind::DNA(aln) = kind {
        Ok(aln.consensus_sequence())
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>("Consensus is only available for DNA alignments"))
    }
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[pyfunction]
pub fn parse_tabular<'py>(py: Python<'py>, path: String) -> PyResult<Bound<'py, PyAny>> {
    let rt = tokio::runtime::Builder::new_current_thread()
         .enable_all()
	 .build()
	 .unwrap();
    
    rt.block_on(async {
        let reader = open_async_reader(&path).await
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
        let mut results = Vec::new();
        let mut lines = reader.lines();
        
        while let Some(line) = lines.next_line().await.unwrap_or(None) {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            let cols: Vec<&str> = t.split('\t').collect();
            if cols.len() >= 12 {
                results.push(BlastTabRecord {
                    qseqid: cols[0].to_string(),
                    sseqid: cols[1].to_string(),
                    pident: cols[2].parse().unwrap_or(0.0),
		    length: cols[3].parse().unwrap_or(0),
		    mismatch: cols[4].parse::<u32>().ok(),
		    gapopen: cols[5].parse::<u32>().ok(),
		    qstart: cols[6].parse::<u32>().ok(),
		    qend: cols[7].parse::<u32>().ok(),
		    sstart: cols[8].parse::<u32>().ok(),
		    send: cols[9].parse::<u32>().ok(),
		    evalue: cols[10].parse().unwrap_or(0.0),
                    bitscore: cols[11].parse().unwrap_or(0.0),
                });
            }
        }
        
        //convert the Vec<BlastTabRecord> into a Python List of Dicts
        pythonize(py, &results).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}
#[allow(unused_imports)]
#[allow(unused_variables)]
#[pyfunction]
pub fn parse_xml<'py>(py: Python<'py>, path: String) -> PyResult<Bound<'py, PyAny>> {
    let rt = tokio::runtime::Builder::new_current_thread()
            .enable_all()
	    .build()
	    .unwrap();
    
    rt.block_on(async {
        let reader = open_async_reader(&path).await
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
        let mut iterations = Vec::new();
        let mut iter = AsyncBlastXmlIter::from_reader(reader);
        
        while let Some(res) = iter.next_iteration().await {
            if let Ok(val) = res {
                iterations.push(val);
            }
        }
        
        //convert nested XML structs into Python Dicts
        pythonize(py, &iterations).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}
#[pyfunction]
pub fn register_seqmetrics(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "seqmetrics")?;
    register_functions!(m, hydrophobicity, amino_counts, amino_percentage);
    parent.add_submodule(&m)?;
    py.import("sys")?
      .getattr("modules")?
      .set_item("microbiorust.seqmetrics", &m)?;
    
    Ok(())
}
#[pyfunction]
pub fn register_align(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "align")?;
    register_functions!(m, subset_msa_alignment, get_consensus);
    parent.add_submodule(&m)?;
    py.import("sys")?
      .getattr("modules")?
      .set_item("microbiorust.align", &m)?;
    
    Ok(())
}
#[pyfunction]
pub fn register_blast(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "blast")?;
    register_functions!(m, parse_tabular, parse_xml);
    parent.add_submodule(&m)?;
    py.import("sys")?
      .getattr("modules")?
      .set_item("microbiorust.blast", &m)?;
    
    Ok(())
}
#[pyfunction]
pub fn register_gbk(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "gbk")?;
    register_functions!(m, gbk_to_faa, gbk_to_fna, gbk_to_ffn, gbk_to_faa_count, gbk_to_gff);
    parent.add_submodule(&m)?;
    py.import("sys")?
      .getattr("modules")?
      .set_item("microbiorust.gbk", &m)?;
    
    Ok(())
}
#[pyfunction]
pub fn register_embl(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "embl")?;
    //register_functions!(m, embl_to_gff, embl_to_faa, embl_to_fna);
    m.add_function(pyo3::wrap_pyfunction!(embl_to_gff, &m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(embl_to_faa, &m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(embl_to_fna, &m)?)?;
    parent.add_submodule(&m)?;
    py.import("sys")?
      .getattr("modules")?
      .set_item("microbiorust.embl", &m)?;
    
    Ok(())
}


#[pymodule]
fn microbiorust(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_gbk(py, m)?;
    register_embl(py, m)?;
    register_align(py, m)?;
    register_seqmetrics(py, m)?;

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::microbiorust;
    use pyo3::prelude::*;
    use pyo3::types::PyModule;
    #[test]
    #[allow(unused_imports)]
    fn test_functions_are_registered() {
        pyo3::prepare_freethreaded_python();
        Python::with_gil(|py| {
            let m = PyModule::new(py, "microbiorust").unwrap();
            microbiorust(py, &m).unwrap();
 
            // GBK
            let gbk = m.getattr("gbk").expect("gbk submodule missing");
            for func in &["gbk_to_faa", "gbk_to_fna", "gbk_to_faa_count", "gbk_to_gff"] {
                assert!(gbk.getattr(func).is_ok(), "Function gbk.{} not found", func);
            }

            // EMBL
            let embl = m.getattr("embl").expect("embl submodule missing");
            for func in &["embl_to_faa", "embl_to_fna", "embl_to_gff"] {
                assert!(embl.getattr(func).is_ok(), "Function embl.{} not found", func);
            }

            // seqmetrics
            let seqmetrics = m.getattr("seqmetrics").expect("seqmetrics submodule missing");
            for func in &["hydrophobicity", "amino_counts", "amino_percentage"] {
                assert!(
                    seqmetrics.getattr(func).is_ok(),
                    "Function seqmetrics.{} not found",
                    func
                );
            }

            // align
            let align = m.getattr("align").expect("align submodule missing");
            assert!(
                align.getattr("subset_msa_alignment").is_ok(),
                "Function align.subset_msa_alignment not found"
            );
        });
    }

}
