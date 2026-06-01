//! # microbiorust-py
//!
//! High-performance Python bindings for the microBioRust bioinformatics library using PyO3.
//! This crate provides a bridge for rapid genomic data processing, sequence analysis, and
//! alignment manipulation.
//!
//! ## Setup
//!
//! To install the module into your Python environment:
//! ```bash
//! pip install microbiorust
//! ```
//! If you are on Mac OSX and experience a cc linker error when importing, please add the libpython to the LIBRARY_PATH or DYLD_LIBRARY_PATH
//! for example on conda:
//! ```bash
//! export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH
//! ```!
//! //! ## Collections
//!
//! microbiorust uses typed collections — the return type tells you exactly what data is available:
//!
//! | Type | Returned by | Fields |
//! |---|---|---|
//! | `RecordCollection` | `parse_gbk`, `parse_embl` | contig-level records, `source`, `fna` |
//! | `SequenceCollection` | `record.sequences()` | `locus_tag`, `has_faa`, `has_ffn` |
//! | `FeatureCollection` | `record.features()` | `locus_tag`, `gene`, `product`, `start`, `stop`, `strand` |
//! | `FaaCollection` | `gbk_to_faa`, `embl_to_faa`, `record.faa()` | `locus_tag`, `faa` |
//! | `FfnCollection` | `gbk_to_ffn`, `embl_to_ffn`, `record.ffn()` | `locus_tag`, `ffn` |
//!
//! Reminder: .faa (protein fasta) and .ffn (gene nucleotide fasta), .fna (contig nucleotide fasta)
//! All the above collections support `keys()`, `values()`, `items()`, `__len__`, `__contains__`, and iteration.
//!
//! ## Usage
//! To interact with data in Python access the collection like a dictionary
//!
//! ## Genome Data Access (conversions of filetypes)
//! ### Writing directly to file (most efficient)
//! You can easily directly convert filetypes using write_faa, write_ffn, write_fna
//!
//! import microbiorust as mb
//! collection = mb.parse_gbk("test_input.gbk")
//! collection.write_faa("test_output.faa")
//! collection.write_ffn("test_output.ffn")
//!
//! ### Flat access across the whole genomefile (multiple contigs)
//!
//! ```python
//! import microbiorust
//!
//! # returns FaaCollection
//! faa_collection = microbiorust.gbk_to_faa("genome.gbk")
//!
//! #prints protein fasta format using .values()
//! for info in faa_collection.values():
//!     print(f">{info.locus_tag}\n{info.faa}")
//!
//! # or with items() - useful if you are building a dictionary, doing a lookup or need the key as a separate variable
//! for tag, info in faa_collection.items():
//!     print(f">{tag}\n{info.faa}")
//!
//! # the same pattern for ffn
//! ffn_collection = microbiorust.gbk_to_ffn("genome.gbk")
//! for info in ffn_collection.values():
//!     print(f">{info.locus_tag}\n{info.ffn}")
//! ```
//!
//! parse_gbk return a RecordCollection, one PyRecord per contig which has both sequences() and features()
//!
//! ```python
//! import microbiorust
//!
//! #parse records from a GenBank file
//! collection = microbiorust.parse_gbk("input.gbk")
//!
//! #iterating over the RecordCollection yields keys (usually as 'RecordID|LocusTag')
//! for record in collection.values():
//!     print(f"Processing {record.id()}...")
//!
//!     # per-gene protein sequences — FaaCollection
//!     for info in record.faa().values():
//!         print(f">{info.locus_tag}\n{info.faa}")
//!
//!     # per-gene nucleotide sequences — FfnCollection
//!     for info in record.ffn().values():
//!         print(f">{info.locus_tag}\n{info.ffn}")
//!
//!     # both sequences and features together — SequenceCollection/FeatureCollection
//!     # faa and ffn are Optional here in case not every locus tag has both
//!     sequences = record.sequences()
//!     features = record.features()
//!
//!     if "b3304" in features:
//!         feat = features["b3304"]
//!         print(f"Gene: {feat.gene}, Product: {feat.product}")
//!         print(f"Location: {feat.start}..{feat.stop}, Strand: {feat.strand}")
//!
//!     # whole genome nucleotide sequence
//!     print(f"Genome length: {len(record.sequence())}")
//! ```
//!
//! ### Choosing the best method
//!
//! - **"Write everything to file"** → `collection.write_faa()` / `write_ffn()` / `write_fna()`
//! - **"Get all proteins across a whole genome file"** → `gbk_to_faa()` — flat, one shot
//! - **"Work per genome record (contig-wise)"** → `parse_gbk()` then `record.faa()` or `record.ffn()`
//! - **"Need features and sequences together"** → `parse_gbk()` then `record.sequences()` + `record.features()`
//! - **"Count proteins directly without loading sequences"** → `gbk_to_faa_count()`
//! - **"Export data to web apps, databases, or Pandas"** → `collection.to_json()`
//!
//! //! ## JSON Serialization
//!
//! All core collections (`RecordCollection`, `SequenceCollection`, `FeatureCollection`, etc.)
//! expose a `.to_json()` method to Python. This allows users to instantly export the
//! highly-optimized Rust memory space into standard JSON strings for use in web applications,
//! databases, or native Python data pipelines.
//!
//! ### Example Usage in Python
//!
//! ```python
//! import microbiorust as mb
//! import json
//!
//! # 1. Parse a file into a RecordCollection (lives in Rust memory)
//! collection = mb.parse_gbk("genome.gbk")
//!
//! # 2. Export the entire collection to a JSON string
//! json_str = collection.to_json()
//!
//! # 3. Load into a native Python dictionary for downstream tasks
//! data = json.loads(json_str)
//! ```
//! ## sequencemetrics
//!
//!functions for protein and nucleotide analysis are available at the top level
//!or via the seqmetrics submodule.
//!
//! ```python
//! from microbiorust import amino_percentage
//!
//! seq = "MSNTQKKNVPELRFPGFEGEWEEKKLGDLTTKIGSGKTPKGGSENYTNKG"
//! stats = amino_percentage(seq)
//! print(f"Alanine content: {stats.get('A', 0)}%")
//! ```
//!
//! ## multiple sequence alignment functionality
//!
//!work with Multiple Sequence Alignments (MSA) through the align submodule
//!
//! ```python
//! from microbiorust.align import subset_msa_alignment, get_consensus
//!
//! #subset an alignment (row_start, row_end, col_start, col_end)
//! sub_aln = subset_msa_alignment("input.fasta", 0, 10, 0, 100)
//!
//! #view the consensus
//! sub_aln.get_consensus()
//! ```
//!
//! ## BLAST
//!
//! ```python
//! import microbiorust
//!
//! results = microbiorust.parse_tabular("blast_results.tab")
//! for hit in results:
//!     print(hit["qseqid"], hit["pident"], hit["bitscore"])
//! ```
//! ## Organization
//!
//! Functions are also organized into the following submodules for access:
//! - `microbiorust.gbk`: GenBank loaders and exporters.
//! - `microbiorust.embl`: EMBL loaders and exporters.
//! - `microbiorust.align`: MSA tools and gap purging.
//! - `microbiorust.seqmetrics`: Hydrophobicity, amino acid counts, and percentages.
//! - `microbiorust.blast`: Tabular and XML BLAST parsers.
#![allow(unused_imports)]
#[macro_use]
mod macros;

use microBioRust::align::{Alignment, AlignmentKind, load_msa_auto};
use microBioRust::blast::*;
use microBioRust::embl;
use microBioRust::embl::gff_write as embl_gff_write;
use microBioRust::gbk::{
    FeatureAttributeBuilder, FeatureAttributes, RangeValue, Reader, Record, SequenceAttributes,
    gff_write,
};
use microBioRust::genbank;
use microBioRust_seqmetrics::metrics::amino_counts as rust_amino_counts;
use microBioRust_seqmetrics::metrics::amino_percentage as rust_amino_percentage;
use microBioRust_seqmetrics::metrics::hydrophobicity as rust_hydrophobicity;
use pyo3::exceptions::{PyIOError, PyKeyError};
use pyo3::types::PyList;
use pyo3::{prelude::*, types::PyModule};
use pyo3_stub_gen::define_stub_info_gatherer;
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};
use pyo3_stub_gen::{Result, StubInfo};
use pythonize::pythonize;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::io::BufWriter;
use std::ops::Range;
use std::{
    collections::BTreeMap,
    fs::OpenOptions,
    io::{self, Write},
};
use tokio::io::AsyncBufReadExt;

#[derive(Serialize, Clone)]
pub enum InternalRecord {
    Gbk(microBioRust::gbk::Record),
    Embl(microBioRust::embl::Record),
}

//PyFeatureInfo contains all the feature attributes, there is an extras field for any uncaptured field (e.g. db_xref, EC)
#[gen_stub_pyclass]
#[pyclass]
#[derive(Serialize, Clone)]
pub struct PyFeatureInfo {
    #[pyo3(get)]
    pub locus_tag: String,
    #[pyo3(get)]
    pub gene: Option<String>,
    #[pyo3(get)]
    pub product: Option<String>,
    #[pyo3(get)]
    pub start: Option<u32>,
    #[pyo3(get)]
    pub stop: Option<u32>,
    #[pyo3(get)]
    pub strand: Option<i8>,
    #[pyo3(get)]
    pub codon_start: Option<u8>,
    #[pyo3(get)]
    pub extras: Vec<String>,
}
impl PyFeatureInfo {
    pub fn new(id: &str) -> Self {
        Self {
            locus_tag: id.to_string(),
            gene: None,
            product: None,
            start: None,
            stop: None,
            strand: None,
            codon_start: None,
            extras: Vec::new(),
        }
    }
}
#[gen_stub_pymethods]
#[pymethods]
impl PyFeatureInfo {
    pub fn __repr__(&self) -> String {
        let start = self
            .start
            .map(|v| v.to_string())
            .unwrap_or_else(|| "None".into());
        let stop = self
            .stop
            .map(|v| v.to_string())
            .unwrap_or_else(|| "None".into());
        let strand = self
            .strand
            .map(|v| v.to_string())
            .unwrap_or_else(|| "None".into());
        format!(
            "FeatureInfo(locus_tag='{}', gene={:?}, product={:?}, {}..{}, strand={}, extras={:?})",
            self.locus_tag, self.gene, self.product, start, stop, strand, self.extras
        )
    }
}

//PySequenceInfo contains the sequenceattributes, with an extras field as for PyFeatureInfo
#[gen_stub_pyclass]
#[pyclass]
#[derive(Serialize, Clone)]
pub struct PySequenceInfo {
    #[pyo3(get)]
    pub locus_tag: String,
    #[pyo3(get)]
    pub faa: Option<String>,
    #[pyo3(get)]
    pub ffn: Option<String>,
    #[pyo3(get)]
    pub extras: Vec<String>,
}
impl PySequenceInfo {
    pub fn new(id: &str) -> Self {
        Self {
            locus_tag: id.to_string(),
            faa: None,
            ffn: None,
            extras: Vec::new(),
        }
    }
}
#[gen_stub_pymethods]
#[pymethods]
impl PySequenceInfo {
    fn __repr__(&self) -> String {
        format!(
            "SequenceInfo(locus_tag='{}', has_faa={}, has_ffn={}, extras={:?})",
            self.locus_tag,
            self.faa.is_some(),
            self.ffn.is_some(),
            self.extras,
        )
    }
}
//create an iter struct for Locus_tag key for both FeatureCollection and SequenceCollection
#[gen_stub_pyclass]
#[pyclass]
pub struct LocusTagIterator {
    keys: Vec<String>,
    index: usize,
}
#[gen_stub_pymethods]
#[pymethods]
impl LocusTagIterator {
    fn __iter__(slf: ::pyo3::PyRef<'_, Self>) -> ::pyo3::PyRef<'_, Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<String> {
        if slf.index < slf.keys.len() {
            let res = slf.keys[slf.index].clone();
            slf.index += 1;
            Some(res)
        } else {
            None
        }
    }
}
//create the SequenceCollection and the FeatureCollection to hold the features using macro from macro.rs
//SequenceCollection is a HashMap of Id, Py<PySequenceInfo>, upfront cost - preventing allocations for each entry
//FeatureCollection is a HashMap of Id, Py<PyFeatureInfo>,
//since they are similar except in the data and types contained, both are built here from macro in macro.rs
crate::create_collection!(FeatureCollection, PyFeatureInfo, "FeatureCollection");
crate::create_collection!(SequenceCollection, PySequenceInfo, "SequenceCollection");
//generate the RecordCollection separately with the write functions
#[gen_stub_pyclass]
#[pyclass]
pub struct RecordCollection {
    pub inner: HashMap<String, Py<PyRecord>>,
}
#[gen_stub_pymethods]
#[pymethods]
impl RecordCollection {
    #[new]
    fn new() -> Self {
        RecordCollection {
            inner: HashMap::new(),
        }
    }
    fn __len__(&self) -> usize {
        self.inner.len()
    }
    fn __contains__(&self, id: &str) -> bool {
        self.inner.contains_key(id)
    }
    fn __repr__(&self) -> String {
        format!("RecordCollection({} entries)", self.inner.len())
    }
    fn __getitem__<'py>(&self, record_id: &str, py: Python<'py>) -> PyResult<Bound<'py, PyRecord>> {
        self.inner
            .get(record_id)
            .map(|obj| obj.bind(py).clone())
            .ok_or_else(|| PyKeyError::new_err(record_id.to_string()))
    }
    fn keys<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, self.inner.keys())
    }
    fn __iter__(slf: PyRef<'_, Self>) -> LocusTagIterator {
        LocusTagIterator {
            keys: slf.inner.keys().cloned().collect(),
            index: 0,
        }
    }
    fn insert(&mut self, record_id: String, record: Py<PyRecord>) {
        //TODO: decide if duplicate keys should error on a per line basis
        self.inner.insert(record_id, record);
    }
    //write all the faa format to a file, no conversion to Python overhead
    fn write_faa(&self, filename: &str, py: Python<'_>) -> PyResult<(String, usize)> {
        let mut count = 0;
        let f = std::fs::File::create(filename).map_err(|e| PyIOError::new_err(e.to_string()))?;
        let mut w = BufWriter::new(f);
        for obj in self.inner.values() {
            let record = obj.bind(py).borrow();
            match &record.inner {
                InternalRecord::Gbk(r) => {
                    for (tag, _) in &r.cds.attributes {
                        if let Some(seq) = r.seq_features.get_sequence_faa(tag) {
                            writeln!(w, ">{}|{}\n{}", r.id, tag, seq)
                                .map_err(|e| PyIOError::new_err(e.to_string()))?;
                            count += 1;
                        }
                    }
                }
                InternalRecord::Embl(r) => {
                    for (tag, _) in &r.cds.attributes {
                        if let Some(seq) = r.seq_features.get_sequence_faa(tag) {
                            writeln!(w, ">{}|{}\n{}", r.id, tag, seq)
                                .map_err(|e| PyIOError::new_err(e.to_string()))?;
                            count += 1;
                        }
                    }
                }
            }
        }
        Ok((filename.to_string(), count))
    }
    //write all the ffn to a file, no conversion to python overhead
    fn write_ffn(&self, filename: &str, py: Python<'_>) -> PyResult<(String, usize)> {
        let mut count = 0;
        let f = std::fs::File::create(filename).map_err(|e| PyIOError::new_err(e.to_string()))?;
        let mut w = BufWriter::new(f);
        for obj in self.inner.values() {
            let record = obj.bind(py).borrow();
            match &record.inner {
                InternalRecord::Gbk(r) => {
                    for (tag, _) in &r.cds.attributes {
                        if let Some(seq) = r.seq_features.get_sequence_ffn(tag) {
                            writeln!(w, ">{}|{}\n{}", r.id, tag, seq)
                                .map_err(|e| PyIOError::new_err(e.to_string()))?;
                            count += 1;
                        }
                    }
                }
                InternalRecord::Embl(r) => {
                    for (tag, _) in &r.cds.attributes {
                        if let Some(seq) = r.seq_features.get_sequence_ffn(tag) {
                            writeln!(w, ">{}|{}\n{}", r.id, tag, seq)
                                .map_err(|e| PyIOError::new_err(e.to_string()))?;
                            count += 1;
                        }
                    }
                }
            }
        }
        Ok((filename.to_string(), count))
    }
    //write fna of records to a file, no conversion to python overhead
    fn write_fna(&self, filename: &str, py: Python<'_>) -> PyResult<(String, usize)> {
        let mut count = 0;
        let f = std::fs::File::create(filename).map_err(|e| PyIOError::new_err(e.to_string()))?;
        let mut w = BufWriter::new(f);
        for obj in self.inner.values() {
            let record = obj.bind(py).borrow();
            match &record.inner {
                InternalRecord::Gbk(r) => writeln!(w, ">{}\n{}", r.id, r.sequence)
                    .map_err(|e| PyIOError::new_err(e.to_string()))?,
                InternalRecord::Embl(r) => writeln!(w, ">{}\n{}", r.id, r.sequence)
                    .map_err(|e| PyIOError::new_err(e.to_string()))?,
            }
            count += 1;
        }
        Ok((filename.to_string(), count))
    }
    fn values<'py>(
        &self,
        py: ::pyo3::Python<'py>,
    ) -> ::pyo3::PyResult<::pyo3::Bound<'py, ::pyo3::types::PyList>> {
        let values: Vec<_> = self
            .inner
            .values()
            .map(|obj| obj.bind(py).clone())
            .collect();
        ::pyo3::types::PyList::new(py, values)
    }

    fn items<'py>(
        &self,
        py: ::pyo3::Python<'py>,
    ) -> ::pyo3::PyResult<::pyo3::Bound<'py, ::pyo3::types::PyList>> {
        let items: Vec<::pyo3::Bound<'py, ::pyo3::types::PyTuple>> = self
            .inner
            .iter()
            .map(|(k, v)| {
                ::pyo3::types::PyTuple::new(
                    py,
                    &[
                        k.clone().into_pyobject(py).unwrap().into_any(),
                        v.bind(py).clone().into_any(),
                    ],
                )
            })
            .collect::<::pyo3::PyResult<Vec<_>>>()?;
        ::pyo3::types::PyList::new(py, items)
    }
    /// Serializes the RecordCollection to a formatted JSON string
    pub fn to_json<'py>(&self, py: Python<'py>) -> PyResult<String> {
        let mut export_map = std::collections::HashMap::new();

        for (key, py_ptr) in &self.inner {
            let rust_struct = py_ptr.bind(py).borrow().clone();
            export_map.insert(key.clone(), rust_struct);
        }
        serde_json::to_string_pretty(&export_map).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("JSON serialization failed: {}", e))
        })
    }
}
//ability to return a HashMap to python containing all the data
//create build sequences and build features methods for both gbk and embl types
//TODO: replace with a single method for each if we move from gbk and embl types to a generic
// Generate the two functions
crate::impl_build_sequences!(
    build_sequences_gbk,
    microBioRust::gbk::Record,
    microBioRust::gbk::SequenceAttributes
);
crate::impl_build_sequences!(
    build_sequences_embl,
    microBioRust::embl::Record,
    microBioRust::embl::SequenceAttributes
);
crate::impl_build_features!(
    build_features_gbk,
    microBioRust::gbk::Record,
    microBioRust::gbk::FeatureAttributes
);
crate::impl_build_features!(
    build_features_embl,
    microBioRust::embl::Record,
    microBioRust::embl::FeatureAttributes
);
create_collection!(FaaCollection, PyFaaInfo, "FaaCollection");
create_collection!(FfnCollection, PyFfnInfo, "FfnCollection");
//create a PyRecord with inner of InternalRecord for gbk and embl, contains .sequences() and .features() methods
//TODO: replace the InternalRecord with a top level Record if we move from gbk/embl to a generic type
//and all match arms collapse to single call
#[gen_stub_pyclass]
#[pyclass]
#[derive(Serialize, Clone)]
pub struct PyRecord {
    inner: InternalRecord,
}
#[gen_stub_pymethods]
#[pymethods]
impl PyRecord {
    fn id(&self) -> &str {
        match &self.inner {
            InternalRecord::Gbk(r) => &r.id,
            InternalRecord::Embl(r) => &r.id,
        }
    }
    //whole-genome nucleotide — distinct from per-CDS ffn
    fn sequence(&self) -> &str {
        match &self.inner {
            InternalRecord::Gbk(r) => &r.sequence,
            InternalRecord::Embl(r) => &r.sequence,
        }
    }
    fn locus_tag(&self) -> Vec<String> {
        match &self.inner {
            InternalRecord::Gbk(r) => r.cds.attributes.keys().cloned().collect(),
            InternalRecord::Embl(r) => r.cds.attributes.keys().cloned().collect(),
        }
    }
    // py token required — Python objects allocated inside build helpers
    fn sequences(&self, py: Python<'_>) -> SequenceCollection {
        let inner = match &self.inner {
            InternalRecord::Gbk(r) => build_sequences_gbk(r, py),
            InternalRecord::Embl(r) => build_sequences_embl(r, py),
        };
        SequenceCollection { inner }
    }
    fn features(&self, py: Python<'_>) -> FeatureCollection {
        let inner = match &self.inner {
            InternalRecord::Gbk(r) => build_features_gbk(r, py),
            InternalRecord::Embl(r) => build_features_embl(r, py),
        };
        FeatureCollection { inner }
    }
    fn faa(&self, py: Python<'_>) -> FaaCollection {
        let inner = match &self.inner {
            InternalRecord::Gbk(r) => build_faa_from_gbk(r, py),
            InternalRecord::Embl(r) => build_faa_from_embl(r, py),
        };
        FaaCollection { inner }
    }
    fn ffn(&self, py: Python<'_>) -> FfnCollection {
        let inner = match &self.inner {
            InternalRecord::Gbk(r) => build_ffn_from_gbk(r, py),
            InternalRecord::Embl(r) => build_ffn_from_embl(r, py),
        };
        FfnCollection { inner }
    }
    fn __repr__(&self) -> String {
        match &self.inner {
            InternalRecord::Gbk(r) => format!(
                "PyRecord(id='{}', format=GBK,  seq_len={})",
                r.id,
                r.sequence.len()
            ),
            InternalRecord::Embl(r) => format!(
                "PyRecord(id='{}', format=EMBL, seq_len={})",
                r.id,
                r.sequence.len()
            ),
        }
    }
    fn __getitem__(&self, py: Python<'_>, tag: &str) -> PyResult<Py<PyAny>> {
        let features = self.features(py);
        if features.inner.contains_key(tag) {
            // 1. Get the Bound<'_, PyFeatureInfo>
            // 2. Use .into_any() to get Bound<'_, PyAny>
            // 3. Use .unbind() to get PyObject
            return features.__getitem__(tag, py).map(|b| b.into_any().unbind());
        }

        let sequences = self.sequences(py);
        if sequences.inner.contains_key(tag) {
            return sequences
                .__getitem__(tag, py)
                .map(|b| b.into_any().unbind());
        }

        Err(PyKeyError::new_err(tag.to_string()))
    }
}
//parse genbank function - returns the RecordCollection, good for fna and metadata from source
#[gen_stub_pyfunction]
#[pyfunction]
pub fn parse_gbk(filename: &str, py: Python<'_>) -> PyResult<RecordCollection> {
    let mut inner = HashMap::new();
    for r in genbank!(filename) {
        let id = r.id.clone();
        //TODO: STRICT -Duplicate Check and line number error
        //if inner.contains_key(&id) {
        //    return Err(PyValueError::new_err(format!(
        //       "Integrity Error: Duplicate record ID '{}' encountered in file '{}'. \
        //         Check line {} (approx). Bioinformatic formats require unique identifiers.",
        //        id,
        //        filename,
        //        r.start_line // Assumes microBioRust tracks the line number
        //    )));
        //}
        let obj = Py::new(
            py,
            PyRecord {
                inner: InternalRecord::Gbk(r),
            },
        )
        .expect("Failed to allocate PyRecord");
        inner.insert(id, obj);
    }
    Ok(RecordCollection { inner })
}
//same function as above but for embl
#[gen_stub_pyfunction]
#[pyfunction]
pub fn parse_embl(filename: &str, py: Python<'_>) -> PyResult<RecordCollection> {
    let mut inner = HashMap::new();
    for r in embl!(filename) {
        let id = r.id.clone();
        //TODO: STRICT -Duplicate Check and line number error
        //if inner.contains_key(&id) {
        //    return Err(PyValueError::new_err(format!(
        //       "Integrity Error: Duplicate record ID '{}' encountered in file '{}'. \
        //         Check line {} (approx). Bioinformatic formats require unique identifiers.",
        //        id,
        //        filename,
        //        r.start_line // Assumes microBioRust tracks the line number
        //    )));
        //}
        let obj = Py::new(
            py,
            PyRecord {
                inner: InternalRecord::Embl(r),
            },
        )
        .expect("Failed to allocate PyRecord");
        inner.insert(id, obj);
    }
    Ok(RecordCollection { inner })
}

//the fna is the whole record nucleotide sequence
//python callers use record.sequence() to get the string.
#[gen_stub_pyfunction]
#[pyfunction]
pub fn gbk_to_fna(filename: &str, py: Python<'_>) -> PyResult<RecordCollection> {
    parse_gbk(filename, py)
}
#[gen_stub_pyfunction]
#[pyfunction]
pub fn embl_to_fna(filename: &str, py: Python<'_>) -> PyResult<RecordCollection> {
    parse_embl(filename, py)
}

//to return a flat SequenceCollection keyed by "record_id|locus_tag".
//need to use parse_gbk/parse_embl to access both sequences and features
#[gen_stub_pyfunction]
#[pyfunction]
pub fn gbk_to_faa(filename: &str, py: Python<'_>) -> PyResult<FaaCollection> {
    let mut inner = HashMap::new();
    for r in genbank!(filename) {
        for (tag, _) in &r.cds.attributes {
            if let Some(seq) = r.seq_features.get_sequence_faa(tag) {
                let key = format!("{}|{}", r.id, tag);
                let info = PyFaaInfo {
                    locus_tag: key.clone(),
                    faa: seq.to_string(),
                };
                let obj = Py::new(py, info).expect("failed to allocate PyFaaInfo");
                inner.insert(key, obj);
            }
        }
    }
    Ok(FaaCollection { inner })
}
#[gen_stub_pyfunction]
#[pyfunction]
pub fn gbk_to_ffn(filename: &str, py: Python<'_>) -> PyResult<FfnCollection> {
    let mut inner = HashMap::new();
    for r in genbank!(filename) {
        for (tag, _) in &r.cds.attributes {
            if let Some(seq) = r.seq_features.get_sequence_ffn(tag) {
                let key = format!("{}|{}", r.id, tag);
                let info = PyFfnInfo {
                    locus_tag: key.clone(),
                    ffn: seq.to_string(),
                };
                let obj = Py::new(py, info).expect("failed to allocate PyFfnInfo");
                inner.insert(key, obj);
            }
        }
    }
    Ok(FfnCollection { inner })
}
#[gen_stub_pyclass]
#[pyclass]
#[derive(Serialize, Clone)]
pub struct PyFaaInfo {
    #[pyo3(get)]
    pub locus_tag: String,
    #[pyo3(get)]
    pub faa: String,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyFaaInfo {
    fn __repr__(&self) -> String {
        format!(
            "FaaInfo(locus_tag='{}', faa='{}...')",
            self.locus_tag,
            &self.faa[..20.min(self.faa.len())]
        )
    }
}
crate::impl_build_faa!(
    build_faa_from_gbk,
    microBioRust::gbk::Record,
    microBioRust::gbk::SequenceAttributes
);
crate::impl_build_faa!(
    build_faa_from_embl,
    microBioRust::embl::Record,
    microBioRust::embl::SequenceAttributes
);

#[gen_stub_pyclass]
#[pyclass]
#[derive(Serialize, Clone)]
pub struct PyFfnInfo {
    #[pyo3(get)]
    pub locus_tag: String,
    #[pyo3(get)]
    pub ffn: String,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyFfnInfo {
    fn __repr__(&self) -> String {
        format!(
            "FfnInfo(locus_tag='{}', ffn='{}...')",
            self.locus_tag,
            &self.ffn[..20.min(self.ffn.len())]
        )
    }
}
crate::impl_build_ffn!(
    build_ffn_from_gbk,
    microBioRust::gbk::Record,
    microBioRust::gbk::SequenceAttributes
);
crate::impl_build_ffn!(
    build_ffn_from_embl,
    microBioRust::embl::Record,
    microBioRust::embl::SequenceAttributes
);
#[gen_stub_pyfunction]
#[pyfunction]
pub fn embl_to_faa(filename: &str, py: Python<'_>) -> PyResult<FaaCollection> {
    let mut inner = HashMap::new();
    for r in embl!(filename) {
        for (tag, _) in &r.cds.attributes {
            if let Some(seq) = r.seq_features.get_sequence_faa(tag) {
                let key = format!("{}|{}", r.id, tag);
                let info = PyFaaInfo {
                    locus_tag: key.clone(),
                    faa: seq.to_string(),
                };
                let obj = Py::new(py, info).expect("failed to allocate PyFaaInfo");
                inner.insert(key, obj);
            }
        }
    }
    Ok(FaaCollection { inner })
}
#[gen_stub_pyfunction]
#[pyfunction]
pub fn embl_to_ffn(filename: &str, py: Python<'_>) -> PyResult<FfnCollection> {
    let mut inner = HashMap::new();
    for r in embl!(filename) {
        for (tag, _) in &r.cds.attributes {
            if let Some(seq) = r.seq_features.get_sequence_ffn(tag) {
                let key = format!("{}|{}", r.id, tag);
                let info = PyFfnInfo {
                    locus_tag: key.clone(),
                    ffn: seq.to_string(),
                };
                let obj = Py::new(py, info).expect("failed to allocate PyFfnInfo");
                inner.insert(key, obj);
            }
        }
    }
    Ok(FfnCollection { inner })
}

//count function
#[gen_stub_pyfunction]
#[pyfunction]
pub fn gbk_to_faa_count(filename: &str) -> PyResult<usize> {
    Ok(genbank!(filename)
        .iter()
        .flat_map(|r| r.seq_features.seq_attributes.values())
        .flat_map(|set| set.iter())
        .filter(|attr| matches!(attr, SequenceAttributes::SequenceFaa { .. }))
        .count())
}

pub fn process_alignment<T>(
    aln: Alignment<T>,
    rows: Range<usize>,
    cols: Range<usize>,
) -> Vec<String> {
    let sub = aln.subset(rows.clone(), cols.clone());
    let sequences: Vec<String> = sub
        .sequence_data
        .chunks_exact(sub.cols)
        .map(|chunk| String::from_utf8_lossy(chunk).into_owned())
        .collect();

    sub.ids
        .iter()
        .zip(sequences.iter())
        .map(|(id, seq)| format!(">{}\n{}", id, seq))
        .collect()
}
#[gen_stub_pyfunction]
#[pyfunction]
pub fn subset_msa_alignment(
    filename: &str,
    row_index: (usize, usize),
    col_index: (usize, usize),
) -> PyResult<Vec<String>> {
    let records = load_msa_auto(&filename)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;
    let subrows = row_index.0..row_index.1;
    let subcols = col_index.0..col_index.1;

    let result = match records {
        AlignmentKind::DNA(aln) => process_alignment(aln, subrows, subcols),
        AlignmentKind::Protein(aln) => process_alignment(aln, subrows, subcols),
    };

    Ok(result)
}

#[allow(unused_variables)]
#[gen_stub_pyfunction]
#[pyfunction]
pub fn gbk_to_gff(filename: &str, dna: bool) -> PyResult<()> {
    let records = genbank!(&filename);
    //let prev_start: u32 = 0;
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
            prev_end += ending; //this is to create the joined record if there are multiple
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
    Ok(())
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[gen_stub_pyfunction]
#[pyfunction]
pub fn embl_to_gff(filename: &str, dna: bool) -> PyResult<()> {
    let records = embl!(&filename);
    //let prev_start: u32 = 0;
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
            prev_end += ending; //this is to create the joined record if there are multiple
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
    Ok(())
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[gen_stub_pyfunction]
#[pyfunction]
fn hydrophobicity(seq: &str, window_size: usize) -> Vec<f64> {
    rust_hydrophobicity(seq, window_size)
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[gen_stub_pyfunction]
#[pyfunction]
fn amino_percentage(seq: &str) -> HashMap<char, f64> {
    rust_amino_percentage(seq)
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[gen_stub_pyfunction]
#[pyfunction]
fn amino_counts(seq: &str) -> HashMap<char, u64> {
    rust_amino_counts(seq)
}

//bridge function for GAP PURGING
#[gen_stub_pyfunction]
#[pyfunction]
pub fn purge_gaps(input: &str, output: &str, threshold: f32) -> PyResult<()> {
    let kind = load_msa_auto(input)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;

    match kind {
        AlignmentKind::DNA(mut aln) => {
            aln.purge_gappy_columns(threshold);
            aln.write_fasta(output, Some(60))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;
        }
        AlignmentKind::Protein(mut aln) => {
            aln.purge_gappy_columns(threshold);
            aln.write_fasta(output, Some(60))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;
        }
    }
    Ok(())
}
#[gen_stub_pyfunction]
#[pyfunction]
pub fn get_consensus(filename: &str) -> PyResult<String> {
    let kind = load_msa_auto(filename)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;

    if let AlignmentKind::DNA(aln) = kind {
        Ok(aln.consensus_sequence())
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
            "Consensus is only available for DNA alignments",
        ))
    }
}

#[allow(unused_imports)]
#[allow(unused_variables)]
#[gen_stub_pyfunction]
#[pyfunction]
pub fn parse_tabular<'py>(py: Python<'py>, path: String) -> PyResult<Bound<'py, PyAny>> {
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;
    rt.block_on(async {
        let reader = open_async_reader(&path)
            .await
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        let mut results = Vec::new();
        let mut lines = reader.lines();
        while let Some(line) = lines.next_line().await.unwrap_or(None) {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') {
                continue;
            }
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
        pythonize(py, &results)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}
#[allow(unused_imports)]
#[allow(unused_variables)]
#[gen_stub_pyfunction]
#[pyfunction]
pub fn parse_xml<'py>(py: Python<'py>, path: String) -> PyResult<Bound<'py, PyAny>> {
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;

    rt.block_on(async {
        let reader = open_async_reader(&path)
            .await
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

        let mut iterations = Vec::new();
        let mut iter = AsyncBlastXmlIter::from_reader(reader);

        while let Some(res) = iter.next_iteration().await {
            if let Ok(val) = res {
                iterations.push(val);
            }
        }
        //convert nested XML structs into Python Dicts
        pythonize(py, &iterations)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

#[pyfunction]
pub fn register_seqmetrics(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "seqmetrics")?;
    register_all!(m, parent, hydrophobicity, amino_counts, amino_percentage);
    parent.add_submodule(&m)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("microbiorust.seqmetrics", &m)?;

    Ok(())
}

#[pyfunction]
pub fn register_align(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "align")?;
    register_all!(m, parent, purge_gaps, subset_msa_alignment, get_consensus);
    parent.add_submodule(&m)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("microbiorust.align", &m)?;

    Ok(())
}

#[pyfunction]
pub fn register_blast(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "blast")?;
    register_all!(m, parent, parse_tabular, parse_xml);
    parent.add_submodule(&m)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("microbiorust.blast", &m)?;

    Ok(())
}

#[pyfunction]
pub fn register_gbk(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "gbk")?;
    register_all!(
        m,
        parent,
        parse_gbk,
        gbk_to_faa,
        gbk_to_fna,
        gbk_to_ffn,
        gbk_to_faa_count,
        gbk_to_gff
    );
    parent.add_submodule(&m)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("microbiorust.gbk", &m)?;

    Ok(())
}

#[pyfunction]
pub fn register_embl(py: Python<'_>, parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(py, "embl")?;
    register_all!(
        m,
        parent,
        parse_embl,
        embl_to_faa,
        embl_to_fna,
        embl_to_ffn,
        embl_to_gff
    );
    parent.add_submodule(&m)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("microbiorust.embl", &m)?;

    Ok(())
}

#[pymodule]
fn _microbiorust(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_gbk(py, m)?;
    register_embl(py, m)?;
    register_align(py, m)?;
    register_seqmetrics(py, m)?;
    register_blast(py, m)?;
    //register classes
    m.add_class::<PyRecord>()?;
    m.add_class::<PyFeatureInfo>()?;
    m.add_class::<PySequenceInfo>()?;

    //register Collections and Iterators!
    m.add_class::<RecordCollection>()?;
    m.add_class::<FeatureCollection>()?;
    m.add_class::<SequenceCollection>()?;
    m.add_class::<LocusTagIterator>()?;
    Ok(())
}
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);
pub fn get_stub_info() -> pyo3_stub_gen::Result<StubInfo> {
    stub_info()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::_microbiorust;
    use pyo3::prelude::*;
    use pyo3::types::PyModule;
    #[test]
    #[allow(unused_imports)]
    fn test_functions_are_registered() {
        pyo3::Python::initialize();
        Python::attach(|py| {
            let m = PyModule::new(py, "microbiorust").unwrap();
            _microbiorust(py, &m).unwrap();

            //GBK
            let gbk = m.getattr("gbk").expect("gbk submodule missing");
            for func in &[
                "gbk_to_faa",
                "parse_gbk",
                "gbk_to_fna",
                "gbk_to_ffn",
                "gbk_to_faa_count",
                "gbk_to_gff",
            ] {
                assert!(gbk.getattr(func).is_ok(), "Function gbk.{} not found", func);
            }

            //EMBL
            let embl = m.getattr("embl").expect("embl submodule missing");
            for func in &[
                "embl_to_faa",
                "parse_embl",
                "embl_to_ffn",
                "embl_to_fna",
                "embl_to_gff",
            ] {
                assert!(
                    embl.getattr(func).is_ok(),
                    "Function embl.{} not found",
                    func
                );
            }

            //seqmetrics
            let seqmetrics = m
                .getattr("seqmetrics")
                .expect("seqmetrics submodule missing");
            for func in &["hydrophobicity", "amino_counts", "amino_percentage"] {
                assert!(
                    seqmetrics.getattr(func).is_ok(),
                    "Function seqmetrics.{} not found",
                    func
                );
            }

            //align
            let align = m.getattr("align").expect("align submodule missing");
            assert!(
                align.getattr("subset_msa_alignment").is_ok(),
                "Function align.subset_msa_alignment not found"
            );
            //blast
            let blast = m.getattr("blast").expect("blast submodule missing");
            assert!(
                blast.getattr("parse_tabular").is_ok(),
                "Function blast.parse_tabular not found"
            );
        });
    }
}
