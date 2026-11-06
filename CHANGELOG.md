# Changelog #
---

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project should adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

## [Unreleased] 06/11/26

### Added
- A `FeatureKind` enum for GenBank/EMBL feature keys to hold the feature type (`gene`, `CDS`, `tRNA`, `rRNA`,
  `tmRNA`, `ncRNA`, `misc_RNA`, `regulatory`, `repeat_region`, `mobile_element`,
  `misc_feature`, other), with `from_key`, and `is_non_coding` to access the data
- Non-coding features (tRNA, rRNA, ncRNA and similar) are now parsed from GenBank and
  EMBL records into `Record::nc_features` (`NonCodingAttributes`, and
  `NonCodingAttributeBuilder`). Joined locations (for joined feature) are now stored as one entry per segment.
- `feature_type` attribute is now added on CDS features.
- `feature_key` and `embl_feature_key` helpers for detecting feature header lines, cleaner than before listing each 
  possibel option in turn
- Test added for joined feature locations

### Fixed
- Parser hang on records with no CDS after first implementing the non coding feature parsing.  
  The source qualifier loop only stopped at a CDS line, so a feature-less record could loop forever in a multi-record file,
  or run on into the next record.
- The CDS qualifier loop now does not read through a following non-coding feature, which
  could possibly overwrite the CDS's product or locus tag. It stops at any feature header.
- Records with a `CONTIG` line followed by `ORIGIN` are no longer rejected as master
  records, this type of output comes from PGAP annotation pipeline of NCBI.

### Changed
- Accessors on `Record` (`cds()`, `source_map()`, `non_coding()`, `seq_features()`)
  return references instead of clones

### Internal
- Clippy: removed the inherent `default()` from `create_builder!` (builders derive
  `Default`), derived `Default` for `FeatureKind`, replaced with `or_default()`, and 
   tidied `Option` handling in the sequence loop.

## Unreleased ##
---

### Added ###
- docs
- windows install description for docs
- this changelog

### Changed ###
- Ran cargo update and removed unused crates futures and thiserror
- Replaced AttributeBuilder<K, HashSet<V>> in record.rs with a single-value AttributeBuilder<K, V> 
  enforcing one entry per locus tag; duplicate keys now halt the parser with a diagnostic 
  message instead of silently dropping data
- SourceAttributes, FeatureAttributes, and SequenceAttributes converted from enums 
  to structs; fields are now accessed directly instead of via pattern matching
- All optional fields on attribute structs are now typed as Option<T> reflecting their 
  true presence/absence semantics in GenBank/EMBL files
- We had two types of Records (gbk and embl) after parsing, they have been converted to one type of records after parsing - several structures and functions have been moved from gbk.rs and embl.rs into record.rs allowing a generic records type to be used for both
- Heatmap path fix
- Moved images folder to assets in docs windows install section

### Removed ###


