# Changelog

## [0.1.7] - 2026-05-30

### Added
- `PyFaaInfo` struct — dedicated return type for FAA sequences with non-Optional `faa: String` field
- `PyFfnInfo` struct — dedicated return type for FFN sequences with non-Optional `ffn: String` field
- `FaaCollection` — typed collection of `PyFaaInfo`, instantiated via `create_collection!` macro
- `FfnCollection` — typed collection of `PyFfnInfo`, instantiated via `create_collection!` macro
- `record.faa()` method on `PyRecord` — returns `FaaCollection` for per-genome protein access
- `record.ffn()` method on `PyRecord` — returns `FfnCollection` for per-genome nucleotide access
- `impl_build_faa!` macro — builds `FaaCollection` from gbk/embl record types
- `impl_build_ffn!` macro — builds `FfnCollection` from gbk/embl record types
- `build_faa_from_gbk`, `build_faa_from_embl` — instantiated from `impl_build_faa!`
- `build_ffn_from_gbk`, `build_ffn_from_embl` — instantiated from `impl_build_ffn!`
- `values()`, `items()`, `keys()` methods on all collections — enables Pythonic dict-like access
- `__len__` and `__contains__` on all collections
- `LocusTagIterator` — shared iterator type for all collections, fixes `iter() returned non-iterator` error
- `embl_to_ffn` function — was missing from embl module
- `to_json()` method function - to allow creation of a json string on all Collections

### Changed
- gbk_to_faa now returns FaaCollection of PyFaaInfo instead of SequenceCollection of generic PySequenceInfo
- gbk_to_ffn now returns FfnCollection of PyFfnInfo instead of SequenceCollection of generic PySequenceInfo
- embl_to_faa now returns FaaCollection of PyFaaInfo instead of SequenceCollection of generic PySequenceInfo
- embl_to_ffn now returns FfnCollection of PyFfnInfo instead of SequenceCollection of generic PySequenceInfo
- __getitem__ PyError message now includes the requested tag name and provides a hint for list(collection.keys())
- LocusTagIterator.__next__ changed from &mut self to PyRefMut<'_, Self> to fix iterator protocol registration
- create_collection! macro __iter__ now returns PyResult<Py<LocusTagIterator>> instead of bare struct
- RecordCollection.__iter__ fixed to return PyResult<Py<LocusTagIterator>>
- Deprecated prepare_freethreaded_python replaced with Python::initialize
- Deprecated Python::with_gil replaced with Python::attach
- Deprecated PyObject replaced with Py<PyAny>
- Added serde flags to existing structs of PyFaaInfo, PyFfnInfo, PySequenceInfo, PyFeatureInfo and Collections
- Cargo.toml default features changed, extension-module removed from defaults to allow cargo test to link libpython correctly

### Fixed
- `iter() returned non-iterator of type 'builtins.LocusTagIterator'` — fixed by returning `Py<LocusTagIterator>` from `__iter__`
- `dyld: symbol not found '_PyBaseObject_Type'` on macOS in cargo test only, was resolved by removing `extension-module` from default features

### Design Decisions
- `PySequenceInfo` reserved exclusively for `parse_gbk`/`parse_embl` — retains `Option<String>` fields for faa/ffn since not all locus tags have both
- `PyFaaInfo`/`PyFfnInfo` use non-Optional fields — if a record is in the collection it by definition has the sequence

### New Tests Added
