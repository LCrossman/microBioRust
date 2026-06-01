#[macro_export]
macro_rules! impl_build_features {
    ($fn_name:ident, $record_type:ty, $enum_path:path) => {
        fn $fn_name(
            r: &$record_type,
            py: ::pyo3::Python<'_>,
        ) -> ::std::collections::HashMap<String, ::pyo3::Py<$crate::PyFeatureInfo>> {
            use $enum_path as Attrs;
            let mut map = ::std::collections::HashMap::new();
            for (tag, attrs) in &r.cds.attributes {
                let tag_owned = tag.clone();
                let mut info = $crate::PyFeatureInfo::new(&tag_owned);
                for attr in attrs {
                    #[allow(unreachable_patterns)]
                    match attr {
                        Attrs::Gene { value } => info.gene = Some(value.clone()),
                        Attrs::Product { value } => info.product = Some(value.clone()),
                        Attrs::Start { value } => info.start = Some(value.get_value()),
                        Attrs::Stop { value } => info.stop = Some(value.get_value()),
                        Attrs::Strand { value } => info.strand = Some(*value),
                        Attrs::CodonStart { value } => info.codon_start = Some(*value),
                        other => {
                            info.extras.push(format!("{:?}", other));
                        }
                    }
                }
                let obj = ::pyo3::Py::new(py, info).expect("failed to allocate the PyFeatureInfo");
                map.insert(tag_owned, obj);
            }
            map
        }
    };
}
#[macro_export]
macro_rules! impl_build_sequences {
    ($fn_name:ident, $record_type:ty, $enum_path:path) => {
        fn $fn_name(
            r: &$record_type,
            py: ::pyo3::Python<'_>,
        ) -> ::std::collections::HashMap<String, ::pyo3::Py<$crate::PySequenceInfo>> {
            use $enum_path as Attrs;
            let mut map = ::std::collections::HashMap::new();
            for (tag, attrs) in &r.seq_features.seq_attributes {
                let tag_owned = tag.clone();
                let mut info = $crate::PySequenceInfo::new(&tag_owned);
                for attr in attrs {
                    #[allow(unreachable_patterns)]
                    match attr {
                        Attrs::SequenceFaa { value } => info.faa = Some(value.clone()),
                        Attrs::SequenceFfn { value } => info.ffn = Some(value.clone()),
                        other => {
                            info.extras.push(format!("{:?}", other));
                        }
                    }
                }
                if info.faa.is_some() || info.ffn.is_some() {
                    let obj =
                        ::pyo3::Py::new(py, info).expect("failed to allocate the PySequenceInfo");
                    map.insert(tag_owned, obj);
                }
            }
            map
        }
    };
}
#[macro_export]
macro_rules! impl_build_faa {
    ($fn_name:ident, $record_type:ty, $enum_path:path) => {
        fn $fn_name(
            r: &$record_type,
            py: ::pyo3::Python<'_>,
        ) -> ::std::collections::HashMap<String, ::pyo3::Py<$crate::PyFaaInfo>> {
            use $enum_path as Attrs;
            let mut map = ::std::collections::HashMap::new();
            for (tag, attrs) in &r.seq_features.seq_attributes {
                let tag_owned = tag.clone();
                let mut faa: Option<String> = None;
                for attr in attrs {
                    #[allow(unreachable_patterns)]
                    match attr {
                        Attrs::SequenceFaa { value } => faa = Some(value.clone()),
                        _ => {}
                    }
                }
                if let Some(faa_seq) = faa {
                    let info = $crate::PyFaaInfo {
                        locus_tag: tag_owned.clone(),
                        faa: faa_seq,
                    };
                    let obj = ::pyo3::Py::new(py, info).expect("failed to allocate a PyFaaInfo");
                    map.insert(tag_owned, obj);
                }
            }
            map
        }
    };
}

#[macro_export]
macro_rules! impl_build_ffn {
    ($fn_name:ident, $record_type:ty, $enum_path:path) => {
        fn $fn_name(
            r: &$record_type,
            py: ::pyo3::Python<'_>,
        ) -> ::std::collections::HashMap<String, ::pyo3::Py<$crate::PyFfnInfo>> {
            use $enum_path as Attrs;
            let mut map = ::std::collections::HashMap::new();
            for (tag, attrs) in &r.seq_features.seq_attributes {
                let tag_owned = tag.clone();
                let mut ffn: Option<String> = None;
                for attr in attrs {
                    #[allow(unreachable_patterns)]
                    match attr {
                        Attrs::SequenceFfn { value } => ffn = Some(value.clone()),
                        _ => {}
                    }
                }
                if let Some(ffn_seq) = ffn {
                    let info = $crate::PyFfnInfo {
                        locus_tag: tag_owned.clone(),
                        ffn: ffn_seq,
                    };
                    let obj = ::pyo3::Py::new(py, info).expect("failed to allocate a PyFfnInfo");
                    map.insert(tag_owned, obj);
                }
            }
            map
        }
    };
}
#[macro_export]
macro_rules! create_collection {
    ($struct_name:ident, $item_type:ty, $label:expr) => {
        #[::pyo3_stub_gen::derive::gen_stub_pyclass]
        #[pyclass]
        pub struct $struct_name {
            // Use absolute paths so it works in any module
            pub inner: ::std::collections::HashMap<String, ::pyo3::Py<$item_type>>,
        }
        #[::pyo3_stub_gen::derive::gen_stub_pymethods]
        #[pymethods]
        impl $struct_name {
            fn __len__(&self) -> usize {
                self.inner.len()
            }

            fn __contains__(&self, tag: &str) -> bool {
                self.inner.contains_key(tag)
            }

            fn __repr__(&self) -> String {
                format!("{}({} entries)", $label, self.inner.len())
            }

            fn __getitem__<'py>(
                &self,
                tag: &str,
                py: ::pyo3::Python<'py>,
            ) -> ::pyo3::PyResult<::pyo3::Bound<'py, $item_type>> {
                self.inner
                    .get(tag)
                    .map(|obj| obj.bind(py).clone())
                    .ok_or_else(|| {
                        let available = self
                            .inner
                            .keys()
                            .take(5)
                            .cloned()
                            .collect::<Vec<_>>()
                            .join(", ");
                        ::pyo3::exceptions::PyKeyError::new_err(format!(
                            "locus tag {} is not found. The first 5 available tags are: {},\n
                            Use list(collection.keys()) to see all tags.",
                            tag, available
                        ))
                    })
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
            fn keys<'py>(
                &self,
                py: ::pyo3::Python<'py>,
            ) -> ::pyo3::PyResult<::pyo3::Bound<'py, ::pyo3::types::PyList>> {
                ::pyo3::types::PyList::new(py, self.inner.keys())
            }

            fn __iter__(
                slf: ::pyo3::PyRef<'_, Self>,
            ) -> ::pyo3::PyResult<::pyo3::Py<$crate::LocusTagIterator>> {
                let iter = $crate::LocusTagIterator {
                    keys: slf.inner.keys().cloned().collect(),
                    index: 0,
                };
                ::pyo3::Py::new(slf.py(), iter).map_err(|e| e.into())
            }
            /// Serializes the collection to a formatted JSON string
            fn to_json<'py>(&self, py: ::pyo3::Python<'py>) -> ::pyo3::PyResult<String> {
                // 1. Create a pure Rust HashMap
                let mut export_map = ::std::collections::HashMap::new();

                // 2. Extract the pure Rust data from the Py wrappers
                for (key, py_ptr) in &self.inner {
                    let rust_struct = py_ptr.bind(py).borrow().clone();
                    export_map.insert(key.clone(), rust_struct);
                }
                // 3. Serialize using serde_json
                ::serde_json::to_string_pretty(&export_map).map_err(|e| {
                    ::pyo3::exceptions::PyValueError::new_err(format!(
                        "JSON serialization failed: {}",
                        e
                    ))
                })
            }
        }
    };
}
//register the functions on submodules and parent
#[macro_export]
macro_rules! register_all {
    ($sub:expr, $parent:expr, $($func:ident),+ $(,)?) => {
        $(
            $sub.add_function(::pyo3::wrap_pyfunction!($func, &$sub)?)?;
            $parent.add_function(::pyo3::wrap_pyfunction!($func, $parent)?)?;
        )+
    };
}
//only register the functions on the submodules
//#[macro_export]
//macro_rules! register_functions {
//    ($module:expr, $($func:ident),*) => {
//        $(
// Using the full path pyo3::wrap_pyfunction ensures it always resolves
//            $module.add_function(pyo3::wrap_pyfunction!($func, &$module)?)?;
//        )*
//    };
//}
