#[macro_export]
macro_rules! impl_build_features {
    ($fn_name:ident, $record_type:ty, $enum_path:path) => {
        fn $fn_name(
            r: &$record_type,
            py: ::pyo3::Python<'_>
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
                let obj = ::pyo3::Py::new(py, info).expect("failed to allocate PyFeatureInfo");
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
            py: ::pyo3::Python<'_>
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
                    let obj = ::pyo3::Py::new(py, info).expect("failed to allocate PySequenceInfo");
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
        #[pyclass]
        pub struct $struct_name {
            // Use absolute paths so it works in any module
            pub inner: ::std::collections::HashMap<String, ::pyo3::Py<$item_type>>,
        }

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
                    .ok_or_else(|| ::pyo3::exceptions::PyKeyError::new_err(tag.to_string()))
            }

            fn keys<'py>(&self, py: ::pyo3::Python<'py>) -> ::pyo3::PyResult<::pyo3::Bound<'py, ::pyo3::types::PyList>> {
                ::pyo3::types::PyList::new(py, self.inner.keys())
            }

            fn __iter__(slf: ::pyo3::PyRef<'_, Self>) -> $crate::LocusTagIterator {
                $crate::LocusTagIterator {
                    keys: slf.inner.keys().cloned().collect(),
                    index: 0,
                }
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
