///Shared generic record types to reduce duplication between gbk and embl
///Minimal initial introduction: defines generic containers and builders that mirror the existing API where possible

use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum RangeValue {
    Exact(u32),
    LessThan(u32),
    GreaterThan(u32),
}

impl RangeValue {
    pub fn get_value(&self) -> u32 {
        match self {
            RangeValue::Exact(v) => *v,
            RangeValue::LessThan(v) => *v,
            RangeValue::GreaterThan(v) => *v,
        }
    }
}

///Traits to unify attribute enums across formats. Existing enums can implement Into these trait views if needed
pub trait HasStartStopStrand {
    fn start(&self) -> Option<RangeValue> { None }
    fn stop(&self) -> Option<RangeValue> { None }
    fn strand(&self) -> Option<i8> { None }
}

///Generic attribute builders
#[derive(Clone, Debug, Default)]
pub struct AttributeBuilder<K, V> {
    pub name: Option<String>,
    pub attributes: HashMap<K, HashSet<V>>,
}

impl<K, V> AttributeBuilder<K, V>
where
    K: Eq + std::hash::Hash,
    V: Eq + std::hash::Hash,
{
    pub fn set_name(&mut self, name: String) { self.name = Some(name); }
    pub fn get_name(&self) -> Option<&String> { self.name.as_ref() }
    pub fn add(&mut self, key: K, value: V) {
        self.attributes.entry(key).or_insert_with(HashSet::new).insert(value);
    }
    pub fn get(&self, key: &K) -> Option<&HashSet<V>> { self.attributes.get(key) }
}

///Generic record and records container
#[derive(Clone, Debug, Default)]
pub struct GenericRecord<S, F, Q> {
    pub id: String,
    pub seq: String,
    pub seqid: String,
    pub start: u32,
    pub end: u32,
    pub strand: i32,
    pub source: S,
    pub cds: F,
    pub seq_features: Q,
}

impl<S, F, Q> GenericRecord<S, F, Q> {
    pub fn is_empty(&self) -> bool { self.id.is_empty() && self.seq.is_empty() }
}

pub struct GenericRecords<R> {
    inner: R,
}

impl<R> GenericRecords<R> {
    pub fn new(reader: R) -> Self { Self { inner: reader } }
}
