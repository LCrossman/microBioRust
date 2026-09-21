use paste::paste;
///Shared generic record types to reduce duplication between gbk and embl
///Minimal initial introduction: defines generic containers and builders that mirror the existing API where possible
use serde::Serialize;
use std::collections::{BTreeMap, HashSet};
use std::fmt;
use std::io;

#[macro_export]
macro_rules! create_getters {
    ($struct_name:ident,$attributes:ident, $enum_name:ident,$( $field:ident { value:$type:ty } ),* ) => {
        impl $struct_name {$(
                paste! {
                  pub fn [<get_$field:snake>](&self, key: &str) -> Option<&$type> {
                      // Access the inner `attributes` BTreeMap on the AttributeBuilder
                      self.$attributes.attributes.get(key).and_then(|set| {
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

#[macro_export]
macro_rules! create_builder {
    ($builder_name:ident, $attributes:ident,$enum_name:ident, $counter_name:ident,$( $field:ident { value:$type:ty } ),* ) => {
        impl $builder_name {
            pub fn new() -> Self {
                $builder_name {
                    $attributes: AttributeBuilder::new(),$counter_name: None,
                }
            }
            pub fn set_counter(&mut self, counter: String) -> &mut Self {
                self.$counter_name = Some(counter);
                self
            }
            pub fn insert_to(&mut self, value: $enum_name) {
                if let Some(counter) = &self.$counter_name {
                    // Access the inner `attributes` BTreeMap
                    self.$attributes.attributes
                        .entry(counter.to_string())
                        .or_insert_with(HashSet::new)
                        .insert(value);
                } else {
                    panic!("Counter key not set");
                }
            }
            $(
              paste! {
                pub fn [<set_$field:snake>](&mut self, value:$type) -> &mut Self {
                   self.insert_to($enum_name::$field { value });
                   self
                }
              }
            )*
            pub fn build(self) -> BTreeMap<String, HashSet<$enum_name>> {
                // Call build on the underlying AttributeBuilder
                self.$attributes.build()
            }
            pub fn iter_sorted(&'_ self) -> std::collections::btree_map::Iter<'_, String, HashSet<$enum_name>> {
                // Call iter_sorted on the underlying AttributeBuilder
                self.$attributes.iter_sorted()
            }
            pub fn default() -> Self {
                Self::new()
            }
        }
    };
}

//stores the details of the source features in genbank (contigs)
#[derive(Debug, Serialize, Eq, PartialEq, Hash, Clone)]
pub enum SourceAttributes {
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Organism { value: String },
    MolType { value: String },
    Strain { value: String },
    CultureCollection { value: String },
    TypeMaterial { value: String },
    DbXref { value: String },
}

#[derive(Clone, Serialize, Debug, Default)]
pub struct SourceAttributeBuilder {
    pub source_attributes: AttributeBuilder<String, SourceAttributes>,
    pub source_name: Option<String>, // source-level only, not in the generic
}

#[derive(Clone, Serialize, Debug, Default)]
pub struct FeatureAttributeBuilder {
    pub attributes: AttributeBuilder<String, FeatureAttributes>,
    pub locus_tag: Option<String>,
}

#[derive(Clone, Serialize, Debug, Default)]
pub struct SequenceAttributeBuilder {
    pub seq_attributes: AttributeBuilder<String, SequenceAttributes>,
    pub locus_tag: Option<String>,
}
//macro for creating the getters
create_getters!(
    SourceAttributeBuilder,
    source_attributes,
    SourceAttributes,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Organism { value: String },
    MolType { value: String },
    Strain { value: String },
    // CultureCollection { value: String},
    TypeMaterial { value: String },
    DbXref { value: String }
);

create_builder!(
    SourceAttributeBuilder,
    source_attributes,
    SourceAttributes,
    source_name,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    Organism { value: String },
    MolType { value: String },
    Strain { value: String },
    // CultureCollection { value: String},
    TypeMaterial { value: String },
    DbXref { value: String }
);

///attributes for each feature, cds or gene
#[derive(Debug, Eq, Serialize, Hash, PartialEq, Clone)]
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
#[derive(Debug, Eq, Serialize, PartialEq, Hash, Clone)]
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
    SequenceFfn { value: String },
    SequenceFaa { value: String },
    CodonStart { value: u8 },
    Strand { value: i8 }
);

create_builder!(
    SequenceAttributeBuilder,
    seq_attributes,
    SequenceAttributes,
    locus_tag,
    Start { value: RangeValue },
    Stop { value: RangeValue },
    SequenceFfn { value: String },
    SequenceFaa { value: String },
    CodonStart { value: u8 },
    Strand { value: i8 }
);

#[derive(Clone, Serialize, Debug, PartialEq, Eq, Hash)]
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
    fn start(&self) -> Option<RangeValue> {
        None
    }
    fn stop(&self) -> Option<RangeValue> {
        None
    }
    fn strand(&self) -> Option<i8> {
        None
    }
}

///Generic attribute builders
#[derive(Clone, Serialize, Debug)]
pub struct AttributeBuilder<K, V> {
    pub name: Option<String>,
    pub counter: Option<K>,
    pub attributes: BTreeMap<K, HashSet<V>>,
}
impl<K, V> Default for AttributeBuilder<K, V> {
    fn default() -> Self {
        Self {
            name: None,
            counter: None,
            attributes: BTreeMap::new(),
        }
    }
}
impl<K, V> AttributeBuilder<K, V>
where
    K: Eq + std::hash::Hash + Ord + Clone,
    V: Eq + std::hash::Hash,
{
    pub fn new() -> Self {
        Self {
            name: None,
            counter: None,
            attributes: BTreeMap::new(),
        }
    }
    pub fn set_counter(&mut self, key: K) -> &mut Self {
        self.counter = Some(key);
        self
    }
    pub fn set_name(&mut self, name: String) {
        self.name = Some(name);
    }
    pub fn get_name(&self) -> Option<&String> {
        self.name.as_ref()
    }
    pub fn add(&mut self, key: K, value: V) {
        self.attributes
            .entry(key)
            .or_insert_with(HashSet::new)
            .insert(value);
    }
    pub fn get_set(&self, key: &K) -> Option<&HashSet<V>> {
        self.attributes.get(key)
    }
    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.attributes.keys()
    }
    pub fn first_key_value(&self) -> Option<(&K, &HashSet<V>)> {
        self.attributes.iter().next()
    }
    pub fn insert(&mut self, value: V) {
        let key = self.counter.clone().expect("Counter key not set");
        self.attributes.entry(key).or_default().insert(value);
    }

    pub fn iter_sorted(&self) -> std::collections::btree_map::Iter<'_, K, HashSet<V>> {
        self.attributes.iter()
    }

    pub fn build(self) -> BTreeMap<K, HashSet<V>> {
        self.attributes
    }
}
impl<'a, K, V> IntoIterator for &'a AttributeBuilder<K, V>
where
    K: Ord + Eq + std::hash::Hash + Clone,
    V: Eq + std::hash::Hash,
{
    type Item = (&'a K, &'a std::collections::HashSet<V>);
    type IntoIter = std::collections::btree_map::Iter<'a, K, std::collections::HashSet<V>>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_sorted()
    }
}

///Generic record and records container
#[derive(Clone, Serialize, Debug)]
pub struct Record {
    pub id: String,
    pub sequence: String,
    pub length: u32,
    pub start: usize,
    pub end: usize,
    pub strand: i32,
    pub source_map: SourceAttributeBuilder,
    pub cds: FeatureAttributeBuilder,
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
    pub(crate) fn rec_clear(&mut self) {
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

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Record")
    }
}

/// Trait implemented by all format-specific readers.
/// Each format (GBK, EMBL, GFF) implements this to parse one Record at a time.
pub trait RecordRead {
    fn read(&mut self, record: &mut Record) -> Result<Record, anyhow::Error>;
}

/// Generic iterator over Records — works for any format that implements RecordRead.
/// Lives in record.rs because it is format-agnostic.
#[derive(Debug)]
pub struct Records<R: RecordRead> {
    reader: R,
    error_has_occurred: bool,
}

impl<R: RecordRead> Records<R> {
    pub fn new(reader: R) -> Self {
        Records {
            reader,
            error_has_occurred: false,
        }
    }
}

impl<R: RecordRead> Iterator for Records<R> {
    type Item = Result<Record, anyhow::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.error_has_occurred {
            eprintln!("error was encountered in iteration");
            return None;
        }
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Ok(_) => {
                if record.is_empty() {
                    None
                } else {
                    Some(Ok(record))
                }
            }
            Err(e) => {
                self.error_has_occurred = true;
                Some(Err(anyhow::anyhow!("record read error: {:?}", e)))
            }
        }
    }
}
