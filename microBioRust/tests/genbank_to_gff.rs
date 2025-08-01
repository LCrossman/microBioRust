use microBioRust::gbk::{gff_write, Reader, Record};
use std::collections::BTreeMap;
use std::fs;
use std::io;
#[test]
pub fn genbank_to_gff() -> io::Result<()> {
    let file_gbk = fs::File::open("test_output.gbk")?;
    let _prev_start: u32 = 0;
    let mut prev_end: u32 = 0;
    let reader = Reader::new(file_gbk);
    let mut records = reader.records();
    let mut read_counter: u32 = 0;
    let mut seq_region: BTreeMap<String, (u32, u32)> = BTreeMap::new();
    let mut record_vec: Vec<Record> = Vec::new();
    loop {
        match records.next() {
            Some(Ok(record)) => {
                //println!("next record");
                //println!("Record id: {:?}", record.id);
                let sour = record
                    .source_map
                    .source_name
                    .clone()
                    .expect("issue collecting source name");
                let beginning = match record.source_map.get_start(&sour) {
                    Some(value) => value.get_value(),
                    _ => 0,
                };
                let ending = match record.source_map.get_stop(&sour) {
                    Some(value) => value.get_value(),
                    _ => 0,
                };
                if ending + prev_end < beginning + prev_end {
                    println!(
                        "debug since the end value smaller is than the start {:?}",
                        beginning
                    );
                }
                seq_region.insert(sour, (beginning + prev_end, ending + prev_end));
                record_vec.push(record);
                // Add additional fields to print if needed
                read_counter += 1;
                prev_end += ending;
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
    let output_file = format!("test_output.gff");
    gff_write(seq_region.clone(), record_vec, &output_file, true)?;
    println!("Total records processed: {}", read_counter);
    return Ok(());
}
