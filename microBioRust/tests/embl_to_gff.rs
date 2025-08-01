use microBioRust::embl::{gff_write, Reader, Record};
use std::collections::BTreeMap;
use std::fs;

#[test]
fn test_embl_to_gff() -> std::io::Result<()> {
    let file_embl = fs::File::open("example.embl")?;
    let reader = Reader::new(file_embl);
    let mut records = reader.records();
    let mut read_counter: u32 = 0;
    let mut prev_end: u32 = 0;
    let mut seq_region: BTreeMap<String, (u32, u32)> = BTreeMap::new();
    let mut record_vec: Vec<Record> = Vec::new();

    while let Some(record_result) = records.next() {
        match record_result {
            Ok(record) => {
                let sour = record
                    .source_map
                    .source_name
                    .clone()
                    .expect("Missing source name");
                let beginning = record
                    .source_map
                    .get_start(&sour)
                    .map(|v| v.get_value())
                    .unwrap_or(0);
                let ending = record
                    .source_map
                    .get_stop(&sour)
                    .map(|v| v.get_value())
                    .unwrap_or(0);

                if ending + prev_end < beginning + prev_end {
                    println!("start > end: {:?}", beginning);
                }

                seq_region.insert(sour, (beginning + prev_end, ending + prev_end));
                record_vec.push(record);
                read_counter += 1;
                prev_end += ending;
            }
            Err(e) => eprintln!("Error: {:?}", e),
        }
    }

    let output_file = "test_output_embl.gff";
    gff_write(seq_region.clone(), record_vec, output_file, true)?;
    println!("Total records processed: {}", read_counter);
    Ok(())
}
