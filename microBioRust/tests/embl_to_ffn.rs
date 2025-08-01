use microBioRust::embl::Reader;
use std::fs;
#[test]
pub fn embl_to_ffn() -> Result<(), anyhow::Error> {
    let file_embl = fs::File::open("example.embl")?;
    let reader = Reader::new(file_embl);
    let mut records = reader.records();
    let mut read_counter: u32 = 0;
    loop {
        match records.next() {
            Some(Ok(record)) => {
                //println!("next record");
                //println!("Record id: {:?}", record.id);
                for (k, _v) in &record.cds.attributes {
                    match record.seq_features.get_sequence_ffn(&k) {
                        Some(value) => {
                            let seq_ffn = value.to_string();
                            println!(">{}|{}\n{}", &record.id, &k, seq_ffn);
                        }
                        _ => (),
                    };
                }
                read_counter += 1;
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
    println!("Total records processed: {}", read_counter);
    return Ok(());
}
