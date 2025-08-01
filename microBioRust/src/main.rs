use clap::Parser;
use microBioRust::genbank;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Arguments {
    #[clap(short, long)]
    filename: String,
}

fn main() -> Result<(), anyhow::Error> {
    let args = Arguments::parse();
    let records = genbank!(&args.filename);
    for record in records.iter() {
        for (k, _v) in &record.cds.attributes {
            if let Some(seq) = record.seq_features.get_sequence_faa(k) {
                println!(">{}|{}\n{}", &record.id, &k, seq);
            }
        }
    }
    return Ok(());
}
