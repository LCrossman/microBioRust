 use clap::Parser;
 use std::{
     fs::File,
     io::{Write, BufWriter},
 };
 use microBioRust::{
     genbank,
 };

 #[derive(Parser, Debug)]
 #[clap(author, version, about)]
 struct Arguments {
 #[clap(short, long)]
 filename: String,
 #[clap(short, long)]
 output: String,
 }

fn main() -> Result<(), anyhow::Error> {
            let args = Arguments::parse();
            let records = genbank!(&args.filename);
            let file = File::create(&args.output)?;
            let mut writer = BufWriter::new(file);
            for record in records {
               for (k, _v) in &record.cds.attributes {
                  if let Some(seq) = record.seq_features.get_sequence_faa(k) {
                     writeln!(writer, ">{}|{}\n{}", &record.id, &k, seq)?;
                     }
                  }
            }
            writer.flush()?;
            Ok(())
}

