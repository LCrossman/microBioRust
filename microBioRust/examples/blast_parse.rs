use anyhow::{Context, Result};
use async_compression::tokio::bufread::GzipDecoder as AsyncGzDecoder;
use clap::Parser;
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use quick_xml::escape::unescape;
use serde::Serialize;
use serde_json::ser::Serializer as JsonSerializer;
use microBioRust::blast::*;
use std::io::Cursor;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncWriteExt, BufReader};

#[derive(Parser, Debug)]
#[command(name = "blast-parsers", author, version, about = "async microBioRust BLAST parsers: for outfmt6 (single line tabular) and outfmt5 (xml)")]
struct Cli {
    ///Use .gz for gzip-compressed files.
    #[arg(short, long, default_value = "-")]
    input: String,
    /// Format: '6' (tabular) or '5' (xml). If omitted we try to infer by file suffix only
    #[arg(short, long)]
    format: Option<String>,
    /// Output newline-delimited JSON (one JSON object per record/iteration)
    #[arg(long)]
    json: bool,
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = Cli::parse();
    let fmt = infer_format(&args.input, &args.format);
    let reader_box = open_async_reader(&args.input).await?;
    if fmt == "6" {
        stream_outfmt6_to_json(reader_box).await?;
    } else {
        // Build AsyncBlastXmlIter from reader_box
        let iter_reader = reader_box;
        let mut iter = AsyncBlastXmlIter::from_reader(iter_reader);
        while let Some(res) = iter.next_iteration().await {
            match res {
                Ok(iter_rec) => {
                    if args.json {
                        let mut buf = Vec::new();
                        serde_json::to_writer(&mut buf, &iter_rec)?;
                        buf.push(b'\n');
                        tokio::io::stdout().write_all(&buf).await?;
                    } else {
                        println!("query {:?} hits {}", iter_rec.query_def, iter_rec.hits.len());
                    }
                }
                Err(e) => eprintln!("xml parse error: {}", e),
            }
        }
    }

    Ok(())
}
