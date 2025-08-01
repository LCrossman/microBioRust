use anyhow::Result;
use tokio::fs::File;
use tokio::io::{AsyncWriteExt, BufWriter};

pub async fn write_distances_csv(
    ids: Vec<String>,
    distances: Vec<Vec<usize>>,
    filename: &str,
) -> Result<(), anyhow::Error> {
    let file = File::create(filename).await?;
    let mut writer = BufWriter::new(file);
    let mut header = String::from("id");
    for id in &ids {
        header.push(',');
        header.push_str(id);
    }
    header.push('\n');
    writer.write_all(header.as_bytes()).await?;
    // Write each row
    for (i, row) in distances.iter().enumerate() {
        let mut line = format!("{}", ids[i]);
        for val in row {
            line.push_str(&format!(",{:.3}", val));
        }
        line.push('\n');
        writer.write_all(line.as_bytes()).await?;
    }
    writer.flush().await?;
    Ok(())
}
