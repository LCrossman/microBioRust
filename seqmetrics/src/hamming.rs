#[allow(unused_imports)]

//function to calculate the hamming distance between 2 sequences
pub fn hamming_distance(seq1: &str, seq2: &str) -> usize {
    //sequences must be of equal length, i.e. a sequence alignment
    assert_eq!(
        seq1.len(),
        seq2.len(),
        "Sequences must be of equal length for Hamming distance"
    );
    seq1.as_bytes()
        .iter()
        .zip(seq2.as_bytes())
        .filter(|(a, b)| a != b)
        .count()
}

pub async fn hamming_matrix(sequences: &Vec<String>) -> Result<Vec<Vec<usize>>, anyhow::Error> {
    assert!(
        sequences.windows(2).all(|w| w[0].len() == w[1].len()),
        "all sequences must be the same length, i.e an alignment"
    );
    println!("all sequences are the same length, proceeding");
    let n = sequences.len();
    let mut distances = vec![vec![0usize; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let d = hamming_distance(&sequences[i], &sequences[j]);
            distances[i][j] = d;
            distances[j][i] = d;
        }
    }
    Ok(distances)
}
