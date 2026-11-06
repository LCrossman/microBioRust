use std::io::{BufRead, Error, ErrorKind};

use phylo::prelude::*;

struct Scanner<R: BufRead> {
    col: usize,
    line: usize,
    reader: R,
    scratch: Vec<u8>,
}

impl<R: BufRead> Scanner<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            line: 1,
            col: 0,
            scratch: vec![],
        }
    }
}

impl<R: BufRead> Iterator for Scanner<R> {
    type Item = std::io::Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let buffer = self.reader.fill_buf().ok()?;
            let mut consumed_amt = 0;
            let mut flush = false;
            let buffer_len = buffer.len();

            if buffer.is_empty() {
                if self.scratch.is_empty() {
                    return None; // reached EOF without any issues
                } else {
                    return Some(Err(Error::new(
                        ErrorKind::InvalidData,
                        format!(
                            "Incomplete newick tree found at line: {}, col:{}",
                            self.line, self.col
                        ),
                    )));
                }
            }

            while consumed_amt < buffer_len {
                let b = buffer[consumed_amt];
                consumed_amt += 1;

                if b == b'\n' {
                    self.line += 1;
                    self.col = 0;
                } else {
                    self.col += 1;
                }

                if b.is_ascii_whitespace() {
                    continue;
                }
                self.scratch.push(b);
                if b == b';' {
                    flush = true;
                    break;
                }
            }
            self.reader.consume(consumed_amt);
            if flush {
                return Some(Ok(std::mem::take(&mut self.scratch)));
            }
        }
    }
}

struct PhyloParser<R: BufRead> {
    scanner: Scanner<R>,
}

impl<R: BufRead> PhyloParser<R> {
    fn new(scanner: Scanner<R>) -> Self {
        Self { scanner }
    }
}

impl<R: BufRead> Iterator for PhyloParser<R> {
    type Item = std::io::Result<PhyloTree>;

    fn next(&mut self) -> Option<Self::Item> {
        let bytes = match self.scanner.next()? {
            Ok(b) => b,
            Err(e) => return Some(Err(e)),
        };
        Some(
            PhyloTree::from_newick(&bytes)
                .map_err(|e| Error::new(ErrorKind::InvalidData, e.to_string())),
        )
    }
}

#[cfg(test)]
mod tests {
    use std::{fs::File, io::BufReader};

    use super::*;

    fn get_scanner(file_name: &str) -> Scanner<BufReader<File>> {
        let path = format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), file_name);
        let file_handle = File::open(path).unwrap();
        Scanner::new(BufReader::new(file_handle))
    }

    #[test]
    fn test_scanner_newick_multitree_file() {
        let mut scanner = get_scanner("basic_trees.nw");
        let first_tree = scanner.next().unwrap().unwrap();
        assert_eq!(String::from_utf8(first_tree).unwrap(), "(,,(,));")
    }

    #[test]
    fn test_basic_parsing() {
        let scanner = get_scanner("basic_trees.nw");
        let mut parser = PhyloParser::new(scanner);
        let first_tree = parser.next().unwrap().unwrap();
        assert_eq!(first_tree.to_newick().to_string(), "(,,(,));");
    }
}
