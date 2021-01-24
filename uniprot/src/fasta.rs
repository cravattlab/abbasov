//! Utilities for loading genomic information
use memchr::{memchr_iter, Memchr};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, prelude::*};
use std::path::Path;
use std::str;

struct Pitchfork<'a> {
    pos: usize,
    haystack: &'a [u8],
    inner: Memchr<'a>,
}

impl<'a> Pitchfork<'a> {
    pub fn new(needle: u8, haystack: &'a [u8]) -> Self {
        Self {
            pos: 0,
            haystack,
            inner: memchr_iter(needle, haystack),
        }
    }
}

impl<'a> Iterator for Pitchfork<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let end = match self.inner.next() {
            Some(e) => e,
            None => {
                if self.pos < self.haystack.len() {
                    self.haystack.len()
                } else {
                    return None;
                }
            }
        };
        let slice = &self.haystack[self.pos..end];
        self.pos = end + 1;
        Some(slice)
    }
}

#[derive(Debug, Clone)]
pub struct Fasta {
    pub map: HashMap<String, String>,
}

impl Fasta {
    /// Build a fasta database
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Fasta> {
        let mut buf = String::new();
        File::open(path)?.read_to_string(&mut buf)?;

        let mut map = HashMap::new();
        let mut iter = Pitchfork::new('\n' as u8, buf.as_bytes());
        let mut last_id = iter.next().unwrap();
        let mut s = String::new();

        for line in iter {
            if line.len() == 0 {
                continue;
            }
            // dbg!(str::from_utf8(line).unwrap());
            if line[0] == '>' as u8 {
                if s != "" {
                    let id = str::from_utf8(last_id).unwrap();
                    if id.contains("Reverse") {
                        s.clear();
                        continue;
                    }
                    let acc = id.split('|').skip(1).next().unwrap().into();
                    map.insert(acc, std::mem::replace(&mut s, String::new()));
                    last_id = line;
                // s.clear();
                } else {
                    last_id = line;
                }
            } else {
                s.push_str(str::from_utf8(line).unwrap());
            }
        }
        Ok(Fasta { map })
    }

    pub fn sequence(&self, acc: &str) -> Option<&String> {
        self.map.get(acc)
    }

    pub fn assign(&self, acc: &str, seq: &str) -> Option<usize> {
        // Handle tryptic cleavage sites
        let peptide = if seq.contains(".") {
            seq.split(".").skip(1).next()?
        } else {
            seq
        };
        let primary = self.map.get(acc)?;
        match peptide.find('*') {
            Some(offset) => {
                let needle = peptide.chars().filter(|&c| c != '*').collect::<String>();
                primary.find(&needle).map(|idx| idx + offset)
            }
            None => primary.find(peptide),
        }
    }
}
