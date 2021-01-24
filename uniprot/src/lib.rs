//! Retrieve basic data about a protein from a local UniprotKB file
//!
//! # File format
//!
//! Files should be tab delimited, with 4 fields: Uniprot Accession, Gene Name,
//! Molecular Weight, and Protein Sequence. Each protein appears on it's own
//! line in the file, and no header
//! should be present
//!
//! ```text
//! $ cat uniprot.txt
//! ...
//! Q13526	PIN1	18243	MADEEKLPPGWEKRMSRSSGRVYYF ... GEMSGPVFTDSGIHIILRTE
//! Q13547	HDAC1	55103	MAQTQGTRRKVCYYYDGDVGNYYYG ... EKTKEEKPEAKGVKEEVKLA
//! ...
//! ```
//!
//! # Example
//!
//! ```rust,ignore
//! # use uniprot::Uniprot;
//! let db: Uniprot = match Uniprot::load("uniprot.txt") {
//!     Ok(db) => db,
//!     Err(e) => panic!("Error loading local Uniprot database: {}", e),
//! };
//! ```
//!

pub mod kw;

use memchr::{memchr_iter, Memchr};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, prelude::*};
use std::path::Path;
use std::str;

pub mod fasta;

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

#[derive(Debug, PartialEq, Clone)]
/// Represents an entry in the locally-stored UniprotKB database
pub struct Entry {
    /// Uniprot accession identifier
    pub accession: String,
    /// Gene name
    pub identifier: String,
    /// Protein sequence
    pub sequence: String,
}

#[derive(Debug, PartialEq, Clone)]
/// Wraps a hashtable
pub struct Uniprot {
    inner: HashMap<String, Entry>,
}

impl Uniprot {
    /// Load a UniprotKB file to build a [`Uniprot`] object
    ///
    /// # File format
    ///
    /// Files should be tab delimited, with 4 fields: Uniprot Accession, Gene
    /// Name, Molecular Weight, and Protein Sequence. Each protein appears
    /// on it's own line in the file, and no header
    /// should be present
    ///
    /// ```text
    /// $ cat uniprot.txt
    /// ...
    /// Q13526	PIN1	18243	MADEEKLPPGWEKRMSRSSGRVYYF ... GEMSGPVFTDSGIHIILRTE
    /// Q13547	HDAC1	55103	MAQTQGTRRKVCYYYDGDVGNYYYG ... EKTKEEKPEAKGVKEEVKLA
    /// ...
    /// ```
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// # use uniprot::Uniprot;
    ///
    /// let db: Uniprot = match Uniprot::load("uniprot.txt") {
    ///     Ok(db) => db,
    ///     Err(e) => panic!("Error loading local Uniprot database: {}", e),
    /// };
    /// ```
    // pub fn load<P: AsRef<Path>>(path: P) -> io::Result<Uniprot> {
    //     let mut buf = String::new();
    //     File::open(path)?.read_to_string(&mut buf)?;

    //     let mut inner = HashMap::new();

    //     let mut iter = Pitchfork::new('\n' as u8, buf.as_bytes()).filter(|x| x.len() != 0);
    //     let mut last = iter.next().unwrap();

    //     let mut id = "";
    //     let mut seq = String::new();
    //     for line in iter {
    //         if line[0] == '>' as u8 {
    //             let id = str::from_utf8(last).unwrap();
    //             last = line;

    //             if id != "" {
    //                 let deets = id.split('|').skip(1).take(2).collect::<Vec<&str>>();
    //                 if deets.len() != 2 {
    //                     panic!("{}", id);
    //                 }
    //                 inner.insert(deets[0].into(), Entry {
    //                     accession: deets[0].into(),
    //                     identifier: deets[1].into(),
    //                     sequence: seq,
    //                 });
    //                 seq = String::new();
    //             }
    //         } else {
    //             seq.push_str(str::from_utf8(line).unwrap().trim());
    //         }
    //     }

    //     Ok(Uniprot {
    //         inner
    //     })

    // }

    #[deprecated(since = "0.2.0", note = "Use fasta::Fasta interface instead")]
    pub fn load<T: AsRef<Path>>(path: T) -> io::Result<Self> {
        let mut buffer = String::new();
        let mut inner = HashMap::new();

        File::open(path)?.read_to_string(&mut buffer)?;

        for line in buffer.lines() {
            let s = line.split('\t').collect::<Vec<_>>();

            let entry = Entry {
                accession: s[0].into(),
                identifier: s[1].into(),
                sequence: s[3].into(),
                // molecular_weight: s[2]
                //     .parse::<u32>()
                //     .map_err(|_| io::Error::from(io::ErrorKind::InvalidData))?,
            };

            inner.insert(entry.accession.clone(), entry);
        }
        Ok(Uniprot { inner })
    }

    /// Search a [`Uniprot`] database object by Uniprot accession
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// # use uniprot::Uniprot;
    /// # let uni = Uniprot::load("./data/uniprot/db").unwrap();
    ///
    /// let seq: &String = match uni.lookup("Q13526") {
    ///     Some(entry) => &entry.sequence,
    ///     None => 0,
    /// };
    /// ```
    pub fn lookup<T: AsRef<str>>(&self, acc: T) -> Option<&Entry> {
        self.inner.get(acc.as_ref())
    }
}

impl Entry {
    pub fn assign_residue(&self, seq: &str) -> Option<usize> {
        let peptide = if seq.contains(".") {
            seq.split(".").skip(1).next()?
        } else {
            seq
        };
        match peptide.find('*') {
            Some(offset) => {
                let needle = peptide.chars().filter(|&c| c != '*').collect::<String>();
                self.sequence.find(&needle).map(|idx| idx + offset)
            }
            None => self.sequence.find(peptide),
        }
    }
}
