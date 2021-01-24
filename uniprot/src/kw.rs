//! Utilities for parsing Keywords from uniprot_sprot.dat database
//! The Keywords represent molecular function ontologies
//!

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

#[derive(Clone, Default, Debug, PartialEq)]
pub struct Annotation {
    pub kw: HashMap<String, String>,
    pub go: HashMap<String, String>,
    pub enz: HashSet<String>,
}

impl Annotation {
    pub fn keyword(&self, accession: &str) -> Option<&str> {
        self.kw.get(accession).map(|s| s as &str)
    }
    pub fn go(&self, accession: &str) -> Option<&str> {
        self.go.get(accession).map(|s| s as &str)
    }
    pub fn enzyme(&self, accession: &str) -> bool {
        self.enz.contains(accession)
    }
}

pub fn load<T: AsRef<Path>>(path: T) -> io::Result<Annotation> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);

    let mut current_ac = String::default();
    let mut kw = String::default();
    let mut go = String::default();

    let mut ann = Annotation::default();

    for line in reader.lines() {
        let line = line?;

        if line.starts_with("AC") {
            ann.kw.insert(current_ac.clone(), kw);
            ann.go.insert(current_ac, go);

            kw = String::default();
            go = String::default();
            current_ac = line[5..11].into();
        } else if line.starts_with("KW") {
            kw.push_str(line.trim_start_matches("KW   "));
        } else if line.starts_with("DR   GO") {
            let s = line.split(';');
            go.push_str(s.skip(2).next().unwrap_or_default());
        } else if line.starts_with("DE            EC=") {
            ann.enz.insert(current_ac.clone());
        }
    }
    Ok(ann)
}
