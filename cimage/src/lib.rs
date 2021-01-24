//! A high-performance Rust library for working with data from the
//! Cravatt lab's cimage proteomics program.
//!
//! This library's API is based around several types that form
//! a data analysis pipeline.
//!
//! Raw cimage_combined output is transformed into a [`Raw`] object,
//! which can then be grouped by protein into a [`Grouped`] struct.
//!
//! ```rust,ignore
//! #  use cimage::*;
//! let g = Raw::load("./data/exp1_combined_dta.txt").unwrap().group();
//! ```
//!
//! [`Grouped`] objects can be transformed into a [`Filtered`] struct
//! by applying a collection of [`Filter`]'s on both the protein and
//! peptide level, including destructively filtering out spurious ratios
//!
//! ```rust,ignore
//! # use cimage::*;
//!
//! let g1 = Raw::load("./data/exp1_combined_dta.txt").unwrap().group();
//! let g2 = Raw::load("./data/exp2_combined_dta.txt").unwrap().group();
//!
//! let filters = Filter::new()
//!     .add_protein_filter(ProteinFilter::Reverse)
//!     .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::CV(0.6)));
//!
//! let f1 = g1.filter(&filters);
//! let f2 = g2.filter(&filters);
//!
//! ```
//!
//! [`Filtered`] datasets can then be joined together into an [`Aggregate`]
//! dataset, where peptide values from each constituent experiment are mereged
//! by median value. An [`Aggregate`] dataset can further be condensed back
//! into a [`Grouped`] struct for further filtering and reaggregation into
//! larger and more complex aggregate of aggregate datasets, or it can be
//! written to a file
//!
//! ```rust,ignore
//! # use cimage::*;
//!
//! let g1 = Raw::load("./data/exp1_combined_dta.txt").unwrap().group();
//! let g2 = Raw::load("./data/exp2_combined_dta.txt").unwrap().group();
//!
//! let filters = Filter::new()
//!     .add_protein_filter(ProteinFilter::Reverse)
//!     .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::CV(0.6)));
//!
//! let f1 = g1.filter(&filters);
//! let f2 = g2.filter(&filters);
//!
//! // Aggregate implements FromIterator<Filtered>
//! let agg = vec![f1, f2].into_iter().collect::<Aggregate>();
//!
//! // Now condense it back into a grouped object, combining individual
//! // experiments on the peptide level
//! let grouped: Grouped = agg.condense("path");
//!
//! ```

use std::char;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{self, BufRead, BufReader};
use std::iter::FromIterator;
use std::path::{Path, PathBuf};

mod aggregate;
mod aggregation;
mod filter;
mod parser;
mod peptide;
mod protein;
mod stats;

pub use aggregate::*;
pub use aggregation::*;
pub use filter::{Filter, PeptideFilter, ProteinFilter, RatioFilter, SecondPassFilter};
pub use parser::Raw;
pub use peptide::Peptide;
pub use protein::Protein;

pub type Residue = u16;

/// Cimage data that has been grouped by UniProt ID and peptide sequence
pub struct Grouped {
    pub proteins: HashMap<String, Protein>,
    pub path: String,
}

pub struct Filtered {
    pub proteins: HashMap<String, Protein>,
    pub path: String,
}

impl Grouped {
    pub fn filter<'a>(self, filters: &Filter<'a>) -> Filtered {
        let proteins = self
            .proteins
            .into_iter()
            .filter_map(|(acc, mut protein)| {
                // protein.collapse_redundant_sites();
                Some((acc, filters.filter(protein)?))
            })
            .collect::<HashMap<String, Protein>>();

        Filtered {
            proteins,
            path: self.path,
        }
    }
}

impl Filtered {
    pub fn write<P: AsRef<std::path::Path>>(&self, p: P) -> std::io::Result<()> {
        use std::io::prelude::*;
        let mut f = std::fs::File::create(p)?;
        writeln!(f, "identifier\tms2\tratio")?;
        for (acc, prot) in &self.proteins {
            for peptide in &prot.peptides {
                if let Some(r) = peptide.median_ratio() {
                    writeln!(f, "{}_{}\t{}\t{}", acc, peptide.residue, peptide.ms2, r)?;
                }
                // for ratio in peptide.ratios.iter().copied().filter_map(|f| f) {
                //     writeln!(f, "{}_{}\t{}\t{}", acc, peptide.residue, peptide.ms2, ratio)?;
                // }
            }
        }
        Ok(())
    }
}

unsafe impl Send for Filtered {}
unsafe impl Sync for Filtered {}
