use super::*;

/// A collection of [`Filtered`] datasets, aggregated
/// on the peptide level
pub struct Aggregate {
    pub proteins: HashMap<String, Vec<Option<Protein>>>,
    pub experiments: Vec<String>,
}

impl Default for Aggregate {
    fn default() -> Aggregate {
        Aggregate {
            proteins: HashMap::default(),
            experiments: Vec::default(),
        }
    }
}

impl FromIterator<Filtered> for Aggregate {
    fn from_iter<I: IntoIterator<Item = Filtered>>(iter: I) -> Self {
        let mut agg = Aggregate::default();
        for i in iter {
            agg.merge(&i);
        }
        agg
    }
}

impl Aggregate {
    /// Add a [`Filtered`] dataset into the aggregate dataset
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// # use cimage::*;
    /// # let g1 = Raw::load("./data/exp1_combined_dta.txt").unwrap().group();
    /// # let g2 = Raw::load("./data/exp2_combined_dta.txt").unwrap().group();
    ///
    /// let filters = Filter::new()
    ///     .add_protein_filter(ProteinFilter::Reverse)
    ///     .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::CV(0.6)));
    ///
    /// let f1 = g1.filter(&filters);
    /// let f2 = g2.filter(&filters);
    ///
    /// let mut agg = Aggregate::default();
    ///
    /// agg.merge(&f1);
    /// agg.merge(&f2);
    /// // agg now contains the full outer join of rep1 and rep2
    ///
    /// ```
    pub fn merge(&mut self, dataset: &Filtered) {
        // Union of current and new uniprot accessions
        let accs = self
            .proteins
            .keys()
            .chain(dataset.proteins.keys())
            .cloned()
            .collect::<HashSet<_>>();

        let n = self.experiments.len();

        for acc in accs {
            self.proteins
                .entry(acc.clone())
                .or_insert_with(|| (0..n).map(|_| None).collect())
                .push(dataset.proteins.get(acc.as_ref() as &str).cloned());
        }

        self.experiments.push(dataset.path.clone());
    }

    /// Condense an aggregate dataset, taking the median ratio for each peptide
    /// in each experiment.
    ///
    /// At the initial [`Grouped`] stage, this should be equivalent to collapsing by peptide
    ///
    /// Each [`Peptide`] struct in the returned [`Grouped`] will contain a
    /// Vec<Option<f32>> for ratios of length `self.experiments.len()`; each
    /// entry in this vector represents the median ratio for the given peptide
    /// sequence in the corresponding experiment, where the index of the item
    /// in `ratios` == index of experiment in `self.experiments`
    pub fn condense<P: AsRef<str>>(self, path: P, collapse: bool) -> Grouped {
        let mut proteins = HashMap::new();
        for (acc, mut experiments) in self.proteins {
            // Hashmap for mapping peptide residue to a constructed Peptide,
            // where the ratios contained within the peptide represent the
            // median ratio across the experiments contained in `self`
            let mut map: HashMap<Residue, Peptide> = HashMap::new();
            let mut desc = "";
            // Iterate over the Vec of Option<Protein>, where each entry
            // in the vector represents data from one of the constituent
            // experiments
            //
            // If ex == None, it means the protein was not detected in the
            // experiment, or was later filtered out
            let n = experiments.len();
            for (idx, ex) in experiments.iter_mut().enumerate() {
                // Protein was detected in experiment #idx
                if let Some(protein) = ex.as_mut() {
                    if collapse {
                        // protein.collapse_redundant_sites();
                    }
                    desc = &protein.description;
                    for peptide in &protein.peptides {
                        // Insert a vector of Nones with length equal to number of experiments
                        // if this is the first time we've encountered this
                        // sequence.
                        let data = map.entry(peptide.residue).or_insert_with(|| Peptide {
                            sequence: peptide.sequence.clone(),
                            residue: peptide.residue,
                            ms2: 0,
                            ratios: (0..n).map(|_| None).collect::<Vec<_>>(),
                        });
                        // Set the value for this experiment
                        data.ratios[idx] = peptide.median_ratio();

                        // Add in more ms2 events
                        data.ms2 += peptide.ms2;
                    }
                }
            }

            // Collect all of the constructed peptides into a new `Protein`
            // struct
            let mut prot = map
                .into_iter()
                .map(|(_, peptide)| peptide)
                .collect::<Protein>();

            // FromIterator doesn't carry over acc/description, so manually move
            prot.accession = acc.clone();
            prot.description = desc.into();
            proteins.insert(acc, prot);
        }

        Grouped {
            proteins,
            path: path.as_ref().into(),
        }
    }
}
