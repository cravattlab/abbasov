use super::*;

#[derive(Clone, Debug, Default, PartialEq)]
/// Container for [`Peptide`] sequences and their quantified ratios
///
/// Allows constant-time lookup by amino-acid sequence
pub struct Protein {
    pub accession: String,
    pub description: String,
    pub peptides: Vec<Peptide>,

    /// Constant-time indexing into peptides array - do we even need the Vec?
    ///
    /// Should look into average length of peptide Vec to determine
    /// if the extra memory usage here is worth the speedup.
    pub map: HashMap<Residue, usize>,
}

impl Protein {
    pub fn new(accession: String, description: String) -> Protein {
        Protein {
            accession,
            description,
            peptides: Vec::new(),
            map: HashMap::new(),
        }
    }

    /// Return a [`HashMap`] mapping peptide sequences to a [`Peptide`] struct
    pub fn as_map(&self) -> HashMap<&str, &Peptide> {
        self.peptides
            .iter()
            .map(|pep| (pep.sequence.as_ref(), pep))
            .collect()
    }

    /// Add a ratio for a given peptide sequence to the [`Protein`]
    ///
    /// If the protein already contains a peptide with the same sequence as
    /// the peptide being added, the new ratio will be appended to
    /// the existing peptide match's ratios.
    pub fn add_ratio(&mut self, residue: Residue, seq: &str, ratio: Option<f32>, ms2: usize) {
        match self.map.get(&residue) {
            Some(idx) => {
                self.peptides[*idx].ratios.push(ratio);
                self.peptides[*idx].ms2 += ms2;

                // if self.peptides[*idx].ms2 == 0 {
                //     dbg!(&self.peptides[*idx]);
                // }
            }
            None => {
                let idx = self.peptides.len();
                self.peptides.push(Peptide {
                    sequence: seq.into(),
                    residue: residue,
                    ms2,
                    ratios: vec![ratio],
                });
                self.map.insert(residue, idx);
            }
        }
    }

    pub fn collapse_redundant_sites(&mut self) {
        let mut count: HashMap<String, Vec<_>> = HashMap::new();
        for (idx, p) in self.peptides.iter().enumerate() {
            count
                .entry(p.sequence.replace('*', ""))
                .or_insert_with(Vec::new)
                .push((idx, p.residue));
        }

        for (_, mut indices) in count {
            indices.sort_by(|a, b| a.1.cmp(&b.1));
            if indices.len() > 1 {
                let (keep, resi) = indices.pop().unwrap();
                // let resi = self.peptides[keep].residue;
                for (remove_idx, _) in indices {
                    let mut removed: Peptide =
                        std::mem::replace(&mut self.peptides[remove_idx], Peptide::default());
                    removed.residue = resi;
                    self.add_peptide(removed);
                }
            }
        }

        // Final step, remove the default peptide that we inserted to keep
        // indices in place, and then remake our hashmap
        let def = Peptide::default();
        self.peptides = self.peptides.drain(..).filter(|p| *p != def).collect();
        self.map = self
            .peptides
            .iter()
            .enumerate()
            .map(|(idx, pep)| (pep.residue, idx))
            .collect();
    }

    /// Add a new [`Peptide`] struct to the [`Protein`]
    ///
    /// If the protein already contains a peptide with the same sequence as
    /// the peptide being added, the new peptide's ratios will be appended to
    /// the existing peptide match's ratios.
    pub fn add_peptide(&mut self, peptide: Peptide) {
        match self.map.get(&peptide.residue) {
            Some(idx) => {
                self.peptides[*idx].ratios.extend(peptide.ratios);
                self.peptides[*idx].ms2 += peptide.ms2;
            }
            None => {
                let idx = self.peptides.len();
                self.map.insert(peptide.residue, idx);
                self.peptides.push(peptide);
            }
        }
    }

    pub fn get(&self, residue: Residue) -> Option<&Peptide> {
        self.peptides.get(*self.map.get(&residue)?)
    }

    pub fn get_mut(&mut self, residue: Residue) -> Option<&mut Peptide> {
        self.peptides.get_mut(*self.map.get(&residue)?)
    }

    /// Similar to a [`HashMap`]'s `Entry`, return mutable reference to
    /// existing peptide, or create a new peptide entry for the sequence and
    /// return a mutable reference to it
    // pub fn get_or_insert(&mut self, residue: Residue) -> &mut Peptide {
    //     match self.map.get(&residue) {
    //         Some(idx) => &mut self.peptides[*idx],
    //         None => {
    //             let idx = self.peptides.len();
    //             self.peptides.push(Peptide {
    //                 sequence: seq.into(),
    //                 residue: None,
    //                 ms2: 1,
    //                 ratios: Vec::new(),
    //             });
    //             self.map.insert(seq.into(), idx);
    //             &mut self.peptides[idx]
    //         }
    //     }
    // }

    pub fn median_ratio(&self) -> Option<f32> {
        let medians = self
            .peptides
            .iter()
            .map(Peptide::median_ratio)
            .collect::<Vec<Option<f32>>>();
        let mock = Peptide {
            sequence: String::default(),
            residue: 0,
            ms2: 0,
            ratios: medians,
        };
        mock.median_ratio()
    }

    pub fn spectral_counts(&self) -> usize {
        self.peptides.iter().fold(0, |acc, x| {
            acc + x.ratios.iter().filter(|r| r.is_some()).count()
        })
    }
}

impl FromIterator<Peptide> for Protein {
    fn from_iter<I: IntoIterator<Item = Peptide>>(iter: I) -> Self {
        let mut protein = Protein::default();
        for peptide in iter {
            protein.add_peptide(peptide);
        }
        protein
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn add_ratios() {
        let mut prot = Protein::default();
        macro_rules! add {
            ($s:expr, $site:expr, $f:expr) => {
                prot.add_ratio($site, $s.into(), Some($f), 0)
            };
        }
        add!("MRL", 0, 1.);
        add!("MRL", 0, 2.);
        add!("MRL", 0, 3.);
        add!("MEHQLL", 1, 20.0);
        assert_eq!(prot.peptides.len(), 2);
        assert_eq!(
            prot.get(0).unwrap().ratios,
            vec![Some(1.), Some(2.), Some(3.)]
        );
        assert_eq!(prot.get(1).unwrap().ratios, vec![Some(20.0)]);
        assert_eq!(prot.spectral_counts(), 4);
    }

    #[test]
    fn removed_redundant_sites() {
        let mut prot = Protein::default();
        macro_rules! add {
            ($s:expr, $site:expr, $f:expr) => {
                prot.add_ratio($site, $s.into(), Some($f), 0)
            };
        }
        add!("K.LQFGSQPQVYNDFLDIMKEFK*SQSIDTPGVISR.V", 155, 1.);
        add!("K.LQFGSQPQVYNDFLDIMKEFK*SQSIDTPGVISR.V", 155, 5.);
        add!("K.LQFGSQPQVYNDFLDIMKEFK*SQSIDTPGVISR.V", 152, 2.);
        add!("K.LQFGSQPQVYNDFLDIMKEFK*SQSIDTPGVISR.V", 152, 3.5);
        add!("R.LK*VEDALSYLDQVK.L", 122, 20.0);

        assert_eq!(prot.peptides.len(), 3);
        assert_eq!(prot.get(152).unwrap().ratios, vec![Some(2.), Some(3.5)]);

        prot.collapse_redundant_sites();
        assert_eq!(prot.peptides.len(), 2);
        assert_eq!(
            prot.get(152).unwrap().ratios,
            vec![Some(2.), Some(3.5), Some(1.), Some(5.)]
        );
        assert_eq!(prot.get(155), None);
    }
}
