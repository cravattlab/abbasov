use super::*;

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub enum ProteinFilter {
    Reverse,
    Spectra(usize),
    ExcludeMatch(String),
}

/// A collection of information about a set of ratios
/// INVARIANT: The sum of twenties + liganded + not_liganded must equal
/// the total number of quantified ratios
#[derive(Copy, Clone, Default, Debug)]
pub struct SecondPassFilter {
    /// How many ratios are exactly 20
    pub twenties: usize,
    /// How many ratios are >= 4, and < 20
    pub liganded: usize,
    /// How many ratios are < 4
    pub not_liganded: usize,
}

impl SecondPassFilter {
    pub fn from_ratios(ratios: &[Option<f32>]) -> Self {
        let mut s = SecondPassFilter::default();
        for r in ratios {
            if let Some(x) = *r {
                if x == 20.0 {
                    s.twenties += 1;
                } else if x >= 4.0 {
                    s.liganded += 1;
                } else {
                    s.not_liganded += 1;
                }
            }
        }
        s
    }

    pub fn total(&self) -> usize {
        self.twenties + self.liganded + self.not_liganded
    }
}

#[derive(Clone)]
pub enum PeptideFilter<'s> {
    HalfTryptic,
    ExcludeMatch(&'s str),
    Ratios(RatioFilter),
    Ms2(usize),
    // SPF2(&'s HashMap<String, HashMap<Residue, SecondPassFilter>>),
    SecondPassFilter(&'s PeptideCollection, &'s HashMap<String, String>, String),
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub enum RatioFilter {
    /// Filter out 20s from sets of ratios with high stdev
    CV(f32),
    /// Number of quantified ratios must be >= usize
    Count(usize),
    Spurious,
}

#[derive(Clone)]
pub struct Filter<'a> {
    peptide_filters: Vec<PeptideFilter<'a>>,
    protein_filters: Vec<ProteinFilter>,
}

impl<'a> Filter<'a> {
    pub fn new() -> Self {
        Filter {
            peptide_filters: Vec::new(),
            protein_filters: Vec::new(),
        }
    }

    pub fn add_peptide_filter(mut self, f: PeptideFilter<'a>) -> Self {
        self.peptide_filters.push(f);
        self
    }

    pub fn add_protein_filter(mut self, f: ProteinFilter) -> Self {
        self.protein_filters.push(f);
        self
    }

    pub fn filter(&self, protein: Protein) -> Option<Protein> {
        use PeptideFilter::*;
        use RatioFilter::*;

        let mut filtered = Vec::new();

        // Check for reverse protein ID filter first - no need to spend
        // time filtering peptides if we have reverse accession
        for filter in &self.protein_filters {
            if let ProteinFilter::Reverse = filter {}
            match filter {
                ProteinFilter::Reverse => {
                    if protein.accession.contains("Reverse") {
                        return None;
                    }
                }
                ProteinFilter::ExcludeMatch(s) => {
                    if protein.description.contains(s) {
                        return None;
                    }
                }
                _ => {}
            }
        }
        let mut trigger = false;
        if !self.peptide_filters.is_empty() {
            for mut peptide in protein.peptides {
                let mut pass = true;
                if peptide.ratios.iter().all(|x| x.is_none()) {
                    // Don't need to filter an empty peptide
                    continue;
                }
                // if trigger &&  peptide.non_zeroes().len() < 2 {
                //     pass = false;
                //     break;
                // }
                for filter in &self.peptide_filters {
                    match filter {
                        HalfTryptic => {
                            if !peptide.is_not_half_tryptic() {
                                pass = false;
                                break;
                            }
                        }
                        ExcludeMatch(seq) => {
                            if peptide.sequence.contains(seq) {
                                pass = false;
                                break;
                            }
                        }
                        Ratios(rf) => match rf {
                            CV(cv) => peptide.cv_filter(*cv),
                            Count(count) => {
                                if peptide.non_zeroes().len() < *count {
                                    pass = false;
                                    break;
                                }
                            }

                            Spurious => peptide.spurious_filter(),
                        },
                        Ms2(cutoff) => {
                            if peptide.ms2 < *cutoff {
                                // For now, just remove twenties
                                peptide.remove_twenties();
                            }
                        }
                        SecondPassFilter(pc, chemotype_map, expt) => {
                            // Ok this is funky.. but the second pass filter should contain
                            // all of the information about the count of liganded, twenties,
                            // and non-liganded ratios for all other compounds within the
                            // same chemotype as [`expt`].
                            let spf = pc.second_pass_filter(
                                expt,
                                chemotype_map,
                                &protein.accession,
                                peptide.residue,
                            );
                            let this = filter::SecondPassFilter::from_ratios(&peptide.ratios);

                            // remove cases where for this compound, all of the reported
                            // replicates have a median ratio of 20, but there are no non-
                            // 20's in any other compound aggregate set
                            //
                            // Modified to be as stated:
                            //
                            // If all ratios within a set of replicates for a
                            // given compound ([`expt`]) are 20s, and no other
                            // compound within our chemotype has a single
                            // ligand event for this peptide, then we remove
                            // 20s
                            if spf.liganded == 0 {
                                peptide.remove_twenties();
                            }

                            // if this.twenties == 1 && this.liganded == 0 {
                            //     peptide.remove_twenties();
                            // }

                            if spf.liganded == 1 && this.twenties == 1 && this.liganded == 0 {
                                peptide.remove_twenties();
                            }

                            // Remake our second pass filter, this time
                            // collecting ratios from the *entire* dataset, not
                            // just our chemotype group. We use this to ensure
                            // that every ligandable ratio is quantified somewhere
                            // else in the dataset
                            // let spf = pc.second_pass_filter(
                            //     expt,
                            //     chemotype_map,
                            //     &protein.accession,
                            //     peptide.residue,
                            //     FilterScope::Dataset,
                            // );

                            // "I’d also recommend that we require for a site to
                            // be quantified in at least 5 distinct data sets for
                            // interpretation – this will remove some events with
                            // very sparse coverage where the single ligandability
                            // event is also borderline quality (we could remove these
                            // events manually too, but I didn’t notice many convincing
                            // liganding events for sites with very sparse coverage; e.g.,
                            // fewer than 5 quantification events across the entire data set"
                            if spf.total() < 5 && !trigger {
                                pass = false;
                                if trigger {
                                    // eprintln!("5 {:?}", peptide);
                                }
                                break;
                            }
                        }
                    }
                }
                if pass {
                    filtered.push(peptide);
                } else {
                    // eprintln!("dropping peptide");
                    if trigger {
                        // eprintln!("drop {:?}", peptide);
                    }
                }
            }
        } else {
            filtered = protein.peptides;
        }

        // No peptides met the filters, so the protein is removed
        if filtered.is_empty() {
            return None;
        }

        // FromIterator doesn't carry over acc/description, so manually move
        let mut constructed = filtered.into_iter().collect::<Protein>();
        constructed.accession = protein.accession;
        constructed.description = protein.description;
        // constructed.map = constructed.peptides.iter().enumerate().map(|(idx, p)| (p.residue, idx)).collect();

        for filter in &self.protein_filters {
            match filter {
                ProteinFilter::Spectra(count) => {
                    if constructed.spectral_counts() < *count {
                        return None;
                    }
                }
                _ => {}
            }
        }

        Some(constructed)
    }
}

#[cfg(test)]
mod test {
    use super::PeptideFilter::*;
    use super::RatioFilter::*;
    use super::*;

    #[test]
    fn filter_20s() {
        let mut pepa = Peptide::default();
        pepa.ratios
            .extend((0..10).map(|i| Some(1.0 + (i as f32 / 10.0))));
        pepa.ratios.push(Some(20.0));
        pepa.sequence = String::from("A");
        pepa.residue = 0;

        let mut pepb = Peptide::default();
        pepb.ratios.extend((0..5).map(|_| Some(8.0)));
        pepb.ratios.push(Some(20.0));
        pepb.sequence = String::from("B");
        pepb.residue = 1;

        let mut prot = Protein::default();
        prot.add_peptide(pepa);
        prot.add_peptide(pepb);

        let a = prot.get(0).unwrap();
        let b = prot.get(1).unwrap();

        let filter = Filter::new()
            .add_peptide_filter(Ratios(CV(0.6)))
            .add_peptide_filter(Ratios(Spurious));

        let f = filter.filter(prot).unwrap();
        let a = f.get(0).unwrap();
        let b = f.get(1).unwrap();

        assert_eq!(a.ratios.iter().filter(|&r| *r == Some(20.0)).count(), 0);
        assert_eq!(b.ratios.iter().filter(|&r| *r == Some(20.0)).count(), 1);
    }

    #[test]
    fn reverse() {
        let filter = Filter::new().add_protein_filter(ProteinFilter::Reverse);

        assert_eq!(
            filter.filter(Protein::new(
                String::from("Reverse_QQQQ"),
                String::default()
            )),
            None
        );
    }
}
