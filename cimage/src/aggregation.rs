use super::*;

use std::collections::HashMap;
use std::fs;
use std::io::{self, prelude::*};
use std::path::Path;

#[derive(Default, Debug)]
struct FilteredPeptide {
    residue: Residue,
    ms2: usize,
    sequence: String,
    desc: String,
    ratios: Vec<Option<f32>>,
}

pub struct PeptideCollection {
    peptides: HashMap<String, Vec<FilteredPeptide>>,
    experiments: Vec<String>,
}

/// The scope of the rule 1 filter that can affect 20s
#[derive(Copy, Clone, PartialEq)]
pub enum FilterScope {
    /// Look at ratios collected from compounds across the *entire* dataset
    Dataset,
    /// Collect only ratios associated with the same chemotype group as the
    /// compound that we are currently filtering
    Chemotype,
}

impl PeptideCollection {
    pub fn counts(&self) -> HashMap<String, HashMap<Residue, SecondPassFilter>> {
        let mut map = HashMap::new();
        for (acc, peptides) in &self.peptides {
            map.insert(
                acc.clone(),
                peptides
                    .iter()
                    .map(|pep| (pep.residue, SecondPassFilter::from_ratios(&pep.ratios)))
                    .collect::<HashMap<Residue, SecondPassFilter>>(),
            );
        }
        map
    }

    // pub fn reverse(&self) -> Vec<
    pub fn second_pass_filter(
        &self,
        expt: &String,
        chemotype_map: &HashMap<String, String>,
        acc: &str,
        residue: Residue,
    ) -> SecondPassFilter {
        if !self.experiments.contains(expt) {
            dbg!(&self.experiments);
            panic!("Experiment not found in aggregate set! {}", expt);
        }

        let mut indices = Vec::new();

        let chemotype = chemotype_map
            .get(expt)
            .expect(&format!("missing {} from chemotype map", expt));
        for (cmpd, _) in chemotype_map.iter().filter(|(_, v)| v == &chemotype) {
            let mut x = None;
            for (idx, name) in self.experiments.iter().enumerate() {
                if name == cmpd {
                    x = Some(idx);
                    break;
                }
            }
            match x {
                Some(idx) => indices.push(idx),
                None => panic!(
                    "This is very bad, there is some eldritch horror lurking in the codebase"
                ),
            }
        }

        // println!("SPF: {}, {:?}, {}", expt,indices,acc);

        // for (i, x) in self.experiments.iter().enumerate() {
        //     if x == expt {
        //         idx = i;
        //         break;
        //     }
        // }

        let mut spf = SecondPassFilter::default();
        if let Some(peptides) = self.peptides.get(acc) {
            for pep in peptides.iter().filter(|p| p.residue == residue) {
                for i in 0..pep.ratios.len() {
                    if &self.experiments[i] == expt {
                        continue;
                    }

                    if let Some(x) = pep.ratios[i] {
                        if x == 20.0 {
                            spf.twenties += 1;
                        } else if x >= 4.0 {
                            if expt.contains("SF-21") && acc == "Q9NR28" && residue == 207 {
                                eprintln!("{} {}", self.experiments[i], x);
                            }
                            spf.liganded += 1;
                        } else {
                            spf.not_liganded += 1;
                        }
                    }
                }
            }
        }
        spf
    }
}

impl PeptideCollection {
    pub fn new(agg: &Aggregate) -> Self {
        let mut map = HashMap::new();

        for (acc, experiments) in &agg.proteins {
            if acc.contains("Reverse") {
                continue;
            }
            let mut desc = "";

            // Issue with duplicate sites being listed, since we're combining by sequence
            let mut site_to_ratios: HashMap<_, Vec<Vec<Option<f32>>>> = HashMap::new();
            let mut site_to_seq: HashMap<Residue, &str> = HashMap::new();
            let mut site_to_ms2: HashMap<Residue, usize> = HashMap::new();

            for (idx, ex) in experiments.iter().enumerate() {
                if let Some(protein) = ex {
                    desc = &protein.description;
                    for peptide in &protein.peptides {
                        let ratios = site_to_ratios.entry(peptide.residue).or_insert_with(|| {
                            (0..experiments.len())
                                .map(|_| Vec::new())
                                .collect::<Vec<_>>()
                        });

                        // Since we may have multiple sequences that map to the
                        // same residue in any given experiment (idx), we keep
                        // a running list of the ratios detected in each exp,
                        // we will calculate a median ratio once we've iterated
                        // over all of the experiments and peptides for this protein
                        ratios[idx].extend(peptide.ratios.iter());
                        site_to_seq.insert(peptide.residue, &peptide.sequence);
                        *site_to_ms2.entry(peptide.residue).or_insert(0) += peptide.ms2;
                    }
                }
            }

            let mut seq_to_sites: HashMap<String, Vec<Residue>> = HashMap::new();

            for (site, seq) in &site_to_seq {
                seq_to_sites
                    .entry(seq.replace('*', ""))
                    .or_insert_with(Vec::new)
                    .push(*site);
            }

            for (_, mut sites) in seq_to_sites {
                sites.sort();
                if sites.len() > 1 {
                    let keep_residue = sites.pop().unwrap();
                    for remove in sites {
                        // All of these unwraps should be gucci
                        let ms2 = site_to_ms2.remove(&remove).unwrap();
                        *site_to_ms2.entry(keep_residue).or_insert(0) += ms2;

                        let ratios = site_to_ratios.remove(&remove).unwrap();

                        let mutref = site_to_ratios.get_mut(&keep_residue).unwrap();

                        for (idx, v) in ratios.into_iter().enumerate() {
                            mutref[idx].extend(v.iter());
                        }
                    }
                }
            }

            // At this stage we've rearranged the data into a matrix of peptide
            // sequences and their ratios, which is essentially the final form
            // that will be presented in an Excel file, etc. It is at this stage
            // that we can apply filters across multiple chemotypes/compounds

            let mut peptides = Vec::new();

            for (residue, ratios) in site_to_ratios {
                // Take our matrix of ratios (experiment by ratios per peptide),
                // and condense it down into a vector of median ratios.
                // This is essentially one row in the CSV, if output on peptide level
                let condensed: Vec<Option<f32>> = ratios
                    .into_iter()
                    .map(|exp| {
                        let mut pep = Peptide::default();
                        pep.ratios = exp;
                        pep.median_ratio()
                    })
                    .collect();

                let mut twenties = 0;
                let mut liganded = 0;
                // let mut not_liganded = 0;

                for r in &condensed {
                    if let Some(x) = r {
                        if *x == 20.0 {
                            twenties += 1;
                        } else if *x >= 4.0 {
                            liganded += 1;
                        } else {
                            // not_liganded += 1;
                        }
                    }
                }

                // if not_liganded + twenties + liganded >= 1 {
                peptides.push(FilteredPeptide {
                    residue,
                    sequence: site_to_seq
                        .get(&residue)
                        .copied()
                        .unwrap_or_default()
                        .into(),
                    desc: desc.into(),
                    ms2: site_to_ms2.get(&residue).copied().unwrap_or_default(),
                    ratios: condensed,
                });
            }
            map.insert(acc.clone(), peptides);
        }
        PeptideCollection {
            peptides: map,
            experiments: agg.experiments.clone(),
        }
    }

    pub fn write_functions<P: AsRef<Path>>(&self, path: P, keywords: P) -> io::Result<()> {
        let mut f = fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)?;

        let ann = uniprot::kw::load(keywords)?;
        writeln!(f, "accession\tdescription\tkeywords\tgo term\tenzyme")?;

        for (acc, peptides) in &self.peptides {
            let mut desc = "";
            if peptides.len() == 0 {
                continue;
            }
            for p in peptides {
                desc = &p.desc;

                if desc != "" {
                    break;
                }
            }

            let kw = ann.keyword(acc).unwrap_or_default();
            let go = ann.go(acc).unwrap_or_default();
            let enz = ann.enzyme(acc);
            writeln!(f, "{}\t{}\t{}\t{}\t{}", acc, desc, kw, go, enz)?;
        }

        Ok(())
    }

    pub fn write_peptides<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut f = fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(&path)
            .expect(&format!("failed to open {}", path.as_ref().display()));

        writeln!(
            f,
            "accession\tgene_name\tdescription\tsequence\tsite\tmax_ratio\taverage_ratio\t{}",
            self.experiments
                .iter()
                .map(|p| p.clone())
                .collect::<Vec<String>>()
                .join("\t")
        )?;

        for (acc, peptides) in &self.peptides {
            for peptide in peptides {
                let mut max = 0.;
                let mut sum = 0.;
                let mut n = 0;
                for r in peptide.ratios.iter().filter_map(|x| *x) {
                    sum += r;
                    if r > max {
                        max = r;
                    }
                    n += 1;
                }

                if n == 0 {
                    // continue;
                }
                writeln!(
                    f,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    acc,
                    peptide.desc.split(' ').next().unwrap_or_default(),
                    peptide.desc,
                    peptide.sequence,
                    peptide.residue,
                    max,
                    sum / (n as f32),
                    peptide
                        .ratios
                        .iter()
                        .map(|x| x.map(|y| y.to_string()).unwrap_or_default())
                        .collect::<Vec<String>>()
                        .join("\t")
                )?;
            }
        }
        Ok(())
    }

    ///
    /// Args:
    /// * chemotype_map: map a compound name to a chemotype, e.g. (HA17, Scouts) or (DAPG1, DAPG)
    pub fn combine_into_chemotypes(
        &self,
        chemotype_map: HashMap<String, String>,
    ) -> PeptideCollection {
        let mut reorg = PeptideCollection {
            peptides: HashMap::new(),
            experiments: chemotype_map
                .values()
                .cloned()
                .collect::<HashSet<_>>()
                .into_iter()
                .collect(),
        };
        reorg.experiments.sort();

        // dbg!(&chemotype_map);
        let compounds = self
            .experiments
            .iter()
            .map(|x| x.clone())
            .collect::<Vec<String>>();

        // dbg!(&compounds);
        let chemotypes = reorg
            .experiments
            .iter()
            .enumerate()
            .map(|(idx, chemotype)| (chemotype.clone(), idx))
            .collect::<HashMap<String, usize>>();
        // dbg!(&chemotypes);
        for (acc, peptides) in &self.peptides {
            let mut aggregated_peptides = Vec::new();

            for peptide in peptides.iter() {
                let mut new_peptide = FilteredPeptide::default();
                new_peptide.ms2 = peptide.ms2;
                new_peptide.desc = peptide.desc.clone();
                new_peptide.sequence = peptide.sequence.clone();
                new_peptide.residue = peptide.residue;
                new_peptide.ratios = (0..chemotypes.len()).map(|_| None).collect::<Vec<_>>();
                for (idx, ratio) in peptide.ratios.iter().enumerate() {
                    if let Some(r) = ratio {
                        let chemotype_name =
                            chemotype_map.get(&compounds[idx]).unwrap_or_else(|| {
                                panic!("Chemotype for {} not found!", &compounds[idx])
                            });
                        let &chemotype_idx = chemotypes.get(chemotype_name).unwrap();
                        new_peptide.ratios[chemotype_idx] = match new_peptide.ratios[chemotype_idx]
                        {
                            Some(v) => {
                                if *r > v {
                                    Some(*r)
                                } else {
                                    Some(v)
                                }
                            }
                            None => Some(*r),
                        }
                    }
                }
                aggregated_peptides.push(new_peptide);
            }
            reorg.peptides.insert(acc.clone(), aggregated_peptides);
        }
        reorg
    }

    /// This function needs to be called before the [`PeptideCollection`]
    /// has been combined by chemotype
    pub fn extract_chemotype(
        &self,
        chemotype: &str,
        chemotype_map: &HashMap<String, String>,
    ) -> PeptideCollection {
        // We want to have all of the individuals compounds sorted alphabetically
        let mut indices = self
            .experiments
            .iter()
            .enumerate()
            .filter(|(_, cmpd)| match chemotype_map.get(cmpd.clone()) {
                Some(chem) if chem == chemotype => true,
                _ => false,
            })
            .collect::<Vec<(usize, &String)>>();

        indices.sort_by_key(|t| t.1);

        let mut pc = PeptideCollection {
            experiments: indices.iter().map(|(_, cmpd)| *cmpd).cloned().collect(),
            peptides: HashMap::new(),
        };

        for (acc, peptides) in &self.peptides {
            let mut pass_thru = Vec::new();
            for peptide in peptides.iter() {
                let mut new_peptide = FilteredPeptide::default();
                new_peptide.ms2 = peptide.ms2;
                new_peptide.desc = peptide.desc.clone();
                new_peptide.sequence = peptide.sequence.clone();
                new_peptide.residue = peptide.residue;
                new_peptide.ratios = (0..indices.len()).map(|_| None).collect::<Vec<_>>();

                // i represents the index into the `experiments` field of the new PeptideCollection,
                // and j represents the index into the `experiments` field of `self`
                for (i, j) in indices.iter().map(|(idx, _)| idx).enumerate() {
                    new_peptide.ratios[i] = peptide.ratios[*j];
                }
                pass_thru.push(new_peptide);
            }
            pc.peptides.insert(acc.clone(), pass_thru);
        }
        pc
    }

    pub fn write_proteins<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut f = fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)?;

        writeln!(
            f,
            "accession\tdescription\tmax_ratio\taverage_ratio\t{}",
            self.experiments
                .iter()
                .map(|p| p.clone())
                .collect::<Vec<String>>()
                .join("\t")
        )?;

        for (acc, peptides) in &self.peptides {
            let mut selected_ratios: Vec<f32> = (0..self.experiments.len()).map(|_| 0f32).collect();
            let mut desc = "";
            let mut n = 0;
            let mut max = 0.;
            let mut sum = 0.;
            for peptide in peptides {
                desc = &peptide.desc;
                for (idx, r) in peptide.ratios.iter().enumerate() {
                    // take max ratio of all the peptides
                    if let Some(rat) = r {
                        n += 1;
                        if *rat > selected_ratios[idx] {
                            selected_ratios[idx] = *rat;
                        }
                        if *rat > max {
                            max = *rat;
                        }
                        sum += *rat;
                    }
                }
            }

            // not detected, or not ligandend
            if n == 0
            /* || selected_ratios.iter().filter(|&&x| x >= 4.0).count() == 0 */
            {
                // continue;
            }

            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                acc,
                desc,
                max,
                sum / (n as f32),
                selected_ratios
                    .into_iter()
                    .map(Self::eliminate_zeros)
                    .collect::<Vec<String>>()
                    .join("\t")
            )?;
        }
        Ok(())
    }

    fn eliminate_zeros(s: f32) -> String {
        if s == 0.0 {
            String::default()
        } else {
            s.to_string()
        }
    }
}
