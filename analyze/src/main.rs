use cimage;
use cimage::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use uniprot::fasta::Fasta;

fn scan_directory<P: AsRef<Path>>(path: P, fasta: &Fasta) -> io::Result<Vec<Grouped>> {
    let mut v = Vec::new();
    for entry in fs::read_dir(path)? {
        let path = entry?.path();
        if path.is_file() {
            if let Some(ext) = path.extension() {
                if ext == "combined" {
                    let mut dta = path.clone();
                    dta.set_extension("dtaselect");
                    v.push((path, dta));
                }
            }
        }
    }
    Ok(v.into_par_iter()
        .map(|p| Raw::load(p.0, p.1, fasta).unwrap().group())
        .collect())
}

struct Details {
    cmpd: String,
    conc: String,
}

/// Given a file name, return the name of the compound
fn extract_meta_data(s: &str) -> Option<Details> {
    let mut x = if s.starts_with("MIKA") {
        s[5..].split('_')
    } else {
        s.split('_')
    };
    let cmpd = x.next()?.into();
    let conc = x.next()?.into();
    Some(Details { cmpd, conc })
}

/// Generate a mapping of `CompoungName_Concentration` -> `Chemotype` names
fn generate_chemotype_map<P: AsRef<Path>>(path: P) -> io::Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    for entry in fs::read_dir(path.as_ref().clone())? {
        let path = entry?.path();
        if path.is_dir() {
            let chemotype_name = path.file_name().unwrap().to_str().unwrap().to_string();
            for entry in fs::read_dir(&path)? {
                let path = entry?.path();
                if path.is_file() {
                    if let Some(ext) = path.extension() {
                        if ext == "combined" {
                            let details =
                                extract_meta_data(path.file_name().unwrap().to_str().unwrap())
                                    .unwrap();

                            map.insert(
                                format!("{}_{}", details.cmpd, details.conc),
                                chemotype_name.clone(),
                            );
                        }
                    }
                }
            }
        }
    }
    Ok(map)
}

/// Search the parent directory, calling `scan_directory` on all subdirectories
///
/// Each subdirectory (chemotype) is fltered using `filters`, aggregated, and
/// saved into it's own file for further chemotype-specific analysis. The
/// aggregated datasets are then condensed back into a `Grouped` dataset,
/// as if each folder was one experiment.
fn generate_funnel<P: AsRef<Path>>(
    path: P,
    ind_filters: &Filter,
    fasta: &Fasta,
) -> io::Result<PeptideCollection> {
    let mut v = Vec::new();
    for entry in fs::read_dir(path.as_ref().clone())? {
        let path = entry?.path();
        if path.is_dir() {
            let intermediate = scan_directory(&path, fasta)?
                .into_iter()
                .map(|g| g.filter(ind_filters))
                .collect::<Vec<_>>();

            let mut map = HashMap::new();
            for x in intermediate {
                let details = extract_meta_data(&x.path).unwrap();
                map.entry(format!("{}_{}", details.cmpd, details.conc))
                    .or_insert_with(Vec::new)
                    .push(x);
            }
            for (cmpd, data) in map {
                let agg = data.into_iter().collect::<Aggregate>();
                v.push((agg, cmpd));
            }
        }
    }
    v.sort_by(|a, b| a.1.cmp(&b.1));

    let mut a = Vec::new();
    for (agg, path) in v {
        a.push(agg.condense(&path, false).filter(
            &Filter::new().add_peptide_filter(PeptideFilter::Ratios(RatioFilter::CV(0.6))),
        ));
    }
    let a = a.into_iter().collect::<Aggregate>();
    Ok(PeptideCollection::new(&a))
}

/// Search the parent directory, calling `scan_directory` on all subdirectories
///
/// Each subdirectory (chemotype) is fltered using `filters`, aggregated, and
/// saved into it's own file for further chemotype-specific analysis. The
/// aggregated datasets are then condensed back into a `Grouped` dataset,
/// as if each folder was one experiment.
fn aggregate_compounds<P: AsRef<Path>>(
    path: P,
    ind_filters: &Filter,
    agg_filters: &PeptideCollection,
    fasta: &Fasta,
    chemotype_map: &HashMap<String, String>,
) -> io::Result<PeptideCollection> {
    let mut v = Vec::new();
    for entry in fs::read_dir(path.as_ref().clone())? {
        let path2 = entry?.path();
        if path2.is_dir() {
            dbg!(&path2.display());
            let intermediate = scan_directory(&path2, fasta)?
                .into_iter()
                .map(|g| g.filter(&ind_filters))
                .collect::<Vec<_>>();

            let mut map = HashMap::new();
            for x in intermediate {
                // println!("{:?}\t{:?}", x.path, count_quantified_lysines(&x));

                let details = extract_meta_data(&x.path).unwrap();
                map.entry(format!("{}_{}", details.cmpd, details.conc))
                    .or_insert_with(Vec::new)
                    .push(x);
            }

            // This is the part where individual replicates are combined
            for (cmpd, data) in map {
                let agg = data.into_iter().collect::<Aggregate>();
                v.push((agg, cmpd));
            }
        } else {
            let intermediate = scan_directory(path.as_ref(), fasta)?
                .into_iter()
                .map(|g| g.filter(&ind_filters))
                .collect::<Vec<_>>();

            let mut map = HashMap::new();
            for x in intermediate {
                // println!("{:?}\t{:?}", x.path, count_quantified_lysines(&x));

                let details = extract_meta_data(&x.path).unwrap();
                map.entry(format!("{}_{}", details.cmpd, details.conc))
                    .or_insert_with(Vec::new)
                    .push(x);
            }

            // This is the part where individual replicates are combined
            for (cmpd, data) in map {
                let agg = data.into_iter().collect::<Aggregate>();
                v.push((agg, cmpd));
            }
            break;
        }
    }

    v.sort_by(|a, b| a.1.cmp(&b.1));

    let mut a = Vec::new();
    for (agg, path) in v {
        a.push(
            agg.condense(path.clone(), false).filter(
                &Filter::new()
                    .add_peptide_filter(PeptideFilter::SecondPassFilter(
                        agg_filters,
                        chemotype_map,
                        path,
                    ))
                    .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::CV(0.6))),
            ),
        );
    }
    let a = a.into_iter().collect::<Aggregate>();

    Ok(PeptideCollection::new(&a))
}

fn aggregate_compounds2<P: AsRef<Path>>(
    path: P,
    ind_filters: &Filter,
    agg_filters: &PeptideCollection,
    fasta: &Fasta,
    chemotype_map: &HashMap<String, String>,
) -> io::Result<PeptideCollection> {
    let mut v = Vec::new();
    for entry in fs::read_dir(path.as_ref().clone())? {
        let path2 = entry?.path();
        if path2.is_dir() {
            dbg!(&path2.display());
            let intermediate = scan_directory(&path2, fasta)?
                .into_iter()
                .map(|g| g.filter(&ind_filters))
                .collect::<Vec<_>>();

            let mut map = HashMap::new();
            for x in intermediate {
                let details = extract_meta_data(&x.path).unwrap();
                map.entry(format!("{}_{}", details.cmpd, details.conc))
                    .or_insert_with(Vec::new)
                    .push(x);
            }

            // This is the part where individual replicates are combined
            for (cmpd, data) in map {
                let agg = data.into_iter().collect::<Aggregate>();
                v.push((agg, cmpd));
            }
        } else {
            let intermediate = scan_directory(path.as_ref(), fasta)?
                .into_iter()
                .map(|g| g.filter(&ind_filters))
                .collect::<Vec<_>>();

            let mut map = HashMap::new();
            for x in intermediate {
                // println!("{:?}\t{:?}", x.path, count_quantified_lysines(&x));

                let details = extract_meta_data(&x.path).unwrap();
                map.entry(format!("{}_{}", details.cmpd, details.conc))
                    .or_insert_with(Vec::new)
                    .push(x);
            }

            // This is the part where individual replicates are combined
            for (cmpd, data) in map {
                let agg = data.into_iter().collect::<Aggregate>();
                v.push((agg, cmpd));
            }
            break;
        }
    }

    v.sort_by(|a, b| a.1.cmp(&b.1));

    let mut a = Vec::new();
    for (agg, path) in v {
        a.push(
            agg.condense(path.clone(), false).filter(
                &Filter::new()
                    .add_peptide_filter(PeptideFilter::SecondPassFilter(
                        agg_filters,
                        chemotype_map,
                        path,
                    ))
                    .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::CV(0.6)))
                    .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::Count(2))),
            ),
        );
    }
    let a = a.into_iter().collect::<Aggregate>();

    Ok(PeptideCollection::new(&a))
}

/// Search the parent directory, calling `scan_directory` on all subdirectories
///
/// Each subdirectory (chemotype) is fltered using `filters`, aggregated, and
/// saved into it's own file for further chemotype-specific analysis. The
/// aggregated datasets are then condensed back into a `Grouped` dataset,
/// as if each folder was one experiment.
fn aggregate_scouts<P: AsRef<Path>>(
    path: P,
    ind_filters: &Filter,
    agg_filters: &PeptideCollection,
    fasta: &Fasta,
    chemotype_map: &HashMap<String, String>,
) -> io::Result<PeptideCollection> {
    let mut v = Vec::new();
    for entry in fs::read_dir(path.as_ref().clone())? {
        let path2 = entry?.path();
        dbg!(&path2.display());
        let intermediate = scan_directory(&path2, fasta)?
            .into_iter()
            .map(|g| g.filter(&ind_filters))
            .collect::<Vec<_>>();

        let mut map = HashMap::new();
        for x in intermediate {
            // println!("{:?}\t{:?}", x.path, count_quantified_lysines(&x));

            let details = extract_meta_data(&x.path).unwrap();
            map.entry(format!("{}_{}", details.cmpd, details.conc))
                .or_insert_with(Vec::new)
                .push(x);
        }

        // This is the part where individual replicates are combined
        for (cmpd, data) in map {
            let agg = data.into_iter().collect::<Aggregate>();
            let name = format!("{}_{}", path2.file_name().unwrap().to_str().unwrap(), cmpd);
            v.push((agg, cmpd, name));
        }
    }

    let mut a = Vec::new();
    for (agg, path, name) in v {
        a.push(
            agg.condense(name.clone(), true).filter(
                &Filter::new()
                    .add_peptide_filter(PeptideFilter::SecondPassFilter(
                        agg_filters,
                        chemotype_map,
                        path,
                    ))
                    .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::CV(0.6))),
            ),
        );
    }
    let a = a.into_iter().collect::<Aggregate>();

    Ok(PeptideCollection::new(&a))
}

/// For use with the scout fragments. Instead of aggregating
/// by compound, we just output all of the replicates as single
/// columns
fn replicates_only<P: AsRef<Path>>(
    path: P,
    ind_filters: &Filter,
    agg_filters: &PeptideCollection,
    scope: FilterScope,
    fasta: &Fasta,
    chemotype_map: &HashMap<String, String>,
) -> io::Result<PeptideCollection> {
    let mut v = Vec::new();
    for entry in fs::read_dir(path.as_ref().clone())? {
        let path = entry?.path();
        if path.is_dir() {
            let intermediate = scan_directory(&path, fasta)?
                .into_iter()
                .map(|g| {
                    // Only execute this part for the funky replicate filtering
                    // for scout fragments
                    let p = g.path.clone();
                    let details = extract_meta_data(&p).unwrap();
                    let filters =
                        ind_filters
                            .clone()
                            .add_peptide_filter(PeptideFilter::SecondPassFilter(
                                agg_filters,
                                chemotype_map,
                                format!("{}_{}", details.cmpd, details.conc),
                            ));
                    g.filter(&filters)
                })
                .collect::<Vec<_>>();
            v.extend(intermediate)
        }
    }

    let a = v.into_iter().collect::<Aggregate>();
    dbg!(&a.experiments);
    Ok(PeptideCollection::new(&a))
}

fn write_all_chemotypes(
    pc: &PeptideCollection,
    chemotype_map: &HashMap<String, String>,
) -> io::Result<()> {
    for chemotype in chemotype_map.values().collect::<HashSet<_>>() {
        pc.extract_chemotype(chemotype, chemotype_map)
            .write_peptides(format!("./data/processed/{}.tsv", chemotype))?;
    }
    Ok(())
}

fn main() -> io::Result<()> {
    let ind_filters = Filter::new()
        .add_protein_filter(ProteinFilter::Reverse)
        .add_protein_filter(ProteinFilter::ExcludeMatch("KRT".into()))
        .add_protein_filter(ProteinFilter::ExcludeMatch("Keratin".into()))
        .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::Spurious))
        .add_peptide_filter(PeptideFilter::Ratios(RatioFilter::CV(0.6)))
        .add_peptide_filter(PeptideFilter::Ms2(2));

    let fasta = Fasta::open("./data/db.fasta")?;
    // println!("funnel");

    // {
    //     let chemotype_map = generate_chemotype_map("./data/abbasov_all_insitu")?;
    //     let first_pass = generate_funnel("./data/abbasov_all_insitu", &ind_filters, &fasta)?;
    //     let compounds = aggregate_compounds(
    //         "./data/Chemotypes_insitu_231",
    //         &ind_filters,
    //         &first_pass,
    //         FilterScope::Dataset,
    //         &fasta,
    //         &chemotype_map,
    //     )?;
    //     compounds
    //         .write_peptides("./data/processed/Chemotypes_insitu_231peptides_compound_filt.txt")?;
    //     compounds
    //         .write_proteins("./data/processed/Chemotypes_insitu_231proteins_compound_filt.txt")?;
    // }

    let chemotype_map = generate_chemotype_map("./data/abbasov_noscouts")?;
    let first_pass = generate_funnel("./data/abbasov_noscouts", &ind_filters, &fasta)?;

    let compounds = aggregate_compounds(
        "./data/abbasov_noscouts",
        &ind_filters,
        &first_pass,
        &fasta,
        &chemotype_map,
    )?;
    compounds.write_peptides("peptides_compounds.tsv")?;
    compounds.write_proteins("proteins_compounds.tsv")?;
    let chemo = compounds.combine_into_chemotypes(chemotype_map.clone());
    chemo.write_peptides("peptides_chemotypes.tsv")?;
    chemo.write_proteins("proteins_chemotypes.tsv")?;

    write_all_chemotypes(&compounds, &chemotype_map)?;

    let cas = aggregate_compounds2(
        "./data/abbasov_noscouts/CAS1_IFIT5",
        &ind_filters,
        &first_pass,
        &fasta,
        &chemotype_map,
    )?;
    cas.write_peptides("cas1.tsv")?;

    let first_pass = generate_funnel("./data/abbasov_all", &ind_filters, &fasta)?;
    let chemotype_map = generate_chemotype_map("./data/scouts")?;
    let s2 = aggregate_scouts(
        "./data/scouts",
        &ind_filters,
        &first_pass,
        &fasta,
        &chemotype_map,
    )?;
    s2.write_peptides("scouts.tsv")?;

    Ok(())
}
