//! Structs and methods for working with raw CIMAGE data with
//! the output_excel_rt format output by cimage
use super::*;
use memchr::{memchr_iter, Memchr};
use std::io::prelude::*;

use uniprot::fasta::Fasta;

/// Generalized wrapper around [`Memchr`] iterator for splitting `&[u8]` slices
/// by a byte.
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

/// Raw cimage data taken directly from a combined_dta file
pub struct Raw {
    pub events: Vec<Event>,
    pub ms2: HashMap<String, usize>,
    pub path: PathBuf,
}

/// Quantified MS1 event
#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct Event {
    acc: String,
    desc: String,
    seq: String,
    residue: Residue,
    ratio: f32,
    r_squared: f32,
}

pub fn is_not_half_tryptic(sequence: &str) -> bool {
    let cterm = sequence.ends_with('-');
    let front = sequence.starts_with(|c| match c {
        'K' | 'R' | '-' => true,
        _ => false,
    });
    let end = sequence
        .split('.')
        .skip(1)
        .next()
        .map(|s| {
            s.ends_with(|c| match c {
                'K' | 'R' => true,
                _ => cterm,
            })
        })
        .unwrap_or(false);
    front && end
}

pub fn liganded_tryptic_end(sequence: &str) -> bool {
    // .add_peptide_filter(PeptideFilter::ExcludeMatch("K.K*"))
    // .add_peptide_filter(PeptideFilter::ExcludeMatch("R.K*"))
    // .add_peptide_filter(PeptideFilter::ExcludeMatch("K*.K"))
    // .add_peptide_filter(PeptideFilter::ExcludeMatch("K*.R"))
    sequence.contains("K.K*")
        || sequence.contains("R.K*")
        || sequence.contains("K*.K")
        || sequence.contains("K*.R")
}

#[inline]
/// parse a line of the output_rt_to_excel output containing information
/// relevant to us
unsafe fn read_event(buf: &[u8], fasta: &Fasta) -> Option<Event> {
    // Skip index field
    let mut pitchfork = Pitchfork::new('\t' as u8, buf).skip(1);
    let acc = String::from_utf8_unchecked(pitchfork.next()?.to_owned());
    let desc = String::from_utf8_unchecked(pitchfork.next()?.to_owned());
    // skip symbol
    let _ = pitchfork.next();
    let seq = String::from_utf8_unchecked(pitchfork.next()?.to_owned());
    // skip mass, charge, segment
    let _ = pitchfork.next()?;
    let _ = pitchfork.next()?;
    let _ = pitchfork.next()?;
    let ratio = std::str::from_utf8_unchecked(pitchfork.next()?)
        .parse::<f32>()
        .ok()?;

    // skip int, np fields
    let _ = pitchfork.next()?;
    let _ = pitchfork.next()?;
    let r_squared = std::str::from_utf8_unchecked(pitchfork.next()?)
        .parse::<f32>()
        .ok()?;

    if !is_not_half_tryptic(&seq) || liganded_tryptic_end(&seq) {
        return None;
    }

    let normalized_seq = seq.split('.').skip(1).next()?;
    let residue = fasta.assign(&acc, &normalized_seq)? as Residue;

    Some(Event {
        acc,
        desc,
        seq,
        residue,
        ratio,
        r_squared,
    })
}

impl Raw {
    /// Parse output_rt_to_excel.txt and DTASelect output files to generate spectral counts
    /// and ratiometric information for peptide events
    pub fn load<P: AsRef<Path>>(output_rt: P, dtaselect: P, fasta: &Fasta) -> io::Result<Raw> {
        let mut ms2 = HashMap::new();
        let mut buffer = Vec::new();
        let mut file = fs::File::open(output_rt.as_ref())?;
        file.read_to_end(&mut buffer)?;

        let events = Pitchfork::new('\n' as u8, &buffer)
            .filter_map(|buf| unsafe { read_event(buf, fasta) })
            .filter(|ev| ev.r_squared >= 0.80)
            .collect::<Vec<Event>>();

        // Parse DTASelect files to extract spectral count data
        buffer.clear();
        file = fs::File::open(dtaselect.as_ref())?;
        file.read_to_end(&mut buffer)?;
        for line in Pitchfork::new('\n' as u8, &buffer).skip(28) {
            if line.is_empty() {
                continue;
            }
            if line[0] == '\t' as u8 || line[0] == '*' as u8 {
                let fields = String::from_utf8(line.to_vec()).unwrap();
                let fields = fields.split('\t').collect::<Vec<_>>();
                if fields.len() >= 15 {
                    let mut seq: String = fields[14].into();
                    seq = seq.replace("(464.24957)", "*");
                    seq = seq.replace("(470.26338)", "*");
                    *ms2.entry(seq.into()).or_insert(0) += 1;
                }
            }
        }

        Ok(Raw {
            events,
            ms2,
            path: PathBuf::from(output_rt.as_ref()),
        })
    }

    /// Group the dataset by protein
    pub fn group(self) -> Grouped {
        let mut table = HashMap::new();
        let mut used: HashSet<(String, String)> = HashSet::new();

        for ev in self.events {
            // only add ms2 once.
            let ms2 = if !used.contains(&(ev.acc.clone(), ev.seq.clone())) {
                used.insert((ev.acc.clone(), ev.seq.clone()));
                self.ms2.get(&ev.seq).copied().unwrap_or(0)
            } else {
                0
            };

            table
                .entry(ev.acc.clone())
                .or_insert_with(|| Protein::new(ev.acc.clone(), ev.desc.clone()))
                .add_ratio(ev.residue, &ev.seq, Some(ev.ratio), ms2);
        }

        Grouped {
            proteins: table,
            path: self.path.file_name().unwrap().to_str().unwrap().to_string(),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn pitchfork() {
        let input = "hello\tworld\tfield\t1111";
        let mut pitch = Pitchfork::new('\t' as u8, input.as_bytes());
        assert_eq!(pitch.next().unwrap(), "hello".as_bytes());
        assert_eq!(pitch.next().unwrap(), "world".as_bytes());
        assert_eq!(pitch.next().unwrap(), "field".as_bytes());
        assert_eq!(pitch.next().unwrap(), "1111".as_bytes());
        assert_eq!(pitch.next(), None);
    }

    #[test]
    /// Basic functionality test for grouping by protein
    fn group() {
        macro_rules! ev {
            ($acc:expr, $seq:expr, $site:expr, $r:expr) => {
                Event {
                    acc: $acc.into(),
                    seq: $seq.into(),
                    residue: $site,
                    ratio: $r,
                    desc: String::new(),
                    r_squared: 1.0,
                }
            };
        }

        macro_rules! pep {
            ($seq: expr, $site:expr, $($r:expr),+) => {
                Peptide { sequence: $seq.into(), ms2: 0, residue: $site, ratios: vec![$(Some($r)),+]}
            };
        }

        let events = vec![
            ev!("Q1", "-.MEHK*QLL.K", 0, 10.0),
            ev!("Q1", "-.MEHK*QLL.K", 0, 20.0),
            ev!("Q1", "K.RYK*TQC.K", 1, 2.0),
            ev!("Q1", "K.RYK*TQC.K", 1, 2.5),
            ev!("Q2", "K.ITYAK*AG.K", 2, 2.0),
            ev!("Q2", "K.ITYAK*AG.K", 2, 1.3),
        ];

        let raw = Raw {
            events,
            ms2: HashMap::default(),
            path: PathBuf::from("test.txt"),
        };

        let grouped = raw.group();

        let mut q1 = Protein::from_iter(
            vec![
                pep!("-.MEHK*QLL.K", 0, 10.0, 20.0),
                pep!("K.RYK*TQC.K", 1, 2.0, 2.5),
            ]
            .into_iter(),
        );
        q1.accession = String::from("Q1");

        let mut q2 = Protein::from_iter(vec![pep!("K.ITYAK*AG.K", 2, 2.0, 1.3)].into_iter());
        q2.accession = String::from("Q2");

        let mut table = HashMap::new();
        table.insert(q1.accession.clone(), q1);
        table.insert(q2.accession.clone(), q2);
        assert_eq!(grouped.proteins, table);
    }
}
