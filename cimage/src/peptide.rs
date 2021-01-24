use super::*;

#[derive(Clone, Debug, Default, PartialEq)]
/// A container of MS1 ratios, and associated sequence and residue information
pub struct Peptide {
    pub sequence: String,
    pub residue: Residue,
    pub ms2: usize,
    /// A ratio can be set to None if it's filtered out, or if it's
    /// not quantified.
    pub ratios: Vec<Option<f32>>,
}

#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct Identifier {
    pub accession: String,
    pub residue: Residue,
}

impl Peptide {
    /// Turn all `Some(20.0)` into `None`
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// # use cimage::*;
    /// let mut peptide = Peptide {
    ///     sequence: String::default(),
    ///     residue: None,
    ///     ratios: vec![Some(0.), Some(0.), Some(20.), Some(1.2), Some(1.3)],
    /// };
    ///
    /// peptide.remove_twenties();
    ///
    /// assert_eq!(peptide.non_zeroes(), vec![1.2, 1.3]);
    /// ```
    pub fn remove_twenties(&mut self) {
        self.ratios.iter_mut().for_each(|r| {
            if *r == Some(20.0) {
                *r = None
            }
        })
    }

    /// Return a `Vec` of all non-zero, non-None ratios
    ///     
    /// # Example
    ///
    /// ```rust,ignore
    /// # use cimage::*;
    /// let mut peptide = Peptide {
    ///     sequence: String::default(),
    ///     residue: None,
    ///     ratios: vec![Some(0.), Some(0.), Some(20.), Some(1.2), Some(1.3)],
    /// };
    ///
    /// assert_eq!(peptide.non_zeroes(), vec![20.0, 1.2, 1.3]);
    /// ```
    pub fn non_zeroes(&self) -> Vec<f32> {
        self.ratios
            .iter()
            .cloned()
            .filter_map(|r| r)
            .filter(|&r| r != 0.)
            .collect()
    }

    /// Return the median ratio, excluding zeroes and None values
    pub fn median_ratio(&self) -> Option<f32> {
        let mut v = self.non_zeroes();
        v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        if v.len() % 2 == 1 {
            v.get(v.len() / 2).cloned()
        } else {
            let s = v.get(v.len() / 2)? + v.get((v.len() / 2) - 1)?;
            Some(s / 2.0)
        }
    }

    // pub fn max_ratio(&self) -> Option<f32> {
    //     let mut max = None;
    //     for r in self.ratios.iter() {
    //         if let Some(x) = r {
    //             if max.is_none() {
    //                 max = Some(*x);
    //             } else {
    //                 if *x > max.unwrap() {
    //                     max = Some(*x);
    //                 }
    //             }
    //         }
    //     }
    //     max
    // }

    /// Return the number of fully tryptic ends
    pub fn tryptic_ends(&self) -> usize {
        self.sequence.matches("K.").count() + self.sequence.matches("R.").count()
        // + self.sequence.matches('-').count()
    }

    pub fn is_not_half_tryptic(&self) -> bool {
        let cterm = self.sequence.ends_with('-');
        let front = self.sequence.starts_with(|c| match c {
            'K' | 'R' | '-' => true,
            _ => false,
        });
        let end = self
            .sequence
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

    /// Remove ratios equal to 20 if the coefficient of variation for non-zero
    /// ratios is greater than or equal to `cutoff`
    pub fn cv_filter(&mut self, cutoff: f32) {
        let r = self.non_zeroes();
        let cv = stats::coefficient_variation(&r);
        if cv >= cutoff {
            self.remove_twenties();
        }
    }

    /// From Radu's filter_20s method
    /// 20s are stripped if the following conditons are met:
    /// - the set of ratios is not just composed of 0s and 20s
    /// - there is only one 20
    /// - the lowest non-zero non-20 value is below a cutoff
    pub fn spurious_filter(&mut self) {
        let spur = self
            .ratios
            .iter()
            .cloned()
            .filter_map(|r| r)
            .collect::<Spurious>();

        if spur.twenties + spur.zeroes != self.ratios.len() && spur.twenties == 1 && spur.min < 4.0
        {
            self.remove_twenties();
        }
    }
}

#[derive(Copy, Clone, Debug)]
struct Spurious {
    min: f32,
    twenties: usize,
    zeroes: usize,
}

impl FromIterator<f32> for Spurious {
    fn from_iter<I: IntoIterator<Item = f32>>(iter: I) -> Self {
        let mut spur = Spurious {
            min: 20.0,
            twenties: 0,
            zeroes: 0,
        };
        for ratio in iter {
            if ratio == 0.0 {
                spur.zeroes += 1;
            } else if (ratio - 20.0).abs() < 0.001 {
                spur.twenties += 1;
            } else {
                spur.min = ratio.min(spur.min)
            }
        }
        spur
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn peptide(ratios: Vec<f32>) -> Peptide {
        Peptide {
            ratios: ratios.into_iter().map(|r| Some(r)).collect(),
            sequence: String::new(),
            residue: 0,
            ms2: 0,
        }
    }

    fn sequence(s: &str) -> Peptide {
        Peptide {
            ratios: Vec::new(),
            residue: 0,
            ms2: 0,
            sequence: String::from(s),
        }
    }

    #[test]
    fn remove_twenties() {
        let mut peptide = peptide(vec![0., 0., 20., 1.2, 1.3]);
        peptide.remove_twenties();
        assert_eq!(peptide.non_zeroes(), vec![1.2, 1.3])
    }

    #[test]
    fn tryptic_ends() {
        assert_eq!(sequence("R.FGTKGLAITFVSDENDAK.I").tryptic_ends(), 2);
        assert_eq!(sequence("-.MQIFVKTLTGK.T").tryptic_ends(), 1);
    }

    #[test]
    fn non_zeros() {
        let mut peptide = peptide(vec![0., 0., 20., 1.2, 1.3]);
        peptide.ratios.insert(1, None);
        peptide.ratios.insert(4, None);
        assert_eq!(peptide.non_zeroes(), vec![20., 1.2, 1.3])
    }

    #[test]
    fn median_ratio() {
        let mut peptide = peptide(vec![0., 0., 20., 1.2, 1.3]);
        assert_eq!(peptide.median_ratio(), Some(1.3));
        peptide.ratios.remove(4);
        assert_eq!(peptide.median_ratio(), Some(10.6));
    }

    #[test]
    fn cv_filter() {
        let mut peptide = peptide(vec![0., 0., 20., 1.2, 1.3]);

        peptide.cv_filter(0.6);
        assert_eq!(peptide.ratios.iter().position(|&r| r == Some(20.0)), None);
        assert_eq!(
            peptide.ratios,
            vec![Some(0.0), Some(0.0), None, Some(1.2), Some(1.3)]
        )
    }

    #[test]
    fn spurious_filter() {
        let mut pep = peptide(vec![1., 1., 20.]);
        pep.spurious_filter();
        assert_eq!(pep.ratios, vec![Some(1.), Some(1.), None]);

        let mut pep = peptide(vec![4.1, 4.2, 20.]);
        pep.spurious_filter();
        assert_eq!(pep.ratios, vec![Some(4.1), Some(4.2), Some(20.)]);

        let mut pep = peptide(vec![0., 0., 20., 20.]);
        pep.spurious_filter();
        assert_eq!(pep.ratios, vec![Some(0.), Some(0.), Some(20.), Some(20.)]);
    }

    #[test]
    fn spurious_collect() {
        let ratios = vec![19.0, 19.5, 1.2, 1.3, 0.67, 0., 0., 20.0, 8.0];
        let spur = ratios.into_iter().collect::<Spurious>();
        assert!((spur.min - 0.67).abs() < 0.0001);
        assert_eq!(spur.twenties, 1);
        assert_eq!(spur.zeroes, 2);
    }

    #[test]
    fn unwrap() {
        let x = Some(10.);
        assert_eq!(x.unwrap(), 10.0);
        assert_eq!(x.unwrap(), 10.0);
    }
}
