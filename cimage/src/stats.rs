#![allow(dead_code)]

/// Calculate the standard deviation (population) of a slice
#[inline]
pub fn stddev(slice: &[f32]) -> f32 {
    let (sum, len) = slice.iter().fold((0.0, 0), |acc, x| (acc.0 + x, acc.1 + 1));
    let n = len as f32;
    let mean = sum / n;

    (slice.iter().fold(0.0f32, |acc, x| acc + (x - mean).powi(2)) / n).sqrt()
}

/// Calculate the coefficient of variation (population) of a slice
#[inline]
pub fn coefficient_variation(slice: &[f32]) -> f32 {
    let (sum, len) = slice.iter().fold((0.0, 0), |acc, x| (acc.0 + x, acc.1 + 1));
    let n = len as f32;
    let mean = sum / n;

    (slice.iter().fold(0.0f32, |acc, x| acc + (x - mean).powi(2)) / n).sqrt() / mean
}
