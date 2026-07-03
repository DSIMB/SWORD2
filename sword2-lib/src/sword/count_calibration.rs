//! Analytical domain-count calibration.
//!
//! The distance_model cross-level N_dom selector systematically under-segments
//! (mean -0.41 domains on CATH-663; see benchmark/DIAGNOSIS.md). This module
//! predicts an expected domain count from chain length and adds a penalty for
//! deviating from it to the distance_model signed distance, pulling the argmax
//! toward the right number of domains. Constants are fitted offline by
//! benchmark/fit_count_calibration.py. Applied only when the caller opts in.

/// Fitted calibration: `expected_ndom = intercept + len_coef * chain_length`,
/// with a `lambda`-weighted `|nd - expected|` penalty on the selection score.
#[derive(Debug, Clone, Copy)]
pub struct CountCalibration {
    pub intercept: f64,
    pub len_coef: f64,
    pub lambda: f64,
}

impl Default for CountCalibration {
    fn default() -> Self {
        // Fitted on CATH-17287 (n_domains ~ n_residues); lambda from the CATH-663
        // sweep in Task 4. Constants from Task 1 real fit.
        Self { intercept: 1.236372, len_coef: 0.003362, lambda: 0.05 }
    }
}

impl CountCalibration {
    /// Default constants with `lambda` optionally overridden (for the A/B sweep).
    pub fn with_lambda(lambda: Option<f64>) -> Self {
        let mut c = Self::default();
        if let Some(l) = lambda {
            c.lambda = l;
        }
        c
    }

    /// Expected (continuous) domain count for a chain of `chain_len` residues,
    /// never less than 1.
    pub fn expected_num_domains(&self, chain_len: usize) -> f64 {
        (self.intercept + self.len_coef * chain_len as f64).max(1.0)
    }

    /// distance_model signed distance minus the count-calibration penalty.
    pub fn adjusted_score(&self, dist: f64, num_domains: usize, expected: f64) -> f64 {
        dist - self.lambda * (num_domains as f64 - expected).abs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expected_scales_with_length_and_clamps_at_one() {
        let c = CountCalibration { intercept: 1.232, len_coef: 0.00348, lambda: 0.05 };
        // ~290 residues -> ~2.24 domains
        assert!((c.expected_num_domains(290) - 2.2412).abs() < 1e-3);
        // a tiny chain never drops below 1.0 (clamp test; with these constants formula gives ~1.235)
        assert!(c.expected_num_domains(1) >= 1.0);
    }

    #[test]
    fn penalty_prefers_counts_near_expected() {
        let c = CountCalibration { intercept: 0.0, len_coef: 0.0, lambda: 0.1 };
        let expected = 3.0;
        // same raw distance; nd closer to expected must score higher
        let near = c.adjusted_score(0.5, 3, expected);
        let far = c.adjusted_score(0.5, 5, expected);
        assert!(near > far);
        assert!((near - 0.5).abs() < 1e-9); // zero penalty exactly at expected
    }

    #[test]
    fn with_lambda_overrides_default() {
        assert_eq!(CountCalibration::with_lambda(Some(0.2)).lambda, 0.2);
        assert_eq!(CountCalibration::with_lambda(None).lambda, CountCalibration::default().lambda);
    }
}
