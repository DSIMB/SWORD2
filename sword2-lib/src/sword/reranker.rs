//! Logistic re-ranker for domain-partitioning candidate selection (Phase B.3).
//!
//! Trained weights come from `benchmark/data/reranker_weights.json`.
//! At inference time, features are chain-local z-score normalised before scoring.

use crate::sword::compute_measure::MeasureLine;

/// 9-feature vector for one domain-partitioning candidate.
pub struct CandidateFeatures {
    pub num_domains: f64,
    pub min_size: f64,
    pub max_cr: f64,
    pub density_min: f64,
    pub mean_density: f64,
    pub n_discontinuous: f64,
    pub size_balance: f64,
    pub largest_domain_frac: f64,
    pub mean_junction_support: f64,
}

impl CandidateFeatures {
    fn as_arr(&self) -> [f64; 9] {
        [
            self.num_domains,
            self.min_size,
            self.max_cr,
            self.density_min,
            self.mean_density,
            self.n_discontinuous,
            self.size_balance,
            self.largest_domain_frac,
            self.mean_junction_support,
        ]
    }
}

const WEIGHTS: [f64; 9] = [
    -1.3329183384825485,   // num_domains
    -0.35472829959224755,  // min_size
    -1.0920904735895476,   // max_cr
     0.11649184545463633,  // density_min
    -0.026792688866579928, // mean_density
    -0.9183576148068141,   // n_discontinuous
     0.6269091008096437,   // size_balance
    -0.4374425483499017,   // largest_domain_frac
     0.020104180409427255, // mean_junction_support
];
const BIAS: f64 = -8.645465737029765;

/// Extract the 9 re-ranker features from a MeasureLine and junction-support map.
pub fn extract_features(
    ml: &MeasureLine,
    jsup_map: &std::collections::HashMap<i32, f64>,
) -> CandidateFeatures {
    let dom_tokens: Vec<&str> = ml.delineation.trim().split_whitespace().collect();
    let n_discontinuous = dom_tokens.iter().filter(|d| d.contains(';')).count() as f64;

    let sizes: Vec<f64> = dom_tokens.iter().map(|d| {
        d.split(';').map(|seg| {
            let parts: Vec<&str> = seg.split('-').collect();
            if parts.len() == 2 {
                parts[1].parse::<f64>().unwrap_or(0.0)
                    - parts[0].parse::<f64>().unwrap_or(0.0) + 1.0
            } else { 0.0 }
        }).sum::<f64>()
    }).collect();
    let total_res: f64 = sizes.iter().sum();
    let n = sizes.len() as f64;
    let mean_sz = if n > 0.0 { total_res / n } else { 1.0 };
    let min_sz = sizes.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_sz = sizes.iter().cloned().fold(0.0_f64, f64::max);
    let size_balance = if mean_sz > 0.0 { (min_sz / mean_sz).min(1.0) } else { 0.0 };
    let largest_domain_frac = if total_res > 0.0 { max_sz / total_res } else { 0.0 };

    let junctions: Vec<i32> = dom_tokens.iter().flat_map(|d| {
        let mut v = Vec::new();
        if let Some(p) = d.find('-') {
            if let Ok(j) = d[..p].parse::<i32>() { v.push(j); }
        }
        if let Some(p) = d.rfind('-') {
            if let Ok(j) = d[p + 1..].parse::<i32>() { v.push(j); }
        }
        v
    }).collect();
    let mean_junction_support = if junctions.is_empty() {
        0.0
    } else {
        let s: f64 = junctions.iter()
            .map(|j| jsup_map.get(j).copied().unwrap_or(0.0))
            .sum();
        s / junctions.len() as f64
    };

    CandidateFeatures {
        num_domains: ml.num_domains as f64,
        min_size: ml.min_size as f64,
        max_cr: ml.max_cr,
        density_min: ml.density_min,
        mean_density: ml.mean_density,
        n_discontinuous,
        size_balance,
        largest_domain_frac,
        mean_junction_support,
    }
}

/// Select the best candidate by logistic score with chain-local z-score normalization.
/// Returns the index of the highest-scoring candidate in `candidates`.
/// Returns 0 if `candidates` is empty.
pub fn rerank(candidates: &[CandidateFeatures]) -> usize {
    if candidates.is_empty() {
        return 0;
    }
    if candidates.len() == 1 {
        return 0;
    }

    let n = candidates.len();
    let nf = n as f64;

    // Compute chain-level mean and std per feature
    let mut means = [0.0f64; 9];
    for c in candidates {
        let arr = c.as_arr();
        for (i, &x) in arr.iter().enumerate() {
            means[i] += x;
        }
    }
    for m in &mut means {
        *m /= nf;
    }

    let mut stds = [0.0f64; 9];
    for c in candidates {
        let arr = c.as_arr();
        for (i, &x) in arr.iter().enumerate() {
            let d = x - means[i];
            stds[i] += d * d;
        }
    }
    for s in &mut stds {
        *s = (*s / nf).sqrt() + 1e-8;
    }

    // Score each candidate: linear(z-scored features) — pick argmax
    let mut best_idx = 0;
    let mut best_score = f64::NEG_INFINITY;
    for (i, c) in candidates.iter().enumerate() {
        let arr = c.as_arr();
        let score: f64 = arr.iter()
            .zip(means.iter())
            .zip(stds.iter())
            .zip(WEIGHTS.iter())
            .map(|(((x, m), s), w)| w * (x - m) / s)
            .sum::<f64>()
            + BIAS;
        if score > best_score {
            best_score = score;
            best_idx = i;
        }
    }
    best_idx
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_candidate(
        nd: f64, ms: f64, cr: f64, dm: f64, md: f64,
        nd2: f64, sb: f64, ldf: f64, mjs: f64,
    ) -> CandidateFeatures {
        CandidateFeatures {
            num_domains: nd,
            min_size: ms,
            max_cr: cr,
            density_min: dm,
            mean_density: md,
            n_discontinuous: nd2,
            size_balance: sb,
            largest_domain_frac: ldf,
            mean_junction_support: mjs,
        }
    }

    #[test]
    fn rerank_empty() {
        assert_eq!(rerank(&[]), 0);
    }

    #[test]
    fn rerank_single() {
        let features = make_candidate(3.0, 50.0, 0.5, 2.0, 2.5, 0.0, 0.8, 0.4, 0.6);
        assert_eq!(rerank(&[features]), 0);
    }

    #[test]
    fn rerank_prefers_fewer_domains() {
        // With negative weight on num_domains, candidate with fewer domains wins
        // (assuming all other features are identical)
        let c1 = make_candidate(5.0, 30.0, 0.7, 2.0, 2.5, 0.0, 0.8, 0.4, 0.5);
        let c2 = make_candidate(2.0, 30.0, 0.7, 2.0, 2.5, 0.0, 0.8, 0.4, 0.5);
        let c3 = make_candidate(8.0, 30.0, 0.7, 2.0, 2.5, 0.0, 0.8, 0.4, 0.5);
        // c2 has fewest domains, should win
        let result = rerank(&[c1, c2, c3]);
        assert_eq!(result, 1);  // c2 is index 1
    }
}
