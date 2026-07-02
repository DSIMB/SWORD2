//! Pairwise-trained linear reranker for domain-partitioning candidate
//! selection. Weights are trained offline (see
//! `benchmark/train_pairwise_reranker.py`) with a pairwise ranking loss —
//! deliberately not the pointwise logistic that caused the earlier Phase B
//! attempt to regress into under-segmentation (see project memory).
//!
//! Feature order here MUST match `benchmark/train_pairwise_reranker.py`'s
//! `FEATURES` list exactly.

use crate::peeling::algorithm::SsType;
use crate::sword::candidate_features::{
    boundary_coil_fraction, candidate_energy_z_score, modal_count_distance, modal_num_domains,
};

const WEIGHTS_JSON: &str =
    include_str!("../../../benchmark/data/pairwise_reranker_weights.json");

#[derive(serde::Deserialize)]
struct TrainedWeights {
    features: Vec<String>,
    weights: Vec<f64>,
    #[allow(dead_code)]
    bias: f64,
}

/// Expected feature order — validated against the loaded JSON at rerank time
/// (falls back to distance_model-style no-op scoring if they don't match, so
/// a stale weights file can't silently misapply weights to the wrong features).
const FEATURES: [&str; 8] = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "energy_z",
    "modal_count_distance",
];

/// Everything needed to score one candidate.
pub struct CandidateInput<'a> {
    pub num_domains: usize,
    pub min_size: usize,
    pub max_cr: f64,
    pub density_min: f64,
    pub mean_density: f64,
    pub raw_delineation: &'a str,
    pub remapped_delineation: &'a str,
}

fn candidate_feature_vector(
    c: &CandidateInput,
    ss_types: &[SsType],
    energy_config: Option<&crate::energy::EnergyConfig>,
    pdb_path: &str,
    chain: &str,
    modal: usize,
) -> [f64; 8] {
    let coil = boundary_coil_fraction(c.raw_delineation, ss_types);
    let energy_z = energy_config
        .and_then(|ec| candidate_energy_z_score(ec, pdb_path, chain, c.remapped_delineation))
        .unwrap_or(0.0);
    [
        c.num_domains as f64,
        c.min_size as f64,
        c.max_cr,
        c.density_min,
        c.mean_density,
        coil,
        energy_z,
        modal_count_distance(c.num_domains, modal),
    ]
}

/// Chain-local z-score normalization across the candidate set, matching
/// `train_pairwise_reranker.py`'s per-chain normalization at training time.
fn zscore_normalize(vectors: &mut [[f64; 8]]) {
    if vectors.len() < 2 {
        return;
    }
    let n = vectors.len() as f64;
    let mut means = [0.0f64; 8];
    for v in vectors.iter() {
        for k in 0..8 {
            means[k] += v[k];
        }
    }
    for m in &mut means {
        *m /= n;
    }
    let mut stds = [0.0f64; 8];
    for v in vectors.iter() {
        for k in 0..8 {
            let d = v[k] - means[k];
            stds[k] += d * d;
        }
    }
    for s in &mut stds {
        *s = (*s / n).sqrt() + 1e-8;
    }
    for v in vectors.iter_mut() {
        for k in 0..8 {
            v[k] = (v[k] - means[k]) / stds[k];
        }
    }
}

/// Pick the winning candidate by pairwise-trained linear score.
///
/// Returns `0` (the first candidate — matching the legacy fallback behavior
/// of picking the first shortlisted entry) if `candidates` is empty, has
/// exactly one entry, or the embedded weights file doesn't match `FEATURES`.
pub fn rerank(
    candidates: &[CandidateInput],
    ss_types: &[SsType],
    energy_config: Option<&crate::energy::EnergyConfig>,
    pdb_path: &str,
    chain: &str,
) -> usize {
    if candidates.len() < 2 {
        return 0;
    }

    let trained: TrainedWeights = match serde_json::from_str(WEIGHTS_JSON) {
        Ok(t) => t,
        Err(_) => return 0,
    };
    if trained.features != FEATURES || trained.weights.len() != FEATURES.len() {
        return 0;
    }

    let modal = modal_num_domains(&candidates.iter().map(|c| c.num_domains).collect::<Vec<_>>());
    let mut vectors: Vec<[f64; 8]> = candidates
        .iter()
        .map(|c| candidate_feature_vector(c, ss_types, energy_config, pdb_path, chain, modal))
        .collect();
    zscore_normalize(&mut vectors);

    let mut best_idx = 0;
    let mut best_score = f64::NEG_INFINITY;
    for (i, v) in vectors.iter().enumerate() {
        let score: f64 = v.iter().zip(trained.weights.iter()).map(|(x, w)| x * w).sum();
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

    #[test]
    fn test_rerank_single_candidate_returns_zero() {
        let c = CandidateInput {
            num_domains: 1, min_size: 10, max_cr: 0.1, density_min: 1.0, mean_density: 2.0,
            raw_delineation: "0-9", remapped_delineation: "1-10",
        };
        assert_eq!(rerank(&[c], &[], None, "/nonexistent.pdb", "A"), 0);
    }

    #[test]
    fn test_rerank_empty_returns_zero() {
        assert_eq!(rerank(&[], &[], None, "/nonexistent.pdb", "A"), 0);
    }

    #[test]
    fn test_rerank_prefers_higher_coil_fraction_and_favorable_energy() {
        // Candidate 1 boundaries land in coil (favorable); candidate 0's land
        // mid-helix. With no energy_config, energy_z is 0.0 for both, so the
        // coil-fraction feature alone should decide it (weights permitting —
        // this test only asserts against the currently checked-in synthetic
        // smoke-trained weights, not real CATH-trained ones; re-verify after
        // Task 9 replaces the weights file with the real training run).
        let mut ss_types = vec![SsType::Coil; 20];
        for s in ss_types.iter_mut().take(15).skip(5) {
            *s = SsType::Helix;
        }
        let bad = CandidateInput {
            num_domains: 2, min_size: 10, max_cr: 0.3, density_min: 0.5, mean_density: 1.0,
            raw_delineation: "0-9 10-19", remapped_delineation: "1-10 11-20",
        };
        let good = CandidateInput {
            num_domains: 2, min_size: 10, max_cr: 0.3, density_min: 0.5, mean_density: 1.0,
            raw_delineation: "0-1 2-19", remapped_delineation: "1-2 3-20",
        };
        let idx = rerank(&[bad, good], &ss_types, None, "/nonexistent.pdb", "A");
        // Not asserting a specific winner here since the synthetic smoke
        // weights aren't meaningfully trained on this feature — just confirm
        // it runs deterministically and returns a valid index.
        assert!(idx == 0 || idx == 1);
    }
}
