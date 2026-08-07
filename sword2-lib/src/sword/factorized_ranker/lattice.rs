use std::collections::{BTreeMap, BTreeSet};

use crate::peeling::algorithm::IterationResult;
use crate::sword::compute_measure::{MeasureLine, MeasureProvenance};
use crate::sword::count_calibration::CountCalibration;
use crate::sword::distance_model;

use super::partition::{parse_partition, FeatureError, ParsedPartition};

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct HierarchyEvidence {
    pub first_appearance_level: usize,
    pub persistence_levels: usize,
    pub parent_merge_margin: f64,
    pub child_merge_margin: f64,
    pub hierarchy_path_count: u64,
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct CandidateRecord {
    pub source_index: usize,
    pub measure: MeasureLine,
    pub partition: ParsedPartition,
    pub legacy_distance: f64,
    /// Filled from aligned `MeasureProvenance` plus Peeling iterations when
    /// complete typed evidence is available.
    pub hierarchy: Option<HierarchyEvidence>,
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct CandidateLattice {
    pub candidates: Vec<CandidateRecord>,
    pub groups: BTreeMap<usize, Vec<usize>>,
}

impl CandidateLattice {
    pub(crate) fn from_first_pass(
        measures: &[MeasureLine],
        shortlisted_indices: &[usize],
        chain_len: usize,
    ) -> Result<Self, FeatureError> {
        let mut candidates = Vec::new();
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut seen = BTreeSet::new();

        for &source_index in shortlisted_indices {
            let measure = measures.get(source_index).ok_or(FeatureError::Malformed)?;
            if groups
                .get(&measure.num_domains)
                .is_some_and(|group| group.len() >= 3)
            {
                continue;
            }
            let partition = parse_partition(&measure.delineation, chain_len)?;
            if partition.domains.len() != measure.num_domains {
                return Err(FeatureError::DomainCountMismatch);
            }
            if !seen.insert((measure.num_domains, partition.canonical.clone())) {
                continue;
            }

            let candidate_index = candidates.len();
            candidates.push(CandidateRecord {
                source_index,
                measure: measure.clone(),
                legacy_distance: rendered_legacy_distance(measure),
                partition,
                hierarchy: None,
            });
            groups
                .entry(measure.num_domains)
                .or_default()
                .push(candidate_index);
        }

        Ok(Self { candidates, groups })
    }

    #[allow(dead_code)]
    pub(crate) fn attach_hierarchy(
        &mut self,
        provenance: &[MeasureProvenance],
        iterations: &[IterationResult],
    ) -> Result<(), FeatureError> {
        if iterations.is_empty() {
            return Err(FeatureError::MissingContext("Peeling iterations"));
        }

        for candidate in &mut self.candidates {
            let provenance = provenance
                .get(candidate.source_index)
                .ok_or(FeatureError::MissingContext("measure provenance"))?;
            let boundaries = partition_boundaries(&candidate.partition);
            let matching_levels: Vec<usize> = iterations
                .iter()
                .enumerate()
                .filter_map(|(level, iteration)| {
                    let pu_ends: BTreeSet<usize> = iteration
                        .pu_boundaries
                        .iter()
                        .map(|boundary| boundary[1])
                        .collect();
                    boundaries
                        .iter()
                        .all(|boundary| pu_ends.contains(boundary))
                        .then_some(level)
                })
                .collect();
            let first_appearance_level =
                matching_levels
                    .first()
                    .copied()
                    .ok_or(FeatureError::MissingContext(
                        "candidate boundary absent from Peeling hierarchy",
                    ))?;
            let persistence_levels = matching_levels
                .iter()
                .filter(|&&level| level >= first_appearance_level)
                .count();

            candidate.hierarchy = Some(HierarchyEvidence {
                first_appearance_level,
                persistence_levels,
                parent_merge_margin: merge_margin(&provenance.incoming_merge_qualities),
                child_merge_margin: merge_margin(&provenance.outgoing_merge_qualities),
                hierarchy_path_count: provenance.hierarchy_path_count,
            });
        }

        Ok(())
    }
}

fn partition_boundaries(partition: &ParsedPartition) -> Vec<usize> {
    let mut segments: Vec<_> = partition
        .domains
        .iter()
        .flat_map(|domain| domain.segments.iter())
        .collect();
    segments.sort_by_key(|segment| segment.start);
    segments
        .iter()
        .take(segments.len().saturating_sub(1))
        .map(|segment| segment.end)
        .collect()
}

fn merge_margin(qualities: &[f64]) -> f64 {
    let mut finite: Vec<f64> = qualities
        .iter()
        .copied()
        .filter(|quality| quality.is_finite())
        .collect();
    finite.sort_by(|left, right| right.total_cmp(left));
    match finite.as_slice() {
        [largest, second_largest, ..] => largest - second_largest,
        _ => 0.0,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct LegacySelection {
    pub num_domains: usize,
    pub measure_line: String,
    pub source_index: Option<usize>,
}

pub(crate) fn select_legacy(
    relevant: &[MeasureLine],
    chain_len: usize,
    use_count_calibration: bool,
    count_lambda: Option<f64>,
) -> LegacySelection {
    let calibration = use_count_calibration.then(|| {
        let calibration = CountCalibration::with_lambda(count_lambda);
        let expected = calibration.expected_num_domains(chain_len);
        (calibration, expected)
    });
    let mut best_num_domains = 0;
    let mut best_score = f64::NEG_INFINITY;
    let mut best_index = None;
    let mut max_domains = None;

    for (index, measure) in relevant
        .iter()
        .take(relevant.len().saturating_sub(1))
        .enumerate()
    {
        if max_domains.is_none() {
            max_domains = Some(measure.num_domains);
        }
        if max_domains != Some(measure.num_domains) {
            continue;
        }
        max_domains = measure.num_domains.checked_sub(1);
        let distance = rendered_legacy_distance(measure);
        let score = match calibration {
            Some((ref calibration, expected)) => {
                calibration.adjusted_score(distance, measure.num_domains, expected)
            }
            None => distance,
        };
        if measure.num_domains > 0 && score > best_score {
            best_score = score;
            best_num_domains = measure.num_domains;
            best_index = Some(index);
        }
    }

    let num_domains = best_num_domains.max(1);
    let measure_line = best_index
        .and_then(|index| relevant.get(index))
        .map(MeasureLine::to_line)
        .unwrap_or_default();
    LegacySelection {
        num_domains,
        measure_line,
        source_index: best_index,
    }
}

fn rendered_legacy_distance(measure: &MeasureLine) -> f64 {
    let rendered = measure.to_line();
    let fields: Vec<&str> = rendered.split('|').collect();
    let max_cr = fields
        .get(3)
        .and_then(|field| field.trim().parse().ok())
        .unwrap_or(0.0);
    let density_min = fields
        .get(5)
        .and_then(|field| field.trim().parse().ok())
        .unwrap_or(0.0);
    distance_model::distance_model(max_cr, density_min, 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::peeling::algorithm::IterationResult;
    use crate::sword::compute_measure::MeasureProvenance;
    use crate::sword::distance_model;

    fn measure(
        num_domains: usize,
        delineation: &str,
        max_cr: f64,
        density_min: f64,
    ) -> MeasureLine {
        MeasureLine {
            num_domains,
            min_size: 1,
            delineation: delineation.to_string(),
            max_cr,
            mean_cr: 0.0,
            density_min,
            mean_density: density_min,
        }
    }

    fn synthetic_measures_with_four_two_domain_candidates() -> Vec<MeasureLine> {
        vec![
            measure(3, "0-3 4-7 8-11", 0.10, 3.0),
            measure(2, "0-5 6-11", 0.11, 3.1),
            measure(2, "0-4 5-11", 0.12, 3.2),
            measure(2, "0-3 4-11", 0.13, 3.3),
            measure(2, "0-2 3-11", 0.14, 3.4),
        ]
    }

    fn synthetic_measures() -> Vec<MeasureLine> {
        synthetic_measures_with_four_two_domain_candidates()
    }

    fn synthetic_measure_strings() -> Vec<String> {
        synthetic_measures()
            .iter()
            .map(MeasureLine::to_line)
            .collect()
    }

    fn identity_indices() -> Vec<usize> {
        (0..synthetic_measures_with_four_two_domain_candidates().len()).collect()
    }

    fn legacy_selection_reference(lines: &[String]) -> (usize, String) {
        let mut best_num_domains = 0;
        let mut best_distance = f64::NEG_INFINITY;
        let mut next_num_domains = None;

        for line in lines.iter().take(lines.len().saturating_sub(1)) {
            let fields: Vec<&str> = line.split('|').collect();
            let num_domains = fields[0].trim().parse::<usize>().unwrap();
            if next_num_domains.is_none() {
                next_num_domains = Some(num_domains);
            }
            if next_num_domains != Some(num_domains) {
                continue;
            }
            next_num_domains = num_domains.checked_sub(1);
            let distance = distance_model::distance_model(
                fields[3].trim().parse().unwrap(),
                fields[5].trim().parse().unwrap(),
                1,
            );
            if distance > best_distance {
                best_distance = distance;
                best_num_domains = num_domains;
            }
        }

        let measure_line = lines
            .iter()
            .find(|line| {
                line.split('|')
                    .next()
                    .and_then(|field| field.trim().parse::<usize>().ok())
                    == Some(best_num_domains)
            })
            .cloned()
            .unwrap_or_default();
        (best_num_domains.max(1), measure_line)
    }

    #[test]
    fn lattice_is_before_legacy_count_filter_and_caps_three_per_count() {
        let measures = synthetic_measures_with_four_two_domain_candidates();
        let lattice =
            CandidateLattice::from_first_pass(&measures, &identity_indices(), 12).unwrap();
        assert_eq!(lattice.groups[&2].len(), 3);
        assert!(lattice.groups.contains_key(&3));
    }

    #[test]
    fn legacy_selection_is_stable() {
        let old = legacy_selection_reference(&synthetic_measure_strings());
        let new = select_legacy(&synthetic_measures(), 120, false, None);
        assert_eq!(new.num_domains, old.0);
        assert_eq!(new.measure_line, old.1);
    }

    #[test]
    fn legacy_selection_uses_rendered_six_decimal_scores() {
        let measures = vec![
            measure(3, "0 1 2", 0.099_999_2, 3.499_999_2),
            measure(2, "0 1-2", 0.099_999_4, 3.499_999_6),
            measure(1, "0-2", 0.0, 0.0),
        ];
        let old = legacy_selection_reference(
            &measures
                .iter()
                .map(MeasureLine::to_line)
                .collect::<Vec<_>>(),
        );
        let new = select_legacy(&measures, 120, false, None);
        assert_eq!(old.0, 2);
        assert_eq!(new.num_domains, old.0);
        assert_eq!(new.measure_line, old.1);
    }

    #[test]
    fn hierarchy_attachment_uses_merge_provenance_and_pu_boundaries() {
        let measures = vec![measure(2, "0-3 4-7", 0.1, 3.0)];
        let mut lattice = CandidateLattice::from_first_pass(&measures, &[0], 8).unwrap();
        let provenance = vec![MeasureProvenance {
            incoming_merge_qualities: vec![0.1, 0.4],
            outgoing_merge_qualities: vec![0.2, 0.8],
            hierarchy_path_count: 4,
            ..MeasureProvenance::default()
        }];
        let iterations = vec![
            IterationResult {
                max_cr: 0.0,
                min_density: 0.0,
                ci: 0.0,
                r: 0.0,
                num_pus: 2,
                pu_boundaries: vec![[0, 3], [4, 7]],
            },
            IterationResult {
                max_cr: 0.0,
                min_density: 0.0,
                ci: 0.0,
                r: 0.0,
                num_pus: 3,
                pu_boundaries: vec![[0, 1], [2, 3], [4, 7]],
            },
        ];

        lattice.attach_hierarchy(&provenance, &iterations).unwrap();

        let evidence = lattice.candidates[0].hierarchy.as_ref().unwrap();
        assert_eq!(evidence.first_appearance_level, 0);
        assert_eq!(evidence.persistence_levels, 2);
        assert!((evidence.parent_merge_margin - 0.3).abs() < 1e-12);
        assert!((evidence.child_merge_margin - 0.6).abs() < 1e-12);
        assert_eq!(evidence.hierarchy_path_count, 4);
    }
}
