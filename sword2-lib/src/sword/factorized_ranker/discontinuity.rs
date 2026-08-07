use super::boundary::{dssp_inverse, retained_bridge_edges, validate_partition_shape};
use super::partition::{FeatureError, ParsedPartition, Segment};
use super::schema::{FeatureMask, DISCONTINUITY_FEATURE_NAMES};
use super::StructuralContext;

const NONLOCAL_SEQUENCE_CUTOFF: usize = 8;
const MEASURE_COUNT: usize = 5;

#[derive(Debug, Clone, Copy)]
struct SegmentEvidence {
    values: [f64; MEASURE_COUNT],
}

pub(crate) fn extract_discontinuity_features(
    partition: &ParsedPartition,
    context: &StructuralContext<'_>,
) -> Result<Vec<f64>, FeatureError> {
    context.validate(FeatureMask {
        global_count: false,
        domain_conditioned: false,
        boundary_local: false,
        relative_hierarchy: false,
        discontinuity: true,
    })?;
    validate_partition_shape(partition, context.ca_coords.len())?;
    let inverse = dssp_inverse(context)?;
    let bridge_edges = retained_bridge_edges(context, &inverse)?;

    let discontinuous_domains: Vec<usize> = partition
        .domains
        .iter()
        .enumerate()
        .filter_map(|(index, domain)| (domain.segments.len() >= 2).then_some(index))
        .collect();
    if discontinuous_domains.is_empty() {
        return Ok(vec![0.0; DISCONTINUITY_FEATURE_NAMES.len()]);
    }

    let mut segment_of = vec![usize::MAX; context.ca_coords.len()];
    let mut segments = Vec::new();
    for (domain_index, domain) in partition.domains.iter().enumerate() {
        for &segment in &domain.segments {
            let segment_index = segments.len();
            for owner in &mut segment_of[segment.start..=segment.end] {
                *owner = segment_index;
            }
            segments.push((domain_index, segment_index, segment));
        }
    }

    let sheet_edges: Vec<(usize, usize)> = bridge_edges
        .into_iter()
        .filter(|&(left, right)| {
            partition.residue_to_domain[left] == partition.residue_to_domain[right]
                && segment_of[left] != segment_of[right]
        })
        .collect();
    let measurements: Vec<SegmentEvidence> = segments
        .iter()
        .filter(|(domain_index, _, _)| discontinuous_domains.contains(domain_index))
        .map(|&(domain_index, segment_index, segment)| {
            segment_evidence(
                partition,
                context,
                &segment_of,
                &sheet_edges,
                domain_index,
                segment_index,
                segment,
            )
        })
        .collect();

    let mut values = Vec::with_capacity(DISCONTINUITY_FEATURE_NAMES.len());
    values.push(1.0);
    for measure in 0..MEASURE_COUNT {
        let measure_values: Vec<f64> = measurements
            .iter()
            .map(|evidence| evidence.values[measure])
            .collect();
        values.extend(summary(&measure_values));
    }
    validate_final(&values)?;
    Ok(values)
}

fn segment_evidence(
    partition: &ParsedPartition,
    context: &StructuralContext<'_>,
    segment_of: &[usize],
    sheet_edges: &[(usize, usize)],
    domain_index: usize,
    segment_index: usize,
    segment: Segment,
) -> SegmentEvidence {
    let current: Vec<usize> = (segment.start..=segment.end).collect();
    let same_domain_other: Vec<usize> = partition.domains[domain_index]
        .residues
        .iter()
        .copied()
        .filter(|&residue| residue < segment.start || residue > segment.end)
        .collect();
    let same_domain_affinity = contact_mean(context, &current, &same_domain_other);
    let competing_affinity = partition
        .domains
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != domain_index)
        .map(|(_, domain)| contact_mean(context, &current, &domain.residues))
        .fold(None, |largest, value| {
            Some(largest.map_or(value, |current: f64| current.max(value)))
        });
    let affinity_margin = competing_affinity
        .map(|largest| same_domain_affinity - largest)
        .unwrap_or(0.0);

    let mut same_domain_long_range_mass = 0.0;
    let mut all_incident_long_range_mass = 0.0;
    let mut bin_masses = [0.0; 4];
    for &left in &current {
        for right in 0..context.ca_coords.len() {
            if (segment.start..=segment.end).contains(&right)
                || left.abs_diff(right) < NONLOCAL_SEQUENCE_CUTOFF
            {
                continue;
            }
            let mass = context.contacts.get(left, right);
            all_incident_long_range_mass += mass;
            if partition.residue_to_domain[right] == domain_index {
                same_domain_long_range_mass += mass;
                bin_masses[separation_bin(left.abs_diff(right))] += mass;
            }
        }
    }
    let long_range_internal_capture = if all_incident_long_range_mass == 0.0 {
        0.0
    } else {
        same_domain_long_range_mass / all_incident_long_range_mass
    };
    let interface_span_entropy = if same_domain_long_range_mass == 0.0 {
        0.0
    } else {
        -bin_masses
            .iter()
            .filter(|&&mass| mass > 0.0)
            .map(|&mass| {
                let probability = mass / same_domain_long_range_mass;
                probability * probability.ln()
            })
            .sum::<f64>()
            / 4.0_f64.ln()
    };
    let same_domain_sheet_links = sheet_edges
        .iter()
        .filter(|&&(left, right)| {
            segment_of[left] == segment_index || segment_of[right] == segment_index
        })
        .count() as f64;

    SegmentEvidence {
        values: [
            same_domain_affinity,
            affinity_margin,
            long_range_internal_capture,
            interface_span_entropy,
            same_domain_sheet_links,
        ],
    }
}

fn contact_mean(context: &StructuralContext<'_>, left: &[usize], right: &[usize]) -> f64 {
    let mut sum = 0.0;
    for &left_index in left {
        for &right_index in right {
            sum += context.contacts.get(left_index, right_index);
        }
    }
    sum / (left.len() * right.len()) as f64
}

fn separation_bin(separation: usize) -> usize {
    match separation {
        0..=15 => 0,
        16..=31 => 1,
        32..=63 => 2,
        _ => 3,
    }
}

fn summary(values: &[f64]) -> [f64; 3] {
    [
        values.iter().copied().fold(f64::INFINITY, f64::min),
        values.iter().sum::<f64>() / values.len() as f64,
        values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    ]
}

fn validate_final(values: &[f64]) -> Result<(), FeatureError> {
    if values.len() != DISCONTINUITY_FEATURE_NAMES.len() {
        return Err(FeatureError::SchemaMismatch);
    }
    if let Some((name, _)) = DISCONTINUITY_FEATURE_NAMES
        .iter()
        .zip(values)
        .find(|(_, value)| !value.is_finite())
    {
        return Err(FeatureError::NonFinite(name));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::sync::OnceLock;

    use crate::dssp::types::BackboneResidue;
    use crate::dssp::DsspChain;
    use crate::peeling::contact_matrix::ContactMatrix;

    use super::extract_discontinuity_features;
    use crate::sword::factorized_ranker::partition::parse_partition;
    use crate::sword::factorized_ranker::schema::DISCONTINUITY_FEATURE_NAMES;
    use crate::sword::factorized_ranker::StructuralContext;

    fn named(values: &[f64], name: &str) -> f64 {
        let index = DISCONTINUITY_FEATURE_NAMES
            .iter()
            .position(|candidate| *candidate == name)
            .unwrap();
        values[index]
    }

    fn dssp_chain(n_residues: usize) -> DsspChain {
        let mut chain = DsspChain::new();
        for _ in 0..n_residues {
            chain.push(BackboneResidue {
                aa: 'A',
                ..BackboneResidue::default()
            });
        }
        chain
    }

    fn context(
        dssp: DsspChain,
        contact_coords: Vec<[f64; 3]>,
        d0: f64,
        delta: f64,
    ) -> StructuralContext<'static> {
        let n_residues = contact_coords.len();
        StructuralContext {
            ca_coords: Box::leak(vec![[0.0, 0.0, 0.0]; n_residues].into_boxed_slice()),
            dssp: Box::leak(Box::new(dssp)),
            contacts: Box::leak(Box::new(ContactMatrix::from_ca_coords(
                &contact_coords,
                d0,
                delta,
            ))),
            iterations: &[],
            measure_provenance: &[],
            dssp_index_for_residue: (1..=n_residues).collect(),
            contact_feature_cache: OnceLock::new(),
        }
    }

    #[test]
    fn continuous_partition_has_explicit_zero_conditional_vector() {
        let values = extract_discontinuity_features(
            &parse_partition("0-3 4-7", 8).unwrap(),
            &context(dssp_chain(8), vec![[0.0, 0.0, 0.0]; 8], 0.0, 1.0),
        )
        .unwrap();
        assert_eq!(values[0], 0.0);
        assert!(values[1..].iter().all(|value| *value == 0.0));
    }

    #[test]
    fn discontinuous_partition_matches_reference_affinity_capture_entropy_and_sheet_links() {
        let partition = parse_partition("0-1;6-7 2-5", 8).unwrap();
        let mut dssp = dssp_chain(8);
        dssp.get_mut(1).partner[0] = 7;
        dssp.get_mut(7).partner[0] = 1;
        let competitor_x = 16.0_f64.ln();
        let contact_coords = vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [competitor_x, 0.0, 0.0],
            [competitor_x, 0.0, 0.0],
            [competitor_x, 0.0, 0.0],
            [competitor_x, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let values = extract_discontinuity_features(
            &partition,
            &context(dssp, contact_coords, 4.0_f64.ln(), 1.0),
        )
        .unwrap();

        assert_eq!(values[0], 1.0);
        for stat in ["min", "mean", "max"] {
            assert!(
                (named(&values, &format!("segment_same_domain_affinity_{stat}")) - 0.8).abs()
                    < 1e-12
            );
            assert!(
                (named(&values, &format!("segment_affinity_margin_{stat}")) - 0.6).abs() < 1e-12
            );
            assert_eq!(
                named(
                    &values,
                    &format!("segment_long_range_internal_capture_{stat}")
                ),
                0.0
            );
            assert_eq!(
                named(&values, &format!("segment_interface_span_entropy_{stat}")),
                0.0
            );
            assert_eq!(
                named(&values, &format!("segment_same_domain_sheet_links_{stat}")),
                1.0
            );
        }
    }

    #[test]
    fn discontinuity_sheet_links_use_bridge_partners_not_equal_labels() {
        let partition = parse_partition("0;3 1-2", 4).unwrap();
        let mut labels_only = dssp_chain(4);
        for index in 1..=4 {
            labels_only.get_mut(index).sheet_label = 'A';
        }
        let label_values = extract_discontinuity_features(
            &partition,
            &context(labels_only, vec![[0.0, 0.0, 0.0]; 4], 0.0, 1.0),
        )
        .unwrap();
        assert_eq!(
            named(&label_values, "segment_same_domain_sheet_links_mean"),
            0.0
        );

        let mut partners_only = dssp_chain(4);
        partners_only.get_mut(1).partner[0] = 4;
        partners_only.get_mut(4).partner[0] = 1;
        let partner_values = extract_discontinuity_features(
            &partition,
            &context(partners_only, vec![[0.0, 0.0, 0.0]; 4], 0.0, 1.0),
        )
        .unwrap();
        assert_eq!(
            named(&partner_values, "segment_same_domain_sheet_links_mean"),
            1.0
        );
    }

    #[test]
    fn long_range_capture_excludes_contacts_within_the_current_segment() {
        let partition = parse_partition("0-9;16-17 10-15", 18).unwrap();
        let mut contact_coords = vec![[0.0, 0.0, 0.0]; 18];
        for coordinate in &mut contact_coords[10..16] {
            coordinate[0] = 1_000.0;
        }
        let values = extract_discontinuity_features(
            &partition,
            &context(dssp_chain(18), contact_coords, 0.0, 1.0),
        )
        .unwrap();

        assert_eq!(
            named(&values, "segment_long_range_internal_capture_min"),
            1.0
        );
        assert_eq!(
            named(&values, "segment_long_range_internal_capture_mean"),
            1.0
        );
        assert_eq!(
            named(&values, "segment_long_range_internal_capture_max"),
            1.0
        );
        assert!(named(&values, "segment_interface_span_entropy_mean") > 0.0);
    }
}
