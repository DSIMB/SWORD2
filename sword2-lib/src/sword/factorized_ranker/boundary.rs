use std::collections::{BTreeMap, BTreeSet};

use super::features::ContactFeatureCache;
use super::partition::{FeatureError, ParsedPartition, Segment};
use super::schema::{FeatureMask, BOUNDARY_LOCAL_FEATURE_NAMES};
use super::StructuralContext;

const MEASURE_COUNT: usize = 14;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SsClass {
    Helix,
    Strand,
    Coil,
}

#[derive(Debug, Clone, Copy)]
struct BoundaryEvidence {
    values: [f64; MEASURE_COUNT],
}

pub(crate) fn extract_boundary_features(
    partition: &ParsedPartition,
    context: &StructuralContext<'_>,
) -> Result<Vec<f64>, FeatureError> {
    let contact_cache = context
        .contact_feature_cache
        .get_or_init(|| ContactFeatureCache::new(context.contacts));
    context.validate(FeatureMask {
        global_count: false,
        domain_conditioned: false,
        boundary_local: true,
        relative_hierarchy: false,
        discontinuity: false,
    })?;
    validate_partition_shape(partition, context.ca_coords.len())?;

    let inverse = dssp_inverse(context)?;
    let hydrogen_bonds = retained_hydrogen_bonds(context, &inverse)?;
    let bridges = retained_bridge_edges(context, &inverse)?;
    let boundaries: Vec<usize> = partition
        .residue_to_domain
        .windows(2)
        .enumerate()
        .filter_map(|(index, owners)| (owners[0] != owners[1]).then_some(index))
        .collect();
    if boundaries.is_empty() {
        return if partition.domains.len() == 1 {
            Ok(vec![0.0; BOUNDARY_LOCAL_FEATURE_NAMES.len()])
        } else {
            Err(FeatureError::SchemaMismatch)
        };
    }

    let classes: Vec<SsClass> = context
        .dssp_index_for_residue
        .iter()
        .map(|&index| ss_class(context.dssp.get(index).ss[0]))
        .collect();
    let sheet_labels: Vec<char> = context
        .dssp_index_for_residue
        .iter()
        .map(|&index| context.dssp.get(index).sheet_label)
        .collect();
    let measurements: Vec<BoundaryEvidence> = boundaries
        .iter()
        .map(|&boundary| {
            boundary_evidence(
                boundary,
                context,
                contact_cache,
                &classes,
                &sheet_labels,
                &hydrogen_bonds,
                &bridges,
            )
        })
        .collect();

    let mut values = Vec::with_capacity(BOUNDARY_LOCAL_FEATURE_NAMES.len());
    for measure in 0..MEASURE_COUNT {
        let measure_values: Vec<f64> = measurements
            .iter()
            .map(|evidence| evidence.values[measure])
            .collect();
        values.extend(summary(&measure_values));
    }
    validate_final(BOUNDARY_LOCAL_FEATURE_NAMES, &values)?;
    Ok(values)
}

fn boundary_evidence(
    boundary: usize,
    context: &StructuralContext<'_>,
    contact_cache: &ContactFeatureCache,
    classes: &[SsClass],
    sheet_labels: &[char],
    hydrogen_bonds: &BTreeMap<(usize, usize), i64>,
    bridges: &BTreeSet<(usize, usize)>,
) -> BoundaryEvidence {
    let n_residues = context.ca_coords.len();
    let same_class = classes[boundary] == classes[boundary + 1];
    let inside = if same_class {
        classes[boundary]
    } else {
        SsClass::Coil
    };
    let sse_terminus_distance = if same_class && inside != SsClass::Coil {
        let mut start = boundary;
        while start > 0 && classes[start - 1] == inside {
            start -= 1;
        }
        let mut end = boundary + 1;
        while end + 1 < n_residues && classes[end + 1] == inside {
            end += 1;
        }
        (boundary - start + 1).min(end - boundary) as f64
    } else {
        0.0
    };

    let crossing_bonds: Vec<(&(usize, usize), &i64)> = hydrogen_bonds
        .iter()
        .filter(|((left, right), _)| crosses(*left, *right, boundary))
        .collect();
    let hbond_energy_kcal = crossing_bonds
        .iter()
        .map(|(_, energy)| **energy)
        .sum::<i64>() as f64
        / 1_000.0;
    let bridge_count = bridges
        .iter()
        .filter(|&&(left, right)| crosses(left, right, boundary))
        .count();
    let sheet_link_count = same_sheet_pair_count(sheet_labels, boundary);
    let (long_range_sum, long_range_count) = contact_cache.long_rectangle(
        Segment {
            start: 0,
            end: boundary,
        },
        Segment {
            start: boundary + 1,
            end: n_residues - 1,
        },
    );
    let left_residue = context.dssp.get(context.dssp_index_for_residue[boundary]);
    let right_residue = context
        .dssp
        .get(context.dssp_index_for_residue[boundary + 1]);
    let virtual_dihedral_change = if left_residue.alpha == 360.0 || right_residue.alpha == 360.0 {
        0.0
    } else {
        ((left_residue.alpha - right_residue.alpha + 180.0).rem_euclid(360.0) - 180.0).abs() / 180.0
    };

    BoundaryEvidence {
        values: [
            (same_class && inside == SsClass::Helix) as u8 as f64,
            (same_class && inside == SsClass::Strand) as u8 as f64,
            (inside == SsClass::Coil) as u8 as f64,
            sse_terminus_distance,
            crossing_bonds.len() as f64,
            hbond_energy_kcal,
            bridge_count as f64,
            sheet_link_count,
            insulation(context, boundary, 8),
            insulation(context, boundary, 16),
            insulation(context, boundary, 32),
            if long_range_count == 0.0 {
                0.0
            } else {
                long_range_sum / long_range_count
            },
            (left_residue.kappa - right_residue.kappa).abs() / 180.0,
            virtual_dihedral_change,
        ],
    }
}

fn same_sheet_pair_count(sheet_labels: &[char], boundary: usize) -> f64 {
    let mut left_counts = BTreeMap::new();
    let mut right_counts = BTreeMap::new();
    for &label in &sheet_labels[..=boundary] {
        if !label.is_whitespace() {
            *left_counts.entry(label).or_insert(0usize) += 1;
        }
    }
    for &label in &sheet_labels[boundary + 1..] {
        if !label.is_whitespace() {
            *right_counts.entry(label).or_insert(0usize) += 1;
        }
    }
    left_counts
        .into_iter()
        .map(|(label, left)| left as f64 * right_counts.get(&label).copied().unwrap_or(0) as f64)
        .sum()
}

fn insulation(context: &StructuralContext<'_>, boundary: usize, window: usize) -> f64 {
    let left_start = boundary.saturating_sub(window - 1);
    let right_end = boundary
        .saturating_add(window)
        .min(context.ca_coords.len() - 1);
    let count = (boundary - left_start + 1) * (right_end - boundary);
    1.0 - context
        .contacts
        .rectangle_sum(left_start, boundary + 1, boundary, right_end)
        / count as f64
}

fn crosses(left: usize, right: usize, boundary: usize) -> bool {
    left.min(right) <= boundary && boundary < left.max(right)
}

fn ss_class(code: char) -> SsClass {
    match code {
        'H' | 'G' | 'I' => SsClass::Helix,
        'E' | 'B' => SsClass::Strand,
        _ => SsClass::Coil,
    }
}

pub(super) fn validate_partition_shape(
    partition: &ParsedPartition,
    n_residues: usize,
) -> Result<(), FeatureError> {
    if partition.domains.is_empty() || partition.residue_to_domain.len() != n_residues {
        return Err(FeatureError::SchemaMismatch);
    }
    let mut covered = vec![false; n_residues];
    for (domain_index, domain) in partition.domains.iter().enumerate() {
        if domain.segments.is_empty() || domain.residues.is_empty() {
            return Err(FeatureError::SchemaMismatch);
        }
        let mut residues = Vec::new();
        for segment in &domain.segments {
            if segment.start > segment.end || segment.end >= n_residues {
                return Err(FeatureError::SchemaMismatch);
            }
            for residue in segment.start..=segment.end {
                if covered[residue] || partition.residue_to_domain[residue] != domain_index {
                    return Err(FeatureError::SchemaMismatch);
                }
                covered[residue] = true;
                residues.push(residue);
            }
        }
        if residues != domain.residues {
            return Err(FeatureError::SchemaMismatch);
        }
    }
    if covered.iter().any(|covered| !covered)
        || partition
            .residue_to_domain
            .iter()
            .any(|&owner| owner >= partition.domains.len())
    {
        return Err(FeatureError::SchemaMismatch);
    }
    Ok(())
}

pub(super) fn dssp_inverse(
    context: &StructuralContext<'_>,
) -> Result<Vec<Option<usize>>, FeatureError> {
    if context.dssp_index_for_residue.len() != context.ca_coords.len() {
        return Err(FeatureError::MissingContext("DSSP residue mapping"));
    }
    let mut inverse = vec![None; context.dssp.len + 1];
    for (clean_index, &dssp_index) in context.dssp_index_for_residue.iter().enumerate() {
        if dssp_index == 0 || dssp_index > context.dssp.len || inverse[dssp_index].is_some() {
            return Err(FeatureError::MissingContext("DSSP residue mapping"));
        }
        inverse[dssp_index] = Some(clean_index);
    }
    Ok(inverse)
}

fn mapped_partner(
    inverse: &[Option<usize>],
    dssp_index: usize,
    evidence_name: &'static str,
) -> Result<usize, FeatureError> {
    inverse
        .get(dssp_index)
        .and_then(|mapped| *mapped)
        .ok_or(FeatureError::MissingContext(evidence_name))
}

fn retained_hydrogen_bonds(
    context: &StructuralContext<'_>,
    inverse: &[Option<usize>],
) -> Result<BTreeMap<(usize, usize), i64>, FeatureError> {
    let mut bonds = BTreeMap::new();
    for (current_clean, &dssp_index) in context.dssp_index_for_residue.iter().enumerate() {
        let residue = context.dssp.get(dssp_index);
        for bond in residue.acceptor.iter().filter(|bond| bond.residue != 0) {
            let partner_clean =
                mapped_partner(inverse, bond.residue, "DSSP hydrogen bond partner")?;
            insert_bond(&mut bonds, (current_clean, partner_clean), bond.energy)?;
        }
        for bond in residue.donor.iter().filter(|bond| bond.residue != 0) {
            let partner_clean =
                mapped_partner(inverse, bond.residue, "DSSP hydrogen bond partner")?;
            insert_bond(&mut bonds, (partner_clean, current_clean), bond.energy)?;
        }
    }
    Ok(bonds)
}

fn insert_bond(
    bonds: &mut BTreeMap<(usize, usize), i64>,
    key: (usize, usize),
    energy: i64,
) -> Result<(), FeatureError> {
    match bonds.get(&key) {
        Some(existing) if *existing != energy => Err(FeatureError::MissingContext(
            "inconsistent mirrored DSSP hydrogen bond energy",
        )),
        Some(_) => Ok(()),
        None => {
            bonds.insert(key, energy);
            Ok(())
        }
    }
}

pub(super) fn retained_bridge_edges(
    context: &StructuralContext<'_>,
    inverse: &[Option<usize>],
) -> Result<BTreeSet<(usize, usize)>, FeatureError> {
    let mut edges = BTreeSet::new();
    for (current_clean, &dssp_index) in context.dssp_index_for_residue.iter().enumerate() {
        for partner in context
            .dssp
            .get(dssp_index)
            .partner
            .into_iter()
            .filter(|partner| *partner != 0)
        {
            let partner_clean = mapped_partner(inverse, partner, "DSSP bridge partner")?;
            if current_clean != partner_clean {
                edges.insert((
                    current_clean.min(partner_clean),
                    current_clean.max(partner_clean),
                ));
            }
        }
    }
    Ok(edges)
}

fn summary(values: &[f64]) -> [f64; 3] {
    [
        values.iter().copied().fold(f64::INFINITY, f64::min),
        values.iter().sum::<f64>() / values.len() as f64,
        values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    ]
}

fn validate_final(names: &[&'static str], values: &[f64]) -> Result<(), FeatureError> {
    if names.len() != values.len() {
        return Err(FeatureError::SchemaMismatch);
    }
    if let Some((name, _)) = names
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

    use crate::dssp::types::{BackboneResidue, HydrogenBond};
    use crate::dssp::DsspChain;
    use crate::peeling::contact_matrix::ContactMatrix;

    use super::extract_boundary_features;
    use crate::sword::factorized_ranker::partition::{
        parse_partition, FeatureError, ParsedPartition,
    };
    use crate::sword::factorized_ranker::schema::BOUNDARY_LOCAL_FEATURE_NAMES;
    use crate::sword::factorized_ranker::StructuralContext;

    fn named(values: &[f64], names: &[&str], name: &str) -> f64 {
        let index = names
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
        dssp_index_for_residue: Vec<usize>,
        contact_coords: Vec<[f64; 3]>,
        d0: f64,
        delta: f64,
    ) -> StructuralContext<'static> {
        let n_residues = dssp_index_for_residue.len();
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
            dssp_index_for_residue,
            contact_feature_cache: OnceLock::new(),
        }
    }

    fn boundary_fixture_with_mirrored_donor_acceptor_records(
        mirrored_energy: i64,
    ) -> (ParsedPartition, StructuralContext<'static>) {
        let mut dssp = dssp_chain(4);
        dssp.get_mut(1).acceptor[0] = HydrogenBond {
            residue: 4,
            energy: -1_200,
        };
        dssp.get_mut(4).donor[0] = HydrogenBond {
            residue: 1,
            energy: mirrored_energy,
        };
        dssp.get_mut(3).acceptor[0] = HydrogenBond {
            residue: 2,
            energy: -1_200,
        };
        dssp.get_mut(2).donor[0] = HydrogenBond {
            residue: 3,
            energy: -1_200,
        };
        (
            parse_partition("0-1 2-3", 4).unwrap(),
            context(dssp, (1..=4).collect(), vec![[0.0, 0.0, 0.0]; 4], 0.0, 1.0),
        )
    }

    #[test]
    fn crossing_bonds_are_deduplicated_and_converted_to_kcal() {
        let (partition, context) = boundary_fixture_with_mirrored_donor_acceptor_records(-1_200);
        let values = extract_boundary_features(&partition, &context).unwrap();
        assert_eq!(
            named(
                &values,
                BOUNDARY_LOCAL_FEATURE_NAMES,
                "boundary_hbond_count_mean",
            ),
            2.0,
        );
        assert!(
            (named(
                &values,
                BOUNDARY_LOCAL_FEATURE_NAMES,
                "boundary_hbond_energy_kcal_mean",
            ) - -2.4)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn insulation_uses_truncated_windows_and_contact_mean() {
        let partition = parse_partition("0 1-3", 4).unwrap();
        let context = context(
            dssp_chain(4),
            (1..=4).collect(),
            vec![[0.0, 0.0, 0.0]; 4],
            -3.0_f64.ln(),
            1.0,
        );
        let values = extract_boundary_features(&partition, &context).unwrap();
        assert!(
            (named(
                &values,
                BOUNDARY_LOCAL_FEATURE_NAMES,
                "boundary_insulation_w8_mean",
            ) - 0.75)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn mirrored_bond_energy_mismatch_fails_closed() {
        let (partition, context) = boundary_fixture_with_mirrored_donor_acceptor_records(-1_300);

        assert_eq!(
            extract_boundary_features(&partition, &context),
            Err(FeatureError::MissingContext(
                "inconsistent mirrored DSSP hydrogen bond energy",
            ))
        );
    }

    #[test]
    fn nonzero_unmapped_dssp_partner_fails_closed() {
        let mut dssp = dssp_chain(5);
        dssp.residues[2] = BackboneResidue::chain_break();
        dssp.get_mut(1).partner[0] = 2;
        let context = context(dssp, vec![1, 3, 4, 5], vec![[0.0, 0.0, 0.0]; 4], 0.0, 1.0);

        assert_eq!(
            extract_boundary_features(&parse_partition("0-1 2-3", 4).unwrap(), &context),
            Err(FeatureError::MissingContext("DSSP bridge partner"))
        );
    }

    #[test]
    fn one_domain_is_zero_but_multiple_domains_require_an_owner_transition() {
        let context = context(
            dssp_chain(4),
            (1..=4).collect(),
            vec![[0.0, 0.0, 0.0]; 4],
            0.0,
            1.0,
        );
        let one_domain =
            extract_boundary_features(&parse_partition("0-3", 4).unwrap(), &context).unwrap();
        assert_eq!(one_domain, vec![0.0; 42]);
        assert!(one_domain.iter().all(|value| value.is_finite()));

        let mut malformed = parse_partition("0-1 2-3", 4).unwrap();
        malformed.residue_to_domain.fill(0);
        assert_eq!(
            extract_boundary_features(&malformed, &context),
            Err(FeatureError::SchemaMismatch)
        );
    }

    #[test]
    fn boundary_sheet_links_use_labels_not_bridge_partners() {
        let mut dssp = dssp_chain(4);
        for index in 1..=4 {
            dssp.get_mut(index).sheet_label = 'A';
        }
        let context = context(dssp, (1..=4).collect(), vec![[0.0, 0.0, 0.0]; 4], 0.0, 1.0);
        let values =
            extract_boundary_features(&parse_partition("0-1 2-3", 4).unwrap(), &context).unwrap();

        assert_eq!(
            named(
                &values,
                BOUNDARY_LOCAL_FEATURE_NAMES,
                "boundary_sheet_link_count_mean",
            ),
            4.0
        );
        assert_eq!(
            named(
                &values,
                BOUNDARY_LOCAL_FEATURE_NAMES,
                "boundary_bridge_count_mean",
            ),
            0.0
        );
    }
}
