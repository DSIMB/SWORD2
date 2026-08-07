#![allow(dead_code)]

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

use crate::peeling::contact_matrix::ContactMatrix;
use crate::sword::geometry_metrics::{
    gyration_tensor, principal_radii, radius_of_gyration, sym3x3_eigenvalues,
};

use super::boundary::extract_boundary_features;
use super::discontinuity::extract_discontinuity_features;
use super::lattice::CandidateRecord;
use super::partition::{Domain, FeatureError, ParsedPartition, Segment};
use super::schema::{
    CandidateFeatures, FeatureMask, GlobalFeatures, BASE_CANDIDATE_FEATURE_NAMES,
    BOUNDARY_LOCAL_END, BOUNDARY_LOCAL_FEATURE_NAMES, BOUNDARY_LOCAL_START,
    CANDIDATE_FEATURE_NAMES, DISCONTINUITY_END, DISCONTINUITY_FEATURE_NAMES, DISCONTINUITY_START,
    DOMAIN_CONDITIONED_FEATURE_NAMES, GLOBAL_FEATURE_NAMES,
};
use super::StructuralContext;

const NONLOCAL_SEQUENCE_CUTOFF: usize = 8;
const V0: f64 = 141.0;
const RHO_IDEAL: f64 = 1.0 / V0;
const I3: f64 = 780.69552;
const I4: f64 = 6137.0618;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SsClass {
    Helix,
    Strand,
    Coil,
}

#[derive(Debug, Clone, Copy)]
struct ShapeMeasures {
    q1: f64,
    q2: f64,
    q3: f64,
    volume_ratio: f64,
    relative_density: f64,
}

#[derive(Debug, Clone, Copy)]
struct DomainMeasures {
    size_fraction: f64,
    shape: ShapeMeasures,
    internal_contact_density: f64,
    contact_order: f64,
    internal_contact_fraction: f64,
    conductance: f64,
}

pub(crate) struct ContactFeatureCache {
    n: usize,
    long_range_sum: Vec<f64>,
    long_range_count: Vec<f64>,
    weighted_separation_sum: Vec<f64>,
    validation: ContactValidation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ContactValidation {
    Valid,
    NonFinite,
    OutOfRange,
    Asymmetric,
}

impl ContactFeatureCache {
    pub(super) fn new(contacts: &ContactMatrix) -> Self {
        let n = contacts.len();
        let side = n + 1;
        let mut long_range_sum = vec![0.0; side * side];
        let mut long_range_count = vec![0.0; side * side];
        let mut weighted_separation_sum = vec![0.0; side * side];
        let mut validation = ContactValidation::Valid;
        for row in 1..=n {
            for column in 1..=n {
                let i = row - 1;
                let j = column - 1;
                let separation = i.abs_diff(j);
                let probability = contacts.get(i, j);
                if !probability.is_finite() {
                    validation = ContactValidation::NonFinite;
                } else if !(0.0..=1.0).contains(&probability)
                    && validation == ContactValidation::Valid
                {
                    validation = ContactValidation::OutOfRange;
                } else if (probability - contacts.get(j, i)).abs() > 1e-12
                    && validation == ContactValidation::Valid
                {
                    validation = ContactValidation::Asymmetric;
                }
                let index = row * side + column;
                let above = (row - 1) * side + column;
                let left = row * side + column - 1;
                let diagonal = (row - 1) * side + column - 1;
                let eligible = separation >= NONLOCAL_SEQUENCE_CUTOFF;
                long_range_sum[index] = if eligible { probability } else { 0.0 }
                    + long_range_sum[above]
                    + long_range_sum[left]
                    - long_range_sum[diagonal];
                long_range_count[index] = if eligible { 1.0 } else { 0.0 }
                    + long_range_count[above]
                    + long_range_count[left]
                    - long_range_count[diagonal];
                weighted_separation_sum[index] = probability * separation as f64
                    + weighted_separation_sum[above]
                    + weighted_separation_sum[left]
                    - weighted_separation_sum[diagonal];
            }
        }
        Self {
            n,
            long_range_sum,
            long_range_count,
            weighted_separation_sum,
            validation,
        }
    }

    pub(super) fn validate(&self) -> Result<(), FeatureError> {
        match self.validation {
            ContactValidation::Valid => Ok(()),
            ContactValidation::NonFinite => Err(FeatureError::NonFinite("contact probability")),
            ContactValidation::OutOfRange => Err(FeatureError::MissingContext(
                "contact probability outside [0, 1]",
            )),
            ContactValidation::Asymmetric => {
                Err(FeatureError::MissingContext("asymmetric contact matrix"))
            }
        }
    }

    fn prefix_rectangle(
        &self,
        prefix: &[f64],
        row_start: usize,
        col_start: usize,
        row_end: usize,
        col_end: usize,
    ) -> f64 {
        let side = self.n + 1;
        let row_after = row_end + 1;
        let col_after = col_end + 1;
        prefix[row_after * side + col_after]
            - prefix[row_start * side + col_after]
            - prefix[row_after * side + col_start]
            + prefix[row_start * side + col_start]
    }

    fn all_rectangle(&self, contacts: &ContactMatrix, left: Segment, right: Segment) -> f64 {
        contacts.rectangle_sum(left.start, right.start, left.end, right.end)
    }

    pub(crate) fn long_rectangle(&self, left: Segment, right: Segment) -> (f64, f64) {
        (
            self.prefix_rectangle(
                &self.long_range_sum,
                left.start,
                right.start,
                left.end,
                right.end,
            ),
            self.prefix_rectangle(
                &self.long_range_count,
                left.start,
                right.start,
                left.end,
                right.end,
            ),
        )
    }

    fn weighted_rectangle(&self, left: Segment, right: Segment) -> f64 {
        self.prefix_rectangle(
            &self.weighted_separation_sum,
            left.start,
            right.start,
            left.end,
            right.end,
        )
    }

    fn segment_internal(&self, contacts: &ContactMatrix, segment: Segment) -> (f64, f64, f64, f64) {
        let diagonal_mass = (segment.start..=segment.end)
            .map(|index| contacts.get(index, index))
            .sum::<f64>();
        let all_mass = (self.all_rectangle(contacts, segment, segment) - diagonal_mass) / 2.0;
        let (long_mass, long_count) = self.long_rectangle(segment, segment);
        (
            all_mass,
            long_mass / 2.0,
            long_count / 2.0,
            self.weighted_rectangle(segment, segment) / 2.0,
        )
    }

    fn domain_internal(&self, contacts: &ContactMatrix, domain: &Domain) -> (f64, f64, f64, f64) {
        let mut all_mass = 0.0;
        let mut long_mass = 0.0;
        let mut long_count = 0.0;
        let mut weighted_separation = 0.0;
        for (index, &left) in domain.segments.iter().enumerate() {
            let (segment_all, segment_long, segment_count, segment_weighted) =
                self.segment_internal(contacts, left);
            all_mass += segment_all;
            long_mass += segment_long;
            long_count += segment_count;
            weighted_separation += segment_weighted;
            for &right in &domain.segments[index + 1..] {
                all_mass += self.all_rectangle(contacts, left, right);
                let (pair_long, pair_count) = self.long_rectangle(left, right);
                long_mass += pair_long;
                long_count += pair_count;
                weighted_separation += self.weighted_rectangle(left, right);
            }
        }
        (all_mass, long_mass, long_count, weighted_separation)
    }

    fn domain_pair_mass(&self, contacts: &ContactMatrix, left: &Domain, right: &Domain) -> f64 {
        left.segments
            .iter()
            .flat_map(|&left_segment| {
                right.segments.iter().map(move |&right_segment| {
                    self.all_rectangle(contacts, left_segment, right_segment)
                })
            })
            .sum()
    }

    fn chain_totals(&self, contacts: &ContactMatrix) -> (f64, f64, f64, f64) {
        let whole = Segment {
            start: 0,
            end: self.n - 1,
        };
        self.segment_internal(contacts, whole)
    }
}

pub(crate) fn extract_global_features(
    context: &StructuralContext<'_>,
    candidate_counts: &[usize],
) -> Result<GlobalFeatures, FeatureError> {
    let contact_cache = context
        .contact_feature_cache
        .get_or_init(|| ContactFeatureCache::new(context.contacts));
    context.validate(FeatureMask {
        global_count: true,
        domain_conditioned: false,
        boundary_local: false,
        relative_hierarchy: false,
        discontinuity: false,
    })?;
    if candidate_counts.contains(&0) {
        return Err(FeatureError::SchemaMismatch);
    }

    let n = context.ca_coords.len();
    let indices: Vec<usize> = (0..n).collect();
    let eigenvalues =
        checked_gyration_eigenvalues(context.ca_coords, &indices, "chain_rg_normalized")?;
    let radii = principal_radii(eigenvalues);
    let (contact_mass, nonlocal_mass, nonlocal_count, weighted_separation) =
        contact_cache.chain_totals(context.contacts);
    let classes: Vec<SsClass> = (0..n).map(|index| ss_class(context, index)).collect();
    let helix_count = classes
        .iter()
        .filter(|&&class| class == SsClass::Helix)
        .count();
    let strand_count = classes
        .iter()
        .filter(|&&class| class == SsClass::Strand)
        .count();
    let coil_count = n - helix_count - strand_count;

    let mut count_histogram = [0.0; 21];
    let mut frequencies = BTreeMap::new();
    for &count in candidate_counts {
        let bin = if count <= 20 { count - 1 } else { 20 };
        count_histogram[bin] += 1.0;
        *frequencies.entry(count).or_insert(0usize) += 1;
    }
    if !candidate_counts.is_empty() {
        for value in &mut count_histogram {
            *value /= candidate_counts.len() as f64;
        }
    }
    let modal_count = frequencies
        .into_iter()
        .max_by_key(|&(count, frequency)| (frequency, std::cmp::Reverse(count)))
        .map(|(count, _)| count)
        .unwrap_or(0);

    let values = GlobalFeatures {
        n_residues: n as f64,
        rg_normalized: radius_of_gyration(eigenvalues) / (n as f64).cbrt(),
        inertia_ratio_21: ratio(radii[1], radii[0]),
        inertia_ratio_31: ratio(radii[2], radii[0]),
        nonlocal_contact_density: ratio(nonlocal_mass, nonlocal_count),
        contact_order: ratio(weighted_separation, contact_mass) / n as f64,
        helix_fraction: helix_count as f64 / n as f64,
        strand_fraction: strand_count as f64 / n as f64,
        coil_fraction: coil_count as f64 / n as f64,
        helix_blocks: block_count(&classes, SsClass::Helix) as f64,
        strand_blocks: block_count(&classes, SsClass::Strand) as f64,
        peeling_levels: context.iterations.len() as f64,
        finest_pus: context
            .iterations
            .iter()
            .map(|iteration| iteration.num_pus)
            .max()
            .unwrap_or(0) as f64,
        candidate_total: candidate_counts.len() as f64,
        available_count_total: candidate_counts
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .len() as f64,
        count_histogram,
        modal_count: modal_count as f64,
    };
    validate_named(GLOBAL_FEATURE_NAMES, &values.to_vec())?;
    Ok(values)
}

pub(crate) fn extract_candidate_base_and_domain(
    candidate: &CandidateRecord,
    context: &StructuralContext<'_>,
    modal_count: usize,
) -> Result<CandidateFeatures, FeatureError> {
    extract_candidate_base_and_domain_with_mask(
        candidate,
        context,
        modal_count,
        FeatureMask {
            global_count: false,
            domain_conditioned: true,
            boundary_local: false,
            relative_hierarchy: false,
            discontinuity: false,
        },
    )
}

pub(crate) fn extract_candidate_base_and_domain_with_mask(
    candidate: &CandidateRecord,
    context: &StructuralContext<'_>,
    modal_count: usize,
    mask: FeatureMask,
) -> Result<CandidateFeatures, FeatureError> {
    let cache = context
        .contact_feature_cache
        .get_or_init(|| ContactFeatureCache::new(context.contacts));
    context.validate(FeatureMask {
        global_count: false,
        domain_conditioned: mask.domain_conditioned,
        boundary_local: false,
        relative_hierarchy: false,
        discontinuity: false,
    })?;
    validate_partition(candidate, context.ca_coords.len())?;

    let shapes: Vec<ShapeMeasures> = candidate
        .partition
        .domains
        .iter()
        .map(|domain| domain_shape(context.ca_coords, &domain.residues, "domain_q1_mean"))
        .collect::<Result<_, _>>()?;
    let sizes: Vec<usize> = candidate
        .partition
        .domains
        .iter()
        .map(|domain| domain.residues.len())
        .collect();
    let mut external_masses = mask.domain_conditioned.then(|| vec![0.0; sizes.len()]);
    let mut contact_q_values = Vec::new();
    for left in 0..sizes.len() {
        for right in left + 1..sizes.len() {
            let mass = cache.domain_pair_mass(
                context.contacts,
                &candidate.partition.domains[left],
                &candidate.partition.domains[right],
            );
            if let Some(external_masses) = &mut external_masses {
                external_masses[left] += mass;
                external_masses[right] += mass;
            }
            contact_q_values.push(contact_q(mass, sizes[left], sizes[right]));
        }
    }

    let segment_sizes: Vec<usize> = candidate
        .partition
        .domains
        .iter()
        .flat_map(|domain| domain.segments.iter())
        .map(|segment| segment.end - segment.start + 1)
        .collect();
    let total_size = sizes.iter().sum::<usize>();
    let smallest_size = *sizes.iter().min().ok_or(FeatureError::Malformed)?;
    let largest_size = *sizes.iter().max().ok_or(FeatureError::Malformed)?;
    let mut values = vec![
        candidate.measure.num_domains as f64,
        candidate.measure.min_size as f64,
        candidate.measure.max_cr,
        candidate.measure.density_min,
        candidate.measure.mean_density,
        boundary_coil_fraction(&candidate.partition, context),
        candidate.measure.num_domains.abs_diff(modal_count) as f64,
        mean(&shapes.iter().map(|shape| shape.q1).collect::<Vec<_>>()),
        mean(&shapes.iter().map(|shape| shape.q2).collect::<Vec<_>>()),
        mean(&shapes.iter().map(|shape| shape.q3).collect::<Vec<_>>()),
        shapes
            .iter()
            .map(|shape| shape.q3)
            .fold(f64::INFINITY, f64::min),
        mean(
            &shapes
                .iter()
                .map(|shape| shape.volume_ratio)
                .collect::<Vec<_>>(),
        ),
        mean(
            &shapes
                .iter()
                .map(|shape| shape.relative_density)
                .collect::<Vec<_>>(),
        ),
        shapes
            .iter()
            .map(|shape| shape.relative_density)
            .fold(f64::INFINITY, f64::min),
        mean(&contact_q_values),
        contact_q_values.iter().copied().fold(0.0, f64::max),
        segment_sizes.len() as f64,
        segment_sizes.len().saturating_sub(sizes.len()) as f64,
        smallest_size as f64 / largest_size as f64,
        largest_size as f64 / total_size as f64,
        *segment_sizes.iter().min().ok_or(FeatureError::Malformed)? as f64,
        mean(
            &segment_sizes
                .iter()
                .map(|&size| size as f64)
                .collect::<Vec<_>>(),
        ),
    ];
    debug_assert_eq!(values.len(), BASE_CANDIDATE_FEATURE_NAMES.len());
    validate_named(BASE_CANDIDATE_FEATURE_NAMES, &values)?;
    values.resize(CANDIDATE_FEATURE_NAMES.len(), 0.0);

    if let Some(external_masses) = external_masses {
        let mut domain_measures = Vec::with_capacity(sizes.len());
        for (index, domain) in candidate.partition.domains.iter().enumerate() {
            let (internal_mass, nonlocal_mass, nonlocal_count, weighted_separation) =
                cache.domain_internal(context.contacts, domain);
            let external_mass = external_masses[index];
            domain_measures.push(DomainMeasures {
                size_fraction: sizes[index] as f64 / context.ca_coords.len() as f64,
                shape: shapes[index],
                internal_contact_density: ratio(nonlocal_mass, nonlocal_count),
                contact_order: ratio(weighted_separation, internal_mass)
                    / context.ca_coords.len() as f64,
                internal_contact_fraction: ratio(internal_mass, internal_mass + external_mass),
                conductance: ratio(external_mass, 2.0 * internal_mass + external_mass),
            });
        }
        let domain_values = domain_conditioned_values(&candidate.partition, &domain_measures)?;
        let domain_start = BASE_CANDIDATE_FEATURE_NAMES.len();
        let domain_end = domain_start + DOMAIN_CONDITIONED_FEATURE_NAMES.len();
        values[domain_start..domain_end].copy_from_slice(&domain_values);
    }

    validate_named(CANDIDATE_FEATURE_NAMES, &values)?;

    Ok(CandidateFeatures {
        source_index: candidate.source_index,
        canonical: candidate.partition.canonical.clone(),
        num_domains: candidate.measure.num_domains,
        values,
    })
}

pub(crate) fn populate_candidate_conditional_features(
    candidate: &CandidateRecord,
    features: &mut CandidateFeatures,
    context: &StructuralContext<'_>,
    mask: FeatureMask,
) -> Result<(), FeatureError> {
    if features.source_index != candidate.source_index
        || features.canonical != candidate.partition.canonical
        || features.num_domains != candidate.measure.num_domains
        || features.values.len() != CANDIDATE_FEATURE_NAMES.len()
        || BOUNDARY_LOCAL_END - BOUNDARY_LOCAL_START != BOUNDARY_LOCAL_FEATURE_NAMES.len()
        || DISCONTINUITY_END - DISCONTINUITY_START != DISCONTINUITY_FEATURE_NAMES.len()
        || DISCONTINUITY_END != CANDIDATE_FEATURE_NAMES.len()
    {
        return Err(FeatureError::SchemaMismatch);
    }

    let discontinuity_values = if mask.discontinuity {
        let values = extract_discontinuity_features(&candidate.partition, context)?;
        if values.len() != DISCONTINUITY_FEATURE_NAMES.len() {
            return Err(FeatureError::SchemaMismatch);
        }
        values
    } else {
        vec![0.0; DISCONTINUITY_FEATURE_NAMES.len()]
    };
    let boundary_values = if mask.boundary_local {
        let values = extract_boundary_features(&candidate.partition, context)?;
        if values.len() != BOUNDARY_LOCAL_FEATURE_NAMES.len() {
            return Err(FeatureError::SchemaMismatch);
        }
        values
    } else {
        vec![0.0; BOUNDARY_LOCAL_FEATURE_NAMES.len()]
    };

    features.values[BOUNDARY_LOCAL_START..BOUNDARY_LOCAL_END].copy_from_slice(&boundary_values);
    features.values[DISCONTINUITY_START..DISCONTINUITY_END].copy_from_slice(&discontinuity_values);
    Ok(())
}

fn validate_partition(candidate: &CandidateRecord, chain_len: usize) -> Result<(), FeatureError> {
    if candidate.partition.residue_to_domain.len() != chain_len
        || candidate.partition.domains.is_empty()
        || candidate
            .partition
            .domains
            .iter()
            .any(|domain| domain.residues.is_empty() || domain.segments.is_empty())
    {
        return Err(FeatureError::Malformed);
    }
    if candidate.partition.domains.len() != candidate.measure.num_domains {
        return Err(FeatureError::DomainCountMismatch);
    }
    Ok(())
}

fn domain_conditioned_values(
    partition: &ParsedPartition,
    measures: &[DomainMeasures],
) -> Result<Vec<f64>, FeatureError> {
    if measures.len() != partition.domains.len() || measures.is_empty() {
        return Err(FeatureError::SchemaMismatch);
    }
    let smallest_index = (0..measures.len())
        .min_by_key(|&index| {
            (
                partition.domains[index].residues.len(),
                partition.domains[index].residues[0],
            )
        })
        .ok_or(FeatureError::Malformed)?;
    let largest_index = (0..measures.len())
        .min_by_key(|&index| {
            (
                std::cmp::Reverse(partition.domains[index].residues.len()),
                partition.domains[index].residues[0],
            )
        })
        .ok_or(FeatureError::Malformed)?;
    let smallest = measures[smallest_index];
    let largest = measures[largest_index];
    let internal_fractions: Vec<f64> = measures
        .iter()
        .map(|measure| measure.internal_contact_fraction)
        .collect();
    let conductances: Vec<f64> = measures.iter().map(|measure| measure.conductance).collect();
    let (internal_min, internal_mean, internal_max) = summary(&internal_fractions);
    let (conductance_min, conductance_mean, conductance_max) = summary(&conductances);
    let values = vec![
        smallest.size_fraction,
        smallest.shape.q1,
        smallest.shape.q2,
        smallest.shape.q3,
        smallest.shape.relative_density,
        smallest.internal_contact_density,
        smallest.contact_order,
        smallest.internal_contact_fraction,
        largest.size_fraction,
        largest.shape.q1,
        largest.shape.q2,
        largest.shape.q3,
        largest.shape.relative_density,
        largest.internal_contact_density,
        largest.contact_order,
        largest.internal_contact_fraction,
        ratio(
            smallest.shape.relative_density,
            largest.shape.relative_density,
        ),
        ratio(
            smallest.internal_contact_density,
            largest.internal_contact_density,
        ),
        ratio(smallest.shape.q1, largest.shape.q1),
        ratio(smallest.shape.q2, largest.shape.q2),
        ratio(smallest.shape.q3, largest.shape.q3),
        internal_min,
        internal_mean,
        internal_max,
        conductance_min,
        conductance_mean,
        conductance_max,
    ];
    validate_named(DOMAIN_CONDITIONED_FEATURE_NAMES, &values)?;
    Ok(values)
}

fn checked_gyration_eigenvalues(
    coords: &[[f64; 3]],
    indices: &[usize],
    feature_name: &'static str,
) -> Result<[f64; 3], FeatureError> {
    let tensor = gyration_tensor(coords, indices);
    if tensor.iter().flatten().any(|value| !value.is_finite()) {
        return Err(FeatureError::NonFinite(feature_name));
    }

    // Cardano's closed form squares tensor elements before normalizing them.
    // Reject magnitudes that can overflow those intermediates even though the
    // input coordinates and tensor elements themselves are finite.
    let max_component = tensor
        .iter()
        .flatten()
        .map(|value| value.abs())
        .fold(0.0, f64::max);
    if max_component > f64::MAX.sqrt() / 8.0 {
        return Err(FeatureError::NonFinite(feature_name));
    }

    let eigenvalues = sym3x3_eigenvalues(tensor);
    if eigenvalues.iter().any(|value| !value.is_finite()) {
        return Err(FeatureError::NonFinite(feature_name));
    }
    Ok(eigenvalues)
}

fn domain_shape(
    coords: &[[f64; 3]],
    indices: &[usize],
    feature_name: &'static str,
) -> Result<ShapeMeasures, FeatureError> {
    if indices.len() < 4 {
        return Ok(ShapeMeasures {
            q1: 0.0,
            q2: 0.0,
            q3: 0.0,
            volume_ratio: 0.0,
            relative_density: 0.0,
        });
    }
    let radii = principal_radii(checked_gyration_eigenvalues(coords, indices, feature_name)?);
    let c_ideal = (3.0 * V0 / (4.0 * PI)).cbrt();
    let ideal_radius = c_ideal / 5.0f64.sqrt() * (indices.len() as f64).cbrt();
    let q1 = radii[0] / ideal_radius.max(1e-12);
    let q2 = radii[1] / ideal_radius.max(1e-12);
    let q3 = radii[2] / ideal_radius.max(1e-12);
    let ellipsoid_volume = (20.0 * 5.0f64.sqrt() / 3.0)
        * PI
        * radii[0].max(1e-3)
        * radii[1].max(1e-3)
        * radii[2].max(1e-3);
    Ok(ShapeMeasures {
        q1,
        q2,
        q3,
        volume_ratio: q1 * q2 * q3,
        relative_density: (indices.len() as f64 / ellipsoid_volume) / RHO_IDEAL,
    })
}

fn contact_q(probability: f64, size_a: usize, size_b: usize) -> f64 {
    if size_a < 4 || size_b < 4 {
        return 0.0;
    }
    let c_ideal = (3.0 * V0 / (4.0 * PI)).cbrt();
    let k_min = (2.0 / 3.0) * PI.powi(2) * RHO_IDEAL.powi(2) * I4;
    let k_max = PI.powi(2) * RHO_IDEAL.powi(2) * I3;
    let radius_a = c_ideal * (size_a as f64).cbrt();
    let radius_b = c_ideal * (size_b as f64).cbrt();
    let min_contacts = k_min * (radius_a * radius_b) / (radius_a + radius_b);
    let radius_ab = c_ideal * ((size_a + size_b) as f64).cbrt();
    let fraction = size_a.min(size_b) as f64 / (size_a + size_b) as f64;
    let angle = (2.0 * fraction - 1.0).clamp(-1.0, 1.0).acos() / 3.0;
    let x_value = 2.0 * (angle + 4.0 * PI / 3.0).cos();
    let max_contacts = k_max * radius_ab.powi(2) * (1.0 - x_value.powi(2));
    let denominator = max_contacts - min_contacts;
    if denominator.abs() <= 1e-5 {
        0.0
    } else {
        ((probability - min_contacts) / denominator).clamp(0.0, 1.0)
    }
}

fn boundary_coil_fraction(partition: &ParsedPartition, context: &StructuralContext<'_>) -> f64 {
    const WINDOW: usize = 2;
    let mut total = 0usize;
    let mut coil = 0usize;
    for segment in partition
        .domains
        .iter()
        .flat_map(|domain| domain.segments.iter())
    {
        for junction in [segment.start, segment.end] {
            let start = junction.saturating_sub(WINDOW);
            let end = (junction + WINDOW).min(context.ca_coords.len() - 1);
            for index in start..=end {
                total += 1;
                if ss_class(context, index) == SsClass::Coil {
                    coil += 1;
                }
            }
        }
    }
    ratio(coil as f64, total as f64)
}

fn ss_class(context: &StructuralContext<'_>, clean_index: usize) -> SsClass {
    match context
        .dssp
        .get(context.dssp_index_for_residue[clean_index])
        .ss[0]
    {
        'H' | 'G' | 'I' => SsClass::Helix,
        'E' | 'B' => SsClass::Strand,
        _ => SsClass::Coil,
    }
}

fn block_count(classes: &[SsClass], target: SsClass) -> usize {
    classes
        .iter()
        .enumerate()
        .filter(|&(index, &class)| class == target && (index == 0 || classes[index - 1] != target))
        .count()
}

fn ratio(numerator: f64, denominator: f64) -> f64 {
    if denominator == 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

fn summary(values: &[f64]) -> (f64, f64, f64) {
    if values.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    (
        values.iter().copied().fold(f64::INFINITY, f64::min),
        mean(values),
        values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    )
}

fn validate_named(names: &[&'static str], values: &[f64]) -> Result<(), FeatureError> {
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
    use std::collections::BTreeSet;

    use super::{
        domain_conditioned_values, extract_candidate_base_and_domain,
        extract_candidate_base_and_domain_with_mask, extract_global_features,
        populate_candidate_conditional_features, DomainMeasures, ShapeMeasures,
    };
    use crate::dssp::types::BackboneResidue;
    use crate::dssp::DsspChain;
    use crate::peeling::algorithm::IterationResult;
    use crate::peeling::contact_matrix::ContactMatrix;
    use crate::sword::compute_measure::{MeasureLine, MeasureProvenance};
    use crate::sword::factorized_ranker::lattice::CandidateRecord;
    use crate::sword::factorized_ranker::partition::{parse_partition, FeatureError};
    use crate::sword::factorized_ranker::schema::{
        candidate_pair_vector, CandidateFeatures, FeatureMask, BASE_AND_DOMAIN_FEATURE_NAMES,
        BASE_CANDIDATE_FEATURE_NAMES, BOUNDARY_LOCAL_END, BOUNDARY_LOCAL_FEATURE_NAMES,
        BOUNDARY_LOCAL_START, CANDIDATE_FEATURE_NAMES, COUNT_ITEM_FEATURE_NAMES, DISCONTINUITY_END,
        DISCONTINUITY_FEATURE_NAMES, DISCONTINUITY_START, DOMAIN_CONDITIONED_FEATURE_NAMES,
        FEATURE_SCHEMA_VERSION, GLOBAL_FEATURE_NAMES, RELATIVE_HIERARCHY_FEATURE_NAMES,
    };
    use crate::sword::factorized_ranker::StructuralContext;

    fn dssp_chain(sequence: &str) -> DsspChain {
        let mut chain = DsspChain::new();
        for code in sequence.chars() {
            let mut residue = BackboneResidue {
                aa: 'A',
                ..BackboneResidue::default()
            };
            residue.ss[0] = code;
            chain.push(residue);
        }
        chain
    }

    fn dssp(sequence: &str) -> &'static DsspChain {
        Box::leak(Box::new(dssp_chain(sequence)))
    }

    fn iteration(chain_len: usize, num_pus: usize) -> IterationResult {
        IterationResult {
            max_cr: 0.0,
            min_density: 0.0,
            ci: 0.0,
            r: 0.0,
            num_pus,
            pu_boundaries: vec![[0, chain_len - 1]],
        }
    }

    fn context(
        ca_coords: Vec<[f64; 3]>,
        contact_coords: Vec<[f64; 3]>,
        d0: f64,
        delta: f64,
        sequence: &str,
        iterations: Vec<IterationResult>,
        provenance: Vec<MeasureProvenance>,
    ) -> StructuralContext<'static> {
        let chain_len = ca_coords.len();
        StructuralContext {
            ca_coords: Box::leak(ca_coords.into_boxed_slice()),
            dssp: dssp(sequence),
            contacts: Box::leak(Box::new(ContactMatrix::from_ca_coords(
                &contact_coords,
                d0,
                delta,
            ))),
            iterations: Box::leak(iterations.into_boxed_slice()),
            measure_provenance: Box::leak(provenance.into_boxed_slice()),
            dssp_index_for_residue: (1..=chain_len).collect(),
            contact_feature_cache: std::sync::OnceLock::new(),
        }
    }

    fn contact_fixture() -> StructuralContext<'static> {
        let coords: Vec<[f64; 3]> = (0..9).map(|index| [index as f64, 0.0, 0.0]).collect();
        context(
            coords.clone(),
            coords,
            11.605_671_529_869_176,
            -3.282_023_664_818_524_5,
            "CCCCCCCCC",
            vec![iteration(9, 3)],
            vec![],
        )
    }

    fn fixture_context() -> StructuralContext<'static> {
        let coords: Vec<[f64; 3]> = (0..14).map(|index| [index as f64, 0.0, 0.0]).collect();
        context(
            coords,
            vec![[0.0, 0.0, 0.0]; 14],
            0.0,
            1.0,
            "HHHCCEEECCCCCC",
            vec![iteration(14, 5)],
            vec![],
        )
    }

    fn fixture_candidate() -> CandidateRecord {
        CandidateRecord {
            source_index: 7,
            measure: MeasureLine {
                num_domains: 2,
                min_size: 4,
                delineation: "0-3 4-13".to_string(),
                max_cr: 0.25,
                mean_cr: 0.0,
                density_min: 1.5,
                mean_density: 2.5,
            },
            partition: parse_partition("0-3 4-13", 14).unwrap(),
            legacy_distance: 0.0,
            hierarchy: None,
        }
    }

    fn base_only_mask() -> FeatureMask {
        FeatureMask {
            global_count: false,
            domain_conditioned: false,
            boundary_local: false,
            relative_hierarchy: false,
            discontinuity: false,
        }
    }

    fn extract_fixture_features() -> Result<
        (
            crate::sword::factorized_ranker::schema::GlobalFeatures,
            CandidateFeatures,
        ),
        FeatureError,
    > {
        let context = fixture_context();
        Ok((
            extract_global_features(&context, &[2, 2, 4, 21])?,
            extract_candidate_base_and_domain(&fixture_candidate(), &context, 2)?,
        ))
    }

    fn named(candidate: &CandidateFeatures, name: &str) -> f64 {
        candidate.base_and_domain_value(name).unwrap()
    }

    #[test]
    fn global_contact_features_exclude_local_pairs() {
        let context = contact_fixture();
        let values = extract_global_features(&context, &[2, 2, 3]).unwrap();
        assert!((values.nonlocal_contact_density - 0.25).abs() < 1e-12);
        assert!((values.contact_order - 0.5).abs() < 1e-12);
    }

    #[test]
    fn global_geometry_dssp_and_histogram_match_reference_formulas() {
        let coords = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let context = context(
            coords,
            vec![[0.0, 0.0, 0.0]; 6],
            0.0,
            1.0,
            "HHHCEE",
            vec![iteration(6, 3), iteration(6, 5), iteration(6, 4)],
            vec![],
        );

        let values = extract_global_features(&context, &[2, 2, 4, 21]).unwrap();

        assert!((values.rg_normalized - 1.0 / 6.0f64.cbrt()).abs() < 1e-12);
        assert!((values.inertia_ratio_21 - 1.0).abs() < 1e-12);
        assert!((values.inertia_ratio_31 - 1.0).abs() < 1e-12);
        assert!((values.helix_fraction - 0.5).abs() < 1e-12);
        assert!((values.strand_fraction - 2.0 / 6.0).abs() < 1e-12);
        assert!((values.coil_fraction - 1.0 / 6.0).abs() < 1e-12);
        assert_eq!(values.helix_blocks, 1.0);
        assert_eq!(values.strand_blocks, 1.0);
        assert_eq!(values.peeling_levels, 3.0);
        assert_eq!(values.finest_pus, 5.0);
        assert_eq!(values.count_histogram[1], 0.5);
        assert_eq!(values.count_histogram[3], 0.25);
        assert_eq!(values.count_histogram[20], 0.25);
        assert!((values.count_histogram.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert_eq!(values.modal_count, 2.0);

        let tied = extract_global_features(&context, &[3, 2]).unwrap();
        assert_eq!(tied.modal_count, 2.0);
        let empty = extract_global_features(&context, &[]).unwrap();
        assert_eq!(empty.modal_count, 0.0);
        assert!(empty.count_histogram.iter().all(|value| *value == 0.0));
    }

    #[test]
    fn asymmetric_domain_features_preserve_smallest_and_largest() {
        let values =
            extract_candidate_base_and_domain(&fixture_candidate(), &fixture_context(), 2).unwrap();
        let smallest = named(&values, "smallest_size_fraction");
        let largest = named(&values, "largest_size_fraction");
        let q1_ratio = named(&values, "smallest_to_largest_q1_ratio");
        let smallest_q1 = named(&values, "smallest_q1");
        let largest_q1 = named(&values, "largest_q1");
        assert!(smallest < largest);
        assert!((q1_ratio - smallest_q1 / largest_q1).abs() < 1e-12);
        assert!((smallest_q1 - 0.487_764_112_577_698).abs() < 1e-12);
        assert!((largest_q1 - 0.923_283_643_796_329_1).abs() < 1e-12);
        assert!((named(&values, "domain_q1_mean") - 0.705_523_878_187_013_6).abs() < 1e-12);
        assert!((named(&values, "contact_q_mean") - 0.906_850_466_360_395_3).abs() < 1e-12);
        assert!((named(&values, "boundary_coil_fraction") - 7.0 / 16.0).abs() < 1e-12);
        assert_eq!(named(&values, "smallest_to_largest_q2_ratio"), 0.0);
    }

    #[test]
    fn domain_contact_masses_use_unordered_internal_and_directed_external_pairs() {
        let values =
            extract_candidate_base_and_domain(&fixture_candidate(), &fixture_context(), 2).unwrap();

        assert_eq!(named(&values, "smallest_internal_contact_density"), 0.0);
        assert!((named(&values, "largest_internal_contact_density") - 0.5).abs() < 1e-12);
        assert!((named(&values, "smallest_contact_order") - 5.0 / 42.0).abs() < 1e-12);
        assert!((named(&values, "largest_contact_order") - 11.0 / 42.0).abs() < 1e-12);
        assert!(
            (named(&values, "domain_internal_contact_fraction_min") - 3.0 / 23.0).abs() < 1e-12
        );
        assert!(
            (named(&values, "domain_internal_contact_fraction_max") - 9.0 / 17.0).abs() < 1e-12
        );
        assert!((named(&values, "domain_conductance_min") - 4.0 / 13.0).abs() < 1e-12);
        assert!((named(&values, "domain_conductance_max") - 10.0 / 13.0).abs() < 1e-12);
    }

    #[test]
    fn equal_size_domain_ties_choose_lower_minimum_residue() {
        let coords = vec![
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [10.0, 0.0, 0.0],
            [11.0, 0.0, 0.0],
            [12.0, 0.0, 0.0],
            [13.0, 0.0, 0.0],
        ];
        let context = context(
            coords,
            vec![[0.0, 0.0, 0.0]; 8],
            0.0,
            1.0,
            "CCCCCCCC",
            vec![iteration(8, 2)],
            vec![],
        );
        let candidate = CandidateRecord {
            source_index: 0,
            measure: MeasureLine {
                num_domains: 2,
                min_size: 4,
                delineation: "0-3 4-7".into(),
                max_cr: 0.0,
                mean_cr: 0.0,
                density_min: 0.0,
                mean_density: 0.0,
            },
            partition: parse_partition("0-3 4-7", 8).unwrap(),
            legacy_distance: 0.0,
            hierarchy: None,
        };

        let values = extract_candidate_base_and_domain(&candidate, &context, 2).unwrap();
        assert!(named(&values, "smallest_q3") > 0.0);
        assert_eq!(named(&values, "smallest_q3"), named(&values, "largest_q3"));
    }

    #[test]
    fn every_emitted_value_is_finite_and_schema_ordered() {
        let (global, candidate) = extract_fixture_features().unwrap();
        assert_eq!(global.to_vec().len(), GLOBAL_FEATURE_NAMES.len());
        assert_eq!(
            candidate.base_and_domain_vec().len(),
            BASE_AND_DOMAIN_FEATURE_NAMES.len()
        );
        assert_eq!(candidate.values.len(), CANDIDATE_FEATURE_NAMES.len());
        assert!(candidate.values[BASE_AND_DOMAIN_FEATURE_NAMES.len()..]
            .iter()
            .all(|value| *value == 0.0));
        assert!(global
            .to_vec()
            .iter()
            .chain(candidate.base_and_domain_vec().iter())
            .all(|value| value.is_finite()));
    }

    #[test]
    fn masked_candidate_extraction_retains_base_and_zeros_domain_slots() {
        let context = fixture_context();
        let candidate = fixture_candidate();
        let full = extract_candidate_base_and_domain(&candidate, &context, 2).unwrap();
        let masked =
            extract_candidate_base_and_domain_with_mask(&candidate, &context, 2, base_only_mask())
                .unwrap();

        assert_eq!(masked.values.len(), CANDIDATE_FEATURE_NAMES.len());
        assert_eq!(
            &masked.values[..BASE_CANDIDATE_FEATURE_NAMES.len()],
            &full.values[..BASE_CANDIDATE_FEATURE_NAMES.len()]
        );
        assert!(masked.values
            [BASE_CANDIDATE_FEATURE_NAMES.len()..BASE_AND_DOMAIN_FEATURE_NAMES.len()]
            .iter()
            .all(|value| *value == 0.0));
        assert!(masked.values[BASE_AND_DOMAIN_FEATURE_NAMES.len()..]
            .iter()
            .all(|value| *value == 0.0));
        assert!(full.values
            [BASE_CANDIDATE_FEATURE_NAMES.len()..BASE_AND_DOMAIN_FEATURE_NAMES.len()]
            .iter()
            .any(|value| *value != 0.0));
    }

    #[test]
    fn conditional_population_writes_only_fixed_boundary_and_discontinuity_slots() {
        let context = fixture_context();
        let candidate = fixture_candidate();
        let mut boundary = extract_candidate_base_and_domain(&candidate, &context, 2).unwrap();
        let base = boundary.values[..BASE_AND_DOMAIN_FEATURE_NAMES.len()].to_vec();
        populate_candidate_conditional_features(
            &candidate,
            &mut boundary,
            &context,
            FeatureMask {
                boundary_local: true,
                ..base_only_mask()
            },
        )
        .unwrap();
        assert_eq!(boundary.values.len(), CANDIDATE_FEATURE_NAMES.len());
        assert_eq!(CANDIDATE_FEATURE_NAMES.len(), 154);
        assert_eq!(
            &boundary.values[..BASE_AND_DOMAIN_FEATURE_NAMES.len()],
            base.as_slice()
        );
        assert!(boundary.values[BOUNDARY_LOCAL_START..BOUNDARY_LOCAL_END]
            .iter()
            .any(|value| *value != 0.0));
        assert!(boundary.values[BOUNDARY_LOCAL_END..DISCONTINUITY_START]
            .iter()
            .all(|value| *value == 0.0));
        assert!(boundary.values[DISCONTINUITY_START..DISCONTINUITY_END]
            .iter()
            .all(|value| *value == 0.0));

        let discontinuous = CandidateRecord {
            source_index: 8,
            measure: MeasureLine {
                num_domains: 2,
                min_size: 6,
                delineation: "0-3;10-13 4-9".into(),
                max_cr: 0.25,
                mean_cr: 0.0,
                density_min: 1.5,
                mean_density: 2.5,
            },
            partition: parse_partition("0-3;10-13 4-9", 14).unwrap(),
            legacy_distance: 0.0,
            hierarchy: None,
        };
        let mut segment = extract_candidate_base_and_domain(&discontinuous, &context, 2).unwrap();
        populate_candidate_conditional_features(
            &discontinuous,
            &mut segment,
            &context,
            FeatureMask {
                discontinuity: true,
                ..base_only_mask()
            },
        )
        .unwrap();
        assert!(segment.values[BOUNDARY_LOCAL_START..BOUNDARY_LOCAL_END]
            .iter()
            .all(|value| *value == 0.0));
        assert!(segment.values[BOUNDARY_LOCAL_END..DISCONTINUITY_START]
            .iter()
            .all(|value| *value == 0.0));
        assert_eq!(segment.values[DISCONTINUITY_START], 1.0);
        assert_eq!(DISCONTINUITY_END, segment.values.len());
    }

    #[test]
    fn disabled_conditional_families_do_not_inspect_their_malformed_evidence() {
        let candidate = fixture_candidate();
        let mut boundary_malformed_context = fixture_context();
        let mut boundary_malformed_dssp = dssp_chain("HHHCCEEECCCCCC");
        boundary_malformed_dssp.get_mut(1).kappa = f64::NAN;
        boundary_malformed_dssp.get_mut(1).acceptor[0].residue = usize::MAX;
        boundary_malformed_context.dssp = Box::leak(Box::new(boundary_malformed_dssp));
        let mut boundary_masked =
            extract_candidate_base_and_domain(&candidate, &boundary_malformed_context, 2).unwrap();
        populate_candidate_conditional_features(
            &candidate,
            &mut boundary_masked,
            &boundary_malformed_context,
            FeatureMask {
                discontinuity: true,
                ..base_only_mask()
            },
        )
        .unwrap();
        assert!(
            boundary_masked.values[BOUNDARY_LOCAL_START..BOUNDARY_LOCAL_END]
                .iter()
                .all(|value| *value == 0.0)
        );

        let mut discontinuity_malformed_context = fixture_context();
        let mut discontinuity_malformed_dssp = dssp_chain("HHHCCEEECCCCCC");
        discontinuity_malformed_dssp.get_mut(1).partner[0] = usize::MAX;
        discontinuity_malformed_context.dssp = Box::leak(Box::new(discontinuity_malformed_dssp));
        let mut discontinuity_masked =
            extract_candidate_base_and_domain(&candidate, &discontinuity_malformed_context, 2)
                .unwrap();
        populate_candidate_conditional_features(
            &candidate,
            &mut discontinuity_masked,
            &discontinuity_malformed_context,
            base_only_mask(),
        )
        .unwrap();
        assert!(
            discontinuity_masked.values[DISCONTINUITY_START..DISCONTINUITY_END]
                .iter()
                .all(|value| *value == 0.0)
        );
    }

    #[test]
    fn conditional_population_validates_candidate_identity_and_vector_length() {
        let context = fixture_context();
        let candidate = fixture_candidate();
        let mut wrong_identity =
            extract_candidate_base_and_domain(&candidate, &context, 2).unwrap();
        wrong_identity.source_index += 1;
        assert_eq!(
            populate_candidate_conditional_features(
                &candidate,
                &mut wrong_identity,
                &context,
                FeatureMask::all(),
            ),
            Err(FeatureError::SchemaMismatch)
        );

        let mut wrong_length = extract_candidate_base_and_domain(&candidate, &context, 2).unwrap();
        wrong_length.values.pop();
        assert_eq!(
            populate_candidate_conditional_features(
                &candidate,
                &mut wrong_length,
                &context,
                FeatureMask::all(),
            ),
            Err(FeatureError::SchemaMismatch)
        );
    }

    #[test]
    fn conditional_population_rolls_back_when_the_second_extractor_fails() {
        let candidate = fixture_candidate();
        let mut context = fixture_context();
        let mut invalid_boundary_dssp = dssp_chain("HHHCCEEECCCCCC");
        invalid_boundary_dssp.get_mut(1).kappa = f64::NAN;
        context.dssp = Box::leak(Box::new(invalid_boundary_dssp));
        assert!(
            crate::sword::factorized_ranker::discontinuity::extract_discontinuity_features(
                &candidate.partition,
                &context,
            )
            .is_ok()
        );

        let mut features = extract_candidate_base_and_domain(&candidate, &context, 2).unwrap();
        for (index, value) in features.values[BOUNDARY_LOCAL_START..DISCONTINUITY_END]
            .iter_mut()
            .enumerate()
        {
            *value = index as f64 + 0.25;
        }
        let before = features.values.clone();

        assert_eq!(
            populate_candidate_conditional_features(
                &candidate,
                &mut features,
                &context,
                FeatureMask {
                    boundary_local: true,
                    discontinuity: true,
                    ..base_only_mask()
                },
            ),
            Err(FeatureError::NonFinite("DSSP angle"))
        );
        assert_eq!(features.values, before);
    }

    #[test]
    fn disabled_domain_family_bypasses_domain_only_validation_and_overflow() {
        let candidate = fixture_candidate();
        assert_eq!(
            domain_conditioned_values(&candidate.partition, &[]),
            Err(FeatureError::SchemaMismatch)
        );

        let finite_shape = ShapeMeasures {
            q1: 1.0,
            q2: 1.0,
            q3: 1.0,
            volume_ratio: 1.0,
            relative_density: 1.0,
        };
        let mut measures = vec![
            DomainMeasures {
                size_fraction: 0.25,
                shape: finite_shape,
                internal_contact_density: 1.0,
                contact_order: 1.0,
                internal_contact_fraction: 1.0,
                conductance: 1.0,
            },
            DomainMeasures {
                size_fraction: 0.75,
                shape: finite_shape,
                internal_contact_density: 1.0,
                contact_order: 1.0,
                internal_contact_fraction: 1.0,
                conductance: 1.0,
            },
        ];
        measures[0].shape.relative_density = f64::MAX;
        measures[1].shape.relative_density = f64::MIN_POSITIVE;
        assert_eq!(
            domain_conditioned_values(&candidate.partition, &measures),
            Err(FeatureError::NonFinite("smallest_to_largest_density_ratio"))
        );

        let masked = extract_candidate_base_and_domain_with_mask(
            &candidate,
            &fixture_context(),
            2,
            base_only_mask(),
        )
        .unwrap();
        assert!(masked.values
            [BASE_CANDIDATE_FEATURE_NAMES.len()..BASE_AND_DOMAIN_FEATURE_NAMES.len()]
            .iter()
            .all(|value| *value == 0.0));
    }

    #[test]
    fn frozen_schema_lengths_and_family_order_are_exact() {
        assert_eq!(FEATURE_SCHEMA_VERSION, 1);
        assert_eq!(GLOBAL_FEATURE_NAMES.len(), 37);
        assert_eq!(COUNT_ITEM_FEATURE_NAMES.len(), 98);
        assert_eq!(BASE_CANDIDATE_FEATURE_NAMES.len(), 22);
        assert_eq!(DOMAIN_CONDITIONED_FEATURE_NAMES.len(), 27);
        assert_eq!(BOUNDARY_LOCAL_FEATURE_NAMES.len(), 42);
        assert_eq!(RELATIVE_HIERARCHY_FEATURE_NAMES.len(), 47);
        assert_eq!(DISCONTINUITY_FEATURE_NAMES.len(), 16);
        assert_eq!(CANDIDATE_FEATURE_NAMES.len(), 154);
        assert_eq!(
            &CANDIDATE_FEATURE_NAMES[..BASE_CANDIDATE_FEATURE_NAMES.len()],
            BASE_CANDIDATE_FEATURE_NAMES
        );
        assert_eq!(
            &CANDIDATE_FEATURE_NAMES
                [BASE_CANDIDATE_FEATURE_NAMES.len()..BASE_AND_DOMAIN_FEATURE_NAMES.len()],
            DOMAIN_CONDITIONED_FEATURE_NAMES
        );
        assert_eq!(
            &CANDIDATE_FEATURE_NAMES[BASE_AND_DOMAIN_FEATURE_NAMES.len()
                ..BASE_AND_DOMAIN_FEATURE_NAMES.len() + BOUNDARY_LOCAL_FEATURE_NAMES.len()],
            BOUNDARY_LOCAL_FEATURE_NAMES
        );
        let unique: BTreeSet<_> = GLOBAL_FEATURE_NAMES
            .iter()
            .chain(COUNT_ITEM_FEATURE_NAMES)
            .chain(CANDIDATE_FEATURE_NAMES)
            .collect();
        assert_eq!(
            unique.len(),
            GLOBAL_FEATURE_NAMES.len()
                + COUNT_ITEM_FEATURE_NAMES.len()
                + CANDIDATE_FEATURE_NAMES.len()
        );
    }

    #[test]
    fn pair_vectors_preserve_shared_diff_and_absolute_diff_order() {
        let (global, mut left) = extract_fixture_features().unwrap();
        let mut right = left.clone();
        left.values[0] = 4.0;
        right.values[0] = 2.0;
        let names = [
            "chain_n_residues",
            "diff__num_domains",
            "abs_diff__num_domains",
        ];

        assert_eq!(
            candidate_pair_vector(&names, &global, &left, &right).unwrap(),
            vec![14.0, 2.0, 2.0]
        );
        right.values.pop();
        assert_eq!(
            candidate_pair_vector(&names, &global, &left, &right),
            Err(FeatureError::SchemaMismatch)
        );
        right.values.push(0.0);
        left.values[0] = f64::MAX;
        right.values[0] = -f64::MAX;
        assert_eq!(
            candidate_pair_vector(&names, &global, &left, &right),
            Err(FeatureError::NonFinite("pair_input"))
        );
    }

    #[test]
    fn masked_validation_only_requires_retained_conditional_evidence() {
        let mut context = fixture_context();
        context.iterations = &[];
        context.measure_provenance = &[];
        let mut invalid_dssp = dssp_chain("HHHCCEEECCCCCC");
        invalid_dssp.get_mut(1).kappa = f64::NAN;
        context.dssp = Box::leak(Box::new(invalid_dssp));
        let base_only = FeatureMask {
            global_count: false,
            domain_conditioned: true,
            boundary_local: false,
            relative_hierarchy: false,
            discontinuity: false,
        };
        assert_eq!(context.validate(base_only), Ok(()));
        assert!(matches!(
            context.validate(FeatureMask {
                global_count: true,
                ..base_only
            }),
            Err(FeatureError::MissingContext("Peeling iterations"))
        ));
        assert!(matches!(
            context.validate(FeatureMask {
                boundary_local: true,
                ..base_only
            }),
            Err(FeatureError::NonFinite("DSSP angle"))
        ));
        context.iterations = Box::leak(Box::new([iteration(14, 5)]));
        assert!(matches!(
            context.validate(FeatureMask {
                relative_hierarchy: true,
                ..base_only
            }),
            Err(FeatureError::MissingContext("measure provenance"))
        ));
    }

    #[test]
    fn candidate_extraction_names_the_first_non_finite_schema_field() {
        let mut candidate = fixture_candidate();
        candidate.measure.max_cr = f64::INFINITY;
        assert!(matches!(
            extract_candidate_base_and_domain(&candidate, &fixture_context(), 2),
            Err(FeatureError::NonFinite("max_cr"))
        ));
    }

    #[test]
    fn global_huge_finite_coordinates_fail_closed_without_panicking() {
        let coords: Vec<[f64; 3]> = (0..8)
            .map(|index| {
                if index % 2 == 0 {
                    [1.0e308, 0.0, 0.0]
                } else {
                    [-1.0e308, 0.0, 0.0]
                }
            })
            .collect();
        let context = context(
            coords,
            vec![[0.0, 0.0, 0.0]; 8],
            0.0,
            1.0,
            "CCCCCCCC",
            vec![iteration(8, 2)],
            vec![],
        );

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            extract_global_features(&context, &[2])
        }));
        assert!(result.is_ok(), "global extraction panicked");
        assert!(matches!(
            result.unwrap(),
            Err(FeatureError::NonFinite("chain_rg_normalized"))
        ));
    }

    #[test]
    fn candidate_huge_finite_coordinates_fail_closed_without_panicking() {
        let coords: Vec<[f64; 3]> = (0..14)
            .map(|index| {
                if index % 2 == 0 {
                    [1.0e308, 0.0, 0.0]
                } else {
                    [-1.0e308, 0.0, 0.0]
                }
            })
            .collect();
        let context = context(
            coords,
            vec![[0.0, 0.0, 0.0]; 14],
            0.0,
            1.0,
            "CCCCCCCCCCCCCC",
            vec![iteration(14, 2)],
            vec![],
        );

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            extract_candidate_base_and_domain(&fixture_candidate(), &context, 2)
        }));
        assert!(result.is_ok(), "candidate extraction panicked");
        assert!(matches!(
            result.unwrap(),
            Err(FeatureError::NonFinite("domain_q1_mean"))
        ));
    }
}
