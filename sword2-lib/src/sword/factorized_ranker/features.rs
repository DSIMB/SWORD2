#![allow(dead_code)]

use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

use crate::peeling::contact_matrix::ContactMatrix;
use crate::sword::geometry_metrics::{
    gyration_tensor, principal_radii, radius_of_gyration, sym3x3_eigenvalues,
};

use super::boundary::extract_boundary_features;
use super::discontinuity::extract_discontinuity_features;
use super::lattice::{sequential_boundaries, CandidateLattice, CandidateRecord};
use super::partition::{Domain, FeatureError, ParsedPartition, Segment};
use super::schema::{
    CandidateFeatures, CountFeatures, FeatureMask, GlobalFeatures, BASE_AND_DOMAIN_FEATURE_NAMES,
    BASE_CANDIDATE_FEATURE_NAMES, BOUNDARY_LOCAL_END, BOUNDARY_LOCAL_FEATURE_NAMES,
    BOUNDARY_LOCAL_START, CANDIDATE_FEATURE_NAMES, COUNT_ITEM_FEATURE_NAMES,
    COUNT_SUMMARY_SOURCE_FEATURE_NAMES, DISCONTINUITY_END, DISCONTINUITY_FEATURE_NAMES,
    DISCONTINUITY_START, DOMAIN_CONDITIONED_FEATURE_NAMES, GLOBAL_FEATURE_NAMES,
    RELATIVE_CORE_FEATURE_NAMES, RELATIVE_HIERARCHY_FEATURE_NAMES,
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
    if !candidate.legacy_distance.is_finite() {
        return Err(FeatureError::NonFinite("legacy_distance"));
    }

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
        legacy_distance: candidate.legacy_distance,
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
        || !features.legacy_distance.is_finite()
        || !candidate.legacy_distance.is_finite()
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

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct CandidateIdentity {
    source_index: usize,
    num_domains: usize,
    canonical: String,
}

impl CandidateIdentity {
    fn from_record(candidate: &CandidateRecord) -> Self {
        Self {
            source_index: candidate.source_index,
            num_domains: candidate.measure.num_domains,
            canonical: candidate.partition.canonical.clone(),
        }
    }

    fn from_features(features: &CandidateFeatures) -> Self {
        Self {
            source_index: features.source_index,
            num_domains: features.num_domains,
            canonical: features.canonical.clone(),
        }
    }
}

pub(crate) fn add_sibling_and_hierarchy_features(
    lattice: &CandidateLattice,
    features: &mut [CandidateFeatures],
) -> Result<(), FeatureError> {
    let relative_start = BASE_AND_DOMAIN_FEATURE_NAMES
        .len()
        .checked_add(BOUNDARY_LOCAL_FEATURE_NAMES.len())
        .ok_or(FeatureError::SchemaMismatch)?;
    let relative_end = relative_start
        .checked_add(RELATIVE_HIERARCHY_FEATURE_NAMES.len())
        .ok_or(FeatureError::SchemaMismatch)?;
    if relative_start != BOUNDARY_LOCAL_END
        || relative_end != DISCONTINUITY_START
        || relative_end.checked_add(DISCONTINUITY_FEATURE_NAMES.len())
            != Some(CANDIDATE_FEATURE_NAMES.len())
        || RELATIVE_HIERARCHY_FEATURE_NAMES.len() != 2 * RELATIVE_CORE_FEATURE_NAMES.len() + 7
    {
        return Err(FeatureError::SchemaMismatch);
    }

    let mut lattice_by_identity = BTreeMap::new();
    let mut lattice_count_canonicals = BTreeSet::new();
    for candidate in &lattice.candidates {
        if candidate.measure.num_domains == 0 {
            return Err(FeatureError::SchemaMismatch);
        }
        if !candidate.legacy_distance.is_finite() {
            return Err(FeatureError::NonFinite("legacy_distance"));
        }
        if lattice_by_identity
            .insert(CandidateIdentity::from_record(candidate), candidate)
            .is_some()
        {
            return Err(FeatureError::SchemaMismatch);
        }
        if !lattice_count_canonicals.insert((
            candidate.measure.num_domains,
            candidate.partition.canonical.clone(),
        )) {
            return Err(FeatureError::SchemaMismatch);
        }
    }
    let mut feature_indices = BTreeMap::new();
    let mut feature_count_canonicals = BTreeSet::new();
    for (index, row) in features.iter().enumerate() {
        if row.num_domains == 0 || row.values.len() != CANDIDATE_FEATURE_NAMES.len() {
            return Err(FeatureError::SchemaMismatch);
        }
        if !row.legacy_distance.is_finite() {
            return Err(FeatureError::NonFinite("legacy_distance"));
        }
        if feature_indices
            .insert(CandidateIdentity::from_features(row), index)
            .is_some()
        {
            return Err(FeatureError::SchemaMismatch);
        }
        if !feature_count_canonicals.insert((row.num_domains, row.canonical.clone())) {
            return Err(FeatureError::SchemaMismatch);
        }
    }
    if lattice_by_identity.len() != feature_indices.len()
        || lattice_by_identity.keys().ne(feature_indices.keys())
    {
        return Err(FeatureError::SchemaMismatch);
    }

    let mut groups: BTreeMap<usize, Vec<CandidateIdentity>> = BTreeMap::new();
    let mut source_values: BTreeMap<CandidateIdentity, Vec<f64>> = BTreeMap::new();
    let mut boundary_sets: BTreeMap<CandidateIdentity, Vec<usize>> = BTreeMap::new();
    for (identity, candidate) in &lattice_by_identity {
        let row = &features[feature_indices[identity]];
        let mut values = Vec::with_capacity(RELATIVE_CORE_FEATURE_NAMES.len());
        for &name in RELATIVE_CORE_FEATURE_NAMES {
            let value = row
                .base_and_domain_value(name)
                .ok_or(FeatureError::SchemaMismatch)?;
            if !value.is_finite() {
                return Err(FeatureError::NonFinite(name));
            }
            values.push(value);
        }
        let hierarchy = candidate
            .hierarchy
            .as_ref()
            .ok_or(FeatureError::MissingContext("candidate hierarchy"))?;
        if !hierarchy.parent_merge_margin.is_finite() {
            return Err(FeatureError::NonFinite("hierarchy_parent_merge_margin"));
        }
        if !hierarchy.child_merge_margin.is_finite() {
            return Err(FeatureError::NonFinite("hierarchy_child_merge_margin"));
        }
        groups
            .entry(identity.num_domains)
            .or_default()
            .push(identity.clone());
        source_values.insert(identity.clone(), values);
        boundary_sets.insert(
            identity.clone(),
            sequential_boundaries(&candidate.partition),
        );
    }

    let mut updates = BTreeMap::new();
    for (identity, candidate) in &lattice_by_identity {
        let group = groups
            .get(&identity.num_domains)
            .ok_or(FeatureError::SchemaMismatch)?;
        let current = source_values
            .get(identity)
            .ok_or(FeatureError::SchemaMismatch)?;
        let mut values = Vec::with_capacity(RELATIVE_HIERARCHY_FEATURE_NAMES.len());
        for source_index in 0..RELATIVE_CORE_FEATURE_NAMES.len() {
            let value = current[source_index];
            let mut less = 0usize;
            let mut equal = 0usize;
            for sibling in group {
                let sibling_value = source_values[sibling][source_index];
                less += usize::from(sibling_value < value);
                equal += usize::from(sibling_value == value);
            }
            values.push((less as f64 + 0.5 * equal as f64) / group.len() as f64);
        }
        for source_index in 0..RELATIVE_CORE_FEATURE_NAMES.len() {
            let mut ordered = group
                .iter()
                .map(|sibling| source_values[sibling][source_index])
                .collect::<Vec<_>>();
            ordered.sort_by(f64::total_cmp);
            let middle = ordered.len() / 2;
            let median = if ordered.len() % 2 == 1 {
                ordered[middle]
            } else {
                (ordered[middle - 1] + ordered[middle]) / 2.0
            };
            values.push(current[source_index] - median);
        }

        let hierarchy = candidate
            .hierarchy
            .as_ref()
            .ok_or(FeatureError::MissingContext("candidate hierarchy"))?;
        values.extend([
            hierarchy.first_appearance_level as f64,
            hierarchy.persistence_levels as f64,
            hierarchy.parent_merge_margin,
            hierarchy.child_merge_margin,
        ]);

        let mut nearest: Option<(&CandidateIdentity, f64)> = None;
        for sibling in group.iter().filter(|sibling| *sibling != identity) {
            let distance =
                symmetric_boundary_distance(&boundary_sets[identity], &boundary_sets[sibling])?;
            if nearest.is_none_or(|(best, best_distance)| {
                distance.total_cmp(&best_distance).is_lt()
                    || (distance == best_distance && sibling.canonical < best.canonical)
            }) {
                nearest = Some((sibling, distance));
            }
        }
        if let Some((nearest, _)) = nearest {
            let sibling = &source_values[nearest];
            let max_cr_index = RELATIVE_CORE_FEATURE_NAMES
                .iter()
                .position(|name| *name == "max_cr")
                .ok_or(FeatureError::SchemaMismatch)?;
            let density_min_index = RELATIVE_CORE_FEATURE_NAMES
                .iter()
                .position(|name| *name == "density_min")
                .ok_or(FeatureError::SchemaMismatch)?;
            values.push(current[max_cr_index] - sibling[max_cr_index]);
            values.push(current[density_min_index] - sibling[density_min_index]);
        } else {
            values.extend([0.0, 0.0]);
        }
        values.push(hierarchy.hierarchy_path_count as f64);

        if values.len() != RELATIVE_HIERARCHY_FEATURE_NAMES.len() {
            return Err(FeatureError::SchemaMismatch);
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(FeatureError::NonFinite("sibling relative feature"));
        }
        updates.insert(identity.clone(), values);
    }

    for (identity, values) in updates {
        let index = feature_indices[&identity];
        features[index].values[relative_start..relative_end].copy_from_slice(&values);
    }
    features.sort_by(|left, right| {
        (left.num_domains, &left.canonical).cmp(&(right.num_domains, &right.canonical))
    });
    Ok(())
}

pub(crate) fn symmetric_boundary_distance(
    left: &[usize],
    right: &[usize],
) -> Result<f64, FeatureError> {
    if left.is_empty() && right.is_empty() {
        return Ok(0.0);
    }
    if left.is_empty() || right.is_empty() {
        return Err(FeatureError::SchemaMismatch);
    }
    let directed = |source: &[usize], target: &[usize]| {
        source
            .iter()
            .map(|value| {
                target
                    .iter()
                    .map(|other| value.abs_diff(*other) as f64)
                    .fold(f64::INFINITY, f64::min)
            })
            .sum::<f64>()
            / source.len() as f64
    };
    let distance = 0.5 * (directed(left, right) + directed(right, left));
    if !distance.is_finite() {
        return Err(FeatureError::NonFinite("sibling boundary distance"));
    }
    Ok(distance)
}

pub(crate) fn extract_count_features(
    global: &GlobalFeatures,
    candidates: &[CandidateFeatures],
) -> Result<Vec<CountFeatures>, FeatureError> {
    if candidates.is_empty() || COUNT_ITEM_FEATURE_NAMES.len() != 98 {
        return Err(FeatureError::SchemaMismatch);
    }
    let candidate_total = finite_nonnegative_integer(global.candidate_total)?;
    let available_count_total = finite_nonnegative_integer(global.available_count_total)?;
    if candidate_total != candidates.len() {
        return Err(FeatureError::SchemaMismatch);
    }

    let mut identities = BTreeSet::new();
    let mut count_canonicals = BTreeSet::new();
    let mut groups: BTreeMap<usize, Vec<&CandidateFeatures>> = BTreeMap::new();
    for candidate in candidates {
        if candidate.num_domains == 0 || candidate.values.len() != CANDIDATE_FEATURE_NAMES.len() {
            return Err(FeatureError::SchemaMismatch);
        }
        if !identities.insert(CandidateIdentity::from_features(candidate)) {
            return Err(FeatureError::SchemaMismatch);
        }
        if !count_canonicals.insert((candidate.num_domains, candidate.canonical.clone())) {
            return Err(FeatureError::SchemaMismatch);
        }
        if !candidate.legacy_distance.is_finite() {
            return Err(FeatureError::NonFinite("legacy_distance"));
        }
        for &name in &COUNT_SUMMARY_SOURCE_FEATURE_NAMES[1..] {
            let value = candidate
                .base_and_domain_value(name)
                .ok_or(FeatureError::SchemaMismatch)?;
            if !value.is_finite() {
                return Err(FeatureError::NonFinite(name));
            }
        }
        groups
            .entry(candidate.num_domains)
            .or_default()
            .push(candidate);
    }
    if available_count_total != groups.len() {
        return Err(FeatureError::SchemaMismatch);
    }

    let modal_count = groups
        .iter()
        .max_by_key(|(count, group)| (group.len(), std::cmp::Reverse(**count)))
        .map(|(count, _)| *count)
        .ok_or(FeatureError::SchemaMismatch)?;
    if !global.modal_count.is_finite() || global.modal_count != modal_count as f64 {
        return Err(FeatureError::SchemaMismatch);
    }
    let mut histogram = [0.0; 21];
    for candidate in candidates {
        histogram[if candidate.num_domains <= 20 {
            candidate.num_domains - 1
        } else {
            20
        }] += 1.0;
    }
    for value in &mut histogram {
        *value /= candidates.len() as f64;
    }
    if global
        .count_histogram
        .iter()
        .zip(histogram)
        .any(|(actual, expected)| !actual.is_finite() || *actual != expected)
    {
        return Err(FeatureError::SchemaMismatch);
    }

    let mut count_summaries = BTreeMap::new();
    for (&count, group) in &groups {
        let mut summaries = Vec::with_capacity(COUNT_SUMMARY_SOURCE_FEATURE_NAMES.len() * 3);
        for &name in COUNT_SUMMARY_SOURCE_FEATURE_NAMES {
            let source = group
                .iter()
                .map(|candidate| {
                    if name == "legacy_distance" {
                        Some(candidate.legacy_distance)
                    } else {
                        candidate.base_and_domain_value(name)
                    }
                    .ok_or(FeatureError::SchemaMismatch)
                })
                .collect::<Result<Vec<_>, _>>()?;
            let minimum = source.iter().copied().fold(f64::INFINITY, f64::min);
            let average = source.iter().sum::<f64>() / source.len() as f64;
            let maximum = source.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            summaries.extend([minimum, average, maximum]);
        }
        if summaries.iter().any(|value| !value.is_finite()) {
            return Err(FeatureError::NonFinite("count summary"));
        }
        count_summaries.insert(count, summaries);
    }

    let counts = groups.keys().copied().collect::<Vec<_>>();
    let mut result = Vec::with_capacity(counts.len());
    for (index, count) in counts.iter().copied().enumerate() {
        let lower = index.checked_sub(1).map(|lower| counts[lower]);
        let higher = counts.get(index + 1).copied();
        let summaries = &count_summaries[&count];
        let mut values = vec![
            count as f64,
            groups[&count].len() as f64,
            groups[&count].len() as f64 / candidates.len() as f64,
            count.abs_diff(modal_count) as f64,
            if lower.is_some() { 1.0 } else { 0.0 },
            if higher.is_some() { 1.0 } else { 0.0 },
            lower.map_or(0.0, |lower| count.abs_diff(lower) as f64),
            higher.map_or(0.0, |higher| count.abs_diff(higher) as f64),
        ];
        values.extend(summaries);
        if let Some(lower) = lower {
            values.extend(
                summaries
                    .iter()
                    .zip(&count_summaries[&lower])
                    .map(|(current, neighbor)| current - neighbor),
            );
        } else {
            values.resize(values.len() + summaries.len(), 0.0);
        }
        if let Some(higher) = higher {
            values.extend(
                summaries
                    .iter()
                    .zip(&count_summaries[&higher])
                    .map(|(current, neighbor)| current - neighbor),
            );
        } else {
            values.resize(values.len() + summaries.len(), 0.0);
        }
        if values.len() != COUNT_ITEM_FEATURE_NAMES.len() {
            return Err(FeatureError::SchemaMismatch);
        }
        validate_named(COUNT_ITEM_FEATURE_NAMES, &values)?;
        result.push(CountFeatures {
            num_domains: count,
            values,
        });
    }
    Ok(result)
}

fn finite_nonnegative_integer(value: f64) -> Result<usize, FeatureError> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > usize::MAX as f64 {
        return Err(FeatureError::SchemaMismatch);
    }
    Ok(value as usize)
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
        add_sibling_and_hierarchy_features, domain_conditioned_values,
        extract_candidate_base_and_domain, extract_candidate_base_and_domain_with_mask,
        extract_count_features, extract_global_features, populate_candidate_conditional_features,
        symmetric_boundary_distance, DomainMeasures, ShapeMeasures,
    };
    use crate::dssp::types::BackboneResidue;
    use crate::dssp::DsspChain;
    use crate::peeling::algorithm::IterationResult;
    use crate::peeling::contact_matrix::ContactMatrix;
    use crate::sword::compute_measure::{MeasureLine, MeasureProvenance};
    use crate::sword::factorized_ranker::lattice::{
        sequential_boundaries, CandidateLattice, CandidateRecord, HierarchyEvidence,
    };
    use crate::sword::factorized_ranker::partition::{parse_partition, FeatureError};
    use crate::sword::factorized_ranker::schema::{
        candidate_pair_vector, CandidateFeatures, FeatureMask, GlobalFeatures,
        BASE_AND_DOMAIN_FEATURE_NAMES, BASE_CANDIDATE_FEATURE_NAMES, BOUNDARY_LOCAL_END,
        BOUNDARY_LOCAL_FEATURE_NAMES, BOUNDARY_LOCAL_START, CANDIDATE_FEATURE_NAMES,
        COUNT_ITEM_FEATURE_NAMES, DISCONTINUITY_END, DISCONTINUITY_FEATURE_NAMES,
        DISCONTINUITY_START, DOMAIN_CONDITIONED_FEATURE_NAMES, FEATURE_SCHEMA_VERSION,
        GLOBAL_FEATURE_NAMES, RELATIVE_CORE_FEATURE_NAMES, RELATIVE_HIERARCHY_FEATURE_NAMES,
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

    fn lattice_candidate(
        source_index: usize,
        num_domains: usize,
        delineation: &str,
        min_size: usize,
        max_cr: f64,
        density_min: f64,
    ) -> CandidateRecord {
        CandidateRecord {
            source_index,
            measure: MeasureLine {
                num_domains,
                min_size,
                delineation: delineation.into(),
                max_cr,
                mean_cr: 0.0,
                density_min,
                mean_density: source_index as f64 + 1.5,
            },
            partition: parse_partition(delineation, 12).unwrap(),
            legacy_distance: source_index as f64 / 10.0 + 0.1,
            hierarchy: Some(HierarchyEvidence {
                first_appearance_level: source_index % 3,
                persistence_levels: source_index % 4 + 1,
                parent_merge_margin: source_index as f64 / 20.0 + 0.05,
                child_merge_margin: source_index as f64 / 25.0 + 0.04,
                hierarchy_path_count: source_index as u64 + 2,
            }),
        }
    }

    fn feature_row(candidate: &CandidateRecord, min_size: f64) -> CandidateFeatures {
        let mut values = (0..CANDIDATE_FEATURE_NAMES.len())
            .map(|index| candidate.source_index as f64 * 10.0 + index as f64 / 100.0)
            .collect::<Vec<_>>();
        values[0] = candidate.measure.num_domains as f64;
        values[1] = min_size;
        values[2] = candidate.measure.max_cr;
        values[3] = candidate.measure.density_min;
        CandidateFeatures {
            source_index: candidate.source_index,
            canonical: candidate.partition.canonical.clone(),
            num_domains: candidate.measure.num_domains,
            legacy_distance: candidate.legacy_distance,
            values,
        }
    }

    fn lattice(candidates: Vec<CandidateRecord>) -> CandidateLattice {
        let mut groups = std::collections::BTreeMap::new();
        for (index, candidate) in candidates.iter().enumerate() {
            groups
                .entry(candidate.measure.num_domains)
                .or_insert_with(Vec::new)
                .push(index);
        }
        CandidateLattice { candidates, groups }
    }

    fn value_bits(rows: &[CandidateFeatures]) -> Vec<Vec<u64>> {
        rows.iter()
            .map(|row| row.values.iter().map(|value| value.to_bits()).collect())
            .collect()
    }

    #[derive(serde::Deserialize)]
    struct SyntheticFixture {
        chain_length: usize,
        candidate_feature_names: Vec<String>,
        count_feature_names: Vec<String>,
        global_feature_names: Vec<String>,
        global_values: Vec<f64>,
        candidate_inputs: Vec<SyntheticCandidateInput>,
        candidate_rows: Vec<SyntheticCandidateRow>,
        count_rows: Vec<SyntheticCountRow>,
    }

    #[derive(serde::Deserialize)]
    struct SyntheticCandidateInput {
        source_index: usize,
        num_domains: usize,
        delineation: String,
        canonical_delineation: String,
        legacy_distance: f64,
        sequential_boundaries: Vec<usize>,
        hierarchy: SyntheticHierarchy,
        values_before: Vec<f64>,
    }

    #[derive(serde::Deserialize)]
    struct SyntheticHierarchy {
        first_appearance_level: usize,
        persistence_levels: usize,
        parent_merge_margin: f64,
        child_merge_margin: f64,
        hierarchy_path_count: u64,
    }

    #[derive(serde::Deserialize)]
    struct SyntheticCandidateRow {
        source_index: usize,
        num_domains: usize,
        canonical_delineation: String,
        legacy_distance: f64,
        values: Vec<f64>,
    }

    #[derive(serde::Deserialize)]
    struct SyntheticCountRow {
        num_domains: usize,
        values: Vec<f64>,
    }

    fn global_for(candidates: &[CandidateFeatures]) -> GlobalFeatures {
        let counts: Vec<usize> = candidates.iter().map(|row| row.num_domains).collect();
        let mut frequencies = std::collections::BTreeMap::new();
        let mut histogram = [0.0; 21];
        for count in counts.iter().copied() {
            *frequencies.entry(count).or_insert(0usize) += 1;
            histogram[if count <= 20 { count - 1 } else { 20 }] += 1.0;
        }
        for value in &mut histogram {
            *value /= counts.len() as f64;
        }
        let modal_count = frequencies
            .iter()
            .max_by_key(|(count, frequency)| (**frequency, std::cmp::Reverse(**count)))
            .map(|(count, _)| *count)
            .unwrap();
        GlobalFeatures {
            n_residues: 12.0,
            rg_normalized: 0.0,
            inertia_ratio_21: 0.0,
            inertia_ratio_31: 0.0,
            nonlocal_contact_density: 0.0,
            contact_order: 0.0,
            helix_fraction: 0.0,
            strand_fraction: 0.0,
            coil_fraction: 1.0,
            helix_blocks: 0.0,
            strand_blocks: 0.0,
            peeling_levels: 1.0,
            finest_pus: 7.0,
            candidate_total: candidates.len() as f64,
            available_count_total: frequencies.len() as f64,
            count_histogram: histogram,
            modal_count: modal_count as f64,
        }
    }

    #[test]
    fn sibling_population_is_key_joined_permutation_invariant_and_slot_scoped() {
        let candidates = vec![
            lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 10.0),
            lattice_candidate(1, 2, "0-3 4-11", 1, 4.0, 9.0),
            lattice_candidate(2, 2, "0-4 5-11", 3, 8.0, 7.0),
        ];
        let original_features = vec![
            feature_row(&candidates[2], 3.0),
            feature_row(&candidates[0], 1.0),
            feature_row(&candidates[1], 1.0),
        ];
        let before = original_features.clone();
        let mut first = original_features.clone();
        add_sibling_and_hierarchy_features(&lattice(candidates.clone()), &mut first).unwrap();

        let mut reversed_lattice = candidates;
        reversed_lattice.reverse();
        let mut second = original_features;
        second.reverse();
        add_sibling_and_hierarchy_features(&lattice(reversed_lattice), &mut second).unwrap();

        assert_eq!(
            first
                .iter()
                .map(|row| (&row.canonical, &row.values))
                .collect::<Vec<_>>(),
            second
                .iter()
                .map(|row| (&row.canonical, &row.values))
                .collect::<Vec<_>>()
        );
        assert_eq!(first.len(), 3);
        assert_eq!(first[0].values.len(), 154);
        for row in &first {
            let old = before
                .iter()
                .find(|candidate| candidate.source_index == row.source_index)
                .unwrap();
            assert_eq!(
                row.values[..BOUNDARY_LOCAL_END]
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                old.values[..BOUNDARY_LOCAL_END]
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>()
            );
            let relative_end = BOUNDARY_LOCAL_END + RELATIVE_HIERARCHY_FEATURE_NAMES.len();
            assert_eq!(
                row.values[relative_end..]
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                old.values[relative_end..]
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>()
            );
        }

        let percentiles = first
            .iter()
            .map(|row| row.values[BOUNDARY_LOCAL_END])
            .collect::<Vec<_>>();
        assert_eq!(percentiles, vec![1.0 / 3.0, 1.0 / 3.0, 5.0 / 6.0]);
        let last = BOUNDARY_LOCAL_END + RELATIVE_HIERARCHY_FEATURE_NAMES.len();
        let first_row = &first[0];
        assert_eq!(
            &first_row.values[last - 7..last],
            &[0.0, 1.0, 0.05, 0.04, -2.0, 1.0, 2.0]
        );
    }

    #[test]
    fn sibling_medians_nearest_ties_and_singletons_are_exact() {
        let candidates = vec![
            lattice_candidate(0, 2, "0-3 4-11", 1, 10.0, 3.0),
            lattice_candidate(1, 2, "0-2 3-11", 2, 4.0, 1.0),
            lattice_candidate(2, 2, "0-4 5-11", 8, 7.0, 2.0),
            lattice_candidate(3, 2, "0-5 6-11", 9, 12.0, 4.0),
            lattice_candidate(4, 7, "0 1 2 3-4 5-6 7-8 9-11", 1, 5.0, 5.0),
        ];
        let mut features = candidates
            .iter()
            .zip([1.0, 2.0, 8.0, 9.0, 1.0])
            .map(|(candidate, min_size)| feature_row(candidate, min_size))
            .collect::<Vec<_>>();

        add_sibling_and_hierarchy_features(&lattice(candidates), &mut features).unwrap();
        let median_delta = BOUNDARY_LOCAL_END + RELATIVE_CORE_FEATURE_NAMES.len();
        let two_domain = features
            .iter()
            .filter(|row| row.num_domains == 2)
            .map(|row| row.values[median_delta])
            .collect::<Vec<_>>();
        assert_eq!(two_domain, vec![-3.0, -4.0, 3.0, 4.0]);

        let middle = features.iter().find(|row| row.source_index == 0).unwrap();
        let sibling_delta_start = BOUNDARY_LOCAL_END + RELATIVE_HIERARCHY_FEATURE_NAMES.len() - 3;
        assert_eq!(middle.values[sibling_delta_start], 6.0);
        assert_eq!(middle.values[sibling_delta_start + 1], 2.0);

        let singleton = features.iter().find(|row| row.num_domains == 7).unwrap();
        assert_eq!(singleton.values[sibling_delta_start], 0.0);
        assert_eq!(singleton.values[sibling_delta_start + 1], 0.0);
        assert_eq!(singleton.values[BOUNDARY_LOCAL_END], 0.5);
    }

    #[test]
    fn sibling_median_constants_and_asymmetric_boundary_distance_are_exact() {
        assert_eq!(
            symmetric_boundary_distance(&[1], &[0, 2, 10]).unwrap(),
            0.5 * (1.0 + 11.0 / 3.0)
        );

        let candidates = vec![
            lattice_candidate(0, 2, "0-1 2-11", 1, 2.0, 10.0),
            lattice_candidate(1, 2, "0-2 3-11", 3, 4.0, 9.0),
            lattice_candidate(2, 3, "0-1 2-5 6-11", 1, 6.0, 8.0),
            lattice_candidate(3, 3, "0-2 3-6 7-11", 2, 8.0, 7.0),
            lattice_candidate(4, 3, "0-3 4-7 8-11", 100, 10.0, 6.0),
        ];
        let mut features = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();
        add_sibling_and_hierarchy_features(&lattice(candidates), &mut features).unwrap();
        let median_delta = BOUNDARY_LOCAL_END + RELATIVE_CORE_FEATURE_NAMES.len();
        assert_eq!(
            features
                .iter()
                .filter(|row| row.num_domains == 2)
                .map(|row| row.values[median_delta])
                .collect::<Vec<_>>(),
            vec![-1.0, 1.0]
        );
        assert_eq!(
            features
                .iter()
                .filter(|row| row.num_domains == 3)
                .map(|row| row.values[median_delta])
                .collect::<Vec<_>>(),
            vec![-1.0, 0.0, 98.0]
        );
    }

    #[test]
    fn nearest_sibling_deltas_read_feature_slots_not_measure_lines() {
        let candidates = vec![
            lattice_candidate(0, 2, "0-3 4-11", 1, -100.0, -200.0),
            lattice_candidate(1, 2, "0-2 3-11", 2, 500.0, 600.0),
            lattice_candidate(2, 2, "0-4 5-11", 3, 700.0, 800.0),
        ];
        let mut features = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();
        features[0].values[2] = 10.0;
        features[0].values[3] = 3.0;
        features[1].values[2] = 4.0;
        features[1].values[3] = 1.0;
        features[2].values[2] = 7.0;
        features[2].values[3] = 2.0;

        add_sibling_and_hierarchy_features(&lattice(candidates), &mut features).unwrap();
        let middle = features.iter().find(|row| row.source_index == 0).unwrap();
        let delta_start = BOUNDARY_LOCAL_END + RELATIVE_HIERARCHY_FEATURE_NAMES.len() - 3;
        assert_eq!(middle.values[delta_start], 6.0);
        assert_eq!(middle.values[delta_start + 1], 2.0);
    }

    #[test]
    fn sibling_population_rejects_identity_defects_transactionally() {
        let candidates = vec![
            lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 10.0),
            lattice_candidate(1, 2, "0-3 4-11", 2, 4.0, 9.0),
        ];
        let candidate_lattice = lattice(candidates.clone());
        let valid = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();

        let mut duplicate = vec![valid[0].clone(), valid[0].clone()];
        let before = value_bits(&duplicate);
        assert_eq!(
            add_sibling_and_hierarchy_features(&candidate_lattice, &mut duplicate),
            Err(FeatureError::SchemaMismatch)
        );
        assert_eq!(value_bits(&duplicate), before);

        let mut missing = vec![valid[0].clone()];
        let before = value_bits(&missing);
        assert_eq!(
            add_sibling_and_hierarchy_features(&candidate_lattice, &mut missing),
            Err(FeatureError::SchemaMismatch)
        );
        assert_eq!(value_bits(&missing), before);

        let mut mismatched = valid.clone();
        mismatched[1].canonical.push_str("-wrong");
        let before = value_bits(&mismatched);
        assert_eq!(
            add_sibling_and_hierarchy_features(&candidate_lattice, &mut mismatched),
            Err(FeatureError::SchemaMismatch)
        );
        assert_eq!(value_bits(&mismatched), before);

        let mut duplicate_lattice_candidates = candidates;
        duplicate_lattice_candidates[1] = duplicate_lattice_candidates[0].clone();
        let mut untouched = valid;
        let before = value_bits(&untouched);
        assert_eq!(
            add_sibling_and_hierarchy_features(
                &lattice(duplicate_lattice_candidates),
                &mut untouched,
            ),
            Err(FeatureError::SchemaMismatch)
        );
        assert_eq!(value_bits(&untouched), before);
    }

    #[test]
    fn duplicate_count_canonical_is_rejected_independent_of_input_order() {
        let first = lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 10.0);
        let mut second = first.clone();
        second.source_index = 1;
        second.legacy_distance = 0.2;
        let candidates = vec![first, second];
        let features = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();

        for reverse in [false, true] {
            let mut ordered_candidates = candidates.clone();
            let mut ordered_features = features.clone();
            if reverse {
                ordered_candidates.reverse();
                ordered_features.reverse();
            }
            let before = value_bits(&ordered_features);
            assert_eq!(
                add_sibling_and_hierarchy_features(
                    &lattice(ordered_candidates),
                    &mut ordered_features,
                ),
                Err(FeatureError::SchemaMismatch)
            );
            assert_eq!(value_bits(&ordered_features), before);
            assert!(matches!(
                extract_count_features(&global_for(&ordered_features), &ordered_features),
                Err(FeatureError::SchemaMismatch)
            ));
        }
    }

    #[test]
    fn sibling_population_rejects_missing_hierarchy_and_nonfinite_sources_transactionally() {
        let mut candidates = vec![
            lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 10.0),
            lattice_candidate(1, 2, "0-3 4-11", 2, 4.0, 9.0),
        ];
        let valid_features = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();
        candidates[1].hierarchy = None;
        let mut features = valid_features.clone();
        let before = value_bits(&features);
        assert_eq!(
            add_sibling_and_hierarchy_features(&lattice(candidates), &mut features),
            Err(FeatureError::MissingContext("candidate hierarchy"))
        );
        assert_eq!(value_bits(&features), before);

        let candidates = vec![
            lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 10.0),
            lattice_candidate(1, 2, "0-3 4-11", 2, 4.0, 9.0),
        ];
        let mut features = valid_features;
        features[1].values[1] = f64::NAN;
        let before = value_bits(&features);
        assert_eq!(
            add_sibling_and_hierarchy_features(&lattice(candidates), &mut features),
            Err(FeatureError::NonFinite("min_size"))
        );
        assert_eq!(value_bits(&features), before);
    }

    #[test]
    fn one_empty_nearest_boundary_set_is_rejected_transactionally() {
        let mut candidates = vec![
            lattice_candidate(0, 2, "0-2 3-11", 3, 2.0, 10.0),
            lattice_candidate(1, 2, "0-5 6-11", 6, 4.0, 9.0),
        ];
        candidates[0].partition.residue_to_domain.fill(0);
        let mut features = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();
        let before = value_bits(&features);

        assert_eq!(
            add_sibling_and_hierarchy_features(&lattice(candidates), &mut features),
            Err(FeatureError::SchemaMismatch)
        );
        assert_eq!(value_bits(&features), before);
    }

    #[test]
    fn count_features_use_exact_summaries_and_nonconsecutive_neighbors() {
        let candidates = vec![
            lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 8.0),
            lattice_candidate(1, 2, "0-3 4-11", 3, 6.0, 4.0),
            lattice_candidate(2, 4, "0-1 2-4 5-8 9-11", 2, 10.0, 2.0),
            lattice_candidate(3, 7, "0 1 2 3-4 5-6 7-8 9-11", 1, 14.0, 1.0),
        ];
        let mut features = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();
        features[0].legacy_distance = 0.5;
        features[1].legacy_distance = 1.5;
        features[2].legacy_distance = 2.5;
        features[3].legacy_distance = 3.5;
        let global = global_for(&features);

        let counts = extract_count_features(&global, &features).unwrap();
        assert_eq!(counts.len(), 3);
        assert!(counts.iter().all(|row| row.values.len() == 98));
        let two = counts.iter().find(|row| row.num_domains == 2).unwrap();
        let four = counts.iter().find(|row| row.num_domains == 4).unwrap();
        assert_eq!(&two.values[..8], &[2.0, 2.0, 0.5, 0.0, 0.0, 1.0, 0.0, 2.0]);
        assert_eq!(&two.values[8..11], &[0.5, 1.0, 1.5]);
        assert_eq!(&two.values[14..17], &[2.0, 4.0, 6.0]);
        assert_eq!(four.values[6], 2.0);
        assert_eq!(four.values[7], 3.0);
        assert_eq!(four.values[8 + 2 * 3 + 1], 10.0);
        assert_eq!(four.values[38 + 2 * 3 + 1], 6.0);
        assert_eq!(four.values.len(), COUNT_ITEM_FEATURE_NAMES.len());
        assert!(counts
            .iter()
            .flat_map(|row| &row.values)
            .all(|value| value.is_finite()));
    }

    #[test]
    fn count_features_reject_mixed_global_population_and_bad_candidates() {
        let candidates = vec![
            lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 8.0),
            lattice_candidate(1, 4, "0-1 2-4 5-8 9-11", 2, 10.0, 2.0),
        ];
        let features = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();

        let mut global = global_for(&features);
        global.candidate_total = 3.0;
        assert!(matches!(
            extract_count_features(&global, &features),
            Err(FeatureError::SchemaMismatch)
        ));
        let mut global = global_for(&features);
        global.modal_count = 4.0;
        assert!(matches!(
            extract_count_features(&global, &features),
            Err(FeatureError::SchemaMismatch)
        ));
        let mut global = global_for(&features);
        global.count_histogram[1] = 0.0;
        assert!(matches!(
            extract_count_features(&global, &features),
            Err(FeatureError::SchemaMismatch)
        ));

        let mut duplicate = features.clone();
        duplicate[1] = duplicate[0].clone();
        assert!(matches!(
            extract_count_features(&global_for(&duplicate), &duplicate),
            Err(FeatureError::SchemaMismatch)
        ));
        let mut nonfinite = features;
        nonfinite[1].legacy_distance = f64::INFINITY;
        assert!(matches!(
            extract_count_features(&global_for(&nonfinite), &nonfinite),
            Err(FeatureError::NonFinite("legacy_distance"))
        ));
    }

    #[test]
    fn every_global_population_field_and_histogram_bin_is_validated() {
        let candidates = vec![
            lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 8.0),
            lattice_candidate(1, 4, "0-1 2-4 5-8 9-11", 2, 10.0, 2.0),
        ];
        let features = candidates
            .iter()
            .map(|candidate| feature_row(candidate, candidate.measure.min_size as f64))
            .collect::<Vec<_>>();

        for field in 0..3 {
            let mut global = global_for(&features);
            match field {
                0 => global.candidate_total += 1.0,
                1 => global.available_count_total += 1.0,
                _ => global.modal_count += 1.0,
            }
            assert!(matches!(
                extract_count_features(&global, &features),
                Err(FeatureError::SchemaMismatch)
            ));
        }
        for bin in 0..21 {
            let mut global = global_for(&features);
            global.count_histogram[bin] += 0.125;
            assert!(
                matches!(
                    extract_count_features(&global, &features),
                    Err(FeatureError::SchemaMismatch)
                ),
                "histogram bin {bin}"
            );
        }
        for invalid in [-1.0, 0.5, f64::NAN, f64::INFINITY] {
            let mut global = global_for(&features);
            global.candidate_total = invalid;
            assert!(extract_count_features(&global, &features).is_err());
            let mut global = global_for(&features);
            global.available_count_total = invalid;
            assert!(extract_count_features(&global, &features).is_err());
            let mut global = global_for(&features);
            global.modal_count = invalid;
            assert!(extract_count_features(&global, &features).is_err());
        }
    }

    #[test]
    fn count_summary_values_are_source_major_with_boundary_coil_fraction_last() {
        let candidate = lattice_candidate(0, 2, "0-2 3-11", 1, 2.0, 8.0);
        let mut features = feature_row(&candidate, 11.0);
        features.legacy_distance = 10.0;
        for (name, value) in [
            ("max_cr", 12.0),
            ("density_min", 13.0),
            ("mean_density", 14.0),
            ("contact_q_mean", 15.0),
            ("contact_q_max", 16.0),
            ("n_segments", 17.0),
            ("n_discontinuous", 18.0),
            ("boundary_coil_fraction", 19.0),
        ] {
            let index = CANDIDATE_FEATURE_NAMES
                .iter()
                .position(|candidate_name| *candidate_name == name)
                .unwrap();
            features.values[index] = value;
        }
        let [count] = extract_count_features(&global_for(&[features.clone()]), &[features])
            .unwrap()
            .try_into()
            .unwrap();
        let expected = (10..=19)
            .flat_map(|value| std::iter::repeat_n(value as f64, 3))
            .collect::<Vec<_>>();
        assert_eq!(&count.values[8..38], expected);
        assert_eq!(
            &COUNT_ITEM_FEATURE_NAMES[35..38],
            &[
                "count_boundary_coil_fraction_min",
                "count_boundary_coil_fraction_mean",
                "count_boundary_coil_fraction_max",
            ]
        );
    }

    #[test]
    fn legacy_distance_is_finite_metadata_not_a_candidate_value_or_pair_feature() {
        let mut candidate = fixture_candidate();
        candidate.legacy_distance = f64::NAN;
        assert!(matches!(
            extract_candidate_base_and_domain(&candidate, &fixture_context(), 2),
            Err(FeatureError::NonFinite("legacy_distance"))
        ));

        let (global, left) = extract_fixture_features().unwrap();
        let mut right = left.clone();
        right.legacy_distance = left.legacy_distance + 123.0;
        assert_eq!(left.values, right.values);
        let names = [
            "chain_n_residues",
            "diff__num_domains",
            "abs_diff__num_domains",
        ];
        assert_eq!(
            candidate_pair_vector(&names, &global, &left, &right).unwrap(),
            candidate_pair_vector(&names, &global, &left, &left).unwrap()
        );
        assert!(!CANDIDATE_FEATURE_NAMES.contains(&"legacy_distance"));
        assert_eq!(left.values.len(), 154);
    }

    #[test]
    fn synthetic_python_fixture_matches_rust_candidate_and_count_vectors() {
        let fixture: SyntheticFixture = serde_json::from_str(include_str!(
            "../../../../benchmark/fixtures/factorized_features_synthetic.json"
        ))
        .unwrap();
        assert_eq!(
            fixture.candidate_feature_names,
            CANDIDATE_FEATURE_NAMES
                .iter()
                .map(|name| name.to_string())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            fixture.count_feature_names,
            COUNT_ITEM_FEATURE_NAMES
                .iter()
                .map(|name| name.to_string())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            fixture.global_feature_names,
            GLOBAL_FEATURE_NAMES
                .iter()
                .map(|name| name.to_string())
                .collect::<Vec<_>>()
        );
        assert_eq!(fixture.global_values.len(), GLOBAL_FEATURE_NAMES.len());
        let global_value = |name: &str| {
            fixture.global_values[GLOBAL_FEATURE_NAMES
                .iter()
                .position(|candidate| *candidate == name)
                .unwrap()]
        };
        let global = GlobalFeatures {
            n_residues: global_value("chain_n_residues"),
            rg_normalized: global_value("chain_rg_normalized"),
            inertia_ratio_21: global_value("chain_inertia_ratio_21"),
            inertia_ratio_31: global_value("chain_inertia_ratio_31"),
            nonlocal_contact_density: global_value("chain_nonlocal_contact_density"),
            contact_order: global_value("chain_contact_order"),
            helix_fraction: global_value("chain_helix_fraction"),
            strand_fraction: global_value("chain_strand_fraction"),
            coil_fraction: global_value("chain_coil_fraction"),
            helix_blocks: global_value("chain_helix_blocks"),
            strand_blocks: global_value("chain_strand_blocks"),
            peeling_levels: global_value("chain_peeling_levels"),
            finest_pus: global_value("chain_finest_pus"),
            candidate_total: global_value("chain_candidate_total"),
            available_count_total: global_value("chain_available_count_total"),
            count_histogram: std::array::from_fn(|index| {
                if index < 20 {
                    global_value(&format!("chain_count_hist_{}", index + 1))
                } else {
                    global_value("chain_count_hist_21_plus")
                }
            }),
            modal_count: global_value("chain_modal_count"),
        };

        let mut candidates = Vec::new();
        let mut features = Vec::new();
        for input in &fixture.candidate_inputs {
            let partition = parse_partition(&input.delineation, fixture.chain_length).unwrap();
            assert_eq!(partition.canonical, input.canonical_delineation);
            assert_eq!(
                sequential_boundaries(&partition),
                input.sequential_boundaries
            );
            assert_eq!(input.values_before.len(), 154);
            let record = CandidateRecord {
                source_index: input.source_index,
                measure: MeasureLine {
                    num_domains: input.num_domains,
                    min_size: input.values_before[1] as usize,
                    delineation: input.delineation.clone(),
                    max_cr: input.values_before[2],
                    mean_cr: 0.0,
                    density_min: input.values_before[3],
                    mean_density: input.values_before[4],
                },
                partition,
                legacy_distance: input.legacy_distance,
                hierarchy: Some(HierarchyEvidence {
                    first_appearance_level: input.hierarchy.first_appearance_level,
                    persistence_levels: input.hierarchy.persistence_levels,
                    parent_merge_margin: input.hierarchy.parent_merge_margin,
                    child_merge_margin: input.hierarchy.child_merge_margin,
                    hierarchy_path_count: input.hierarchy.hierarchy_path_count,
                }),
            };
            features.push(CandidateFeatures {
                source_index: input.source_index,
                canonical: input.canonical_delineation.clone(),
                num_domains: input.num_domains,
                legacy_distance: input.legacy_distance,
                values: input.values_before.clone(),
            });
            candidates.push(record);
        }
        candidates.reverse();
        features.reverse();
        add_sibling_and_hierarchy_features(&lattice(candidates), &mut features).unwrap();
        assert_eq!(features.len(), fixture.candidate_rows.len());
        for (actual, expected) in features.iter().zip(&fixture.candidate_rows) {
            assert_eq!(actual.source_index, expected.source_index);
            assert_eq!(actual.num_domains, expected.num_domains);
            assert_eq!(actual.canonical, expected.canonical_delineation);
            assert_eq!(actual.legacy_distance, expected.legacy_distance);
            assert_eq!(actual.values.len(), expected.values.len());
            for (name, (actual, expected)) in CANDIDATE_FEATURE_NAMES
                .iter()
                .zip(actual.values.iter().zip(&expected.values))
            {
                assert!(
                    (actual - expected).abs() <= 1e-12,
                    "candidate feature {name}: {actual} != {expected}"
                );
            }
        }

        let counts = extract_count_features(&global, &features).unwrap();
        assert_eq!(counts.len(), fixture.count_rows.len());
        for (actual, expected) in counts.iter().zip(&fixture.count_rows) {
            assert_eq!(actual.num_domains, expected.num_domains);
            assert_eq!(actual.values.len(), expected.values.len());
            for (name, (actual, expected)) in COUNT_ITEM_FEATURE_NAMES
                .iter()
                .zip(actual.values.iter().zip(&expected.values))
            {
                assert!(
                    (actual - expected).abs() <= 1e-12,
                    "count feature {name}: {actual} != {expected}"
                );
            }
        }
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
