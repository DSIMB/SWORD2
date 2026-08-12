pub(crate) mod boundary;
pub(crate) mod discontinuity;
pub(crate) mod features;
pub(crate) mod generated_model;
pub(crate) mod lattice;
pub(crate) mod model;
pub(crate) mod partition;
pub(crate) mod schema;

use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use crate::dssp::DsspChain;
use crate::peeling::algorithm::IterationResult;
use crate::peeling::contact_matrix::ContactMatrix;
use crate::sword::compute_measure::MeasureProvenance;

use self::features::ContactFeatureCache;
use self::generated_model::{
    CANDIDATE_FEATURE_NAMES as MODEL_CANDIDATE_FEATURE_NAMES, CANDIDATE_MODEL,
    COUNT_FEATURE_NAMES as MODEL_COUNT_FEATURE_NAMES, COUNT_MODEL, MODEL_INPUT_DTYPE,
    MODEL_THRESHOLD_POLICY, RETAINED_FEATURE_FAMILIES,
};
use self::lattice::CandidateLattice;
use self::model::{normalized_borda, predict_probability_validated, validate_model, ModelError};
use self::partition::FeatureError;
use self::schema::{
    candidate_pair_vector, count_pair_vector, CandidateFeatures, CountFeatures, FeatureMask,
    GlobalFeatures, BASE_CANDIDATE_FEATURE_NAMES, BOUNDARY_LOCAL_FEATURE_NAMES,
    CANDIDATE_FEATURE_NAMES, COUNT_ITEM_FEATURE_NAMES, DISCONTINUITY_FEATURE_NAMES,
    DOMAIN_CONDITIONED_FEATURE_NAMES, GLOBAL_FEATURE_NAMES, RELATIVE_HIERARCHY_FEATURE_NAMES,
};

#[derive(Debug, thiserror::Error)]
pub(crate) enum FactorizedError {
    #[error(transparent)]
    Feature(#[from] FeatureError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error("factorized feature schema mismatch")]
    SchemaMismatch,
    #[error("factorized candidate identity mismatch")]
    IdentityMismatch,
    #[error("factorized count-group mismatch")]
    CountGroupMismatch,
}

const FAMILY_ORDER: &[&str] = &[
    "base",
    "global_count",
    "domain_conditioned",
    "boundary_local",
    "relative_hierarchy",
    "discontinuity",
];

const COUNT_BASE_SOURCES: &[&str] = &[
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "contact_q_mean",
    "contact_q_max",
    "n_segments",
    "n_discontinuous",
    "boundary_coil_fraction",
];

#[derive(Debug)]
struct HeadSpec {
    shared: Vec<&'static str>,
    items: Vec<&'static str>,
    pair_names: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FactorizedSelection {
    pub measure_index: usize,
    pub num_domains: usize,
    pub canonical: String,
}

fn validate_retained_families(retained: &[&str]) -> Result<(), ModelError> {
    if retained.first() != Some(&"base") {
        return Err(ModelError::FeatureSchema {
            head: "retained_families",
        });
    }
    let mut previous = None;
    for family in retained {
        let position = FAMILY_ORDER
            .iter()
            .position(|candidate| candidate == family)
            .ok_or(ModelError::FeatureSchema {
                head: "retained_families",
            })?;
        if previous.is_some_and(|previous| position <= previous) {
            return Err(ModelError::FeatureSchema {
                head: "retained_families",
            });
        }
        previous = Some(position);
    }
    Ok(())
}

fn named_count_item(name: &str) -> Result<&'static str, ModelError> {
    COUNT_ITEM_FEATURE_NAMES
        .iter()
        .copied()
        .find(|candidate| *candidate == name)
        .ok_or(ModelError::FeatureSchema { head: "count" })
}

fn count_base_items() -> Result<Vec<&'static str>, ModelError> {
    let mut items = [
        "count_num_domains",
        "count_n_candidates",
        "count_candidate_fraction",
        "count_modal_distance",
    ]
    .into_iter()
    .map(named_count_item)
    .collect::<Result<Vec<_>, _>>()?;
    for source in COUNT_BASE_SOURCES {
        for summary in ["min", "mean", "max"] {
            items.push(named_count_item(&format!("count_{source}_{summary}"))?);
        }
    }
    Ok(items)
}

fn pair_names(shared: &[&str], items: &[&str]) -> Vec<String> {
    shared
        .iter()
        .map(|name| (*name).to_string())
        .chain(items.iter().map(|name| format!("diff__{name}")))
        .chain(items.iter().map(|name| format!("abs_diff__{name}")))
        .collect()
}

fn derive_head_specs(retained: &[&str]) -> Result<(HeadSpec, HeadSpec), ModelError> {
    validate_retained_families(retained)?;
    let has = |family| retained.contains(&family);

    let count_shared = if has("global_count") {
        GLOBAL_FEATURE_NAMES.to_vec()
    } else {
        Vec::new()
    };
    let count_items = if has("global_count") {
        COUNT_ITEM_FEATURE_NAMES.to_vec()
    } else {
        count_base_items()?
    };

    let candidate_shared = if has("global_count") {
        GLOBAL_FEATURE_NAMES.to_vec()
    } else {
        Vec::new()
    };
    let mut candidate_items = BASE_CANDIDATE_FEATURE_NAMES.to_vec();
    for (family, names) in [
        ("domain_conditioned", DOMAIN_CONDITIONED_FEATURE_NAMES),
        ("boundary_local", BOUNDARY_LOCAL_FEATURE_NAMES),
        ("relative_hierarchy", RELATIVE_HIERARCHY_FEATURE_NAMES),
        ("discontinuity", DISCONTINUITY_FEATURE_NAMES),
    ] {
        if has(family) {
            candidate_items.extend_from_slice(names);
        }
    }

    Ok((
        HeadSpec {
            pair_names: pair_names(&count_shared, &count_items),
            shared: count_shared,
            items: count_items,
        },
        HeadSpec {
            pair_names: pair_names(&candidate_shared, &candidate_items),
            shared: candidate_shared,
            items: candidate_items,
        },
    ))
}

fn exact_name_match(expected: &[String], actual: &[&str]) -> bool {
    expected.len() == actual.len()
        && expected
            .iter()
            .zip(actual)
            .all(|(expected, actual)| expected == actual)
}

fn family_for_candidate_item(name: &str) -> Result<Option<&'static str>, ModelError> {
    let mut matches = Vec::new();
    for (family, names) in [
        ("base", BASE_CANDIDATE_FEATURE_NAMES),
        ("domain_conditioned", DOMAIN_CONDITIONED_FEATURE_NAMES),
        ("boundary_local", BOUNDARY_LOCAL_FEATURE_NAMES),
        ("relative_hierarchy", RELATIVE_HIERARCHY_FEATURE_NAMES),
        ("discontinuity", DISCONTINUITY_FEATURE_NAMES),
    ] {
        if names.contains(&name) {
            matches.push(family);
        }
    }
    match matches.as_slice() {
        [] => Ok(None),
        [family] => Ok(Some(*family)),
        _ => Err(ModelError::FeatureSchema { head: "candidate" }),
    }
}

fn derive_mask_from_names(
    count_names: &[&str],
    candidate_names: &[&str],
    retained: &[&str],
) -> Result<FeatureMask, ModelError> {
    let mut mask = FeatureMask {
        global_count: false,
        domain_conditioned: false,
        boundary_local: false,
        relative_hierarchy: false,
        discontinuity: false,
    };
    for (head, names) in [("count", count_names), ("candidate", candidate_names)] {
        let mut seen = BTreeSet::new();
        for name in names {
            if !seen.insert(*name) {
                return Err(ModelError::FeatureSchema { head });
            }
            let item = name
                .strip_prefix("abs_diff__")
                .or_else(|| name.strip_prefix("diff__"));
            if let Some(item) = item {
                if item == "legacy_distance" {
                    return Err(ModelError::FeatureSchema { head });
                }
                if head == "count" {
                    if !COUNT_ITEM_FEATURE_NAMES.contains(&item) {
                        return Err(ModelError::FeatureSchema { head });
                    }
                } else {
                    match family_for_candidate_item(item)? {
                        Some("base") => {}
                        Some("domain_conditioned") => mask.domain_conditioned = true,
                        Some("boundary_local") => mask.boundary_local = true,
                        Some("relative_hierarchy") => mask.relative_hierarchy = true,
                        Some("discontinuity") => mask.discontinuity = true,
                        _ => return Err(ModelError::FeatureSchema { head }),
                    }
                }
            } else if GLOBAL_FEATURE_NAMES.contains(name) {
                mask.global_count = true;
            } else {
                return Err(ModelError::FeatureSchema { head });
            }
        }
    }
    let declared = FeatureMask {
        global_count: retained.contains(&"global_count"),
        domain_conditioned: retained.contains(&"domain_conditioned"),
        boundary_local: retained.contains(&"boundary_local"),
        relative_hierarchy: retained.contains(&"relative_hierarchy"),
        discontinuity: retained.contains(&"discontinuity"),
    };
    if !same_mask(mask, declared) {
        return Err(ModelError::FeatureSchema { head: "mask" });
    }
    Ok(mask)
}

fn same_mask(left: FeatureMask, right: FeatureMask) -> bool {
    left.global_count == right.global_count
        && left.domain_conditioned == right.domain_conditioned
        && left.boundary_local == right.boundary_local
        && left.relative_hierarchy == right.relative_hierarchy
        && left.discontinuity == right.discontinuity
}

fn validate_embedded_models_uncached() -> Result<FeatureMask, ModelError> {
    validate_policies(MODEL_INPUT_DTYPE, MODEL_THRESHOLD_POLICY)?;
    let (count, candidate) = derive_head_specs(RETAINED_FEATURE_FAMILIES)?;
    if !exact_name_match(&count.pair_names, &MODEL_COUNT_FEATURE_NAMES) {
        return Err(ModelError::FeatureSchema { head: "count" });
    }
    if !exact_name_match(&candidate.pair_names, &MODEL_CANDIDATE_FEATURE_NAMES) {
        return Err(ModelError::FeatureSchema { head: "candidate" });
    }
    validate_model(&COUNT_MODEL)?;
    validate_model(&CANDIDATE_MODEL)?;
    let tree_total = COUNT_MODEL.trees.len() + CANDIDATE_MODEL.trees.len();
    if tree_total > 192 {
        return Err(ModelError::TreeLimit {
            actual: tree_total,
            maximum: 192,
        });
    }
    let node_total = COUNT_MODEL.nodes.len() + CANDIDATE_MODEL.nodes.len();
    if node_total > 2_880 {
        return Err(ModelError::NodeLimit {
            actual: node_total,
            maximum: 2_880,
        });
    }
    derive_mask_from_names(
        &MODEL_COUNT_FEATURE_NAMES,
        &MODEL_CANDIDATE_FEATURE_NAMES,
        RETAINED_FEATURE_FAMILIES,
    )
}

fn validate_policies(
    input_dtype: &'static str,
    threshold_policy: &'static str,
) -> Result<(), ModelError> {
    if input_dtype != "float32" {
        return Err(ModelError::Policy {
            policy: "input_dtype",
            expected: "float32",
            actual: input_dtype,
        });
    }
    if threshold_policy != "floor_to_f32" {
        return Err(ModelError::Policy {
            policy: "threshold_policy",
            expected: "floor_to_f32",
            actual: threshold_policy,
        });
    }
    Ok(())
}

static EMBEDDED_FEATURE_MASK: OnceLock<FeatureMask> = OnceLock::new();

pub(crate) fn embedded_feature_mask() -> Result<FeatureMask, FactorizedError> {
    if let Some(mask) = EMBEDDED_FEATURE_MASK.get() {
        return Ok(*mask);
    }
    let mask = validate_embedded_models_uncached()?;
    let _ = EMBEDDED_FEATURE_MASK.set(mask);
    Ok(*EMBEDDED_FEATURE_MASK.get().unwrap_or(&mask))
}

fn same_float_vectors(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

fn candidate_unused_families_are_zero(values: &[f64], mask: FeatureMask) -> bool {
    let base_end = BASE_CANDIDATE_FEATURE_NAMES.len();
    let domain_end = base_end + DOMAIN_CONDITIONED_FEATURE_NAMES.len();
    let boundary_end = domain_end + BOUNDARY_LOCAL_FEATURE_NAMES.len();
    let relative_end = boundary_end + RELATIVE_HIERARCHY_FEATURE_NAMES.len();
    let discontinuity_end = relative_end + DISCONTINUITY_FEATURE_NAMES.len();
    let zero = |start: usize, end: usize| {
        values
            .get(start..end)
            .is_some_and(|slice| slice.iter().all(|value| value.to_bits() == 0))
    };
    (mask.domain_conditioned || zero(base_end, domain_end))
        && (mask.boundary_local || zero(domain_end, boundary_end))
        && (mask.relative_hierarchy || zero(boundary_end, relative_end))
        && (mask.discontinuity || zero(relative_end, discontinuity_end))
        && discontinuity_end == CANDIDATE_FEATURE_NAMES.len()
}

fn choose_count(scores: &[(usize, f64)], legacy_count: usize) -> Result<usize, ModelError> {
    let maximum = scores
        .iter()
        .map(|(_, score)| *score)
        .reduce(f64::max)
        .ok_or(ModelError::EmptyGroup)?;
    if scores
        .iter()
        .any(|(count, score)| *count == legacy_count && *score == maximum)
    {
        return Ok(legacy_count);
    }
    scores
        .iter()
        .filter(|(_, score)| *score == maximum)
        .map(|(count, _)| *count)
        .min()
        .ok_or(ModelError::EmptyGroup)
}

fn choose_candidate(scores: &[(String, f64)]) -> Result<String, ModelError> {
    let maximum = scores
        .iter()
        .map(|(_, score)| *score)
        .reduce(f64::max)
        .ok_or(ModelError::EmptyGroup)?;
    scores
        .iter()
        .filter(|(_, score)| *score == maximum)
        .map(|(canonical, _)| canonical.clone())
        .min()
        .ok_or(ModelError::EmptyGroup)
}

pub(crate) fn select_factorized(
    lattice: &CandidateLattice,
    global: &GlobalFeatures,
    counts: &[CountFeatures],
    candidates: &[CandidateFeatures],
    legacy_count: usize,
) -> Result<FactorizedSelection, FactorizedError> {
    let mask = embedded_feature_mask()?;
    let global_values = global.to_vec();
    if global_values.len() != GLOBAL_FEATURE_NAMES.len()
        || global_values.iter().any(|value| !value.is_finite())
    {
        return Err(FactorizedError::SchemaMismatch);
    }

    let mut lattice_keys = BTreeMap::new();
    let mut lattice_canonicals = BTreeSet::new();
    let mut grouped_indices = BTreeSet::new();
    let mut chain_len = None;
    if lattice.candidates.is_empty() || lattice.groups.is_empty() {
        return Err(ModelError::LatticeGroupMismatch.into());
    }
    for (&count, group) in &lattice.groups {
        if count == 0 || group.is_empty() || group.len() > 3 {
            return Err(ModelError::LatticeGroupMismatch.into());
        }
        for &index in group {
            let record = lattice
                .candidates
                .get(index)
                .ok_or(ModelError::LatticeGroupMismatch)?;
            let record_chain_len = record.partition.residue_to_domain.len();
            let reparsed =
                partition::parse_partition(&record.partition.canonical, record_chain_len)
                    .map_err(|_| ModelError::LatticeGroupMismatch)?;
            if !grouped_indices.insert(index)
                || record.measure.num_domains != count
                || record.partition.domains.len() != count
                || record.partition.canonical.is_empty()
                || reparsed != record.partition
                || chain_len.is_some_and(|expected| expected != record_chain_len)
            {
                return Err(ModelError::LatticeGroupMismatch.into());
            }
            chain_len.get_or_insert(record_chain_len);
        }
    }
    if grouped_indices.len() != lattice.candidates.len()
        || global.n_residues != chain_len.unwrap_or(0) as f64
    {
        return Err(ModelError::LatticeGroupMismatch.into());
    }
    for (index, record) in lattice.candidates.iter().enumerate() {
        let key = (
            record.source_index,
            record.measure.num_domains,
            record.partition.canonical.clone(),
        );
        if !lattice_canonicals.insert((
            record.measure.num_domains,
            record.partition.canonical.clone(),
        )) || lattice_keys.insert(key, index).is_some()
        {
            return Err(ModelError::LatticeIdentityMismatch.into());
        }
    }

    let mut feature_keys = BTreeMap::new();
    let mut feature_canonicals = BTreeSet::new();
    for (index, candidate) in candidates.iter().enumerate() {
        if candidate.num_domains == 0
            || candidate.canonical.is_empty()
            || !candidate.legacy_distance.is_finite()
            || candidate.values.len() != CANDIDATE_FEATURE_NAMES.len()
            || candidate.values.iter().any(|value| !value.is_finite())
            || candidate.values[0] != candidate.num_domains as f64
            || !candidate_unused_families_are_zero(&candidate.values, mask)
        {
            return Err(FactorizedError::SchemaMismatch);
        }
        let key = (
            candidate.source_index,
            candidate.num_domains,
            candidate.canonical.clone(),
        );
        if !feature_canonicals.insert((candidate.num_domains, candidate.canonical.clone()))
            || feature_keys.insert(key, index).is_some()
        {
            return Err(ModelError::LatticeIdentityMismatch.into());
        }
    }
    if lattice_keys.keys().ne(feature_keys.keys()) {
        return Err(ModelError::LatticeIdentityMismatch.into());
    }
    for (key, lattice_index) in &lattice_keys {
        let record = &lattice.candidates[*lattice_index];
        let candidate = &candidates[feature_keys[key]];
        let expected_base = [
            record.measure.num_domains as f64,
            record.measure.min_size as f64,
            record.measure.max_cr,
            record.measure.density_min,
            record.measure.mean_density,
        ];
        if record.legacy_distance.to_bits() != candidate.legacy_distance.to_bits()
            || !same_float_vectors(&candidate.values[..expected_base.len()], &expected_base)
        {
            return Err(ModelError::LatticeIdentityMismatch.into());
        }
    }

    let expected_counts = features::extract_count_features(global, candidates)?;
    let modal_count = global.modal_count as usize;
    if candidates.iter().any(|candidate| {
        candidate.values[6].to_bits()
            != (candidate.num_domains.abs_diff(modal_count) as f64).to_bits()
    }) {
        return Err(ModelError::LatticeIdentityMismatch.into());
    }
    let mut expected_count_map = BTreeMap::new();
    for count in &expected_counts {
        expected_count_map.insert(count.num_domains, count);
    }
    let mut count_map = BTreeMap::new();
    for count in counts {
        if count.num_domains == 0
            || count.values.len() != COUNT_ITEM_FEATURE_NAMES.len()
            || count.values.iter().any(|value| !value.is_finite())
            || count.values[0] != count.num_domains as f64
            || count_map.insert(count.num_domains, count).is_some()
        {
            return Err(FactorizedError::CountGroupMismatch);
        }
    }
    let lattice_counts = lattice.groups.keys().copied().collect::<BTreeSet<_>>();
    if lattice_counts != count_map.keys().copied().collect()
        || lattice_counts != expected_count_map.keys().copied().collect()
        || expected_count_map.iter().any(|(count, expected)| {
            count_map
                .get(count)
                .is_none_or(|actual| !same_float_vectors(&actual.values, &expected.values))
        })
    {
        return Err(FactorizedError::CountGroupMismatch);
    }

    let count_identities = count_map.keys().copied().collect::<Vec<_>>();
    let count_scores = normalized_borda(&count_identities, |left, right| {
        let pair = count_pair_vector(
            &MODEL_COUNT_FEATURE_NAMES,
            global,
            count_map[left],
            count_map[right],
        )
        .map_err(|_| ModelError::FeatureSchema { head: "count" })?;
        predict_probability_validated(&COUNT_MODEL, &pair)
    })?;
    let selected_count = choose_count(&count_scores, legacy_count)?;

    let mut selected_candidates = BTreeMap::new();
    for candidate in candidates
        .iter()
        .filter(|candidate| candidate.num_domains == selected_count)
    {
        if selected_candidates
            .insert(candidate.canonical.clone(), candidate)
            .is_some()
        {
            return Err(ModelError::LatticeIdentityMismatch.into());
        }
    }
    let candidate_identities = selected_candidates.keys().cloned().collect::<Vec<_>>();
    let candidate_scores = normalized_borda(&candidate_identities, |left, right| {
        let pair = candidate_pair_vector(
            &MODEL_CANDIDATE_FEATURE_NAMES,
            global,
            selected_candidates[left],
            selected_candidates[right],
        )
        .map_err(|_| ModelError::FeatureSchema { head: "candidate" })?;
        predict_probability_validated(&CANDIDATE_MODEL, &pair)
    })?;
    let canonical = choose_candidate(&candidate_scores)?;
    let winner = selected_candidates
        .get(&canonical)
        .ok_or(ModelError::LatticeIdentityMismatch)?;
    let key = (
        winner.source_index,
        winner.num_domains,
        winner.canonical.clone(),
    );
    let lattice_index = lattice_keys
        .get(&key)
        .ok_or(ModelError::LatticeIdentityMismatch)?;
    let record = &lattice.candidates[*lattice_index];
    Ok(FactorizedSelection {
        measure_index: record.source_index,
        num_domains: selected_count,
        canonical,
    })
}

fn dump_header() -> Vec<&'static str> {
    std::iter::once("chain_id")
        .chain(std::iter::once("canonical_delineation"))
        .chain(std::iter::once("source_index"))
        .chain(std::iter::once("legacy_distance"))
        .chain(GLOBAL_FEATURE_NAMES.iter().copied())
        .chain(COUNT_ITEM_FEATURE_NAMES.iter().copied())
        .chain(CANDIDATE_FEATURE_NAMES.iter().copied())
        .collect()
}

fn temporary_dump_path(path: &Path) -> PathBuf {
    let file_name = path
        .file_name()
        .map(|name| name.to_string_lossy())
        .unwrap_or_default();
    path.with_file_name(format!(".{file_name}.tmp-{}", std::process::id()))
}

fn atomically_write_dump<F>(path: &Path, write_rows: F) -> Result<(), FactorizedError>
where
    F: FnOnce(&mut BufWriter<std::fs::File>) -> Result<(), FactorizedError>,
{
    let temporary = temporary_dump_path(path);
    let result = (|| {
        let file = std::fs::File::create(&temporary)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "{}", dump_header().join(","))?;
        write_rows(&mut writer)?;
        writer.flush()?;
        writer.get_ref().sync_all()?;
        std::fs::rename(&temporary, path)?;
        Ok(())
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    result
}

pub(crate) fn write_empty_feature_dump(path: &Path) -> Result<(), FactorizedError> {
    atomically_write_dump(path, |_| Ok(()))
}

pub(crate) fn install_failure_dump_header(
    path: &Path,
    result: Result<(), FactorizedError>,
) -> Result<Option<FactorizedError>, FactorizedError> {
    match result {
        Ok(()) => Ok(None),
        Err(error) => {
            write_empty_feature_dump(path)?;
            Ok(Some(error))
        }
    }
}

fn csv_quoted(value: &str) -> String {
    format!("\"{}\"", value.replace('"', "\"\""))
}

pub(crate) fn write_feature_dump(
    path: &Path,
    chain_id: &str,
    global: &GlobalFeatures,
    counts: &[CountFeatures],
    candidates: &[CandidateFeatures],
) -> Result<(), FactorizedError> {
    let global_values = global.to_vec();
    if global_values.len() != GLOBAL_FEATURE_NAMES.len()
        || global_values.iter().any(|value| !value.is_finite())
        || candidates.is_empty()
    {
        return Err(FactorizedError::SchemaMismatch);
    }

    let mut identities = BTreeSet::new();
    for candidate in candidates {
        if candidate.num_domains == 0
            || candidate.values.len() != CANDIDATE_FEATURE_NAMES.len()
            || candidate.values.iter().any(|value| !value.is_finite())
            || !candidate.legacy_distance.is_finite()
            || candidate.values[0] != candidate.num_domains as f64
        {
            return Err(FactorizedError::SchemaMismatch);
        }
        if !identities.insert((
            candidate.source_index,
            candidate.num_domains,
            candidate.canonical.as_str(),
        )) {
            return Err(FactorizedError::IdentityMismatch);
        }
    }

    let expected_counts = features::extract_count_features(global, candidates)?;
    let mut count_lookup = BTreeMap::new();
    for count in counts {
        if count.num_domains == 0
            || count.values.len() != COUNT_ITEM_FEATURE_NAMES.len()
            || count.values.iter().any(|value| !value.is_finite())
            || count.values[0] != count.num_domains as f64
            || count_lookup.insert(count.num_domains, count).is_some()
        {
            return Err(FactorizedError::CountGroupMismatch);
        }
    }
    if count_lookup.len() != expected_counts.len()
        || expected_counts.iter().any(|expected| {
            count_lookup
                .get(&expected.num_domains)
                .is_none_or(|actual| !same_float_vectors(&actual.values, &expected.values))
        })
    {
        return Err(FactorizedError::CountGroupMismatch);
    }

    let mut ordered_candidates = candidates.iter().collect::<Vec<_>>();
    ordered_candidates.sort_by(|left, right| {
        (left.num_domains, left.canonical.as_str(), left.source_index).cmp(&(
            right.num_domains,
            right.canonical.as_str(),
            right.source_index,
        ))
    });
    atomically_write_dump(path, |writer| {
        for candidate in ordered_candidates {
            let count = count_lookup
                .get(&candidate.num_domains)
                .ok_or(FactorizedError::CountGroupMismatch)?;
            let mut fields = Vec::with_capacity(dump_header().len());
            fields.push(csv_quoted(chain_id));
            fields.push(csv_quoted(&candidate.canonical));
            fields.push(candidate.source_index.to_string());
            fields.push(candidate.legacy_distance.to_string());
            fields.extend(global_values.iter().map(ToString::to_string));
            fields.extend(count.values.iter().map(ToString::to_string));
            fields.extend(candidate.values.iter().map(ToString::to_string));
            writeln!(writer, "{}", fields.join(","))?;
        }
        Ok(())
    })
}

#[allow(dead_code)]
pub(crate) struct StructuralContext<'a> {
    pub ca_coords: &'a [[f64; 3]],
    pub dssp: &'a DsspChain,
    pub contacts: &'a ContactMatrix,
    pub iterations: &'a [IterationResult],
    pub measure_provenance: &'a [MeasureProvenance],
    pub dssp_index_for_residue: Vec<usize>,
    pub contact_feature_cache: OnceLock<ContactFeatureCache>,
}

impl StructuralContext<'_> {
    #[allow(dead_code)]
    pub(crate) fn validate(&self, mask: FeatureMask) -> Result<(), FeatureError> {
        let mapped_dssp_indices: Vec<usize> = (1..=self.dssp.len)
            .filter(|&index| self.dssp.get(index).aa != '!')
            .collect();
        let chain_len = self.ca_coords.len();

        if chain_len == 0 {
            return Err(FeatureError::MissingContext("empty chain"));
        }
        if mapped_dssp_indices.len() != chain_len || self.contacts.len() != chain_len {
            return Err(FeatureError::MissingContext("DSSP/chain length mismatch"));
        }
        if self.dssp_index_for_residue != mapped_dssp_indices {
            return Err(FeatureError::MissingContext(
                "DSSP residue mapping mismatch",
            ));
        }
        if (mask.global_count || mask.relative_hierarchy) && self.iterations.is_empty() {
            return Err(FeatureError::MissingContext("Peeling iterations"));
        }
        if mask.relative_hierarchy && self.measure_provenance.is_empty() {
            return Err(FeatureError::MissingContext("measure provenance"));
        }
        if self
            .ca_coords
            .iter()
            .flatten()
            .any(|coordinate| !coordinate.is_finite())
        {
            return Err(FeatureError::NonFinite("CA coordinates"));
        }
        if let Some(cache) = self.contact_feature_cache.get() {
            cache.validate()?;
        } else {
            for row in 0..chain_len {
                for column in 0..chain_len {
                    let probability = self.contacts.get(row, column);
                    if !probability.is_finite() {
                        return Err(FeatureError::NonFinite("contact probability"));
                    }
                    if !(0.0..=1.0).contains(&probability) {
                        return Err(FeatureError::MissingContext(
                            "contact probability outside [0, 1]",
                        ));
                    }
                    if (probability - self.contacts.get(column, row)).abs() > 1e-12 {
                        return Err(FeatureError::MissingContext("asymmetric contact matrix"));
                    }
                }
            }
        }
        if mask.boundary_local || mask.discontinuity {
            for &dssp_index in &self.dssp_index_for_residue {
                let residue = self.dssp.get(dssp_index);
                if mask.boundary_local && (!residue.kappa.is_finite() || !residue.alpha.is_finite())
                {
                    return Err(FeatureError::NonFinite("DSSP angle"));
                }
                if mask.boundary_local {
                    for bond in residue.acceptor.iter().chain(residue.donor.iter()) {
                        if bond.residue > self.dssp.len {
                            return Err(FeatureError::MissingContext("DSSP hydrogen bond partner"));
                        }
                    }
                }
                for partner in residue.partner {
                    if partner > self.dssp.len {
                        return Err(FeatureError::MissingContext("DSSP bridge partner"));
                    }
                }
            }
        }

        Ok(())
    }
}

#[allow(dead_code)]
pub(crate) fn prepare_factorized_context<'a>(
    ca_coords: Option<&'a [[f64; 3]]>,
    dssp: Option<&'a DsspChain>,
    iterations: &'a [IterationResult],
    typed_evidence: Option<(&'a ContactMatrix, &'a [MeasureProvenance])>,
    mask: FeatureMask,
) -> Result<StructuralContext<'a>, FeatureError> {
    let (ca_coords, dssp, (contacts, measure_provenance)) = match (ca_coords, dssp, typed_evidence)
    {
        (Some(ca_coords), Some(dssp), Some((contacts, measure_provenance))) => {
            (ca_coords, dssp, (contacts, measure_provenance))
        }
        _ => {
            return Err(FeatureError::MissingContext(
                "complete typed factorized context",
            ))
        }
    };
    let dssp_index_for_residue = (1..=dssp.len)
        .filter(|&index| dssp.get(index).aa != '!')
        .collect();
    let context = StructuralContext {
        ca_coords,
        dssp,
        contacts,
        iterations,
        measure_provenance,
        dssp_index_for_residue,
        contact_feature_cache: OnceLock::new(),
    };
    context
        .contact_feature_cache
        .get_or_init(|| ContactFeatureCache::new(context.contacts));
    context.validate(mask)?;
    Ok(context)
}

#[cfg(test)]
mod tests {
    use super::{
        choose_candidate, choose_count, derive_head_specs, derive_mask_from_names, dump_header,
        embedded_feature_mask, exact_name_match, install_failure_dump_header,
        prepare_factorized_context, select_factorized, validate_policies, write_empty_feature_dump,
        write_feature_dump, FeatureMask, StructuralContext,
    };
    use crate::dssp::DsspChain;
    use crate::peeling::algorithm::IterationResult;
    use crate::peeling::contact_matrix::ContactMatrix;
    use crate::sword::factorized_ranker::partition::FeatureError;

    #[test]
    fn embedded_schema_and_mask_are_derived_from_exact_model_names() {
        use crate::sword::factorized_ranker::generated_model::{
            CANDIDATE_FEATURE_NAMES, COUNT_FEATURE_NAMES, RETAINED_FEATURE_FAMILIES,
        };

        let (count, candidate) = derive_head_specs(RETAINED_FEATURE_FAMILIES).unwrap();
        assert_eq!(count.shared.len(), 37);
        assert_eq!(count.items.len(), 98);
        assert_eq!(candidate.shared.len(), 37);
        assert_eq!(candidate.items.len(), 111);
        assert!(exact_name_match(&count.pair_names, &COUNT_FEATURE_NAMES));
        assert!(exact_name_match(
            &candidate.pair_names,
            &CANDIDATE_FEATURE_NAMES
        ));

        let mask = embedded_feature_mask().unwrap();
        assert!(mask.global_count);
        assert!(!mask.domain_conditioned);
        assert!(mask.boundary_local);
        assert!(mask.relative_hierarchy);
        assert!(!mask.discontinuity);

        let (base_count, base_candidate) = derive_head_specs(&["base"]).unwrap();
        assert!(base_count.shared.is_empty());
        assert_eq!(base_count.items.len(), 31);
        assert!(!base_count.items.contains(&"count_legacy_distance_min"));
        assert!(base_candidate.shared.is_empty());
        assert_eq!(base_candidate.items, super::BASE_CANDIDATE_FEATURE_NAMES);
    }

    #[test]
    fn schema_validation_rejects_family_and_pair_name_drift() {
        use crate::sword::factorized_ranker::generated_model::{
            CANDIDATE_FEATURE_NAMES, COUNT_FEATURE_NAMES,
        };
        use crate::sword::factorized_ranker::model::ModelError;

        assert!(validate_policies("float32", "floor_to_f32").is_ok());
        assert!(matches!(
            validate_policies("float64", "floor_to_f32"),
            Err(ModelError::Policy {
                policy: "input_dtype",
                ..
            })
        ));
        assert!(matches!(
            validate_policies("float32", "raw_f64"),
            Err(ModelError::Policy {
                policy: "threshold_policy",
                ..
            })
        ));

        for retained in [
            vec!["global_count"],
            vec!["base", "base"],
            vec!["base", "boundary_local", "global_count"],
            vec!["base", "unknown"],
        ] {
            assert!(derive_head_specs(&retained).is_err());
        }

        let (count, _) = derive_head_specs(&["base", "global_count"]).unwrap();
        let mut missing = count.pair_names.clone();
        missing.pop();
        assert!(!exact_name_match(&missing, &COUNT_FEATURE_NAMES));
        let mut reordered = count.pair_names.clone();
        reordered.swap(0, 1);
        assert!(!exact_name_match(&reordered, &COUNT_FEATURE_NAMES));
        let mut extra = count.pair_names.clone();
        extra.push("diff__count_num_domains".to_string());
        assert!(!exact_name_match(&extra, &COUNT_FEATURE_NAMES));

        assert!(matches!(
            derive_mask_from_names(
                &["diff__count_num_domains", "diff__count_num_domains"],
                &CANDIDATE_FEATURE_NAMES,
                &[
                    "base",
                    "global_count",
                    "boundary_local",
                    "relative_hierarchy"
                ]
            ),
            Err(ModelError::FeatureSchema { head: "count" })
        ));
        assert!(derive_mask_from_names(
            &["diff__num_domains"],
            &["diff__legacy_distance"],
            &["base"]
        )
        .is_err());
        assert!(derive_mask_from_names(&[], &["diff__count_num_domains"], &["base"]).is_err());
    }

    #[test]
    fn exact_ties_use_frozen_count_and_candidate_rules() {
        assert_eq!(choose_count(&[(1, 0.5), (2, 0.5), (3, 0.5)], 2).unwrap(), 2);
        assert_eq!(choose_count(&[(1, 0.5), (2, 0.5)], 9).unwrap(), 1);
        assert_eq!(
            choose_count(&[(1, 0.5), (2, 0.5 + f64::EPSILON)], 1).unwrap(),
            2
        );
        assert_eq!(
            choose_candidate(&[("z".to_string(), 0.5), ("a".to_string(), 0.5)]).unwrap(),
            "a"
        );
        assert_eq!(
            choose_candidate(&[
                ("a".to_string(), 0.5),
                ("z".to_string(), 0.5 + f64::EPSILON),
            ])
            .unwrap(),
            "z"
        );
    }

    #[test]
    fn selection_joins_by_full_identity_and_returns_original_source_index() {
        use std::collections::BTreeMap;

        use crate::sword::compute_measure::MeasureLine;
        use crate::sword::factorized_ranker::features::extract_count_features;
        use crate::sword::factorized_ranker::lattice::{CandidateLattice, CandidateRecord};
        use crate::sword::factorized_ranker::model::ModelError;
        use crate::sword::factorized_ranker::partition::parse_partition;
        use crate::sword::factorized_ranker::schema::{
            CandidateFeatures, GlobalFeatures, BASE_CANDIDATE_FEATURE_NAMES,
            CANDIDATE_FEATURE_NAMES,
        };

        fn record(source_index: usize, delineation: &str, count: usize) -> CandidateRecord {
            CandidateRecord {
                source_index,
                measure: MeasureLine {
                    num_domains: count,
                    min_size: 0,
                    delineation: delineation.to_string(),
                    max_cr: 0.0,
                    mean_cr: 0.0,
                    density_min: 0.0,
                    mean_density: 0.0,
                },
                partition: parse_partition(delineation, 6).unwrap(),
                legacy_distance: 0.0,
                hierarchy: None,
            }
        }

        fn features(source_index: usize, canonical: &str, count: usize) -> CandidateFeatures {
            let mut values = vec![0.0; CANDIDATE_FEATURE_NAMES.len()];
            values[0] = count as f64;
            values[6] = count.abs_diff(2) as f64;
            CandidateFeatures {
                source_index,
                canonical: canonical.to_string(),
                num_domains: count,
                legacy_distance: 0.0,
                values,
            }
        }

        let lattice = CandidateLattice {
            candidates: vec![
                record(77, "0-2 3-5", 2),
                record(41, "0-5", 1),
                record(12, "0-1 2-5", 2),
            ],
            groups: BTreeMap::from([(1, vec![1]), (2, vec![0, 2])]),
        };
        let candidates = vec![
            features(12, "0-1 2-5", 2),
            features(41, "0-5", 1),
            features(77, "0-2 3-5", 2),
        ];
        let mut histogram = [0.0; 21];
        histogram[0] = 1.0 / 3.0;
        histogram[1] = 2.0 / 3.0;
        let global = GlobalFeatures {
            n_residues: 6.0,
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
            finest_pus: 1.0,
            candidate_total: 3.0,
            available_count_total: 2.0,
            count_histogram: histogram,
            modal_count: 2.0,
        };
        let counts = extract_count_features(&global, &candidates).unwrap();
        let selected = select_factorized(&lattice, &global, &counts, &candidates, 2).unwrap();
        let mut reversed_candidates = candidates.clone();
        reversed_candidates.reverse();
        let mut reversed_counts = counts.clone();
        reversed_counts.reverse();
        assert_eq!(
            select_factorized(&lattice, &global, &reversed_counts, &reversed_candidates, 2)
                .unwrap(),
            selected
        );
        match selected.num_domains {
            1 => {
                assert_eq!(selected.measure_index, 41);
                assert_eq!(selected.canonical, "0-5");
            }
            2 => {
                assert_eq!(selected.measure_index, 12);
                assert_eq!(selected.canonical, "0-1 2-5");
            }
            _ => panic!("selected unavailable count"),
        }

        let mut mismatched = candidates.clone();
        mismatched[0].source_index = 999;
        assert!(matches!(
            select_factorized(&lattice, &global, &counts, &mismatched, 2),
            Err(super::FactorizedError::Model(
                ModelError::LatticeIdentityMismatch
            ))
        ));

        let mut duplicate_candidates = candidates.clone();
        duplicate_candidates.push(candidates[0].clone());
        assert!(matches!(
            select_factorized(&lattice, &global, &counts, &duplicate_candidates, 2),
            Err(super::FactorizedError::Model(
                ModelError::LatticeIdentityMismatch
            ))
        ));
        assert!(matches!(
            select_factorized(&lattice, &global, &counts, &candidates[..2], 2),
            Err(super::FactorizedError::Model(
                ModelError::LatticeIdentityMismatch
            ))
        ));
        assert!(matches!(
            select_factorized(&lattice, &global, &counts[..1], &candidates, 2),
            Err(super::FactorizedError::CountGroupMismatch)
        ));

        let mut nonfinite_candidate = candidates.clone();
        nonfinite_candidate[0].values[10] = f64::NAN;
        assert!(matches!(
            select_factorized(&lattice, &global, &counts, &nonfinite_candidate, 2),
            Err(super::FactorizedError::SchemaMismatch)
        ));

        let mut wrong_modal_distance = candidates.clone();
        wrong_modal_distance[0].values[6] = 99.0;
        assert!(matches!(
            select_factorized(&lattice, &global, &counts, &wrong_modal_distance, 2),
            Err(super::FactorizedError::Model(
                ModelError::LatticeIdentityMismatch
            ))
        ));

        let mut unmasked = candidates.clone();
        unmasked[0].values[BASE_CANDIDATE_FEATURE_NAMES.len()] = 1.0;
        assert!(matches!(
            select_factorized(&lattice, &global, &counts, &unmasked, 2),
            Err(super::FactorizedError::SchemaMismatch)
        ));
        let mut bad_groups = lattice.clone();
        bad_groups.groups.get_mut(&2).unwrap().push(99);
        assert!(matches!(
            select_factorized(&bad_groups, &global, &counts, &candidates, 2),
            Err(super::FactorizedError::Model(
                ModelError::LatticeGroupMismatch
            ))
        ));

        let mut bad_partition = lattice.clone();
        bad_partition.candidates[0].partition.residue_to_domain[0] = usize::MAX;
        assert!(matches!(
            select_factorized(&bad_partition, &global, &counts, &candidates, 2),
            Err(super::FactorizedError::Model(
                ModelError::LatticeGroupMismatch
            ))
        ));

        let mut bad_legacy = lattice.clone();
        bad_legacy.candidates[0].legacy_distance = 1.0;
        assert!(matches!(
            select_factorized(&bad_legacy, &global, &counts, &candidates, 2),
            Err(super::FactorizedError::Model(
                ModelError::LatticeIdentityMismatch
            ))
        ));

        let mut wrong_length = global.clone();
        wrong_length.n_residues = 7.0;
        assert!(matches!(
            select_factorized(&lattice, &wrong_length, &counts, &candidates, 2),
            Err(super::FactorizedError::Model(
                ModelError::LatticeGroupMismatch
            ))
        ));
    }

    fn synthetic_context_with_lengths(
        ca_len: usize,
        dssp_len: usize,
        contact_len: usize,
    ) -> StructuralContext<'static> {
        let ca_coords = Box::leak(
            std::iter::repeat_n([0.0, 0.0, 0.0], ca_len)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        );
        let mut dssp = DsspChain::new();
        for _ in 0..dssp_len {
            dssp.push(Default::default());
        }
        let dssp = Box::leak(Box::new(dssp));
        let contacts = Box::leak(Box::new(ContactMatrix::from_ca_coords(
            &vec![[0.0, 0.0, 0.0]; contact_len],
            6.0,
            1.5,
        )));
        let iterations = Box::leak(Box::new([IterationResult {
            max_cr: 0.0,
            min_density: 0.0,
            ci: 0.0,
            r: 0.0,
            num_pus: 1,
            pu_boundaries: vec![[0, ca_len.saturating_sub(1)]],
        }]));
        StructuralContext {
            ca_coords,
            dssp,
            contacts,
            iterations,
            measure_provenance: &[],
            dssp_index_for_residue: Vec::new(),
            contact_feature_cache: std::sync::OnceLock::new(),
        }
    }

    #[test]
    fn context_rejects_dimension_mismatch() {
        let context = synthetic_context_with_lengths(8, 7, 8);
        assert!(matches!(
            context.validate(FeatureMask::all()),
            Err(FeatureError::MissingContext("DSSP/chain length mismatch"))
        ));
    }

    #[test]
    fn unavailable_typed_cache_requests_whole_chain_fallback() {
        let result = prepare_factorized_context(None, None, &[], None, FeatureMask::all());
        assert!(matches!(result, Err(FeatureError::MissingContext(_))));
    }

    #[test]
    fn preparation_uses_caller_feature_mask() {
        let ca_coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut dssp = DsspChain::new();
        dssp.push(Default::default());
        dssp.push(Default::default());
        let contacts = ContactMatrix::from_ca_coords(&ca_coords, 6.0, 1.5);
        let mask = FeatureMask {
            global_count: false,
            domain_conditioned: true,
            boundary_local: false,
            relative_hierarchy: false,
            discontinuity: false,
        };

        assert!(prepare_factorized_context(
            Some(&ca_coords),
            Some(&dssp),
            &[],
            Some((&contacts, &[])),
            mask,
        )
        .is_ok());
    }

    #[test]
    fn exact_dump_round_trips_quotes_and_rust_float_strings() {
        use crate::sword::factorized_ranker::features::extract_count_features;
        use crate::sword::factorized_ranker::schema::{
            CandidateFeatures, GlobalFeatures, CANDIDATE_FEATURE_NAMES,
        };

        let candidate = CandidateFeatures {
            source_index: 7,
            canonical: "0-2 \"3-5\"".to_string(),
            num_domains: 2,
            legacy_distance: -0.0,
            values: (0..CANDIDATE_FEATURE_NAMES.len())
                .map(|index| index as f64 + 0.25)
                .collect(),
        };
        let mut candidate = candidate;
        candidate.values[0] = 2.0;
        let mut histogram = [0.0; 21];
        histogram[1] = 1.0;
        let global = GlobalFeatures {
            n_residues: 6.0,
            rg_normalized: 0.1,
            inertia_ratio_21: 0.2,
            inertia_ratio_31: 0.3,
            nonlocal_contact_density: 0.4,
            contact_order: 0.5,
            helix_fraction: 0.6,
            strand_fraction: 0.2,
            coil_fraction: 0.2,
            helix_blocks: 1.0,
            strand_blocks: 1.0,
            peeling_levels: 2.0,
            finest_pus: 2.0,
            candidate_total: 1.0,
            available_count_total: 1.0,
            count_histogram: histogram,
            modal_count: 2.0,
        };
        let counts = extract_count_features(&global, std::slice::from_ref(&candidate)).unwrap();
        let path = std::env::temp_dir().join(format!(
            "sword2-feature-dump-{}-{}.csv",
            std::process::id(),
            candidate.source_index
        ));
        write_feature_dump(
            &path,
            "chain,\n\"id",
            &global,
            &counts,
            std::slice::from_ref(&candidate),
        )
        .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert!(text.contains("\"chain,\n\"\"id\""));
        assert!(text.contains("\"0-2 \"\"3-5\"\"\""));
        assert!(text.contains(&(-0.0_f64).to_string()));

        let protected = std::env::temp_dir().join(format!(
            "sword2-feature-protected-{}.csv",
            std::process::id()
        ));
        std::fs::write(&protected, "previous\n").unwrap();
        let mut wrong_counts = counts.clone();
        wrong_counts[0].values[1] += 1.0;
        assert!(matches!(
            write_feature_dump(
                &protected,
                "chain",
                &global,
                &wrong_counts,
                std::slice::from_ref(&candidate),
            ),
            Err(super::FactorizedError::CountGroupMismatch)
        ));
        assert_eq!(std::fs::read_to_string(&protected).unwrap(), "previous\n");

        let mut wrong_candidate = candidate.clone();
        wrong_candidate.values[0] = 3.0;
        assert!(matches!(
            write_feature_dump(&protected, "chain", &global, &counts, &[wrong_candidate]),
            Err(super::FactorizedError::SchemaMismatch)
        ));
        assert_eq!(std::fs::read_to_string(&protected).unwrap(), "previous\n");
        std::fs::remove_file(protected).unwrap();
    }

    #[test]
    fn failure_artifact_is_exact_header_only() {
        use crate::sword::factorized_ranker::schema::{
            CANDIDATE_FEATURE_NAMES, COUNT_ITEM_FEATURE_NAMES, GLOBAL_FEATURE_NAMES,
        };

        let path =
            std::env::temp_dir().join(format!("sword2-feature-header-{}.csv", std::process::id()));
        write_empty_feature_dump(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let expected = std::iter::once("chain_id")
            .chain(std::iter::once("canonical_delineation"))
            .chain(std::iter::once("source_index"))
            .chain(std::iter::once("legacy_distance"))
            .chain(GLOBAL_FEATURE_NAMES.iter().copied())
            .chain(COUNT_ITEM_FEATURE_NAMES.iter().copied())
            .chain(CANDIDATE_FEATURE_NAMES.iter().copied())
            .collect::<Vec<_>>()
            .join(",")
            + "\n";
        assert_eq!(text, expected);
        assert_eq!(dump_header().len(), 293);
        assert_eq!(text.matches("num_domains").count(), 2); // frozen candidate + count name
        assert_eq!(text.matches(",num_domains").count(), 1);

        std::fs::write(&path, "partial,row\n").unwrap();
        let caught =
            install_failure_dump_header(&path, Err(super::FactorizedError::SchemaMismatch))
                .unwrap();
        assert!(matches!(
            caught,
            Some(super::FactorizedError::SchemaMismatch)
        ));
        let replaced = std::fs::read_to_string(&path).unwrap();
        assert_eq!(replaced, expected);
        std::fs::remove_file(&path).unwrap();

        let missing_parent = std::env::temp_dir()
            .join(format!("sword2-missing-parent-{}", std::process::id()))
            .join("dump.csv");
        assert!(matches!(
            write_empty_feature_dump(&missing_parent),
            Err(super::FactorizedError::Io(_))
        ));
    }
}
