pub(crate) mod boundary;
pub(crate) mod discontinuity;
pub(crate) mod features;
pub(crate) mod lattice;
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
use self::partition::FeatureError;
use self::schema::{
    CandidateFeatures, CountFeatures, FeatureMask, GlobalFeatures, CANDIDATE_FEATURE_NAMES,
    COUNT_ITEM_FEATURE_NAMES, GLOBAL_FEATURE_NAMES,
};

#[derive(Debug, thiserror::Error)]
pub(crate) enum FactorizedError {
    #[error(transparent)]
    Feature(#[from] FeatureError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("factorized feature schema mismatch")]
    SchemaMismatch,
    #[error("factorized candidate identity mismatch")]
    IdentityMismatch,
    #[error("factorized count-group mismatch")]
    CountGroupMismatch,
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

fn same_float_vectors(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| left.to_bits() == right.to_bits())
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
        dump_header, install_failure_dump_header, prepare_factorized_context,
        write_empty_feature_dump, write_feature_dump, FeatureMask, StructuralContext,
    };
    use crate::dssp::DsspChain;
    use crate::peeling::algorithm::IterationResult;
    use crate::peeling::contact_matrix::ContactMatrix;
    use crate::sword::factorized_ranker::partition::FeatureError;

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
