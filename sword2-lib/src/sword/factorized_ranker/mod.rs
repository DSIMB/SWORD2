pub(crate) mod lattice;
pub(crate) mod partition;

use crate::dssp::DsspChain;
use crate::peeling::algorithm::IterationResult;
use crate::peeling::contact_matrix::ContactMatrix;
use crate::sword::compute_measure::MeasureProvenance;

use self::partition::FeatureError;

#[allow(dead_code)]
pub(crate) struct StructuralContext<'a> {
    pub ca_coords: &'a [[f64; 3]],
    pub dssp: &'a DsspChain,
    pub contacts: &'a ContactMatrix,
    pub iterations: &'a [IterationResult],
    pub measure_provenance: &'a [MeasureProvenance],
    pub dssp_index_for_residue: Vec<usize>,
}

impl StructuralContext<'_> {
    #[allow(dead_code)]
    pub(crate) fn validate(&self) -> Result<(), FeatureError> {
        let mapped_dssp_indices: Vec<usize> = (1..=self.dssp.len)
            .filter(|&index| self.dssp.get(index).aa != '!')
            .collect();
        let chain_len = self.ca_coords.len();

        if mapped_dssp_indices.len() != chain_len || self.contacts.len() != chain_len {
            return Err(FeatureError::MissingContext("DSSP/chain length mismatch"));
        }
        if self.dssp_index_for_residue != mapped_dssp_indices {
            return Err(FeatureError::MissingContext(
                "DSSP residue mapping mismatch",
            ));
        }
        if self.iterations.is_empty() {
            return Err(FeatureError::MissingContext("Peeling iterations"));
        }
        if self
            .ca_coords
            .iter()
            .flatten()
            .any(|coordinate| !coordinate.is_finite())
        {
            return Err(FeatureError::NonFinite("CA coordinates"));
        }
        if (0..chain_len)
            .any(|row| (0..chain_len).any(|column| !self.contacts.get(row, column).is_finite()))
        {
            return Err(FeatureError::NonFinite("contact probability"));
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
    };
    context.validate()?;
    Ok(context)
}

#[cfg(test)]
mod tests {
    use super::{prepare_factorized_context, StructuralContext};
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
        }
    }

    #[test]
    fn context_rejects_dimension_mismatch() {
        let context = synthetic_context_with_lengths(8, 7, 8);
        assert!(matches!(
            context.validate(),
            Err(FeatureError::MissingContext("DSSP/chain length mismatch"))
        ));
    }

    #[test]
    fn unavailable_typed_cache_requests_whole_chain_fallback() {
        let result = prepare_factorized_context(None, None, &[], None);
        assert!(matches!(result, Err(FeatureError::MissingContext(_))));
    }
}
