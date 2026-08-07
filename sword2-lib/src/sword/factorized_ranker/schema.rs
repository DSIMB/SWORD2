#![allow(dead_code)]

use super::partition::FeatureError;

pub(crate) const FEATURE_SCHEMA_VERSION: u32 = 1;

pub(crate) const GLOBAL_FEATURE_NAMES: &[&str] = &[
    "chain_n_residues",
    "chain_rg_normalized",
    "chain_inertia_ratio_21",
    "chain_inertia_ratio_31",
    "chain_nonlocal_contact_density",
    "chain_contact_order",
    "chain_helix_fraction",
    "chain_strand_fraction",
    "chain_coil_fraction",
    "chain_helix_blocks",
    "chain_strand_blocks",
    "chain_peeling_levels",
    "chain_finest_pus",
    "chain_candidate_total",
    "chain_available_count_total",
    "chain_count_hist_1",
    "chain_count_hist_2",
    "chain_count_hist_3",
    "chain_count_hist_4",
    "chain_count_hist_5",
    "chain_count_hist_6",
    "chain_count_hist_7",
    "chain_count_hist_8",
    "chain_count_hist_9",
    "chain_count_hist_10",
    "chain_count_hist_11",
    "chain_count_hist_12",
    "chain_count_hist_13",
    "chain_count_hist_14",
    "chain_count_hist_15",
    "chain_count_hist_16",
    "chain_count_hist_17",
    "chain_count_hist_18",
    "chain_count_hist_19",
    "chain_count_hist_20",
    "chain_count_hist_21_plus",
    "chain_modal_count",
];

pub(crate) const COUNT_ITEM_FEATURE_NAMES: &[&str] = &[
    "count_num_domains",
    "count_n_candidates",
    "count_candidate_fraction",
    "count_modal_distance",
    "count_has_lower",
    "count_has_higher",
    "count_lower_gap",
    "count_higher_gap",
    "count_legacy_distance_min",
    "count_legacy_distance_mean",
    "count_legacy_distance_max",
    "count_min_size_min",
    "count_min_size_mean",
    "count_min_size_max",
    "count_max_cr_min",
    "count_max_cr_mean",
    "count_max_cr_max",
    "count_density_min_min",
    "count_density_min_mean",
    "count_density_min_max",
    "count_mean_density_min",
    "count_mean_density_mean",
    "count_mean_density_max",
    "count_contact_q_mean_min",
    "count_contact_q_mean_mean",
    "count_contact_q_mean_max",
    "count_contact_q_max_min",
    "count_contact_q_max_mean",
    "count_contact_q_max_max",
    "count_n_segments_min",
    "count_n_segments_mean",
    "count_n_segments_max",
    "count_n_discontinuous_min",
    "count_n_discontinuous_mean",
    "count_n_discontinuous_max",
    "count_boundary_coil_fraction_min",
    "count_boundary_coil_fraction_mean",
    "count_boundary_coil_fraction_max",
    "count_legacy_distance_min_delta_lower",
    "count_legacy_distance_mean_delta_lower",
    "count_legacy_distance_max_delta_lower",
    "count_min_size_min_delta_lower",
    "count_min_size_mean_delta_lower",
    "count_min_size_max_delta_lower",
    "count_max_cr_min_delta_lower",
    "count_max_cr_mean_delta_lower",
    "count_max_cr_max_delta_lower",
    "count_density_min_min_delta_lower",
    "count_density_min_mean_delta_lower",
    "count_density_min_max_delta_lower",
    "count_mean_density_min_delta_lower",
    "count_mean_density_mean_delta_lower",
    "count_mean_density_max_delta_lower",
    "count_contact_q_mean_min_delta_lower",
    "count_contact_q_mean_mean_delta_lower",
    "count_contact_q_mean_max_delta_lower",
    "count_contact_q_max_min_delta_lower",
    "count_contact_q_max_mean_delta_lower",
    "count_contact_q_max_max_delta_lower",
    "count_n_segments_min_delta_lower",
    "count_n_segments_mean_delta_lower",
    "count_n_segments_max_delta_lower",
    "count_n_discontinuous_min_delta_lower",
    "count_n_discontinuous_mean_delta_lower",
    "count_n_discontinuous_max_delta_lower",
    "count_boundary_coil_fraction_min_delta_lower",
    "count_boundary_coil_fraction_mean_delta_lower",
    "count_boundary_coil_fraction_max_delta_lower",
    "count_legacy_distance_min_delta_higher",
    "count_legacy_distance_mean_delta_higher",
    "count_legacy_distance_max_delta_higher",
    "count_min_size_min_delta_higher",
    "count_min_size_mean_delta_higher",
    "count_min_size_max_delta_higher",
    "count_max_cr_min_delta_higher",
    "count_max_cr_mean_delta_higher",
    "count_max_cr_max_delta_higher",
    "count_density_min_min_delta_higher",
    "count_density_min_mean_delta_higher",
    "count_density_min_max_delta_higher",
    "count_mean_density_min_delta_higher",
    "count_mean_density_mean_delta_higher",
    "count_mean_density_max_delta_higher",
    "count_contact_q_mean_min_delta_higher",
    "count_contact_q_mean_mean_delta_higher",
    "count_contact_q_mean_max_delta_higher",
    "count_contact_q_max_min_delta_higher",
    "count_contact_q_max_mean_delta_higher",
    "count_contact_q_max_max_delta_higher",
    "count_n_segments_min_delta_higher",
    "count_n_segments_mean_delta_higher",
    "count_n_segments_max_delta_higher",
    "count_n_discontinuous_min_delta_higher",
    "count_n_discontinuous_mean_delta_higher",
    "count_n_discontinuous_max_delta_higher",
    "count_boundary_coil_fraction_min_delta_higher",
    "count_boundary_coil_fraction_mean_delta_higher",
    "count_boundary_coil_fraction_max_delta_higher",
];

pub(crate) const BASE_CANDIDATE_FEATURE_NAMES: &[&str] = &[
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "modal_count_distance",
    "domain_q1_mean",
    "domain_q2_mean",
    "domain_q3_mean",
    "domain_q3_min",
    "domain_vol_ratio_mean",
    "domain_density_mean",
    "domain_density_min",
    "contact_q_mean",
    "contact_q_max",
    "n_segments",
    "n_discontinuous",
    "size_balance",
    "largest_domain_fraction",
    "min_segment_size",
    "mean_segment_size",
];

pub(crate) const DOMAIN_CONDITIONED_FEATURE_NAMES: &[&str] = &[
    "smallest_size_fraction",
    "smallest_q1",
    "smallest_q2",
    "smallest_q3",
    "smallest_relative_density",
    "smallest_internal_contact_density",
    "smallest_contact_order",
    "smallest_internal_contact_fraction",
    "largest_size_fraction",
    "largest_q1",
    "largest_q2",
    "largest_q3",
    "largest_relative_density",
    "largest_internal_contact_density",
    "largest_contact_order",
    "largest_internal_contact_fraction",
    "smallest_to_largest_density_ratio",
    "smallest_to_largest_internal_contact_density_ratio",
    "smallest_to_largest_q1_ratio",
    "smallest_to_largest_q2_ratio",
    "smallest_to_largest_q3_ratio",
    "domain_internal_contact_fraction_min",
    "domain_internal_contact_fraction_mean",
    "domain_internal_contact_fraction_max",
    "domain_conductance_min",
    "domain_conductance_mean",
    "domain_conductance_max",
];

pub(crate) const BASE_AND_DOMAIN_FEATURE_NAMES: &[&str] = &[
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "modal_count_distance",
    "domain_q1_mean",
    "domain_q2_mean",
    "domain_q3_mean",
    "domain_q3_min",
    "domain_vol_ratio_mean",
    "domain_density_mean",
    "domain_density_min",
    "contact_q_mean",
    "contact_q_max",
    "n_segments",
    "n_discontinuous",
    "size_balance",
    "largest_domain_fraction",
    "min_segment_size",
    "mean_segment_size",
    "smallest_size_fraction",
    "smallest_q1",
    "smallest_q2",
    "smallest_q3",
    "smallest_relative_density",
    "smallest_internal_contact_density",
    "smallest_contact_order",
    "smallest_internal_contact_fraction",
    "largest_size_fraction",
    "largest_q1",
    "largest_q2",
    "largest_q3",
    "largest_relative_density",
    "largest_internal_contact_density",
    "largest_contact_order",
    "largest_internal_contact_fraction",
    "smallest_to_largest_density_ratio",
    "smallest_to_largest_internal_contact_density_ratio",
    "smallest_to_largest_q1_ratio",
    "smallest_to_largest_q2_ratio",
    "smallest_to_largest_q3_ratio",
    "domain_internal_contact_fraction_min",
    "domain_internal_contact_fraction_mean",
    "domain_internal_contact_fraction_max",
    "domain_conductance_min",
    "domain_conductance_mean",
    "domain_conductance_max",
];

pub(crate) const BOUNDARY_LOCAL_FEATURE_NAMES: &[&str] = &[
    "boundary_inside_helix_min",
    "boundary_inside_helix_mean",
    "boundary_inside_helix_max",
    "boundary_inside_strand_min",
    "boundary_inside_strand_mean",
    "boundary_inside_strand_max",
    "boundary_inside_coil_min",
    "boundary_inside_coil_mean",
    "boundary_inside_coil_max",
    "boundary_sse_terminus_distance_min",
    "boundary_sse_terminus_distance_mean",
    "boundary_sse_terminus_distance_max",
    "boundary_hbond_count_min",
    "boundary_hbond_count_mean",
    "boundary_hbond_count_max",
    "boundary_hbond_energy_kcal_min",
    "boundary_hbond_energy_kcal_mean",
    "boundary_hbond_energy_kcal_max",
    "boundary_bridge_count_min",
    "boundary_bridge_count_mean",
    "boundary_bridge_count_max",
    "boundary_sheet_link_count_min",
    "boundary_sheet_link_count_mean",
    "boundary_sheet_link_count_max",
    "boundary_insulation_w8_min",
    "boundary_insulation_w8_mean",
    "boundary_insulation_w8_max",
    "boundary_insulation_w16_min",
    "boundary_insulation_w16_mean",
    "boundary_insulation_w16_max",
    "boundary_insulation_w32_min",
    "boundary_insulation_w32_mean",
    "boundary_insulation_w32_max",
    "boundary_long_range_contact_density_min",
    "boundary_long_range_contact_density_mean",
    "boundary_long_range_contact_density_max",
    "boundary_bend_change_min",
    "boundary_bend_change_mean",
    "boundary_bend_change_max",
    "boundary_virtual_dihedral_change_min",
    "boundary_virtual_dihedral_change_mean",
    "boundary_virtual_dihedral_change_max",
];

pub(crate) const RELATIVE_HIERARCHY_FEATURE_NAMES: &[&str] = &[
    "sibling_percentile_min_size",
    "sibling_percentile_max_cr",
    "sibling_percentile_density_min",
    "sibling_percentile_mean_density",
    "sibling_percentile_boundary_coil_fraction",
    "sibling_percentile_domain_q1_mean",
    "sibling_percentile_domain_q2_mean",
    "sibling_percentile_domain_q3_mean",
    "sibling_percentile_domain_q3_min",
    "sibling_percentile_domain_vol_ratio_mean",
    "sibling_percentile_domain_density_mean",
    "sibling_percentile_domain_density_min",
    "sibling_percentile_contact_q_mean",
    "sibling_percentile_contact_q_max",
    "sibling_percentile_n_segments",
    "sibling_percentile_n_discontinuous",
    "sibling_percentile_size_balance",
    "sibling_percentile_largest_domain_fraction",
    "sibling_percentile_min_segment_size",
    "sibling_percentile_mean_segment_size",
    "sibling_median_delta_min_size",
    "sibling_median_delta_max_cr",
    "sibling_median_delta_density_min",
    "sibling_median_delta_mean_density",
    "sibling_median_delta_boundary_coil_fraction",
    "sibling_median_delta_domain_q1_mean",
    "sibling_median_delta_domain_q2_mean",
    "sibling_median_delta_domain_q3_mean",
    "sibling_median_delta_domain_q3_min",
    "sibling_median_delta_domain_vol_ratio_mean",
    "sibling_median_delta_domain_density_mean",
    "sibling_median_delta_domain_density_min",
    "sibling_median_delta_contact_q_mean",
    "sibling_median_delta_contact_q_max",
    "sibling_median_delta_n_segments",
    "sibling_median_delta_n_discontinuous",
    "sibling_median_delta_size_balance",
    "sibling_median_delta_largest_domain_fraction",
    "sibling_median_delta_min_segment_size",
    "sibling_median_delta_mean_segment_size",
    "hierarchy_first_appearance_level",
    "hierarchy_persistence_levels",
    "hierarchy_parent_merge_margin",
    "hierarchy_child_merge_margin",
    "sibling_nearest_cr_delta",
    "sibling_nearest_density_delta",
    "hierarchy_path_count",
];

pub(crate) const DISCONTINUITY_FEATURE_NAMES: &[&str] = &[
    "has_discontinuity",
    "segment_same_domain_affinity_min",
    "segment_same_domain_affinity_mean",
    "segment_same_domain_affinity_max",
    "segment_affinity_margin_min",
    "segment_affinity_margin_mean",
    "segment_affinity_margin_max",
    "segment_long_range_internal_capture_min",
    "segment_long_range_internal_capture_mean",
    "segment_long_range_internal_capture_max",
    "segment_interface_span_entropy_min",
    "segment_interface_span_entropy_mean",
    "segment_interface_span_entropy_max",
    "segment_same_domain_sheet_links_min",
    "segment_same_domain_sheet_links_mean",
    "segment_same_domain_sheet_links_max",
];

pub(crate) const BOUNDARY_LOCAL_START: usize = BASE_AND_DOMAIN_FEATURE_NAMES.len();
pub(crate) const BOUNDARY_LOCAL_END: usize =
    BOUNDARY_LOCAL_START + BOUNDARY_LOCAL_FEATURE_NAMES.len();
pub(crate) const DISCONTINUITY_START: usize =
    BOUNDARY_LOCAL_END + RELATIVE_HIERARCHY_FEATURE_NAMES.len();
pub(crate) const DISCONTINUITY_END: usize = DISCONTINUITY_START + DISCONTINUITY_FEATURE_NAMES.len();

pub(crate) const CANDIDATE_FEATURE_NAMES: &[&str] = &[
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "modal_count_distance",
    "domain_q1_mean",
    "domain_q2_mean",
    "domain_q3_mean",
    "domain_q3_min",
    "domain_vol_ratio_mean",
    "domain_density_mean",
    "domain_density_min",
    "contact_q_mean",
    "contact_q_max",
    "n_segments",
    "n_discontinuous",
    "size_balance",
    "largest_domain_fraction",
    "min_segment_size",
    "mean_segment_size",
    "smallest_size_fraction",
    "smallest_q1",
    "smallest_q2",
    "smallest_q3",
    "smallest_relative_density",
    "smallest_internal_contact_density",
    "smallest_contact_order",
    "smallest_internal_contact_fraction",
    "largest_size_fraction",
    "largest_q1",
    "largest_q2",
    "largest_q3",
    "largest_relative_density",
    "largest_internal_contact_density",
    "largest_contact_order",
    "largest_internal_contact_fraction",
    "smallest_to_largest_density_ratio",
    "smallest_to_largest_internal_contact_density_ratio",
    "smallest_to_largest_q1_ratio",
    "smallest_to_largest_q2_ratio",
    "smallest_to_largest_q3_ratio",
    "domain_internal_contact_fraction_min",
    "domain_internal_contact_fraction_mean",
    "domain_internal_contact_fraction_max",
    "domain_conductance_min",
    "domain_conductance_mean",
    "domain_conductance_max",
    "boundary_inside_helix_min",
    "boundary_inside_helix_mean",
    "boundary_inside_helix_max",
    "boundary_inside_strand_min",
    "boundary_inside_strand_mean",
    "boundary_inside_strand_max",
    "boundary_inside_coil_min",
    "boundary_inside_coil_mean",
    "boundary_inside_coil_max",
    "boundary_sse_terminus_distance_min",
    "boundary_sse_terminus_distance_mean",
    "boundary_sse_terminus_distance_max",
    "boundary_hbond_count_min",
    "boundary_hbond_count_mean",
    "boundary_hbond_count_max",
    "boundary_hbond_energy_kcal_min",
    "boundary_hbond_energy_kcal_mean",
    "boundary_hbond_energy_kcal_max",
    "boundary_bridge_count_min",
    "boundary_bridge_count_mean",
    "boundary_bridge_count_max",
    "boundary_sheet_link_count_min",
    "boundary_sheet_link_count_mean",
    "boundary_sheet_link_count_max",
    "boundary_insulation_w8_min",
    "boundary_insulation_w8_mean",
    "boundary_insulation_w8_max",
    "boundary_insulation_w16_min",
    "boundary_insulation_w16_mean",
    "boundary_insulation_w16_max",
    "boundary_insulation_w32_min",
    "boundary_insulation_w32_mean",
    "boundary_insulation_w32_max",
    "boundary_long_range_contact_density_min",
    "boundary_long_range_contact_density_mean",
    "boundary_long_range_contact_density_max",
    "boundary_bend_change_min",
    "boundary_bend_change_mean",
    "boundary_bend_change_max",
    "boundary_virtual_dihedral_change_min",
    "boundary_virtual_dihedral_change_mean",
    "boundary_virtual_dihedral_change_max",
    "sibling_percentile_min_size",
    "sibling_percentile_max_cr",
    "sibling_percentile_density_min",
    "sibling_percentile_mean_density",
    "sibling_percentile_boundary_coil_fraction",
    "sibling_percentile_domain_q1_mean",
    "sibling_percentile_domain_q2_mean",
    "sibling_percentile_domain_q3_mean",
    "sibling_percentile_domain_q3_min",
    "sibling_percentile_domain_vol_ratio_mean",
    "sibling_percentile_domain_density_mean",
    "sibling_percentile_domain_density_min",
    "sibling_percentile_contact_q_mean",
    "sibling_percentile_contact_q_max",
    "sibling_percentile_n_segments",
    "sibling_percentile_n_discontinuous",
    "sibling_percentile_size_balance",
    "sibling_percentile_largest_domain_fraction",
    "sibling_percentile_min_segment_size",
    "sibling_percentile_mean_segment_size",
    "sibling_median_delta_min_size",
    "sibling_median_delta_max_cr",
    "sibling_median_delta_density_min",
    "sibling_median_delta_mean_density",
    "sibling_median_delta_boundary_coil_fraction",
    "sibling_median_delta_domain_q1_mean",
    "sibling_median_delta_domain_q2_mean",
    "sibling_median_delta_domain_q3_mean",
    "sibling_median_delta_domain_q3_min",
    "sibling_median_delta_domain_vol_ratio_mean",
    "sibling_median_delta_domain_density_mean",
    "sibling_median_delta_domain_density_min",
    "sibling_median_delta_contact_q_mean",
    "sibling_median_delta_contact_q_max",
    "sibling_median_delta_n_segments",
    "sibling_median_delta_n_discontinuous",
    "sibling_median_delta_size_balance",
    "sibling_median_delta_largest_domain_fraction",
    "sibling_median_delta_min_segment_size",
    "sibling_median_delta_mean_segment_size",
    "hierarchy_first_appearance_level",
    "hierarchy_persistence_levels",
    "hierarchy_parent_merge_margin",
    "hierarchy_child_merge_margin",
    "sibling_nearest_cr_delta",
    "sibling_nearest_density_delta",
    "hierarchy_path_count",
    "has_discontinuity",
    "segment_same_domain_affinity_min",
    "segment_same_domain_affinity_mean",
    "segment_same_domain_affinity_max",
    "segment_affinity_margin_min",
    "segment_affinity_margin_mean",
    "segment_affinity_margin_max",
    "segment_long_range_internal_capture_min",
    "segment_long_range_internal_capture_mean",
    "segment_long_range_internal_capture_max",
    "segment_interface_span_entropy_min",
    "segment_interface_span_entropy_mean",
    "segment_interface_span_entropy_max",
    "segment_same_domain_sheet_links_min",
    "segment_same_domain_sheet_links_mean",
    "segment_same_domain_sheet_links_max",
];

#[derive(Debug, Clone)]
pub(crate) struct GlobalFeatures {
    pub n_residues: f64,
    pub rg_normalized: f64,
    pub inertia_ratio_21: f64,
    pub inertia_ratio_31: f64,
    pub nonlocal_contact_density: f64,
    pub contact_order: f64,
    pub helix_fraction: f64,
    pub strand_fraction: f64,
    pub coil_fraction: f64,
    pub helix_blocks: f64,
    pub strand_blocks: f64,
    pub peeling_levels: f64,
    pub finest_pus: f64,
    pub candidate_total: f64,
    pub available_count_total: f64,
    /// Bins 0..=19 represent counts 1..=20; bin 20 represents count >=21.
    pub count_histogram: [f64; 21],
    pub modal_count: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct CountFeatures {
    pub num_domains: usize,
    pub values: Vec<f64>,
}

#[derive(Debug, Clone)]
pub(crate) struct CandidateFeatures {
    pub source_index: usize,
    pub canonical: String,
    pub num_domains: usize,
    pub values: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct FeatureMask {
    pub global_count: bool,
    pub domain_conditioned: bool,
    pub boundary_local: bool,
    pub relative_hierarchy: bool,
    pub discontinuity: bool,
}

impl FeatureMask {
    pub(crate) const fn all() -> Self {
        Self {
            global_count: true,
            domain_conditioned: true,
            boundary_local: true,
            relative_hierarchy: true,
            discontinuity: true,
        }
    }
}

impl GlobalFeatures {
    pub(crate) fn to_vec(&self) -> Vec<f64> {
        let mut values = vec![
            self.n_residues,
            self.rg_normalized,
            self.inertia_ratio_21,
            self.inertia_ratio_31,
            self.nonlocal_contact_density,
            self.contact_order,
            self.helix_fraction,
            self.strand_fraction,
            self.coil_fraction,
            self.helix_blocks,
            self.strand_blocks,
            self.peeling_levels,
            self.finest_pus,
            self.candidate_total,
            self.available_count_total,
        ];
        values.extend(self.count_histogram);
        values.push(self.modal_count);
        values
    }
}

impl CandidateFeatures {
    pub(crate) fn base_and_domain_value(&self, name: &str) -> Option<f64> {
        BASE_AND_DOMAIN_FEATURE_NAMES
            .iter()
            .position(|candidate| *candidate == name)
            .and_then(|index| self.values.get(index).copied())
    }

    pub(crate) fn base_and_domain_vec(&self) -> Vec<f64> {
        self.values
            .iter()
            .take(BASE_AND_DOMAIN_FEATURE_NAMES.len())
            .copied()
            .collect()
    }
}

pub(crate) fn count_pair_vector(
    feature_names: &[&str],
    global: &GlobalFeatures,
    left: &CountFeatures,
    right: &CountFeatures,
) -> Result<Vec<f64>, FeatureError> {
    assemble_pair(
        feature_names,
        global,
        COUNT_ITEM_FEATURE_NAMES,
        &left.values,
        &right.values,
    )
}

pub(crate) fn candidate_pair_vector(
    feature_names: &[&str],
    global: &GlobalFeatures,
    left: &CandidateFeatures,
    right: &CandidateFeatures,
) -> Result<Vec<f64>, FeatureError> {
    assemble_pair(
        feature_names,
        global,
        CANDIDATE_FEATURE_NAMES,
        &left.values,
        &right.values,
    )
}

fn assemble_pair(
    feature_names: &[&str],
    global: &GlobalFeatures,
    item_names: &[&str],
    left: &[f64],
    right: &[f64],
) -> Result<Vec<f64>, FeatureError> {
    if left.len() != right.len() || left.len() != item_names.len() {
        return Err(FeatureError::SchemaMismatch);
    }
    let global_values = global.to_vec();
    if global_values.len() != GLOBAL_FEATURE_NAMES.len() {
        return Err(FeatureError::SchemaMismatch);
    }
    let mut values = Vec::with_capacity(feature_names.len());
    for name in feature_names {
        let value = if let Some(raw) = name.strip_prefix("diff__") {
            let index = item_names
                .iter()
                .position(|candidate| *candidate == raw)
                .ok_or(FeatureError::SchemaMismatch)?;
            left[index] - right[index]
        } else if let Some(raw) = name.strip_prefix("abs_diff__") {
            let index = item_names
                .iter()
                .position(|candidate| *candidate == raw)
                .ok_or(FeatureError::SchemaMismatch)?;
            (left[index] - right[index]).abs()
        } else {
            let index = GLOBAL_FEATURE_NAMES
                .iter()
                .position(|candidate| *candidate == *name)
                .ok_or(FeatureError::SchemaMismatch)?;
            global_values[index]
        };
        values.push(value);
    }
    if values.iter().any(|value| !value.is_finite()) {
        return Err(FeatureError::NonFinite("pair_input"));
    }
    Ok(values)
}
