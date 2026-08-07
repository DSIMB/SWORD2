"""Ordered feature contract shared by factorized-ranker tooling."""

SCHEMA_VERSION = 1
SEED = 37
MAX_COUNT_HISTOGRAM_BIN = 20

GLOBAL_FEATURES = (
    "chain_n_residues", "chain_rg_normalized",
    "chain_inertia_ratio_21", "chain_inertia_ratio_31",
    "chain_nonlocal_contact_density", "chain_contact_order",
    "chain_helix_fraction", "chain_strand_fraction", "chain_coil_fraction",
    "chain_helix_blocks", "chain_strand_blocks",
    "chain_peeling_levels", "chain_finest_pus",
    "chain_candidate_total", "chain_available_count_total",
    *(f"chain_count_hist_{value}" for value in range(1, 21)),
    "chain_count_hist_21_plus", "chain_modal_count",
)

COUNT_SUMMARY_SOURCES = (
    "legacy_distance", "min_size", "max_cr", "density_min", "mean_density",
    "contact_q_mean", "contact_q_max", "n_segments", "n_discontinuous",
    "boundary_coil_fraction",
)
COUNT_SUMMARIES = tuple(
    f"count_{source}_{stat}"
    for source in COUNT_SUMMARY_SOURCES for stat in ("min", "mean", "max")
)
COUNT_ITEM_FEATURES = (
    "count_num_domains", "count_n_candidates", "count_candidate_fraction",
    "count_modal_distance", "count_has_lower", "count_has_higher",
    "count_lower_gap", "count_higher_gap", *COUNT_SUMMARIES,
    *(f"{name}_delta_lower" for name in COUNT_SUMMARIES),
    *(f"{name}_delta_higher" for name in COUNT_SUMMARIES),
)

BASE_CANDIDATE_FEATURES = (
    "num_domains", "min_size", "max_cr", "density_min", "mean_density",
    "boundary_coil_fraction", "modal_count_distance",
    "domain_q1_mean", "domain_q2_mean", "domain_q3_mean", "domain_q3_min",
    "domain_vol_ratio_mean", "domain_density_mean", "domain_density_min",
    "contact_q_mean", "contact_q_max", "n_segments", "n_discontinuous",
    "size_balance", "largest_domain_fraction", "min_segment_size",
    "mean_segment_size",
)

DOMAIN_SIDE_MEASURES = (
    "size_fraction", "q1", "q2", "q3", "relative_density",
    "internal_contact_density", "contact_order", "internal_contact_fraction",
)
DOMAIN_CONDITIONED_FEATURES = (
    *(f"smallest_{name}" for name in DOMAIN_SIDE_MEASURES),
    *(f"largest_{name}" for name in DOMAIN_SIDE_MEASURES),
    "smallest_to_largest_density_ratio",
    "smallest_to_largest_internal_contact_density_ratio",
    "smallest_to_largest_q1_ratio", "smallest_to_largest_q2_ratio",
    "smallest_to_largest_q3_ratio",
    *(f"domain_internal_contact_fraction_{stat}" for stat in ("min", "mean", "max")),
    *(f"domain_conductance_{stat}" for stat in ("min", "mean", "max")),
)

BOUNDARY_MEASURES = (
    "inside_helix", "inside_strand", "inside_coil", "sse_terminus_distance",
    "hbond_count", "hbond_energy_kcal", "bridge_count", "sheet_link_count",
    "insulation_w8", "insulation_w16", "insulation_w32",
    "long_range_contact_density", "bend_change", "virtual_dihedral_change",
)
BOUNDARY_LOCAL_FEATURES = tuple(
    f"boundary_{name}_{stat}"
    for name in BOUNDARY_MEASURES for stat in ("min", "mean", "max")
)

RELATIVE_CORE_FEATURES = tuple(
    name for name in BASE_CANDIDATE_FEATURES
    if name not in {"num_domains", "modal_count_distance"}
)
RELATIVE_HIERARCHY_FEATURES = (
    *(f"sibling_percentile_{name}" for name in RELATIVE_CORE_FEATURES),
    *(f"sibling_median_delta_{name}" for name in RELATIVE_CORE_FEATURES),
    "hierarchy_first_appearance_level", "hierarchy_persistence_levels",
    "hierarchy_parent_merge_margin", "hierarchy_child_merge_margin",
    "sibling_nearest_cr_delta", "sibling_nearest_density_delta",
    "hierarchy_path_count",
)

DISCONTINUITY_MEASURES = (
    "same_domain_affinity", "affinity_margin", "long_range_internal_capture",
    "interface_span_entropy", "same_domain_sheet_links",
)
DISCONTINUITY_FEATURES = (
    "has_discontinuity",
    *(f"segment_{name}_{stat}" for name in DISCONTINUITY_MEASURES
      for stat in ("min", "mean", "max")),
)
CANDIDATE_FEATURES = (
    *BASE_CANDIDATE_FEATURES, *DOMAIN_CONDITIONED_FEATURES,
    *BOUNDARY_LOCAL_FEATURES, *RELATIVE_HIERARCHY_FEATURES,
    *DISCONTINUITY_FEATURES,
)


def pair_feature_names(shared, item):
    return (*shared, *(f"diff__{name}" for name in item),
            *(f"abs_diff__{name}" for name in item))
