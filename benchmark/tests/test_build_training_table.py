from benchmark.build_training_table import FIELDNAMES, _score_candidates


def test_fieldnames_include_new_reranker_features():
    for col in ("boundary_coil_fraction", "energy_z", "modal_count_distance"):
        assert col in FIELDNAMES, f"missing column: {col}"


def test_score_candidates_passes_through_new_features(tmp_path):
    # score_choppings needs real chopping strings; use a trivial 1-domain case
    # so the ground-truth scoring path is exercised without needing a real PDB.
    reference = {"testchain": ("A", "1-20", 20)}
    candidates = [
        {
            "output_dir": str(tmp_path),
            "num_domains": "1",
            "min_size": "20",
            "max_cr": "0.1",
            "density_min": "1.0",
            "mean_density": "2.0",
            "delineation": "0-19",
            "boundary_coil_fraction": "0.75",
            "energy_z": "-2.5",
            "modal_count_distance": "0.0",
        }
    ]
    rows = _score_candidates("testchain", candidates, reference)
    assert len(rows) == 1
    assert rows[0]["boundary_coil_fraction"] == "0.75"
    assert rows[0]["energy_z"] == "-2.5"
    assert rows[0]["modal_count_distance"] == "0.0"
