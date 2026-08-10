from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
import pytest

import benchmark.factorized_ranker.pairs as pairs
from benchmark.factorized_ranker.corpus import CANDIDATE_FIELDS, CHAIN_FIELDS, COUNT_FIELDS
from benchmark.factorized_ranker.pairs import (
    PairBatch,
    build_candidate_pairs,
    build_count_pairs,
    pair_vector,
)
from benchmark.factorized_ranker.schema import (
    BASE_CANDIDATE_FEATURES,
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    GLOBAL_FEATURES,
    pair_feature_names,
)


def _candidate_id(chain_id: str, delineation: str) -> str:
    return hashlib.sha256((chain_id + "\0" + delineation).encode()).hexdigest()


def _chain_row(chain_id: str, n_true: int, offset: float = 0.0) -> dict[str, object]:
    return {
        "chain_id": chain_id,
        "n_true_domains": n_true,
        **{
            name: offset + index / 100.0
            for index, name in enumerate(GLOBAL_FEATURES, start=1)
        },
    }


def _count_row(chain_id: str, count: int, offset: float = 0.0) -> dict[str, object]:
    row: dict[str, object] = {
        "chain_id": chain_id,
        **{
            name: offset + count * 10.0 + index / 100.0
            for index, name in enumerate(COUNT_ITEM_FEATURES, start=1)
        },
    }
    row["count_num_domains"] = count
    return row


def _candidate_row(
    chain_id: str,
    n_true: int,
    count: int,
    source_index: int,
    ndo: float,
    *,
    delineation: str | None = None,
    offset: float = 0.0,
) -> dict[str, object]:
    canonical = delineation or f"0-{source_index} {source_index + 1}-{source_index + 2}"
    row: dict[str, object] = {
        "candidate_id": _candidate_id(chain_id, canonical),
        "chain_id": chain_id,
        "canonical_delineation": canonical,
        "source_index": source_index,
        "legacy_distance": offset + source_index / 10.0,
        **{
            name: offset + count * 10.0 + source_index + index / 1000.0
            for index, name in enumerate(CANDIDATE_FEATURES, start=1)
        },
        "n_true_domains": n_true,
        "n_pred_domains": count,
        "ndo": ndo,
        "iou": 0.5,
        "boundary_f1_10": 0.5,
        "matched_dice": 0.5,
        "d_count_acc": float(count == n_true),
        "S": ndo,
        "is_oracle_s": 0,
    }
    row["num_domains"] = count
    return row


def _frame_hash(parts_: tuple[str, ...]) -> str:
    framed = b"".join(len(part.encode()).to_bytes(8, "big") + part.encode() for part in parts_)
    return hashlib.sha256(framed).hexdigest()


def _two_candidate_fixture():
    chain_id = "chainA"
    rows = [
        _candidate_row(chain_id, 2, 2, 0, 0.0, delineation="0-4 5-9"),
        _candidate_row(chain_id, 2, 2, 1, 0.0, delineation="0-5 6-9"),
    ]
    rows.sort(key=lambda row: row["candidate_id"])
    rows[0]["ndo"] = 0.9
    rows[1]["ndo"] = 0.1
    return (
        pd.DataFrame([_chain_row(chain_id, 2)], columns=CHAIN_FIELDS),
        pd.DataFrame(rows, columns=CANDIDATE_FIELDS),
    )


def _batch_bytes(batch: PairBatch) -> tuple[object, ...]:
    return (
        batch.feature_names,
        batch.x.shape,
        batch.x.dtype.str,
        batch.x.tobytes(),
        batch.y.dtype.str,
        batch.y.tobytes(),
        batch.sample_weight.dtype.str,
        batch.sample_weight.tobytes(),
        batch.chain_ids.dtype.str,
        batch.chain_ids.tobytes(),
        batch.left_ids.dtype.str,
        batch.left_ids.tobytes(),
        batch.right_ids.dtype.str,
        batch.right_ids.tobytes(),
    )


def test_mirrored_pair_has_reversed_diff_and_same_absolute_part() -> None:
    batch = build_candidate_pairs(
        *_two_candidate_fixture(),
        item_features=BASE_CANDIDATE_FEATURES,
    )
    assert batch.y.tolist() == [1, 0]
    shared = len(GLOBAL_FEATURES)
    item = len(BASE_CANDIDATE_FEATURES)
    np.testing.assert_allclose(batch.x[0, :shared], batch.x[1, :shared])
    np.testing.assert_allclose(
        batch.x[0, shared : shared + item],
        -batch.x[1, shared : shared + item],
    )
    np.testing.assert_allclose(
        batch.x[0, shared + item :],
        batch.x[1, shared + item :],
    )
    assert batch.feature_names == pair_feature_names(
        GLOBAL_FEATURES, BASE_CANDIDATE_FEATURES
    )


def test_pair_cap_and_chain_weights_are_exact() -> None:
    chain_id = "large"
    candidates = [
        _candidate_row(chain_id, 3, 3, index, index / 100.0)
        for index in range(12)
    ]
    batch = build_candidate_pairs(
        pd.DataFrame([_chain_row(chain_id, 3)], columns=CHAIN_FIELDS),
        pd.DataFrame(candidates, columns=CANDIDATE_FIELDS),
        max_unordered_pairs=64,
        seed=37,
    )
    assert len(batch.y) == 128
    for observed_chain in np.unique(batch.chain_ids):
        mask = batch.chain_ids == observed_chain
        assert batch.sample_weight[mask].sum() == pytest.approx(1.0)
        assert np.all(batch.sample_weight[mask] == np.float64(1.0 / 128.0))


def test_count_labels_skip_equal_error_and_require_true_count_row() -> None:
    chains = pd.DataFrame(
        [_chain_row("usable", 3), _chain_row("missing", 3)], columns=CHAIN_FIELDS
    )
    counts = pd.DataFrame(
        [
            _count_row("usable", 1),
            _count_row("usable", 2),
            _count_row("usable", 3),
            _count_row("usable", 4),
            _count_row("missing", 1),
            _count_row("missing", 2),
        ],
        columns=COUNT_FIELDS,
    )
    batch = build_count_pairs(chains, counts)
    assert set(batch.chain_ids.tolist()) == {"usable"}
    # Six unordered pairs minus the equal-error 2-vs-4 pair.
    assert len(batch.y) == 10
    assert np.array_equal(batch.y[0::2], 1 - batch.y[1::2])


def test_candidate_pairs_never_cross_chain_or_count_and_pool_one_cap() -> None:
    chains = pd.DataFrame(
        [_chain_row("a", 2), _chain_row("b", 4)], columns=CHAIN_FIELDS
    )
    candidates = []
    for chain_id, true_count in (("a", 2), ("b", 4)):
        source = 0
        for count in (2, 3):
            for index in range(9):
                candidates.append(
                    _candidate_row(
                        chain_id,
                        true_count,
                        count,
                        source,
                        source / 100.0,
                    )
                )
                source += 1
    batch = build_candidate_pairs(
        chains,
        pd.DataFrame(candidates, columns=CANDIDATE_FIELDS),
        max_unordered_pairs=64,
    )
    assert len(batch.y) == 2 * 64 * 2
    for index in range(0, len(batch.y), 2):
        assert batch.chain_ids[index] == batch.chain_ids[index + 1]
        assert batch.left_ids[index] == batch.right_ids[index + 1]
        assert batch.right_ids[index] == batch.left_ids[index + 1]


@pytest.mark.parametrize("table", ["chains", "counts", "candidates"])
@pytest.mark.parametrize("defect", ["missing", "extra", "reordered", "duplicate"])
def test_exact_table_headers_are_required(table: str, defect: str) -> None:
    chains, candidates = _two_candidate_fixture()
    counts = pd.DataFrame([_count_row("chainA", 2)], columns=COUNT_FIELDS)
    frames = {"chains": chains, "counts": counts, "candidates": candidates}
    frame = frames[table].copy()
    columns = list(frame.columns)
    if defect == "missing":
        frame = frame.drop(columns=columns[-1])
    elif defect == "extra":
        frame["merizo_prediction"] = 1.0
    elif defect == "reordered":
        frame = frame[[columns[1], columns[0], *columns[2:]]]
    else:
        frame.columns = [*columns[:-1], columns[-2]]
    with pytest.raises(ValueError):
        if table == "candidates":
            build_candidate_pairs(frames["chains"], frame)
        elif table == "counts":
            build_count_pairs(frames["chains"], frame)
        else:
            build_candidate_pairs(frame, frames["candidates"])


@pytest.mark.parametrize(
    "defect",
    [
        "duplicate_chain",
        "unknown_chain",
        "duplicate_candidate_id",
        "duplicate_canonical",
        "duplicate_source",
        "bad_candidate_id",
        "count_disagreement",
    ],
)
def test_candidate_join_identity_defects_fail_closed(defect: str) -> None:
    chains, candidates = _two_candidate_fixture()
    if defect == "duplicate_chain":
        chains = pd.concat([chains, chains.iloc[[0]]], ignore_index=True)
    elif defect == "unknown_chain":
        candidates.loc[0, "chain_id"] = "other"
    elif defect == "duplicate_candidate_id":
        candidates.loc[1, "candidate_id"] = candidates.loc[0, "candidate_id"]
    elif defect == "duplicate_canonical":
        candidates.loc[1, "canonical_delineation"] = candidates.loc[0, "canonical_delineation"]
    elif defect == "duplicate_source":
        candidates.loc[1, "source_index"] = candidates.loc[0, "source_index"]
    elif defect == "bad_candidate_id":
        candidates.loc[0, "candidate_id"] = "0" * 64
    else:
        candidates.loc[0, "n_pred_domains"] = 3
    with pytest.raises(ValueError):
        build_candidate_pairs(chains, candidates)


def test_count_duplicates_and_unknown_chains_fail_closed() -> None:
    chains = pd.DataFrame([_chain_row("a", 2)], columns=CHAIN_FIELDS)
    counts = pd.DataFrame([_count_row("a", 2), _count_row("a", 3)], columns=COUNT_FIELDS)
    with pytest.raises(ValueError, match="duplicate"):
        build_count_pairs(chains, pd.concat([counts, counts.iloc[[0]]], ignore_index=True))
    counts.loc[0, "chain_id"] = "missing"
    with pytest.raises(ValueError):
        build_count_pairs(chains, counts)


@pytest.mark.parametrize(
    "shared,item",
    [
        ((GLOBAL_FEATURES[1], GLOBAL_FEATURES[0]), BASE_CANDIDATE_FEATURES),
        ((GLOBAL_FEATURES[0], GLOBAL_FEATURES[0]), BASE_CANDIDATE_FEATURES),
        (("n_true_domains",), BASE_CANDIDATE_FEATURES),
        (GLOBAL_FEATURES, ("ndo",)),
        ((), ()),
    ],
)
def test_feature_names_must_be_unique_ordered_schema_subsequences(shared, item) -> None:
    with pytest.raises(ValueError):
        build_candidate_pairs(
            *_two_candidate_fixture(), shared_features=shared, item_features=item
        )


def test_input_permutation_produces_byte_identical_batches() -> None:
    chain_rows = [_chain_row("b", 3, 1.0), _chain_row("a", 2, 2.0)]
    count_rows = [
        _count_row("b", 2),
        _count_row("b", 3),
        _count_row("b", 4),
        _count_row("a", 1),
        _count_row("a", 2),
        _count_row("a", 3),
    ]
    forward = build_count_pairs(
        pd.DataFrame(chain_rows, columns=CHAIN_FIELDS),
        pd.DataFrame(count_rows, columns=COUNT_FIELDS),
    )
    reverse = build_count_pairs(
        pd.DataFrame(list(reversed(chain_rows)), columns=CHAIN_FIELDS),
        pd.DataFrame(list(reversed(count_rows)), columns=COUNT_FIELDS),
    )
    assert _batch_bytes(forward) == _batch_bytes(reverse)


def test_empty_batch_has_exact_shapes_and_nonobject_dtypes() -> None:
    chains, candidates = _two_candidate_fixture()
    candidates["ndo"] = 0.5
    batch = build_candidate_pairs(chains, candidates)
    assert batch.x.shape == (0, len(batch.feature_names))
    assert batch.y.shape == batch.sample_weight.shape == batch.chain_ids.shape == (0,)
    assert batch.x.dtype == batch.sample_weight.dtype == np.float64
    assert batch.y.dtype == np.int8
    assert batch.chain_ids.dtype == np.dtype("<U1")
    assert batch.left_ids.dtype == batch.right_ids.dtype == np.dtype("<U64")
    assert not batch.x.flags.writeable
    assert not batch.y.flags.writeable


@pytest.mark.parametrize("field,value", [("ndo", np.nan), ("num_domains", 0), ("source_index", 1.5)])
def test_malformed_or_nonfinite_candidate_values_fail(field: str, value: object) -> None:
    chains, candidates = _two_candidate_fixture()
    candidates[field] = candidates[field].astype(object)
    candidates.loc[0, field] = value
    with pytest.raises(ValueError):
        build_candidate_pairs(chains, candidates)


@pytest.mark.parametrize(
    "seed,cap",
    [
        (True, 64),
        (37.0, 64),
        (-1, 64),
        (2**64, 64),
        (37, True),
        (37, 64.0),
        (37, 0),
    ],
)
def test_invalid_sampling_settings_fail(seed: object, cap: object) -> None:
    with pytest.raises(ValueError):
        build_candidate_pairs(
            *_two_candidate_fixture(), seed=seed, max_unordered_pairs=cap
        )


def test_pair_vector_is_exact_finite_float64() -> None:
    vector = pair_vector([10.0], [3.0, 8.0], [1.0, 10.0])
    assert vector.dtype == np.float64
    np.testing.assert_array_equal(vector, np.asarray([10.0, 2.0, -2.0, 2.0, 2.0]))
    with pytest.raises(ValueError):
        pair_vector([], [1.0], [1.0, 2.0])
    with pytest.raises(ValueError):
        pair_vector([np.inf], [1.0], [2.0])


def test_framed_count_and_pair_ids_are_unambiguous() -> None:
    assert _frame_hash(("count", "a:2", "3")) != _frame_hash(("count", "a", "2:3"))
    assert _frame_hash(("pair", "ab", "c")) != _frame_hash(("pair", "a", "bc"))


def test_different_seed_changes_only_over_cap_selection() -> None:
    chain_id = "large"
    chains = pd.DataFrame([_chain_row(chain_id, 3)], columns=CHAIN_FIELDS)
    candidates = pd.DataFrame(
        [_candidate_row(chain_id, 3, 3, index, index / 100.0) for index in range(13)],
        columns=CANDIDATE_FIELDS,
    )
    first = build_candidate_pairs(chains, candidates, seed=37, max_unordered_pairs=16)
    second = build_candidate_pairs(chains, candidates, seed=38, max_unordered_pairs=16)
    assert first.feature_names == second.feature_names
    assert set(zip(first.left_ids[0::2], first.right_ids[0::2])) != set(
        zip(second.left_ids[0::2], second.right_ids[0::2])
    )
    uncapped_a = build_candidate_pairs(chains, candidates, seed=37, max_unordered_pairs=100)
    uncapped_b = build_candidate_pairs(chains, candidates, seed=38, max_unordered_pairs=100)
    assert _batch_bytes(uncapped_a) == _batch_bytes(uncapped_b)


def _sampling_records(magnitudes: list[float]) -> list[pairs._PairRecord]:
    return [
        pairs._PairRecord(
            chain_id="chain",
            pair_id=f"{index:064x}",
            left_id=f"{2 * index:064x}",
            right_id=f"{2 * index + 1:064x}",
            shared=(),
            left=(),
            right=(),
            label=1,
            magnitude=magnitude,
        )
        for index, magnitude in enumerate(magnitudes)
    ]


def test_over_cap_sampling_represents_each_nonempty_quartile() -> None:
    selected = pairs._select_records(
        _sampling_records([float(value) for value in range(1, 101)]),
        seed=37,
        cap=4,
    )
    selected_magnitudes = sorted(record.magnitude for record in selected)
    assert len(selected_magnitudes) == 4
    assert sum(value <= 25 for value in selected_magnitudes) == 1
    assert sum(25 < value <= 50 for value in selected_magnitudes) == 1
    assert sum(50 < value <= 75 for value in selected_magnitudes) == 1
    assert sum(value > 75 for value in selected_magnitudes) == 1


def test_one_unique_sampling_magnitude_is_one_deterministic_bin() -> None:
    records = _sampling_records([1.0] * 20)
    selected = pairs._select_records(records, seed=37, cap=5)
    expected = sorted(
        records,
        key=lambda record: (
            pairs._stable_hash(("sample", "37", "chain", record.pair_id)),
            record.pair_id,
        ),
    )[:5]
    assert {record.pair_id for record in selected} == {
        record.pair_id for record in expected
    }
