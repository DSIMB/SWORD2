"""Deterministic validated pair construction for the factorized ranker."""

from __future__ import annotations

import hashlib
import math
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations
from numbers import Integral, Real
from typing import Any

import numpy as np
import pandas as pd

from benchmark.factorized_ranker.corpus import CANDIDATE_FIELDS, CHAIN_FIELDS, COUNT_FIELDS
from benchmark.factorized_ranker.schema import (
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    GLOBAL_FEATURES,
    pair_feature_names,
)


_PROHIBITED_COLUMN_TOKENS = ("merizo", "chainsaw")
_U64_MAX = 2**64 - 1


@dataclass(frozen=True)
class PairBatch:
    feature_names: tuple[str, ...]
    x: np.ndarray
    y: np.ndarray
    sample_weight: np.ndarray
    chain_ids: np.ndarray
    left_ids: np.ndarray
    right_ids: np.ndarray


@dataclass(frozen=True)
class _ChainRow:
    chain_id: str
    n_true_domains: int
    values: dict[str, float]


@dataclass(frozen=True)
class _CountRow:
    chain_id: str
    count: int
    item_id: str
    values: dict[str, float]


@dataclass(frozen=True)
class _CandidateRow:
    chain_id: str
    count: int
    item_id: str
    source_index: int
    canonical_delineation: str
    ndo: float
    values: dict[str, float]


@dataclass(frozen=True)
class _PairRecord:
    chain_id: str
    pair_id: str
    left_id: str
    right_id: str
    shared: tuple[float, ...]
    left: tuple[float, ...]
    right: tuple[float, ...]
    label: int
    magnitude: float


def _validate_header(
    frame: pd.DataFrame, expected: Sequence[str], description: str
) -> None:
    if not isinstance(frame, pd.DataFrame):
        raise ValueError(f"{description} must be a pandas DataFrame")
    columns = list(frame.columns)
    prohibited = [
        str(column)
        for column in columns
        if any(token in str(column).casefold() for token in _PROHIBITED_COLUMN_TOKENS)
    ]
    if prohibited:
        raise ValueError(f"{description} contains prohibited predictor columns")
    if len(set(columns)) != len(columns):
        raise ValueError(f"{description} contains duplicate columns")
    if tuple(columns) != tuple(expected):
        raise ValueError(f"{description} header does not match the exact Task 8 schema")


def _rows(frame: pd.DataFrame) -> list[dict[str, Any]]:
    columns = list(frame.columns)
    return [
        dict(zip(columns, values, strict=True))
        for values in frame.itertuples(index=False, name=None)
    ]


def _text(value: object, description: str) -> str:
    if not isinstance(value, str) or not value or "\0" in value:
        raise ValueError(f"{description} must be nonempty text without NUL")
    return value


def _finite(value: object, description: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{description} must be a finite real number")
    if isinstance(value, str):
        if not value:
            raise ValueError(f"{description} is empty")
        try:
            number = float(value)
        except ValueError as error:
            raise ValueError(f"{description} is not numeric") from error
        if not math.isfinite(number) or format(number, ".17g") != value:
            raise ValueError(f"{description} is not canonical finite numeric text")
        return number
    if not isinstance(value, Real):
        raise ValueError(f"{description} must be a finite real number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{description} must be finite")
    return number


def _integer(value: object, description: str, *, minimum: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{description} must be an integer")
    if isinstance(value, str):
        try:
            number = int(value)
        except ValueError as error:
            raise ValueError(f"{description} must be an ordinary integer") from error
        if str(number) != value:
            raise ValueError(f"{description} must be an ordinary integer")
    elif isinstance(value, Integral):
        number = int(value)
    elif isinstance(value, Real):
        numeric = float(value)
        if not math.isfinite(numeric) or not numeric.is_integer():
            raise ValueError(f"{description} must be an exact integer")
        number = int(numeric)
    else:
        raise ValueError(f"{description} must be an integer")
    if number < minimum:
        raise ValueError(f"{description} is out of range")
    if number > _U64_MAX:
        raise ValueError(f"{description} exceeds u64")
    return number


def _feature_tuple(
    requested: Sequence[str], frozen: Sequence[str], description: str
) -> tuple[str, ...]:
    try:
        names = tuple(requested)
    except TypeError as error:
        raise ValueError(f"{description} must be a feature-name sequence") from error
    if any(not isinstance(name, str) for name in names):
        raise ValueError(f"{description} contains a non-string name")
    if len(set(names)) != len(names):
        raise ValueError(f"{description} contains duplicate names")
    positions = {name: index for index, name in enumerate(frozen)}
    if any(name not in positions for name in names):
        raise ValueError(f"{description} contains a non-feature column")
    indices = [positions[name] for name in names]
    if indices != sorted(indices):
        raise ValueError(f"{description} is not an ordered frozen-schema subsequence")
    return names


def _stable_hash(parts: Sequence[str]) -> str:
    digest = hashlib.sha256()
    for part in parts:
        if not isinstance(part, str):
            raise ValueError("stable ID parts must be strings")
        encoded = part.encode("utf-8")
        if len(encoded) > _U64_MAX:
            raise ValueError("stable ID part is too large for u64 framing")
        digest.update(len(encoded).to_bytes(8, "big", signed=False))
        digest.update(encoded)
    return digest.hexdigest()


def _candidate_id(chain_id: str, canonical_delineation: str) -> str:
    return hashlib.sha256((chain_id + "\0" + canonical_delineation).encode("utf-8")).hexdigest()


def _hash_text(value: object, description: str) -> str:
    text = _text(value, description)
    if len(text) != 64 or any(character not in "0123456789abcdef" for character in text):
        raise ValueError(f"{description} must be a lowercase SHA-256")
    return text


def _validate_chains(chains: pd.DataFrame) -> dict[str, _ChainRow]:
    _validate_header(chains, CHAIN_FIELDS, "chains")
    validated: dict[str, _ChainRow] = {}
    for raw in _rows(chains):
        chain_id = _text(raw["chain_id"], "chain_id")
        if chain_id in validated:
            raise ValueError("duplicate chains.chain_id")
        true_count = _integer(raw["n_true_domains"], "n_true_domains", minimum=1)
        values = {
            name: _finite(raw[name], f"chains.{name}") for name in GLOBAL_FEATURES
        }
        validated[chain_id] = _ChainRow(chain_id, true_count, values)
    if not validated:
        raise ValueError("chains table is empty")
    return validated


def _validate_counts(
    counts: pd.DataFrame, chains: dict[str, _ChainRow]
) -> dict[str, list[_CountRow]]:
    _validate_header(counts, COUNT_FIELDS, "counts")
    grouped: dict[str, list[_CountRow]] = {chain_id: [] for chain_id in chains}
    seen: set[tuple[str, int]] = set()
    item_ids: set[str] = set()
    for raw in _rows(counts):
        chain_id = _text(raw["chain_id"], "counts.chain_id")
        if chain_id not in chains:
            raise ValueError("counts row does not join to chains")
        count = _integer(raw["count_num_domains"], "count_num_domains", minimum=1)
        key = (chain_id, count)
        if key in seen:
            raise ValueError("duplicate count row")
        seen.add(key)
        values = {
            name: _finite(raw[name], f"counts.{name}") for name in COUNT_ITEM_FEATURES
        }
        if values["count_num_domains"] != float(count):
            raise ValueError("count_num_domains is inconsistent")
        item_id = _stable_hash(("count", chain_id, str(count)))
        if item_id in item_ids:
            raise ValueError("count item ID collision")
        item_ids.add(item_id)
        grouped[chain_id].append(_CountRow(chain_id, count, item_id, values))
    if any(not rows for rows in grouped.values()):
        raise ValueError("every accepted chain must have at least one count row")
    for rows in grouped.values():
        rows.sort(key=lambda row: row.item_id)
    return grouped


def _validate_candidates(
    candidates: pd.DataFrame, chains: dict[str, _ChainRow]
) -> dict[str, list[_CandidateRow]]:
    _validate_header(candidates, CANDIDATE_FIELDS, "candidates")
    grouped: dict[str, list[_CandidateRow]] = {chain_id: [] for chain_id in chains}
    candidate_ids: set[str] = set()
    canonicals: set[tuple[str, str]] = set()
    source_indices: set[tuple[str, int]] = set()
    for raw in _rows(candidates):
        item_id = _hash_text(raw["candidate_id"], "candidate_id")
        chain_id = _text(raw["chain_id"], "candidates.chain_id")
        canonical = _text(raw["canonical_delineation"], "canonical_delineation")
        if chain_id not in chains:
            raise ValueError("candidate row does not join to chains")
        source_index = _integer(raw["source_index"], "source_index", minimum=0)
        if item_id in candidate_ids:
            raise ValueError("duplicate candidate_id")
        if (chain_id, canonical) in canonicals:
            raise ValueError("duplicate candidate canonical identity")
        if (chain_id, source_index) in source_indices:
            raise ValueError("duplicate candidate source identity")
        candidate_ids.add(item_id)
        canonicals.add((chain_id, canonical))
        source_indices.add((chain_id, source_index))
        if item_id != _candidate_id(chain_id, canonical):
            raise ValueError("candidate_id does not match Task 8 identity")

        _finite(raw["legacy_distance"], "legacy_distance")
        values = {
            name: _finite(raw[name], f"candidates.{name}")
            for name in CANDIDATE_FEATURES
        }
        count = _integer(raw["num_domains"], "num_domains", minimum=1)
        if values["num_domains"] != float(count):
            raise ValueError("candidate num_domains is inconsistent")
        true_count = _integer(raw["n_true_domains"], "n_true_domains", minimum=1)
        predicted_count = _integer(
            raw["n_pred_domains"], "n_pred_domains", minimum=1
        )
        if true_count != chains[chain_id].n_true_domains:
            raise ValueError("candidate n_true_domains disagrees with chains")
        if predicted_count != count:
            raise ValueError("candidate num_domains disagrees with n_pred_domains")
        ndo = _finite(raw["ndo"], "ndo")
        for name in (
            "iou",
            "boundary_f1_10",
            "matched_dice",
            "d_count_acc",
            "S",
        ):
            _finite(raw[name], name)
        oracle = _integer(raw["is_oracle_s"], "is_oracle_s", minimum=0)
        if oracle not in {0, 1}:
            raise ValueError("is_oracle_s must be 0 or 1")
        grouped[chain_id].append(
            _CandidateRow(
                chain_id,
                count,
                item_id,
                source_index,
                canonical,
                ndo,
                values,
            )
        )
    if any(not rows for rows in grouped.values()):
        raise ValueError("every accepted chain must have at least one candidate row")
    for rows in grouped.values():
        rows.sort(key=lambda row: (row.count, row.item_id))
    return grouped


def _vector(values: Sequence[float], description: str) -> np.ndarray:
    try:
        materialized = list(values)
    except TypeError as error:
        raise ValueError(f"{description} must be a one-dimensional sequence") from error
    converted = np.asarray(
        [_finite(value, f"{description} value") for value in materialized],
        dtype=np.float64,
    )
    if converted.ndim != 1 or not np.isfinite(converted).all():
        raise ValueError(f"{description} must be a finite one-dimensional vector")
    return converted


def pair_vector(
    shared: Sequence[float],
    left: Sequence[float],
    right: Sequence[float],
) -> np.ndarray:
    shared_vector = _vector(shared, "shared")
    left_vector = _vector(left, "left")
    right_vector = _vector(right, "right")
    if left_vector.shape != right_vector.shape:
        raise ValueError("left/right item vectors have different lengths")
    difference = left_vector - right_vector
    result = np.concatenate((shared_vector, difference, np.abs(difference))).astype(
        np.float64, copy=False
    )
    if result.ndim != 1 or not np.isfinite(result).all():
        raise ValueError("pair vector is non-finite")
    return result


def _sampling_settings(seed: object, cap: object) -> tuple[int, int]:
    if not isinstance(seed, Integral) or isinstance(seed, (bool, np.bool_)):
        raise ValueError("seed must be a non-boolean integer")
    seed_value = int(seed)
    if seed_value < 0:
        raise ValueError("seed must be nonnegative")
    if seed_value > _U64_MAX:
        raise ValueError("seed exceeds u64")
    if not isinstance(cap, Integral) or isinstance(cap, (bool, np.bool_)):
        raise ValueError("max_unordered_pairs must be a non-boolean integer")
    cap_value = int(cap)
    if cap_value <= 0:
        raise ValueError("max_unordered_pairs must be positive")
    return seed_value, cap_value


def _sample_key(seed: int, chain_id: str, pair_id: str) -> str:
    return _stable_hash(("sample", str(seed), chain_id, pair_id))


def _select_records(
    records: Sequence[_PairRecord], seed: int, cap: int
) -> list[_PairRecord]:
    ordered = sorted(records, key=lambda record: record.pair_id)
    if len({record.pair_id for record in ordered}) != len(ordered):
        raise ValueError("eligible pair IDs are not unique")
    if len(ordered) <= cap:
        return ordered

    magnitudes = pd.Series(
        [record.magnitude for record in ordered], dtype=np.float64
    )
    try:
        raw_codes = pd.qcut(
            magnitudes,
            q=4,
            labels=False,
            duplicates="drop",
        ).to_numpy()
    except (TypeError, ValueError) as error:
        raise ValueError("pair magnitudes could not be quantile-binned") from error
    if len(raw_codes) != len(ordered):
        raise ValueError("qcut returned the wrong number of codes")
    if pd.isna(raw_codes).all():
        codes = np.zeros(len(ordered), dtype=np.int64)
    elif pd.isna(raw_codes).any():
        raise ValueError("qcut returned partially missing codes")
    else:
        observed_codes = sorted({int(code) for code in raw_codes})
        remap = {code: index for index, code in enumerate(observed_codes)}
        codes = np.asarray([remap[int(code)] for code in raw_codes], dtype=np.int64)

    bins: dict[int, list[_PairRecord]] = {}
    for record, code in zip(ordered, codes, strict=True):
        bins.setdefault(int(code), []).append(record)
    if not bins:
        bins = {0: list(ordered)}
    for code in bins:
        bins[code].sort(
            key=lambda record: (
                _sample_key(seed, record.chain_id, record.pair_id),
                record.pair_id,
            )
        )

    n_nonempty = len(bins)
    base_quota, remainder = divmod(cap, n_nonempty)
    selected_by_bin: dict[int, list[_PairRecord]] = {
        code: records_in_bin[: min(base_quota, len(records_in_bin))]
        for code, records_in_bin in bins.items()
    }
    nonexhausted = [
        code
        for code, records_in_bin in bins.items()
        if len(selected_by_bin[code]) < len(records_in_bin)
    ]
    nonexhausted.sort(
        key=lambda code: (
            -len(bins[code]),
            tuple(sorted(record.pair_id for record in bins[code])),
            code,
        )
    )
    for code in nonexhausted[:remainder]:
        selected_by_bin[code].append(bins[code][len(selected_by_bin[code])])

    selected = [record for code in sorted(selected_by_bin) for record in selected_by_bin[code]]
    selected_ids = {record.pair_id for record in selected}
    if len(selected) < cap:
        remaining = [record for record in ordered if record.pair_id not in selected_ids]
        remaining.sort(
            key=lambda record: (
                _sample_key(seed, record.chain_id, record.pair_id),
                record.pair_id,
            )
        )
        selected.extend(remaining[: cap - len(selected)])
    if len(selected) != cap or len({record.pair_id for record in selected}) != cap:
        raise ValueError("pair cap selection did not produce the exact unique capacity")
    return sorted(selected, key=lambda record: record.pair_id)


def _record(
    chain_id: str,
    left_id: str,
    right_id: str,
    shared: tuple[float, ...],
    left: tuple[float, ...],
    right: tuple[float, ...],
    label: int,
    magnitude: float,
) -> _PairRecord:
    if not left_id < right_id:
        raise ValueError("eligible pair orientation is not canonical")
    pair_id = _stable_hash(("pair", left_id, right_id))
    if label not in {0, 1} or not math.isfinite(magnitude) or magnitude <= 0:
        raise ValueError("eligible pair label/magnitude is invalid")
    return _PairRecord(
        chain_id,
        pair_id,
        left_id,
        right_id,
        shared,
        left,
        right,
        label,
        magnitude,
    )


def _pair_batch(
    records_by_chain: dict[str, list[_PairRecord]],
    feature_names: tuple[str, ...],
    seed: int,
    cap: int,
) -> PairBatch:
    selected: list[_PairRecord] = []
    selected_counts: dict[str, int] = {}
    for chain_id in sorted(records_by_chain):
        chain_records = _select_records(records_by_chain[chain_id], seed, cap)
        if chain_records:
            selected_counts[chain_id] = len(chain_records)
            selected.extend(chain_records)
    selected.sort(key=lambda record: (record.chain_id, record.pair_id))

    vectors: list[np.ndarray] = []
    labels: list[int] = []
    weights: list[float] = []
    chain_ids: list[str] = []
    left_ids: list[str] = []
    right_ids: list[str] = []
    for record in selected:
        weight = np.float64(1.0 / (2 * selected_counts[record.chain_id]))
        for orientation in (0, 1):
            if orientation == 0:
                left_id, right_id = record.left_id, record.right_id
                left, right = record.left, record.right
                label = record.label
            else:
                left_id, right_id = record.right_id, record.left_id
                left, right = record.right, record.left
                label = 1 - record.label
            vectors.append(pair_vector(record.shared, left, right))
            labels.append(label)
            weights.append(float(weight))
            chain_ids.append(record.chain_id)
            left_ids.append(left_id)
            right_ids.append(right_id)

    n_features = len(feature_names)
    if vectors:
        x = np.vstack(vectors).astype(np.float64, copy=False)
    else:
        x = np.empty((0, n_features), dtype=np.float64)
    y = np.asarray(labels, dtype=np.int8)
    sample_weight = np.asarray(weights, dtype=np.float64)
    chain_width = max((len(chain_id) for chain_id in chain_ids), default=1)
    chain_array = np.asarray(chain_ids, dtype=f"<U{max(1, chain_width)}")
    left_array = np.asarray(left_ids, dtype="<U64")
    right_array = np.asarray(right_ids, dtype="<U64")
    expected_rows = len(labels)
    if x.shape != (expected_rows, n_features):
        raise ValueError("pair feature matrix has the wrong shape")
    if any(array.shape != (expected_rows,) for array in (y, sample_weight, chain_array, left_array, right_array)):
        raise ValueError("pair metadata arrays have the wrong shape")
    if x.dtype != np.float64 or sample_weight.dtype != np.float64 or y.dtype != np.int8:
        raise ValueError("pair arrays have the wrong dtype")
    if not np.isfinite(x).all() or not np.isfinite(sample_weight).all():
        raise ValueError("pair arrays contain non-finite values")
    if not set(y.tolist()).issubset({0, 1}):
        raise ValueError("pair labels are not binary")
    shared_count = next(
        (index for index, name in enumerate(feature_names) if name.startswith("diff__")),
        n_features,
    )
    item_count = (n_features - shared_count) // 2
    if shared_count + 2 * item_count != n_features:
        raise ValueError("pair feature-name blocks have inconsistent lengths")
    for row in range(0, expected_rows, 2):
        if (
            chain_array[row] != chain_array[row + 1]
            or left_array[row] != right_array[row + 1]
            or right_array[row] != left_array[row + 1]
            or int(y[row]) + int(y[row + 1]) != 1
            or sample_weight[row] != sample_weight[row + 1]
            or not np.array_equal(x[row, :shared_count], x[row + 1, :shared_count])
            or not np.array_equal(
                x[row, shared_count : shared_count + item_count],
                -x[row + 1, shared_count : shared_count + item_count],
            )
            or not np.array_equal(
                x[row, shared_count + item_count :],
                x[row + 1, shared_count + item_count :],
            )
        ):
            raise ValueError("mirrored pair rows are inconsistent")
    for chain_id, unordered_count in selected_counts.items():
        mask = chain_array == chain_id
        if int(mask.sum()) != 2 * unordered_count:
            raise ValueError("pair chain row count is inconsistent")
        expected_weight = np.float64(1.0 / (2 * unordered_count))
        if not np.all(sample_weight[mask] == expected_weight):
            raise ValueError("pair chain weights are inconsistent")
    for array in (x, y, sample_weight, chain_array, left_array, right_array):
        array.flags.writeable = False
    return PairBatch(
        feature_names,
        x,
        y,
        sample_weight,
        chain_array,
        left_array,
        right_array,
    )


def build_count_pairs(
    chains: pd.DataFrame,
    counts: pd.DataFrame,
    *,
    shared_features: Sequence[str] = GLOBAL_FEATURES,
    item_features: Sequence[str] = COUNT_ITEM_FEATURES,
    seed: int = 37,
    max_unordered_pairs: int = 64,
) -> PairBatch:
    shared_names = _feature_tuple(shared_features, GLOBAL_FEATURES, "shared_features")
    item_names = _feature_tuple(item_features, COUNT_ITEM_FEATURES, "item_features")
    if not shared_names and not item_names:
        raise ValueError("pair feature vector cannot have zero columns")
    feature_names = tuple(pair_feature_names(shared_names, item_names))
    seed_value, cap = _sampling_settings(seed, max_unordered_pairs)
    chain_rows = _validate_chains(chains)
    count_rows = _validate_counts(counts, chain_rows)

    records_by_chain: dict[str, list[_PairRecord]] = {}
    for chain_id in sorted(chain_rows):
        chain = chain_rows[chain_id]
        items = count_rows[chain_id]
        if chain.n_true_domains not in {item.count for item in items}:
            records_by_chain[chain_id] = []
            continue
        shared = tuple(chain.values[name] for name in shared_names)
        records: list[_PairRecord] = []
        for first, second in combinations(items, 2):
            left, right = (first, second) if first.item_id < second.item_id else (second, first)
            left_error = abs(left.count - chain.n_true_domains)
            right_error = abs(right.count - chain.n_true_domains)
            if left_error == right_error:
                continue
            records.append(
                _record(
                    chain_id,
                    left.item_id,
                    right.item_id,
                    shared,
                    tuple(left.values[name] for name in item_names),
                    tuple(right.values[name] for name in item_names),
                    int(left_error < right_error),
                    float(abs(left_error - right_error)),
                )
            )
        records_by_chain[chain_id] = records
    return _pair_batch(records_by_chain, feature_names, seed_value, cap)


def build_candidate_pairs(
    chains: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    shared_features: Sequence[str] = GLOBAL_FEATURES,
    item_features: Sequence[str] = CANDIDATE_FEATURES,
    seed: int = 37,
    max_unordered_pairs: int = 64,
) -> PairBatch:
    shared_names = _feature_tuple(shared_features, GLOBAL_FEATURES, "shared_features")
    item_names = _feature_tuple(item_features, CANDIDATE_FEATURES, "item_features")
    if not shared_names and not item_names:
        raise ValueError("pair feature vector cannot have zero columns")
    feature_names = tuple(pair_feature_names(shared_names, item_names))
    seed_value, cap = _sampling_settings(seed, max_unordered_pairs)
    chain_rows = _validate_chains(chains)
    candidate_rows = _validate_candidates(candidates, chain_rows)

    records_by_chain: dict[str, list[_PairRecord]] = {}
    for chain_id in sorted(chain_rows):
        chain = chain_rows[chain_id]
        shared = tuple(chain.values[name] for name in shared_names)
        groups: dict[int, list[_CandidateRow]] = {}
        for item in candidate_rows[chain_id]:
            groups.setdefault(item.count, []).append(item)
        records: list[_PairRecord] = []
        for count in sorted(groups):
            items = sorted(groups[count], key=lambda item: item.item_id)
            for left, right in combinations(items, 2):
                if left.ndo == right.ndo:
                    continue
                records.append(
                    _record(
                        chain_id,
                        left.item_id,
                        right.item_id,
                        shared,
                        tuple(left.values[name] for name in item_names),
                        tuple(right.values[name] for name in item_names),
                        int(left.ndo > right.ndo),
                        abs(left.ndo - right.ndo),
                    )
                )
        records_by_chain[chain_id] = records
    return _pair_batch(records_by_chain, feature_names, seed_value, cap)
