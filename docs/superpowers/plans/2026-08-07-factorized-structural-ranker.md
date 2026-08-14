# Factorized Structural Ranker Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build, train, freeze, embed, and validate a two-head structural candidate ranker that selects SWORD2's domain count and then its best partition at that count, without any competitor-derived feature or runtime ML dependency.

**Architecture:** Rust retains a typed, bounded candidate lattice plus DSSP, Peeling, contact, and merge-provenance evidence and emits the exact features used at inference. Python builds fail-closed CATH labels, leakage-resistant five-fold splits, deterministic pairwise training data, grouped ablations, and two compact `GradientBoostingClassifier` artifacts; a deterministic exporter turns the frozen trees into static Rust arrays. The opt-in Rust selector applies symmetrized pair probabilities and normalized Borda twice—first over available counts, then over candidates at the winning count—and falls back as a whole chain to the unchanged legacy selector on any validation or feature error.

**Tech Stack:** Rust 2021 (`sword2-lib`, `sword2-cli`, clap, serde, rayon); Python 3.10+ (`numpy`, `pandas`, `scipy`, `scikit-learn`, `pytest`); existing SWORD2 benchmark runners and CATH-17287/CATH-663 metadata.

## Global Constraints

- Do not use Merizo, Chainsaw, or any other predictor as a training target, feature, pseudo-label, or runtime input. They appear only in the final locked comparison.
- Preserve the current dirty worktree. In particular, finish rather than discard the canonical-numbering edits already present in `benchmark/build_training_table.py` and its tests; stage only files named by each task.
- The factorized selector stays off by default behind `--use-factorized-ranker` until every acceptance gate passes. Flag-off output must remain byte-identical to the current legacy path.
- Define the runtime lattice as the current first `ParseMeasure` shortlist with `option_alt_dist=false`, `n_dom=0`, `alt_b=3`, `alt_l=3`, before any legacy count filtering. Preserve legacy order, deduplicate canonical partitions, and retain at most the first three candidates per available count. Training dumps and inference must use this identical lattice.
- Candidate partitions use zero-based inclusive residue indices and must cover every clean-chain residue exactly once. Empty, malformed, reversed, overlapping, duplicated, out-of-range, or incomplete partitions are invalid.
- On a fresh run, factorized features consume the in-memory `DsspChain`, `ContactMatrix`, Peeling iterations, and `ComputeMeasure` provenance. On a cache hit where complete typed evidence is unavailable, emit one structured warning and fall back; do not reconstruct partial evidence from permissive text parsing.
- Use the existing logistic contact model with `d0=6.0`, `delta=1.5`. “Nonlocal” and “long-range” both mean sequence separation `>= 8` residues.
- DSSP crossing-bond features mean the retained best donor/acceptor bonds present in `DsspChain`, deduplicated by `(donor, acceptor)`; energy is stored in kcal/mol by dividing the retained integer value by `1000.0`. Bridge features use retained partner IDs and sheet labels because ladder objects are not retained.
- Pair vectors are `[shared chain context, left-right item features, abs(left-right) item features]`. Skip target ties, mirror every retained pair, cap each head at 64 unordered pairs per chain with deterministic seed 37, and give each eligible chain total ordered-pair weight exactly `1.0`.
- The offline grid is exactly trees `{64, 96}` × learning rate `{0.03, 0.05}` × minimum leaf samples `{32, 64}`, always depth `3`, `loss="log_loss"`, and `random_state=37`.
- Each frozen head has at most 96 depth-3 trees; together they have at most 192 trees and 2,880 nodes. Rust inference is `f64`, has no model parser or new ML dependency, and must match Python within `1e-10`.
- Development uses five deterministic connected-component folds that keep a four-character PDB ID or an exact sorted CATH family-label multiset in one fold. CATH metadata is never a runtime feature.
- CATH-663 is not accepted by the development training CLI. After the development schema, folds, retained families, and model hashes are committed, run the locked CATH-663 benchmark once and make no model or feature choice from that result.
- Acceptance requires: NDO `> 0.8389`; a strictly positive paired 95% CI versus Merizo on 663 common chains; a strictly positive paired CI versus Chainsaw on its common successful subset; count accuracy `>= 0.745`; BF1@10 `>= 0.620`; contiguous and discontinuous NDO regressions each no worse than `0.005`; median runtime overhead `<= 15%`; maximum per-process peak RSS overhead `<= 10%`; and all parity/tests passing.
- After every Rust checkpoint run `cargo check` and `cargo test`. After every benchmark-tooling checkpoint run focused tests and `benchmark/.venv/bin/python -m pytest benchmark/tests -q`.

---

## File Structure

### Python tooling

- `benchmark/factorized_ranker/__init__.py` — package marker and public schema version.
- `benchmark/factorized_ranker/integrity.py` — canonical CATH mapping, strict partition validation, and rejection records.
- `benchmark/factorized_ranker/schema.py` — the sole ordered Python feature contract and pair-vector names.
- `benchmark/factorized_ranker/structure_features.py` — independent NumPy reference formulas used by synthetic and parity tests.
- `benchmark/factorized_ranker/lattice_features.py` — sibling-relative features, count summaries, and neighbor deltas.
- `benchmark/factorized_ranker/corpus.py` — normalized corpus writers and deterministic manifests.
- `benchmark/factorized_ranker/folds.py` — exact-family/PDB connected components and balanced fold assignment.
- `benchmark/factorized_ranker/pairs.py` — deterministic weighted count/candidate pair construction.
- `benchmark/factorized_ranker/ranking.py` — symmetrized probabilities, normalized Borda, and tie rules.
- `benchmark/factorized_ranker/training.py` — fixed hyperparameter grid, grouped CV, and ordered feature-family ablations.
- `benchmark/factorized_ranker/model_artifact.py` — sklearn tree freezing, artifact validation, reference prediction, hashes, and Rust rendering.
- `benchmark/build_factorized_corpus.py`, `benchmark/build_factorized_folds.py`, `benchmark/train_factorized_ranker.py`, `benchmark/export_factorized_ranker.py` — thin reproducible CLIs.
- `benchmark/evaluate_factorized_acceptance.py` — locked accuracy/resource gates and JSON/Markdown report.
- `benchmark/models/` — committed fold manifest, OOF/ablation report, frozen models, provenance hashes, and golden vectors; bulk CSV corpora remain ignored.

### Rust runtime

- `sword2-lib/src/sword/factorized_ranker/mod.rs` — public-in-crate orchestration, full-chain fallback result, and selection seam.
- `sword2-lib/src/sword/factorized_ranker/partition.rs` — strict delineation parser and canonical partition key.
- `sword2-lib/src/sword/factorized_ranker/lattice.rs` — typed bounded shortlist and merge/Peeling provenance.
- `sword2-lib/src/sword/factorized_ranker/schema.rs` — named global/count/candidate feature records and ordered pair-vector assembly.
- `sword2-lib/src/sword/factorized_ranker/features.rs` — global, count, geometry, contact, and domain-conditioned extraction.
- `sword2-lib/src/sword/factorized_ranker/boundary.rs` — boundary-local DSSP/contact/hinge extraction.
- `sword2-lib/src/sword/factorized_ranker/discontinuity.rs` — discontinuous segment-association extraction.
- `sword2-lib/src/sword/factorized_ranker/model.rs` — static tree validation/traversal and Borda inference.
- `sword2-lib/src/sword/factorized_ranker/generated_model.rs` — generated static arrays only; never hand-edit.
- `sword2-lib/tests/factorized_ranker_golden.rs` — committed Python/Rust parity and integration fixture.

---

### Task 1: Fail-closed training labels and deterministic rejection records

**Files:**
- Create: `benchmark/factorized_ranker/__init__.py`
- Create: `benchmark/factorized_ranker/integrity.py`
- Modify: `benchmark/build_training_table.py`
- Modify: `benchmark/candidate_geometry.py`
- Test: `benchmark/tests/test_factorized_integrity.py`
- Test: `benchmark/tests/test_build_training_table.py`
- Test: `benchmark/tests/test_candidate_geometry.py`

**Interfaces:**
- Consumes: `CathEntry`, `numbering_from_pdb`, `map_author_chopping`, and raw SWORD delineations.
- Produces: `CanonicalReference`, `ValidatedPartition`, `RejectionRecord`, and `ScoredChain`; Tasks 7 and 8 use these exact types.

- [ ] **Step 1: Write failing integrity tests**

```python
# benchmark/tests/test_factorized_integrity.py
import pytest

from benchmark.factorized_ranker.integrity import (
    RejectionCode,
    validate_partition,
)


def test_partition_requires_exact_coverage():
    with pytest.raises(ValueError, match=RejectionCode.CANDIDATE_INCOMPLETE_COVERAGE.value):
        validate_partition("0-2 4-5", n_residues=6, declared_domains=2)


@pytest.mark.parametrize(
    ("delineation", "code"),
    [
        ("0-3 3-5", RejectionCode.CANDIDATE_OVERLAP),
        ("0-2 3-6", RejectionCode.CANDIDATE_OUT_OF_RANGE),
        ("2-0 1-5", RejectionCode.CANDIDATE_PARSE_FAILED),
        ("0-2 3-5", RejectionCode.CANDIDATE_DOMAIN_COUNT_MISMATCH),
    ],
)
def test_partition_rejects_invalid_candidates(delineation, code):
    declared = 3 if code is RejectionCode.CANDIDATE_DOMAIN_COUNT_MISMATCH else 2
    with pytest.raises(ValueError, match=code.value):
        validate_partition(delineation, n_residues=6, declared_domains=declared)


def test_discontinuous_partition_has_stable_canonical_key():
    parsed = validate_partition("4-5 0-1;3-3 2-2", 6, 3)
    assert parsed.canonical_delineation == "0-1;3 2 4-5"
    assert parsed.residue_to_domain == (0, 0, 1, 0, 2, 2)
```

Append tests in `test_build_training_table.py` that pass a missing chain-cache PDB and a PDB whose mapped residue/domain counts disagree with CATH. Assert `_score_candidates(...).rows == []` and rejection codes `missing_chain_pdb`, `residue_count_mismatch`, and `true_domain_count_mismatch`. Add a test in `test_candidate_geometry.py` asserting invalid bounds raise `ValueError` instead of returning an all-zero geometry record.

- [ ] **Step 2: Run the tests and confirm the old permissive behavior fails**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_integrity.py \
  benchmark/tests/test_build_training_table.py \
  benchmark/tests/test_candidate_geometry.py -q
```

Expected: failures because the package/types do not exist, raw author numbering still falls back, and invalid geometry is zero-filled.

- [ ] **Step 3: Add the strict integrity types and parser**

```python
# benchmark/factorized_ranker/integrity.py
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Literal

from benchmark.datasets import CathEntry, strip_cath_labels
from benchmark.numbering import map_author_chopping, numbering_from_pdb


class RejectionCode(str, Enum):
    MISSING_REFERENCE = "missing_reference"
    MISSING_CHAIN_PDB = "missing_chain_pdb"
    AUTHOR_MAPPING_FAILED = "author_mapping_failed"
    RESIDUE_COUNT_MISMATCH = "residue_count_mismatch"
    TRUE_DOMAIN_COUNT_MISMATCH = "true_domain_count_mismatch"
    CANDIDATE_PARSE_FAILED = "candidate_parse_failed"
    CANDIDATE_OUT_OF_RANGE = "candidate_out_of_range"
    CANDIDATE_OVERLAP = "candidate_overlap"
    CANDIDATE_INCOMPLETE_COVERAGE = "candidate_incomplete_coverage"
    CANDIDATE_DOMAIN_COUNT_MISMATCH = "candidate_domain_count_mismatch"
    SCORING_FAILED = "scoring_failed"
    NONFINITE_FEATURE = "nonfinite_feature"
    SCHEMA_MISMATCH = "schema_mismatch"


@dataclass(frozen=True)
class RejectionRecord:
    chain_id: str
    scope: Literal["chain", "candidate"]
    code: RejectionCode
    detail: str
    delineation: str | None = None


@dataclass(frozen=True)
class CanonicalReference:
    chopping: str
    n_residues: int
    n_domains: int


@dataclass(frozen=True)
class ValidatedPartition:
    canonical_delineation: str
    domains: tuple[tuple[tuple[int, int], ...], ...]
    residue_to_domain: tuple[int, ...]


def map_cath_reference(entry: CathEntry, pdb_path: Path) -> CanonicalReference:
    numbering = numbering_from_pdb(pdb_path, chain_id=entry.chain_id)
    if numbering.n_residues != entry.n_residues:
        raise ValueError(RejectionCode.RESIDUE_COUNT_MISMATCH.value)
    chopping = map_author_chopping(
        strip_cath_labels(entry.chopping), numbering, chain_id=entry.chain_id
    )
    domains = tuple(part for part in chopping.split(",") if part)
    if len(domains) != entry.n_domains:
        raise ValueError(RejectionCode.TRUE_DOMAIN_COUNT_MISMATCH.value)
    return CanonicalReference(chopping, numbering.n_residues, len(domains))


def validate_partition(
    delineation: str, n_residues: int, declared_domains: int
) -> ValidatedPartition:
    parsed: list[list[tuple[int, int]]] = []
    for raw_domain in delineation.strip().split():
        segments: list[tuple[int, int]] = []
        for raw_segment in raw_domain.split(";"):
            try:
                pieces = raw_segment.split("-", 1)
                start = int(pieces[0])
                end = int(pieces[1]) if len(pieces) == 2 else start
            except (ValueError, IndexError) as exc:
                raise ValueError(RejectionCode.CANDIDATE_PARSE_FAILED.value) from exc
            if start < 0 or end < start:
                raise ValueError(RejectionCode.CANDIDATE_PARSE_FAILED.value)
            if end >= n_residues:
                raise ValueError(RejectionCode.CANDIDATE_OUT_OF_RANGE.value)
            segments.append((start, end))
        if not segments:
            raise ValueError(RejectionCode.CANDIDATE_PARSE_FAILED.value)
        parsed.append(sorted(segments))
    if len(parsed) != declared_domains:
        raise ValueError(RejectionCode.CANDIDATE_DOMAIN_COUNT_MISMATCH.value)

    parsed.sort(key=lambda domain: min(start for start, _ in domain))
    owner = [-1] * n_residues
    for domain_index, segments in enumerate(parsed):
        for start, end in segments:
            for residue in range(start, end + 1):
                if owner[residue] != -1:
                    raise ValueError(RejectionCode.CANDIDATE_OVERLAP.value)
                owner[residue] = domain_index
    if any(value == -1 for value in owner):
        raise ValueError(RejectionCode.CANDIDATE_INCOMPLETE_COVERAGE.value)

    def render(segment: tuple[int, int]) -> str:
        return str(segment[0]) if segment[0] == segment[1] else f"{segment[0]}-{segment[1]}"

    frozen = tuple(tuple(domain) for domain in parsed)
    canonical = " ".join(";".join(map(render, domain)) for domain in frozen)
    return ValidatedPartition(canonical, frozen, tuple(owner))
```

- [ ] **Step 4: Make training-table construction fail closed**

Change `_load_reference()` to return `dict[str, CathEntry]`; delete `_true_chopping_for_chain()` and all raw-chopping fallback logic. Add:

```python
@dataclass(frozen=True)
class ScoredChain:
    rows: list[dict[str, Any]]
    rejections: list[RejectionRecord]
```

Make `_score_candidates(...) -> ScoredChain`. A missing reference/PDB or canonical mapping exception returns no rows and one chain rejection. Validate every raw candidate before `score_choppings`; candidate failures add candidate-level records and do not become zero features. If a candidate feature cannot parse as a finite float, reject it as `NONFINITE_FEATURE`. Sort rejection output by `(chain_id, scope, code.value, delineation or "", detail)`. Add `--rejections`, default `<output stem>.rejections.csv`, and make both CSV and dump-directory builders write the same five columns: `chain_id,scope,code,detail,delineation`.

In `candidate_geometry_features`, replace its invalid-input `default` returns with `raise ValueError("invalid candidate partition")`; a valid continuous candidate may still legitimately return zeros for conditional quantities.

- [ ] **Step 5: Run focused and full Python tests**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_integrity.py \
  benchmark/tests/test_build_training_table.py \
  benchmark/tests/test_candidate_geometry.py -q
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: all pass; existing callers now unwrap `.rows`, and no test observes raw-author fallback.

- [ ] **Step 6: Commit only the integrity checkpoint**

```bash
git add benchmark/factorized_ranker/__init__.py \
  benchmark/factorized_ranker/integrity.py \
  benchmark/build_training_table.py benchmark/candidate_geometry.py \
  benchmark/tests/test_factorized_integrity.py \
  benchmark/tests/test_build_training_table.py \
  benchmark/tests/test_candidate_geometry.py
git commit -m "fix: fail closed when building structural training labels"
```

---

### Task 2: Typed Rust partitions, bounded lattice, and exact legacy-selection seam

**Files:**
- Create: `sword2-lib/src/sword/factorized_ranker/mod.rs`
- Create: `sword2-lib/src/sword/factorized_ranker/partition.rs`
- Create: `sword2-lib/src/sword/factorized_ranker/lattice.rs`
- Modify: `sword2-lib/src/sword/mod.rs`
- Modify: `sword2-lib/src/sword/parse_measure.rs`
- Modify: `sword2-lib/src/sword/compute_measure.rs`
- Test: inline tests in the three modified/new modules

**Interfaces:**
- Consumes: `MeasureLine` and the existing first-pass `ParseMeasure` behavior.
- Produces: `ParsedPartition`, `CandidateRecord`, `CandidateLattice`, `HierarchyEvidence`, and `LegacySelection`; later Rust feature and inference tasks depend on these names.

- [ ] **Step 1: Write strict parser and lattice regression tests**

```rust
// partition.rs tests
#[test]
fn rejects_overlap_gap_and_out_of_range() {
    assert!(matches!(parse_partition("0-3 3-5", 6), Err(FeatureError::Overlap(3))));
    assert!(matches!(parse_partition("0-2 4-5", 6), Err(FeatureError::IncompleteCoverage)));
    assert!(matches!(parse_partition("0-2 3-6", 6), Err(FeatureError::OutOfRange { .. })));
}

#[test]
fn canonicalizes_domain_and_segment_order() {
    let p = parse_partition("4-5 3;0-1 2", 6).unwrap();
    assert_eq!(p.canonical, "0-1;3 2 4-5");
    assert_eq!(p.residue_to_domain, vec![0, 0, 1, 1, 2, 2]);
}

// lattice.rs tests
#[test]
fn lattice_is_before_legacy_count_filter_and_caps_three_per_count() {
    let measures = synthetic_measures_with_four_two_domain_candidates();
    let lattice = CandidateLattice::from_first_pass(&measures, &identity_indices(), 12).unwrap();
    assert_eq!(lattice.groups[&2].len(), 3);
    assert!(lattice.groups.contains_key(&3));
}

#[test]
fn legacy_selection_is_stable() {
    let old = legacy_selection_reference(&synthetic_measure_strings());
    let new = select_legacy(&synthetic_measures(), 120, false, None);
    assert_eq!(new.num_domains, old.0);
    assert_eq!(new.measure_line, old.1);
}
```

- [ ] **Step 2: Run the tests and verify the new modules are absent**

Run: `cargo test -p sword2-lib factorized_ranker --no-fail-fast`

Expected: compile failure because `factorized_ranker` and its types are not defined.

- [ ] **Step 3: Implement the strict Rust partition contract**

Use these exact public-in-crate types in `partition.rs`:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Segment { pub start: usize, pub end: usize }

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Domain {
    pub segments: Vec<Segment>,
    pub residues: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ParsedPartition {
    pub domains: Vec<Domain>,
    pub residue_to_domain: Vec<usize>,
    pub canonical: String,
}

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub(crate) enum FeatureError {
    #[error("empty or malformed delineation")]
    Malformed,
    #[error("segment {start}-{end} exceeds chain length {chain_len}")]
    OutOfRange { start: usize, end: usize, chain_len: usize },
    #[error("residue {0} occurs in more than one domain")]
    Overlap(usize),
    #[error("partition does not cover every residue")]
    IncompleteCoverage,
    #[error("declared domain count does not match parsed partition")]
    DomainCountMismatch,
    #[error("required structural evidence is unavailable: {0}")]
    MissingContext(&'static str),
    #[error("non-finite feature {0}")]
    NonFinite(&'static str),
    #[error("feature schema or vector length mismatch")]
    SchemaMismatch,
}

pub(crate) fn parse_partition(text: &str, chain_len: usize) -> Result<ParsedPartition, FeatureError>;
```

Parsing must accept `start-end` and singleton `start`, sort segments and then domains by their lowest residue, reject duplicates/overlaps, require exact coverage, and render singletons without `-`. This implementation must mirror Task 1's Python canonicalization byte-for-byte.

- [ ] **Step 4: Preserve merge genealogy instead of inferring it from strings**

In `compute_measure.rs`, add:

```rust
#[derive(Debug, Clone, Default)]
pub struct MeasureProvenance {
    pub canonical_pu_key: String,
    pub incoming_merge_qualities: Vec<f64>,
    pub outgoing_merge_qualities: Vec<f64>,
    pub distinct_parent_keys: Vec<String>,
    pub hierarchy_path_count: u64,
}

#[derive(Debug, Clone)]
pub struct MeasureCorpus {
    pub lines: Vec<MeasureLine>,
    /// Same length/order as `lines`; provenance[i] describes lines[i].
    pub provenance: Vec<MeasureProvenance>,
}

pub fn compute_measure_from_data_with_provenance(
    pu_contacts: &[(usize, usize, f64)],
    pu_delineation: &[(usize, usize, usize)],
) -> MeasureCorpus;
```

Extend `MergeResult` with `parent_key: String`. Its child key is the sorted PU-domain representation already held in `new_domains`; canonicalize each domain's PU IDs numerically and sort domains by their minimum PU ID. During deduplication, retain all distinct parent keys and all finite `ratio_pdp` values even when only one `MeasureLine` is emitted. Record the same ratio as an outgoing quality on the parent. Set the initial partition's path count to one; for each emitted child, sum the already-computed path counts of its distinct parents with `u64::saturating_add`. Keep existing `compute_measure* -> Vec<MeasureLine>` wrappers returning `.lines`, so flag-off behavior and callers remain unchanged.

Derive later margins exactly as `largest - second_largest`, or `0.0` with fewer than two finite values. Use the recursively accumulated `hierarchy_path_count`, not merely the number of immediate parents.

- [ ] **Step 5: Add typed shortlist indices without changing legacy formatting**

Add `parse_measure_indices(...) -> Vec<usize>` beside `parse_measure()`. It must run the same filtering/Jones logic but carry original `MeasureLine` indices; implement the existing string API by mapping the returned indices back to `to_line()`. Define:

```rust
pub(crate) struct HierarchyEvidence {
    pub first_appearance_level: usize,
    pub persistence_levels: usize,
    pub parent_merge_margin: f64,
    pub child_merge_margin: f64,
    pub hierarchy_path_count: u64,
}

pub(crate) struct CandidateRecord {
    pub source_index: usize,
    pub measure: MeasureLine,
    pub partition: ParsedPartition,
    pub legacy_distance: f64,
    /// Filled from aligned `MeasureProvenance` plus Peeling iterations when
    /// complete typed evidence is available.
    pub hierarchy: Option<HierarchyEvidence>,
}

pub(crate) struct CandidateLattice {
    pub candidates: Vec<CandidateRecord>,
    pub groups: std::collections::BTreeMap<usize, Vec<usize>>,
}

impl CandidateLattice {
    pub(crate) fn from_first_pass(
        measures: &[MeasureLine],
        shortlisted_indices: &[usize],
        chain_len: usize,
    ) -> Result<Self, FeatureError>;
}
```

`from_first_pass` preserves the shortlist order, deduplicates by `(num_domains, canonical)`, and ignores entries after the first three valid entries for a count. It does not filter around the legacy selected count.

Extract the current cross-level distance/count-calibration loop into:

```rust
pub(crate) struct LegacySelection {
    pub num_domains: usize,
    pub measure_line: String,
    pub source_index: Option<usize>,
}

pub(crate) fn select_legacy(
    relevant: &[MeasureLine],
    chain_len: usize,
    use_count_calibration: bool,
    count_lambda: Option<f64>,
) -> LegacySelection;
```

Call this function from `run_pipeline` even when the factorized selector is eventually enabled; it is the fallback and supplies the count tie-break.

- [ ] **Step 6: Verify legacy behavior and the workspace**

Run:

```bash
cargo fmt -- --check
cargo check
cargo test -p sword2-lib factorized_ranker --no-fail-fast
cargo test
```

Expected: all pass; the existing output-format tests remain unchanged.

- [ ] **Step 7: Commit the typed lattice checkpoint**

```bash
git add sword2-lib/src/sword/factorized_ranker \
  sword2-lib/src/sword/mod.rs sword2-lib/src/sword/parse_measure.rs \
  sword2-lib/src/sword/compute_measure.rs
git commit -m "refactor: retain typed candidate lattice provenance"
```

---

### Task 3: Retain complete DSSP, contact, and Peeling evidence for opt-in runs

**Files:**
- Modify: `sword2-lib/src/sword/mod.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/mod.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/lattice.rs`
- Test: inline tests in `factorized_ranker/mod.rs`

**Interfaces:**
- Consumes: `DsspResult.chain`, `PeelingOutput.contact_matrix`, Peeling iterations, and `MeasureCorpus`.
- Produces: validated `StructuralContext<'a>` and completed hierarchy evidence for Task 5 onward.

- [ ] **Step 1: Write context-validation and cache-fallback tests**

```rust
#[test]
fn context_rejects_dimension_mismatch() {
    let context = synthetic_context_with_lengths(8, 7, 8);
    assert!(matches!(context.validate(), Err(FeatureError::MissingContext("DSSP/chain length mismatch"))));
}

#[test]
fn unavailable_typed_cache_requests_whole_chain_fallback() {
    let result = prepare_factorized_context(None, None, &[], None);
    assert!(matches!(result, Err(FeatureError::MissingContext(_))));
}
```

- [ ] **Step 2: Confirm the tests fail before context retention**

Run: `cargo test -p sword2-lib factorized_ranker::tests::context -- --nocapture`

Expected: compile failure for missing `StructuralContext`.

- [ ] **Step 3: Add the structural context and explicit residue mapping**

```rust
pub(crate) struct StructuralContext<'a> {
    pub ca_coords: &'a [[f64; 3]],
    pub dssp: &'a crate::dssp::DsspChain,
    pub contacts: &'a crate::peeling::contact_matrix::ContactMatrix,
    pub iterations: &'a [crate::peeling::algorithm::IterationResult],
    pub measure_provenance: &'a [MeasureProvenance],
    pub dssp_index_for_residue: Vec<usize>,
}

impl StructuralContext<'_> {
    pub(crate) fn validate(&self) -> Result<(), FeatureError>;
}
```

Build `dssp_index_for_residue` from 1-based DSSP residues whose `aa != '!'`. Require its length, CA length, and contact-matrix length to match. Require nonempty Peeling iterations and finite coordinates/contact values used by features. Do not silently align by truncation.

- [ ] **Step 4: Retain fresh-run evidence in `run_pipeline`**

Capture `Some(DsspResult)` when `run_dssp` executes instead of discarding it. Capture `MeasureCorpus` from `compute_measure_from_data_with_provenance` on a fresh Peeling run. Keep the current file-based cache path for legacy output, but set the typed factorized context to `None` on cache hits because it lacks DSSP bonds, the contact matrix, and genealogy. No feature work runs unless `config.use_factorized_ranker` is true in Task 14; therefore flag-off CPU, I/O, and output remain unchanged.

For each candidate's partition boundary set, define a boundary as available at a Peeling iteration when every segment end followed by another domain/segment corresponds to the end of a PU in that iteration. `first_appearance_level` is the first zero-based matching iteration; `persistence_levels` is the number of matching iterations from first appearance through the final iteration. If no iteration matches, context validation fails rather than inventing zero provenance.

Add `CandidateLattice::attach_hierarchy(&mut self, provenance: &[MeasureProvenance], iterations: &[IterationResult]) -> Result<(), FeatureError>`. It indexes aligned provenance by each candidate's `source_index`, computes the two merge margins/path count, computes Peeling appearance/persistence, and replaces every `CandidateRecord::hierarchy` `None` with `Some(HierarchyEvidence)`. Feature extraction rejects a remaining `None` whenever the relative/hierarchy family is requested.

- [ ] **Step 5: Run Rust verification**

Run:

```bash
cargo fmt -- --check
cargo check
cargo test -p sword2-lib factorized_ranker --no-fail-fast
cargo test
```

Expected: all pass; a factorized-context cache miss is represented as an error value, not a panic.

- [ ] **Step 6: Commit evidence retention**

```bash
git add sword2-lib/src/sword/mod.rs sword2-lib/src/sword/factorized_ranker
git commit -m "feat: retain structural evidence for factorized selection"
```

---

### Task 4: Freeze the ordered feature schema and Python reference formulas

**Files:**
- Create: `benchmark/factorized_ranker/schema.py`
- Create: `benchmark/factorized_ranker/structure_features.py`
- Modify: `benchmark/requirements.txt`
- Modify: `benchmark/pyproject.toml`
- Test: `benchmark/tests/test_factorized_structure_features.py`

**Interfaces:**
- Consumes: `ValidatedPartition`, NumPy CA coordinates/contact matrix, and synthetic DSSP evidence.
- Produces: the exact ordered raw and pair feature names used by every later Python/Rust task.

- [ ] **Step 1: Write schema and exact-formula tests**

Create fixtures for: a symmetric octahedron; the sequence `HHHCCEEECC`; asymmetric two-domain points; a boundary with two retained crossing H-bonds; one bridge and one shared-sheet pair; and discontinuous partition `0-1;6-7 2-5`. Assert:

```python
def test_schema_is_unique_finite_and_has_no_competitor_fields():
    all_names = (*GLOBAL_FEATURES, *COUNT_ITEM_FEATURES, *CANDIDATE_FEATURES)
    assert len(all_names) == len(set(all_names))
    assert not any("merizo" in name.lower() or "chainsaw" in name.lower() for name in all_names)


def test_boundary_aggregation_is_exact():
    values = aggregate_measurements([{"hbond_count": 2.0}, {"hbond_count": 4.0}])
    assert values["boundary_hbond_count_min"] == 2.0
    assert values["boundary_hbond_count_mean"] == 3.0
    assert values["boundary_hbond_count_max"] == 4.0


def test_continuous_candidate_has_explicit_zero_discontinuity_vector():
    values = compute_discontinuity_features(continuous_fixture())
    assert values["has_discontinuity"] == 0.0
    assert all(values[name] == 0.0 for name in DISCONTINUITY_FEATURES[1:])
```

- [ ] **Step 2: Run tests and confirm the schema modules are missing**

Run: `benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_structure_features.py -q`

Expected: import failure.

- [ ] **Step 3: Add `scikit-learn` as an explicit offline dependency**

Append `scikit-learn` to `benchmark/requirements.txt` and to `[project].dependencies` in `benchmark/pyproject.toml`. Do not add it to either Rust crate.

- [ ] **Step 4: Define the schema programmatically but deterministically**

```python
# benchmark/factorized_ranker/schema.py
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
```

- [ ] **Step 5: Implement the Python reference formulas exactly**

Use these fixture/input records in `structure_features.py`:

```python
@dataclass(frozen=True)
class DsspHydrogenBond:
    partner: int
    energy_kcal: float


@dataclass(frozen=True)
class DsspResidue:
    index: int
    ss: Literal["helix", "strand", "coil"]
    bridge_partners: tuple[int, ...]
    sheet_label: str
    hydrogen_bonds: tuple[DsspHydrogenBond, ...]
    kappa: float
    alpha: float


@dataclass(frozen=True)
class PeelingSummary:
    n_levels: int
    finest_pu_count: int
```

Expose `compute_global_features(coordinates, contacts, dssp, peeling, candidate_counts)`, `compute_domain_conditioned_features(partition, coordinates, contacts)`, `compute_boundary_local_features(partition, coordinates, contacts, dssp)`, and `compute_discontinuity_features(partition, contacts, dssp)`, each returning a dictionary whose keys exactly equal its schema tuple.

Use `R_g / n^(1/3)`; sorted descending gyration eigenvalues; ratios `sqrt(lambda_2/lambda_1)` and `sqrt(lambda_3/lambda_1)` with zero when the denominator is zero. Map DSSP `H/G/I` to helix, `E/B` to strand, and every other code to coil. Compute nonlocal density over unordered pairs `i<j` with `j-i>=8`; divide contact order's contact-weighted separation by `n_residues`. Count maximal contiguous helix/strand runs.

For a domain `D`, define internal contact density as the mean `P(i,j)` over unordered nonlocal pairs in `D`; internal contact fraction as `internal_mass / (internal_mass + external_mass)`; conductance as `external_mass / (2*internal_mass + external_mass)`; and domain contact order as the contact-weighted separation divided by chain length. Use the existing `candidate_geometry.py` constants/formulas for q1/q2/q3, volume ratio, relative density, and contact-Q so the base 22 remain comparable.

For a sequential cut after residue `b`: an “inside” SSE is true only when residues `b` and `b+1` have the same class; otherwise `inside_coil=1`. The terminus distance is zero unless inside helix/strand, then the minimum residue distance from the cut to either end of that maximal run. Deduplicate retained crossing bonds by donor/acceptor; count bridge partners crossing; count nonblank equal sheet labels across the cut. For a half-window `w`, use truncated left `[b-w+1,b]`, right `[b+1,b+w]`, and `insulation = 1 - mean(P(left,right))`. Long-range contact density is the mean crossing contact over eligible pairs with separation `>=8`. `bend_change` is `abs(kappa_b-kappa_b+1)/180`; virtual-dihedral change is the wrapped absolute alpha difference divided by `180`, with DSSP sentinel `360` producing zero.

For a discontinuous segment, affinity is mean contact to other segments in its assigned domain; margin subtracts the largest mean contact to any competing domain. Long-range capture is same-domain long-range mass divided by all long-range mass incident to the segment. Interface-span entropy is normalized Shannon entropy across separation bins `8-15`, `16-31`, `32-63`, `64+`. Sheet links are deduplicated retained partners joining different segments of the same domain. Aggregate every boundary/segment measure as min/mean/max; use exactly zero only for a valid empty conditional set.

- [ ] **Step 6: Run focused and full tests**

Run:

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_structure_features.py -q
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: exact fixture assertions pass and every schema value is finite.

- [ ] **Step 7: Commit the frozen schema/reference checkpoint**

```bash
git add benchmark/factorized_ranker/schema.py \
  benchmark/factorized_ranker/structure_features.py \
  benchmark/requirements.txt benchmark/pyproject.toml \
  benchmark/tests/test_factorized_structure_features.py
git commit -m "feat: define factorized structural feature schema"
```

---

### Task 5: Rust global, count-independent candidate, and domain-conditioned features

**Files:**
- Create: `sword2-lib/src/sword/factorized_ranker/schema.rs`
- Create: `sword2-lib/src/sword/factorized_ranker/features.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/mod.rs`
- Modify: `sword2-lib/src/sword/geometry_metrics.rs`
- Test: inline tests in `factorized_ranker/features.rs`

**Interfaces:**
- Consumes: `CandidateLattice`, `ParsedPartition`, and `StructuralContext`.
- Produces: `GlobalFeatures`, `CandidateFeatures`, schema-ordered `Vec<f64>` values, and the reusable raw feature map consumed by Tasks 6–8.

- [ ] **Step 1: Write hand-computed Rust feature tests**

Port the same octahedron, asymmetric domains, and contact fixtures from Task 4. Include these assertions:

```rust
#[test]
fn global_contact_features_exclude_local_pairs() {
    let context = contact_fixture();
    let values = extract_global_features(&context, &[2, 2, 3]).unwrap();
    assert!((values.nonlocal_contact_density - 0.25).abs() < 1e-12);
    assert!((values.contact_order - 0.5).abs() < 1e-12);
}

#[test]
fn asymmetric_domain_features_preserve_smallest_and_largest() {
    let values = extract_candidate_base_and_domain(
        &fixture_candidate(), &fixture_context(), 2,
    ).unwrap();
    let smallest = values.base_and_domain_value("smallest_size_fraction").unwrap();
    let largest = values.base_and_domain_value("largest_size_fraction").unwrap();
    let q1_ratio = values.base_and_domain_value("smallest_to_largest_q1_ratio").unwrap();
    let smallest_q1 = values.base_and_domain_value("smallest_q1").unwrap();
    let largest_q1 = values.base_and_domain_value("largest_q1").unwrap();
    assert!(smallest < largest);
    assert!((q1_ratio - smallest_q1 / largest_q1).abs() < 1e-12);
}

#[test]
fn every_emitted_value_is_finite_and_schema_ordered() {
    let (global, candidate) = extract_fixture_features().unwrap();
    assert_eq!(global.to_vec().len(), GLOBAL_FEATURE_NAMES.len());
    assert_eq!(candidate.base_and_domain_vec().len(), BASE_AND_DOMAIN_FEATURE_NAMES.len());
    assert!(global.to_vec().iter().chain(candidate.base_and_domain_vec().iter()).all(|v| v.is_finite()));
}
```

- [ ] **Step 2: Run the focused Rust tests and confirm missing modules**

Run: `cargo test -p sword2-lib factorized_ranker::features -- --nocapture`

Expected: compile failure for missing feature structs/functions.

- [ ] **Step 3: Define named Rust records and one ordered conversion seam**

In `schema.rs`, set `pub(crate) const FEATURE_SCHEMA_VERSION: u32 = 1;`. Define named structs rather than using positional vectors internally:

```rust
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
pub(crate) struct CountFeatures { pub num_domains: usize, pub values: Vec<f64> }

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
        Self { global_count: true, domain_conditioned: true,
            boundary_local: true, relative_hierarchy: true,
            discontinuity: true }
    }
}

impl GlobalFeatures {
    pub(crate) fn to_vec(&self) -> Vec<f64> {
        let mut values = vec![
            self.n_residues, self.rg_normalized,
            self.inertia_ratio_21, self.inertia_ratio_31,
            self.nonlocal_contact_density, self.contact_order,
            self.helix_fraction, self.strand_fraction, self.coil_fraction,
            self.helix_blocks, self.strand_blocks,
            self.peeling_levels, self.finest_pus,
            self.candidate_total, self.available_count_total,
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
}

pub(crate) fn count_pair_vector(
    feature_names: &[&str],
    global: &GlobalFeatures,
    left: &CountFeatures,
    right: &CountFeatures,
) -> Result<Vec<f64>, FeatureError> {
    assemble_pair(feature_names, global, COUNT_ITEM_FEATURE_NAMES, &left.values, &right.values)
}

pub(crate) fn candidate_pair_vector(
    feature_names: &[&str],
    global: &GlobalFeatures,
    left: &CandidateFeatures,
    right: &CandidateFeatures,
) -> Result<Vec<f64>, FeatureError> {
    assemble_pair(feature_names, global, CANDIDATE_FEATURE_NAMES, &left.values, &right.values)
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
    let mut values = Vec::with_capacity(feature_names.len());
    for name in feature_names {
        let value = if let Some(raw) = name.strip_prefix("diff__") {
            let index = item_names.iter().position(|candidate| *candidate == raw)
                .ok_or(FeatureError::SchemaMismatch)?;
            left[index] - right[index]
        } else if let Some(raw) = name.strip_prefix("abs_diff__") {
            let index = item_names.iter().position(|candidate| *candidate == raw)
                .ok_or(FeatureError::SchemaMismatch)?;
            (left[index] - right[index]).abs()
        } else {
            let index = GLOBAL_FEATURE_NAMES.iter().position(|candidate| *candidate == *name)
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
```

Copy the generated names from Task 4 into `const` string slices and test them against the committed golden schema in Task 13. Pair assembly appends shared values, signed differences, then absolute differences, and rejects unequal lengths/non-finite values. Do not use a `HashMap` to determine order.

Every raw vector retains the full schema length. `SWORD2_DUMP_CANDIDATES` uses `FeatureMask::all()`. Final inference derives a mask from the generated model feature names; an unretained family is not computed and its raw slots are explicit zeros that the validated model cannot reference. This keeps rejected feature families from adding failure modes or runtime cost.

In this task, change context validation to `StructuralContext::validate(&self, mask: FeatureMask)`: dimensions and contact evidence are always required for the base model; DSSP details, merge provenance, and conditional extractors are required only when a retained feature name actually consumes them.

- [ ] **Step 4: Reuse geometry primitives without changing existing reports**

Make `centroid`, `gyration_tensor`, `sym3x3_eigenvalues`, `principal_radii`, and `radius_of_gyration` in `geometry_metrics.rs` `pub(crate)` where necessary. Leave existing `GeometryReport` formulas/output untouched. Implement schema-compatible q/contact formulas in `features.rs`; do not call permissive `parse_domain_indices` or the hard-cutoff `interdomain_contacts` helper.

- [ ] **Step 5: Implement global and per-domain formulas exactly as Task 4**

Build one `ContactFeatureCache` per chain from `ContactMatrix::get`: a flat 2-D prefix sum whose value is `P(i,j)` only when `|i-j|>=8`, plus a prefix sum of the corresponding eligible-pair indicators. Reuse `ContactMatrix::rectangle_sum` for all-contact rectangles and this cache for long-range rectangles. For one segment's self-rectangle, subtract diagonal mass and divide the symmetric remainder by two; for the nonlocal self-rectangle divide by two because its diagonal is already zero. For two distinct segments/domains, use one directed rectangle exactly once. Assemble domain/segment contact sums from these validated contiguous rectangles; do not recompute distances or scan the full residue matrix once per candidate.

Expose:

```rust
pub(crate) fn extract_global_features(
    context: &StructuralContext<'_>,
    candidate_counts: &[usize],
) -> Result<GlobalFeatures, FeatureError>;

pub(crate) fn extract_candidate_base_and_domain(
    candidate: &CandidateRecord,
    context: &StructuralContext<'_>,
    modal_count: usize,
) -> Result<CandidateFeatures, FeatureError>;
```

Use the Task 4 denominators and sequence cutoff verbatim. The histogram has normalized fractions for counts `1..=20` and one `21+` bin; modal-count ties choose the lower count. Smallest/largest domain ties choose the lower minimum residue. For a ratio with a zero denominator, emit `0.0`, not infinity. Validate every final field with `is_finite()` and return `FeatureError::NonFinite(the_schema_name)` on the first failure.

- [ ] **Step 6: Run parity-oriented Rust verification**

Run:

```bash
cargo fmt -- --check
cargo check
cargo test -p sword2-lib factorized_ranker::features -- --nocapture
cargo test
```

Expected: all pass; existing geometry report snapshots remain unchanged.

- [ ] **Step 7: Commit the core Rust features**

```bash
git add sword2-lib/src/sword/factorized_ranker \
  sword2-lib/src/sword/geometry_metrics.rs
git commit -m "feat: extract factorized chain and domain features"
```

---

### Task 6: Rust boundary-local and discontinuous segment features

**Files:**
- Create: `sword2-lib/src/sword/factorized_ranker/boundary.rs`
- Create: `sword2-lib/src/sword/factorized_ranker/discontinuity.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/features.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/schema.rs`
- Test: inline tests in both new modules

**Interfaces:**
- Consumes: validated partitions and complete `StructuralContext`.
- Produces: all 42 boundary fields and 16 discontinuity fields in Task 4 schema order.

- [ ] **Step 1: Write exact crossing-evidence tests**

```rust
fn named(values: &[f64], names: &[&str], name: &str) -> f64 {
    let index = names.iter().position(|candidate| *candidate == name).unwrap();
    values[index]
}

#[test]
fn crossing_bonds_are_deduplicated_and_converted_to_kcal() {
    let evidence = boundary_fixture_with_mirrored_donor_acceptor_records();
    let values = extract_boundary_features(&evidence.partition, &evidence.context).unwrap();
    assert_eq!(named(&values, BOUNDARY_LOCAL_FEATURE_NAMES, "boundary_hbond_count_mean"), 2.0);
    assert!((named(&values, BOUNDARY_LOCAL_FEATURE_NAMES,
        "boundary_hbond_energy_kcal_mean") - -2.4).abs() < 1e-12);
}

#[test]
fn insulation_uses_truncated_windows_and_contact_mean() {
    let values = extract_boundary_features(&edge_boundary_fixture(), &edge_context()).unwrap();
    assert!((named(&values, BOUNDARY_LOCAL_FEATURE_NAMES,
        "boundary_insulation_w8_mean") - 0.75).abs() < 1e-12);
}

#[test]
fn continuous_partition_has_explicit_zero_conditional_vector() {
    let values = extract_discontinuity_features(&continuous_partition(), &context()).unwrap();
    assert_eq!(values[0], 0.0);
    assert!(values[1..].iter().all(|value| *value == 0.0));
}
```

Also assert the discontinuous `0-1;6-7 2-5` fixture's affinity, competing-domain margin, long-range capture, four-bin entropy, and sheet-link count against Task 4's Python expected constants.

- [ ] **Step 2: Run tests and verify missing extractors**

Run:

```bash
cargo test -p sword2-lib factorized_ranker::boundary -- --nocapture
cargo test -p sword2-lib factorized_ranker::discontinuity -- --nocapture
```

Expected: compile failure.

- [ ] **Step 3: Implement sequential-boundary discovery and aggregation**

Sequential boundaries are cuts `b` for which `residue_to_domain[b] != residue_to_domain[b+1]`; do not treat a gap between two segments of the same discontinuous domain as a domain boundary. Implement a private `BoundaryEvidence` with the 14 Task 4 scalar measurements and:

```rust
pub(crate) fn extract_boundary_features(
    partition: &ParsedPartition,
    context: &StructuralContext<'_>,
) -> Result<Vec<f64>, FeatureError>;
```

Use `BTreeSet<(usize, usize)>` to deduplicate retained H-bonds, where the tuple is `(donor_clean_index, acceptor_clean_index)`. For boundary cut `b`, a retained pair crosses exactly when `min(i,j) <= b < max(i,j)`; the same definition applies to bridge and same-sheet pairs. Bridge partners are converted through the DSSP-to-clean mapping and deduplicated as sorted pairs. For same-sheet links, require equal nonblank `sheet_label` and count each residue pair once. Apply the exact windows, long-range cutoff, bend/dihedral normalization, and min/mean/max order from Task 4.

- [ ] **Step 4: Implement discontinuous segment association**

```rust
pub(crate) fn extract_discontinuity_features(
    partition: &ParsedPartition,
    context: &StructuralContext<'_>,
) -> Result<Vec<f64>, FeatureError>;
```

Only segments in domains with at least two segments enter the conditional aggregates. Contact means exclude self pairs and use all cross-segment pairs for affinity; the capture/entropy calculations use only separation `>=8`. Normalize entropy by `ln(4)` and emit zero when total mass is zero. A valid continuous candidate returns `[0.0; 16]`; a discontinuous candidate with no eligible long-range mass still sets `has_discontinuity=1.0` and returns finite zeros for those individual measures.

- [ ] **Step 5: Append the vectors in schema order and reject missing boundaries**

`extract_candidate_base_and_domain` must append boundary then discontinuity values only after their independent extractors succeed. A multi-domain candidate with no sequential boundary is valid only if all its domains are discontinuous; its boundary aggregate is a finite zero vector. Do not partially score a candidate after one extractor fails.

- [ ] **Step 6: Run Rust verification**

Run:

```bash
cargo fmt -- --check
cargo check
cargo test -p sword2-lib factorized_ranker --no-fail-fast
cargo test
```

Expected: all boundary/discontinuity fixture tests pass.

- [ ] **Step 7: Commit boundary and discontinuity features**

```bash
git add sword2-lib/src/sword/factorized_ranker
git commit -m "feat: add boundary and discontinuity evidence"
```

---

### Task 7: Sibling-relative, hierarchy, and per-count summary features in Python and Rust

**Files:**
- Create: `benchmark/factorized_ranker/lattice_features.py`
- Create: `benchmark/tests/test_factorized_lattice_features.py`
- Generate: `benchmark/fixtures/factorized_features_synthetic.json`
- Modify: `sword2-lib/src/sword/factorized_ranker/lattice.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/features.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/schema.rs`
- Test: inline Rust tests in `lattice.rs` and `features.rs`

**Interfaces:**
- Consumes: base candidate values, provenance retained in Task 2, and same-count sibling groups.
- Produces: complete schema-ordered `CountFeatures` and `CandidateFeatures` used for dumping, training, and inference.

- [ ] **Step 1: Write permutation, ties, and neighbor tests in Python**

```python
def test_sibling_percentile_uses_midrank_and_is_permutation_invariant():
    rows = sibling_rows(values=[1.0, 1.0, 3.0])
    first = add_sibling_relative_features(rows)
    second = add_sibling_relative_features(list(reversed(rows)))
    assert canonical_json(first) == canonical_json(second)
    assert [r["sibling_percentile_min_size"] for r in first] == [1/3, 1/3, 5/6]


def test_nonconsecutive_neighbor_deltas_include_gaps():
    counts = build_count_rows(global_fixture(), candidate_rows_for_counts(2, 4, 7))
    four = next(row for row in counts if row["count_num_domains"] == 4)
    assert four["count_lower_gap"] == 2.0
    assert four["count_higher_gap"] == 3.0
    assert four["count_max_cr_mean_delta_lower"] == (
        four["count_max_cr_mean"] - count_row(counts, 2)["count_max_cr_mean"]
    )
```

Add Rust equivalents and a test that nearest-sibling choice is invariant under lattice order.

- [ ] **Step 2: Run focused tests and confirm the modules/functions are absent**

Run:

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_lattice_features.py -q
cargo test -p sword2-lib factorized_ranker::features -- --nocapture
```

Expected: failures for missing lattice feature helpers.

- [ ] **Step 3: Implement sibling-relative features identically in both languages**

Expose Python:

```python
def add_sibling_relative_features(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]
```

and Rust:

```rust
pub(crate) fn add_sibling_and_hierarchy_features(
    lattice: &CandidateLattice,
    features: &mut [CandidateFeatures],
) -> Result<(), FeatureError>;
```

Group by domain count. For each `RELATIVE_CORE_FEATURE`, define percentile as `(number_less + 0.5*number_equal)/sibling_count` using exact `f64` equality on dump-derived values, and median delta as `value - sorted_midpoint_median`. Sort output candidates by `(num_domains, canonical)` before writing or comparing.

Choose the nearest sibling by minimizing the symmetric mean nearest-boundary distance between the two sorted sequential boundary sets; ties use canonical text. Emit signed `candidate.max_cr - sibling.max_cr` and `candidate.density_min - sibling.density_min`, or zero when there is no sibling. Append the seven hierarchy fields from `CandidateRecord::hierarchy`; no sign or monotonic transform is applied.

- [ ] **Step 4: Build exact per-count summaries and neighbor deltas**

Expose Python `build_count_rows(chain_features, candidate_rows)` and Rust:

```rust
pub(crate) fn extract_count_features(
    global: &GlobalFeatures,
    candidates: &[CandidateFeatures],
) -> Result<Vec<CountFeatures>, FeatureError>;
```

For each available count, compute min/arithmetic mean/max for the ten `COUNT_SUMMARY_SOURCES`; `count_candidate_fraction` divides by total lattice candidates. The closest lower/higher numeric counts are neighbors. Missing neighbor deltas and gaps are zero with the corresponding `has_*` indicator zero; present neighbor deltas are current minus neighbor and gaps are positive numeric differences. Return count rows sorted numerically.

- [ ] **Step 5: Compare Python and Rust fixture values**

Add a small canonical JSON fixture under `benchmark/fixtures/factorized_features_synthetic.json`, produced once from the literal synthetic inputs in the Python reference test and committed. Normal tests read rather than rewrite it: Python recomputes and compares its bytes/values, while the Rust test in `features.rs` loads it with `include_str!("../../../../benchmark/fixtures/factorized_features_synthetic.json")` and asserts name equality plus `abs(delta) <= 1e-12` for every value. This is formula parity; model parity is added after freezing.

- [ ] **Step 6: Run full verification**

Run:

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_lattice_features.py -q
benchmark/.venv/bin/python -m pytest benchmark/tests -q
cargo fmt -- --check
cargo check
cargo test
```

Expected: all tests pass and the synthetic vectors agree within `1e-12`.

- [ ] **Step 7: Commit complete feature assembly**

```bash
git add benchmark/factorized_ranker/lattice_features.py \
  benchmark/tests/test_factorized_lattice_features.py \
  benchmark/fixtures/factorized_features_synthetic.json \
  sword2-lib/src/sword/factorized_ranker
git commit -m "feat: assemble lattice-relative ranker features"
```

---

### Task 8: Exact Rust feature dump and deterministic normalized corpus

**Files:**
- Modify: `sword2-lib/src/sword/mod.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/mod.rs`
- Modify: `benchmark/dump_candidate_corpus.py`
- Create: `benchmark/factorized_ranker/corpus.py`
- Create: `benchmark/build_factorized_corpus.py`
- Create: `benchmark/tests/test_factorized_corpus.py`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: complete Rust raw features and Task 1 labels/rejections.
- Produces: `CorpusPaths` with `chains.csv`, `counts.csv`, `candidates.csv`, `rejections.csv`, and `corpus_manifest.json`; Tasks 9–13 consume these files/hashes.

- [ ] **Step 1: Write deterministic-corpus tests**

```python
def test_shuffled_input_produces_identical_corpus_bytes(tmp_path):
    first = build_fixture_corpus(tmp_path / "first", rows=fixture_rows())
    second = build_fixture_corpus(tmp_path / "second", rows=list(reversed(fixture_rows())))
    for name in ("chains", "counts", "candidates", "rejections", "manifest"):
        assert getattr(first, name).read_bytes() == getattr(second, name).read_bytes()


def test_manifest_changes_when_one_feature_changes(tmp_path):
    a = build_fixture_corpus(tmp_path / "a", rows=fixture_rows())
    b = build_fixture_corpus(tmp_path / "b", rows=mutate_one_feature(fixture_rows()))
    assert manifest_hash(a.manifest) != manifest_hash(b.manifest)


def test_schema_rejects_external_predictor_columns():
    with pytest.raises(ValueError, match="external predictor"):
        validate_dump_header([*required_dump_fields(), "merizo_score"])
```

- [ ] **Step 2: Run tests and confirm corpus tooling is missing**

Run: `benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_corpus.py -q`

Expected: import failure.

- [ ] **Step 3: Replace the old shortlist dump with the exact typed-lattice dump**

Keep `SWORD2_DUMP_CANDIDATES` as the opt-in environment variable, but route it through:

```rust
pub(crate) fn write_feature_dump(
    path: &std::path::Path,
    chain_id: &str,
    global: &GlobalFeatures,
    counts: &[CountFeatures],
    candidates: &[CandidateFeatures],
) -> Result<(), FactorizedError>;
```

Write one CSV row per candidate, repeat the small global/count vectors, and include `chain_id`, `canonical_delineation`, `source_index`, `num_domains`, `legacy_distance`, and every ordered schema value. Use Rust's round-trip `f64::to_string`; quote canonical delineations by doubling embedded quotes. The header comes from schema constants. If a chain cannot produce a complete feature lattice, write no rows and make `dump_candidate_corpus.py` record the chain failure. Do not run energy scoring or Python geometry enrichment on this path.

The dump must occur from the first-pass bounded lattice defined globally, not `relevant_measure2`. Add an integration test whose legacy-selected count is 2 but whose dump contains available 3- and 4-domain groups.

- [ ] **Step 4: Simplify `dump_candidate_corpus.py` to validate and preserve Rust features**

Remove calls that recompute candidate geometry. `_run_one` validates the exact schema header and copies rows into per-chain parts sorted by `(num_domains, canonical_delineation)`. `_merge_parts` remains resumable, but writes a failure CSV sorted by entry ID and reason. The temporary output stays alive until the Rust dump and canonical chain PDB have been consumed.

- [ ] **Step 5: Implement normalized corpus files and manifest**

```python
@dataclass(frozen=True)
class CorpusPaths:
    chains: Path
    counts: Path
    candidates: Path
    rejections: Path
    manifest: Path
```

Implement `write_corpus(paths, chain_rows, count_rows, candidate_rows, rejections, provenance) -> dict[str, object]`. It sorts and writes all four tables, computes their SHA-256 values from the written bytes, builds the canonical manifest dictionary, writes that dictionary, and returns the same dictionary.

Use `chains.csv` for one global row per chain, `counts.csv` for one row per `(chain_id,count)`, and `candidates.csv` for one row per canonical candidate plus CATH metrics/labels. Candidate ID is `sha256(chain_id + "\0" + canonical_delineation).hexdigest()`. Sort keys exactly as described in File Structure; write `\n` newlines and floats with `.17g`.

Before normalization, require all repeated global fields for a chain to be byte-identical after `.17g` formatting and all repeated count fields for a `(chain,count)` group to be identical. Reject the chain as a schema-integrity failure if they disagree; never choose the first conflicting row.

The canonical JSON manifest uses sorted keys, separators `(',', ':')`, and no timestamp or absolute path. Record dataset name/hash, binary hash, Git commit, exact argv, seed, thread/job counts, feature-schema hash, each CSV hash/row count, accepted/rejected chain counts, candidate rejection counts by code, and output hashes. `build_factorized_corpus.py` accepts a required mutually exclusive `--dump` or `--dump-dir`, plus `--dataset cath17287`, `--chain-cache-dir`, `--out-dir`, `--binary`, and `--rejections`; reject any dataset other than CATH-17287 in the training path.

- [ ] **Step 6: Update ignore rules narrowly**

Keep bulk corpora ignored. Add `benchmark/models/` as the committed artifact directory and explicit ignores for `benchmark/factorized-corpus*/` and `benchmark/data/*factorized*.csv`; later tasks stage only small JSON/golden/model files under `benchmark/models`.

- [ ] **Step 7: Run a two-chain smoke corpus and all tests**

Run:

```bash
cargo build --release
benchmark/.venv/bin/python -m benchmark.dump_candidate_corpus \
  --dataset cath17287 --limit 2 --jobs 1 --threads 1 \
  --out /tmp/sword2-factorized-smoke.csv
benchmark/.venv/bin/python -m benchmark.build_factorized_corpus \
  --dataset cath17287 \
  --dump /tmp/sword2-factorized-smoke.csv \
  --out-dir /tmp/sword2-factorized-corpus-smoke
benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_corpus.py -q
benchmark/.venv/bin/python -m pytest benchmark/tests -q
cargo check
cargo test
```

Expected: two-chain manifest is deterministic across two runs; all tests pass. If one of the first two dataset chains is a known malformed entry, the manifest instead deterministically records its rejection and the command still completes.

- [ ] **Step 8: Commit the reproducible corpus pipeline**

```bash
git add .gitignore benchmark/dump_candidate_corpus.py \
  benchmark/build_factorized_corpus.py benchmark/factorized_ranker/corpus.py \
  benchmark/tests/test_factorized_corpus.py \
  sword2-lib/src/sword/mod.rs sword2-lib/src/sword/factorized_ranker
git commit -m "feat: build deterministic factorized feature corpus"
```

---

### Task 9: Leakage-resistant five-fold manifest and paired-bootstrap primitive

**Files:**
- Modify: `benchmark/datasets.py`
- Create: `benchmark/factorized_ranker/folds.py`
- Create: `benchmark/build_factorized_folds.py`
- Modify: `benchmark/stats.py`
- Test: `benchmark/tests/test_factorized_folds.py`
- Create: `benchmark/tests/test_stats.py`

**Interfaces:**
- Consumes: CATH-17287 `CathEntry` metadata and the corpus manifest hash.
- Produces: deterministic `FoldAssignment` records, a canonical fold-manifest hash, seen/unseen individual-label cohorts, and `paired_chain_bootstrap`; Tasks 11, 15, and 16 depend on them.

- [ ] **Step 1: Write family-combination and connected-component tests**

```python
def test_family_combination_is_sorted_multiset_without_sentinel():
    chopping = "1-10:2.40.50.140|11-20:1.10.8.10|21-30:2.40.50.140|31-40:999_999"
    assert cath_family_combination(chopping) == (
        "1.10.8.10", "2.40.50.140", "2.40.50.140"
    )


def test_transitive_pdb_and_combination_links_share_component():
    # A shares PDB with B; B shares exact combination with C.
    components = connected_components(transitive_entries())
    assert component_ids(components, "A", "B", "C") == {component_id(components, "A")}


def test_no_pdb_or_exact_combination_crosses_folds():
    assignments = assign_folds(synthetic_entries(), n_folds=5, seed=37)
    validate_folds(assignments)
    assert len({assignment.fold for assignment in assignments}) == 5
```

Add a shuffle-invariance byte comparison and a sentinel-only test proving unknown labels receive entry-specific component keys instead of joining every unknown chain.

- [ ] **Step 2: Write paired bootstrap tests**

```python
def test_paired_bootstrap_uses_only_common_chain_ids():
    result = paired_chain_bootstrap(
        {"a": 0.7, "b": 0.9, "only_base": 1.0},
        {"a": 0.8, "b": 1.0, "only_new": 0.0},
        n_resamples=1000, seed=37,
    )
    assert result.n == 2
    assert result.mean_delta == pytest.approx(0.1)


def test_paired_bootstrap_is_deterministic():
    assert bootstrap_fixture() == bootstrap_fixture()
```

- [ ] **Step 3: Run tests and confirm the new APIs are absent**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_folds.py benchmark/tests/test_stats.py -q
```

Expected: import failures.

- [ ] **Step 4: Parse CATH family labels without losing multiplicity**

Add to `datasets.py`:

```python
def cath_family_labels(chopping: str) -> tuple[str, ...]:
    labels = []
    for raw_domain in chopping.split("|"):
        if ":" not in raw_domain:
            continue
        label = raw_domain.rsplit(":", 1)[1].strip()
        if label and label != "999_999":
            labels.append(label)
    return tuple(labels)


def cath_family_combination(chopping: str) -> tuple[str, ...]:
    return tuple(sorted(cath_family_labels(chopping)))
```

- [ ] **Step 5: Implement deterministic union-find and fold balancing**

```python
@dataclass(frozen=True)
class FoldAssignment:
    chain_id: str
    pdb_id: str
    family_combination: tuple[str, ...]
    family_labels: tuple[str, ...]
    true_count_bin: Literal["2", "3", "4", "5+"]
    length_bin: Literal["<250", "250-349", "350-449", "450+"]
    component_id: str
    fold: int
```

Implement these exact public functions: `connected_components(entries) -> list[Sequence[CathEntry]]`, `assign_folds(entries, n_folds=5, seed=37) -> Sequence[FoldAssignment]`, `validate_folds(assignments) -> None`, and `write_fold_manifest(path, assignments, dataset_sha256, corpus_sha256, seed=37) -> str`.

Union by normalized four-character `pdb_id` and nonempty exact combination. Component ID is SHA-256 of sorted entry IDs joined by NUL. Order components by descending chain count, descending largest stratum count, then SHA-256 of `f"{seed}:{component_id}"`. For each candidate fold, calculate the sum of squared `(new_count-target)/max(target,1)` deviations for total chains, four count bins, and four length bins. Choose minimum cost, then lower current fold size, then lower fold number. Persist sorted records plus source hashes and the canonical manifest hash. Validate chain uniqueness and assert each PDB/combination maps to one fold.

Expose `individual_label_seen(assignments, validation_fold) -> dict[str, bool]` using the other four folds as the seen-label set.

- [ ] **Step 6: Implement the paired chain bootstrap in `stats.py`**

```python
@dataclass(frozen=True)
class PairedBootstrap:
    n: int
    mean_delta: float
    ci_low: float
    ci_high: float

def paired_chain_bootstrap(
    baseline_by_chain: Mapping[str, float],
    experiment_by_chain: Mapping[str, float],
    *, n_resamples: int = 10_000, seed: int = 37,
    higher_is_better: bool = True,
) -> PairedBootstrap:
    ids = sorted(set(baseline_by_chain) & set(experiment_by_chain))
    if not ids:
        raise ValueError("no common chains for paired bootstrap")
    baseline = np.asarray([baseline_by_chain[key] for key in ids], dtype=float)
    experiment = np.asarray([experiment_by_chain[key] for key in ids], dtype=float)
    delta = experiment - baseline if higher_is_better else baseline - experiment
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(ids), size=(n_resamples, len(ids)))
    means = delta[draws].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return PairedBootstrap(len(ids), float(delta.mean()), float(low), float(high))
```

- [ ] **Step 7: Run tests and build a smoke fold manifest**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_folds.py benchmark/tests/test_stats.py -q
benchmark/.venv/bin/python -m benchmark.build_factorized_folds \
  --dataset cath17287 \
  --corpus-manifest /tmp/sword2-factorized-corpus-smoke/corpus_manifest.json \
  --out /tmp/cath17287-factorized-folds-smoke.json
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: zero leakage assertions and byte-identical output on a repeated command.

- [ ] **Step 8: Commit fold/statistics tooling**

```bash
git add benchmark/datasets.py benchmark/stats.py \
  benchmark/factorized_ranker/folds.py benchmark/build_factorized_folds.py \
  benchmark/tests/test_factorized_folds.py benchmark/tests/test_stats.py
git commit -m "feat: add leakage-resistant structural validation folds"
```

---

### Task 10: Deterministic pair construction and reference Borda selection

**Files:**
- Create: `benchmark/factorized_ranker/pairs.py`
- Create: `benchmark/factorized_ranker/ranking.py`
- Test: `benchmark/tests/test_factorized_pairs.py`
- Test: `benchmark/tests/test_factorized_ranking.py`

**Interfaces:**
- Consumes: normalized chain/count/candidate tables and Task 4 schema names.
- Produces: `PairBatch`, pair feature names, symmetrized probabilities, Borda scores, and exact tie-break functions used by training/evaluation/export parity.

- [ ] **Step 1: Write pair symmetry, weighting, and cap tests**

```python
def test_mirrored_pair_has_reversed_diff_and_same_absolute_part():
    batch = build_candidate_pairs(*two_candidate_fixture())
    assert batch.y.tolist() == [1, 0]
    shared = len(GLOBAL_FEATURES)
    item = len(candidate_item_features(("base",)))
    np.testing.assert_allclose(batch.x[0, :shared], batch.x[1, :shared])
    np.testing.assert_allclose(batch.x[0, shared:shared+item], -batch.x[1, shared:shared+item])
    np.testing.assert_allclose(batch.x[0, shared+item:], batch.x[1, shared+item:])


def test_pair_cap_and_chain_weights_are_exact():
    batch = build_candidate_pairs(*large_fixture(), max_unordered_pairs=64, seed=37)
    assert len(batch.y) == 128
    for chain_id in np.unique(batch.chain_ids):
        assert batch.sample_weight[batch.chain_ids == chain_id].sum() == pytest.approx(1.0)
```

Add tests that count pairs skip equal absolute true-count distance, candidate pairs do not cross count/chain and skip NDO ties, quartiles are represented when nonempty, and row permutation yields byte-equivalent arrays/IDs.

- [ ] **Step 2: Write Borda orientation and tie tests**

```python
def test_symmetrization_calls_both_orientations_once():
    calls = []
    probability = symmetrized_probability(recording_predictor(calls), shared(), left(), right())
    assert len(calls) == 2
    assert probability == pytest.approx(0.5 * (calls[0].prob + 1.0 - calls[1].prob))


def test_ties_follow_approved_rules():
    assert select_count({2: 0.5, 3: 0.5}, legacy_count=3) == 3
    assert select_count({2: 0.5, 3: 0.5}, legacy_count=4) == 2
    assert select_candidate({"0-4 5-9": 0.5, "0-5 6-9": 0.5}) == "0-4 5-9"
```

- [ ] **Step 3: Run tests and verify missing modules**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_pairs.py benchmark/tests/test_factorized_ranking.py -q
```

Expected: import failures.

- [ ] **Step 4: Implement stable, weighted pair construction**

```python
@dataclass(frozen=True)
class PairBatch:
    feature_names: tuple[str, ...]
    x: np.ndarray
    y: np.ndarray
    sample_weight: np.ndarray
    chain_ids: np.ndarray
    left_ids: np.ndarray
    right_ids: np.ndarray
```

Implement `build_count_pairs(chains, counts, *, shared_features=GLOBAL_FEATURES, item_features=COUNT_ITEM_FEATURES, seed=37, max_unordered_pairs=64) -> PairBatch` and the corresponding `build_candidate_pairs(..., item_features=CANDIDATE_FEATURES, ...) -> PairBatch`.

Deduplicate candidates by `(chain_id, canonical_delineation)` before pairing. A chain contributes count-head pairs only when its exact true count exists among its available count rows. Stable pair ID is the sorted concatenation of item IDs. Difference target for count sampling is `abs(abs(k_i-true)-abs(k_j-true))`; for candidates it is `abs(ndo_i-ndo_j)`. Assign quartiles with `pandas.qcut(..., q=4, duplicates="drop")`; allocate `64 // n_nonempty` per bin, distribute the remainder from largest bins then lexical stable IDs, and fill unused capacity from remaining pairs sorted by SHA-256 of `f"{seed}:{chain}:{pair_id}"`. Mirror only after selecting unordered pairs. Each ordered row weight is `1/(2*n_selected_for_chain)`.

- [ ] **Step 5: Implement reference ranking**

Implement `symmetrized_probability(predict_proba, shared, left, right) -> float`, `normalized_borda(item_ids, probability) -> dict[object, float]`, `select_count(scores, legacy_count) -> int`, and `select_candidate(scores) -> str` as the shared reference-selection API.

For one item, Borda is `1.0` without model invocation. Otherwise its score is the arithmetic mean of symmetrized win probability against every other item, always iterating opponents in canonical item-ID order. Validate finite `[0,1]` probabilities. Apply a tie rule only when the resulting `f64` scores are exactly equal; canonical iteration makes this enumeration-independent in both Python and Rust.

- [ ] **Step 6: Run focused and full tests**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_pairs.py benchmark/tests/test_factorized_ranking.py -q
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: all pass.

- [ ] **Step 7: Commit pair/ranking primitives**

```bash
git add benchmark/factorized_ranker/pairs.py \
  benchmark/factorized_ranker/ranking.py \
  benchmark/tests/test_factorized_pairs.py \
  benchmark/tests/test_factorized_ranking.py
git commit -m "feat: add deterministic factorized pair ranking"
```

---

### Task 11: Grouped CV, fixed-grid selection, and ordered feature-family ablations

**Files:**
- Create: `benchmark/factorized_ranker/training.py`
- Create: `benchmark/train_factorized_ranker.py`
- Test: `benchmark/tests/test_factorized_training.py`
- Create: `benchmark/tests/fixtures/factorized_corpus/chains.csv`
- Create: `benchmark/tests/fixtures/factorized_corpus/counts.csv`
- Create: `benchmark/tests/fixtures/factorized_corpus/candidates.csv`
- Create: `benchmark/tests/fixtures/factorized_corpus/rejections.csv`
- Create: `benchmark/tests/fixtures/factorized_corpus/corpus_manifest.json`
- Create: `benchmark/tests/fixtures/factorized_folds.json`

**Interfaces:**
- Consumes: `CorpusPaths`, fold manifest, pair/ranking APIs, and paired bootstrap.
- Produces: OOF chain predictions, fold/cohort metrics, selected head hyperparameters, retained feature families, and final fitted sklearn heads for Task 12.

- [ ] **Step 1: Write grid, leakage, macro-metric, and ablation-gate tests**

```python
def test_model_grid_is_exact_and_bounded():
    assert len(MODEL_GRID) == 8
    assert {(p.n_estimators, p.learning_rate, p.min_samples_leaf, p.max_depth)
            for p in MODEL_GRID} == {
        (n, lr, leaf, 3)
        for n in (64, 96) for lr in (0.03, 0.05) for leaf in (32, 64)
    }


def test_training_never_sees_validation_chain(monkeypatch):
    seen = capture_fit_chain_ids(monkeypatch)
    cross_validate_configuration(synthetic_corpus(), synthetic_folds(), params(), params(), ("base",))
    assert all(fit.isdisjoint(validation) for fit, validation in seen)


def test_ablation_guards_are_exact():
    assert retain_ablation(report(ndo_delta=0.01, ci_low=0.0, count_delta=0.0))
    assert not retain_ablation(report(ndo_delta=0.01, ci_low=-1e-6, count_delta=0.0))
    assert not retain_ablation(report(ndo_delta=0.01, ci_low=0.0, count_delta=-1e-6,
                                      family="global_count"))
    assert not retain_ablation(report(ndo_delta=0.01, ci_low=0.0,
                                      contiguous_delta=-0.0051,
                                      family="discontinuity"))
```

Also test deterministic reports, chain-macro rather than candidate-row metrics, validation cohort reporting, final all-fold refit, and CLI rejection of a `--test`/CATH-663 argument.

- [ ] **Step 2: Run tests and confirm training module is absent**

Run: `benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_training.py -q`

Expected: import failure.

- [ ] **Step 3: Define the only legal model grid and fitting function**

```python
@dataclass(frozen=True, order=True)
class Hyperparameters:
    n_estimators: Literal[64, 96]
    learning_rate: Literal[0.03, 0.05]
    min_samples_leaf: Literal[32, 64]
    max_depth: Literal[3] = 3

MODEL_GRID = tuple(
    Hyperparameters(n, rate, leaf)
    for n in (64, 96) for rate in (0.03, 0.05) for leaf in (32, 64)
)

def fit_head(batch: PairBatch, params: Hyperparameters, seed: int = 37):
    model = GradientBoostingClassifier(
        n_estimators=params.n_estimators,
        learning_rate=params.learning_rate,
        min_samples_leaf=params.min_samples_leaf,
        max_depth=3, random_state=seed, loss="log_loss",
    )
    return model.fit(batch.x, batch.y, sample_weight=batch.sample_weight)
```

No scaling, imputation, probability calibration, early stopping, or alternate model family is permitted.

- [ ] **Step 4: Select head hyperparameters without coupling the grids**

For each fold and grid entry, fit on four folds only. Select count-head parameters by highest chain-macro count accuracy, then lowest chain-macro absolute count error. Select candidate-head parameters by averaging selected within-count NDO across counts inside each chain and then taking the macro mean over chains. Both ties prefer fewer trees, then larger leaf size, then lower learning rate. This avoids a 64-way coupled search and keeps the two factorized objectives separately auditable.

After choosing one configuration per head, run end-to-end count then candidate Borda to produce one OOF winner per chain. Report mean NDO, count accuracy, BF1@10, matched Dice, predicted-domain mean, total/count/within-count regret, each fold, count bins `2/3/4/5+`, length bins, contiguous/discontinuous cohorts, and individual-family seen/unseen cohorts.

- [ ] **Step 5: Implement ordered ablations exactly**

Use these family stages:

```python
FEATURE_FAMILY_ORDER = (
    "base",
    "global_count",
    "domain_conditioned",
    "boundary_local",
    "relative_hierarchy",
    "discontinuity",
)
```

`base` means the candidate head's current 22 energy-free fields and a minimal count head using `count_num_domains`, candidate count/fraction/modal distance, and summaries derived from those 22 fields; it has no new global structural context. `global_count` adds all `GLOBAL_FEATURES`, full per-count summaries, and neighbor deltas to the count head and shared globals to the candidate head. Later stages add their names from `schema.py` only to the candidate head.

Compare the proposed stage's OOF chain NDO to the last retained stage with 10,000 paired resamples, seed 37. Retain only if mean delta `>0` and CI low `>=0`. `global_count` must also have count-accuracy delta `>=0`. `discontinuity` must also have contiguous-chain NDO delta `>=-0.005`; apply the same `-0.005` no-regression diagnostic to every stage's contiguous and discontinuous cohorts and record any failure.

- [ ] **Step 6: Add a development-only CLI**

`benchmark/train_factorized_ranker.py` accepts `--corpus-dir`, `--fold-manifest`, `--out-dir`, and `--seed` (must equal 37 for a freeze). It rejects corpus manifests whose dataset is not `cath17287`, verifies all hashes before reading CSVs, has no CATH-663/test option, writes canonical `oof_predictions.csv`, `cv_report.json`, and `ablation_report.json`, then refits the selected heads on all development rows in memory for artifact freezing.

- [ ] **Step 7: Run synthetic CV twice and all Python tests**

Run:

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_training.py -q
benchmark/.venv/bin/python -m benchmark.train_factorized_ranker \
  --corpus-dir benchmark/tests/fixtures/factorized_corpus \
  --fold-manifest benchmark/tests/fixtures/factorized_folds.json \
  --out-dir /tmp/factorized-training-a --seed 37
benchmark/.venv/bin/python -m benchmark.train_factorized_ranker \
  --corpus-dir benchmark/tests/fixtures/factorized_corpus \
  --fold-manifest benchmark/tests/fixtures/factorized_folds.json \
  --out-dir /tmp/factorized-training-b --seed 37
diff -ru /tmp/factorized-training-a /tmp/factorized-training-b
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: `diff` is empty and all tests pass.

- [ ] **Step 8: Commit grouped training and ablations**

```bash
git add benchmark/factorized_ranker/training.py \
  benchmark/train_factorized_ranker.py \
  benchmark/tests/test_factorized_training.py \
  benchmark/tests/fixtures/factorized_corpus \
  benchmark/tests/fixtures/factorized_folds.json
git commit -m "feat: train factorized heads with grouped ablations"
```

---

### Task 12: Versioned model artifacts, deterministic Rust exporter, and golden vectors

**Files:**
- Create: `benchmark/factorized_ranker/model_artifact.py`
- Create: `benchmark/export_factorized_ranker.py`
- Modify: `benchmark/train_factorized_ranker.py`
- Test: `benchmark/tests/test_factorized_model_artifact.py`
- Test: `benchmark/tests/test_export_factorized_ranker.py`

**Interfaces:**
- Consumes: fitted binary `GradientBoostingClassifier` heads and frozen provenance hashes.
- Produces: validated canonical JSON heads, reference probabilities, generated Rust source text, and golden vectors. Task 13 runs this on the real development models; Task 14 consumes its Rust output.

- [ ] **Step 1: Write sklearn parity and invalid-artifact tests**

```python
def test_frozen_tree_probability_matches_sklearn():
    model, x = fitted_tiny_classifier()
    artifact = freeze_classifier(model, head_name="count", feature_names=("x",), **hash_args())
    expected = model.predict_proba(x)[:, 1]
    actual = predict_artifact(artifact, x)
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-12)


def test_threshold_equality_takes_sklearn_left_branch():
    artifact = one_split_artifact(threshold=1.25, left=-1.0, right=1.0)
    assert raw_tree_value(artifact, np.array([1.25])) == -1.0


@pytest.mark.parametrize("mutation", [
    invalid_child, nan_threshold, nan_leaf, too_many_trees,
    too_deep, too_many_combined_nodes, wrong_feature_count,
])
def test_invalid_artifact_is_rejected(mutation):
    with pytest.raises(ValueError):
        validate_artifacts(*mutation(valid_artifacts()))
```

Add byte-determinism tests for JSON, Rust rendering, and golden vectors.

- [ ] **Step 2: Run tests and confirm artifact modules are absent**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_model_artifact.py \
  benchmark/tests/test_export_factorized_ranker.py -q
```

Expected: import failures.

- [ ] **Step 3: Define and validate the canonical artifact**

```python
@dataclass(frozen=True)
class TreeNode:
    feature_index: int       # -1 for a leaf
    threshold: float
    left_child: int
    right_child: int
    leaf_value: float

@dataclass(frozen=True)
class FrozenTree:
    nodes: tuple[TreeNode, ...]

```

Implement `freeze_classifier(model, head_name, feature_names, corpus_sha256, fold_manifest_sha256, retained_feature_families, training_command, seed=37) -> dict[str, object]`, `validate_artifacts(count_artifact, candidate_artifact) -> None`, `predict_artifact(artifact, x) -> np.ndarray`, and `write_artifact(path, artifact) -> str`.

Extract each `model.estimators_[stage, 0].tree_` in its native node order. A leaf has `feature_index=-1`, children `-1`, threshold `0.0`, and its exact `tree_.value[node,0,0]`. An internal node has leaf value `0.0`; prediction uses `<= threshold` for left. Obtain initial log odds by evaluating `model._raw_predict_init(np.zeros((1,n_features)))[0,0]` and store it explicitly. Probability is `sigmoid(initial_log_odds + learning_rate*sum(tree_leaf))` with a numerically stable sigmoid.

Artifact keys are: `schema_version`, `head`, `feature_names`, `feature_count`, `corpus_sha256`, `fold_manifest_sha256`, `retained_feature_families`, `training_command`, `seed`, `n_estimators`, `max_depth`, `learning_rate`, `initial_log_odds`, and `trees`. Validate head names, exact schema version, finite numeric values, node/child reachability, no cycles, depth/tree/node caps, ordered-feature uniqueness, both hashes, and combined caps.

- [ ] **Step 4: Render static Rust arrays deterministically**

Implement `render_rust_models(count_artifact, candidate_artifact) -> str` and `write_golden_vectors(path, artifacts, vectors) -> str` using the validated representation above.

Flatten each head's nodes globally and offset child indices/tree roots. Render `f64` using `repr(float(value))`, feature indices as `u16`, leaf feature as `u16::MAX`, and one `StaticTree { root }` per estimator. Emit `RETAINED_FEATURE_FAMILIES`, `COUNT_FEATURE_NAMES`, `CANDIDATE_FEATURE_NAMES`, `COUNT_NODES`, `COUNT_TREES`, `COUNT_MODEL`, then candidate equivalents. The generated file imports handwritten `StaticBoostedModel`, `StaticNode`, and `StaticTree` from `super::model`; it contains no traversal logic.

Golden JSON includes schema names, at least four vectors per head (below/equal/above thresholds), both pair orientations, Python raw/probability values, symmetrized values, Borda totals, winners, and malformed/non-finite fallback cases.

- [ ] **Step 5: Add the exporter CLI and formatting check**

Extend `benchmark/train_factorized_ranker.py` with required `--count-model-out` and `--candidate-model-out` arguments. After the all-development refit, call `freeze_classifier` and `write_artifact` for each head; include the exact invoked argv, corpus/fold hashes, retained families, and selected parameters. `benchmark/export_factorized_ranker.py` accepts both model paths, corpus/fold/CV/ablation/OOF/baseline provenance paths, `--rust-out`, `--golden-out`, and `--manifest-out`; it validates and atomically writes all three derived artifacts. A test writes to temp paths twice and compares bytes, then invokes `rustfmt --check` on the generated source. Neither CLI may read CATH-663 except to hash the already-frozen standalone baseline file supplied as opaque provenance.

- [ ] **Step 6: Run focused, full Python, and Rust formatting tests**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_model_artifact.py \
  benchmark/tests/test_export_factorized_ranker.py -q
benchmark/.venv/bin/python -m pytest benchmark/tests -q
cargo fmt -- --check
```

Expected: exact sklearn parity within `1e-12`; byte-identical exports.

- [ ] **Step 7: Commit artifact/export tooling only**

```bash
git add benchmark/factorized_ranker/model_artifact.py \
  benchmark/export_factorized_ranker.py benchmark/train_factorized_ranker.py \
  benchmark/tests/test_factorized_model_artifact.py \
  benchmark/tests/test_export_factorized_ranker.py
git commit -m "feat: export compact factorized tree artifacts"
```

---

### Task 13: Build the full development corpus, run grouped ablations, and freeze the real models

**Files:**
- Generate: `benchmark/models/cath17287_factorized_corpus_v1_manifest.json`
- Generate: `benchmark/models/cath17287_factorized_folds_v1.json`
- Generate: `benchmark/models/factorized_ranker_v1_cv.json`
- Generate: `benchmark/models/factorized_ranker_v1_ablations.json`
- Generate: `benchmark/models/factorized_ranker_v1_oof.csv`
- Generate: `benchmark/models/cath663_standalone_structural_baseline.csv`
- Generate: `benchmark/models/factorized_count_v1.json`
- Generate: `benchmark/models/factorized_candidate_v1.json`
- Generate: `benchmark/models/factorized_ranker_v1_manifest.json`
- Generate: `benchmark/models/factorized_ranker_v1_golden.json`
- Generate: `sword2-lib/src/sword/factorized_ranker/generated_model.rs`

**Interfaces:**
- Consumes: all development tooling from Tasks 1–12 and only CATH-17287/cache structures.
- Produces: immutable real model/provenance hashes and generated arrays. Task 14 must reproduce their goldens; Task 17's locked benchmark is forbidden until this task is committed.

- [ ] **Step 1: Record the clean execution inputs before the long run**

Run:

```bash
git status --short
cargo build --release
sha256sum target/release/sword2
benchmark/.venv/bin/python -c \
  'import sklearn,sys; print(sys.version); print(sklearn.__version__)'
```

Expected: record the exact feature-dump binary hash and versions in the corpus provenance. This is distinct from the final inference binary hash recorded by the locked run. Do not clean or reset unrelated worktree files.

Copy the already-existing, pre-specification standalone structural selection artifact without rerunning or selecting among alternatives:

```bash
cp benchmark/data/cath663_structural_ranker_full_selection.csv \
  benchmark/models/cath663_standalone_structural_baseline.csv
benchmark/.venv/bin/python -c \
  'import pandas as pd; p="benchmark/models/cath663_standalone_structural_baseline.csv"; d=pd.read_csv(p); assert d.chain_id.nunique()==663; print(len(d), d.chain_id.nunique())'
sha256sum benchmark/models/cath663_standalone_structural_baseline.csv
```

This artifact is used only for the predeclared contiguous/discontinuous no-regression gates. It is not a teacher, feature source, or model-selection input.

- [ ] **Step 2: Build/resume the full CATH-17287 feature dump**

Run:

```bash
benchmark/.venv/bin/python -m benchmark.dump_candidate_corpus \
  --dataset cath17287 \
  --out benchmark/data/cath17287_factorized_candidates_v1.csv \
  --parts-dir benchmark/data/cath17287_factorized_candidates_v1_parts \
  --chain-dir benchmark/cache/chains \
  --binary target/release/sword2 \
  --jobs 16 --threads 1 --timeout 300 --seed 37 --resume
```

Expected: every successful part has the exact schema header; CATH-663 IDs are excluded; failures are sorted and written separately. Resume until all non-failed entries have parts. Do not change feature/model choices based on failures—only Task 1 integrity rules may exclude chains.

- [ ] **Step 3: Label and normalize the full corpus**

Run:

```bash
benchmark/.venv/bin/python -m benchmark.build_factorized_corpus \
  --dataset cath17287 \
  --dump benchmark/data/cath17287_factorized_candidates_v1.csv \
  --chain-cache-dir benchmark/cache/chains \
  --out-dir benchmark/data/cath17287_factorized_corpus_v1 \
  --binary target/release/sword2 \
  --rejections benchmark/data/cath17287_factorized_corpus_v1/rejections.csv
cp benchmark/data/cath17287_factorized_corpus_v1/corpus_manifest.json \
  benchmark/models/cath17287_factorized_corpus_v1_manifest.json
```

Expected: the six known malformed chains and 26 corrupt candidate rows are rejected by codes rather than silently scored; record actual counts without hard-coding them into tests. Re-run the builder and verify the manifest hash is unchanged.

- [ ] **Step 4: Build and validate the real five-fold manifest**

Run:

```bash
benchmark/.venv/bin/python -m benchmark.build_factorized_folds \
  --dataset cath17287 \
  --corpus-manifest benchmark/models/cath17287_factorized_corpus_v1_manifest.json \
  --out benchmark/models/cath17287_factorized_folds_v1.json \
  --seed 37
benchmark/.venv/bin/python -c \
  'from pathlib import Path; from benchmark.factorized_ranker.folds import load_fold_manifest,validate_folds; validate_folds(load_fold_manifest(Path("benchmark/models/cath17287_factorized_folds_v1.json")))'
```

Expected: no four-character PDB ID or exact family combination crosses folds; all five folds are populated.

- [ ] **Step 5: Run the fixed-grid grouped CV and ordered ablations once**

Run:

```bash
benchmark/.venv/bin/python -m benchmark.train_factorized_ranker \
  --corpus-dir benchmark/data/cath17287_factorized_corpus_v1 \
  --fold-manifest benchmark/models/cath17287_factorized_folds_v1.json \
  --out-dir benchmark/data/factorized_ranker_v1_training \
  --count-model-out benchmark/models/factorized_count_v1.json \
  --candidate-model-out benchmark/models/factorized_candidate_v1.json \
  --seed 37
cp benchmark/data/factorized_ranker_v1_training/cv_report.json \
  benchmark/models/factorized_ranker_v1_cv.json
cp benchmark/data/factorized_ranker_v1_training/ablation_report.json \
  benchmark/models/factorized_ranker_v1_ablations.json
cp benchmark/data/factorized_ranker_v1_training/oof_predictions.csv \
  benchmark/models/factorized_ranker_v1_oof.csv
```

Expected: the report explicitly accepts/rejects every family in the approved order, includes CIs/guards/cohorts, and contains no CATH-663 metric. Do not override a rejected family manually.

- [ ] **Step 6: Export the frozen arrays and goldens**

Run:

```bash
benchmark/.venv/bin/python -m benchmark.export_factorized_ranker \
  --count-model benchmark/models/factorized_count_v1.json \
  --candidate-model benchmark/models/factorized_candidate_v1.json \
  --corpus-manifest benchmark/models/cath17287_factorized_corpus_v1_manifest.json \
  --fold-manifest benchmark/models/cath17287_factorized_folds_v1.json \
  --cv-report benchmark/models/factorized_ranker_v1_cv.json \
  --ablation-report benchmark/models/factorized_ranker_v1_ablations.json \
  --oof-predictions benchmark/models/factorized_ranker_v1_oof.csv \
  --standalone-baseline benchmark/models/cath663_standalone_structural_baseline.csv \
  --rust-out sword2-lib/src/sword/factorized_ranker/generated_model.rs \
  --golden-out benchmark/models/factorized_ranker_v1_golden.json \
  --manifest-out benchmark/models/factorized_ranker_v1_manifest.json
rustfmt --check sword2-lib/src/sword/factorized_ranker/generated_model.rs
```

Expected: combined trees/nodes remain within 192/2,880 and reference probabilities reproduce sklearn within `1e-12`.

- [ ] **Step 7: Verify the generated top-level immutable manifest**

Require `factorized_ranker_v1_manifest.json` to contain canonical JSON keys for schema version, seed, retained families, exact training/export commands, Git commit, feature-dump binary hash, corpus/fold hashes, both model hashes, CV/ablation/OOF/golden/generated-Rust hashes, the existing standalone CATH-663 baseline hash, Python/sklearn versions, tree/node counts, and ordered feature-schema hashes. Then run:

```bash
sha256sum benchmark/models/factorized_*_v1.json \
  benchmark/models/factorized_ranker_v1_oof.csv \
  benchmark/models/cath663_standalone_structural_baseline.csv \
  sword2-lib/src/sword/factorized_ranker/generated_model.rs
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: all manifest hashes match the files; Python tests pass.

- [ ] **Step 8: Commit the frozen development decision before any locked evaluation**

```bash
git add benchmark/models/cath17287_factorized_corpus_v1_manifest.json \
  benchmark/models/cath17287_factorized_folds_v1.json \
  benchmark/models/factorized_ranker_v1_cv.json \
  benchmark/models/factorized_ranker_v1_ablations.json \
  benchmark/models/factorized_ranker_v1_oof.csv \
  benchmark/models/cath663_standalone_structural_baseline.csv \
  benchmark/models/factorized_count_v1.json \
  benchmark/models/factorized_candidate_v1.json \
  benchmark/models/factorized_ranker_v1_manifest.json \
  benchmark/models/factorized_ranker_v1_golden.json \
  sword2-lib/src/sword/factorized_ranker/generated_model.rs
git commit -m "model: freeze factorized structural ranker v1"
```

Do not stage bulk `benchmark/data` CSVs/parts.

---

### Task 14: Static Rust tree inference and Python/Rust golden parity

**Files:**
- Create: `sword2-lib/src/sword/factorized_ranker/model.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/mod.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/generated_model.rs` only by rerunning the exporter if formatting requires it
- Create: `sword2-lib/tests/factorized_ranker_golden.rs`

**Interfaces:**
- Consumes: generated `COUNT_MODEL`, `CANDIDATE_MODEL`, complete feature vectors, and the golden JSON.
- Produces: validated probability evaluation, symmetrized pair probability, normalized Borda, and the pure `select_factorized` result used by Task 15.

- [ ] **Step 1: Write tree traversal, validation, and Borda tests**

```rust
#[test]
fn equality_takes_left_branch() {
    let model = one_split_model(1.25, -1.0, 1.0);
    assert_eq!(raw_tree_value(&model, &[1.25]).unwrap(), -1.0);
}

#[test]
fn validator_rejects_bad_schema_nonfinite_and_child_indices() {
    assert!(matches!(validate_model(&wrong_schema_model()), Err(ModelError::SchemaVersion { .. })));
    assert!(matches!(validate_model(&nan_model()), Err(ModelError::NonFiniteNode { .. })));
    assert!(matches!(validate_model(&bad_child_model()), Err(ModelError::ChildIndex { .. })));
}

#[test]
fn borda_is_permutation_invariant_and_ties_are_canonical() {
    let first = rank_candidates(&fixture_candidates_in_order(&[2, 0, 1])).unwrap();
    let second = rank_candidates(&fixture_candidates_in_order(&[1, 2, 0])).unwrap();
    assert_eq!(first.canonical, second.canonical);
}
```

Add model depth/tree/node cap tests, both-orientation call tests, count ties preferring legacy then lower count, one-item direct wins, and non-finite pair-vector errors.

- [ ] **Step 2: Run tests and confirm the handwritten model runtime is absent**

Run: `cargo test -p sword2-lib factorized_model --no-fail-fast`

Expected: compile failure.

- [ ] **Step 3: Implement the static types and validator**

```rust
pub(crate) const LEAF_FEATURE: u16 = u16::MAX;

#[derive(Debug, Clone, Copy)]
pub(crate) struct StaticNode {
    pub feature: u16,
    pub threshold: f64,
    pub left: u16,
    pub right: u16,
    pub leaf_value: f64,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct StaticTree { pub root: u16 }

pub(crate) struct StaticBoostedModel {
    pub schema_version: u32,
    pub feature_names: &'static [&'static str],
    pub initial_log_odds: f64,
    pub learning_rate: f64,
    pub trees: &'static [StaticTree],
    pub nodes: &'static [StaticNode],
}

pub(crate) fn validate_model(model: &StaticBoostedModel) -> Result<(), ModelError>;
pub(crate) fn predict_probability(model: &StaticBoostedModel, features: &[f64]) -> Result<f64, ModelError>;
```

Validation walks every root with a three-color DFS, rejects cycles/unreachable invalid children/depth over 3, requires leaf sentinel consistency, finite thresholds/leaves/init/rate, exact schema version and feature count, and at most 96 trees. A separate `validate_embedded_models()` enforces combined 192 trees/2,880 nodes; unique generated feature names must be valid members of the full pair schema and preserve that schema's relative order, allowing ablation-rejected families to be absent.

Traversal takes left on `value <= threshold`, sums one leaf per tree, applies learning rate and stable sigmoid, and rejects any non-finite input/output.

- [ ] **Step 4: Implement symmetrized Borda and full two-head selection**

```rust
pub(crate) struct FactorizedSelection {
    pub measure_index: usize,
    pub num_domains: usize,
    pub canonical: String,
}

pub(crate) fn select_factorized(
    lattice: &CandidateLattice,
    global: &GlobalFeatures,
    counts: &[CountFeatures],
    candidates: &[CandidateFeatures],
    legacy_count: usize,
) -> Result<FactorizedSelection, FactorizedError>;
```

For every unordered pair, build/evaluate both orientations and use `0.5*(P(left,right)+1-P(right,left))`; the reverse win is `1-p`. Iterate opponents in canonical order and average each item's wins. Exactly equal count scores use legacy count when present then lower count; exactly equal candidate scores use canonical text. After choosing count, rank only its candidates. Map the winning candidate back by `source_index`; duplicate/missing mapping is an error.

- [ ] **Step 5: Add the committed cross-language golden integration test**

`factorized_ranker_golden.rs` loads `include_str!("../../benchmark/models/factorized_ranker_v1_golden.json")` using `serde_json` in the test build. For both heads assert generated feature names, raw probabilities, both orientations, symmetrized probability, Borda scores, and final winner within `1e-10`. Assert malformed and NaN cases produce the expected error category without panic.

- [ ] **Step 6: Run full Rust and Python parity verification**

Run:

```bash
cargo fmt -- --check
cargo check
cargo test -p sword2-lib factorized --no-fail-fast
cargo test
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_model_artifact.py \
  benchmark/tests/test_export_factorized_ranker.py -q
```

Expected: Python/Rust deltas are at most `1e-10`; all tests pass.

- [ ] **Step 7: Commit the static inference engine**

```bash
git add sword2-lib/src/sword/factorized_ranker/model.rs \
  sword2-lib/src/sword/factorized_ranker/mod.rs \
  sword2-lib/src/sword/factorized_ranker/generated_model.rs \
  sword2-lib/tests/factorized_ranker_golden.rs
git commit -m "feat: run frozen factorized trees in Rust"
```

---

### Task 15: Pipeline integration, CLI opt-in, and whole-chain fallback

**Files:**
- Modify: `sword2-lib/src/sword/mod.rs`
- Modify: `sword2-lib/src/sword/factorized_ranker/mod.rs`
- Modify: `sword2-cli/src/main.rs`
- Test: inline tests in `sword2-lib/src/sword/mod.rs`
- Test: CLI tests in `sword2-cli/src/main.rs`
- Test: `sword2-lib/tests/factorized_ranker_golden.rs`

**Interfaces:**
- Consumes: `select_legacy`, exact typed feature lattice/context, and static selector.
- Produces: `SwordConfig::use_factorized_ranker`, `--use-factorized-ranker`, selected-count propagation, one-warning fallback, and an end-to-end opt-in binary ready for the locked run.

- [ ] **Step 1: Write default-off, conflict, propagation, and fallback tests**

```rust
#[test]
fn factorized_ranker_defaults_off() {
    assert!(!SwordConfig::default().use_factorized_ranker);
}

#[test]
fn selected_count_and_candidate_propagate_together() {
    let outcome = apply_factorized_fixture(legacy_count_2_factorized_count_3());
    assert_eq!(outcome.n_dom, 3);
    assert_eq!(outcome.to_print, expected_three_domain_line());
}

#[test]
fn any_required_count_group_loss_falls_back_to_exact_legacy_result() {
    let legacy = legacy_fixture_result();
    let outcome = apply_factorized_fixture(malformed_only_candidate_at_one_count());
    assert_eq!(outcome.selection, legacy);
    assert_eq!(outcome.warning_count, 1);
}
```

Add CLI parser tests showing `--use-factorized-ranker` conflicts with `--use-pairwise-reranker`, `--use-count-calibration`, `--count-lambda`, `--use-geometry-metrics`, and `--geometry-lambda`. Add a flag-off integration snapshot comparing raw output lines byte-for-byte to the existing fixture.

- [ ] **Step 2: Run tests and verify the flag/config field do not exist**

Run:

```bash
cargo test -p sword factorized -- --nocapture
cargo test -p sword2-lib factorized_pipeline -- --nocapture
```

Expected: compile/parser failures.

- [ ] **Step 3: Add the opt-in config and clap flag**

Add `pub use_factorized_ranker: bool` to `SwordConfig` and `false` in `Default`. Add:

```rust
/// Select domain count and partition with the embedded factorized structural ranker.
/// Experimental and off by default; falls back to the legacy selector on incomplete evidence.
#[arg(long, conflicts_with_all = [
    "use_pairwise_reranker", "use_count_calibration", "count_lambda",
    "use_geometry_metrics", "geometry_lambda"
])]
use_factorized_ranker: bool,
```

Wire it into the sole `SwordConfig` construction. Do not add `--legacy-selector` yet because the new model is not default.

- [ ] **Step 4: Integrate factorized selection before old experimental branches**

Always compute `LegacySelection` first. If the new flag is off, execute the existing branches exactly. If on:

1. Require complete fresh-run context and validate embedded models.
2. Build the first-pass bounded lattice.
3. Derive `FeatureMask` from `RETAINED_FEATURE_FAMILIES`, then extract all required global/base/domain/boundary/discontinuity/sibling/hierarchy/count features; unused family slots remain zero and are not allowed in the model's feature-name arrays.
4. Apply the malformed-candidate rule: exclude one invalid candidate only if every original available count retains at least one valid candidate; otherwise return an error for whole-chain fallback.
5. Run `select_factorized` and set both `n_dom` and `to_print` from its result.
6. Run second `ParseMeasure` around the factorized count for display alternatives; if the factorized winner is absent, prepend it, then canonical-deduplicate alternatives.

On any error, restore both legacy count and legacy line and emit exactly one event:

```rust
tracing::warn!(
    target: "sword2",
    error = %error,
    "factorized selector unavailable; using legacy selector"
);
```

Do not call the old pairwise/count/geometry selectors after a successful factorized choice.

- [ ] **Step 5: Verify representative fresh and cache-hit behavior**

Use the committed 1JX4 fixture or a benchmark cache chain. Run once into a fresh output and once into the same intermediate directory. The fresh opt-in run must select deterministically; the typed-cache-missing run must complete with the legacy result and one warning, not panic. Also run without the flag twice and compare summaries with `cmp`.

```bash
cargo build --release
tmp_dir=$(mktemp -d)
./target/release/sword2 -i benchmark/cache/chains/1jx4A.pdb \
  -o "$tmp_dir/legacy-a" --threads 1
./target/release/sword2 -i benchmark/cache/chains/1jx4A.pdb \
  -o "$tmp_dir/legacy-b" --threads 1
cmp "$tmp_dir/legacy-a/SWORD2_summary.txt" "$tmp_dir/legacy-b/SWORD2_summary.txt"
./target/release/sword2 -i benchmark/cache/chains/1jx4A.pdb \
  -o "$tmp_dir/factorized" --threads 1 --use-factorized-ranker \
  >"$tmp_dir/factorized-fresh.stdout" 2>"$tmp_dir/factorized-fresh.stderr"
./target/release/sword2 -i benchmark/cache/chains/1jx4A.pdb \
  -o "$tmp_dir/factorized" --threads 1 --use-factorized-ranker \
  >"$tmp_dir/factorized-cache.stdout" 2>"$tmp_dir/factorized-cache.stderr"
test "$(rg -c 'factorized selector unavailable; using legacy selector' \
  "$tmp_dir/factorized-cache.stderr")" -eq 1
```

- [ ] **Step 6: Run repository verification**

Run:

```bash
cargo fmt -- --check
cargo check
cargo test
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: all pass; flag-off snapshot is unchanged.

- [ ] **Step 7: Commit opt-in integration**

```bash
git add sword2-lib/src/sword/mod.rs \
  sword2-lib/src/sword/factorized_ranker/mod.rs \
  sword2-lib/tests/factorized_ranker_golden.rs \
  sword2-cli/src/main.rs
git commit -m "feat: add opt-in factorized structural selector"
```

---

### Task 16: Locked benchmark manifest and acceptance-report tooling

**Files:**
- Create: `benchmark/evaluate_factorized_acceptance.py`
- Modify: `benchmark/run_benchmark.py`
- Modify: `benchmark/compare_sword2_experiments.py`
- Test: `benchmark/tests/test_factorized_acceptance.py`
- Test: `benchmark/tests/test_compare_sword2_experiments.py`
- Create: `benchmark/tests/test_run_benchmark.py`

**Interfaces:**
- Consumes: benchmark `scores.csv`/`runs.csv`, frozen model manifest, Task 9 paired bootstrap, and current competitor outputs only for comparison.
- Produces: gate-by-gate JSON/Markdown results with exact common-chain denominators and resource comparisons. Task 17 invokes it without modifying model code.

- [ ] **Step 1: Write accuracy/resource gate tests**

```python
def test_accuracy_gates_use_declared_common_subsets():
    result = evaluate_accuracy_gates(synthetic_scores(), frozen_baseline_predictions())
    assert result["merizo_ndo_delta"]["n"] == 663
    assert result["chainsaw_ndo_delta"]["n"] == 623
    assert result["count_accuracy"]["threshold"] == 0.745
    assert result["boundary_f1_10"]["threshold"] == 0.620


def test_resource_gates_use_paired_median_runtime_and_max_peak_rss():
    result = evaluate_resource_gates(resource_runs())
    assert result["runtime"]["statistic"] == "paired median per-chain ratio"
    assert result["peak_rss"]["statistic"] == "maximum per-process RSS ratio"


def test_report_fails_closed_on_missing_denominator_or_hash():
    with pytest.raises(ValueError):
        evaluate_acceptance(missing_chain_scores(), mismatched_model_manifest())
```

- [ ] **Step 2: Run tests and confirm acceptance tooling is absent**

Run: `benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_acceptance.py -q`

Expected: import failure.

- [ ] **Step 3: Expand benchmark manifests before collecting locked evidence**

Extend `run_benchmark.py::write_manifest` to record executable SHA-256, factorized top-level/count/candidate/schema/fold/corpus hashes, exact command/extra args, Git commit, OS/kernel/CPU/RAM, Rust version, Python and package versions, thread count, and cache/reuse mode. Refuse to label a run “locked” if a hash differs from `factorized_ranker_v1_manifest.json` or `--skip-existing`/reused SWORD outputs are enabled.

- [ ] **Step 4: Implement the gate evaluator**

Implement `evaluate_accuracy_gates(scores, baseline_predictions, competitors) -> dict[str, object]`, `evaluate_resource_gates(baseline_runs, factorized_runs) -> dict[str, object]`, and `write_acceptance_report(json_path, markdown_path, results) -> None`.

Use 10,000 paired chain bootstraps with seed 37. Merizo must have exactly 663 common successful chains; Chainsaw uses its recorded successful intersection and reports the actual denominator. Compute factorized NDO/count/BF1 from rank-1 optimal rows. Compare contiguous/discontinuous CATH-663 NDO to `benchmark/models/cath663_standalone_structural_baseline.csv` on exact common chain IDs, with thresholds `-0.005`; its hash must match the frozen model manifest. Runtime uses paired per-chain factorized/baseline ratios and gates their median at `1.15`; report p95 diagnostically. RSS gates `max(factorized_peak_rss)/max(baseline_peak_rss) <= 1.10`, with median/p95 diagnostic. Every gate record contains value, threshold, denominator, method, pass boolean, and source hashes.

Use strict comparisons exactly where approved: overall NDO must be `>0.8389`; each competitor CI lower bound must be `>0.0`; count accuracy and BF1 are inclusive `>=0.745` and `>=0.620`; each cohort delta is inclusive `>=-0.005`; runtime/RSS ratios are inclusive `<=1.15` and `<=1.10`.

- [ ] **Step 5: Produce deterministic reports without external writes**

Canonical JSON has no timestamp; Markdown is rendered from JSON in a fixed gate order. Exit code is 0 only when all gates pass, 1 when any measured gate fails, and 2 for invalid/missing evidence. Extend `compare_sword2_experiments.py` to include paired CI columns but preserve existing outputs when no acceptance arguments are given.

- [ ] **Step 6: Run tooling tests and full verification**

Run:

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_acceptance.py \
  benchmark/tests/test_compare_sword2_experiments.py \
  benchmark/tests/test_run_benchmark.py -q
benchmark/.venv/bin/python -m pytest benchmark/tests -q
cargo check
cargo test
```

Expected: all pass.

- [ ] **Step 7: Commit acceptance tooling before the locked run**

```bash
git add benchmark/evaluate_factorized_acceptance.py \
  benchmark/run_benchmark.py benchmark/compare_sword2_experiments.py \
  benchmark/tests/test_factorized_acceptance.py \
  benchmark/tests/test_compare_sword2_experiments.py \
  benchmark/tests/test_run_benchmark.py
git commit -m "feat: certify factorized benchmark acceptance gates"
```

---

### Task 17: Execute the single locked CATH-663 benchmark and publish the engineering decision

**Files:**
- Generate: `benchmark/results_factorized_locked_legacy/benchmark_manifest.json`
- Generate: `benchmark/results_factorized_locked_legacy/scores.csv`
- Generate: `benchmark/results_factorized_locked_legacy/runs.csv`
- Generate: `benchmark/results_factorized_locked_model/benchmark_manifest.json`
- Generate: `benchmark/results_factorized_locked_model/scores.csv`
- Generate: `benchmark/results_factorized_locked_model/runs.csv`
- Generate: `benchmark/models/factorized_ranker_v1_acceptance.json`
- Generate: `benchmark/FACTORIZED_RANKER_ACCEPTANCE.md`
- Modify: `benchmark/REPORT.md`

**Interfaces:**
- Consumes: the committed frozen binary/model and acceptance tooling.
- Produces: the final go/no-go decision. No feature, model, threshold, or hyperparameter change is authorized in this task.

- [ ] **Step 1: Verify the freeze and create a release binary from it**

Run:

```bash
git status --short
benchmark/.venv/bin/python -c \
  'from pathlib import Path; from benchmark.factorized_ranker.model_artifact import verify_top_level_manifest; verify_top_level_manifest(Path("benchmark/models/factorized_ranker_v1_manifest.json"))'
cargo build --release
cargo check
cargo test
benchmark/.venv/bin/python -m pytest benchmark/tests -q
sha256sum target/release/sword2
```

Expected: every frozen hash and test passes. Unrelated dirty files may exist, but no model/runtime source named in the manifest may differ. If it differs, stop; do not run CATH-663.

- [ ] **Step 2: Run fresh baseline and factorized SWORD2 on identical inputs**

Use the existing CATH-663 raw structure cache, the same machine and thread count, and no SWORD output reuse. Run the legacy/competitor benchmark first and the factorized SWORD2 benchmark immediately afterward; record this fixed order because the current harness runs one SWORD2 configuration at a time.

```bash
benchmark/.venv/bin/python -m benchmark.run_benchmark \
  --dataset cath663 \
  --cache-dir benchmark/cache \
  --results-dir benchmark/results_factorized_locked_legacy \
  --tools sword2-rust,merizo,chainsaw \
  --sword2-threads 1 \
  --no-download
benchmark/.venv/bin/python -m benchmark.run_benchmark \
  --dataset cath663 \
  --cache-dir benchmark/cache \
  --results-dir benchmark/results_factorized_locked_model \
  --tools sword2-rust \
  --sword2-threads 1 \
  --sword2-extra-args=--use-factorized-ranker \
  --no-download --strict
```

Do not pass `--skip-existing` or `--reuse-tool-results-dir`. Competitor runs are comparisons only; their outputs never enter SWORD2 training or inference.

- [ ] **Step 3: Validate coverage before looking at scores**

Run a manifest-only validator that checks 663 factorized and baseline successes, exact input IDs, factorized warning/fallback counts, Merizo 663 successes, the recorded Chainsaw successful subset, complete per-chain runtime/RSS, matching executable/model hashes, and no reused SWORD output. Invalid coverage exits 2 and blocks evaluation.

- [ ] **Step 4: Produce the immutable acceptance report**

Run:

```bash
benchmark/.venv/bin/python -m benchmark.evaluate_factorized_acceptance \
  --baseline-scores benchmark/results_factorized_locked_legacy/scores.csv \
  --baseline-runs benchmark/results_factorized_locked_legacy/runs.csv \
  --baseline-manifest benchmark/results_factorized_locked_legacy/benchmark_manifest.json \
  --factorized-scores benchmark/results_factorized_locked_model/scores.csv \
  --factorized-runs benchmark/results_factorized_locked_model/runs.csv \
  --factorized-manifest benchmark/results_factorized_locked_model/benchmark_manifest.json \
  --model-manifest benchmark/models/factorized_ranker_v1_manifest.json \
  --standalone-baseline benchmark/models/cath663_standalone_structural_baseline.csv \
  --json-out benchmark/models/factorized_ranker_v1_acceptance.json \
  --markdown-out benchmark/FACTORIZED_RANKER_ACCEPTANCE.md \
  --bootstrap-replicates 10000 --seed 37
```

Expected: a complete pass/fail record for all nine approved gates, with common-subset IDs/denominators and hashes. Do not rerun after seeing a failure unless the run is invalidated by a documented infrastructure fault; a valid model failure remains the engineering result.

- [ ] **Step 5: Document the result without tuning language**

Append to `benchmark/REPORT.md`: development protocol, retained/rejected families, corpus/fold/model hashes, exact locked command, every metric/CI, count and structural cohorts, fallback count, runtime/RSS statistics, competitor denominators, and the note that publication-quality claims require another untouched benchmark because CATH-663 informed diagnosis/design.

- [ ] **Step 6: Commit the locked evidence and decision**

Stage only small manifests/reports and the locked scores/runs necessary to audit the result; do not stage raw tool output directories.

```bash
git add benchmark/results_factorized_locked_legacy/benchmark_manifest.json \
  benchmark/results_factorized_locked_legacy/scores.csv \
  benchmark/results_factorized_locked_legacy/runs.csv \
  benchmark/results_factorized_locked_model/benchmark_manifest.json \
  benchmark/results_factorized_locked_model/scores.csv \
  benchmark/results_factorized_locked_model/runs.csv \
  benchmark/models/factorized_ranker_v1_acceptance.json \
  benchmark/FACTORIZED_RANKER_ACCEPTANCE.md benchmark/REPORT.md
git commit -m "bench: record locked factorized ranker evaluation"
```

---

### Task 18: Conditional default-on rollout or documented opt-in retention

**Files:**
- Modify only if all gates pass: `sword2-cli/src/main.rs`
- Modify only if all gates pass: `sword2-lib/src/sword/mod.rs`
- Modify: `README.md`
- Modify: `benchmark/REPORT.md`
- Test: CLI/unit tests beside modified Rust code

**Interfaces:**
- Consumes: `factorized_ranker_v1_acceptance.json` from Task 17.
- Produces: either a default-on selector with explicit `--legacy-selector`, or an accurately documented opt-in experiment with no runtime default change.

- [ ] **Step 1: Gate the branch mechanically**

Run:

```bash
benchmark/.venv/bin/python -c \
  'import json; r=json.load(open("benchmark/models/factorized_ranker_v1_acceptance.json")); print("PASS" if r["all_gates_pass"] else "KEEP_OPT_IN")'
```

Expected: exactly one of `PASS` or `KEEP_OPT_IN`.

- [ ] **Step 2A: If any gate failed, retain opt-in behavior**

Do not edit selector defaults or add `--legacy-selector`. Update README/report with the opt-in command, failed gates, cohort diagnosis, hashes, fallback behavior, and rollback/removal instructions. Run documentation link checks if present, then commit:

```bash
git add README.md benchmark/REPORT.md
git commit -m "docs: retain factorized ranker as opt-in"
```

This completes the implementation even though promotion failed; the valid locked result must not trigger a new feature attempt in this plan.

- [ ] **Step 2B: If and only if every gate passed, write failing default/escape-hatch tests**

```rust
#[test]
fn factorized_is_default_after_accepted_rollout() {
    assert!(SwordConfig::default().use_factorized_ranker);
}

#[test]
fn legacy_selector_disables_factorized_path() {
    let cli = Cli::try_parse_from(["sword2", "--legacy-selector", "-i", "x.pdb"]).unwrap();
    assert!(cli.legacy_selector);
}
```

Run the tests and confirm they fail before changing defaults.

- [ ] **Step 3B: Implement the accepted default and compatibility switches**

Set `SwordConfig::default().use_factorized_ranker=true`. Add `--legacy-selector`, conflicting with all experimental selector flags, and construct the CLI config with `use_factorized_ranker: !cli.legacy_selector`; do not continue wiring it from the opt-in boolean. Keep `--use-factorized-ranker` as a hidden/deprecated compatibility alias that is accepted but redundant and conflicts with `--legacy-selector`. Fallback always remains the exact legacy path. Update CLI help, README, and report with model hash, training protocol, locked metrics, runtime cost, fallback conditions, and escape hatch.

- [ ] **Step 4B: Verify default output, explicit legacy output, and all tests**

Run:

```bash
cargo fmt -- --check
cargo check
cargo test
benchmark/.venv/bin/python -m pytest benchmark/tests -q
./target/release/sword2 --help | rg 'legacy-selector'
```

Expected: all pass; default selects factorized, `--legacy-selector` reproduces the pre-rollout snapshot byte-for-byte.

- [ ] **Step 5B: Commit the accepted rollout**

```bash
git add sword2-cli/src/main.rs sword2-lib/src/sword/mod.rs \
  README.md benchmark/REPORT.md
git commit -m "feat: make accepted factorized selector the default"
```

---

## Final Verification Checklist

- [ ] `benchmark/.venv/bin/python -m pytest benchmark/tests -q` passes.
- [ ] `cargo fmt -- --check`, `cargo check`, and `cargo test` pass.
- [ ] Strict Python/Rust partition parsing and synthetic feature vectors agree.
- [ ] Frozen Python and Rust head probabilities/Borda winners agree within `1e-10`.
- [ ] Generated model caps are at most 192 trees and 2,880 nodes combined.
- [ ] The corpus, fold, schema, model, binary, golden, and locked-run hashes form one verified chain in the manifests.
- [ ] No feature/artifact name contains or imports competitor predictions.
- [ ] Default behavior changes only on the mechanically verified all-gates-pass branch; otherwise the factorized selector remains opt-in.
