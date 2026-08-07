# Factorized Compact Structural Ranker Design

**Date:** 2026-08-07
**Status:** Approved for implementation planning
**Scope:** Standalone SWORD2 candidate selection; no external predictor at training or inference time

## Context

SWORD2 already generates a strong candidate lattice on CATH-663. The best
candidate available per chain reaches NDO `0.9019`, while the legacy rank-1
selection reaches `0.7772`. A pointwise ExtraTrees model trained on 352,586
candidates from 17,283 CATH-17287 chains improved rank-1 NDO to `0.8142`, but
remained below Merizo's `0.8389`.

The apparent CATH-17287 validation score of `0.8665` was optimistic. A random
chain split allowed 91.6% of validation chains to reuse an exact CATH family
combination present in the fit set; only 24.0% of CATH-663 chains have a family
combination present in CATH-17287. Holding out exact CATH combinations lowers
the same 100-tree model from `0.8665` to `0.8362`, while its candidate oracle
remains near the CATH-663 oracle.

The CATH-663 selector loss decomposes into two actionable parts:

- `0.0538` mean NDO is lost by selecting the wrong domain count.
- `0.0338` mean NDO is lost while selecting among candidates at the already
  chosen count.

The best candidate at the model's selected count reaches NDO `0.8480`, already
above Merizo. Candidate generation is therefore not the immediate bottleneck.
The new selector must improve count selection and within-count ranking without
using Merizo, Chainsaw, or any other predictor as a feature or teacher.

## Goals

1. Factor domain-count selection from candidate selection within that count.
2. Add global chain context, boundary-local evidence, domain-conditioned
   geometry, and candidate-hierarchy context missing from the current ranker.
3. Train compact models with validation that holds out structural-family
   combinations rather than only individual chain identifiers.
4. Embed frozen model weights in the Rust binary with deterministic,
   dependency-free inference.
5. Beat Merizo and Chainsaw independently on the locked CATH-663 benchmark.

## Non-goals

- Do not use competitor predictions, embeddings, scores, or pseudo-labels.
- Do not replace SWORD2 candidate generation or Peeling in this project.
- Do not add a Python runtime, model server, GPU dependency, or native ML
  library to the released binary.
- Do not optimize a conventional weakest-contact cut as the primary rule.
  Correct CATH boundaries in the failure cohort often retain more interface
  contact than the wrong SWORD2 selection.
- Do not enable the model by default until every acceptance gate passes.

## Architecture

The selector has two independent pairwise boosted-tree heads. Both operate on
the candidates already produced by SWORD2.

### 1. Count head

The count head compares every domain count represented in the candidate
lattice. It receives global chain context plus summaries of the candidates
available at each count. It is trained to order counts by absolute distance
from the CATH reference count, whose distance is zero.

For each pair of available counts `(k_i, k_j)`, the model predicts which count
is more likely to be correct. Pairwise probabilities are symmetrized and
aggregated with normalized Borda scoring. The highest-scoring count wins.
This representation supports arbitrary counts present in the lattice and does
not require a fixed multiclass output dimension.

### 2. Within-count candidate head

The candidate head compares candidates only with siblings having the selected
count. It learns the sign of the NDO difference for each non-tied pair. The
same symmetrized Borda rule selects the winning partition.

Keeping the two heads separate prevents a corpus-wide domain-count bias from
overwhelming the finer boundary decision. It also makes count and within-count
regret independently measurable.

### Runtime data flow

1. Generate the normal SWORD2 candidate lattice.
2. Compute global chain features once.
3. Compute candidate and hierarchy features once per candidate.
4. Group candidates by their number of domains.
5. Apply the count head to the available count groups.
6. Apply the candidate head to candidates in the winning group.
7. Return the Borda winner through the existing output pipeline.
8. Fall back to the legacy distance-model selector if validation of the model
   artifact, chain coverage, or any required feature fails.

No absolute score is mixed across the two heads.
If only one count or one candidate is available at its respective stage, that
entry wins directly without invoking the corresponding pairwise model.

## Feature contract

All runtime features come from the input structure, the pure-Rust DSSP and
Peeling results, and SWORD2's own candidate lattice.

### Global chain features

The count head receives:

- Number of residues.
- Radius of gyration normalized by `n_residues^(1/3)`.
- Ratios of the three principal inertia axes.
- Nonlocal C-alpha contact density using the existing logistic contact model.
- Contact order: contact-weighted sequence separation divided by chain length.
- DSSP helix, strand, and coil fractions.
- Number of contiguous helix and strand blocks.
- Number of Peeling levels and finest-level Protein Units.
- Histogram of candidate counts represented in the lattice.
- Modal candidate count and distance from each proposed count to the mode.

### Per-count summary features

For each available count, the count head receives:

- Proposed count and number of candidates at that count.
- Minimum, mean, and maximum legacy distance-model score.
- Minimum, mean, and maximum `max_cr`, `density_min`, `mean_density`, minimum
  domain size, contact-Q summaries, fragmentation, and boundary integrity.
- Difference between those summaries and the corresponding summaries at the
  neighboring available counts.

### Existing per-candidate features

The candidate head retains the tested energy-free features:

- `num_domains`, `min_size`, `max_cr`, `density_min`, `mean_density`.
- Boundary coil fraction and modal-count distance.
- Domain `q1/q2/q3`, volume-ratio, and density summaries.
- Inter-domain contact-Q summaries.
- Segment count, discontinuity count, size balance, largest-domain fraction,
  and segment-size summaries.

### Domain-conditioned features

Aggregate means currently hide the asymmetric-domain signature seen in the
largest failures. The candidate head adds:

- Size fraction, `q1/q2/q3`, relative density, internal contact density, and
  contact order for the smallest and largest predicted domains.
- Smallest-to-largest ratios for density, internal contact density, and each
  gyration axis.
- Minimum, mean, and maximum per-domain internal-contact fraction and graph
  conductance. These are learned descriptors, not monotonic cut penalties.

### Boundary-local features

For every sequential boundary introduced by a candidate, compute:

- Whether it falls inside helix, strand, or coil and its distance to the
  nearest secondary-structure-element terminus.
- Counts and summed energies of DSSP hydrogen bonds crossing the boundary.
- Counts of beta-bridge and sheet relationships crossing the boundary.
- Logistic-contact insulation in residue windows of 8, 16, and 32 residues.
- Long-range contacts crossing the boundary, using sequence separation >= 8.
- Local C-alpha bend and virtual-dihedral changes as hinge descriptors.

Candidate features contain the minimum, mean, and maximum over boundaries.
The model may learn that a valid boundary can retain substantial interface
contact; no boundary-contact feature is assigned a fixed favorable direction.

### Candidate-relative and hierarchy features

For each candidate among siblings of the same count, compute:

- Percentile rank and delta to the sibling median for every core geometry,
  density, contact, fragmentation, and boundary-integrity feature.
- Peeling level of first appearance.
- Number of levels for which the exact partition remains available.
- Parent and child merge-quality margins.
- Difference in CR and density score from the nearest competing sibling.
- Number of distinct hierarchy paths producing the partition.

Persistence is contextual only. It is not rewarded monotonically because the
wrong boundary was more persistent than the oracle in the diagnosed
correct-count failures.

### Discontinuous segment-association features

For candidates containing discontinuous domains, compute:

- Contact affinity of each segment to the other segments assigned to its
  domain.
- Affinity margin over the best competing domain.
- Long-range internal-contact capture and interface-span entropy.
- DSSP beta-sheet links joining segments assigned to the same domain.

Continuous candidates receive an explicit `has_discontinuity = 0` indicator
and zero values for the conditional summaries.

## Pairwise training

Both heads use `sklearn.ensemble.GradientBoostingClassifier` offline. Each
head is limited to at most 96 trees of maximum depth 3. The bounded model grid
is:

- Trees: `{64, 96}`.
- Learning rate: `{0.03, 0.05}`.
- Minimum samples per leaf: `{32, 64}`.

No other model family or hyperparameter is selected using CATH-663.

### Pair construction

- Count head: for each chain whose true count exists in the lattice, pair all
  available counts whose absolute distances from the true count differ. The
  count closer to truth wins; equal-distance pairs are skipped. Add both
  orientations for every retained pair.
- Candidate head: within every chain and count, pair candidates whose NDO
  values differ and add both orientations.
- Exact NDO ties are skipped.
- Every chain has total training weight 1.0 in each head, divided equally over
  its retained ordered pairs.
- When a chain yields more than 64 unordered pairs for one head, retain exactly
  64 using deterministic seed 37, stratified over target-difference quartiles
  before mirroring.

Pair inputs contain shared chain context, signed feature differences, and
absolute feature differences. Mirrored training makes the comparison
direction explicit. At inference, the probability that `i` beats `j` is:

`0.5 * (P(i, j) + 1 - P(j, i))`.

The normalized Borda score is the mean win probability against all siblings.
A count-head tie prefers the count chosen by the legacy selector when present,
then the lower numeric count. A candidate-head tie is broken by canonical
delineation text. These rules make output independent of enumeration order.

## Training-data integrity

Residue-coordinate mapping must fail closed. The training-table builder must
not fall back from an unavailable canonical mapping to raw author-numbered
CATH chopping. A chain is excluded and reported if:

- Author numbering cannot be mapped to canonical zero-based indices.
- Mapped residue count disagrees with the cached chain.
- Mapped true-domain count disagrees with CATH metadata.
- Any candidate references an out-of-range residue or overlaps another domain.

The six known malformed CATH-17287 chains and their 26 corrupted candidate
rows are excluded by these rules. A machine-readable rejection report is part
of every corpus build.

## Validation protocol

### Development folds

CATH metadata is used only to form validation groups and labels, never as a
runtime feature.

1. Construct connected groups that keep the same four-character PDB ID or the
   same exact sorted CATH family combination in one fold. The combination is
   the multiset of `topology_family` suffixes after `:` in the reference
   chopping, excluding the sentinel label `999_999`.
2. Assign groups deterministically to five folds with greedy balancing over:
   chain count, true-domain count bins `{2, 3, 4, 5+}`, and chain-length bins
   `<250`, `250-349`, `350-449`, and `450+`.
3. Verify that no PDB ID or exact family combination crosses folds.
4. Report macro means over chains, never over candidate rows.

Individual CATH family labels cannot be made fully disjoint without joining
nearly the complete multi-domain corpus into one connected component. Exact
family combinations plus PDB grouping are therefore the enforceable primary
split. Results are additionally reported for validation chains containing
seen versus unseen individual family labels.

### Feature-family ablations

Start from the existing 22-feature model and add feature families in this
order:

1. Global count context.
2. Domain-conditioned features.
3. Boundary-local features.
4. Candidate-relative and hierarchy features.
5. Discontinuous segment association.

A family is retained only when its five-fold mean NDO improves and its paired
chain-bootstrap 95% confidence interval is non-negative. Count-head features
must also avoid reducing domain-count accuracy. The discontinuity family must
not lower contiguous-chain NDO by more than `0.005`.

### Locked benchmark

The existing CATH-663 results informed the problem, this approved feature
design, and the acceptance thresholds. After this specification is approved,
no additional model choice, feature choice, or hyperparameter choice is made
from a new CATH-663 result. After grouped CATH-17287 development is complete:

1. Freeze the feature schema, model artifacts, and their SHA-256 hashes.
2. Commit the frozen artifacts before benchmark execution.
3. Run SWORD2, Merizo, and Chainsaw once on CATH-663.
4. Produce paired chain-bootstrap confidence intervals with 10,000 replicates
   and deterministic seed 37.

Because CATH-663 was inspected during diagnosis, publication-quality claims
will require an additional untouched benchmark. It remains the engineering
gate requested for this project.

## Model artifact and Rust inference

Offline training writes a versioned JSON artifact containing:

- Schema version and head name.
- Ordered feature names.
- Training-corpus and fold-manifest hashes.
- Tree nodes: feature index, threshold, left child, right child, and leaf
  value.
- Learning rate and initial log odds. No post-hoc probability calibration is
  applied.
- Training command, seed, and retained feature-family list.

A deterministic exporter converts the JSON into static Rust arrays committed
under `sword2-lib/src/sword/`. The release binary does not parse a model file
at runtime.

Rust validates the compiled schema version, ordered feature count, tree child
indices, finite thresholds, and finite leaf values in tests. Inference uses
`f64` throughout. Python and Rust probabilities must agree within `1e-10` on
committed golden vectors.

The combined two-head model may contain no more than 192 depth-3 trees and
2,880 tree nodes.

## Error handling and fallback

The factorized selector is skipped for the entire chain, rather than partially
applied, when any of these conditions occurs:

- Candidate residue coverage is invalid or inconsistent.
- Required DSSP, Peeling, coordinate, or hierarchy data is unavailable.
- Any required feature is non-finite.
- The candidate lattice has no valid count group.
- Compiled model validation fails.

Fallback uses the current legacy distance-model path and emits one structured
warning. A single malformed candidate is excluded only when at least one other
complete candidate remains at every available count; otherwise the whole
factorized selection falls back. No panic is permitted for model or feature
errors.

## Testing

### Python tests

- Canonical numbering fails closed and produces rejection records.
- Every new feature has synthetic geometric and DSSP fixtures with exact
  expected values.
- Pair construction mirrors orientations, skips ties, and gives every chain
  total weight 1.0.
- Fold manifests have zero PDB-ID and exact-family-combination leakage.
- Borda selection is invariant under candidate row permutation.
- Exported tree predictions reproduce scikit-learn probabilities.
- Corpus rebuilds and model training are deterministic for seed 37.

### Rust tests

- Chain, domain, boundary, hierarchy, and discontinuity feature extraction.
- Tree traversal at thresholds and leaves.
- Symmetrized pair probability and normalized Borda aggregation.
- Candidate-order permutation invariance and canonical tie-breaking.
- Python/Rust golden-vector parity within `1e-10`.
- Model-schema and non-finite-feature fallback.
- Representative contiguous, discontinuous, two-domain, and four-domain
  integration fixtures.

### Repository verification

Every Rust checkpoint runs:

```bash
cargo check
cargo test
```

Every benchmark-tooling checkpoint runs the focused tests plus:

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

## Acceptance gates

The factorized selector is eligible to become the default only when all gates
pass on the locked CATH-663 run:

1. Mean NDO is greater than Merizo's established `0.8389`.
2. The paired-bootstrap 95% confidence interval for NDO delta versus Merizo is
   strictly positive on their 663 common chains.
3. NDO also exceeds Chainsaw with a strictly positive paired interval on the
   common successful-chain subset.
4. Domain-count accuracy is at least `0.745`.
5. Boundary F1@10 is at least `0.620`.
6. Contiguous-chain and discontinuous-chain NDO each regress by no more than
   `0.005` versus the best standalone SWORD2 development baseline.
7. Median end-to-end runtime is no more than 15% above current SWORD2 on the
   same hardware and thread count.
8. Peak RSS is no more than 10% above current SWORD2.
9. All Python and Rust tests pass and Python/Rust model parity holds.

If any gate fails, the new selector remains opt-in and the failure is reported
by count and structural cohort before another feature family is attempted.

## Rollout

1. Implement behind `--use-factorized-ranker`; legacy selection remains the
   default and fallback during development.
2. Commit each independently verified tooling, feature, training, export, Rust
   inference, and benchmark checkpoint. Keep generated candidate/training
   corpora ignored; commit only their reproducibility manifests and hashes.
3. Freeze and commit the model artifacts before the locked benchmark.
4. If every acceptance gate passes, make the factorized selector the default
   and retain `--legacy-selector` as an explicit compatibility escape hatch.
5. Document the training corpus, model hash, validation protocol, feature
   schema, benchmark result, and runtime cost in the repository report.
