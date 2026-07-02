# Reranking Metrics Design

## Goal

Substantially raise the quality of SWORD2's rank-1 domain partition pick, closing the gap between what the pipeline already generates and what it selects.

## Evidence

`benchmark/results_*/figures/sword2_topk_ndo.png` shows the rank-1 ("optimal") pick averages **NDO ≈ 0.77–0.78**, while the best NDO already present among the top ~10 alternatives the pipeline generates for the same entries averages **≈ 0.89–0.90**. `sword2_oracle_gap.png` confirms a long right tail of entries where a discarded alternative scores far above the pick across NDO, boundary F1, matched-domain Dice, and pairwise F1. The candidates are already good; selection is the bottleneck.

## Current Selection Logic

`sword2-lib/src/sword/distance_model.rs` ranks candidates using a single hand-fit linear decision boundary over exactly two features: `max_cr` (contact ratio) and `mean_density`, both purely geometric properties of the PU contact matrix, ported unchanged from 2003-era SWORD1 Perl (`DistanceModel.pm`). It is used both to pick `n_dom` (predicted domain count) and to pick the representative candidate within that count, in `sword2-lib/src/sword/mod.rs`.

The Phase B logistic re-ranker (commits `c6a691f`..`08909be`, reverted) added 7 more features — `num_domains`, `min_size`, `density_min`, `n_discontinuous`, `size_balance`, `largest_domain_frac`, `mean_junction_support` — all variations on PU-merge geometry/topology. It regressed CATH-663 NDO from 0.777 to 0.719 because a pointwise logistic model learned a strong negative weight on `num_domains`, systematically under-segmenting.

Two signals the pipeline already computes are **never used for selection**:

- **Pseudo-energy Z-score** (`sword2-lib/src/energy/score.rs`): a statistical-potential-based energy score with a decoy-shuffle Z-score, computed only for the already-chosen winner, gated behind the optional `-E` display flag.
- **DSSP secondary structure** (`sword2-lib/src/dssp/`): computed for every residue before selection runs, but not consulted by it.

## New Metric 1: Energy Z-score shortlist rescoring

For the candidate set the pipeline already produces per protein — the `alt_b`/`alt_l`-bounded alternatives from the two-pass `parse_measure` (`relevant_measure`/`relevant_measure2` in `sword2-lib/src/sword/mod.rs`), typically 10–20 candidates, matching the x-axis of `sword2_topk_ndo.png` — compute a pseudo-energy Z-score for each candidate's full domain assignment via the existing `energy::get_energy_and_z_score`, at a reduced **200 shuffles** (vs. the 2000 used for display).

This reuses fully-implemented, tested infrastructure. No new candidate generation. The Z-score is information-theoretically distinct from `max_cr`/`mean_density`: those measure contact-geometry compactness of the PU merge, the Z-score measures whether the resulting fold looks real according to a potential trained on known domains. It is the one signal source Phase B never touched.

Cost is bounded to ~10–20 rescoring calls per protein, not the full peeling-level candidate tree. Because this makes energy computation a mandatory pipeline step rather than the current `-E`-gated opt-in, its overhead must be measured against `benchmark/results_*/figures/runtime_s_vs_length.png` before shipping default-on.

## New Metric 2: Boundary secondary-structure (coil) fraction

For each candidate, compute the fraction of its domain-boundary junctions (±2 residues) that fall in Coil/Turn vs. Helix/Strand, using the DSSP secondary structure already computed in Step 1 of `run_pipeline`. Domain linkers are structurally coil far more often than not; cutting a boundary through a helix is a strong tell of a bad partition.

This is effectively free: a lookup against data already in memory/on disk, no new computation.

## Combination: pairwise-ranking-trained reranker

Train a small linear model on **feature differences between candidate pairs within the same protein**:

```
Δenergy_z, Δcoil_fraction, Δmax_cr, Δmean_density,
Δ|num_domains − chain's modal candidate count|, ...
```

predicting `sign(NDO_i − NDO_j)` for pairs `(i, j)` from the same chain — a pairwise (RankNet-style) loss, not Phase B's pointwise logistic. This directly targets the recorded Phase B failure: a pointwise model compared candidates across the whole corpus and learned a global domain-count bias; a pairwise loss only ever compares candidates from the *same* chain, which cannot structurally collapse to that bias. The `|num_domains − modal count|` feature anchors each candidate to what's typical for that specific chain's own candidate spread (memory's suggested fix), not a corpus-wide average.

Training data reuses the existing `SWORD2_DUMP_CANDIDATES` dump mechanism (built for Phase B), extended with the two new per-candidate features, over the CATH benchmark set.

## Evaluation Gate

Must clear the current baseline (CATH-663: NDO 0.777, d_count_acc 0.670 at rank-1) by a real margin, measured with the existing `benchmark/run_benchmark.py`/`figures.py` harness. Specifically, rank-1 mean NDO in `sword2_topk_ndo.png` must move meaningfully toward the k≈10 ceiling (~0.89), not stay pinned near 0.77. Cross-validate by CATH superfamily, not a random split, to avoid homolog leakage inflating the result — Phase B's training/eval split discipline should be re-checked for this before trusting new numbers.

## Risks

- **Runtime**: always-on energy computation changes baseline runtime for every user; needs explicit before/after measurement, not an assumption.
- **Z-score sampling noise**: 200 shuffles has more variance than the 2000-shuffle display version; check whether that noise ever flips a ranking decision between two close candidates before trusting it in production.
- **Overfitting to the benchmark set**: same risk class as Phase B; the superfamily-level split is the mitigation, not optional.

## Rollout

Same discipline as Phase B: build behind a flag, benchmark before flipping default, keep a clean revert path if it regresses. The project already has working infrastructure for this (`SWORD2_DUMP_CANDIDATES`, `run_benchmark.py`, the `results_bypass`/`results_rust_original` comparison pattern) — reuse it rather than rebuilding it.

## Out of Scope

- Retrofitting energy Z-score into the full peeling-level candidate tree (only the already-shortlisted `alt_b`/`alt_l` candidates are rescored).
- Any change to candidate *generation* (Peeling, ComputeMeasure) — this design only changes how already-generated candidates are scored and picked.
- Making `-E`/display-time energy computation reuse the 200-shuffle reranking result — display keeps its own 2000-shuffle run for now; unifying them is a possible future simplification, not part of this design.
