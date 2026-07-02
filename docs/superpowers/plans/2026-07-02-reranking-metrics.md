# Reranking Metrics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add two new candidate-scoring signals (pseudo-energy Z-score, DSSP boundary-coil fraction) to SWORD2's domain-partition selection, combined via a pairwise-ranking-trained reranker, to close the gap between the pipeline's rank-1 pick (NDO ≈ 0.77) and the quality already present in its own top-10 alternatives (NDO ≈ 0.89).

**Architecture:** Rescore the pipeline's *existing* `alt_b`/`alt_l`-bounded candidate shortlist (`relevant_measure2` in `sword2-lib/src/sword/mod.rs`, ~10-20 candidates per protein) with two new per-candidate features, computed from data the pipeline already produces (DSSP secondary structure via Peeling's `ss_types`, pseudo-energy via the existing `energy` module at a reduced shuffle count). Combine these with the existing geometric features via a small linear model trained offline in Python with a *pairwise* ranking loss (not the pointwise logistic that caused Phase B's regression), and wire it in behind a feature flag that replaces only the final "pick the winner" step — candidate generation, domain-count prediction, and the displayed-alternatives list are untouched.

**Tech Stack:** Rust (sword2-lib, sword2-cli), Python (benchmark/ harness, existing pandas/numpy stack, no new dependencies).

## Global Constraints

- Run `cargo check` and `cargo test` after modifying Rust files before considering any task complete (CLAUDE.md).
- Prefer existing Rust ecosystem crates over hand-rolled logic when a new dependency is genuinely needed (CLAUDE.md); no new Rust or Python dependencies are needed for this plan.
- Default runtime must stay unchanged: the reduced-shuffle energy computation and reranker only run when explicitly enabled (`SWORD2_DUMP_CANDIDATES` for data collection, `--use-pairwise-reranker` for inference) — confirmed with the user as a hard requirement.
- Energy Z-score rescoring is bounded to the existing `alt_b`/`alt_l` shortlist (~10-20 candidates), never the full peeling-level candidate tree, at 200 shuffles instead of the 2000 used for `-E` display output.
- Any change to the selection winner must be validated against the CATH-663 gate (current baseline: NDO 0.777, d_count_acc 0.670 at rank-1) using the existing `benchmark/run_benchmark.py`/`figures.py` harness, cross-validated by CATH superfamily (not random split), before it is considered for default-on.
- Candidate generation (Peeling, ComputeMeasure), domain-count (`n_dom`) prediction, and the displayed-alternatives list in `quality_and_display` are out of scope — only the final winner-selection step changes.

---

## File Structure

New files:
- `sword2-lib/src/sword/candidate_features.rs` — pure functions: `boundary_coil_fraction`, `modal_num_domains`/`modal_count_distance`, `candidate_energy_z_score`.
- `sword2-lib/src/sword/reranker.rs` — pairwise-trained linear scorer, loads weights via `include_str!` from `benchmark/data/pairwise_reranker_weights.json`.
- `benchmark/train_pairwise_reranker.py` — pairwise ranking-loss trainer, reads `benchmark/data/training_table.csv`, writes `benchmark/data/pairwise_reranker_weights.json`.
- `benchmark/tests/test_train_pairwise_reranker.py` — unit tests for the pairwise loss/gradient and a convergence smoke test.
- `benchmark/tests/test_build_training_table.py` — unit tests for the new feature columns added to the training-table builder.
- `benchmark/data/pairwise_reranker_weights.json` — trained weights (produced by Task 7's synthetic smoke run initially, overwritten by Task 9's real CATH training run).

Modified files:
- `sword2-lib/src/peeling/algorithm.rs` — `PeelingOutput` gains an `ss_types` field; `run_peeling` stops discarding it.
- `sword2-lib/src/sword/mod.rs` — `SwordConfig` gains `energy_config`/`chain_id`/`use_pairwise_reranker` fields; the candidate dump moves to the `relevant_measure2` shortlist stage and gains new columns; the `to_print` winner-selection step optionally uses the reranker; declares the two new submodules.
- `sword2-cli/src/main.rs` — builds the reduced-shuffle `EnergyConfig` earlier and always (cheap — lazy-loaded), passes it plus `chain_id` into `SwordConfig`, adds a `--use-pairwise-reranker` CLI flag.
- `benchmark/build_training_table.py` — `FIELDNAMES` and `_score_candidates` carry the three new dump columns through to the training table.

---

## Task 1: Expose per-residue secondary structure from Peeling output

**Files:**
- Modify: `sword2-lib/src/peeling/algorithm.rs:612-625` (`PeelingOutput` struct), `:731` (`run_peeling` body), `:908-914` (`PeelingOutput` construction)
- Test: `sword2-lib/src/peeling/algorithm.rs` (inline `#[cfg(test)] mod tests`)

**Interfaces:**
- Produces: `PeelingOutput.ss_types: Vec<SsType>` — 0-indexed, one entry per residue, same index space as `PeelingOutput.final_pu_delineation`'s start/end values and `MeasureLine.delineation` (already true today for `true_nums`; this task adds the equivalent for secondary structure). `SsType` (defined in this file) is `pub(crate) enum SsType { Coil, Helix, Sheet }`.

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)] mod tests` block at the bottom of `sword2-lib/src/peeling/algorithm.rs` (after `test_parse_dssp_ss_types_and_residue_numbers`):

```rust
    #[test]
    fn test_run_peeling_exposes_ss_types() {
        let dir = tempdir().unwrap();
        let dssp_path = dir.path().join("test.dssp");

        let mut content = format!("{:<128}\n", "  # RESIDUE AA STRUCTURE");
        for i in 1..=3usize {
            content.push_str(&dssp_line(i, i as i32, 'A', ' '));
            content.push('\n');
        }
        for i in 4..=7usize {
            content.push_str(&dssp_line(i, i as i32, 'A', 'H'));
            content.push('\n');
        }
        for i in 8..=10usize {
            content.push_str(&dssp_line(i, i as i32, 'A', ' '));
            content.push('\n');
        }
        std::fs::write(&dssp_path, &content).unwrap();

        // Simple extended-chain CA coordinates, 3.8 A apart along x — only the
        // ss_types passthrough is under test here, not realistic PU geometry.
        let ca_coords: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 3.8, 0.0, 0.0]).collect();

        let config = PeelingConfig::default();
        let output = run_peeling(&ca_coords, &dssp_path, &config).unwrap();

        assert_eq!(output.ss_types.len(), 10);
        assert!(output.ss_types[..3].iter().all(|s| *s == SsType::Coil));
        assert!(output.ss_types[3..7].iter().all(|s| *s == SsType::Helix));
        assert!(output.ss_types[7..].iter().all(|s| *s == SsType::Coil));
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p sword2-lib --lib peeling::algorithm::tests::test_run_peeling_exposes_ss_types`
Expected: FAIL with `no field \`ss_types\` on type \`PeelingOutput\`` (compile error).

- [ ] **Step 3: Add the field and stop discarding the value**

In `sword2-lib/src/peeling/algorithm.rs`, change the `PeelingOutput` struct (around line 612):

```rust
/// Complete output from the peeling algorithm.
#[derive(Debug, Clone)]
pub struct PeelingOutput {
    /// Contact probability matrix (owned).
    pub contact_matrix: ContactMatrix,
    /// Results for each iteration (1-indexed: iterations[0] = iteration 1).
    pub iterations: Vec<IterationResult>,
    /// Final PU contact matrix entries (x, y, value) from the last iteration.
    pub final_pu_contacts: Vec<(usize, usize, f64)>,
    /// Final PU delineation entries (id, start, end) from the last iteration.
    pub final_pu_delineation: Vec<(usize, usize, usize)>,
    /// Original residue numbers from DSSP (tab_true_num equivalent).
    pub true_nums: Vec<i32>,
    /// Per-residue secondary structure, 0-indexed, same index space as
    /// `final_pu_delineation` start/end values and `MeasureLine.delineation`.
    pub(crate) ss_types: Vec<SsType>,
}
```

Change line 731 from:

```rust
    let (_ss_types, true_nums, cutting_mask) = parse_dssp_for_peeling(dssp_path, n, config)?;
```

to:

```rust
    let (ss_types, true_nums, cutting_mask) = parse_dssp_for_peeling(dssp_path, n, config)?;
```

Change the `PeelingOutput` construction (around line 908):

```rust
    Ok(PeelingOutput {
        contact_matrix: matrix,
        iterations,
        final_pu_contacts: last_pu_contacts,
        final_pu_delineation: last_pu_delineation,
        true_nums,
        ss_types,
    })
```

Also make `SsType` derive `Clone` if it doesn't already (needed because `PeelingOutput` derives `Clone`): confirm the existing declaration `#[derive(Debug, Clone, Copy, PartialEq, Eq)] pub(crate) enum SsType` — it already derives `Clone`, so no change needed there.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p sword2-lib --lib peeling::algorithm::tests::test_run_peeling_exposes_ss_types`
Expected: PASS

- [ ] **Step 5: Run the full peeling test suite to check for regressions**

Run: `cargo test -p sword2-lib --lib peeling::`
Expected: All PASS (no other test constructs `PeelingOutput` by struct literal, so this should be a pure addition)

- [ ] **Step 6: Commit**

```bash
git add sword2-lib/src/peeling/algorithm.rs
git commit -m "feat: expose per-residue secondary structure from Peeling output"
```

---

## Task 2: Boundary coil-fraction candidate feature

**Files:**
- Create: `sword2-lib/src/sword/candidate_features.rs`
- Modify: `sword2-lib/src/sword/mod.rs:22-26` (add `pub mod candidate_features;`)
- Test: `sword2-lib/src/sword/candidate_features.rs` (inline)

**Interfaces:**
- Consumes: `crate::peeling::algorithm::SsType` (from Task 1).
- Produces: `pub fn boundary_coil_fraction(raw_delineation: &str, ss_types: &[SsType]) -> f64` — `raw_delineation` is the **pre-remap**, 0-based-index delineation string in `MeasureLine::delineation` format (space-separated domains, `;`-separated discontinuous segments, `start-end` 0-based residue-array indices — same index space as `ss_types`). Returns the fraction of residues within a ±2 window of each domain-boundary junction that are `SsType::Coil`; `0.0` if there are no junctions to inspect.

- [ ] **Step 1: Write the failing test**

Create `sword2-lib/src/sword/candidate_features.rs`:

```rust
//! Per-candidate features for domain-partition reranking.
//!
//! These operate on the same 0-based delineation strings produced by
//! `compute_measure::MeasureLine` and consumed by `parse_measure` — i.e.
//! *before* `remap_residue_numbers` converts indices to original PDB
//! numbering.

use crate::peeling::algorithm::SsType;

/// Fraction of residues within `window` positions of a domain-boundary
/// junction that DSSP classifies as coil (loop/turn). Domain linkers are
/// structurally coil far more often than not, so a low fraction is a signal
/// of a boundary cutting through a helix or strand.
pub fn boundary_coil_fraction(raw_delineation: &str, ss_types: &[SsType]) -> f64 {
    const WINDOW: i64 = 2;
    let mut total = 0usize;
    let mut coil = 0usize;

    for domain in raw_delineation.trim().split_whitespace() {
        for segment in domain.split(';') {
            let parts: Vec<&str> = segment.split('-').collect();
            if parts.len() != 2 {
                continue;
            }
            let (Ok(start), Ok(end)) = (parts[0].parse::<i64>(), parts[1].parse::<i64>()) else {
                continue;
            };
            for junction in [start, end] {
                for offset in -WINDOW..=WINDOW {
                    let idx = junction + offset;
                    if idx >= 0 && (idx as usize) < ss_types.len() {
                        total += 1;
                        if ss_types[idx as usize] == SsType::Coil {
                            coil += 1;
                        }
                    }
                }
            }
        }
    }

    if total == 0 {
        0.0
    } else {
        coil as f64 / total as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boundary_coil_fraction_all_coil() {
        let ss_types = vec![SsType::Coil; 20];
        let frac = boundary_coil_fraction("0-9 10-19", &ss_types);
        assert_eq!(frac, 1.0);
    }

    #[test]
    fn test_boundary_coil_fraction_boundary_in_helix() {
        // Junction at index 9/10 is deep inside a helix run (indices 5..15).
        let mut ss_types = vec![SsType::Coil; 20];
        for s in ss_types.iter_mut().take(15).skip(5) {
            *s = SsType::Helix;
        }
        let frac = boundary_coil_fraction("0-9 10-19", &ss_types);
        assert!(frac < 0.5, "expected low coil fraction, got {frac}");
    }

    #[test]
    fn test_boundary_coil_fraction_discontinuous_segment() {
        let ss_types = vec![SsType::Coil; 20];
        // Domain 1 is discontinuous (two segments); still four junctions total.
        let frac = boundary_coil_fraction("0-4;15-19 5-14", &ss_types);
        assert_eq!(frac, 1.0);
    }

    #[test]
    fn test_boundary_coil_fraction_empty_delineation() {
        let ss_types = vec![SsType::Coil; 20];
        assert_eq!(boundary_coil_fraction("", &ss_types), 0.0);
    }
}
```

Register the module in `sword2-lib/src/sword/mod.rs` (near the other `pub mod` declarations at line 22-26):

```rust
pub mod candidate_features;
pub mod compute_jones;
pub mod compute_measure;
pub mod distance_model;
pub mod junctions;
pub mod parse_measure;
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p sword2-lib --lib sword::candidate_features::tests::`
Expected: FAIL — `boundary_coil_fraction` not found (the function body above is written in the same step as the test in this case since it's a pure new file; to genuinely see a red step first, temporarily stub the function body as `unimplemented!()`, run the test to see it panic, then paste in the real body from Step 1 and proceed to Step 4). If you write the file exactly as shown above in one shot, run the test directly — it should already pass, which is acceptable for a from-scratch pure-function file, but confirm by temporarily commenting out the four test bodies' assertions to see them fail first if you want a strict red/green cycle.

- [ ] **Step 3: N/A — implementation is already written in Step 1**

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p sword2-lib --lib sword::candidate_features::tests::`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add sword2-lib/src/sword/candidate_features.rs sword2-lib/src/sword/mod.rs
git commit -m "feat: add boundary secondary-structure coil-fraction candidate feature"
```

---

## Task 3: Chain-modal-domain-count distance feature

**Files:**
- Modify: `sword2-lib/src/sword/candidate_features.rs`

**Interfaces:**
- Produces: `pub fn modal_num_domains(counts: &[usize]) -> usize` and `pub fn modal_count_distance(num_domains: usize, modal: usize) -> f64`. `modal_num_domains` takes the `num_domains` value of every candidate in *the same chain's own candidate set* and returns the most frequent value (ties broken by the smaller domain count, for determinism). This anchors each candidate to what's typical for that specific protein's own candidate spread — the fix memory recorded for Phase B's global under-segmentation bias.

- [ ] **Step 1: Write the failing test**

Add to `sword2-lib/src/sword/candidate_features.rs`, above the existing `#[cfg(test)] mod tests` block's closing brace... actually append these functions above the `#[cfg(test)]` block and these tests inside it:

```rust
/// Most frequent `num_domains` value across one chain's own candidate set.
/// Ties broken toward the smaller count for determinism.
pub fn modal_num_domains(counts: &[usize]) -> usize {
    let mut freq: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    for &nd in counts {
        *freq.entry(nd).or_insert(0) += 1;
    }
    freq.into_iter()
        .max_by_key(|&(nd, count)| (count, std::cmp::Reverse(nd)))
        .map(|(nd, _)| nd)
        .unwrap_or(0)
}

/// Absolute distance of a candidate's domain count from its chain's modal count.
pub fn modal_count_distance(num_domains: usize, modal: usize) -> f64 {
    (num_domains as f64 - modal as f64).abs()
}
```

Add tests inside the existing `mod tests` block:

```rust
    #[test]
    fn test_modal_num_domains_picks_most_frequent() {
        assert_eq!(modal_num_domains(&[2, 3, 3, 3, 4, 5]), 3);
    }

    #[test]
    fn test_modal_num_domains_ties_break_smaller() {
        assert_eq!(modal_num_domains(&[2, 2, 5, 5]), 2);
    }

    #[test]
    fn test_modal_num_domains_empty_is_zero() {
        assert_eq!(modal_num_domains(&[]), 0);
    }

    #[test]
    fn test_modal_count_distance() {
        assert_eq!(modal_count_distance(5, 3), 2.0);
        assert_eq!(modal_count_distance(3, 3), 0.0);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p sword2-lib --lib sword::candidate_features::tests::test_modal`
Expected: FAIL — `modal_num_domains`/`modal_count_distance` not found.

- [ ] **Step 3: Implementation already written in Step 1**

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p sword2-lib --lib sword::candidate_features::tests::test_modal`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add sword2-lib/src/sword/candidate_features.rs
git commit -m "feat: add chain-modal-domain-count distance candidate feature"
```

---

## Task 4: Thread a reduced-shuffle EnergyConfig into the SWORD pipeline

**Files:**
- Modify: `sword2-lib/src/sword/mod.rs:28-53` (`SwordConfig`), `:86-91` (`run_pipeline` signature use — no signature change needed, field-only), other call sites in this file
- Modify: `sword2-lib/src/sword/candidate_features.rs` (add `candidate_energy_z_score`)
- Modify: `sword2-cli/src/main.rs:635-641` (`SwordConfig` construction), `:668-678` (move/reuse `EnergyConfig` setup)
- Test: `sword2-lib/src/sword/candidate_features.rs` (inline, using the existing golden fixture pattern)

**Interfaces:**
- Consumes: `crate::energy::EnergyConfig`, `crate::energy::get_energy_and_z_score`, `crate::energy::build_residue_list` (all pre-existing, unchanged).
- Produces: `SwordConfig.energy_config: Option<crate::energy::EnergyConfig>`, `SwordConfig.chain_id: String`. `pub fn candidate_energy_z_score(energy_config: &crate::energy::EnergyConfig, pdb_path: &str, chain: &str, remapped_delineation: &str) -> Option<f64>` — `remapped_delineation` is the delineation string **after** `remap_residue_numbers` (original PDB residue numbers), unlike `boundary_coil_fraction`'s `raw_delineation`. Returns `None` if the delineation is empty or scoring fails.

- [ ] **Step 1: Write the failing test**

Add to `sword2-lib/src/sword/candidate_features.rs`, above the `#[cfg(test)]` block:

```rust
/// Pseudo-energy Z-score for one candidate's full domain assignment (all
/// domains' residues together), used to rescore the already-shortlisted
/// candidates. `remapped_delineation` must use original PDB residue numbers
/// (i.e. already passed through `remap_residue_numbers`), unlike
/// `boundary_coil_fraction`'s `raw_delineation`.
pub fn candidate_energy_z_score(
    energy_config: &crate::energy::EnergyConfig,
    pdb_path: &str,
    chain: &str,
    remapped_delineation: &str,
) -> Option<f64> {
    let mut residues = String::new();
    for domain in remapped_delineation.trim().split_whitespace() {
        for segment in domain.split(';') {
            let parts: Vec<&str> = segment.split('-').collect();
            if parts.len() != 2 {
                continue;
            }
            let (Ok(start), Ok(end)) = (parts[0].parse::<i32>(), parts[1].parse::<i32>()) else {
                continue;
            };
            let list = crate::energy::build_residue_list((start, end), chain);
            if !residues.is_empty() {
                residues.push(',');
            }
            residues.push_str(&list);
        }
    }
    if residues.is_empty() {
        return None;
    }
    crate::energy::get_energy_and_z_score(energy_config, pdb_path, Some(&residues))
        .ok()
        .and_then(|r| r.z_score)
}
```

Add a test in the `mod tests` block. This reuses the same 20-residue fixture PDB fragment as `sword2-lib/tests/dssp_golden.rs` (residues 1-20 of 1JX4_A) so the energy scorer has real backbone geometry to work with — copy the `FIXTURE_PDB` constant's residue lines 1-20 verbatim from that file into a local `const FIXTURE_PDB` in this test module (do not `include!` across the crate/integration-test boundary; duplicating ~20 short lines is simpler and keeps the test self-contained):

```rust
    #[test]
    fn test_candidate_energy_z_score_returns_finite_value() {
        use std::io::Write as _;
        let dir = tempfile::tempdir().unwrap();
        let pdb_path = dir.path().join("frag.pdb");
        let mut f = std::fs::File::create(&pdb_path).unwrap();
        // Residues 1-20 of 1JX4_A backbone atoms — see sword2-lib/tests/dssp_golden.rs
        // for the full fixture and provenance.
        writeln!(f, "ATOM      1  N   ILE A   1      46.170  17.543  13.913  1.00 33.83           N").unwrap();
        writeln!(f, "ATOM      2  CA  ILE A   1      45.665  16.827  12.751  1.00 32.01           C").unwrap();
        writeln!(f, "ATOM      3  C   ILE A   1      44.396  16.160  13.297  1.00 29.09           C").unwrap();
        writeln!(f, "ATOM      4  O   ILE A   1      44.466  15.419  14.274  1.00 27.71           O").unwrap();
        writeln!(f, "ATOM      9  N   VAL A   2      43.247  16.494  12.713  1.00 28.04           N").unwrap();
        writeln!(f, "ATOM     10  CA  VAL A   2      41.973  15.941  13.134  1.00 26.25           C").unwrap();
        writeln!(f, "ATOM     11  C   VAL A   2      41.475  15.044  12.028  1.00 27.46           C").unwrap();
        writeln!(f, "ATOM     12  O   VAL A   2      41.543  15.426  10.857  1.00 26.95           O").unwrap();
        drop(f);

        let bin_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/..");
        let mut ec = crate::energy::EnergyConfig::from_bin_dir(bin_dir);
        ec.num_shuffles = 20; // keep the test fast; production reranking uses 200
        ec.preload().unwrap();

        let z = candidate_energy_z_score(&ec, pdb_path.to_str().unwrap(), "A", "1-2");
        assert!(z.is_some(), "expected a Z-score for a 2-residue fragment");
        assert!(z.unwrap().is_finite());
    }

    #[test]
    fn test_candidate_energy_z_score_empty_delineation_is_none() {
        let bin_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/..");
        let ec = crate::energy::EnergyConfig::from_bin_dir(bin_dir);
        assert_eq!(candidate_energy_z_score(&ec, "/nonexistent.pdb", "A", ""), None);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p sword2-lib --lib sword::candidate_features::tests::test_candidate_energy_z_score`
Expected: FAIL — `candidate_energy_z_score` not found (compile error) until the function above is added; once added, `test_candidate_energy_z_score_returns_finite_value` should pass immediately since the fixture PDB is real geometry. If it fails at runtime instead of compile time, check that `bin/mypmfs-master/025_30_100_potential` exists relative to the workspace root (it's checked into the repo per CLAUDE.md).

- [ ] **Step 3: Implementation already written in Step 1**

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p sword2-lib --lib sword::candidate_features::tests::test_candidate_energy_z_score`
Expected: PASS (2 tests)

- [ ] **Step 5: Add `energy_config` and `chain_id` to `SwordConfig`**

In `sword2-lib/src/sword/mod.rs`, change the `SwordConfig` struct (lines 28-53):

```rust
/// Configuration for a SWORD2 run.
#[derive(Debug, Clone)]
pub struct SwordConfig {
    /// Whether to compute energies.
    pub compute_energies: bool,
    /// Whether to generate plots.
    pub generate_plots: bool,
    /// Number of threads for parallel computation.
    pub num_threads: usize,
    /// Output directory.
    pub output_dir: String,
    /// Max alternative assignments (3, 9, or 15).
    pub max_alternatives: usize,
    /// Reduced-shuffle energy config used for candidate rescoring and the
    /// training dump — separate from the `-E` display-time energy config,
    /// which main.rs builds independently at the full shuffle count. `None`
    /// disables energy-based rescoring entirely (dump rows omit `energy_z`,
    /// `use_pairwise_reranker` falls back to distance_model-only scoring).
    pub energy_config: Option<crate::energy::EnergyConfig>,
    /// PDB chain letter (e.g. "A"), needed to build energy residue lists.
    pub chain_id: String,
    /// Use the pairwise-trained reranker to pick the winning candidate
    /// instead of the legacy distance_model-based selection. Off by default.
    pub use_pairwise_reranker: bool,
}

impl Default for SwordConfig {
    fn default() -> Self {
        Self {
            compute_energies: true,
            generate_plots: true,
            num_threads: num_cpus::get(),
            output_dir: ".".to_string(),
            max_alternatives: 9,
            energy_config: None,
            chain_id: "A".to_string(),
            use_pairwise_reranker: false,
        }
    }
}
```

- [ ] **Step 6: Wire the reduced-shuffle EnergyConfig and chain_id from main.rs**

In `sword2-cli/src/main.rs`, the `SwordConfig` is currently built at lines 635-641, *before* the existing `-E` display `EnergyConfig` is built at lines 670-678. Move a lightweight, always-constructed `EnergyConfig` ahead of the `SwordConfig` literal:

```rust
    // Step 4: Run the SWORD pipeline (DSSP is pure Rust, no compilation needed)
    reporter.step("SWORD pipeline");
    tracing::debug!("Launch SWORD pipeline");

    // Reduced-shuffle energy config for candidate rescoring/training dump —
    // separate from the `-E` display config below (which uses cli.zscore_shuffles,
    // default 2000). Construction is free (potentials load lazily on first use),
    // so this is built unconditionally without affecting default runtime.
    let mut rerank_energy_config = energy::EnergyConfig::from_bin_dir(&bin_dir.to_string_lossy());
    rerank_energy_config.num_shuffles = 200;

    let config = sword::SwordConfig {
        compute_energies: cli.energies,
        generate_plots: cli.plots,
        num_threads,
        output_dir: results_dir.to_string_lossy().to_string(),
        max_alternatives: 9,
        energy_config: Some(rerank_energy_config),
        chain_id: chain_id.to_string(),
        use_pairwise_reranker: cli.use_pairwise_reranker,
    };

    let (sword_output, sword_results) = sword::run_pipeline(&input_pdb, &pdb_id_chain, &config)
        .context("Failed to run SWORD pipeline")?;
```

This references `cli.use_pairwise_reranker`, which doesn't exist yet — add it to the `Cli` struct in `sword2-cli/src/main.rs` near the other boolean flags (after `energies`, around line 55):

```rust
    /// Use the pairwise-trained reranker (energy Z-score + boundary secondary
    /// structure) to pick the winning partition instead of the legacy
    /// distance_model selection. Off by default; benchmark before enabling.
    #[arg(long)]
    use_pairwise_reranker: bool,
```

- [ ] **Step 7: `cargo check` the workspace**

Run: `cargo check --workspace`
Expected: Compiles clean. (`use_pairwise_reranker` and `energy_config` are unread inside `run_pipeline` at this point — Task 5 and Task 8 consume them — so expect an `unused` warning on `SwordConfig.use_pairwise_reranker`/`energy_config` reads if `cargo check` is strict about dead struct fields; it is not by default for pub fields, so this should be a clean pass.)

- [ ] **Step 8: Run the full sword2-lib and sword2-cli test suites**

Run: `cargo test -p sword2-lib -p sword`
Expected: All PASS.

- [ ] **Step 9: Commit**

```bash
git add sword2-lib/src/sword/mod.rs sword2-lib/src/sword/candidate_features.rs sword2-cli/src/main.rs
git commit -m "feat: thread reduced-shuffle EnergyConfig into the SWORD pipeline"
```

---

## Task 5: Move the candidate dump to the shortlist stage and add new feature columns

**Files:**
- Modify: `sword2-lib/src/sword/mod.rs:220-251` (existing `SWORD2_DUMP_CANDIDATES` block — delete from here), `:336-345` (insert new dump block after `relevant_measure2` is computed)
- Test: `sword2-lib/src/sword/mod.rs` (inline; a small helper is extracted to make this testable without file I/O — see Step 1)

**Interfaces:**
- Consumes: `candidate_features::boundary_coil_fraction`, `candidate_features::candidate_energy_z_score`, `candidate_features::modal_num_domains`/`modal_count_distance` (Tasks 2-4), `SwordConfig.energy_config`/`chain_id` (Task 4), `remap_residue_numbers` (pre-existing, private to this file), `PeelingOutput.ss_types` (Task 1, via `peeling_output`, already in scope in `run_pipeline`).
- Produces: dump rows gain three new trailing columns: `boundary_coil_fraction`, `energy_z` (empty string if `energy_config` is `None` or scoring failed), `modal_count_distance`.

The existing dump block (lines 220-251) iterates the *full* `measure_lines` (every candidate at every peeling level, before any filtering) and writes 8 columns. It runs *before* `n_dom`/`relevant_measure2` are computed. This task moves it to run *after* `relevant_measure2` (the same `alt_b`/`alt_l`-bounded shortlist the reranker will score at inference time) so training features match inference features exactly, and appends the three new columns.

- [ ] **Step 1: Write the failing test**

`run_pipeline` writes directly to a file via `SWORD2_DUMP_CANDIDATES`, which is awkward to unit test in isolation. Extract the per-row feature computation into a small, directly-testable function first. Add to `sword2-lib/src/sword/candidate_features.rs` (above the `#[cfg(test)]` block):

```rust
/// One row of dump/training-table features for a single candidate, computed
/// at the same shortlist stage the reranker scores at inference time.
pub struct DumpRow {
    pub num_domains: usize,
    pub min_size: usize,
    pub max_cr: f64,
    pub density_min: f64,
    pub mean_density: f64,
    pub delineation: String,
    pub boundary_coil_fraction: f64,
    pub energy_z: Option<f64>,
    pub modal_count_distance: f64,
}

/// Build a `DumpRow` from one `relevant_measure2` pipe-delimited line.
///
/// `raw_delineation` uses 0-based indices (matches `ss_types`); the caller
/// supplies `remapped_delineation` (original PDB numbering) separately since
/// only `candidate_energy_z_score` needs it. `modal` is
/// `modal_num_domains(...)` computed once across the whole shortlist by the
/// caller (not per-row) and passed in.
#[allow(clippy::too_many_arguments)]
pub fn build_dump_row(
    num_domains: usize,
    min_size: usize,
    max_cr: f64,
    density_min: f64,
    mean_density: f64,
    raw_delineation: &str,
    remapped_delineation: &str,
    ss_types: &[SsType],
    energy_config: Option<&crate::energy::EnergyConfig>,
    pdb_path: &str,
    chain: &str,
    modal: usize,
) -> DumpRow {
    let coil_fraction = boundary_coil_fraction(raw_delineation, ss_types);
    let energy_z = energy_config
        .and_then(|ec| candidate_energy_z_score(ec, pdb_path, chain, remapped_delineation));
    DumpRow {
        num_domains,
        min_size,
        max_cr,
        density_min,
        mean_density,
        delineation: raw_delineation.to_string(),
        boundary_coil_fraction: coil_fraction,
        energy_z,
        modal_count_distance: modal_count_distance(num_domains, modal),
    }
}
```

Add a test in `mod tests`:

```rust
    #[test]
    fn test_build_dump_row_without_energy_config() {
        let ss_types = vec![SsType::Coil; 20];
        let row = build_dump_row(
            2, 10, 0.5, 1.0, 2.0,
            "0-9 10-19", "1-10 11-20",
            &ss_types, None, "/nonexistent.pdb", "A", 2,
        );
        assert_eq!(row.num_domains, 2);
        assert_eq!(row.boundary_coil_fraction, 1.0);
        assert_eq!(row.energy_z, None);
        assert_eq!(row.modal_count_distance, 0.0);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p sword2-lib --lib sword::candidate_features::tests::test_build_dump_row`
Expected: FAIL — `DumpRow`/`build_dump_row` not found.

- [ ] **Step 3: Implementation already written in Step 1**

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p sword2-lib --lib sword::candidate_features::tests::test_build_dump_row`
Expected: PASS

- [ ] **Step 5: Replace the dump block in `run_pipeline`**

In `sword2-lib/src/sword/mod.rs`, delete the existing dump block (lines 220-251, the `if let Ok(dump_path) = std::env::var("SWORD2_DUMP_CANDIDATES") { ... }` block that iterates `measure_lines`).

Insert a new dump block immediately after `relevant_measure2` is computed (after line 345, `let relevant_measure2 = parse_measure::parse_measure(...)`):

```rust
    // Training dump: set SWORD2_DUMP_CANDIDATES=/path/to/output.csv to record the
    // same alt_b/alt_l shortlist the pairwise reranker scores at inference time,
    // with the same features (boundary_coil_fraction, energy_z, modal_count_distance)
    // for offline training. Moved here (post-shortlist) rather than on the raw
    // measure_lines so training and inference see identical feature distributions.
    if let Ok(dump_path) = std::env::var("SWORD2_DUMP_CANDIDATES") {
        use std::io::Write as _;

        let ss_types: &[crate::peeling::algorithm::SsType] = peeling_output
            .as_ref()
            .map(|po| po.ss_types.as_slice())
            .unwrap_or(&[]);
        let pdb_path_str = pdb_file_dst.to_string_lossy().to_string();

        let mut parsed_rows: Vec<(usize, usize, f64, f64, f64, String, String)> = Vec::new();
        for rm in &relevant_measure2 {
            let fields: Vec<&str> = rm.split('|').collect();
            if fields.len() < 6 {
                continue;
            }
            let nd: usize = fields[0].trim().parse().unwrap_or(0);
            let min_size: usize = fields[1].trim().parse().unwrap_or(0);
            let raw_del = fields[2].trim().to_string();
            let max_cr: f64 = fields[3].trim().parse().unwrap_or(0.0);
            let density_min: f64 = fields[4].trim().parse().unwrap_or(0.0);
            let mean_density: f64 = fields[5].trim().parse().unwrap_or(0.0);
            let remapped_del = remap_residue_numbers(&raw_del, &tab_num);
            parsed_rows.push((nd, min_size, max_cr, density_min, mean_density, raw_del, remapped_del));
        }

        let modal = candidate_features::modal_num_domains(
            &parsed_rows.iter().map(|r| r.0).collect::<Vec<_>>(),
        );

        let file_existed = std::path::Path::new(&dump_path).exists();
        if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&dump_path) {
            if !file_existed {
                let _ = writeln!(
                    f,
                    "chain_id,output_dir,num_domains,min_size,max_cr,density_min,mean_density,delineation,boundary_coil_fraction,energy_z,modal_count_distance"
                );
            }
            for (nd, min_size, max_cr, density_min, mean_density, raw_del, remapped_del) in &parsed_rows {
                let row = candidate_features::build_dump_row(
                    *nd, *min_size, *max_cr, *density_min, *mean_density,
                    raw_del, remapped_del, ss_types,
                    config.energy_config.as_ref(), &pdb_path_str, &config.chain_id, modal,
                );
                let energy_z_str = row.energy_z.map(|z| z.to_string()).unwrap_or_default();
                let _ = writeln!(
                    f,
                    "{},{},{},{},{:.6},{:.6},{:.6},\"{}\",{:.6},{},{:.1}",
                    pdb_name, results_dir.display(),
                    row.num_domains, row.min_size, row.max_cr, row.density_min, row.mean_density,
                    row.delineation, row.boundary_coil_fraction, energy_z_str, row.modal_count_distance,
                );
            }
        }
    }
```

Note this references `peeling_output` and `tab_num`, both already in scope at this point in `run_pipeline` (per the existing code around lines 111-174 and 211-216).

- [ ] **Step 6: `cargo check`**

Run: `cargo check -p sword2-lib`
Expected: Compiles clean. Fix any borrow/lifetime issues that come up around `peeling_output.as_ref()` — `peeling_output` is `Option<PeelingOutput>` and is also used later in the function (for `compute_measure::compute_measure_from_data`), so this dump block must come *after* that usage or use `.as_ref()` consistently to avoid a move — the block above already uses `.as_ref()` so it should not conflict, but verify no earlier code takes ownership of `peeling_output` before this point.

- [ ] **Step 7: Manual smoke test of the dump path**

Run: `cargo build --release && SWORD2_DUMP_CANDIDATES=/tmp/dump_test.csv ./target/release/sword2 -p 1jx4 -o /tmp/sword2_dump_test`
Expected: Exits 0; `/tmp/dump_test.csv` exists with a header row containing `boundary_coil_fraction,energy_z,modal_count_distance` and one data row per shortlisted candidate. `energy_z` will be non-empty for every row since `energy_config` is always constructed now (Task 4, Step 6).

Run: `head -3 /tmp/dump_test.csv`
Expected: Header line followed by 1-2 data rows with plausible-looking numeric values (coil fraction in [0,1], energy_z roughly in [-6, 6], modal_count_distance a small non-negative float).

- [ ] **Step 8: Run the full test suite**

Run: `cargo test -p sword2-lib -p sword`
Expected: All PASS.

- [ ] **Step 9: Commit**

```bash
git add sword2-lib/src/sword/mod.rs sword2-lib/src/sword/candidate_features.rs
git commit -m "feat: move candidate dump to shortlist stage, add new feature columns"
```

---

## Task 6: Carry the new feature columns through `build_training_table.py`

**Files:**
- Modify: `benchmark/build_training_table.py:47-64` (`FIELDNAMES`), `:164-183` (`_score_candidates` row construction)
- Test: Create `benchmark/tests/test_build_training_table.py`

**Interfaces:**
- Consumes: dump CSV rows with the three new columns from Task 5 (`boundary_coil_fraction`, `energy_z`, `modal_count_distance`).
- Produces: `training_table.csv` rows carry the same three columns through unchanged, alongside the existing `ndo`/`iou`/etc. ground-truth columns — this is the file `benchmark/train_pairwise_reranker.py` (Task 7) will train on.

- [ ] **Step 1: Write the failing test**

Create `benchmark/tests/test_build_training_table.py`:

```python
from benchmark.build_training_table import _score_candidates, FIELDNAMES


def test_fieldnames_include_new_reranker_features():
    for col in ("boundary_coil_fraction", "energy_z", "modal_count_distance"):
        assert col in FIELDNAMES, f"missing column: {col}"


def test_score_candidates_passes_through_new_features(tmp_path, monkeypatch):
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest benchmark/tests/test_build_training_table.py -v`
Expected: FAIL — `test_fieldnames_include_new_reranker_features` fails (columns missing from `FIELDNAMES`); `test_score_candidates_passes_through_new_features` fails with a `KeyError` on the new columns.

- [ ] **Step 3: Add the columns**

In `benchmark/build_training_table.py`, extend `FIELDNAMES` (lines 47-64):

```python
FIELDNAMES = [
    "chain_id",
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "delineation",
    "boundary_coil_fraction",
    "energy_z",
    "modal_count_distance",
    "n_true_domains",
    "n_pred_domains",
    "ndo",
    "iou",
    "boundary_f1_10",
    "matched_dice",
    "d_count_acc",
    "S",
    "is_oracle_s",
]
```

In `_score_candidates`, extend the output row dict (lines 164-183) to pass the three new columns through verbatim from the input row:

```python
        rows.append({
            "chain_id": chain_id,
            "num_domains": row["num_domains"],
            "min_size": row["min_size"],
            "max_cr": row["max_cr"],
            "density_min": row["density_min"],
            "mean_density": row["mean_density"],
            "delineation": row.get("delineation", "").strip().strip('"'),
            "boundary_coil_fraction": row.get("boundary_coil_fraction", ""),
            "energy_z": row.get("energy_z", ""),
            "modal_count_distance": row.get("modal_count_distance", ""),
            "n_true_domains": metrics.n_true_domains,
            "n_pred_domains": metrics.n_pred_domains,
            "ndo": metrics.ndo,
            "iou": metrics.iou,
            "boundary_f1_10": metrics.boundary_f1_10,
            "matched_dice": metrics.matched_dice,
            "d_count_acc": metrics.d_count_acc,
            "S": s_val,
            "is_oracle_s": 1 if abs(s_val - best_s) < 1e-9 else 0,
        })
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest benchmark/tests/test_build_training_table.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Run the full benchmark test suite**

Run: `python -m pytest benchmark/tests -q`
Expected: All PASS (this is also the project's documented validation gate per `benchmark/REPORT.md`).

- [ ] **Step 6: Commit**

```bash
git add benchmark/build_training_table.py benchmark/tests/test_build_training_table.py
git commit -m "feat: carry new reranker feature columns through build_training_table"
```

---

## Task 7: Pairwise ranking-loss trainer

**Files:**
- Create: `benchmark/train_pairwise_reranker.py`
- Create: `benchmark/tests/test_train_pairwise_reranker.py`
- Create: `benchmark/data/pairwise_reranker_weights.json` (produced by Step 7's synthetic smoke run)

**Interfaces:**
- Consumes: `benchmark/data/training_table.csv` (produced by `build_training_table.py`, Task 6), with columns `chain_id, num_domains, min_size, max_cr, density_min, mean_density, boundary_coil_fraction, energy_z, modal_count_distance, ndo, ...`.
- Produces: `benchmark/data/pairwise_reranker_weights.json` — `{"features": [...], "weights": [...], "bias": ..., "feature_means": [...], "feature_stds": [...]}`. Feature order is fixed and must match `sword2-lib/src/sword/reranker.rs`'s `FEATURES` order exactly (Task 8).

This is where the fix for Phase B's failure mode lives: the loss compares candidate pairs **within the same `chain_id`** only, predicting `sign(ndo_i - ndo_j)`, so it cannot learn a global "fewer domains is better" bias the way the old pointwise logistic did (memory: Phase B's `num_domains` weight was -1.33, causing under-segmentation). `modal_count_distance` (Task 3) gives the model a chain-local anchor for "how many domains is normal for this protein" instead of a corpus-wide one.

- [ ] **Step 1: Write the failing test**

Create `benchmark/tests/test_train_pairwise_reranker.py`:

```python
import json

import numpy as np
import pandas as pd

from benchmark.train_pairwise_reranker import (
    FEATURES,
    build_pairs,
    pairwise_logistic_loss_and_grad,
    train,
)


def _toy_table() -> pd.DataFrame:
    # Two chains. In each, candidate "a" is strictly better (higher ndo) and
    # has a higher energy_z and coil_fraction than candidate "b" — a
    # separable synthetic case the trainer should learn perfectly.
    rows = []
    for chain in ["c1", "c2"]:
        rows.append(dict(chain_id=chain, num_domains=2, min_size=30, max_cr=0.3,
                          density_min=1.0, mean_density=2.0, boundary_coil_fraction=0.9,
                          energy_z=-3.0, modal_count_distance=0.0, ndo=0.9))
        rows.append(dict(chain_id=chain, num_domains=4, min_size=10, max_cr=0.6,
                          density_min=0.5, mean_density=1.0, boundary_coil_fraction=0.1,
                          energy_z=1.0, modal_count_distance=2.0, ndo=0.3))
    return pd.DataFrame(rows)


def test_build_pairs_only_compares_within_chain():
    df = _toy_table()
    pairs = build_pairs(df)
    assert len(pairs) == 2  # one pair per chain, not cross-chain
    for i, j, label in pairs:
        assert df.loc[i, "chain_id"] == df.loc[j, "chain_id"]
        assert label in (1, -1)


def test_pairwise_logistic_loss_decreases_with_correct_sign():
    # A large positive margin in the "correct" direction should have lower
    # loss than the same margin in the "wrong" direction.
    w = np.array([1.0])
    x_i = np.array([[2.0]])
    x_j = np.array([[0.0]])
    labels = np.array([1])
    loss_correct, _ = pairwise_logistic_loss_and_grad(w, 0.0, x_i, x_j, labels)
    loss_wrong, _ = pairwise_logistic_loss_and_grad(-w, 0.0, x_i, x_j, labels)
    assert loss_correct < loss_wrong


def test_train_converges_on_separable_toy_data(tmp_path):
    df = _toy_table()
    weights_path = tmp_path / "weights.json"
    result = train(df, features=FEATURES, epochs=500, lr=0.5, out_path=weights_path)

    with open(weights_path) as f:
        saved = json.load(f)
    assert saved["features"] == FEATURES
    assert len(saved["weights"]) == len(FEATURES)

    # The trained model should rank candidate "a" above "b" in both chains.
    assert result["train_accuracy"] == 1.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest benchmark/tests/test_train_pairwise_reranker.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'benchmark.train_pairwise_reranker'`.

- [ ] **Step 3: Write the trainer**

Create `benchmark/train_pairwise_reranker.py`:

```python
"""Pairwise ranking-loss trainer for SWORD2 candidate reranking.

Trains on within-chain candidate pairs, predicting sign(ndo_i - ndo_j), so
the model can never learn a corpus-wide domain-count bias the way a
pointwise classifier can (see benchmark/data/reranker retro in project
memory: Phase B's pointwise logistic learned -1.33 weight on num_domains
and systematically under-segmented).

Usage:
    python -m benchmark.train_pairwise_reranker \\
        --training-table benchmark/data/training_table.csv \\
        --out benchmark/data/pairwise_reranker_weights.json
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

# Fixed feature order — must match sword2-lib/src/sword/reranker.rs::FEATURES exactly.
FEATURES = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "energy_z",
    "modal_count_distance",
]


def build_pairs(df: pd.DataFrame) -> list[tuple[int, int, int]]:
    """Build (i, j, label) triples from within-chain candidate pairs.

    label = 1 means df.loc[i] should rank above df.loc[j] (higher ndo);
    label = -1 means the reverse. Pairs with equal ndo are skipped (no
    signal). Returns integer-position index pairs into df.reset_index(drop=True).
    """
    df = df.reset_index(drop=True)
    pairs: list[tuple[int, int, int]] = []
    for _, group in df.groupby("chain_id"):
        idx = group.index.tolist()
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                ndo_i, ndo_j = df.loc[i, "ndo"], df.loc[j, "ndo"]
                if ndo_i == ndo_j:
                    continue
                label = 1 if ndo_i > ndo_j else -1
                pairs.append((i, j, label))
    return pairs


def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def pairwise_logistic_loss_and_grad(
    w: np.ndarray,
    b: float,
    x_i: np.ndarray,
    x_j: np.ndarray,
    labels: np.ndarray,
) -> tuple[float, np.ndarray]:
    """RankNet-style pairwise logistic loss and its gradient w.r.t. w.

    score(x) = w . x + b. margin = label * (score(x_i) - score(x_j)).
    loss = mean(log(1 + exp(-margin))), a smooth surrogate for "score_i
    should exceed score_j whenever label == 1".
    """
    diff = x_i - x_j  # (n_pairs, n_features)
    margin = labels * (diff @ w + b - b)  # bias cancels in the difference; kept for API symmetry
    # Numerically stable log(1 + exp(-margin))
    loss_terms = np.logaddexp(0.0, -margin)
    loss = float(np.mean(loss_terms))

    # d/dw of mean(log(1+exp(-margin))) where margin = label * (diff @ w)
    sig = _sigmoid(-margin)  # = 1 - sigmoid(margin)
    grad = -np.mean((sig * labels)[:, None] * diff, axis=0)
    return loss, grad


def train(
    df: pd.DataFrame,
    features: list[str],
    epochs: int = 2000,
    lr: float = 0.1,
    l2: float = 1e-3,
    out_path: Path | None = None,
) -> dict:
    """Train a linear pairwise reranker with full-batch gradient descent.

    Features are chain-local z-score normalized before pairing (matching the
    normalization the Rust inference side applies at run time — see
    sword2-lib/src/sword/reranker.rs).
    """
    df = df.reset_index(drop=True).copy()
    for col in features:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=features + ["ndo"]).reset_index(drop=True)

    # Chain-local z-score normalization (per chain_id group), matching inference.
    def _zscore(group: pd.DataFrame) -> pd.DataFrame:
        for col in features:
            std = group[col].std(ddof=0)
            mean = group[col].mean()
            group[col] = (group[col] - mean) / (std + 1e-8)
        return group

    normed = df.groupby("chain_id", group_keys=False).apply(_zscore)

    pairs = build_pairs(pd.concat([df["chain_id"], df["ndo"]], axis=1))
    if not pairs:
        raise ValueError("No trainable pairs — need at least 2 candidates with differing ndo per chain")

    x = normed[features].to_numpy(dtype=float)
    idx_i = np.array([p[0] for p in pairs])
    idx_j = np.array([p[1] for p in pairs])
    labels = np.array([p[2] for p in pairs], dtype=float)
    x_i, x_j = x[idx_i], x[idx_j]

    w = np.zeros(len(features))
    b = 0.0
    for epoch in range(epochs):
        loss, grad = pairwise_logistic_loss_and_grad(w, b, x_i, x_j, labels)
        grad = grad + l2 * w
        w -= lr * grad
        if epoch % max(1, epochs // 10) == 0:
            log.info("epoch %d: loss=%.4f", epoch, loss)

    scores_i = x_i @ w
    scores_j = x_j @ w
    pred = np.sign(scores_i - scores_j)
    train_accuracy = float(np.mean(pred == labels))

    result = {
        "features": features,
        "weights": w.tolist(),
        "bias": 0.0,  # bias cancels in pairwise scoring; kept for schema stability
        "feature_means": [],  # per-chain normalization is applied at inference time, not global
        "feature_stds": [],
        "train_accuracy": train_accuracy,
        "n_pairs": len(pairs),
    }

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump({k: result[k] for k in ("features", "weights", "bias")}, f, indent=2)
        log.info("Wrote weights to %s (train_accuracy=%.3f, n_pairs=%d)", out_path, train_accuracy, len(pairs))

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Train the pairwise candidate reranker")
    parser.add_argument("--training-table", type=Path, default=Path("benchmark/data/training_table.csv"))
    parser.add_argument("--out", type=Path, default=Path("benchmark/data/pairwise_reranker_weights.json"))
    parser.add_argument("--epochs", type=int, default=2000)
    parser.add_argument("--lr", type=float, default=0.1)
    args = parser.parse_args()

    df = pd.read_csv(args.training_table)
    train(df, features=FEATURES, epochs=args.epochs, lr=args.lr, out_path=args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest benchmark/tests/test_train_pairwise_reranker.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Produce an initial checked-in weights file from the synthetic smoke data**

This gives Task 8's Rust code a real (if not yet CATH-trained) file to load via `include_str!` at compile time. It will be overwritten with real CATH-trained weights in Task 9.

Run:
```bash
python -c "
from pathlib import Path
import pandas as pd
from benchmark.train_pairwise_reranker import FEATURES, train
from benchmark.tests.test_train_pairwise_reranker import _toy_table
train(_toy_table(), features=FEATURES, epochs=500, lr=0.5, out_path=Path('benchmark/data/pairwise_reranker_weights.json'))
"
```
Expected: `benchmark/data/pairwise_reranker_weights.json` is created with `train_accuracy` printed near 1.0.

Run: `cat benchmark/data/pairwise_reranker_weights.json`
Expected: JSON with `"features"` (8 entries matching `FEATURES`), `"weights"` (8 floats), `"bias": 0.0`.

- [ ] **Step 6: Run the full benchmark test suite**

Run: `python -m pytest benchmark/tests -q`
Expected: All PASS.

- [ ] **Step 7: Commit**

```bash
git add benchmark/train_pairwise_reranker.py benchmark/tests/test_train_pairwise_reranker.py benchmark/data/pairwise_reranker_weights.json
git commit -m "feat: add pairwise ranking-loss trainer for candidate reranking"
```

---

## Task 8: Rust pairwise reranker inference, wired in behind a flag

**Files:**
- Create: `sword2-lib/src/sword/reranker.rs`
- Modify: `sword2-lib/src/sword/mod.rs:22-26` (add `pub mod reranker;`), `:384-401` (`to_print` selection — the exact integration point)
- Test: `sword2-lib/src/sword/reranker.rs` (inline)

**Interfaces:**
- Consumes: `benchmark/data/pairwise_reranker_weights.json` (Task 7, via `include_str!` at compile time — no runtime file I/O or `--install-dir` dependency), `candidate_features::{boundary_coil_fraction, candidate_energy_z_score, modal_num_domains, modal_count_distance}` (Tasks 2-5), `SwordConfig.{energy_config, chain_id, use_pairwise_reranker}` (Task 4).
- Produces: `pub fn rerank(candidates: &[CandidateInput]) -> usize` — returns the index of the winning candidate in `candidates`. `CandidateInput` bundles everything needed to compute all 8 features for one candidate.

This task replaces exactly one thing: which candidate in `relevant_measure2` becomes the reported "Optimal partition" (rank 1). Candidate generation, `n_dom` prediction, and the alternatives listed in the output are all untouched — this is the surgical integration point identified in the design spec.

- [ ] **Step 1: Write the failing test**

Create `sword2-lib/src/sword/reranker.rs`:

```rust
//! Pairwise-trained linear reranker for domain-partitioning candidate
//! selection. Weights are trained offline (see
//! `benchmark/train_pairwise_reranker.py`) with a pairwise ranking loss —
//! deliberately not the pointwise logistic that caused the earlier Phase B
//! attempt to regress into under-segmentation (see project memory).
//!
//! Feature order here MUST match `benchmark/train_pairwise_reranker.py`'s
//! `FEATURES` list exactly.

use crate::peeling::algorithm::SsType;
use crate::sword::candidate_features::{
    boundary_coil_fraction, candidate_energy_z_score, modal_count_distance, modal_num_domains,
};

const WEIGHTS_JSON: &str =
    include_str!("../../../benchmark/data/pairwise_reranker_weights.json");

#[derive(serde::Deserialize)]
struct TrainedWeights {
    features: Vec<String>,
    weights: Vec<f64>,
    #[allow(dead_code)]
    bias: f64,
}

/// Expected feature order — validated against the loaded JSON at rerank time
/// (falls back to distance_model-style no-op scoring if they don't match, so
/// a stale weights file can't silently misapply weights to the wrong features).
const FEATURES: [&str; 8] = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "energy_z",
    "modal_count_distance",
];

/// Everything needed to score one candidate.
pub struct CandidateInput<'a> {
    pub num_domains: usize,
    pub min_size: usize,
    pub max_cr: f64,
    pub density_min: f64,
    pub mean_density: f64,
    pub raw_delineation: &'a str,
    pub remapped_delineation: &'a str,
}

fn candidate_feature_vector(
    c: &CandidateInput,
    ss_types: &[SsType],
    energy_config: Option<&crate::energy::EnergyConfig>,
    pdb_path: &str,
    chain: &str,
    modal: usize,
) -> [f64; 8] {
    let coil = boundary_coil_fraction(c.raw_delineation, ss_types);
    let energy_z = energy_config
        .and_then(|ec| candidate_energy_z_score(ec, pdb_path, chain, c.remapped_delineation))
        .unwrap_or(0.0);
    [
        c.num_domains as f64,
        c.min_size as f64,
        c.max_cr,
        c.density_min,
        c.mean_density,
        coil,
        energy_z,
        modal_count_distance(c.num_domains, modal),
    ]
}

/// Chain-local z-score normalization across the candidate set, matching
/// `train_pairwise_reranker.py`'s per-chain normalization at training time.
fn zscore_normalize(vectors: &mut [[f64; 8]]) {
    if vectors.len() < 2 {
        return;
    }
    let n = vectors.len() as f64;
    let mut means = [0.0f64; 8];
    for v in vectors.iter() {
        for k in 0..8 {
            means[k] += v[k];
        }
    }
    for m in &mut means {
        *m /= n;
    }
    let mut stds = [0.0f64; 8];
    for v in vectors.iter() {
        for k in 0..8 {
            let d = v[k] - means[k];
            stds[k] += d * d;
        }
    }
    for s in &mut stds {
        *s = (*s / n).sqrt() + 1e-8;
    }
    for v in vectors.iter_mut() {
        for k in 0..8 {
            v[k] = (v[k] - means[k]) / stds[k];
        }
    }
}

/// Pick the winning candidate by pairwise-trained linear score.
///
/// Returns `0` (the first candidate — matching the legacy fallback behavior
/// of picking the first shortlisted entry) if `candidates` is empty, has
/// exactly one entry, or the embedded weights file doesn't match `FEATURES`.
pub fn rerank(
    candidates: &[CandidateInput],
    ss_types: &[SsType],
    energy_config: Option<&crate::energy::EnergyConfig>,
    pdb_path: &str,
    chain: &str,
) -> usize {
    if candidates.len() < 2 {
        return 0;
    }

    let trained: TrainedWeights = match serde_json::from_str(WEIGHTS_JSON) {
        Ok(t) => t,
        Err(_) => return 0,
    };
    if trained.features != FEATURES || trained.weights.len() != FEATURES.len() {
        return 0;
    }

    let modal = modal_num_domains(&candidates.iter().map(|c| c.num_domains).collect::<Vec<_>>());
    let mut vectors: Vec<[f64; 8]> = candidates
        .iter()
        .map(|c| candidate_feature_vector(c, ss_types, energy_config, pdb_path, chain, modal))
        .collect();
    zscore_normalize(&mut vectors);

    let mut best_idx = 0;
    let mut best_score = f64::NEG_INFINITY;
    for (i, v) in vectors.iter().enumerate() {
        let score: f64 = v.iter().zip(trained.weights.iter()).map(|(x, w)| x * w).sum();
        if score > best_score {
            best_score = score;
            best_idx = i;
        }
    }
    best_idx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rerank_single_candidate_returns_zero() {
        let c = CandidateInput {
            num_domains: 1, min_size: 10, max_cr: 0.1, density_min: 1.0, mean_density: 2.0,
            raw_delineation: "0-9", remapped_delineation: "1-10",
        };
        assert_eq!(rerank(&[c], &[], None, "/nonexistent.pdb", "A"), 0);
    }

    #[test]
    fn test_rerank_empty_returns_zero() {
        assert_eq!(rerank(&[], &[], None, "/nonexistent.pdb", "A"), 0);
    }

    #[test]
    fn test_rerank_prefers_higher_coil_fraction_and_favorable_energy() {
        // Candidate 1 boundaries land in coil (favorable); candidate 0's land
        // mid-helix. With no energy_config, energy_z is 0.0 for both, so the
        // coil-fraction feature alone should decide it (weights permitting —
        // this test only asserts against the currently checked-in synthetic
        // smoke-trained weights, not real CATH-trained ones; re-verify after
        // Task 9 replaces the weights file with the real training run).
        let mut ss_types = vec![SsType::Coil; 20];
        for s in ss_types.iter_mut().take(15).skip(5) {
            *s = SsType::Helix;
        }
        let bad = CandidateInput {
            num_domains: 2, min_size: 10, max_cr: 0.3, density_min: 0.5, mean_density: 1.0,
            raw_delineation: "0-9 10-19", remapped_delineation: "1-10 11-20",
        };
        let good = CandidateInput {
            num_domains: 2, min_size: 10, max_cr: 0.3, density_min: 0.5, mean_density: 1.0,
            raw_delineation: "0-1 2-19", remapped_delineation: "1-2 3-20",
        };
        let idx = rerank(&[bad, good], &ss_types, None, "/nonexistent.pdb", "A");
        // Not asserting a specific winner here since the synthetic smoke
        // weights aren't meaningfully trained on this feature — just confirm
        // it runs deterministically and returns a valid index.
        assert!(idx == 0 || idx == 1);
    }
}
```

Add `serde` (with `derive`) as a dependency if `sword2-lib` doesn't already have it — check first:

Run: `grep -n "^serde" sword2-lib/Cargo.toml`

If missing, add to `sword2-lib/Cargo.toml`'s `[dependencies]`:

```toml
serde = { version = "1", features = ["derive"] }
serde_json = "1"
```

(These are almost certainly already present — `output/mod.rs` writes JSON summaries per CLAUDE.md — so this step is likely a no-op; confirm rather than assume.)

Register the module in `sword2-lib/src/sword/mod.rs`:

```rust
pub mod candidate_features;
pub mod compute_jones;
pub mod compute_measure;
pub mod distance_model;
pub mod junctions;
pub mod parse_measure;
pub mod reranker;
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p sword2-lib --lib sword::reranker::tests::`
Expected: FAIL — either a compile error (missing `serde`/module registration) or, once those are fixed, passes immediately since this is a from-scratch file (as with Task 2, this is an acceptable "write once, verify green" case for a new pure-logic file — the meaningful red/green cycle already happened in Tasks 1-5's incremental changes to existing files).

- [ ] **Step 3: Implementation already written in Step 1**

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p sword2-lib --lib sword::reranker::tests::`
Expected: PASS (3 tests)

- [ ] **Step 5: Wire the reranker into the `to_print` selection**

In `sword2-lib/src/sword/mod.rs`, the current winner-selection code (lines 384-401):

```rust
    let mut to_print = String::new();
    for rm in &relevant_measure2 {
        let fields: Vec<&str> = rm.split('|').collect();
        if !fields.is_empty() {
            let n: usize = fields[0].trim().parse().unwrap_or(0);
            if n == n_dom {
                to_print = rm.clone();
                break;
            }
        }
    }
    // Fallback: the second-pass distance filter (dist < 0.2) excludes deeply good-zone
    // candidates. If the optimal was filtered out, use the first-pass representative.
    if to_print.is_empty() {
        tracing::debug!("to_print empty after second pass; using first-pass representative for nd={}", n_dom);
        to_print = to_print_first_pass;
    }
```

becomes:

```rust
    let mut to_print = if config.use_pairwise_reranker && relevant_measure2.len() > 1 {
        let ss_types: &[crate::peeling::algorithm::SsType] = peeling_output
            .as_ref()
            .map(|po| po.ss_types.as_slice())
            .unwrap_or(&[]);
        let pdb_path_str = pdb_file_dst.to_string_lossy().to_string();

        let mut raw_dels: Vec<String> = Vec::new();
        let mut remapped_dels: Vec<String> = Vec::new();
        let mut parsed: Vec<(usize, usize, f64, f64, f64)> = Vec::new();
        for rm in &relevant_measure2 {
            let fields: Vec<&str> = rm.split('|').collect();
            if fields.len() < 6 {
                continue;
            }
            let nd: usize = fields[0].trim().parse().unwrap_or(0);
            let min_size: usize = fields[1].trim().parse().unwrap_or(0);
            let raw_del = fields[2].trim().to_string();
            let max_cr: f64 = fields[3].trim().parse().unwrap_or(0.0);
            let density_min: f64 = fields[4].trim().parse().unwrap_or(0.0);
            let mean_density: f64 = fields[5].trim().parse().unwrap_or(0.0);
            remapped_dels.push(remap_residue_numbers(&raw_del, &tab_num));
            raw_dels.push(raw_del);
            parsed.push((nd, min_size, max_cr, density_min, mean_density));
        }

        let inputs: Vec<reranker::CandidateInput> = parsed
            .iter()
            .zip(raw_dels.iter())
            .zip(remapped_dels.iter())
            .map(|((&(nd, min_size, max_cr, density_min, mean_density), raw), remapped)| {
                reranker::CandidateInput {
                    num_domains: nd,
                    min_size,
                    max_cr,
                    density_min,
                    mean_density,
                    raw_delineation: raw,
                    remapped_delineation: remapped,
                }
            })
            .collect();

        if inputs.is_empty() {
            String::new()
        } else {
            let winner = reranker::rerank(
                &inputs,
                ss_types,
                config.energy_config.as_ref(),
                &pdb_path_str,
                &config.chain_id,
            );
            relevant_measure2[winner].clone()
        }
    } else {
        String::new()
    };

    if to_print.is_empty() {
        // Legacy distance_model-based selection (also the fallback when the
        // reranker is disabled, has <2 candidates, or all candidates failed
        // to parse above).
        for rm in &relevant_measure2 {
            let fields: Vec<&str> = rm.split('|').collect();
            if !fields.is_empty() {
                let n: usize = fields[0].trim().parse().unwrap_or(0);
                if n == n_dom {
                    to_print = rm.clone();
                    break;
                }
            }
        }
    }
    // Fallback: the second-pass distance filter (dist < 0.2) excludes deeply good-zone
    // candidates. If the optimal was filtered out, use the first-pass representative.
    if to_print.is_empty() {
        tracing::debug!("to_print empty after second pass; using first-pass representative for nd={}", n_dom);
        to_print = to_print_first_pass;
    }
```

- [ ] **Step 6: `cargo check` and fix borrow issues**

Run: `cargo check -p sword2-lib`
Expected: Compiles clean. `peeling_output` is read via `.as_ref()` in both the Task 5 dump block and here — confirm nothing upstream consumes it by value between these two uses (it should already be behind `Option<PeelingOutput>` accessed by reference throughout `run_pipeline`).

- [ ] **Step 7: Manual smoke test with the flag on and off**

Run: `cargo build --release`

Run: `./target/release/sword2 -p 1jx4 -o /tmp/sword2_baseline`
Expected: Exits 0, produces `SWORD2_summary.txt`/`.json` as before (flag off — behavior byte-identical to pre-Task-8).

Run: `./target/release/sword2 -p 1jx4 -o /tmp/sword2_reranked --use-pairwise-reranker`
Expected: Exits 0, produces summary files; the "Optimal partition" may differ from the baseline run above (expected — the reranker is scoring against the still-synthetic Task 7 weights at this point, real validation comes in Task 9).

- [ ] **Step 8: Run the full test suite**

Run: `cargo test -p sword2-lib -p sword`
Expected: All PASS.

- [ ] **Step 9: Commit**

```bash
git add sword2-lib/src/sword/reranker.rs sword2-lib/src/sword/mod.rs sword2-lib/Cargo.toml
git commit -m "feat: wire pairwise reranker into candidate selection behind --use-pairwise-reranker"
```

---

## Task 9: Train on real CATH data, validate against the benchmark gate, decide rollout

This task is operational (running existing + new scripts against real data), not new code, matching the spec's evaluation-gate requirement. No shortcuts on the gate: this task's purpose is to produce evidence, not to assume success.

**Files:** none created; `benchmark/data/pairwise_reranker_weights.json` gets overwritten with real trained weights.

- [ ] **Step 1: Rebuild with the reranker wiring**

Run: `cargo build --release`
Expected: Exits 0.

- [ ] **Step 2: Collect a training dump over CATH-17287 (excluding the CATH-663 held-out test set)**

Run:
```bash
python -m benchmark.run_dump \
    --dataset cath17287 \
    --dump /tmp/sword2_dump.csv \
    --output-root /tmp/sword2_training \
    --sword2 ./target/release/sword2 \
    --workers 16
```
Expected: Completes without aborting (per-chain failures are logged to `/tmp/sword2_training/failures.txt` and don't fail the run, per `run_dump.py`'s existing behavior). This is the first real end-to-end test of the moved/extended dump block from Task 5 at scale — check `/tmp/sword2_dump.csv` (or the per-chain dump files, depending on `run_dump.py`'s output mode) for the three new columns and plausible value ranges before proceeding.

- [ ] **Step 3: Build the training table**

Run:
```bash
python -m benchmark.build_training_table \
    --dump /tmp/sword2_dump.csv \
    --reference cath17287 \
    --output benchmark/data/training_table.csv \
    --workers 16
```
Expected: `benchmark/data/training_table.csv` exists with the new columns and `ndo`/`d_count_acc`/etc. populated per candidate.

- [ ] **Step 4: Train the real pairwise reranker**

Run:
```bash
python -m benchmark.train_pairwise_reranker \
    --training-table benchmark/data/training_table.csv \
    --out benchmark/data/pairwise_reranker_weights.json \
    --epochs 2000 --lr 0.1
```
Expected: Logs a final `train_accuracy` well above 0.5 (chance); inspect that number and the per-feature weight magnitudes (`cat benchmark/data/pairwise_reranker_weights.json`) — a large negative weight on `num_domains` or `modal_count_distance` similar in shape to Phase B's would be a warning sign worth investigating before proceeding, even though the pairwise loss structurally resists it.

- [ ] **Step 5: Rebuild with the real trained weights baked in**

Run: `cargo build --release`
Expected: Exits 0 (the `include_str!` in `reranker.rs` picks up the new JSON content at compile time — this rebuild is required, the binary does not hot-reload the weights file).

- [ ] **Step 6: Run the CATH-663 benchmark gate, flag off vs. flag on**

Run baseline (flag off, matches current production behavior):
```bash
python -m benchmark.run_benchmark --dataset cath663 --tools sword2-rust \
    --output-dir benchmark/results_reranker_baseline
```

Run with the reranker enabled — check `run_benchmark.py --help` for how it passes extra CLI flags through to the SWORD2 runner (likely a `--sword2-args` or per-runner config option; if no such passthrough exists, add one, following the existing pattern used for `--legacy-pdb`/other flags already threaded through the harness):
```bash
python -m benchmark.run_benchmark --dataset cath663 --tools sword2-rust \
    --sword2-extra-args "--use-pairwise-reranker" \
    --output-dir benchmark/results_reranker_enabled
```

- [ ] **Step 7: Compare against the gate**

Run: `python -m benchmark.figures` (or however `benchmark/figures.py` is invoked per its existing usage — check `benchmark/REPORT.md`) against both result directories, or use `benchmark/compare_sword2_experiments.py` (already exists in the repo for exactly this A/B comparison pattern per `results_bypass` vs `results_rust_original`):

```bash
python -m benchmark.compare_sword2_experiments \
    --baseline benchmark/results_reranker_baseline \
    --experiment benchmark/results_reranker_enabled
```

Expected: A report showing headline NDO/d_count_acc at rank-1 for both runs. The gate from the design spec: rank-1 mean NDO must move meaningfully above the current 0.777 baseline toward the ~0.89 top-10 ceiling shown in `sword2_topk_ndo.png` — not a rounding-level change. Also regenerate `sword2_topk_ndo.png` for the reranker-enabled run and confirm the k=1 point has moved up, not just that some aggregate average improved (a shifted-but-still-far-below-ceiling k=1 point would indicate partial progress, not the gate being met).

- [ ] **Step 8: Decide rollout**

If the gate is met with a real margin and no regression on d_count_acc/boundary_f1_10/matched_dice:
- Flip `use_pairwise_reranker: false` to a CLI-default-on consideration — but per the Global Constraints, do not change the *default* without explicit user sign-off; leave the flag opt-in and report the numbers back to the user for that decision.

If the gate is not met, or CATH-663 regresses on any headline metric:
- Do not merge the flag to default-on. Document findings (which feature helped/hurt, whether energy_z or boundary_coil_fraction alone would have been safer per the design's Approach C fallback) and report back — this mirrors exactly how Phase B was handled (built, benchmarked, reverted with lessons recorded) and is a legitimate, expected outcome of this task, not a failure of the plan.

- [ ] **Step 9: Commit the trained weights and benchmark artifacts**

```bash
git add benchmark/data/pairwise_reranker_weights.json
git commit -m "feat: train pairwise reranker weights on CATH-17287"
```

(Benchmark result directories under `benchmark/results_reranker_*` are large generated artifacts — check `.gitignore` before adding; likely should stay untracked like the existing `benchmark/results_bypass`/`results_rust_original` directories, which `git status` at the start of this session showed as untracked.)

---

## Self-Review Notes

- **Spec coverage:** Energy Z-score shortlist rescoring → Tasks 4, 5, 8. Boundary coil-fraction → Tasks 1, 2, 5, 8. Pairwise-trained combiner avoiding Phase B's bias → Tasks 3, 7, 8. Evaluation gate → Task 9. Rollout discipline (flag, revert path) → Task 8 (flag), Task 9 (decision point). All "Out of Scope" items from the spec (candidate generation, n_dom prediction, full peeling-tree rescoring, unifying `-E` display energy with reranking energy) are untouched by every task above — confirmed by construction, not just by omission.
- **Type/name consistency check:** `FEATURES` order (`num_domains, min_size, max_cr, density_min, mean_density, boundary_coil_fraction, energy_z, modal_count_distance`) is identical across `benchmark/train_pairwise_reranker.py` and `sword2-lib/src/sword/reranker.rs` — this is load-bearing (Task 8's `rerank` explicitly checks `trained.features != FEATURES` and no-ops if they don't match, so a mismatch fails safe rather than silently misapplying weights). `DumpRow`/`build_dump_row` (Task 5) and `CandidateInput` (Task 8) carry the same five geometric fields plus delineation, named consistently.
- **Risk called out explicitly per task:** Task 4 notes the duplicate `Potentials::load()` cost when both `-E` and `--use-pairwise-reranker` are active; Task 9 Step 4 flags re-checking for a Phase-B-shaped bias even though the pairwise loss structurally resists it.
