# Known Intentional Divergences from the Original SWORD2

This file documents places where the Rust port deliberately produces output that
differs from the original Perl/C SWORD2 pipeline.  Ground truth for tests and
regression checks is the **current Rust output**, not the original tool.

---

## 1. Domain Delineation Output (commit d0c92d6)

**What changed:** The ordering of domain segments in the final output was
corrected.  The original Perl code contained a bug that produced inconsistent
orderings under certain conditions (segments not always listed in ascending
residue-number order).

**Why it was fixed:** The inconsistency was a bug in the original code, not
intentional behaviour.  The Rust port produces a canonical, sorted ordering.

**Impact:** If you diff the Rust output against the original SWORD2 output on
the same structure, the domain boundary strings may appear in a different order.
The actual domain assignments are identical — only the textual representation
differs.

**How to apply:** Do not "fix" the Rust output to match the original.  The
original was wrong.

---

## 2. DSSP Secondary Structure Assignment

The Rust DSSP implementation is an independent port of the Kabsch & Sander
(1983) algorithm.  It produces output that closely matches the original
`dsspcmbi` binary but may differ in edge cases due to:

- Floating-point rounding (Rust `f64` vs. original C `float`/`double` handling).
- The spatial grid H-bond optimization (O(N·k) instead of O(N²)) — both
  strategies search the same neighbour set within the 9 Å CA distance cutoff,
  but floating-point evaluation order differs.

These differences are not considered bugs.  The golden test in
`sword2-lib/tests/dssp_golden.rs` locks the **Rust output** as ground truth.

---

## 3. Peeling Parallelism (rayon)

The Rust Peeling algorithm uses `rayon` to evaluate double-cut candidates in
parallel.  When two cuts score identically, the winner depends on thread
scheduling and may differ from the original serial C implementation.  Results
are deterministic across runs on the same machine (rayon's work-stealing is
deterministic for a fixed thread count and identical inputs), but may differ
from the original tool on ties.

---

## 4. Pseudo-energy / Z-score Scoring

`sword2-lib/src/energy/score.rs` computes pseudo-energy and Z-score in Rust
from the precomputed mypmfs potentials (CA representation, linear interpolation).
Two implementation choices intentionally stabilize edge cases and repeatability:

- **Interpolation boundary.** For an interatomic distance in
  `[last_bin, distmax)` (i.e. `[14.95, 15)` Å with the shipped potentials), the
  scorer clamps to the last valid bin interval `[len-2, len-1]`. This affects
  only the rare pair whose distance falls in that 0.05 Å window.

- **Z-score determinism.** Decoy shuffles use a fixed seed, giving reproducible
  Z-scores across runs.

Ground truth for the golden/regression checks is the current Rust output.
