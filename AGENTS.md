This file provides guidance to Agents when working with code in this repository.

## Project Context

This project primarily involves Rust development (with some Python). When porting Python to Rust, prefer using existing Rust ecosystem crates rather than direct 1:1 translation of Python logic.

## Project Overview

SWORD2 (SWift and Optimized Recognition of protein Domains) is a protein domain partitioning tool implemented as a Rust workspace. It produces multiple alternative domain assignments for a given protein structure via hierarchical clustering of Protein Units (PUs).

## Build Commands

```bash
# Build the Rust binary (release)
cargo build --release
# Binary output: target/release/sword2
# No auxiliary native build step is required.

# Run on a PDB id (from repo root)
./target/release/sword2 -p 1jx4 -o results

# Run on a local file
./target/release/sword2 -i structure.pdb -o results

# Fast mode (skip energy calculation, ~2x faster)
./target/release/sword2 -p 1jx4 -o results --disable-energies

# Run tests
cargo test
```

## Architecture

### Workspace Structure

- **`sword2-cli/`** — Binary crate. CLI entry point (`main.rs`) using clap. Orchestrates the pipeline: input resolution → PDB parsing → chain cleaning → DSSP (pure Rust) → SWORD pipeline → energy calculation → output writing → peeling → plots → junctions → cleanup.
- **`sword2-lib/`** — Library crate. All core logic, organized as modules:
  - `sword/` — Pipeline orchestration (`mod.rs`), PU merging (`compute_measure.rs`), domain selection (`parse_measure.rs`), distance model (`distance_model.rs`), Jones metrics (`compute_jones.rs`), junction analysis (`junctions.rs`)
  - `pdb/` — PDB/mmCIF parsing (`parser.rs`) via `pdbtbx`, type definitions (`types.rs`), PDB writing (`writer.rs`), amino acid definitions (`amino_acids.rs`)
  - `energy/` — **Pure Rust** pseudo-energy and Z-score calculation (`score.rs`) using the precomputed mypmfs potentials (CA representation, linear interpolation). `mod.rs` holds `EnergyConfig`/`EnergyResult` and the parallel batch helpers; potentials are loaded once and cached
  - `peeling/` — **Pure Rust Peeling implementation** (Gelly et al. 2006). Modules: `algorithm.rs` (iterative hierarchical cutting, rayon-parallelized double cuts), `contact_matrix.rs` (contact probability matrix with 2D prefix sums), `mod.rs` (types, result conversion, legacy file parsing)
  - `dssp/` — **Pure Rust DSSP implementation** (Kabsch & Sander 1983 algorithm). Modules: `backbone.rs` (atom extraction, H synthesis), `hbond.rs` (spatial grid H-bond detection), `bridge.rs` (β-sheet assembly), `helix.rs` (helix/turn assignment), `angles.rs` (backbone geometry), `format.rs` (DSSP output format), `types.rs` (data structures)
  - `fetch.rs` — Downloads structures from PDB, AlphaFold, ESM Atlas
  - `output/` — Writes SWORD2_summary.txt/json results
  - `plot/` — SVG plot generation via `plotters`

### Self-contained Rust pipeline

The pipeline implements DSSP, Peeling, and scoring in Rust:
- DSSP: spatial grid optimization for H-bond detection (O(N·k) vs original O(N²))
- Peeling: rayon-parallelized double cutting with 2D prefix-sum contact matrix (O(1) rectangle queries)
- Scoring: `energy/score.rs` (pseudo-energy + rayon-parallelized Z-score decoys)

The only thing kept under `bin/` is the precomputed potential data
(`bin/mypmfs-master/025_30_100_potential/`), read directly by `score.rs`.

### Key Types

- `SwordConfig` / `SwordResults` — Pipeline configuration and results (`sword2-lib/src/sword/mod.rs`)
- `EnergyConfig` / `EnergyResult` — Energy calculation config and results (`sword2-lib/src/energy/mod.rs`)
- `PeelingConfig` / `PeelingOutput` — Peeling algorithm config and raw output (`sword2-lib/src/peeling/algorithm.rs`)
- `ContactMatrix` — Contact probability matrix with prefix sums (`sword2-lib/src/peeling/contact_matrix.rs`)
- `PeelingLevel` / `ProteinUnit` — Peeling decomposition results (`sword2-lib/src/peeling/mod.rs`)
- `Partitioning` — Domain assignment output (`sword2-lib/src/output/mod.rs`)
- PDB types: `Structure`, `Model`, `Chain`, `Residue`, `Atom` (`sword2-lib/src/pdb/types.rs`)
- DSSP types: `DsspChain`, `BackboneResidue`, `HydrogenBond`, `Bridge` (`sword2-lib/src/dssp/types.rs`)

### Pipeline Flow

1. Fetch/load structure → 2. Parse PDB/mmCIF → 3. Clean chain (remove non-standard residues, renumber from 1) → 4. Run DSSP (pure Rust) → 5. Run Peeling (pure Rust, in-memory) → 6. ComputeMeasure (merge PUs, in-memory from peeling) → 7. ParseMeasure + distance model (select domains, pure Rust) → 8. Calculate pseudo-energies (pure Rust, `energy/score.rs`) → 9. Write results (JSON + text) → 10. Junction consistency analysis → 11. Cleanup

### Important Notes

- The binary auto-detects `bin/` directory relative to its location; use `--base-dir` when running from a non-standard location.
- Parallelism uses rayon throughout (DSSP, Peeling, energy decoys); controlled via `-x`/`--cpu` flag.
- Logging via `tracing`; default level is `info` for the `sword2` target. Control with `RUST_LOG` env var.

## Workflow Guidelines

When working on long tasks (porting, migration, research), break work into smaller committed checkpoints so progress isn't lost to rate limits.

## Testing

Always run `cargo check` and `cargo test` after modifying Rust files before considering a task complete.
