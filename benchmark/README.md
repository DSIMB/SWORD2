# scoring backend benchmark

Times the mypmfs scoring backend on proteins of different sizes and reports
wall-clock (median / min / stddev) per structure. Compares the **C++** scorer
(`bin/mypmfs-master/scoring_omp`) against the **pure-Rust** scorer
(`sword2 score`, backed by `sword2-lib/src/energy/score.rs`).

> **Status:** both engines are live. Energy matches the C++ binary essentially
> exactly; the Z-score differs slightly (fixed seed vs C++ `time(NULL)`). On a
> single whole-structure call the Rust scorer is faster across all sizes (raw
> ~2×; full Z-score ~1.5–2×). The Rust win is larger still inside the real
> pipeline, which calls the scorer many times: Rust loads the 210 potential
> files once and caches them, whereas the C++ binary reloads them on every fork —
> a cost this single-call benchmark does **not** measure.

## Prerequisites

```bash
cargo build --release          # builds target/release/sword2 (cleaning + rust engine)
make -C bin/mypmfs-master      # builds scoring_omp (only needed for the cpp engine)
```

## Run

```bash
# Quick smoke run (two small structures, 3 repeats, both engines):
python3 benchmark/run_benchmark.py --pdb-ids 1CRN 1UBQ --repeats 3

# Full default size-spanning set (fetches from RCSB on first run):
python3 benchmark/run_benchmark.py

# One engine only:
python3 benchmark/run_benchmark.py --engines rust
```

Outputs a markdown table to stdout and writes `benchmark/results.csv` +
`benchmark/results.md`. Useful flags: `--repeats`, `--shuffles` (Z-score decoys,
default 2000 to match the pipeline), `--cpu` (threads; 0 = auto), `--pdb-ids`,
`--engines` (default `cpp rust`).

## How it works

For each structure the harness:

1. **Cleans** it by running `sword2 -p <id> -o benchmark/structures/<id>
   --disable-energies`. This reuses the pipeline's own cleaning to produce the
   exact renumbered single-chain `input.pdb` the scorer consumes (cached across
   runs). It does **not** run the energy step.
2. Times the scorer directly on that `input.pdb`, matching the real invocation:
   `scoring_omp -i input.pdb -d 025_30_100_potential -z -s <shuffles>`.
3. Measures **two modes** to separate cost centers:
   - `raw` — no `-z`: deterministic distance + linear-interpolation cost.
   - `full` — adds `-z -s <shuffles>`: the decoy-shuffle Z-score cost (dominant
     at 2000 decoys).

Each measurement does 1 warmup + N timed runs. `energy`/`zscore` columns are a
produced-output sanity check; the authoritative C++-vs-Rust equivalence assertion
lives in `sword2-lib/tests/energy_agreement.rs`.

The two engines are defined by the `Engine` abstraction in `run_benchmark.py`:
- `cpp` → `scoring_omp -i <pdb> -d <dir> [-z -s N]`
- `rust` → `sword2 score --pdb <pdb> --potential-dir <dir> [--shuffles N | --no-zscore]`

Compare the `median_s` columns for matching `(pdb_id, mode)` rows.

### Fairness caveat

Time the **same whole-structure single call** the harness uses here. Do **not**
benchmark at the pipeline level: the production C++ path in `energy/mod.rs` forks
one `scoring_omp` process per PU/domain, so a pipeline comparison would mostly
measure process-startup amortization (which favors an in-process Rust scorer),
not the scoring algorithm itself.

## Notes

- Fetched structures live in `benchmark/structures/` and results in
  `benchmark/results.*` — both git-ignored.
- Default ids span sizes: `1CRN` (~46 res) · `1UBQ` (~76) · `1AKE` (~214) ·
  `1JX4` (335, the test structure) · `1QCF` (~449). Sizes are reported from the
  cleaned PDB, not hard-coded; override with `--pdb-ids`.
