# Deterministic Parallel Factorized Training Design

## Objective

Reduce the wall time of the factorized structural-ranker training run by
executing independent model fits in eight worker processes, while preserving
the existing scientific grid, folds, seed, feature-family order, checkpoint
integrity, output bytes, and provenance guarantees.

The machine exposes 40 logical CPUs (20 physical cores) and 64 GB of physical
RAM. Each sklearn estimator remains single-threaded; concurrency occurs only
across independent `(head, hyperparameters, fold)` or OOF-fold jobs.

## Scope

Modify only the factorized training implementation, its CLI, and focused
tests. The change adds an explicit `--jobs` argument with accepted integer
values `1..=8`.

- `--jobs 1` preserves the existing sequential execution path and results.
- `--jobs 8` is the production Task 13 training command.
- The model grid, objectives, ties, family gates, seed 37, fold assignments,
  input dtype, estimator parameters, and artifact formats do not change.
- No CATH-663 prediction, metric, or competitor output is read or executed.
- No new dependency is added.

## Architecture

### Parent scheduler

The parent process remains the sole authority for task enumeration,
checkpoint reads and writes, result validation, aggregation, stage decisions,
reports, and final artifacts.

For each head and feature-family proposal, it constructs the canonical 40-job
sequence in `MODEL_GRID` order followed by fold order. It first loads and
validates every available fold checkpoint. Missing jobs are submitted to one
fixed-size process pool, up to `jobs` at a time. Returned results are indexed
by their exact task identity and are aggregated only in canonical grid/fold
order, never completion order.

### Worker processes

On Linux, workers use the explicit `fork` multiprocessing context. They
inherit the verified, read-only corpus and folds through copy-on-write rather
than receiving serialized copies of the full data frames. Each submitted grid
task contains only the head, frozen hyperparameters, retained-family tuple,
fold, and seed.

A worker constructs the training subset and pair batch, fits one
`GradientBoostingClassifier`, evaluates the held-out fold, and returns a typed
result containing its identity, aggregate values, and per-chain values. It
does not write a checkpoint, mutate stage state, select hyperparameters, or
inspect another task.

The seven deterministic environment variables remain fixed at one thread per
process:

```text
PYTHONHASHSEED=0
OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
NUMEXPR_NUM_THREADS=1
RAYON_NUM_THREADS=1
```

### OOF jobs

The five folds used for end-to-end OOF generation are independent. They use
the same worker infrastructure and return one typed fold payload. The parent
validates exact fold populations, combines rows and count decisions in fold
and chain-ID order, and produces the existing canonical OOF CSV bytes.

When an ablation reuses count decisions, that exact validated mapping is an
input to each OOF task. It cannot be recomputed or changed by a worker.

### Sequential fallback

`jobs == 1` calls the existing in-process fold logic. It does not construct a
process pool. This is both the compatibility path and the authoritative
equivalence reference in tests.

## Checkpoint and failure contract

The parent loads checkpoints before submission and writes one checkpoint only
after validating a returned result against the expected stage, family set,
head, parameters, fold, train IDs, validation IDs, feature names, and seed.
Checkpoint payloads remain canonical JSON and use the existing atomic writer.

Only the parent writes, so concurrent temporary-file or rename races are
impossible. Completion order cannot alter filenames or bytes.

If any worker raises, exits abnormally, returns a duplicate identity, or
returns a result outside the submitted set, the parent cancels pending work,
shuts down the pool, and fails the stage. Successfully validated checkpoints
already installed remain resumable. No hyperparameter selection, OOF result,
stage decision, report, or final model is written from a partial stage.

`jobs > 1` fails before training on a platform without the `fork` start
method. Worker count must be a non-boolean integer in `1..=8`.

## Determinism

Every estimator keeps its existing independent `random_state=37`. No shared
random-number generator is introduced. Worker completion order is ignored;
all reductions and serialized output follow the existing canonical order.

The parallel implementation is acceptable only if focused fixtures prove:

- sequential and parallel grid selections are exactly equal;
- fold checkpoint payload bytes are exactly equal;
- sequential and parallel OOF CSV bytes are exactly equal;
- repeated parallel runs are byte-identical;
- reversing task completion order does not change aggregation;
- a partial checkpoint population causes only missing jobs to run.

## CLI and provenance

The normalized training command includes `--jobs 8`. The checkpoint context
continues to bind the exact command, source Git commit, training-source hash,
dataset and table hashes, corpus manifest, folds, schema, grid, seed, IDs, and
package versions. Therefore sequential checkpoints are intentionally
incompatible with the parallel run.

After the implementation commit, rebuild the normalized corpus manifest and
fold manifest twice from the unchanged raw acquisition. Require the four
scientific tables and fold assignments to remain byte-identical; only
source/provenance-bound bytes may change. The new training run begins only
after these comparisons pass.

## Resource policy

Launch the production run with eight workers. During its first full eight-job
wave, measure the parent and descendant RSS without inspecting model metrics.
Require total resident memory to remain below 75% of physical RAM
(50,361,477,120 bytes for this machine). A breach is an infrastructure failure:
stop the run before accepting it as the production execution and revise the
explicit worker count. Do not adapt concurrency in response to scores,
selected parameters, or feature-family outcomes.

## Cleanup and restart

Once the implementation plan is approved and before the first source edit:

1. Wait for or identify the current atomic-checkpoint boundary.
2. Stop the old sequential Python process.
3. Confirm the PID and exact command before deletion.
4. Delete `benchmark/data/factorized_ranker_v1_training/`, which contains only
   the now-incompatible sequential reports/checkpoints.
5. Delete the zero-byte
   `benchmark/data/factorized_ranker_v1_run_ledger/training.time.txt`.

Do not archive those obsolete outputs. Do not delete the raw acquisition dump,
parts, normalized scientific tables, chain cache, acquisition provenance,
opaque CATH-663 baseline, or other audit ledgers: they remain required inputs
or evidence.

The replacement run uses a fresh canonical training directory and redirects
stdout, stderr, progress, and GNU-time output to durable ledger files so its
status remains inspectable after a chat/session disconnect.

## Verification

Use test-driven development. Before production code changes, add focused tests
for worker-count validation, sequential/parallel equivalence, checkpoint
reuse, completion-order invariance, OOF equivalence, and worker failure.
Observe the required RED failures, implement the smallest scheduler change,
then run:

```text
focused factorized-training tests
full benchmark Python tests
Python bytecode compilation for changed modules
git diff --check on the exact changed paths
```

After committing the implementation, rebuild and independently reproduce the
corpus/folds, launch the eight-worker run, enforce the memory gate, and verify
that the first wave writes canonical fold checkpoints with unique identities.

## Success criteria

- Eight workers remain simultaneously compute-bound during a full grid wave.
- Peak total RSS stays below the fixed 75% limit.
- Sequential and parallel fixtures are byte-identical.
- The production checkpoint context binds `--jobs 8` and the new source commit.
- No obsolete sequential checkpoint or empty timing artifact remains.
- Scientific data, folds, model choices, and downstream Task 13–18 contracts
  are otherwise unchanged.
