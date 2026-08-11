# Deterministic Parallel Factorized Training Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce the fixed factorized-ranker training wall time by evaluating independent grid folds and OOF folds in eight forked Python worker processes without changing any scientific input, decision, checkpoint byte, or final artifact byte.

**Architecture:** The parent process validates inputs, enumerates tasks in canonical model-grid/fold order, reads and writes every checkpoint, validates worker identities, and aggregates results in canonical order. Linux `fork` workers inherit one verified read-only corpus through copy-on-write, receive only small frozen task identities, fit one fold at a time, and return immutable payloads; `jobs=1` remains an in-process reference path that uses the same fold functions.

**Tech Stack:** Python 3.12, `concurrent.futures.ProcessPoolExecutor`, `multiprocessing` with explicit `fork`, pandas/numpy, scikit-learn `GradientBoostingClassifier`, pytest, canonical JSON/CSV checkpoint formats, GNU `time`, Linux `/proc`/`ps` resource inspection.

---

## Global constraints

- Work in the current checkout. Preserve unrelated user changes, especially `benchmark/data/geometry_reference.json` and untracked files.
- Do not dispatch sub-agents; the user selected inline execution.
- Do not run, read, score, or inspect CATH-663 or competitor outputs before the locked Task 17 runtime freeze.
- Keep `MODEL_GRID`, five folds, feature-family order, seed 37, objectives, tie order, pair construction, estimators, reports, and model formats unchanged.
- Accept only non-boolean integers in `1..=8`; use `jobs=8` for production and `jobs=1` as the exact reference.
- Set `PYTHONHASHSEED=0`, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1`, and `RAYON_NUM_THREADS=1` for parity tests and production.
- The parent is the only checkpoint/stage/report/model writer. A worker may not write files or select a model.
- A parallel worker failure cancels pending tasks and prevents selection/stage/output creation from partial results. Canonical fold checkpoints already validated and installed by the parent remain resumable.
- The first production eight-worker wave must keep total parent-and-descendant RSS below 50,361,477,120 bytes. This resource gate must not examine scores or adapt based on metrics.

## File map

- Modify `benchmark/factorized_ranker/training.py`: frozen task/result types, jobs validation, fork worker contexts, parallel grid scheduler, parallel OOF scheduler, parent-side validation, and progress messages.
- Modify `benchmark/train_factorized_ranker.py`: `--jobs`, propagation, help text, and normalized-command provenance.
- Modify `benchmark/tests/test_factorized_training.py`: RED/GREEN coverage for validation, exact sequential/parallel equality, checkpoint bytes/reuse, canonical completion order, OOF bytes, and worker failure.
- Modify `.superpowers/sdd/2026-08-07-factorized-structural-ranker/task-13-brief.md`: freeze the production command as `--jobs 8` and name the durable logs/memory gate.
- Generate no new Python module and add no dependency.
- Regenerate source-bound files under `benchmark/data/cath17287_factorized_corpus_v1/` and `benchmark/data/factorized_ranker_v1_freeze_stage/` only after the final implementation commit.

### Task 1: Retire the incompatible sequential execution

**Files:**
- Delete: `benchmark/data/factorized_ranker_v1_training/`
- Delete: `benchmark/data/factorized_ranker_v1_run_ledger/training.time.txt`
- Preserve: raw dump, raw parts, acquisition provenance/rejections, chain cache, normalized scientific inputs until their controlled rebuild, and all unrelated ledgers.

- [ ] **Step 1: Resolve the exact live process and targets without mutation**

Run:

```bash
ps -p 2310978 -o pid=,ppid=,etimes=,pcpu=,rss=,stat=,args=
readlink -f /proc/2310978/cwd
tr '\0' ' ' </proc/2310978/cmdline
test "$(readlink -f /proc/2310978/cwd)" = "$PWD"
test -d benchmark/data/factorized_ranker_v1_training
test "$(stat -c %s benchmark/data/factorized_ranker_v1_run_ledger/training.time.txt)" = 0
```

Expected: PID 2310978 is the exact old `benchmark.train_factorized_ranker ... --seed 37 --resume` command rooted in this repository; the timing file is zero bytes.

- [ ] **Step 2: Stop the old process at its current safe computational boundary**

Run:

```bash
kill -INT 2310978
for attempt in $(seq 1 30); do
  if ! kill -0 2310978 2>/dev/null; then break; fi
  sleep 1
done
if kill -0 2310978 2>/dev/null; then
  kill -TERM 2310978
fi
for attempt in $(seq 1 30); do
  if ! kill -0 2310978 2>/dev/null; then break; fi
  sleep 1
done
! kill -0 2310978 2>/dev/null
```

Expected: the obsolete trainer exits. A process still alive after SIGTERM is an infrastructure blocker; do not use SIGKILL without inspecting its exact state again.

- [ ] **Step 3: Delete only the now-invalid outputs the user authorized**

Run:

```bash
test "$PWD" = /home/chili/cretin/PROJECTS/SWORD2
test -d benchmark/data/factorized_ranker_v1_training
rm -r -- benchmark/data/factorized_ranker_v1_training
rm -- benchmark/data/factorized_ranker_v1_run_ledger/training.time.txt
test ! -e benchmark/data/factorized_ranker_v1_training
test ! -e benchmark/data/factorized_ranker_v1_run_ledger/training.time.txt
```

Expected: only the obsolete checkpoint/report directory and empty timing file are gone. Report that this deletion is intentional and not recoverable except by recomputation.

### Task 2: Add RED tests for the jobs contract and grid scheduler

**Files:**
- Modify: `benchmark/tests/test_factorized_training.py`
- Test: `benchmark/tests/test_factorized_training.py`

- [ ] **Step 1: Add imports and a checkpoint-byte helper**

Add `multiprocessing` and `time` imports, then add:

```python
def _checkpoint_bytes(root: Path) -> dict[str, bytes]:
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in sorted(root.rglob("*.json"))
    }
```

- [ ] **Step 2: Add strict jobs and CLI tests**

Add:

```python
@pytest.mark.parametrize("value", [True, False, 0, -1, 9, 1.0, "8", None])
def test_parallel_jobs_reject_invalid_library_values(value: object) -> None:
    with pytest.raises(ValueError, match="jobs"):
        training.validate_training_jobs(value)  # type: ignore[arg-type]


def test_parallel_jobs_accept_exact_supported_range() -> None:
    assert [training.validate_training_jobs(value) for value in range(1, 9)] == list(
        range(1, 9)
    )


def test_training_cli_exposes_bounded_jobs(capsys: pytest.CaptureFixture[str]) -> None:
    from benchmark.train_factorized_ranker import main

    with pytest.raises(SystemExit) as help_exit:
        main(["--help"])
    assert help_exit.value.code == 0
    assert "--jobs" in capsys.readouterr().out
    with pytest.raises(SystemExit) as invalid_exit:
        main([
            "--corpus-dir", "unused", "--fold-manifest", "unused",
            "--out-dir", "unused", "--jobs", "9",
        ])
    assert invalid_exit.value.code == 2
    assert "1..=8" in capsys.readouterr().err
```

- [ ] **Step 3: Add real sequential/parallel grid and checkpoint equivalence**

Add a test that uses the committed 15-chain fixture and two stores with identical contexts:

```python
@pytest.mark.skipif("fork" not in multiprocessing.get_all_start_methods(), reason="fork required")
def test_parallel_grid_matches_sequential_selection_and_checkpoint_bytes(
    tmp_path: Path,
) -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    context = {"source_git_commit": "a" * 40, "seed": 37}
    sequential_root = tmp_path / "sequential"
    parallel_root = tmp_path / "parallel"
    sequential = select_head_hyperparameters(
        data.corpus, data.assignments, "count", ("base",), jobs=1,
        checkpoint_store=TrainingCheckpointStore(sequential_root, context),
        stage_family="base",
    )
    parallel = select_head_hyperparameters(
        data.corpus, data.assignments, "count", ("base",), jobs=8,
        checkpoint_store=TrainingCheckpointStore(parallel_root, context),
        stage_family="base",
    )
    assert parallel == sequential
    assert _checkpoint_bytes(parallel_root) == _checkpoint_bytes(sequential_root)
```

- [ ] **Step 4: Add reversed-completion and partial-resume tests**

The scheduler will expose `_assemble_grid_evaluations(head, results)` as a pure canonical reducer. Add a test that obtains the 40 `_HeadFoldWorkResult` values from the sequential helper, passes them in canonical and reversed order, and asserts identical `HeadGridEvaluation` tuples. Add a partial-resume test that copies 13 canonical grid checkpoint files plus the manifest from a complete sequential directory into a new store, records those 13 bytes and mtimes, runs `jobs=8`, and requires:

```python
assert resumed == complete
assert len(list((partial_root / "grid").glob("*.json"))) == 40
assert {path.name: path.read_bytes() for path in original_cached} == cached_bytes
assert {path.name: path.stat().st_mtime_ns for path in original_cached} == cached_mtimes
```

- [ ] **Step 5: Add worker-failure isolation test**

At module scope define a picklable failing fit function:

```python
def _always_fail_fit(*_args: object, **_kwargs: object) -> object:
    raise RuntimeError("injected worker failure")
```

Monkeypatch `training.fit_head` to it before the fork pool is created, call parallel selection with a fresh checkpoint store, and require `RuntimeError("injected worker failure")`, no `stage_state.json`, and no selection return. Grid fold checkpoints that happened to complete before the failure are allowed, but every file must pass `TrainingCheckpointStore.load_fold` validation.

- [ ] **Step 6: Run the new focused tests and capture RED evidence**

Run:

```bash
PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
RAYON_NUM_THREADS=1 benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_training.py -q
```

Expected: failures name missing `validate_training_jobs`, missing `jobs` parameters, and missing `_assemble_grid_evaluations`; existing tests remain green.

### Task 3: Implement deterministic fold-level grid execution

**Files:**
- Modify: `benchmark/factorized_ranker/training.py`
- Modify: `benchmark/train_factorized_ranker.py`
- Test: `benchmark/tests/test_factorized_training.py`

- [ ] **Step 1: Add standard-library imports and strict validation**

Add:

```python
import multiprocessing
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
```

Define:

```python
def validate_training_jobs(jobs: int) -> int:
    if isinstance(jobs, bool) or not isinstance(jobs, Integral) or not 1 <= int(jobs) <= 8:
        raise ValueError("training jobs must be a non-boolean integer in 1..=8")
    return int(jobs)


def _fork_context() -> multiprocessing.context.BaseContext:
    if "fork" not in multiprocessing.get_all_start_methods():
        raise ValueError("parallel factorized training requires the fork start method")
    return multiprocessing.get_context("fork")
```

Call `validate_training_jobs` before any pool or fit in every public entry point accepting `jobs`.

- [ ] **Step 2: Define immutable grid identities and payloads**

Add beside `HeadFoldResult`:

```python
@dataclass(frozen=True)
class _HeadFoldTask:
    head: Literal["count", "candidate"]
    params: Hyperparameters
    retained_families: tuple[str, ...]
    fold: int
    seed: int


@dataclass(frozen=True)
class _HeadFoldWorkResult:
    task: _HeadFoldTask
    fold_result: HeadFoldResult
    primary_values: tuple[float, ...]
    secondary_values: tuple[float, ...]


@dataclass(frozen=True)
class _GridWorkerContext:
    corpus: CorpusTables
    assignments: tuple[FoldAssignment, ...]


_GRID_WORKER_CONTEXT: _GridWorkerContext | None = None
```

Task identity contains every scientific worker choice; the corpus and assignments remain inherited global read-only state rather than task arguments.

- [ ] **Step 3: Extract one-fold evaluation without changing formulas**

Move the existing loop body into:

```python
def _evaluate_head_fold(
    corpus: CorpusTables,
    assignments: tuple[FoldAssignment, ...],
    task: _HeadFoldTask,
) -> _HeadFoldWorkResult:
```

It must derive `validation_ids` from `task.fold`, derive `training_ids` as the exact complement, build only the selected head's pair batch using `head_feature_spec`, fit with `task.params` and `task.seed`, iterate validation chain IDs in sorted order, preserve the existing count accuracy/absolute-error and candidate mean-NDO formulas, and return tuples. Reject empty folds, nonfinite aggregates, unexpected secondary values, or a fold outside `range(5)`.

Define the picklable worker entry:

```python
def _run_head_fold_task(task: _HeadFoldTask) -> _HeadFoldWorkResult:
    context = _GRID_WORKER_CONTEXT
    if context is None:
        raise RuntimeError("grid worker context is unavailable")
    return _evaluate_head_fold(context.corpus, context.assignments, task)
```

- [ ] **Step 4: Make canonical result validation and aggregation pure**

Define `_validate_head_fold_work_result(expected_task, expected_validation_ids, result)` to require exact task equality, fold equality, sorted exact validation IDs, exact vector lengths, finite values, exact `float(np.mean(...))` aggregates, and the head-specific secondary contract.

Define:

```python
def _assemble_grid_evaluations(
    head: Literal["count", "candidate"],
    results: Mapping[_HeadFoldTask, _HeadFoldWorkResult],
) -> tuple[HeadGridEvaluation, ...]:
```

For each `params` in `MODEL_GRID` and each `fold` in `range(5)`, look up the exact task, extend per-chain values in fold order, and construct `HeadGridEvaluation`. Reject missing or extra tasks. Compute the overall primary/secondary with the same single `np.mean` over concatenated chain values used before; never average fold means.

- [ ] **Step 5: Implement the parent-only 40-task scheduler**

Refactor `select_head_hyperparameters(..., jobs: int = 1)`:

1. Validate and canonicalize corpus/assignments/families once.
2. Enumerate tasks in `MODEL_GRID` then fold order.
3. Build exact checkpoint kwargs in the parent.
4. Load each checkpoint first and turn cached data into validated `_HeadFoldWorkResult` values.
5. For `jobs == 1`, run missing tasks directly in canonical order and write each validated result from the parent.
6. For `jobs > 1`, set `_GRID_WORKER_CONTEXT`, create `ProcessPoolExecutor(max_workers=jobs, mp_context=_fork_context())`, submit only missing tasks, consume `as_completed`, validate exact submitted identity/no duplicates, write the result's checkpoint from the parent, and clear the global in `finally`.
7. On any exception, cancel every pending future, use `shutdown(wait=True, cancel_futures=True)`, clear the global, and re-raise without calling the reducer.
8. Call `_assemble_grid_evaluations`, then retain the existing tie-selection key verbatim.

Emit a flushed parent-side stderr line after each accepted result without metric values:

```python
print(
    f"factorized-training grid head={head} family={stage_family or families[-1]} "
    f"completed={len(results)}/40 cached={cached_count} fitted={fitted_count}",
    file=sys.stderr,
    flush=True,
)
```

- [ ] **Step 6: Preserve the sequential compatibility helper**

Keep `_evaluate_head_configuration` as a jobs=1 helper over `_evaluate_head_fold` so existing callers/tests retain their API. It must load/write checkpoints in the parent and return the same bytes and dataclass equality as before.

- [ ] **Step 7: Add and propagate CLI jobs**

In `benchmark/train_factorized_ranker.py`, define:

```python
def _jobs_argument(value: str) -> int:
    try:
        jobs = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("jobs must be an integer in 1..=8") from error
    if not 1 <= jobs <= 8:
        raise argparse.ArgumentTypeError("jobs must be an integer in 1..=8")
    return jobs
```

Add `parser.add_argument("--jobs", type=_jobs_argument, default=1, help="independent fold workers (1..=8)")`, pass `jobs=args.jobs` to `run_grouped_training`, and let the existing normalized argv retain `--jobs 8` exactly. Do not treat it as a path role.

- [ ] **Step 8: Run focused GREEN and commit the grid checkpoint**

Run the focused command from Task 2. Expected: all focused tests pass, including exact checkpoint bytes. Then run:

```bash
benchmark/.venv/bin/python -m py_compile \
  benchmark/factorized_ranker/training.py benchmark/train_factorized_ranker.py
git diff --check -- benchmark/factorized_ranker/training.py \
  benchmark/train_factorized_ranker.py benchmark/tests/test_factorized_training.py
git add benchmark/factorized_ranker/training.py benchmark/train_factorized_ranker.py \
  benchmark/tests/test_factorized_training.py
git commit -m "feat: parallelize factorized grid folds"
```

Expected: a scoped commit; unrelated dirty files remain unstaged.

### Task 4: Add RED tests and implementation for parallel OOF folds

**Files:**
- Modify: `benchmark/factorized_ranker/training.py`
- Modify: `benchmark/tests/test_factorized_training.py`

- [ ] **Step 1: Add exact sequential/parallel OOF tests**

Using the committed 15-chain fixture, generate base OOF with `jobs=1` and `jobs=8`, then assert complete dataclass equality, exact `csv_bytes`, exact SHA-256, and exact `count_decisions`. Repeat a fresh `jobs=8` call and require identical bytes. Repeat the comparison with `reused_count_decisions=sequential.count_decisions` to cover candidate-only ablations.

- [ ] **Step 2: Add canonical order and failure tests**

Test a pure `_combine_oof_fold_results` with canonical and reversed fold payloads and require equal `OOFResult`. Monkeypatch the module-scope `fit_head` to `_always_fail_fit`, call `generate_oof(..., jobs=8)`, and require the injected error with no returned/serialized OOF object.

- [ ] **Step 3: Run the OOF tests and capture RED evidence**

Run only the new node IDs with the deterministic environment. Expected: `generate_oof` rejects `jobs`, and `_combine_oof_fold_results` is absent.

- [ ] **Step 4: Define immutable OOF tasks/results and inherited context**

Add:

```python
@dataclass(frozen=True)
class _OOFFoldTask:
    fold: int
    retained_families: tuple[str, ...]
    count_params: Hyperparameters
    candidate_params: Hyperparameters
    seed: int
    reuse_count_decisions: bool


@dataclass(frozen=True)
class _OOFFoldWorkResult:
    task: _OOFFoldTask
    rows: tuple[OOFRow, ...]
    count_decisions: tuple[tuple[str, CountDecision], ...]


@dataclass(frozen=True)
class _OOFWorkerContext:
    data: VerifiedTrainingData
    reused_count_decisions: Mapping[str, CountDecision] | None


_OOF_WORKER_CONTEXT: _OOFWorkerContext | None = None
```

- [ ] **Step 5: Extract one OOF fold and validate its exact population**

Move the current `for fold in range(5)` body unchanged into `_generate_oof_fold(data, task, reused_count_decisions)`. Return rows sorted by chain ID and decisions as sorted `(chain_id, CountDecision)` tuples. Require exact task fold population, one row/decision per validation ID, correct row fold values, no duplicate IDs, and reuse-mode agreement.

Add `_run_oof_fold_task(task)` that reads `_OOF_WORKER_CONTEXT`, validates the boolean reuse identity against whether a mapping exists, and calls the pure fold helper.

- [ ] **Step 6: Combine only in canonical fold/chain order**

Define `_combine_oof_fold_results(data, tasks, results)` to require exactly five distinct expected task identities; append each fold's rows and decisions in fold then sorted-chain order; reject extra/missing/duplicate rows or decisions; apply the existing whole-population checks; call `_render_oof`; and return the unchanged `OOFResult`.

- [ ] **Step 7: Implement jobs=1/jobs>1 OOF dispatch**

Extend `generate_oof(..., jobs: int = 1)`. After existing family/population/reuse checks, enumerate five tasks. For jobs=1, call the fold helper in-process. For jobs>1, set `_OOF_WORKER_CONTEXT`, use an explicit fork executor with `min(jobs, 5)` workers, submit all five tasks, accept only exact unique identities, cancel/shutdown/re-raise on failure, and clear the global in `finally`. Emit no metric values; print only:

```python
factorized-training oof family=<last-family> completed=<n>/5
```

Then combine in canonical order.

- [ ] **Step 8: Propagate jobs through every grouped-training call**

Add `jobs: int = 1` to `run_grouped_training`, validate it before loading stage state, and pass it to every `select_head_hyperparameters` and `generate_oof` call, including mandatory base, each proposal, and the final fresh OOF. Keep the two final full-data fits sequential because they are only two dependent output fits.

- [ ] **Step 9: Run GREEN and commit the OOF checkpoint**

Run the full focused file, py_compile, and diff check. Expected: all pass. Commit only the two modified files:

```bash
git add benchmark/factorized_ranker/training.py benchmark/tests/test_factorized_training.py
git commit -m "feat: parallelize factorized OOF folds"
```

### Task 5: Verify determinism, failure behavior, and repository scope

**Files:**
- Modify: `.superpowers/sdd/2026-08-07-factorized-structural-ranker/task-13-brief.md`
- Verify: all Python changes

- [ ] **Step 1: Freeze the Task 13 production invocation and logs**

Change the training command to include `--jobs 8 --resume`, require the seven one-thread environment variables, and name:

```text
benchmark/data/factorized_ranker_v1_run_ledger/training.stdout.log
benchmark/data/factorized_ranker_v1_run_ledger/training.stderr.log
benchmark/data/factorized_ranker_v1_run_ledger/training.time.txt
benchmark/data/factorized_ranker_v1_run_ledger/training.pid
```

Add the fixed 50,361,477,120-byte first-wave RSS gate and explicitly forbid inspecting partial model metrics to tune concurrency.

- [ ] **Step 2: Run focused tests repeatedly**

Run the factorized-training test file twice under the deterministic environment. Expected: identical pass counts and no flaky process-pool failures.

- [ ] **Step 3: Run the full benchmark Python suite**

Run:

```bash
PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
RAYON_NUM_THREADS=1 benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: all tests pass.

- [ ] **Step 4: Run static and scope checks**

Run:

```bash
benchmark/.venv/bin/python -m py_compile \
  benchmark/factorized_ranker/training.py benchmark/train_factorized_ranker.py
git diff --check -- benchmark/factorized_ranker/training.py \
  benchmark/train_factorized_ranker.py benchmark/tests/test_factorized_training.py \
  .superpowers/sdd/2026-08-07-factorized-structural-ranker/task-13-brief.md
git status --short
git diff --stat
```

Expected: only intentional Task 13 parallelism files are changed; unrelated dirt remains untouched.

- [ ] **Step 5: Commit the verified command contract**

```bash
git add .superpowers/sdd/2026-08-07-factorized-structural-ranker/task-13-brief.md
git commit -m "docs: freeze parallel factorized training command"
```

Record the resulting final implementation HEAD; all rebuilt provenance must bind this commit.

### Task 6: Rebuild source-bound corpus/fold artifacts twice

**Files:**
- Regenerate: `benchmark/data/cath17287_factorized_corpus_v1/`
- Regenerate: `benchmark/data/cath17287_factorized_corpus_v1_verify/`
- Regenerate: `benchmark/data/factorized_ranker_v1_freeze_stage/cath17287_factorized_corpus_v1_manifest.json`
- Regenerate: `benchmark/data/factorized_ranker_v1_freeze_stage/cath17287_factorized_folds_v1.json`
- Regenerate: `benchmark/data/factorized_ranker_v1_freeze_stage/cath17287_factorized_folds_v1.verify.json`

- [ ] **Step 1: Snapshot scientific hashes and validate raw acquisition inputs**

Hash `chains.csv`, `counts.csv`, `candidates.csv`, and `rejections.csv` from the current canonical corpus into the run ledger. Verify the raw candidate dump, acquisition rejection/provenance files, release binary, and expected acquisition commit `b384b7d5ba4d63c12153c580a7b20c71a8c7532e` still match the existing accepted manifest. Do not rerun acquisition.

- [ ] **Step 2: Remove only obsolete source-bound normalized/fold outputs**

Resolve every path explicitly, ensure each is under `$PWD/benchmark/data`, then remove the two normalized corpus directories and the three freeze-stage manifest/fold files. Preserve the raw dump, raw parts, acquisition files, cache, and model/baseline evidence.

- [ ] **Step 3: Build the canonical and independent verify corpora**

Run twice, changing only `--out-dir`:

```bash
benchmark/.venv/bin/python -m benchmark.build_factorized_corpus \
  --dataset cath17287 \
  --dump benchmark/data/cath17287_factorized_candidates_v1.csv \
  --chain-cache-dir benchmark/cache/chains \
  --out-dir benchmark/data/cath17287_factorized_corpus_v1 \
  --binary target/release/sword2 \
  --expected-acquisition-commit b384b7d5ba4d63c12153c580a7b20c71a8c7532e \
  --rejections benchmark/data/cath17287_factorized_candidates_v1.rejections.csv
```

Expected: 10,573 accepted chains, 84,617 count rows, 211,446 candidate rows, and complete accounting for all 17,286 metadata rows.

- [ ] **Step 4: Require exact independent corpus reproduction and scientific identity**

Run `cmp` for `chains.csv`, `counts.csv`, `candidates.csv`, `rejections.csv`, and `corpus_manifest.json` between canonical and verify directories. Compare the four new scientific-table hashes to the Task 6 Step 1 snapshot; require exact equality. Only provenance/source-bound manifest fields may legitimately differ from the pre-implementation manifest.

- [ ] **Step 5: Atomically stage the canonical manifest and build folds twice**

Copy the canonical manifest through a `.tmp` plus `mv`, then run:

```bash
benchmark/.venv/bin/python -m benchmark.build_factorized_folds \
  --dataset cath17287 \
  --corpus-dir benchmark/data/cath17287_factorized_corpus_v1 \
  --out benchmark/data/factorized_ranker_v1_freeze_stage/cath17287_factorized_folds_v1.json \
  --n-folds 5 --seed 37
benchmark/.venv/bin/python -m benchmark.build_factorized_folds \
  --dataset cath17287 \
  --corpus-dir benchmark/data/cath17287_factorized_corpus_v1 \
  --out benchmark/data/factorized_ranker_v1_freeze_stage/cath17287_factorized_folds_v1.verify.json \
  --n-folds 5 --seed 37
cmp benchmark/data/factorized_ranker_v1_freeze_stage/cath17287_factorized_folds_v1.json \
  benchmark/data/factorized_ranker_v1_freeze_stage/cath17287_factorized_folds_v1.verify.json
```

Expected: exact bytes, all 10,573 accepted chains once, all five folds populated, and no PDB/family/component leakage. Compare assignment content to the prior fold snapshot and require scientific identity.

### Task 7: Launch and validate the eight-worker production run

**Files:**
- Generate: `benchmark/data/factorized_ranker_v1_training/`
- Generate: `benchmark/data/factorized_ranker_v1_run_ledger/training.stdout.log`
- Generate: `benchmark/data/factorized_ranker_v1_run_ledger/training.stderr.log`
- Generate: `benchmark/data/factorized_ranker_v1_run_ledger/training.time.txt`
- Generate: `benchmark/data/factorized_ranker_v1_run_ledger/training.pid`

- [ ] **Step 1: Validate a clean production target and committed source**

Require no training directory, no staged count/candidate model output from the obsolete run, a clean diff for every training source in `_TRAINING_SOURCE_FILES`, and the final committed HEAD. Ensure the run-ledger directory exists and the four new log/PID targets do not.

- [ ] **Step 2: Launch one durable exact command**

Run with inherited deterministic environment:

```bash
mkdir -p benchmark/data/factorized_ranker_v1_run_ledger
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export RAYON_NUM_THREADS=1
nohup /usr/bin/time -v \
  -o benchmark/data/factorized_ranker_v1_run_ledger/training.time.txt \
  benchmark/.venv/bin/python -m benchmark.train_factorized_ranker \
  --corpus-dir benchmark/data/cath17287_factorized_corpus_v1 \
  --fold-manifest benchmark/data/factorized_ranker_v1_freeze_stage/cath17287_factorized_folds_v1.json \
  --out-dir benchmark/data/factorized_ranker_v1_training \
  --count-model-out benchmark/data/factorized_ranker_v1_freeze_stage/factorized_count_v1.json \
  --candidate-model-out benchmark/data/factorized_ranker_v1_freeze_stage/factorized_candidate_v1.json \
  --seed 37 --jobs 8 --resume \
  >benchmark/data/factorized_ranker_v1_run_ledger/training.stdout.log \
  2>benchmark/data/factorized_ranker_v1_run_ledger/training.stderr.log \
  </dev/null &
printf '%s\n' "$!" >benchmark/data/factorized_ranker_v1_run_ledger/training.pid
```

The user can inspect progress with:

```bash
tail -f benchmark/data/factorized_ranker_v1_run_ledger/training.stderr.log
```

- [ ] **Step 3: Verify process-level parallelism without reading metrics**

Resolve the wrapper PID from `training.pid`, inspect its descendant tree with `ps`, and require one Python parent plus eight compute-bound Python workers during a full grid wave. Confirm all seven thread-limit variables through `/proc/<python-parent>/environ` without printing unrelated environment values.

- [ ] **Step 4: Enforce the first-wave memory gate**

Use `ps -eo pid=,ppid=,rss=,args=` to build the wrapper descendant set and sum RSS in bytes. Sample throughout the first full eight-job wave. Require every total below 50,361,477,120 bytes. If breached, stop the exact process group as an infrastructure failure before accepting checkpoints; do not inspect scores or silently lower jobs.

- [ ] **Step 5: Validate first-wave checkpoint health**

After the first eight fitted completions, require at least eight canonical grid JSON files, unique context-derived names, exact context containing `--jobs`, `8`, the final source commit and hashes, no `.tmp` remnants, and successful `TrainingCheckpointStore.load_fold` validation. Progress logs may expose task identity/count only, never selection metrics.

- [ ] **Step 6: Hand off the durable run**

Report the wrapper and Python parent PIDs, current worker count/CPU/RSS, checkpoint count, fixed memory limit, exact log paths, and `tail` command. If training is still running, do not claim scientific completion; explain that disconnecting chat is safe but powering off requires later resume with the identical command.
