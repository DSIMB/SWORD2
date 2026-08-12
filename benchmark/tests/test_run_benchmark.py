from __future__ import annotations

import csv
import io
from pathlib import Path
from types import SimpleNamespace

import pytest

import benchmark.run_benchmark as run_benchmark
from benchmark.run_benchmark import (
    LockedBenchmarkError,
    _eligible_resource_assignments,
    _gnu_time_version,
    _normalized_command,
    _split_locked_merizo_stdout,
    _validate_locked_factorized_outcome,
    _validate_locked_jobs,
    _validate_locked_role_evidence,
    counterbalanced_pair_orders,
    parse_selector_status,
)
from benchmark.runners.base import parse_locked_peak_rss_kb
from benchmark.runners.sword2_rust import Sword2RustRunner
from benchmark.score import LockedRunRow, RunRow, write_locked_runs_csv, write_runs_csv


LEGACY_STATUS = b'{"error_code":null,"excluded_candidate_count":0,"fallback":false,"requested_selector":"legacy","schema_version":1,"selector_used":"legacy"}\n'
FACTORIZED_STATUS = b'{"error_code":null,"excluded_candidate_count":2,"fallback":false,"requested_selector":"factorized","schema_version":1,"selector_used":"factorized"}\n'
FALLBACK_STATUS = b'{"error_code":"feature_missing_context","excluded_candidate_count":2,"fallback":true,"requested_selector":"factorized","schema_version":1,"selector_used":"legacy"}\n'
STRUCTURAL_ABSTENTION_STATUS = b'{"error_code":"structural_quality_abstention","excluded_candidate_count":0,"fallback":true,"requested_selector":"factorized","schema_version":1,"selector_used":"legacy"}\n'


def test_installed_gnu_time_version_accepts_distribution_capitalization():
    assert "gnu time" in _gnu_time_version().casefold()


def test_locked_artifact_recheck_paths_preserve_symlink_identity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    target = tmp_path / "python-real"
    target.write_bytes(b"python")
    link = tmp_path / "python-link"
    link.symlink_to(target.name)
    monkeypatch.chdir(tmp_path)

    paths = run_benchmark._artifact_paths_for_recheck({"python": Path("python-link")})

    assert paths == {"python": link}
    assert paths["python"].is_symlink()


def test_counterbalanced_orders_are_seed37_permutation_invariant():
    expected = counterbalanced_pair_orders(["c", "a", "b", "d"])
    assert expected == counterbalanced_pair_orders(["d", "b", "a", "c"])
    assert sorted(expected) == ["a", "b", "c", "d"]
    counts = {selector: list(expected.values()).count(selector) for selector in set(expected.values())}
    assert set(counts) == {"legacy", "factorized"}
    assert abs(counts["legacy"] - counts["factorized"]) <= 1


def test_locked_rss_and_selector_status_are_exact_and_fail_closed():
    assert parse_locked_peak_rss_kb("Maximum resident set size (kbytes): 123\n") == 123
    for invalid in [
        "",
        "Maximum resident set size (kbytes): 0\n",
        "Maximum resident set size (bytes): 123\n",
        "Maximum resident set size (kbytes): 1\nMaximum resident set size (kbytes): 2\n",
    ]:
        with pytest.raises(ValueError):
            parse_locked_peak_rss_kb(invalid)

    assert parse_selector_status(LEGACY_STATUS, "legacy")["selector_used"] == "legacy"
    assert parse_selector_status(FACTORIZED_STATUS, "factorized")["excluded_candidate_count"] == 2
    with pytest.raises(ValueError):
        parse_selector_status(FACTORIZED_STATUS.replace(b'"fallback":false', b'"fallback":true'), "factorized")
    rejected = parse_selector_status(
        FALLBACK_STATUS,
        "factorized",
        require_success=False,
    )
    assert rejected["fallback"] is True
    assert rejected["error_code"] == "feature_missing_context"
    with pytest.raises(ValueError, match="requested selector"):
        parse_selector_status(FALLBACK_STATUS, "factorized")
    with pytest.raises(ValueError):
        parse_selector_status(FACTORIZED_STATUS.replace(b'"schema_version":1', b'"schema_version":1,"schema_version":1'), "factorized")


def test_structural_abstention_status_is_valid_only_when_success_not_required():
    parsed = parse_selector_status(
        STRUCTURAL_ABSTENTION_STATUS,
        "factorized",
        require_success=False,
    )
    assert parsed["error_code"] == "structural_quality_abstention"
    with pytest.raises(ValueError):
        parse_selector_status(
            STRUCTURAL_ABSTENTION_STATUS,
            "factorized",
            require_success=True,
        )


def _factorized_row(
    entry_id: str,
    *,
    success: bool,
    code: str = "structural_quality_abstention",
    exclusions: int = 0,
) -> SimpleNamespace:
    return SimpleNamespace(
        entry_id=entry_id,
        requested_selector="factorized",
        selector_used="factorized" if success else "legacy",
        fallback=not success,
        error_code="" if success else code,
        selector_warning_code="" if success else "factorized_fallback",
        excluded_candidate_count=exclusions,
    )


def test_factorized_outcome_scores_eligible_and_retains_ineligible_run_only():
    eligible = frozenset({"a"})
    assert _validate_locked_factorized_outcome(
        "a", _factorized_row("a", success=True), eligible
    )
    assert not _validate_locked_factorized_outcome(
        "b", _factorized_row("b", success=False), eligible
    )


@pytest.mark.parametrize(
    ("entry_id", "row"),
    [
        ("a", _factorized_row("a", success=False)),
        ("b", _factorized_row("b", success=True)),
        (
            "b",
            _factorized_row(
                "b", success=False, code="feature_missing_context"
            ),
        ),
        ("b", _factorized_row("b", success=False, exclusions=1)),
    ],
)
def test_expected_abstention_rejects_wrong_id_error_or_exclusion(
    entry_id: str, row: SimpleNamespace
):
    with pytest.raises(LockedBenchmarkError):
        _validate_locked_factorized_outcome(entry_id, row, frozenset({"a"}))


def test_resource_assignments_use_only_eligible_ids():
    assignments = _eligible_resource_assignments(
        dataset_ids=("a", "b", "c"),
        eligible_ids=frozenset({"a", "c"}),
    )
    assert set(assignments) == {"a", "c"}
    assert "b" not in assignments


def test_locked_worker_contract_is_role_specific():
    assert _validate_locked_jobs("factorized-accuracy", 32) == 32
    assert _validate_locked_jobs("legacy-accuracy", 32) == 32
    assert _validate_locked_jobs("paired-sword-resources", 1) == 1
    with pytest.raises(LockedBenchmarkError):
        _validate_locked_jobs("factorized-accuracy", 31)
    with pytest.raises(LockedBenchmarkError):
        _validate_locked_jobs("paired-sword-resources", 2)


def test_factorized_accuracy_parallelizes_runs_and_does_not_score_abstentions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    structures = tuple(
        SimpleNamespace(
            entry=SimpleNamespace(entry_id=entry_id),
            numbering=SimpleNamespace(),
        )
        for entry_id in ("a", "b")
    )
    preflight = SimpleNamespace(structures=structures, eligible_ids=("a",))
    args = SimpleNamespace(
        locked_role="factorized-accuracy",
        locked_jobs=32,
        results_dir=tmp_path,
    )
    executor_workers: list[int] = []
    real_executor = run_benchmark.ThreadPoolExecutor

    def executor(*, max_workers: int):
        executor_workers.append(max_workers)
        return real_executor(max_workers=max_workers)

    def locked_run(structure, selector, process_dir, **kwargs):
        entry_id = structure.entry.entry_id
        row = _factorized_row(entry_id, success=entry_id == "a")
        row.tool = "sword2-rust"
        row.selector_variant = selector
        row.pair_position = ""
        row.runtime_s = 1.0
        row.peak_rss_kb = 1024
        return [SimpleNamespace(name="prediction")], row

    scored: list[str] = []

    def score(entry, numbering, prediction, **kwargs):
        scored.append(entry.entry_id)
        return SimpleNamespace(
            entry_id=entry.entry_id,
            tool="sword2-rust",
            variant="optimal",
            partition="Optimal partition",
        )

    monkeypatch.setattr(run_benchmark, "ThreadPoolExecutor", executor)
    monkeypatch.setattr(run_benchmark, "_locked_sword_run", locked_run)
    monkeypatch.setattr(run_benchmark, "score_prediction", score)
    monkeypatch.setattr(run_benchmark, "_atomic_write_rows", lambda *args: None)

    scores, runs, failures, order = run_benchmark._run_locked_accuracy(
        preflight, args
    )

    assert executor_workers == [32]
    assert scored == ["a"]
    assert [row.entry_id for row in runs] == ["a", "b"]
    assert [row.entry_id for row in scores] == ["a"]
    assert failures == []
    assert order is None


def test_locked_factorized_role_evidence_matches_frozen_eligibility_sets():
    eligible = _factorized_row("a", success=True)
    eligible.tool = "sword2-rust"
    eligible.selector_variant = "factorized"
    abstained = _factorized_row("b", success=False)
    abstained.tool = "sword2-rust"
    abstained.selector_variant = "factorized"
    preflight = SimpleNamespace(
        dataset_ids=("a", "b"),
        eligible_ids=("a",),
        ineligible_ids=("b",),
    )
    args = SimpleNamespace(locked_role="factorized-accuracy")
    scores = [SimpleNamespace(entry_id="a", tool="sword2-rust")]

    assert _validate_locked_role_evidence(
        preflight, args, scores=scores, runs=[eligible, abstained]
    ) == frozenset({"b"})

    with pytest.raises(LockedBenchmarkError):
        _validate_locked_role_evidence(
            preflight,
            args,
            scores=[*scores, SimpleNamespace(entry_id="b", tool="sword2-rust")],
            runs=[eligible, abstained],
        )


def test_locked_intent_binds_eligibility_population_and_worker_contract():
    preflight = SimpleNamespace(
        dataset_sha256="d" * 64,
        dataset_ids=("a", "b"),
        dataset_id_set_sha256="i" * 64,
        runtime_manifest_sha256="r" * 64,
        model_manifest_sha256="m" * 64,
        runtime_manifest={"binary_sha256": "b" * 64},
        structure_tree_sha256="s" * 64,
        selected_tools=("sword2-rust",),
        external_artifacts={},
        eligibility_manifest_sha256="e" * 64,
        eligibility_policy="strict_complete_backbone_v1",
        eligible_ids=("a",),
        ineligible_ids=("b",),
        eligibility_manifest={
            "eligible_id_set_sha256": "a" * 64,
            "ineligible_id_set_sha256": "c" * 64,
        },
    )
    args = SimpleNamespace(
        locked_role="factorized-accuracy",
        dataset="cath663",
        locked_jobs=32,
    )

    payload = run_benchmark._locked_intent_payload(
        preflight, args, "2026-08-12T00:00:00+00:00"
    )

    assert payload["schema_version"] == 2
    assert payload["eligibility_manifest_sha256"] == "e" * 64
    assert payload["factorized_eligible_count"] == 1
    assert payload["structural_abstention_count"] == 1
    assert payload["locked_jobs"] == 32


def test_locked_parser_requires_explicit_eligibility_and_worker_inputs():
    args = run_benchmark.parse_args(
        [
            "--locked-factorized-manifest",
            "runtime.json",
            "--locked-eligibility-manifest",
            "eligibility.json",
            "--locked-role",
            "factorized-accuracy",
            "--locked-jobs",
            "32",
        ]
    )
    assert args.locked_eligibility_manifest == Path("eligibility.json")
    assert args.locked_jobs == 32


def test_locked_mode_rejects_a_partial_authority_tuple(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(
        run_benchmark.sys,
        "argv",
        [
            "run_benchmark",
            "--locked-factorized-manifest",
            "runtime.json",
            "--locked-eligibility-manifest",
            "eligibility.json",
            "--locked-role",
            "factorized-accuracy",
        ],
    )

    assert run_benchmark.main() == 2


def test_locked_run_header_is_separate_and_ordinary_header_is_unchanged(tmp_path: Path):
    ordinary_path = tmp_path / "ordinary.csv"
    locked_path = tmp_path / "locked.csv"
    write_runs_csv([], ordinary_path)
    write_locked_runs_csv([], locked_path)
    ordinary_header = next(csv.reader(io.StringIO(ordinary_path.read_text())))
    locked_header = next(csv.reader(io.StringIO(locked_path.read_text())))
    assert ordinary_header == list(RunRow.__dataclass_fields__)
    assert locked_header == list(LockedRunRow.__dataclass_fields__)
    assert "peak_rss_kb" not in ordinary_header
    assert locked_header[0:6] == ["dataset", "entry_id", "pdb_id", "chain_id", "tool", "selector_variant"]


def test_fresh_only_sword_runner_refuses_existing_output(tmp_path: Path):
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "keep"
    marker.write_text("owned")
    runner = Sword2RustRunner(repo_dir=tmp_path, binary=tmp_path / "missing", fresh_only=True)
    with pytest.raises(FileExistsError):
        runner.run(tmp_path / "input.pdb", output)
    assert marker.read_text() == "owned"


def test_locked_top_level_command_normalizes_relative_and_embedded_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.chdir(tmp_path)
    for path in (tmp_path / "cache", tmp_path / "tool"):
        path.mkdir()
    for path in (
        tmp_path / "sword2",
        tmp_path / "runtime.json",
        tmp_path / "eligibility.json",
        tmp_path / "ids.txt",
        tmp_path / "tool/python",
    ):
        path.write_text("x")
    preflight = SimpleNamespace(
        repo_root=tmp_path.resolve(),
        binary=(tmp_path / "sword2").resolve(),
        structures=(),
        external_artifact_paths={"merizo_python": (tmp_path / "tool/python").resolve()},
    )
    args = SimpleNamespace(
        cache_dir=Path("cache"),
        results_dir=(tmp_path / "results").absolute(),
        locked_factorized_manifest=Path("runtime.json"),
        locked_eligibility_manifest=Path("eligibility.json"),
        chainsaw_expected_success_ids=Path("ids.txt"),
        locked_artifact=["merizo_python=tool/python"],
    )
    command = [
        "benchmark.run_benchmark",
        "--cache-dir", "cache",
        "--results-dir", "results",
        "--locked-factorized-manifest", "runtime.json",
        "--locked-eligibility-manifest", "eligibility.json",
        "--chainsaw-expected-success-ids", "ids.txt",
        "--locked-artifact", "merizo_python=tool/python",
    ]
    normalized = _normalized_command(command, preflight=preflight, args=args)
    assert normalized == [
        "benchmark.run_benchmark",
        "--cache-dir", "<cache>",
        "--results-dir", "<results>",
        "--locked-factorized-manifest", "<runtime_manifest>",
        "--locked-eligibility-manifest", "<eligibility_manifest>",
        "--chainsaw-expected-success-ids", "<chainsaw_expected_success_ids>",
        "--locked-artifact", "merizo_python=<locked_artifact:merizo_python>",
    ]
    assert str(tmp_path) not in repr(normalized)


def test_locked_merizo_split_is_absent_only_and_exact(tmp_path: Path):
    inputs = [
        (tmp_path / "a.pdb", tmp_path / "out/a.tsv", "A"),
        (tmp_path / "b.pdb", tmp_path / "out/b.tsv", "A"),
    ]
    stdout = "input\tresult\tndom\n/tmp/a.pdb\t1-10\t1\n/tmp/b.pdb\t1-20\t1\n"
    _split_locked_merizo_stdout(stdout, inputs)
    assert inputs[0][1].read_text() == "input\tresult\tndom\n/tmp/a.pdb\t1-10\t1\n"
    with pytest.raises(ValueError, match="already exists"):
        _split_locked_merizo_stdout(stdout, inputs)
    assert inputs[0][1].read_text().endswith("/tmp/a.pdb\t1-10\t1\n")

    fresh = [(tmp_path / "a.pdb", tmp_path / "fresh/a.tsv", "A")]
    with pytest.raises(ValueError, match="unknown input"):
        _split_locked_merizo_stdout(
            "input\tresult\n/tmp/other.pdb\t1-10\n",
            fresh,
        )


def test_locked_merizo_split_accepts_exact_repeated_batch_headers(tmp_path: Path):
    inputs = [
        (tmp_path / "a.pdb", tmp_path / "out/a.tsv", "A"),
        (tmp_path / "b.pdb", tmp_path / "out/b.tsv", "A"),
    ]
    stdout = (
        "input\tresult\tndom\n"
        "/tmp/a.pdb\t1-10\t1\n"
        "input\tresult\tndom\n"
        "/tmp/b.pdb\t1-20\t1\n"
    )

    _split_locked_merizo_stdout(stdout, inputs)

    assert inputs[0][1].read_text() == (
        "input\tresult\tndom\n/tmp/a.pdb\t1-10\t1\n"
    )
    assert inputs[1][1].read_text() == (
        "input\tresult\tndom\n/tmp/b.pdb\t1-20\t1\n"
    )

    altered = [(tmp_path / "a.pdb", tmp_path / "altered/a.tsv", "A")]
    with pytest.raises(ValueError, match="unknown input"):
        _split_locked_merizo_stdout(
            "input\tresult\tndom\n"
            "/tmp/a.pdb\t1-10\t1\n"
            "input\tndom\tresult\n",
            altered,
        )
