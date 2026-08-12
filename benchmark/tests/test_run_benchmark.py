from __future__ import annotations

import csv
import io
from pathlib import Path
from types import SimpleNamespace

import pytest

import benchmark.run_benchmark as run_benchmark
from benchmark.run_benchmark import (
    _normalized_command,
    _split_locked_merizo_stdout,
    counterbalanced_pair_orders,
    parse_selector_status,
)
from benchmark.runners.base import parse_locked_peak_rss_kb
from benchmark.runners.sword2_rust import Sword2RustRunner
from benchmark.score import LockedRunRow, RunRow, write_locked_runs_csv, write_runs_csv


LEGACY_STATUS = b'{"error_code":null,"excluded_candidate_count":0,"fallback":false,"requested_selector":"legacy","schema_version":1,"selector_used":"legacy"}\n'
FACTORIZED_STATUS = b'{"error_code":null,"excluded_candidate_count":2,"fallback":false,"requested_selector":"factorized","schema_version":1,"selector_used":"factorized"}\n'
FALLBACK_STATUS = b'{"error_code":"feature_missing_context","excluded_candidate_count":2,"fallback":true,"requested_selector":"factorized","schema_version":1,"selector_used":"legacy"}\n'


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
    for path in (tmp_path / "sword2", tmp_path / "runtime.json", tmp_path / "ids.txt", tmp_path / "tool/python"):
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
        chainsaw_expected_success_ids=Path("ids.txt"),
        locked_artifact=["merizo_python=tool/python"],
    )
    command = [
        "benchmark.run_benchmark",
        "--cache-dir", "cache",
        "--results-dir", "results",
        "--locked-factorized-manifest", "runtime.json",
        "--chainsaw-expected-success-ids", "ids.txt",
        "--locked-artifact", "merizo_python=tool/python",
    ]
    normalized = _normalized_command(command, preflight=preflight, args=args)
    assert normalized == [
        "benchmark.run_benchmark",
        "--cache-dir", "<cache>",
        "--results-dir", "<results>",
        "--locked-factorized-manifest", "<runtime_manifest>",
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
