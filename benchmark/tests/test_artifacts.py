from pathlib import Path
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

import benchmark.figures as figures
from benchmark.datasets import CathEntry
from benchmark.numbering import ResidueKey, StructureNumbering
from benchmark.runners.base import PartitionPrediction, ToolRunResult
import benchmark.run_benchmark as run_benchmark
import benchmark.score as score
import benchmark.summaries as summaries


def _score_row(**overrides) -> score.ScoreRow:
    values = {
        "dataset": "cath663",
        "entry_id": "19hcA",
        "pdb_id": "19hc",
        "chain_id": "A",
        "tool": "merizo",
        "variant": "cpu",
        "partition": "Merizo",
        "true_chopping": "0-9",
        "pred_chopping": "0-9",
        "n_residues": 10,
        "n_true_domains": 1,
        "n_pred_domains": 1,
        "ndo": 1.0,
        "boundary_dist_score": 1.0,
        "d_count_acc": 1.0,
        "d_count_dev": 0.0,
        "iou": 1.0,
        "multi_ndo": 1.0,
        "domain_count_bias": 0.0,
        "over_split": 0.0,
        "merge": 0.0,
        "boundary_precision_5": 1.0,
        "boundary_recall_5": 1.0,
        "boundary_f1_5": 1.0,
        "boundary_precision_10": 1.0,
        "boundary_recall_10": 1.0,
        "boundary_f1_10": 1.0,
        "boundary_precision_20": 1.0,
        "boundary_recall_20": 1.0,
        "boundary_f1_20": 1.0,
        "median_boundary_error": 0.0,
        "worst_boundary_error": 0.0,
        "pred_coverage": 1.0,
        "pred_linker_fraction": 0.0,
        "pairwise_precision": 1.0,
        "pairwise_recall": 1.0,
        "pairwise_f1": 1.0,
        "adjusted_rand": 1.0,
        "normalized_mutual_info": 1.0,
        "variation_of_information": 0.0,
        "matched_dice": 1.0,
        "matched_jaccard": 1.0,
        "exact_match": 1.0,
        "runtime_s": None,
        "peak_rss_mb": None,
    }
    values.update(overrides)
    return score.ScoreRow(**values)


def test_parse_args_defaults_merizo_to_cuda(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["run_benchmark.py"])

    args = run_benchmark.parse_args()

    assert args.merizo_device == "cuda"
    assert args.allow_dl_cpu is False


def test_parse_args_accepts_reuse_tool_results_dir(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        sys,
        "argv",
        ["run_benchmark.py", "--reuse-tool-results-dir", str(tmp_path / "reference")],
    )

    args = run_benchmark.parse_args()

    assert args.reuse_tool_results_dir == tmp_path / "reference"


def test_headline_scores_keeps_non_sword_cuda_predictions():
    scores = pd.DataFrame(
        [
            {"entry_id": "19hcA", "tool": "sword2-rust", "variant": "optimal", "ndo": 1.0},
            {"entry_id": "19hcA", "tool": "sword2-rust", "variant": "alternative", "ndo": 0.8},
            {"entry_id": "19hcA", "tool": "merizo", "variant": "cuda", "ndo": 0.7},
            {"entry_id": "19hcA", "tool": "chainsaw", "variant": "single", "ndo": 0.6},
        ]
    )

    headline = summaries.headline_scores(scores)

    assert set(zip(headline["tool"], headline["variant"])) == {
        ("sword2-rust", "optimal"),
        ("merizo", "cuda"),
        ("chainsaw", "single"),
    }


def test_run_rows_are_written_as_one_row_per_tool_invocation(tmp_path: Path):
    assert hasattr(score, "RunRow")
    assert hasattr(score, "write_runs_csv")
    path = tmp_path / "runs.csv"
    rows = [
        score.RunRow(
            dataset="cath663",
            entry_id="19hcA",
            pdb_id="19hc",
            chain_id="A",
            tool="sword2-rust",
            variant="all-partitions",
            reused_output=False,
            returncode=0,
            runtime_s=0.2,
            peak_rss_mb=12.5,
            command="sword2 -i 19hcA.pdb",
            cwd="/repo",
            raw_path="benchmark/results/raw/sword2-rust/19hcA/summary.json",
        )
    ]

    score.write_runs_csv(rows, path)

    written = pd.read_csv(path)
    assert written.to_dict("records")[0]["entry_id"] == "19hcA"
    assert written.to_dict("records")[0]["runtime_s"] == 0.2


def test_failure_rows_capture_non_strict_tool_errors(tmp_path: Path):
    assert hasattr(score, "FailureRow")
    assert hasattr(score, "write_failures_csv")
    path = tmp_path / "failures.csv"
    rows = [
        score.FailureRow(
            dataset="cath663",
            entry_id="19hcA",
            pdb_id="19hc",
            chain_id="A",
            tool="merizo",
            stage="run",
            message="Merizo failed",
            returncode=1,
            command="python predict.py",
            stderr="traceback",
        )
    ]

    score.write_failures_csv(rows, path)

    written = pd.read_csv(path)
    assert written.to_dict("records")[0]["tool"] == "merizo"
    assert written.to_dict("records")[0]["returncode"] == 1


def test_run_benchmark_builds_run_and_failure_rows_from_tool_results():
    assert hasattr(run_benchmark, "_run_row_from_result")
    assert hasattr(run_benchmark, "_failure_row_from_exception")
    entry = CathEntry(
        pdb_id="19hc",
        chain_id="A",
        entry_id="19hcA",
        n_domains=2,
        n_residues=292,
        chopping="1-10|11-20",
        dataset="cath663",
    )
    structure = SimpleNamespace(entry=entry)
    result = ToolRunResult(
        command=["sword2", "-i", "19hcA.pdb"],
        cwd=Path("/repo"),
        returncode=0,
        runtime_s=0.25,
        peak_rss_mb=13.0,
        stdout="",
        stderr="",
    )

    run_row = run_benchmark._run_row_from_result(
        structure,
        tool="sword2-rust",
        variant="all-partitions",
        raw_path=Path("benchmark/results/raw/sword2-rust/19hcA/summary.json"),
        result=result,
        reused_output=False,
    )
    failure_row = run_benchmark._failure_row_from_exception(
        structure,
        tool="sword2-rust",
        exc=RuntimeError("failed"),
    )

    assert run_row.entry_id == "19hcA"
    assert run_row.command == "sword2 -i 19hcA.pdb"
    assert run_row.runtime_s == 0.25
    assert failure_row.message == "failed"


def test_batched_deep_learning_runtime_uses_tool_warm_runtime():
    entry = CathEntry(
        pdb_id="19hc",
        chain_id="A",
        entry_id="19hcA",
        n_domains=1,
        n_residues=10,
        chopping="1-10",
        dataset="cath663",
    )
    numbering = StructureNumbering(
        residues=[ResidueKey("A", str(index), "") for index in range(1, 11)]
    )
    structure = SimpleNamespace(entry=entry, numbering=numbering)
    prediction = PartitionPrediction(
        tool="merizo",
        name="Merizo",
        variant="cpu",
        chopping="0-9",
        n_domains=1,
        raw_path=Path("benchmark/results/raw/merizo-cpu/19hcA.tsv"),
        runtime_s=0.12,
    )
    process_result = ToolRunResult(
        command=["python", "predict.py", "-i", "19hcA.pdb", "1a59A.pdb"],
        cwd=Path("/repo/Merizo"),
        returncode=0,
        runtime_s=5.0,
        peak_rss_mb=512.0,
        stdout="",
        stderr="",
    )

    score_rows, run_row = run_benchmark._rows_from_predictions(
        structure,
        tool="merizo",
        predictions=[prediction],
        result=process_result,
        reused_output=False,
        batch_n_entries=2,
    )

    assert score_rows[0].runtime_s == 0.12
    assert score_rows[0].peak_rss_mb is None
    assert run_row.runtime_s == 0.12
    assert run_row.peak_rss_mb is None
    assert run_row.batch_runtime_s == 5.0
    assert run_row.batch_peak_rss_mb == 512.0
    assert run_row.batch_n_entries == 2


def test_merizo_batch_orchestration_runs_one_process_for_multiple_structures(monkeypatch, tmp_path: Path):
    structures = [
        SimpleNamespace(
            entry=CathEntry("19hc", "A", "19hcA", 1, 10, "1-10", "cath663"),
            numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 11)]),
            pdb_path=tmp_path / "19hcA.pdb",
        ),
        SimpleNamespace(
            entry=CathEntry("1a59", "A", "1a59A", 1, 12, "1-12", "cath663"),
            numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 13)]),
            pdb_path=tmp_path / "1a59A.pdb",
        ),
    ]
    for structure in structures:
        structure.pdb_path.write_text("ATOM\n")

    calls = []

    class FakeMerizoRunner:
        def __init__(self, device, require_cuda=True):
            self.device = device
            self.require_cuda = require_cuda

        def run_batch(self, inputs):
            calls.append(inputs)
            rows = {
                "19hcA": "19hcA.pdb\t10\t10\t0\t1\t0.9\t0.11\t1-10\n",
                "1a59A": "1a59A.pdb\t12\t12\t0\t1\t0.8\t0.22\t1-12\n",
            }
            for structure_file, output_tsv, _chain_id in inputs:
                output_tsv.parent.mkdir(parents=True, exist_ok=True)
                output_tsv.write_text(
                    "input\tnres\tnres_dom\tnres_ndr\tndom\tpIoU\truntime\tresult\n"
                    + rows[structure_file.stem]
                )
            return ToolRunResult(["python", "predict.py"], tmp_path, 0, 7.0, 700.0, "", "")

    monkeypatch.setattr(run_benchmark, "MerizoRunner", FakeMerizoRunner)
    args = SimpleNamespace(results_dir=tmp_path / "results", merizo_device="cpu", skip_existing=False)

    batched = run_benchmark._maybe_run_merizo_batch(structures, args)

    assert len(calls) == 1
    assert len(calls[0]) == 2
    predictions, result, reused_output, batch_n_entries = batched["19hcA"]
    assert predictions[0].runtime_s == 0.11
    assert result.runtime_s == 7.0
    assert reused_output is False
    assert batch_n_entries == 2


def test_merizo_batch_reuses_copied_tsv_from_reference_results_dir(monkeypatch, tmp_path: Path):
    structure = SimpleNamespace(
        entry=CathEntry("19hc", "A", "19hcA", 1, 10, "1-10", "cath663"),
        numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 11)]),
        pdb_path=tmp_path / "19hcA.pdb",
    )
    source_tsv = tmp_path / "reference" / "raw" / "merizo-cpu" / "19hcA.tsv"
    source_tsv.parent.mkdir(parents=True)
    source_tsv.write_text(
        "input\tnres\tnres_dom\tnres_ndr\tndom\tpIoU\truntime\tresult\n"
        "19hcA.pdb\t10\t10\t0\t1\t0.9\t0.11\t1-10\n"
    )

    class UnexpectedMerizoRunner:
        def __init__(self, *_args, **_kwargs):
            raise AssertionError("reference outputs should skip the Merizo runner")

    monkeypatch.setattr(run_benchmark, "MerizoRunner", UnexpectedMerizoRunner)
    args = SimpleNamespace(
        results_dir=tmp_path / "results",
        merizo_device="cpu",
        skip_existing=False,
        reuse_tool_results_dir=tmp_path / "reference",
    )

    batched = run_benchmark._maybe_run_merizo_batch([structure], args)

    copied_tsv = tmp_path / "results" / "raw" / "merizo-cpu" / "19hcA.tsv"
    assert copied_tsv.read_text() == source_tsv.read_text()
    predictions, result, reused_output, batch_n_entries = batched["19hcA"]
    assert predictions[0].raw_path == copied_tsv
    assert predictions[0].runtime_s == 0.11
    assert result is None
    assert reused_output is True
    assert batch_n_entries is None


def test_chainsaw_batch_orchestration_runs_one_process_for_multiple_structures(monkeypatch, tmp_path: Path):
    structures = [
        SimpleNamespace(
            entry=CathEntry("19hc", "A", "19hcA", 1, 10, "1-10", "cath663"),
            numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 11)]),
            pdb_path=tmp_path / "19hcA.pdb",
        ),
        SimpleNamespace(
            entry=CathEntry("1a59", "A", "1a59A", 1, 12, "1-12", "cath663"),
            numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 13)]),
            pdb_path=tmp_path / "1a59A.pdb",
        ),
    ]
    for structure in structures:
        structure.pdb_path.write_text("ATOM\n")

    calls = []

    class FakeChainsawRunner:
        def __init__(self, require_cuda=False):
            self.require_cuda = require_cuda

        def run_batch(self, inputs, stage_dir, combined_output):
            calls.append((self.require_cuda, inputs, stage_dir, combined_output))
            rows = {
                "19hcA": "19hcA\tmd5a\t10\t1\t1-10\t0.9\t0.31\n",
                "1a59A": "1a59A\tmd5b\t12\t1\t1-12\t0.8\t0.42\n",
            }
            for _structure_file, output_tsv in inputs:
                output_tsv.parent.mkdir(parents=True, exist_ok=True)
                output_tsv.write_text(
                    "chain_id\tsequence_md5\tnres\tndom\tchopping\tconfidence\ttime_sec\n"
                    + rows[output_tsv.stem]
                )
            return ToolRunResult(["python", "get_predictions.py"], tmp_path, 0, 8.0, 800.0, "", "")

    monkeypatch.setattr(run_benchmark, "ChainsawRunner", FakeChainsawRunner)
    args = SimpleNamespace(results_dir=tmp_path / "results", skip_existing=False)

    batched = run_benchmark._maybe_run_chainsaw_batch(structures, args)

    assert len(calls) == 1
    require_cuda, inputs, _stage_dir, _combined_output = calls[0]
    assert require_cuda is True
    assert len(inputs) == 2
    predictions, result, reused_output, batch_n_entries = batched["19hcA"]
    assert predictions[0].runtime_s == 0.31
    assert result.runtime_s == 8.0
    assert reused_output is False
    assert batch_n_entries == 2


def test_chainsaw_batch_reuses_copied_tsv_from_reference_results_dir(monkeypatch, tmp_path: Path):
    structure = SimpleNamespace(
        entry=CathEntry("19hc", "A", "19hcA", 1, 10, "1-10", "cath663"),
        numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 11)]),
        pdb_path=tmp_path / "19hcA.pdb",
    )
    source_tsv = tmp_path / "reference" / "raw" / "chainsaw" / "19hcA.tsv"
    source_tsv.parent.mkdir(parents=True)
    source_tsv.write_text(
        "chain_id\tsequence_md5\tnres\tndom\tchopping\tconfidence\ttime_sec\n"
        "19hcA\tmd5a\t10\t1\t1-10\t0.9\t0.31\n"
    )

    class UnexpectedChainsawRunner:
        def __init__(self, *_args, **_kwargs):
            raise AssertionError("reference outputs should skip the Chainsaw runner")

    monkeypatch.setattr(run_benchmark, "ChainsawRunner", UnexpectedChainsawRunner)
    args = SimpleNamespace(
        results_dir=tmp_path / "results",
        skip_existing=False,
        reuse_tool_results_dir=tmp_path / "reference",
    )

    batched = run_benchmark._maybe_run_chainsaw_batch([structure], args)

    copied_tsv = tmp_path / "results" / "raw" / "chainsaw" / "19hcA.tsv"
    assert copied_tsv.read_text() == source_tsv.read_text()
    predictions, result, reused_output, batch_n_entries = batched["19hcA"]
    assert predictions[0].raw_path == copied_tsv
    assert predictions[0].runtime_s == 0.31
    assert result is None
    assert reused_output is True
    assert batch_n_entries is None


def test_reused_tool_output_must_exist_before_batch_run(monkeypatch, tmp_path: Path):
    structure = SimpleNamespace(
        entry=CathEntry("19hc", "A", "19hcA", 1, 10, "1-10", "cath663"),
        numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 11)]),
        pdb_path=tmp_path / "19hcA.pdb",
    )

    class UnexpectedChainsawRunner:
        def __init__(self, *_args, **_kwargs):
            raise AssertionError("missing reference outputs should not fall back to the runner")

    monkeypatch.setattr(run_benchmark, "ChainsawRunner", UnexpectedChainsawRunner)
    args = SimpleNamespace(
        results_dir=tmp_path / "results",
        skip_existing=False,
        reuse_tool_results_dir=tmp_path / "reference",
    )

    with pytest.raises(FileNotFoundError, match="Cannot reuse chainsaw output for 19hcA"):
        run_benchmark._maybe_run_chainsaw_batch([structure], args)


def test_chainsaw_batch_keeps_valid_predictions_when_one_entry_returns_null(monkeypatch, tmp_path: Path):
    structures = [
        SimpleNamespace(
            entry=CathEntry("19hc", "A", "19hcA", 1, 10, "1-10", "cath663"),
            numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 11)]),
            pdb_path=tmp_path / "19hcA.pdb",
        ),
        SimpleNamespace(
            entry=CathEntry("1a59", "A", "1a59A", 1, 12, "1-12", "cath663"),
            numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 13)]),
            pdb_path=tmp_path / "1a59A.pdb",
        ),
    ]
    for structure in structures:
        structure.pdb_path.write_text("ATOM\n")

    class FakeChainsawRunner:
        def __init__(self, require_cuda=False):
            self.require_cuda = require_cuda

        def run_batch(self, inputs, stage_dir, combined_output):
            rows = {
                "19hcA": "19hcA\tmd5a\t10\t1\t1-10\t0.9\t0.31\n",
                "1a59A": "1a59A\tmd5b\t12\t0\tNULL\t0.8\t0.42\n",
            }
            for _structure_file, output_tsv in inputs:
                output_tsv.parent.mkdir(parents=True, exist_ok=True)
                output_tsv.write_text(
                    "chain_id\tsequence_md5\tnres\tndom\tchopping\tconfidence\ttime_sec\n"
                    + rows[output_tsv.stem]
                )
            return ToolRunResult(["python", "get_predictions.py"], tmp_path, 0, 8.0, 800.0, "", "")

    monkeypatch.setattr(run_benchmark, "ChainsawRunner", FakeChainsawRunner)
    args = SimpleNamespace(results_dir=tmp_path / "results", skip_existing=False)

    batched = run_benchmark._maybe_run_chainsaw_batch(structures, args)

    valid_predictions, result, reused_output, batch_n_entries = batched["19hcA"]
    assert valid_predictions[0].chopping == "0-9"
    assert result.runtime_s == 8.0
    assert reused_output is False
    assert batch_n_entries == 2
    assert batched["1a59A"].predictions == []
    assert isinstance(batched["1a59A"].error, ValueError)
    assert "NULL" in str(batched["1a59A"].error)


def test_sword2_rust_skip_existing_finds_nested_summary_without_rerun(monkeypatch, tmp_path: Path):
    structure = SimpleNamespace(
        entry=CathEntry("19hc", "A", "19hcA", 1, 10, "1-10", "cath663"),
        numbering=StructureNumbering([ResidueKey("A", str(index), "") for index in range(1, 11)]),
        pdb_path=tmp_path / "19hcA.pdb",
    )
    nested_summary = tmp_path / "results" / "raw" / "sword2-rust" / "19hcA" / "19hcA_A" / "summary.json"
    nested_summary.parent.mkdir(parents=True)
    nested_summary.write_text("{}")

    class UnexpectedSword2Runner:
        def __init__(self, threads=None):
            self.threads = threads

        def run(self, *_args, **_kwargs):
            raise AssertionError("skip-existing should reuse the nested summary")

    loaded_paths: list[Path] = []

    def fake_load_summary_partitions(summary_path, **_kwargs):
        loaded_paths.append(summary_path)
        return [
            PartitionPrediction(
                tool="sword2-rust",
                name="Optimal partition",
                variant="optimal",
                chopping="0-9",
                n_domains=1,
                raw_path=summary_path,
            )
        ]

    monkeypatch.setattr(run_benchmark, "Sword2RustRunner", UnexpectedSword2Runner)
    monkeypatch.setattr(run_benchmark, "load_summary_partitions", fake_load_summary_partitions)
    args = SimpleNamespace(results_dir=tmp_path / "results", skip_existing=True, sword2_threads=None)

    predictions, result = run_benchmark._maybe_run_sword2_rust(structure, args)

    assert predictions[0].raw_path == nested_summary
    assert loaded_paths == [nested_summary]
    assert result is None


def test_summarize_existing_writes_comparison_analysis_outputs(monkeypatch, tmp_path: Path):
    scores_csv = tmp_path / "scores.csv"
    pd.DataFrame(
        [
            {
                "entry_id": "19hcA",
                "tool": "chainsaw",
                "variant": "single",
                "partition": "Chainsaw",
                "ndo": 0.5,
                "boundary_dist_score": 0.5,
                "boundary_f1_10": 0.5,
                "matched_dice": 0.5,
                "pairwise_f1": 0.5,
                "d_count_acc": 1.0,
                "iou": 0.5,
                "runtime_s": 0.1,
                "peak_rss_mb": 10.0,
                "n_residues": 100,
            }
        ]
    ).to_csv(scores_csv, index=False)

    calls: list[str] = []

    def fake_writer(name):
        def inner(*args, **kwargs):
            calls.append(name)
            return []

        return inner

    monkeypatch.setattr(run_benchmark, "plot_metric_distributions", fake_writer("metric_plots"))
    monkeypatch.setattr(run_benchmark, "plot_runtime_memory", fake_writer("runtime_plots"))
    monkeypatch.setattr(run_benchmark, "plot_paired_deltas", fake_writer("paired_deltas"))
    monkeypatch.setattr(run_benchmark, "write_metric_summaries", fake_writer("metric_summaries"))
    monkeypatch.setattr(run_benchmark, "write_win_rates", fake_writer("win_rates"))
    monkeypatch.setattr(
        run_benchmark,
        "plot_sword2_oracle_gap",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("SWORD2-only plot should not run")),
        raising=False,
    )
    monkeypatch.setattr(
        run_benchmark,
        "plot_sword2_topk_curve",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("SWORD2-only plot should not run")),
        raising=False,
    )
    monkeypatch.setattr(
        run_benchmark,
        "plot_alternatives_per_entry",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("SWORD2-only plot should not run")),
        raising=False,
    )

    run_benchmark.summarize_existing(scores_csv, tmp_path)

    assert calls == [
        "metric_plots",
        "runtime_plots",
        "paired_deltas",
        "metric_summaries",
        "win_rates",
    ]


def test_summarize_existing_writes_sword2_alternative_analysis_outputs(monkeypatch, tmp_path: Path):
    scores_csv = tmp_path / "scores.csv"
    pd.DataFrame(
        [
            {
                "entry_id": "19hcA",
                "tool": "sword2-rust",
                "variant": "optimal",
                "partition": "Optimal partition",
                "ndo": 0.6,
                "boundary_f1_10": 0.5,
                "matched_dice": 0.7,
                "pairwise_f1": 0.8,
                "runtime_s": 0.1,
                "peak_rss_mb": 10.0,
                "n_residues": 100,
            },
            {
                "entry_id": "19hcA",
                "tool": "sword2-rust",
                "variant": "alternative",
                "partition": "Alternative partition 1",
                "ndo": 0.8,
                "boundary_f1_10": 0.7,
                "matched_dice": 0.8,
                "pairwise_f1": 0.9,
                "runtime_s": 0.1,
                "peak_rss_mb": 10.0,
                "n_residues": 100,
            },
            {
                "entry_id": "19hcA",
                "tool": "sword2-rust",
                "variant": "oracle",
                "partition": "Oracle partition",
                "ndo": 0.8,
                "boundary_f1_10": 0.7,
                "matched_dice": 0.8,
                "pairwise_f1": 0.9,
                "runtime_s": 0.1,
                "peak_rss_mb": 10.0,
                "n_residues": 100,
            },
        ]
    ).to_csv(scores_csv, index=False)

    calls: list[str] = []

    def fake_writer(name):
        def inner(*args, **kwargs):
            calls.append(name)
            return []

        return inner

    monkeypatch.setattr(run_benchmark, "plot_metric_distributions", fake_writer("metric_plots"))
    monkeypatch.setattr(run_benchmark, "plot_runtime_memory", fake_writer("runtime_plots"))
    monkeypatch.setattr(run_benchmark, "plot_paired_deltas", fake_writer("paired_deltas"))
    monkeypatch.setattr(run_benchmark, "plot_sword2_oracle_gap", fake_writer("sword2_oracle_gap"))
    monkeypatch.setattr(run_benchmark, "plot_sword2_topk_curve", fake_writer("sword2_topk"))
    monkeypatch.setattr(run_benchmark, "plot_alternatives_per_entry", fake_writer("sword2_alternative_counts"))
    monkeypatch.setattr(run_benchmark, "write_metric_summaries", fake_writer("metric_summaries"))
    monkeypatch.setattr(run_benchmark, "write_win_rates", fake_writer("win_rates"))

    run_benchmark.summarize_existing(scores_csv, tmp_path)

    assert calls == [
        "metric_plots",
        "runtime_plots",
        "paired_deltas",
        "sword2_oracle_gap",
        "sword2_topk",
        "sword2_alternative_counts",
        "metric_summaries",
        "win_rates",
    ]


def test_summarize_existing_removes_stale_figure_files(monkeypatch, tmp_path: Path):
    scores_csv = tmp_path / "scores.csv"
    pd.DataFrame(
        [
            {
                "entry_id": "19hcA",
                "tool": "chainsaw",
                "variant": "single",
                "partition": "Chainsaw",
                "ndo": 0.5,
                "boundary_dist_score": 0.5,
                "boundary_f1_10": 0.5,
                "matched_dice": 0.5,
                "pairwise_f1": 0.5,
                "d_count_acc": 1.0,
                "iou": 0.5,
                "runtime_s": 0.1,
                "peak_rss_mb": 10.0,
                "n_residues": 100,
            }
        ]
    ).to_csv(scores_csv, index=False)
    figure_dir = tmp_path / "figures"
    figure_dir.mkdir()
    stale_svg = figure_dir / "stale.svg"
    stale_png = figure_dir / "stale.png"
    stale_svg.write_text("<svg></svg>")
    stale_png.write_text("png")

    monkeypatch.setattr(run_benchmark, "plot_metric_distributions", lambda *args, **kwargs: [])
    monkeypatch.setattr(run_benchmark, "plot_runtime_memory", lambda *args, **kwargs: [])
    monkeypatch.setattr(run_benchmark, "plot_paired_deltas", lambda *args, **kwargs: [])
    monkeypatch.setattr(run_benchmark, "write_metric_summaries", lambda *args, **kwargs: [])
    monkeypatch.setattr(run_benchmark, "write_win_rates", lambda *args, **kwargs: None)

    run_benchmark.summarize_existing(scores_csv, tmp_path)

    assert not stale_svg.exists()
    assert not stale_png.exists()


def test_progress_messages_include_position_entry_and_tool():
    assert (
        run_benchmark._progress_message(
            current=3,
            total=12,
            entry_id="1abcA",
            tool="sword2-rust",
        )
        == "[3/12] 1abcA sword2-rust"
    )

    assert (
        run_benchmark._progress_message(
            current=2,
            total=663,
            entry_id="2xyzB",
            tool=None,
        )
        == "[2/663] 2xyzB"
    )


def test_figures_are_saved_as_png_only(monkeypatch, tmp_path: Path):
    saved_suffixes: list[str] = []
    monkeypatch.setattr(figures.plt, "tight_layout", lambda: None)
    monkeypatch.setattr(figures.plt, "savefig", lambda path: saved_suffixes.append(Path(path).suffix))
    monkeypatch.setattr(figures.plt, "close", lambda: None)

    paths = figures._savefig(tmp_path / "plot.png")

    assert paths == [tmp_path / "plot.png"]
    assert saved_suffixes == [".png"]


def _stub_figure_io(monkeypatch):
    monkeypatch.setattr(figures, "_setup_theme", lambda: None)
    monkeypatch.setattr(figures.plt, "figure", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        figures.plt,
        "gca",
        lambda: SimpleNamespace(get_legend_handles_labels=lambda: ([], [])),
    )
    monkeypatch.setattr(figures.plt, "xlabel", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "ylabel", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "title", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "xticks", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "legend", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "axhline", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures, "_savefig", lambda path: [path])


@pytest.mark.parametrize(
    ("plotter", "scores", "expected_y"),
    [
        (
            figures.plot_metric_distributions,
            pd.DataFrame(
                [
                    {"entry_id": "a", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.8},
                    {"entry_id": "b", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.6},
                    {"entry_id": "a", "tool": "merizo", "variant": "cpu", "ndo": 0.7},
                    {"entry_id": "b", "tool": "merizo", "variant": "cpu", "ndo": 0.9},
                ]
            ),
            "ndo",
        ),
        (
            figures.plot_sword2_oracle_gap,
            pd.DataFrame(
                [
                    {"entry_id": "a", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.6},
                    {"entry_id": "a", "tool": "sword2-rust", "variant": "oracle", "ndo": 0.8},
                    {"entry_id": "b", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.5},
                    {"entry_id": "b", "tool": "sword2-rust", "variant": "oracle", "ndo": 0.9},
                ]
            ),
            "gap",
        ),
        (
            figures.plot_alternatives_per_entry,
            pd.DataFrame(
                [
                    {"entry_id": "a", "tool": "sword2-rust", "variant": "optimal"},
                    {"entry_id": "a", "tool": "sword2-rust", "variant": "alternative"},
                    {"entry_id": "b", "tool": "sword2-rust", "variant": "optimal"},
                    {"entry_id": "b", "tool": "sword2-rust", "variant": "alternative"},
                    {"entry_id": "b", "tool": "sword2-rust", "variant": "alternative"},
                ]
            ),
            "n_alternatives",
        ),
        (
            figures.plot_paired_deltas,
            pd.DataFrame(
                [
                    {"entry_id": "a", "tool": "merizo", "variant": "cpu", "ndo": 0.7},
                    {"entry_id": "a", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.9},
                    {"entry_id": "b", "tool": "merizo", "variant": "cpu", "ndo": 0.8},
                    {"entry_id": "b", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.6},
                ]
            ),
            "delta",
        ),
    ],
)
def test_distribution_figures_use_violin_plots(monkeypatch, tmp_path: Path, plotter, scores, expected_y):
    box_calls: list[dict[str, object]] = []
    violin_calls: list[dict[str, object]] = []
    strip_calls: list[dict[str, object]] = []

    _stub_figure_io(monkeypatch)
    monkeypatch.setattr(figures.sns, "boxplot", lambda **kwargs: box_calls.append(kwargs))
    monkeypatch.setattr(figures.sns, "violinplot", lambda **kwargs: violin_calls.append(kwargs))
    monkeypatch.setattr(figures.sns, "stripplot", lambda **kwargs: strip_calls.append(kwargs))

    plotter(scores, tmp_path)

    assert box_calls == []
    assert violin_calls
    assert violin_calls[0]["y"] == expected_y
    assert violin_calls[0]["inner"] == "quartile"
    assert violin_calls[0]["cut"] == 0
    assert strip_calls


def test_sword2_topk_curve_keeps_fixed_entry_cohort(monkeypatch, tmp_path: Path):
    scores = pd.DataFrame(
        [
            {
                "entry_id": "easy",
                "tool": "sword2-rust",
                "variant": "optimal",
                "partition": "Optimal partition",
                "ndo": 0.9,
            },
            {
                "entry_id": "hard",
                "tool": "sword2-rust",
                "variant": "optimal",
                "partition": "Optimal partition",
                "ndo": 0.1,
            },
            {
                "entry_id": "hard",
                "tool": "sword2-rust",
                "variant": "alternative",
                "partition": "Alternative partition 1",
                "ndo": 0.2,
            },
            {
                "entry_id": "hard",
                "tool": "sword2-rust",
                "variant": "alternative",
                "partition": "Alternative partition 2",
                "ndo": 0.3,
            },
        ]
    )
    captured: dict[str, pd.DataFrame] = {}

    def fake_lineplot(**kwargs):
        captured["data"] = kwargs["data"].copy()

    monkeypatch.setattr(figures, "_setup_theme", lambda: None)
    monkeypatch.setattr(figures.plt, "figure", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.sns, "lineplot", fake_lineplot)
    monkeypatch.setattr(figures.plt, "xlabel", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "ylabel", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "title", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "legend", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures, "_savefig", lambda path: [path])

    figures.plot_sword2_topk_curve(scores, tmp_path)

    curve = captured["data"].sort_values("k")
    assert curve["k"].tolist() == [1, 2, 3]
    assert curve["best_ndo"].tolist() == pytest.approx([0.5, 0.55, 0.6])


def test_run_level_data_falls_back_to_score_runtime_for_reused_outputs():
    scores = pd.DataFrame(
        [
            {
                "entry_id": "1abcA",
                "tool": "sword2-rust",
                "variant": "optimal",
                "n_residues": 100,
                "runtime_s": 1.0,
                "peak_rss_mb": 100.0,
            },
            {
                "entry_id": "1abcA",
                "tool": "merizo",
                "variant": "cpu",
                "n_residues": 100,
                "runtime_s": 2.0,
                "peak_rss_mb": 200.0,
            },
        ]
    )
    runs = pd.DataFrame(
        [
            {
                "entry_id": "1abcA",
                "tool": "sword2-rust",
                "reused_output": False,
                "runtime_s": 1.5,
                "peak_rss_mb": 150.0,
            },
            {
                "entry_id": "1abcA",
                "tool": "merizo",
                "reused_output": True,
                "runtime_s": None,
                "peak_rss_mb": None,
            },
        ]
    )

    data = figures._run_level_data(scores, runs=runs)

    runtimes = data.set_index("tool")["runtime_s"].to_dict()
    assert runtimes == {"sword2-rust": 1.5, "merizo": 2.0}


def test_runtime_plots_exclude_tools_without_metric_values(monkeypatch, tmp_path: Path):
    scores = pd.DataFrame(
        [
            {
                "entry_id": "1abcA",
                "tool": "sword2-rust",
                "variant": "optimal",
                "n_residues": 100,
                "runtime_s": 1.0,
            },
            {
                "entry_id": "1abcA",
                "tool": "merizo",
                "variant": "cpu",
                "n_residues": 100,
                "runtime_s": None,
            },
        ]
    )
    captured: dict[str, object] = {}

    def fake_scatterplot(**kwargs):
        captured["tools"] = set(kwargs["data"]["tool"])
        captured["hue_order"] = kwargs["hue_order"]

    monkeypatch.setattr(figures, "_setup_theme", lambda: None)
    monkeypatch.setattr(figures.plt, "figure", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.sns, "scatterplot", fake_scatterplot)
    monkeypatch.setattr(figures.plt, "yscale", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "xlabel", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "ylabel", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "title", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures.plt, "legend", lambda *args, **kwargs: None)
    monkeypatch.setattr(figures, "_savefig", lambda path: [path])

    figures.plot_runtime_memory(scores, tmp_path)

    assert captured["tools"] == {"sword2-rust"}
    assert captured["hue_order"] == ["sword2-rust"]


def test_missing_score_runtime_is_filled_from_previous_scores():
    row = _score_row(runtime_s=None, peak_rss_mb=None)
    lookup = {
        ("cath663", "19hcA", "merizo", "cpu", "Merizo"): (2.5, 256.0),
    }

    filled = run_benchmark._fill_missing_score_runtimes([row], lookup)

    assert filled[0].runtime_s == 2.5
    assert filled[0].peak_rss_mb == 256.0
