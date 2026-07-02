from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark.datasets import dataset_path, load_dataset
from benchmark.figures import (
    plot_alternatives_per_entry,
    plot_metric_distributions,
    plot_paired_deltas,
    plot_runtime_memory,
    plot_sword2_oracle_gap,
    plot_sword2_topk_curve,
)
from benchmark.summaries import write_metric_summaries, write_win_rates
from benchmark.runners.base import PartitionPrediction, ToolRunResult
from benchmark.runners.chainsaw import ChainsawRunner, chainsaw_to_common_chopping, parse_chainsaw_tsv
from benchmark.runners.merizo import MerizoRunner, merizo_to_common_chopping, parse_merizo_tsv
from benchmark.runners.sword2_orig import Sword2OriginalRunner, load_original_summary_partitions
from benchmark.runners.sword2_rust import Sword2RustRunner, load_summary_partitions
from benchmark.score import (
    FailureRow,
    RunRow,
    ScoreRow,
    add_sword2_oracle,
    add_sword2_oracle_s,
    score_prediction,
    write_failures_csv,
    write_runs_csv,
    write_scores_csv,
)
from benchmark.structures import canonicalize_entry


@dataclass(frozen=True)
class BatchedToolResult:
    predictions: list[PartitionPrediction]
    result: ToolRunResult | None
    reused_output: bool
    batch_n_entries: int | None
    error: Exception | None = None

    def __iter__(self):
        yield self.predictions
        yield self.result
        yield self.reused_output
        yield self.batch_n_entries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the SWORD2 benchmark harness")
    parser.add_argument("--dataset", default="cath663", help="Dataset key, e.g. cath663")
    parser.add_argument("--limit", type=int, default=None, help="Limit entries for smoke runs")
    parser.add_argument("--cache-dir", type=Path, default=Path("benchmark/cache"))
    parser.add_argument("--results-dir", type=Path, default=Path("benchmark/results"))
    parser.add_argument("--no-download", action="store_true", help="Use only cached raw PDB files")
    parser.add_argument("--skip-existing", action="store_true", help="Reuse existing raw tool outputs")
    parser.add_argument(
        "--reuse-tool-results-dir",
        type=Path,
        default=None,
        help=(
            "Copy reusable Merizo/Chainsaw raw TSV outputs from an existing benchmark "
            "results directory or raw directory instead of recomputing them"
        ),
    )
    parser.add_argument("--strict", action="store_true", help="Fail on first tool error instead of skipping")
    parser.add_argument(
        "--tools",
        default="sword2-rust,merizo,chainsaw,sword2-original",
        help="Comma-separated tools to run",
    )
    parser.add_argument("--merizo-device", default="cuda", choices=["cpu", "cuda", "mps"])
    parser.add_argument(
        "--allow-dl-cpu",
        action="store_true",
        help="Allow DL tools to run when CUDA is unavailable instead of failing before launch",
    )
    parser.add_argument("--sword2-threads", type=int, default=None)
    parser.add_argument(
        "--sword2-experiments",
        default=None,
        help="Comma-separated SWORD2_EXPERIMENTS flags for sword2-rust benchmark runs",
    )
    parser.add_argument(
        "--sword2-extra-args",
        default=None,
        help="Extra CLI args passed verbatim to the sword2-rust binary, e.g. '--use-pairwise-reranker'",
    )
    parser.add_argument(
        "--sword2-original-conda-env",
        type=Path,
        default=Path("/opt/apps/pkgs/anaconda3/2023.09-0/kvhhznl/envs/sword2"),
        help="Conda environment prefix for original SWORD2",
    )
    parser.add_argument(
        "--score-csv",
        type=Path,
        default=None,
        help="Existing scores CSV to summarize without running tools",
    )
    return parser.parse_args()


def _load_runs_for_summary(scores_csv: Path, results_dir: Path) -> pd.DataFrame | None:
    for path in [results_dir / "runs.csv", scores_csv.parent / "runs.csv"]:
        if path.exists():
            return pd.read_csv(path)
    return None


def _clear_figure_dir(figure_dir: Path) -> None:
    if not figure_dir.exists():
        return
    for pattern in ["*.png", "*.svg"]:
        for path in figure_dir.glob(pattern):
            path.unlink()


def _sword2_variants(scores: pd.DataFrame) -> set[str]:
    sword2 = scores[scores["tool"].str.startswith("sword2", na=False)]
    return set(sword2["variant"].dropna())


def summarize_existing(scores_csv: Path, results_dir: Path) -> None:
    scores = pd.read_csv(scores_csv)
    runs = _load_runs_for_summary(scores_csv, results_dir)
    figure_dir = results_dir / "figures"
    table_dir = results_dir / "tables"
    _clear_figure_dir(figure_dir)
    plot_metric_distributions(scores, figure_dir)
    plot_runtime_memory(scores, figure_dir, runs=runs)
    plot_paired_deltas(scores, figure_dir)
    sword2_variants = _sword2_variants(scores)
    if {"optimal", "oracle"}.issubset(sword2_variants):
        plot_sword2_oracle_gap(scores, figure_dir)
    if {"optimal", "alternative"}.issubset(sword2_variants):
        plot_sword2_topk_curve(scores, figure_dir)
        plot_alternatives_per_entry(scores, figure_dir)
    write_metric_summaries(scores, table_dir)
    write_win_rates(scores, table_dir)
    print(f"Wrote figures to {figure_dir}")
    print(f"Wrote tables to {table_dir}")


def _is_missing(value: object) -> bool:
    return value is None or bool(pd.isna(value))


def _score_runtime_key(row: ScoreRow | dict[str, object]) -> tuple[str, str, str, str, str]:
    if isinstance(row, ScoreRow):
        return (row.dataset, row.entry_id, row.tool, row.variant, row.partition)
    return (
        str(row["dataset"]),
        str(row["entry_id"]),
        str(row["tool"]),
        str(row["variant"]),
        str(row["partition"]),
    )


def _score_runtime_lookup(scores_csv: Path) -> dict[tuple[str, str, str, str, str], tuple[float | None, float | None]]:
    if not scores_csv.exists():
        return {}
    scores = pd.read_csv(scores_csv)
    lookup: dict[tuple[str, str, str, str, str], tuple[float | None, float | None]] = {}
    for row in scores.to_dict("records"):
        runtime = row.get("runtime_s")
        peak_rss = row.get("peak_rss_mb")
        if _is_missing(runtime) and _is_missing(peak_rss):
            continue
        lookup[_score_runtime_key(row)] = (
            None if _is_missing(runtime) else float(runtime),
            None if _is_missing(peak_rss) else float(peak_rss),
        )
    return lookup


def _fill_missing_score_runtimes(
    rows: list[ScoreRow],
    runtime_lookup: dict[tuple[str, str, str, str, str], tuple[float | None, float | None]],
) -> list[ScoreRow]:
    filled: list[ScoreRow] = []
    for row in rows:
        historical = runtime_lookup.get(_score_runtime_key(row))
        if historical is None:
            filled.append(row)
            continue
        runtime_s = row.runtime_s
        peak_rss_mb = row.peak_rss_mb
        if _is_missing(runtime_s):
            runtime_s = historical[0]
        if _is_missing(peak_rss_mb):
            peak_rss_mb = historical[1]
        filled.append(replace(row, runtime_s=runtime_s, peak_rss_mb=peak_rss_mb))
    return filled


def _find_one(root: Path, name: str) -> Path:
    matches = sorted(root.rglob(name))
    if not matches:
        raise FileNotFoundError(f"Could not find {name} under {root}")
    return matches[0]


def _check_run(result: ToolRunResult, tool: str) -> None:
    if result.returncode != 0:
        raise RuntimeError(f"{tool} failed with exit {result.returncode}\n{result.stderr}")


def _check_run_or_summary(result: ToolRunResult, tool: str, raw_dir: Path, summary_name: str) -> Path:
    matches = sorted(raw_dir.rglob(summary_name))
    if result.returncode == 0:
        if not matches:
            raise FileNotFoundError(f"Could not find {summary_name} under {raw_dir}")
        return matches[0]
    if matches:
        print(
            f"[{tool}] exited {result.returncode} after writing {matches[0]}; continuing with summary.",
            file=sys.stderr,
        )
        return matches[0]
    raise RuntimeError(f"{tool} failed with exit {result.returncode}\n{result.stderr}")


def _require_dl_cuda(args: argparse.Namespace) -> bool:
    return not getattr(args, "allow_dl_cpu", False)


def _format_command(command: list[str]) -> str:
    return " ".join(str(part) for part in command)


def _progress_message(current: int, total: int, entry_id: str, tool: str | None = None) -> str:
    parts = [f"[{current}/{total}]", entry_id]
    if tool:
        parts.append(tool)
    return " ".join(parts)


def _print_progress(current: int, total: int, entry_id: str, tool: str | None = None) -> None:
    print(_progress_message(current, total, entry_id, tool), file=sys.stderr, flush=True)


def _run_row_from_result(
    structure,
    tool: str,
    variant: str,
    raw_path: Path | None,
    result: ToolRunResult | None,
    reused_output: bool,
    runtime_s: float | None = None,
    peak_rss_mb: float | None = None,
    use_result_runtime: bool = True,
    use_result_peak: bool = True,
    batch_runtime_s: float | None = None,
    batch_peak_rss_mb: float | None = None,
    batch_n_entries: int | None = None,
) -> RunRow:
    row_runtime_s = result.runtime_s if use_result_runtime and result else runtime_s
    row_peak_rss_mb = result.peak_rss_mb if use_result_peak and result else peak_rss_mb
    return RunRow(
        dataset=structure.entry.dataset,
        entry_id=structure.entry.entry_id,
        pdb_id=structure.entry.pdb_id,
        chain_id=structure.entry.chain_id,
        tool=tool,
        variant=variant,
        reused_output=reused_output,
        returncode=result.returncode if result else None,
        runtime_s=row_runtime_s,
        peak_rss_mb=row_peak_rss_mb,
        command=_format_command(result.command) if result else "",
        cwd=str(result.cwd) if result else "",
        raw_path=str(raw_path) if raw_path else "",
        batch_runtime_s=batch_runtime_s,
        batch_peak_rss_mb=batch_peak_rss_mb,
        batch_n_entries=batch_n_entries,
    )


def _failure_row_from_exception(
    structure,
    tool: str,
    exc: Exception,
    stage: str = "run",
    result: ToolRunResult | None = None,
) -> FailureRow:
    return FailureRow(
        dataset=structure.entry.dataset,
        entry_id=structure.entry.entry_id,
        pdb_id=structure.entry.pdb_id,
        chain_id=structure.entry.chain_id,
        tool=tool,
        stage=stage,
        message=str(exc),
        returncode=result.returncode if result else None,
        command=_format_command(result.command) if result else None,
        stderr=result.stderr if result else None,
    )


def _run_variant(tool: str, predictions: list[PartitionPrediction]) -> str:
    if tool.startswith("sword2"):
        return "all-partitions"
    if predictions:
        return predictions[0].variant
    return "single"


def _raw_path(predictions: list[PartitionPrediction]) -> Path | None:
    for prediction in predictions:
        if prediction.raw_path is not None:
            return prediction.raw_path
    return None


def _parsed_runtime(parsed, column: str) -> float | None:
    if parsed.raw is None:
        return None
    value = parsed.raw.get(column)
    if _is_missing(value) or value == "":
        return None
    return float(value)


def _prediction_runtime_s(prediction: PartitionPrediction, result: ToolRunResult | None) -> float | None:
    if prediction.runtime_s is not None:
        return prediction.runtime_s
    return result.runtime_s if result else None


def _prediction_peak_rss_mb(
    prediction: PartitionPrediction,
    result: ToolRunResult | None,
    batched: bool = False,
) -> float | None:
    if batched and prediction.runtime_s is not None:
        return None
    return result.peak_rss_mb if result else None


def _rows_from_predictions(
    structure,
    tool: str,
    predictions: list[PartitionPrediction],
    result: ToolRunResult | None,
    reused_output: bool,
    batch_n_entries: int | None = None,
) -> tuple[list[ScoreRow], RunRow]:
    score_rows = [
        score_prediction(
            structure.entry,
            structure.numbering,
            prediction,
            runtime_s=_prediction_runtime_s(prediction, result),
            peak_rss_mb=_prediction_peak_rss_mb(prediction, result, batched=batch_n_entries is not None),
        )
        for prediction in predictions
    ]
    first_prediction = predictions[0] if predictions else None
    runtime_s = _prediction_runtime_s(first_prediction, result) if first_prediction else None
    peak_rss_mb = (
        _prediction_peak_rss_mb(first_prediction, result, batched=batch_n_entries is not None)
        if first_prediction
        else None
    )
    run_row = _run_row_from_result(
        structure,
        tool=tool,
        variant=_run_variant(tool, predictions),
        raw_path=_raw_path(predictions),
        result=result,
        reused_output=reused_output,
        runtime_s=runtime_s,
        peak_rss_mb=peak_rss_mb,
        use_result_runtime=False,
        use_result_peak=False,
        batch_runtime_s=result.runtime_s if batch_n_entries is not None and result else None,
        batch_peak_rss_mb=result.peak_rss_mb if batch_n_entries is not None and result else None,
        batch_n_entries=batch_n_entries,
    )
    return score_rows, run_row


def _merizo_tsv_path(structure, args: argparse.Namespace) -> Path:
    return args.results_dir / "raw" / f"merizo-{args.merizo_device}" / f"{structure.entry.entry_id}.tsv"


def _reuse_raw_root(args: argparse.Namespace) -> Path | None:
    reuse_dir = getattr(args, "reuse_tool_results_dir", None)
    if reuse_dir is None:
        return None
    raw_dir = reuse_dir / "raw"
    return raw_dir if raw_dir.exists() else reuse_dir


def _copy_reused_output(tool: str, entry_id: str, source_path: Path | None, target_path: Path) -> bool:
    if source_path is None:
        return False
    if not source_path.exists():
        raise FileNotFoundError(f"Cannot reuse {tool} output for {entry_id}: {source_path} does not exist")
    target_path.parent.mkdir(parents=True, exist_ok=True)
    if source_path.resolve() != target_path.resolve():
        shutil.copy2(source_path, target_path)
    return True


def _copy_reused_merizo_output(structure, args: argparse.Namespace) -> bool:
    raw_root = _reuse_raw_root(args)
    source_path = (
        raw_root / f"merizo-{args.merizo_device}" / f"{structure.entry.entry_id}.tsv"
        if raw_root is not None
        else None
    )
    return _copy_reused_output("merizo", structure.entry.entry_id, source_path, _merizo_tsv_path(structure, args))


def _merizo_predictions_from_tsv(structure, args: argparse.Namespace, tsv_path: Path) -> list[PartitionPrediction]:
    parsed = parse_merizo_tsv(tsv_path.read_text())
    chopping = merizo_to_common_chopping(parsed.chopping, structure.numbering, structure.entry.chain_id)
    prediction = PartitionPrediction(
        tool="merizo",
        name="Merizo",
        variant=args.merizo_device,
        chopping=chopping,
        n_domains=parsed.n_domains or len(chopping.split(",")),
        raw_path=tsv_path,
        runtime_s=_parsed_runtime(parsed, "runtime"),
    )
    return [prediction]


def _chainsaw_tsv_path(structure, args: argparse.Namespace) -> Path:
    return args.results_dir / "raw" / "chainsaw" / f"{structure.entry.entry_id}.tsv"


def _copy_reused_chainsaw_output(structure, args: argparse.Namespace) -> bool:
    raw_root = _reuse_raw_root(args)
    source_path = raw_root / "chainsaw" / f"{structure.entry.entry_id}.tsv" if raw_root is not None else None
    return _copy_reused_output("chainsaw", structure.entry.entry_id, source_path, _chainsaw_tsv_path(structure, args))


def _chainsaw_predictions_from_tsv(structure, tsv_path: Path) -> list[PartitionPrediction]:
    parsed = parse_chainsaw_tsv(tsv_path.read_text())
    chopping = chainsaw_to_common_chopping(
        parsed.chopping,
        numbering=structure.numbering,
        chain_id=structure.entry.chain_id,
    )
    prediction = PartitionPrediction(
        tool="chainsaw",
        name="Chainsaw",
        variant="single",
        chopping=chopping,
        n_domains=parsed.n_domains or len(chopping.split(",")),
        raw_path=tsv_path,
        runtime_s=_parsed_runtime(parsed, "time_sec"),
    )
    return [prediction]


def _maybe_run_merizo_batch(
    structures,
    args: argparse.Namespace,
) -> dict[str, BatchedToolResult]:
    missing_by_chain: dict[str, list] = {}
    result_by_entry: dict[str, tuple[ToolRunResult, int]] = {}
    for structure in structures:
        tsv_path = _merizo_tsv_path(structure, args)
        if _copy_reused_merizo_output(structure, args):
            continue
        if args.skip_existing and tsv_path.exists():
            continue
        missing_by_chain.setdefault(structure.entry.chain_id, []).append(structure)

    if missing_by_chain:
        runner = MerizoRunner(device=args.merizo_device, require_cuda=_require_dl_cuda(args))
        for batch_structures in missing_by_chain.values():
            result = runner.run_batch(
                [
                    (structure.pdb_path, _merizo_tsv_path(structure, args), structure.entry.chain_id)
                    for structure in batch_structures
                ]
            )
            _check_run(result, f"merizo-{args.merizo_device}")
            for structure in batch_structures:
                result_by_entry[structure.entry.entry_id] = (result, len(batch_structures))

    batch_results: dict[str, BatchedToolResult] = {}
    for structure in structures:
        result_info = result_by_entry.get(structure.entry.entry_id)
        result = result_info[0] if result_info else None
        batch_n_entries = result_info[1] if result_info else None
        try:
            predictions = _merizo_predictions_from_tsv(structure, args, _merizo_tsv_path(structure, args))
            error = None
        except Exception as exc:  # noqa: BLE001 - preserve other batch entries.
            predictions = []
            error = exc
        batch_results[structure.entry.entry_id] = BatchedToolResult(
            predictions,
            result,
            reused_output=result is None,
            batch_n_entries=batch_n_entries,
            error=error,
        )
    return batch_results


def _maybe_run_chainsaw_batch(
    structures,
    args: argparse.Namespace,
) -> dict[str, BatchedToolResult]:
    missing_structures = []
    for structure in structures:
        if _copy_reused_chainsaw_output(structure, args):
            continue
        if args.skip_existing and _chainsaw_tsv_path(structure, args).exists():
            continue
        missing_structures.append(structure)
    result_by_entry: dict[str, tuple[ToolRunResult, int]] = {}
    if missing_structures:
        raw_dir = args.results_dir / "raw" / "chainsaw"
        result = ChainsawRunner(require_cuda=_require_dl_cuda(args)).run_batch(
            [(structure.pdb_path, _chainsaw_tsv_path(structure, args)) for structure in missing_structures],
            stage_dir=raw_dir / "_batch_stage",
            combined_output=raw_dir / "_batch.tsv",
        )
        _check_run(result, "chainsaw")
        for structure in missing_structures:
            result_by_entry[structure.entry.entry_id] = (result, len(missing_structures))

    batch_results: dict[str, BatchedToolResult] = {}
    for structure in structures:
        result_info = result_by_entry.get(structure.entry.entry_id)
        result = result_info[0] if result_info else None
        batch_n_entries = result_info[1] if result_info else None
        try:
            predictions = _chainsaw_predictions_from_tsv(structure, _chainsaw_tsv_path(structure, args))
            error = None
        except Exception as exc:  # noqa: BLE001 - preserve other batch entries.
            predictions = []
            error = exc
        batch_results[structure.entry.entry_id] = BatchedToolResult(
            predictions,
            result,
            reused_output=result is None,
            batch_n_entries=batch_n_entries,
            error=error,
        )
    return batch_results


def _maybe_run_sword2_rust(structure, args: argparse.Namespace) -> tuple[list[PartitionPrediction], ToolRunResult | None]:
    raw_dir = args.results_dir / "raw" / "sword2-rust" / structure.entry.entry_id
    summary_path = raw_dir / "summary.json"
    if args.skip_existing and not summary_path.exists():
        try:
            summary_path = _find_one(raw_dir, "summary.json")
        except FileNotFoundError:
            pass
    if not args.skip_existing or not summary_path.exists():
        result = Sword2RustRunner(
            threads=args.sword2_threads,
            experiments=args.sword2_experiments,
            extra_args=args.sword2_extra_args,
        ).run(structure.pdb_path, raw_dir)
        _check_run(result, "sword2-rust")
    else:
        result = None
    if not summary_path.exists():
        summary_path = _find_one(raw_dir, "summary.json")
    mapping_path = summary_path.parent / "residue_mapping.txt"
    return (
        load_summary_partitions(
            summary_path,
            numbering=structure.numbering,
            residue_mapping_path=mapping_path,
            chain_id=structure.entry.chain_id,
        ),
        result,
    )


def _maybe_run_merizo(structure, args: argparse.Namespace) -> tuple[list[PartitionPrediction], ToolRunResult | None]:
    tsv_path = _merizo_tsv_path(structure, args)
    if _copy_reused_merizo_output(structure, args):
        result = None
    elif not args.skip_existing or not tsv_path.exists():
        result = MerizoRunner(device=args.merizo_device, require_cuda=_require_dl_cuda(args)).run(
            structure.pdb_path,
            tsv_path,
            chain_id=structure.entry.chain_id,
        )
        _check_run(result, f"merizo-{args.merizo_device}")
    else:
        result = None
    return _merizo_predictions_from_tsv(structure, args, tsv_path), result


def _maybe_run_chainsaw(structure, args: argparse.Namespace) -> tuple[list[PartitionPrediction], ToolRunResult | None]:
    tsv_path = _chainsaw_tsv_path(structure, args)
    if _copy_reused_chainsaw_output(structure, args):
        result = None
    elif not args.skip_existing or not tsv_path.exists():
        result = ChainsawRunner(require_cuda=_require_dl_cuda(args)).run(structure.pdb_path, tsv_path)
        _check_run(result, "chainsaw")
    else:
        result = None
    return _chainsaw_predictions_from_tsv(structure, tsv_path), result


def _maybe_run_sword2_original(structure, args: argparse.Namespace) -> tuple[list[PartitionPrediction], ToolRunResult | None]:
    raw_dir = args.results_dir / "raw" / "sword2-original" / structure.entry.entry_id
    summary_path = raw_dir / "SWORD2_summary.json"
    if not args.skip_existing or not summary_path.exists():
        result = Sword2OriginalRunner(conda_env=args.sword2_original_conda_env).run(
            structure.entry.pdb_id,
            raw_dir,
            structure_file=structure.pdb_path,
            chain_id=structure.entry.chain_id,
        )
        summary_path = _check_run_or_summary(result, "sword2-original", raw_dir, "SWORD2_summary.json")
    else:
        result = None
    if not summary_path.exists():
        summary_path = _find_one(raw_dir, "SWORD2_summary.json")
    mapping_path = summary_path.parent / "residue_mapping.txt"
    return (
        load_original_summary_partitions(
            summary_path,
            numbering=structure.numbering,
            residue_mapping_path=mapping_path,
            chain_id=structure.entry.chain_id,
        ),
        result,
    )


def run_tool_for_structure(structure, tool: str, args: argparse.Namespace) -> tuple[list[ScoreRow], RunRow]:
    runners = {
        "sword2-rust": _maybe_run_sword2_rust,
        "merizo": _maybe_run_merizo,
        "chainsaw": _maybe_run_chainsaw,
        "sword2-original": _maybe_run_sword2_original,
    }
    predictions, result = runners[tool](structure, args)
    return _rows_from_predictions(
        structure,
        tool=tool,
        predictions=predictions,
        result=result,
        reused_output=result is None,
    )


def _sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _package_version(package: str) -> str | None:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return None


def _git_commit() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def write_manifest(
    args: argparse.Namespace,
    selected_tools: list[str],
    n_entries: int,
    n_score_rows: int,
    n_run_rows: int,
    n_failure_rows: int,
) -> None:
    dataset_file = dataset_path(args.dataset)
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": args.dataset,
        "dataset_path": str(dataset_file),
        "dataset_sha256": _sha256(dataset_file),
        "limit": args.limit,
        "n_entries": n_entries,
        "tools": selected_tools,
        "cache_dir": str(args.cache_dir),
        "results_dir": str(args.results_dir),
        "reuse_tool_results_dir": str(args.reuse_tool_results_dir) if args.reuse_tool_results_dir else None,
        "skip_existing": args.skip_existing,
        "strict": args.strict,
        "merizo_device": args.merizo_device,
        "allow_dl_cpu": args.allow_dl_cpu,
        "sword2_experiments": args.sword2_experiments,
        "sword2_extra_args": args.sword2_extra_args,
        "score_csv_input": str(args.score_csv) if args.score_csv else None,
        "n_score_rows": n_score_rows,
        "n_run_rows": n_run_rows,
        "n_failure_rows": n_failure_rows,
        "git_commit": _git_commit(),
        "packages": {
            package: _package_version(package)
            for package in ["numpy", "pandas", "scipy", "matplotlib", "seaborn"]
        },
    }
    (args.results_dir / "benchmark_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))


def main() -> int:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    scores_csv = args.results_dir / "scores.csv"
    previous_runtime_lookup = _score_runtime_lookup(scores_csv)

    if args.score_csv is not None:
        summarize_existing(args.score_csv, args.results_dir)
        return 0

    entries = load_dataset(args.dataset, limit=args.limit)
    structures = []
    for index, entry in enumerate(entries, start=1):
        _print_progress(index, len(entries), entry.entry_id, "prepare")
        structures.append(
            canonicalize_entry(entry, args.cache_dir, download=not args.no_download)
        )

    selected_tools = [tool.strip() for tool in args.tools.split(",") if tool.strip()]
    rows: list[ScoreRow] = []
    run_rows: list[RunRow] = []
    failure_rows: list[FailureRow] = []
    batched_tool_results: dict[
        str,
        dict[str, BatchedToolResult],
    ] = {}
    failed_batched_tools: set[str] = set()
    batch_runners = {
        "merizo": _maybe_run_merizo_batch,
        "chainsaw": _maybe_run_chainsaw_batch,
    }
    for tool, batch_runner in batch_runners.items():
        if tool not in selected_tools:
            continue
        try:
            _print_progress(0, len(structures), "batch", tool)
            batched_tool_results[tool] = batch_runner(structures, args)
        except Exception as exc:  # noqa: BLE001 - report batch failures per entry.
            message = f"[batch][{tool}] skipped: {exc}"
            if args.strict:
                raise RuntimeError(message) from exc
            failed_batched_tools.add(tool)
            for structure in structures:
                failure_rows.append(_failure_row_from_exception(structure, tool, exc, stage="batch"))
            print(message, file=sys.stderr)

    total_tool_runs = len(structures) * len(selected_tools)
    current_tool_run = 0
    for structure in structures:
        for tool in selected_tools:
            current_tool_run += 1
            _print_progress(current_tool_run, total_tool_runs, structure.entry.entry_id, tool)
            if tool in failed_batched_tools:
                continue
            try:
                if tool in batched_tool_results:
                    batched_entry = batched_tool_results[tool][structure.entry.entry_id]
                    if batched_entry.error is not None:
                        raise batched_entry.error
                    predictions, result, reused_output, batch_n_entries = batched_entry
                    score_rows, run_row = _rows_from_predictions(
                        structure,
                        tool=tool,
                        predictions=predictions,
                        result=result,
                        reused_output=reused_output,
                        batch_n_entries=batch_n_entries,
                    )
                else:
                    score_rows, run_row = run_tool_for_structure(structure, tool, args)
                rows.extend(score_rows)
                run_rows.append(run_row)
            except Exception as exc:  # noqa: BLE001 - CLI should report per-tool failures.
                message = f"[{structure.entry.entry_id}][{tool}] skipped: {exc}"
                failure_rows.append(_failure_row_from_exception(structure, tool, exc))
                if args.strict:
                    raise RuntimeError(message) from exc
                print(message, file=sys.stderr)

    rows = _fill_missing_score_runtimes(rows, previous_runtime_lookup)
    rows = add_sword2_oracle(rows)
    rows = add_sword2_oracle_s(rows)
    write_runs_csv(run_rows, args.results_dir / "runs.csv")
    write_failures_csv(failure_rows, args.results_dir / "failures.csv")
    write_manifest(args, selected_tools, len(structures), len(rows), len(run_rows), len(failure_rows))
    if not rows:
        print("No score rows were produced.", file=sys.stderr)
        return 1

    write_scores_csv(rows, scores_csv)
    summarize_existing(scores_csv, args.results_dir)
    print(f"Wrote {len(rows)} score rows to {scores_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
