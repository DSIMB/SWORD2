from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import unicodedata
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
from benchmark.runners.base import PartitionPrediction, ToolRunResult, run_timed
from benchmark.runners.chainsaw import ChainsawRunner, chainsaw_to_common_chopping, parse_chainsaw_tsv
from benchmark.runners.merizo import MerizoRunner, merizo_to_common_chopping, parse_merizo_tsv
from benchmark.runners.sword2_orig import Sword2OriginalRunner, load_original_summary_partitions
from benchmark.runners.sword2_rust import Sword2RustRunner, load_summary_partitions
from benchmark.score import (
    FailureRow,
    LockedRunRow,
    RunRow,
    ScoreRow,
    add_sword2_oracle,
    add_sword2_oracle_s,
    score_prediction,
    write_failures_csv,
    write_locked_runs_csv,
    write_runs_csv,
    write_scores_csv,
)
from benchmark.structures import CanonicalStructure, canonicalize_entry
from benchmark.numbering import numbering_from_pdb
from benchmark.factorized_ranker.runtime_freeze import (
    canonical_id_set_hash,
    canonical_json_bytes,
    hash_file_or_tree,
    sha256_file,
    verify_runtime_freeze,
)


LOCKED_ROLES = (
    "legacy-accuracy",
    "factorized-accuracy",
    "paired-sword-resources",
)
LOCKED_ENVIRONMENT = {
    "PYTHONHASHSEED": "0",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "RAYON_NUM_THREADS": "1",
}


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value}")


def _locked_entry_id(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value.strip() != value
        or unicodedata.normalize("NFC", value) != value
        or "\0" in value
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise ValueError("locked entry_id is not canonical")
    return value


def counterbalanced_pair_orders(entry_ids: list[str]) -> dict[str, str]:
    canonical = [_locked_entry_id(entry_id) for entry_id in entry_ids]
    if len(canonical) != len(set(canonical)):
        raise ValueError("locked entry IDs contain duplicates")
    ordered = sorted(
        canonical,
        key=lambda entry_id: (
            hashlib.sha256(b"37\0" + entry_id.encode("utf-8")).digest(),
            entry_id,
        ),
    )
    return {
        entry_id: "legacy" if index % 2 == 0 else "factorized"
        for index, entry_id in enumerate(ordered)
    }


def parse_selector_status(
    data: bytes,
    expected_selector: str,
    *,
    require_success: bool = True,
) -> dict[str, object]:
    if expected_selector not in {"legacy", "factorized"}:
        raise ValueError("unknown expected selector")
    try:
        payload = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=_unique_json_object,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("selector status is not valid JSON") from error
    expected_keys = {
        "error_code",
        "excluded_candidate_count",
        "fallback",
        "requested_selector",
        "schema_version",
        "selector_used",
    }
    if not isinstance(payload, dict) or set(payload) != expected_keys:
        raise ValueError("selector status schema mismatch")
    if canonical_json_bytes(payload) != data:
        raise ValueError("selector status is not canonically serialized")
    exclusions = payload["excluded_candidate_count"]
    if type(exclusions) is not int or exclusions < 0:
        raise ValueError("selector exclusion count is invalid")
    if type(payload["schema_version"]) is not int or payload["schema_version"] != 1:
        raise ValueError("selector status schema version mismatch")
    requested = payload["requested_selector"]
    used = payload["selector_used"]
    fallback = payload["fallback"]
    error_code = payload["error_code"]
    if requested not in {"legacy", "factorized"} or used not in {
        "legacy",
        "factorized",
    }:
        raise ValueError("selector status has an unknown selector")
    if type(fallback) is not bool:
        raise ValueError("selector fallback field is not a bool")
    if error_code is not None and (
        not isinstance(error_code, str)
        or not error_code
        or error_code.strip() != error_code
    ):
        raise ValueError("selector error code is invalid")
    if fallback:
        if requested != "factorized" or used != "legacy" or error_code is None:
            raise ValueError("selector fallback status is inconsistent")
    elif used != requested or error_code is not None:
        raise ValueError("selector success status is inconsistent")
    if requested == "legacy" and exclusions != 0:
        raise ValueError("legacy selector status records exclusions")
    if require_success and (
        requested != expected_selector or used != expected_selector or fallback
    ):
        raise ValueError("selector status does not prove the requested selector")
    return payload


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
    parser.add_argument("--locked-factorized-manifest", type=Path, default=None)
    parser.add_argument("--locked-role", choices=LOCKED_ROLES, default=None)
    parser.add_argument(
        "--locked-artifact",
        action="append",
        default=None,
        metavar="ROLE=PATH",
    )
    parser.add_argument("--chainsaw-expected-success-ids", type=Path, default=None)
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


class LockedBenchmarkError(ValueError):
    pass


@dataclass(frozen=True)
class LockedPreflight:
    repo_root: Path
    binary: Path
    runtime_manifest: dict[str, object]
    runtime_manifest_sha256: str
    model_manifest_sha256: str
    dataset_file: Path
    dataset_sha256: str
    dataset_ids: tuple[str, ...]
    dataset_id_set_sha256: str
    structures: tuple[CanonicalStructure, ...]
    structure_sha256s: dict[str, str]
    structure_tree_sha256: str
    external_artifacts: dict[str, dict[str, object]]
    external_artifact_paths: dict[str, Path]
    expected_chainsaw_ids: tuple[str, ...]
    expected_chainsaw_path: Path | None
    expected_chainsaw_file_sha256: str | None
    selected_tools: tuple[str, ...]
    gnu_time_version: str


def _length_framed_mapping_hash(values: dict[str, str]) -> str:
    digest = hashlib.sha256()
    for key, value in sorted(values.items()):
        key_bytes = key.encode("utf-8")
        value_bytes = value.encode("ascii")
        digest.update(len(key_bytes).to_bytes(8, "big"))
        digest.update(key_bytes)
        digest.update(len(value_bytes).to_bytes(8, "big"))
        digest.update(value_bytes)
    return digest.hexdigest()


def _path_contains(parent: Path, child: Path) -> bool:
    try:
        child.relative_to(parent)
    except ValueError:
        return False
    return True


def _nonexistent_resolved(path: Path) -> Path:
    path = Path(path).absolute()
    if path.exists() or path.is_symlink():
        raise LockedBenchmarkError(f"locked output already exists: {path}")
    parent = path.parent.resolve(strict=True)
    return parent / path.name


def _reject_output_input_alias(output: Path, inputs: list[tuple[str, Path]]) -> None:
    for role, raw_input in inputs:
        input_path = raw_input.resolve(strict=True)
        if output == input_path or _path_contains(output, input_path):
            raise LockedBenchmarkError(f"locked output aliases input role {role}")
        if input_path.is_dir() and _path_contains(input_path, output):
            raise LockedBenchmarkError(f"locked output is inside input role {role}")


def _parse_locked_artifacts(values: list[str] | None) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values or []:
        if "=" not in value:
            raise LockedBenchmarkError("locked artifacts must use ROLE=PATH")
        role, raw_path = value.split("=", 1)
        if not role or role.strip() != role or not raw_path:
            raise LockedBenchmarkError("locked artifact role/path is empty or noncanonical")
        if role in result:
            raise LockedBenchmarkError(f"duplicate locked artifact role {role}")
        result[role] = Path(raw_path)
    return result


def _artifact_paths_for_recheck(paths: dict[str, Path]) -> dict[str, Path]:
    return {
        role: Path(os.path.abspath(path))
        for role, path in sorted(paths.items())
    }


def _read_expected_success_ids(path: Path) -> tuple[tuple[str, ...], str]:
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise LockedBenchmarkError("Chainsaw expected-success input is not a regular file")
    before_hash = sha256_file(path)
    data = path.read_bytes()
    if sha256_file(path) != before_hash:
        raise LockedBenchmarkError("Chainsaw expected-success input changed while reading")
    if not data or not data.endswith(b"\n") or b"\r" in data:
        raise LockedBenchmarkError("Chainsaw expected-success IDs have noncanonical newlines")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise LockedBenchmarkError("Chainsaw expected-success IDs are not UTF-8") from error
    ids = text.splitlines()
    canonical = [_locked_entry_id(entry_id) for entry_id in ids]
    if canonical != sorted(canonical) or len(canonical) != len(set(canonical)):
        raise LockedBenchmarkError("Chainsaw expected-success IDs are not sorted and unique")
    if hashlib.sha256(data).hexdigest() != before_hash:
        raise LockedBenchmarkError("Chainsaw expected-success bytes changed while reading")
    return tuple(canonical), before_hash


def _read_locked_dataset_ids(path: Path) -> list[str]:
    ids: list[str] = []
    try:
        with Path(path).open(newline="", encoding="utf-8") as handle:
            for row in csv.reader(handle):
                if not row or row[0].startswith("#"):
                    continue
                if len(row) < 7:
                    raise LockedBenchmarkError("locked dataset metadata row is truncated")
                ids.append(_locked_entry_id(row[1]))
    except (OSError, UnicodeError, csv.Error) as error:
        raise LockedBenchmarkError("locked dataset metadata is malformed") from error
    if not ids or len(ids) != len(set(ids)):
        raise LockedBenchmarkError("locked dataset metadata IDs are empty or duplicated")
    return ids


def _gnu_time_version() -> str:
    time_path = Path("/usr/bin/time")
    if not time_path.is_file() or not os.access(time_path, os.X_OK):
        raise LockedBenchmarkError("/usr/bin/time is unavailable")
    completed = subprocess.run(
        [os.fspath(time_path), "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    version = "\n".join(part.strip() for part in (completed.stdout, completed.stderr) if part.strip())
    if completed.returncode != 0 or "GNU time" not in version:
        raise LockedBenchmarkError("/usr/bin/time is not GNU time")
    return version


def _locked_role_contract(args: argparse.Namespace) -> tuple[tuple[str, ...], dict[str, Path]]:
    selected = tuple(tool.strip() for tool in args.tools.split(",") if tool.strip())
    if len(selected) != len(set(selected)) or args.tools != ",".join(selected):
        raise LockedBenchmarkError("locked tool list is duplicated or noncanonical")
    artifacts = _parse_locked_artifacts(args.locked_artifact)
    if args.locked_role == "legacy-accuracy":
        expected_tools = ("sword2-rust", "merizo", "chainsaw")
        expected_artifacts = {
            "merizo_python",
            "merizo_source_tree",
            "merizo_model_tree",
            "chainsaw_python",
            "chainsaw_source_tree",
            "chainsaw_model_tree",
        }
        if args.sword2_extra_args is not None:
            raise LockedBenchmarkError("legacy accuracy forbids SWORD extra arguments")
        if args.chainsaw_expected_success_ids is None:
            raise LockedBenchmarkError("legacy accuracy requires frozen Chainsaw success IDs")
        if set(artifacts) != expected_artifacts:
            raise LockedBenchmarkError("legacy accuracy external artifact roles mismatch")
        if args.merizo_device != "cuda" or args.allow_dl_cpu:
            raise LockedBenchmarkError("legacy accuracy requires frozen CUDA competitor settings")
    elif args.locked_role == "factorized-accuracy":
        expected_tools = ("sword2-rust",)
        if args.sword2_extra_args != "--use-factorized-ranker":
            raise LockedBenchmarkError("factorized accuracy requires the exact factorized flag")
        if artifacts or args.chainsaw_expected_success_ids is not None:
            raise LockedBenchmarkError("factorized accuracy has unexpected competitor inputs")
    elif args.locked_role == "paired-sword-resources":
        expected_tools = ("sword2-rust",)
        if args.sword2_extra_args is not None:
            raise LockedBenchmarkError("paired resources controls both SWORD commands internally")
        if artifacts or args.chainsaw_expected_success_ids is not None:
            raise LockedBenchmarkError("paired resources has unexpected competitor inputs")
    else:
        raise LockedBenchmarkError("unknown locked role")
    if selected != expected_tools:
        raise LockedBenchmarkError("locked role tool set/order mismatch")
    return selected, artifacts


def _locked_preflight(args: argparse.Namespace) -> LockedPreflight:
    if args.results_dir.exists() or args.results_dir.is_symlink():
        raise LockedBenchmarkError("locked results directory must be absent")
    if any(
        [
            args.score_csv is not None,
            args.skip_existing,
            args.reuse_tool_results_dir is not None,
            args.limit is not None,
        ]
    ):
        raise LockedBenchmarkError("locked collection forbids score/reuse/skip/limit")
    if not args.strict or not args.no_download or args.sword2_threads != 1:
        raise LockedBenchmarkError("locked collection requires strict, no-download, and one SWORD thread")
    if args.sword2_experiments is not None:
        raise LockedBenchmarkError("locked collection forbids SWORD2_EXPERIMENTS")
    if os.environ.get("SWORD2_EXPERIMENTS") is not None:
        raise LockedBenchmarkError(
            "locked collection requires ambient SWORD2_EXPERIMENTS to be unset"
        )
    selected_tools, artifact_paths = _locked_role_contract(args)
    for name, expected in LOCKED_ENVIRONMENT.items():
        if os.environ.get(name) != expected:
            raise LockedBenchmarkError(f"locked environment {name} must equal {expected}")

    repo_root = Path(__file__).resolve().parents[1]
    binary = repo_root / "target/release/sword2"
    runtime_path = Path(args.locked_factorized_manifest)
    runtime_manifest = verify_runtime_freeze(runtime_path, binary=binary, repo_root=repo_root)
    for name, expected in runtime_manifest["build_environment"].items():
        if os.environ.get(name) != expected:
            raise LockedBenchmarkError(f"build environment {name} differs from runtime freeze")
    runtime_manifest_sha256 = sha256_file(runtime_path)
    model_manifest_sha256 = str(runtime_manifest["model_manifest_sha256"])

    dataset_file = Path(dataset_path(args.dataset))
    if dataset_file.is_symlink() or not dataset_file.is_file():
        raise LockedBenchmarkError("locked dataset metadata is missing or symlinked")
    dataset_sha256 = sha256_file(dataset_file)
    raw_ids = _read_locked_dataset_ids(dataset_file)
    entries = load_dataset(args.dataset)
    if sha256_file(dataset_file) != dataset_sha256:
        raise LockedBenchmarkError("locked dataset metadata changed while parsing")
    ids = [_locked_entry_id(entry.entry_id) for entry in entries]
    if not ids or len(ids) != len(set(ids)) or ids != raw_ids:
        raise LockedBenchmarkError("locked dataset IDs are empty or duplicated")
    if args.dataset == "cath663" and len(ids) != 663:
        raise LockedBenchmarkError("CATH-663 locked population is not exactly 663")

    cache_root = Path(args.cache_dir).resolve(strict=True)
    chains_root = cache_root / "chains"
    if chains_root.is_symlink() or not chains_root.is_dir():
        raise LockedBenchmarkError("locked canonical chain cache is unavailable")
    structures: list[CanonicalStructure] = []
    structure_hashes: dict[str, str] = {}
    for entry in entries:
        structure_path = chains_root / f"{entry.entry_id}.pdb"
        if structure_path.is_symlink() or not structure_path.is_file():
            raise LockedBenchmarkError(f"locked canonical structure is missing: {entry.entry_id}")
        if structure_path.resolve(strict=True).parent != chains_root.resolve(strict=True):
            raise LockedBenchmarkError("locked canonical structure escapes the cache")
        digest = sha256_file(structure_path)
        structure_hashes[entry.entry_id] = digest
        structures.append(
            CanonicalStructure(
                entry=entry,
                pdb_path=structure_path.resolve(strict=True),
                numbering=numbering_from_pdb(structure_path, chain_id=entry.chain_id),
            )
        )

    external_artifacts: dict[str, dict[str, object]] = {}
    for role, path in sorted(artifact_paths.items()):
        descriptor = hash_file_or_tree(path)
        descriptor["role"] = role
        external_artifacts[role] = descriptor
    if args.locked_role == "legacy-accuracy":
        merizo = MerizoRunner()
        chainsaw = ChainsawRunner()
        if Path(artifact_paths["merizo_python"]).resolve(
            strict=True
        ) != merizo.python.resolve(strict=True):
            raise LockedBenchmarkError("Merizo Python artifact does not match the configured runner")
        if Path(artifact_paths["merizo_source_tree"]).resolve(
            strict=True
        ) != merizo.merizo_dir.resolve(strict=True):
            raise LockedBenchmarkError("Merizo source artifact does not match the configured runner")
        if Path(artifact_paths["chainsaw_python"]).resolve(
            strict=True
        ) != chainsaw.python.resolve(strict=True):
            raise LockedBenchmarkError("Chainsaw Python artifact does not match the configured runner")
        if Path(artifact_paths["chainsaw_source_tree"]).resolve(
            strict=True
        ) != chainsaw.chainsaw_dir.resolve(strict=True):
            raise LockedBenchmarkError("Chainsaw source artifact does not match the configured runner")

    expected_ids: tuple[str, ...] = ()
    expected_file_hash: str | None = None
    if args.chainsaw_expected_success_ids is not None:
        expected_ids, expected_file_hash = _read_expected_success_ids(args.chainsaw_expected_success_ids)
        if not set(expected_ids).issubset(ids):
            raise LockedBenchmarkError("Chainsaw expected-success IDs are not a dataset subset")

    output = _nonexistent_resolved(args.results_dir)
    inputs: list[tuple[str, Path]] = [
        ("cache", cache_root),
        ("dataset", dataset_file),
        ("runtime_manifest", runtime_path),
        ("binary", binary),
    ]
    inputs.extend((role, path) for role, path in artifact_paths.items())
    if args.chainsaw_expected_success_ids is not None:
        inputs.append(("chainsaw_expected_success_ids", args.chainsaw_expected_success_ids))
    _reject_output_input_alias(output, inputs)

    return LockedPreflight(
        repo_root=repo_root,
        binary=binary.resolve(strict=True),
        runtime_manifest=runtime_manifest,
        runtime_manifest_sha256=runtime_manifest_sha256,
        model_manifest_sha256=model_manifest_sha256,
        dataset_file=dataset_file.resolve(strict=True),
        dataset_sha256=dataset_sha256,
        dataset_ids=tuple(sorted(ids)),
        dataset_id_set_sha256=canonical_id_set_hash(ids),
        structures=tuple(structures),
        structure_sha256s=dict(sorted(structure_hashes.items())),
        structure_tree_sha256=_length_framed_mapping_hash(structure_hashes),
        external_artifacts=external_artifacts,
        external_artifact_paths=_artifact_paths_for_recheck(artifact_paths),
        expected_chainsaw_ids=expected_ids,
        expected_chainsaw_path=(
            Path(args.chainsaw_expected_success_ids).resolve(strict=True)
            if args.chainsaw_expected_success_ids is not None
            else None
        ),
        expected_chainsaw_file_sha256=expected_file_hash,
        selected_tools=selected_tools,
        gnu_time_version=_gnu_time_version(),
    )


def _canonical_json_text(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def _atomic_write_bytes(path: Path, data: bytes, *, absent_only: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if absent_only and (path.exists() or path.is_symlink()):
        raise FileExistsError(path)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(temporary)
    try:
        with temporary.open("xb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        if absent_only:
            os.link(temporary, path)
        else:
            os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_write_rows(writer, rows: list, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.rows-{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(temporary)
    writer(rows, temporary)
    os.replace(temporary, path)


def _normalized_command(
    command: list[str],
    *,
    preflight: LockedPreflight,
    args: argparse.Namespace,
    structure: Path | None = None,
    output: Path | None = None,
    status: Path | None = None,
) -> list[str]:
    replacements: dict[str, str] = {}

    def add_path(raw_path: Path, semantic: str, *, must_exist: bool) -> None:
        path = Path(raw_path)
        variants = {os.fspath(path), os.fspath(path.absolute())}
        try:
            variants.add(os.fspath(path.resolve(strict=must_exist)))
        except OSError:
            if must_exist:
                raise
        for variant in variants:
            if variant:
                replacements[variant] = semantic

    add_path(preflight.repo_root, "<repo>", must_exist=True)
    add_path(preflight.binary, "<binary>", must_exist=True)
    add_path(Path(args.cache_dir), "<cache>", must_exist=True)
    add_path(Path(args.results_dir), "<results>", must_exist=False)
    add_path(
        Path(args.locked_factorized_manifest),
        "<runtime_manifest>",
        must_exist=True,
    )
    if args.chainsaw_expected_success_ids is not None:
        add_path(
            Path(args.chainsaw_expected_success_ids),
            "<chainsaw_expected_success_ids>",
            must_exist=True,
        )
    replacements.update(
        {
            os.fspath(structure_input.pdb_path): f"<input_structure:{structure_input.entry.entry_id}>"
            for structure_input in preflight.structures
        }
    )
    raw_artifacts = _parse_locked_artifacts(args.locked_artifact)
    for role, path in preflight.external_artifact_paths.items():
        semantic = f"<locked_artifact:{role}>"
        add_path(path, semantic, must_exist=True)
        add_path(raw_artifacts[role], semantic, must_exist=True)
    if structure is not None:
        replacements[os.fspath(structure.resolve(strict=True))] = "<input_structure>"
    if output is not None:
        replacements[os.fspath(output.absolute())] = "<process_output>"
    if status is not None:
        replacements[os.fspath(status.absolute())] = "<selector_status>"
    normalized: list[str] = []
    for token in command:
        replacement = replacements.get(token)
        if replacement is None:
            replacement = token
            for physical, semantic in sorted(replacements.items(), key=lambda item: -len(item[0])):
                if token.startswith(physical + os.sep):
                    replacement = semantic + token[len(physical) :]
                    break
                marker = "=" + physical
                if marker in token:
                    replacement = token.replace(marker, "=" + semantic)
                    break
        normalized.append(replacement)
    path_flags = {
        "--cache-dir": "<cache>",
        "--results-dir": "<results>",
        "--locked-factorized-manifest": "<runtime_manifest>",
        "--chainsaw-expected-success-ids": "<chainsaw_expected_success_ids>",
    }
    for index, token in enumerate(normalized):
        semantic = path_flags.get(token)
        if semantic is not None:
            if index + 1 >= len(normalized):
                raise LockedBenchmarkError(f"normalized command lacks value after {token}")
            normalized[index + 1] = semantic
            continue
        for flag, flag_semantic in path_flags.items():
            if token.startswith(flag + "="):
                normalized[index] = f"{flag}={flag_semantic}"
                break
        if token == "--locked-artifact":
            if index + 1 >= len(normalized) or "=" not in normalized[index + 1]:
                raise LockedBenchmarkError("normalized locked artifact lacks ROLE=PATH")
            role = normalized[index + 1].split("=", 1)[0]
            if role not in preflight.external_artifact_paths:
                raise LockedBenchmarkError("normalized locked artifact has an unknown role")
            normalized[index + 1] = f"{role}=<locked_artifact:{role}>"
    return normalized


def _write_process_log(path: Path, text: str) -> str:
    data = text.encode("utf-8")
    _atomic_write_bytes(path, data, absent_only=True)
    return hashlib.sha256(data).hexdigest()


def _locked_sword_run(
    structure: CanonicalStructure,
    selector: str,
    process_dir: Path,
    *,
    pair_order: str,
    pair_position: int | str,
    locked_role: str,
    preflight: LockedPreflight,
    args: argparse.Namespace,
) -> tuple[list[PartitionPrediction], LockedRunRow]:
    if selector not in {"legacy", "factorized"}:
        raise LockedBenchmarkError("unknown locked SWORD selector")
    if sha256_file(structure.pdb_path) != preflight.structure_sha256s[structure.entry.entry_id]:
        raise LockedBenchmarkError("locked structure changed before SWORD launch")
    status_path = process_dir / "selector_status.json"
    extra_args = "--use-factorized-ranker" if selector == "factorized" else None
    result = Sword2RustRunner(
        repo_dir=preflight.repo_root,
        binary=preflight.binary,
        threads=1,
        extra_args=extra_args,
        fresh_only=True,
        locked_env=LOCKED_ENVIRONMENT,
        selector_status_path=status_path,
    ).run(structure.pdb_path, process_dir)
    stdout_path = process_dir / "stdout.log"
    stderr_path = process_dir / "stderr.log"
    stdout_sha256 = _write_process_log(stdout_path, result.stdout)
    stderr_sha256 = _write_process_log(stderr_path, result.stderr)
    if result.returncode != 0:
        raise LockedBenchmarkError(f"locked SWORD process exited {result.returncode}")
    status_data = status_path.read_bytes()
    status = parse_selector_status(
        status_data,
        selector,
        require_success=False,
    )
    status_sha256 = hashlib.sha256(status_data).hexdigest()
    summaries = sorted(process_dir.rglob("summary.json"))
    if len(summaries) != 1:
        raise LockedBenchmarkError("locked SWORD output does not contain exactly one summary")
    summary = summaries[0]
    mapping_path = summary.parent / "residue_mapping.txt"
    predictions = load_summary_partitions(
        summary,
        numbering=structure.numbering,
        residue_mapping_path=mapping_path,
        chain_id=structure.entry.chain_id,
    )
    if sum(pred.variant == "optimal" and pred.name == "Optimal partition" for pred in predictions) != 1:
        raise LockedBenchmarkError("locked SWORD summary lacks one rank-1 partition")
    relative_summary = summary.relative_to(Path(args.results_dir)).as_posix()
    row = LockedRunRow(
        dataset=structure.entry.dataset,
        entry_id=structure.entry.entry_id,
        pdb_id=structure.entry.pdb_id,
        chain_id=structure.entry.chain_id,
        tool="sword2-rust",
        selector_variant=selector,
        locked_role=locked_role,
        pair_order=pair_order,
        pair_position=pair_position,
        reused_output=False,
        returncode=result.returncode,
        runtime_s=result.runtime_s,
        peak_rss_kb=result.peak_rss_kb,
        command_json=_canonical_json_text(result.command),
        normalized_command_json=_canonical_json_text(
            _normalized_command(
                result.command,
                preflight=preflight,
                args=args,
                structure=structure.pdb_path,
                output=process_dir,
                status=status_path,
            )
        ),
        cwd_role="repository_root",
        input_structure_sha256=preflight.structure_sha256s[structure.entry.entry_id],
        raw_summary_role=relative_summary,
        raw_summary_sha256=sha256_file(summary),
        stdout_sha256=stdout_sha256,
        stderr_sha256=stderr_sha256,
        selector_status_sha256=status_sha256,
        requested_selector=str(status["requested_selector"]),
        selector_used=str(status["selector_used"]),
        fallback=bool(status["fallback"]),
        error_code=("" if status["error_code"] is None else str(status["error_code"])),
        selector_warning_code=("factorized_fallback" if status["fallback"] is True else ""),
        excluded_candidate_count=int(status["excluded_candidate_count"]),
        binary_sha256=str(preflight.runtime_manifest["binary_sha256"]),
        model_manifest_sha256=preflight.model_manifest_sha256,
        runtime_manifest_sha256=preflight.runtime_manifest_sha256,
    )
    return predictions, row


def _require_locked_sword_success(row: LockedRunRow, expected_selector: str) -> None:
    if (
        row.requested_selector != expected_selector
        or row.selector_used != expected_selector
        or row.fallback is not False
        or row.error_code
        or row.selector_warning_code
        or (expected_selector == "legacy" and row.excluded_candidate_count != 0)
    ):
        raise LockedBenchmarkError(
            "locked SWORD status does not prove the requested selector"
        )


def _locked_competitor_row(
    structure: CanonicalStructure,
    tool: str,
    raw_path: Path,
    result: ToolRunResult,
    *,
    stdout_sha256: str,
    stderr_sha256: str,
    cwd_role: str,
    preflight: LockedPreflight,
    args: argparse.Namespace,
) -> LockedRunRow:
    return LockedRunRow(
        dataset=structure.entry.dataset,
        entry_id=structure.entry.entry_id,
        pdb_id=structure.entry.pdb_id,
        chain_id=structure.entry.chain_id,
        tool=tool,
        selector_variant="",
        locked_role="legacy-accuracy",
        pair_order="",
        pair_position="",
        reused_output=False,
        returncode=result.returncode,
        runtime_s=result.runtime_s,
        peak_rss_kb=result.peak_rss_kb,
        command_json=_canonical_json_text(result.command),
        normalized_command_json=_canonical_json_text(
            _normalized_command(
                result.command,
                preflight=preflight,
                args=args,
            )
        ),
        cwd_role=cwd_role,
        input_structure_sha256=preflight.structure_sha256s[structure.entry.entry_id],
        raw_summary_role=raw_path.relative_to(Path(args.results_dir)).as_posix(),
        raw_summary_sha256=sha256_file(raw_path),
        stdout_sha256=stdout_sha256,
        stderr_sha256=stderr_sha256,
        selector_status_sha256="",
        requested_selector="",
        selector_used="",
        fallback="",
        error_code="",
        selector_warning_code="",
        excluded_candidate_count="",
        binary_sha256="",
        model_manifest_sha256="",
        runtime_manifest_sha256="",
    )


def _locked_merizo_accuracy(
    preflight: LockedPreflight,
    args: argparse.Namespace,
) -> tuple[list[ScoreRow], list[LockedRunRow]]:
    grouped: dict[str, list[CanonicalStructure]] = {}
    for structure in preflight.structures:
        grouped.setdefault(structure.entry.chain_id, []).append(structure)
    score_rows: list[ScoreRow] = []
    run_rows: list[LockedRunRow] = []
    raw_root = Path(args.results_dir) / "raw" / "merizo-cuda"
    log_root = Path(args.results_dir) / "raw" / "_logs"
    runner = MerizoRunner(device="cuda", require_cuda=False)
    for chain_id, structures in sorted(grouped.items()):
        inputs: list[tuple[Path, Path, str]] = []
        for structure in structures:
            output = raw_root / f"{structure.entry.entry_id}.tsv"
            if output.exists() or output.is_symlink():
                raise LockedBenchmarkError("locked Merizo output already exists")
            inputs.append((structure.pdb_path, output, chain_id))
        command = [
            os.fspath(runner.python.absolute()),
            "predict.py",
            "-i",
            *[os.fspath(structure_path.resolve(strict=True)) for structure_path, _, _ in inputs],
            "-d",
            "cuda",
            "--return_indices",
            "--output_headers",
            "--pdb_chain",
            chain_id,
        ]
        result = run_timed(
            command,
            cwd=runner.merizo_dir,
            env=LOCKED_ENVIRONMENT,
            timeout_s=None,
        )
        stdout_path = log_root / f"merizo-{chain_id}.stdout"
        stderr_path = log_root / f"merizo-{chain_id}.stderr"
        stdout_hash = _write_process_log(stdout_path, result.stdout)
        stderr_hash = _write_process_log(stderr_path, result.stderr)
        if result.returncode != 0:
            raise LockedBenchmarkError(f"locked Merizo batch exited {result.returncode}")
        _split_locked_merizo_stdout(result.stdout, inputs)
        for structure, (_, output, _) in zip(structures, inputs, strict=True):
            predictions = _merizo_predictions_from_tsv(structure, args, output)
            if len(predictions) != 1:
                raise LockedBenchmarkError("locked Merizo did not produce one prediction")
            score_rows.append(
                score_prediction(
                    structure.entry,
                    structure.numbering,
                    predictions[0],
                    runtime_s=predictions[0].runtime_s,
                    peak_rss_mb=None,
                )
            )
            run_rows.append(
                _locked_competitor_row(
                    structure,
                    "merizo",
                    output,
                    result,
                    stdout_sha256=stdout_hash,
                    stderr_sha256=stderr_hash,
                    cwd_role="locked_artifact:merizo_source_tree",
                    preflight=preflight,
                    args=args,
                )
            )
    return score_rows, run_rows


def _split_locked_merizo_stdout(
    stdout: str,
    inputs: list[tuple[Path, Path, str]],
) -> None:
    lines = stdout.splitlines()
    if not lines or not lines[0] or any(not line for line in lines):
        raise LockedBenchmarkError("locked Merizo stdout is empty or noncanonical")
    try:
        reader = csv.DictReader(lines, delimiter="\t")
        fieldnames = reader.fieldnames
        if (
            not fieldnames
            or "input" not in fieldnames
            or len(fieldnames) != len(set(fieldnames))
        ):
            raise LockedBenchmarkError("locked Merizo output header is invalid")
        rows_by_input: dict[str, list[dict[str, str]]] = {
            structure_path.name: [] for structure_path, _, _ in inputs
        }
        for row in reader:
            if set(row) != set(fieldnames) or any(value is None for value in row.values()):
                raise LockedBenchmarkError("locked Merizo output row schema is invalid")
            key = Path(row["input"]).name
            if key not in rows_by_input:
                raise LockedBenchmarkError("locked Merizo output contains an unknown input")
            rows_by_input[key].append(row)
    except csv.Error as error:
        raise LockedBenchmarkError("locked Merizo output is malformed") from error

    for structure_path, output, _ in inputs:
        rows = rows_by_input[structure_path.name]
        if len(rows) != 1:
            raise LockedBenchmarkError("locked Merizo output coverage is incomplete or duplicated")
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists() or output.is_symlink():
            raise LockedBenchmarkError("locked Merizo per-entry output already exists")
        temporary = output.with_name(f".{output.name}.tmp-{os.getpid()}")
        if temporary.exists() or temporary.is_symlink():
            raise LockedBenchmarkError("locked Merizo temporary output already exists")
        try:
            with temporary.open("x", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fieldnames,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)
                handle.flush()
                os.fsync(handle.fileno())
            os.link(temporary, output)
        finally:
            temporary.unlink(missing_ok=True)


def _split_locked_chainsaw_output(
    combined_output: Path,
    outputs: dict[str, Path],
) -> set[str]:
    if not combined_output.exists():
        return set()
    with combined_output.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames or "chain_id" not in fieldnames or len(fieldnames) != len(set(fieldnames)):
            raise LockedBenchmarkError("locked Chainsaw output header is invalid")
        rows_by_id: dict[str, list[dict[str, str]]] = {entry_id: [] for entry_id in outputs}
        for row in reader:
            entry_id = row.get("chain_id")
            if entry_id not in rows_by_id:
                raise LockedBenchmarkError("locked Chainsaw output contains an unknown ID")
            rows_by_id[entry_id].append(row)
    successes: set[str] = set()
    for entry_id, rows in rows_by_id.items():
        if len(rows) > 1:
            raise LockedBenchmarkError("locked Chainsaw output contains duplicate predictions")
        if not rows:
            continue
        output = outputs[entry_id]
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists() or output.is_symlink():
            raise LockedBenchmarkError("locked Chainsaw per-entry output already exists")
        temporary = output.with_name(f".{output.name}.tmp-{os.getpid()}")
        with temporary.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, output)
        temporary.unlink()
        successes.add(entry_id)
    return successes


def _locked_chainsaw_accuracy(
    preflight: LockedPreflight,
    args: argparse.Namespace,
) -> tuple[list[ScoreRow], list[LockedRunRow], list[FailureRow]]:
    runner = ChainsawRunner(require_cuda=False)
    raw_root = Path(args.results_dir) / "raw" / "chainsaw"
    stage_dir = raw_root / "_batch_stage"
    combined_output = raw_root / "_batch.tsv"
    if stage_dir.exists() or combined_output.exists():
        raise LockedBenchmarkError("locked Chainsaw batch output already exists")
    stage_dir.mkdir(parents=True)
    outputs: dict[str, Path] = {}
    for structure in preflight.structures:
        staged = stage_dir / f"{structure.entry.entry_id}{structure.pdb_path.suffix}"
        if staged.exists() or staged.is_symlink():
            raise LockedBenchmarkError("locked Chainsaw staged input already exists")
        shutil.copy2(structure.pdb_path, staged)
        outputs[structure.entry.entry_id] = raw_root / f"{structure.entry.entry_id}.tsv"
    command = [
        os.fspath(runner.python.absolute()),
        "get_predictions.py",
        "--structure_directory",
        os.fspath(stage_dir),
        "--output",
        os.fspath(combined_output),
    ]
    result = run_timed(command, cwd=runner.chainsaw_dir, env=LOCKED_ENVIRONMENT, timeout_s=None)
    stdout_hash = _write_process_log(raw_root / "batch.stdout", result.stdout)
    stderr_hash = _write_process_log(raw_root / "batch.stderr", result.stderr)
    if result.returncode != 0:
        raise LockedBenchmarkError(f"locked Chainsaw batch exited {result.returncode}")
    successes = _split_locked_chainsaw_output(combined_output, outputs)
    expected = set(preflight.expected_chainsaw_ids)
    if successes != expected:
        missing = sorted(expected - successes)
        unexpected = sorted(successes - expected)
        raise LockedBenchmarkError(
            f"locked Chainsaw success set mismatch: missing={missing[:3]} unexpected={unexpected[:3]}"
        )
    by_id = {structure.entry.entry_id: structure for structure in preflight.structures}
    score_rows: list[ScoreRow] = []
    run_rows: list[LockedRunRow] = []
    failure_rows: list[FailureRow] = []
    for entry_id in sorted(successes):
        structure = by_id[entry_id]
        output = outputs[entry_id]
        predictions = _chainsaw_predictions_from_tsv(structure, output)
        if len(predictions) != 1:
            raise LockedBenchmarkError("locked Chainsaw did not produce one prediction")
        score_rows.append(
            score_prediction(
                structure.entry,
                structure.numbering,
                predictions[0],
                runtime_s=predictions[0].runtime_s,
                peak_rss_mb=None,
            )
        )
        run_rows.append(
            _locked_competitor_row(
                structure,
                "chainsaw",
                output,
                result,
                stdout_sha256=stdout_hash,
                stderr_sha256=stderr_hash,
                cwd_role="locked_artifact:chainsaw_source_tree",
                preflight=preflight,
                args=args,
            )
        )
    for entry_id in sorted(set(preflight.dataset_ids) - successes):
        structure = by_id[entry_id]
        failure_rows.append(
            FailureRow(
                dataset=structure.entry.dataset,
                entry_id=entry_id,
                pdb_id=structure.entry.pdb_id,
                chain_id=structure.entry.chain_id,
                tool="chainsaw",
                stage="expected_failure",
                message="chainsaw_expected_failure",
                returncode=result.returncode,
                command=_canonical_json_text(result.command),
                stderr=None,
            )
        )
    return score_rows, run_rows, failure_rows


def _run_locked_pairs(
    preflight: LockedPreflight,
    args: argparse.Namespace,
) -> tuple[list[ScoreRow], list[LockedRunRow], list[FailureRow], dict[str, object]]:
    assignments = counterbalanced_pair_orders(list(preflight.dataset_ids))
    structures = {structure.entry.entry_id: structure for structure in preflight.structures}
    rows: list[LockedRunRow] = []
    runs_path = Path(args.results_dir) / "runs.csv"
    failures_path = Path(args.results_dir) / "failures.csv"
    _atomic_write_rows(write_locked_runs_csv, rows, runs_path)
    _atomic_write_rows(write_failures_csv, [], failures_path)
    ordered_assignments: list[list[str]] = []
    for entry_id, first in assignments.items():
        ordered_assignments.append([entry_id, first])
        selectors = (first, "factorized" if first == "legacy" else "legacy")
        pair_order = f"{first}_first"
        for position, selector in enumerate(selectors, start=1):
            process_dir = Path(args.results_dir) / "raw" / "paired" / entry_id / selector
            _, row = _locked_sword_run(
                structures[entry_id],
                selector,
                process_dir,
                pair_order=pair_order,
                pair_position=position,
                locked_role="paired-sword-resources",
                preflight=preflight,
                args=args,
            )
            rows.append(row)
            _atomic_write_rows(write_locked_runs_csv, rows, runs_path)
            _require_locked_sword_success(row, selector)
    order_bytes = canonical_json_bytes(ordered_assignments)
    order = {
        "seed": 37,
        "assignments": ordered_assignments,
        "assignments_sha256": hashlib.sha256(order_bytes).hexdigest(),
        "legacy_first_count": sum(value == "legacy" for value in assignments.values()),
        "factorized_first_count": sum(value == "factorized" for value in assignments.values()),
    }
    rows.sort(
        key=lambda row: (
            row.entry_id,
            row.tool,
            row.selector_variant,
            str(row.pair_position),
        )
    )
    _atomic_write_rows(write_locked_runs_csv, rows, runs_path)
    return [], rows, [], order


def _run_locked_accuracy(
    preflight: LockedPreflight,
    args: argparse.Namespace,
) -> tuple[list[ScoreRow], list[LockedRunRow], list[FailureRow], None]:
    selector = "factorized" if args.locked_role == "factorized-accuracy" else "legacy"
    scores: list[ScoreRow] = []
    runs: list[LockedRunRow] = []
    failures: list[FailureRow] = []
    runs_path = Path(args.results_dir) / "runs.csv"
    _atomic_write_rows(write_locked_runs_csv, runs, runs_path)
    for structure in sorted(preflight.structures, key=lambda item: item.entry.entry_id):
        process_dir = Path(args.results_dir) / "raw" / "sword2-rust" / structure.entry.entry_id
        predictions, row = _locked_sword_run(
            structure,
            selector,
            process_dir,
            pair_order="",
            pair_position="",
            locked_role=args.locked_role,
            preflight=preflight,
            args=args,
        )
        scores.extend(
            score_prediction(
                structure.entry,
                structure.numbering,
                prediction,
                runtime_s=row.runtime_s,
                peak_rss_mb=(row.peak_rss_kb / 1024.0 if row.peak_rss_kb else None),
            )
            for prediction in predictions
        )
        runs.append(row)
        _atomic_write_rows(write_locked_runs_csv, runs, runs_path)
        _require_locked_sword_success(row, selector)
    if args.locked_role == "legacy-accuracy":
        merizo_scores, merizo_runs = _locked_merizo_accuracy(preflight, args)
        scores.extend(merizo_scores)
        runs.extend(merizo_runs)
        _atomic_write_rows(write_locked_runs_csv, runs, runs_path)
        chainsaw_scores, chainsaw_runs, chainsaw_failures = _locked_chainsaw_accuracy(
            preflight,
            args,
        )
        scores.extend(chainsaw_scores)
        runs.extend(chainsaw_runs)
        failures.extend(chainsaw_failures)
        _atomic_write_rows(write_locked_runs_csv, runs, runs_path)
    _atomic_write_rows(
        write_failures_csv,
        sorted(failures, key=lambda row: (row.entry_id, row.tool, row.stage)),
        Path(args.results_dir) / "failures.csv",
    )
    runs.sort(
        key=lambda row: (
            row.entry_id,
            row.tool,
            row.selector_variant,
            str(row.pair_position),
        )
    )
    scores.sort(
        key=lambda row: (row.entry_id, row.tool, row.variant, row.partition)
    )
    _atomic_write_rows(write_locked_runs_csv, runs, runs_path)
    _atomic_write_rows(write_scores_csv, scores, Path(args.results_dir) / "scores.csv")
    return scores, runs, failures, None


def _locked_intent_payload(
    preflight: LockedPreflight,
    args: argparse.Namespace,
    created_at: str,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "status": "started",
        "created_at": created_at,
        "locked_role": args.locked_role,
        "dataset": args.dataset,
        "dataset_sha256": preflight.dataset_sha256,
        "dataset_id_count": len(preflight.dataset_ids),
        "dataset_id_set_sha256": preflight.dataset_id_set_sha256,
        "runtime_manifest_sha256": preflight.runtime_manifest_sha256,
        "model_manifest_sha256": preflight.model_manifest_sha256,
        "binary_sha256": preflight.runtime_manifest["binary_sha256"],
        "structure_tree_sha256": preflight.structure_tree_sha256,
        "tools": list(preflight.selected_tools),
        "environment": LOCKED_ENVIRONMENT,
        "external_artifacts": preflight.external_artifacts,
    }


def _locked_manifest_payload(
    preflight: LockedPreflight,
    args: argparse.Namespace,
    *,
    created_at: str,
    completed_at: str,
    scores: list[ScoreRow],
    runs: list[LockedRunRow],
    failures: list[FailureRow],
    order: dict[str, object] | None,
) -> dict[str, object]:
    runtime_path = Path(args.locked_factorized_manifest)
    current_runtime = verify_runtime_freeze(
        runtime_path,
        binary=preflight.binary,
        repo_root=preflight.repo_root,
    )
    if current_runtime != preflight.runtime_manifest or sha256_file(runtime_path) != preflight.runtime_manifest_sha256:
        raise LockedBenchmarkError("runtime freeze changed during locked collection")
    if sha256_file(preflight.dataset_file) != preflight.dataset_sha256:
        raise LockedBenchmarkError("dataset metadata changed during locked collection")
    if preflight.expected_chainsaw_path is not None:
        if (
            preflight.expected_chainsaw_file_sha256 is None
            or sha256_file(preflight.expected_chainsaw_path)
            != preflight.expected_chainsaw_file_sha256
        ):
            raise LockedBenchmarkError(
                "Chainsaw expected-success input changed during locked collection"
            )
    for entry_id, expected in preflight.structure_sha256s.items():
        path = Path(args.cache_dir).resolve(strict=True) / "chains" / f"{entry_id}.pdb"
        if sha256_file(path) != expected:
            raise LockedBenchmarkError("canonical structure changed during locked collection")
    for role, path in preflight.external_artifact_paths.items():
        current = hash_file_or_tree(path)
        current["role"] = role
        if current != preflight.external_artifacts[role]:
            raise LockedBenchmarkError(f"locked artifact changed during collection: {role}")

    scores_path = Path(args.results_dir) / "scores.csv"
    runs_path = Path(args.results_dir) / "runs.csv"
    failures_path = Path(args.results_dir) / "failures.csv"
    if not runs_path.is_file() or not failures_path.is_file():
        raise LockedBenchmarkError("locked evidence CSV is missing")
    if args.locked_role == "paired-sword-resources":
        if scores_path.exists():
            raise LockedBenchmarkError("resource-only collection unexpectedly produced scores")
        scores_hash: str | None = None
    else:
        if not scores_path.is_file():
            raise LockedBenchmarkError("locked accuracy scores are missing")
        scores_hash = sha256_file(scores_path)

    success_ids: dict[str, set[str]] = {}
    for row in scores:
        success_ids.setdefault(row.tool, set()).add(row.entry_id)
    if args.locked_role == "paired-sword-resources":
        success_ids = {
            "sword2-rust:legacy": {
                row.entry_id for row in runs if row.selector_variant == "legacy"
            },
            "sword2-rust:factorized": {
                row.entry_id for row in runs if row.selector_variant == "factorized"
            },
        }
    failure_ids: dict[str, set[str]] = {}
    for row in failures:
        failure_ids.setdefault(row.tool, set()).add(row.entry_id)
    selector_counts = {
        "legacy": sum(row.selector_variant == "legacy" for row in runs),
        "factorized": sum(row.selector_variant == "factorized" for row in runs),
    }
    fallback_count = sum(row.fallback is True for row in runs)
    exclusion_count = sum(
        row.excluded_candidate_count
        for row in runs
        if type(row.excluded_candidate_count) is int
    )
    expected_chainsaw = None
    if preflight.expected_chainsaw_file_sha256 is not None:
        expected_chainsaw = {
            "file_sha256": preflight.expected_chainsaw_file_sha256,
            "success_count": len(preflight.expected_chainsaw_ids),
            "success_id_set_sha256": canonical_id_set_hash(preflight.expected_chainsaw_ids),
            "failure_count": len(preflight.dataset_ids) - len(preflight.expected_chainsaw_ids),
            "failure_id_set_sha256": canonical_id_set_hash(
                set(preflight.dataset_ids) - set(preflight.expected_chainsaw_ids)
            ),
        }
    raw_root = Path(args.results_dir) / "raw"
    raw_evidence: dict[str, str] = {}
    if raw_root.exists():
        for path in sorted(raw_root.rglob("*")):
            if path.is_symlink():
                raise LockedBenchmarkError("locked raw evidence contains a symlink")
            if path.is_dir():
                continue
            if not path.is_file():
                raise LockedBenchmarkError("locked raw evidence contains a nonregular file")
            raw_evidence[path.relative_to(Path(args.results_dir)).as_posix()] = sha256_file(path)
    return {
        "schema_version": 1,
        "status": "complete",
        "created_at": created_at,
        "completed_at": completed_at,
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "locked_role": args.locked_role,
        "dataset": args.dataset,
        "dataset_sha256": preflight.dataset_sha256,
        "dataset_id_count": len(preflight.dataset_ids),
        "dataset_ids": list(preflight.dataset_ids),
        "dataset_id_set_sha256": preflight.dataset_id_set_sha256,
        "structure_sha256s": preflight.structure_sha256s,
        "structure_tree_sha256": preflight.structure_tree_sha256,
        "raw_evidence_sha256s": raw_evidence,
        "raw_evidence_tree_sha256": _length_framed_mapping_hash(raw_evidence),
        "tools": list(preflight.selected_tools),
        "sword2_extra_args": (
            [] if args.sword2_extra_args is None else args.sword2_extra_args.split()
        ),
        "sword2_threads": 1,
        "environment": LOCKED_ENVIRONMENT,
        "normalized_argv": _normalized_command(
            list(sys.argv),
            preflight=preflight,
            args=args,
        ),
        "gnu_time_version": preflight.gnu_time_version,
        "external_artifacts": preflight.external_artifacts,
        "expected_chainsaw": expected_chainsaw,
        "order_design": order,
        "model_manifest_sha256": preflight.model_manifest_sha256,
        "runtime_manifest_sha256": preflight.runtime_manifest_sha256,
        "binary_sha256": preflight.runtime_manifest["binary_sha256"],
        "runtime_source_git_commit": preflight.runtime_manifest["runtime_source_git_commit"],
        "runtime_input_tree_sha256": preflight.runtime_manifest["runtime_input_tree_sha256"],
        "evidence_tool_tree_sha256": preflight.runtime_manifest["evidence_tool_tree_sha256"],
        "cargo_lock_sha256": preflight.runtime_manifest["cargo_lock_sha256"],
        "row_counts": {
            "scores": len(scores),
            "runs": len(runs),
            "failures": len(failures),
        },
        "success_id_counts": {
            role: len(values) for role, values in sorted(success_ids.items())
        },
        "success_id_set_sha256s": {
            role: canonical_id_set_hash(values)
            for role, values in sorted(success_ids.items())
        },
        "failure_id_counts": {
            role: len(values) for role, values in sorted(failure_ids.items())
        },
        "failure_id_set_sha256s": {
            role: canonical_id_set_hash(values)
            for role, values in sorted(failure_ids.items())
        },
        "selector_counts": selector_counts,
        "fallback_count": fallback_count,
        "excluded_candidate_count": exclusion_count,
        "scores_sha256": scores_hash,
        "runs_sha256": sha256_file(runs_path),
        "failures_sha256": sha256_file(failures_path),
    }


def _run_locked_benchmark(args: argparse.Namespace) -> int:
    preflight = _locked_preflight(args)
    args.results_dir = Path(args.results_dir).absolute()
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=False, exist_ok=False)
    created_at = datetime.now(timezone.utc).isoformat()
    intent = _locked_intent_payload(preflight, args, created_at)
    _atomic_write_bytes(
        results_dir / "benchmark_intent.json",
        canonical_json_bytes(intent),
        absent_only=True,
    )
    try:
        if args.locked_role == "paired-sword-resources":
            scores, runs, failures, order = _run_locked_pairs(preflight, args)
        else:
            scores, runs, failures, order = _run_locked_accuracy(preflight, args)
        completed_at = datetime.now(timezone.utc).isoformat()
        manifest = _locked_manifest_payload(
            preflight,
            args,
            created_at=created_at,
            completed_at=completed_at,
            scores=scores,
            runs=runs,
            failures=failures,
            order=order,
        )
        _atomic_write_bytes(
            results_dir / "benchmark_manifest.json",
            canonical_json_bytes(manifest),
            absent_only=True,
        )
    except Exception as error:
        failure_path = results_dir / "failures.csv"
        if not failure_path.exists():
            failure = FailureRow(
                dataset=args.dataset,
                entry_id="",
                pdb_id="",
                chain_id="",
                tool="locked-harness",
                stage="infrastructure",
                message=type(error).__name__,
            )
            _atomic_write_rows(write_failures_csv, [failure], failure_path)
        raise
    return 0


def main() -> int:
    args = parse_args()
    locked_manifest_set = args.locked_factorized_manifest is not None
    locked_role_set = args.locked_role is not None
    if locked_manifest_set != locked_role_set:
        print(
            "locked-factorized-manifest and locked-role are required together",
            file=sys.stderr,
        )
        return 2
    if locked_manifest_set:
        try:
            return _run_locked_benchmark(args)
        except (OSError, ValueError) as error:
            print(f"locked benchmark evidence is invalid: {error}", file=sys.stderr)
            return 2
    if args.locked_artifact or args.chainsaw_expected_success_ids is not None:
        print("locked-only inputs require locked collection mode", file=sys.stderr)
        return 2
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
