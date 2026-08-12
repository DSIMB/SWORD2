from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

import benchmark.evaluate_factorized_acceptance as acceptance_module
from benchmark.evaluate_factorized_acceptance import (
    BOOTSTRAP,
    GATE_ORDER,
    LOCKED_MANIFEST_KEYS,
    main,
    evaluate_accuracy_gates,
    evaluate_resource_gates,
    write_acceptance_report,
)
from benchmark.factorized_ranker.runtime_freeze import (
    canonical_id_set_hash,
    canonical_json_bytes,
    sha256_file,
)
from benchmark.run_benchmark import counterbalanced_pair_orders
from benchmark.score import (
    FailureRow,
    LockedRunRow,
    write_failures_csv,
    write_locked_runs_csv,
)


def _json_text(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _write_score_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "dataset",
        "entry_id",
        "pdb_id",
        "chain_id",
        "tool",
        "variant",
        "partition",
        "ndo",
        "n_pred_domains",
        "n_true_domains",
        "d_count_acc",
        "boundary_f1_10",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _mapping_hash(mapping: dict[str, str]) -> str:
    digest = hashlib.sha256()
    for name, value in sorted(mapping.items()):
        name_bytes = name.encode()
        value_bytes = value.encode()
        digest.update(len(name_bytes).to_bytes(8, "big"))
        digest.update(name_bytes)
        digest.update(len(value_bytes).to_bytes(8, "big"))
        digest.update(value_bytes)
    return digest.hexdigest()


def _synthetic_evidence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    ids = ["a", "b"]
    metadata = tmp_path / "CATH-test.csv"
    metadata.write_text(
        "1aaa,a,A,2,x,100,1-50:1.1.1.1|51-100:2.2.2.2\n"
        "2bbb,b,B,2,x,100,1-20_80-100:1.1.1.1|21-79:2.2.2.2\n"
    )
    chainsaw_ids = tmp_path / "chainsaw.txt"
    chainsaw_ids.write_text("a\n")
    standalone = tmp_path / "standalone.csv"
    standalone.write_text("chain_id,ndo\na,0.85\nb,0.85\n")
    model = tmp_path / "model.json"
    model.write_bytes(canonical_json_bytes({"model": "test"}))
    runtime = tmp_path / "runtime.json"
    runtime.write_bytes(
        canonical_json_bytes(
            {
                "runtime": "test",
                "selector_cache_contract": acceptance_module.CACHE_CONTRACT,
            }
        )
    )
    model_hash = sha256_file(model)
    runtime_hash = sha256_file(runtime)
    frozen_runtime = {
        "model_manifest_sha256": model_hash,
        "model_artifact_sha256s": {
            "standalone_baseline_sha256": sha256_file(standalone)
        },
        "binary_sha256": "b" * 64,
        "runtime_source_git_commit": "1" * 40,
        "runtime_input_tree_sha256": "c" * 64,
        "evidence_tool_tree_sha256": "d" * 64,
        "cargo_lock_sha256": "e" * 64,
        "selector_cache_contract": acceptance_module.CACHE_CONTRACT,
    }
    monkeypatch.setattr(acceptance_module, "LOCKED_POPULATION_SIZE", 2)
    monkeypatch.setattr(
        acceptance_module,
        "verify_top_level_manifest",
        lambda path: {"standalone_baseline_sha256": sha256_file(standalone)},
    )
    monkeypatch.setattr(
        acceptance_module,
        "verify_runtime_freeze",
        lambda path, **kwargs: frozen_runtime,
    )

    metadata_by_id = {
        "a": ("1aaa", "A"),
        "b": ("2bbb", "B"),
    }
    structures = {entry_id: (entry_id * 64)[:64] for entry_id in ids}
    common = {
        "schema_version": 1,
        "status": "complete",
        "created_at": "2026-01-01T00:00:00+00:00",
        "completed_at": "2026-01-01T00:01:00+00:00",
        "hostname": "test-host",
        "platform": "test-platform",
        "dataset": "cath663",
        "dataset_sha256": sha256_file(metadata),
        "dataset_id_count": 2,
        "dataset_ids": ids,
        "dataset_id_set_sha256": canonical_id_set_hash(ids),
        "structure_sha256s": structures,
        "structure_tree_sha256": _mapping_hash(structures),
        "sword2_threads": 1,
        "environment": {
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "RAYON_NUM_THREADS": "1",
        },
        "normalized_argv": ["benchmark.run_benchmark", "<locked>"],
        "gnu_time_version": "GNU time 1.9",
        "model_manifest_sha256": model_hash,
        "runtime_manifest_sha256": runtime_hash,
        "binary_sha256": "b" * 64,
        "runtime_source_git_commit": "1" * 40,
        "runtime_input_tree_sha256": "c" * 64,
        "evidence_tool_tree_sha256": "d" * 64,
        "cargo_lock_sha256": "e" * 64,
    }

    def add_raw(root: Path, relative: str, data: bytes) -> tuple[str, str]:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return relative, sha256_file(path)

    def sword_row(
        root: Path,
        entry_id: str,
        selector: str,
        role: str,
        pair_order: str = "",
        pair_position: int | str = "",
    ) -> LockedRunRow:
        process = f"raw/sword/{entry_id}/{selector}"
        status = canonical_json_bytes(
            {
                "error_code": None,
                "excluded_candidate_count": 0 if selector == "legacy" else 1,
                "fallback": False,
                "requested_selector": selector,
                "schema_version": 1,
                "selector_used": selector,
            }
        )
        summary_role, summary_hash = add_raw(
            root, f"{process}/output/summary.json", b"{}\n"
        )
        _, status_hash = add_raw(root, f"{process}/selector_status.json", status)
        _, stdout_hash = add_raw(root, f"{process}/stdout.log", b"")
        _, stderr_hash = add_raw(root, f"{process}/stderr.log", b"time\n")
        normalized = [
            "<binary>", "-i", "<input_structure>", "-o", "<process_output>", "-j", "1"
        ]
        if selector == "factorized":
            normalized.append("--use-factorized-ranker")
        pdb_id, chain_id = metadata_by_id[entry_id]
        return LockedRunRow(
            dataset="cath663",
            entry_id=entry_id,
            pdb_id=pdb_id,
            chain_id=chain_id,
            tool="sword2-rust",
            selector_variant=selector,
            locked_role=role,
            pair_order=pair_order,
            pair_position=pair_position,
            reused_output=False,
            returncode=0,
            runtime_s=1.0 if selector == "legacy" else 1.1,
            peak_rss_kb=100 if selector == "legacy" else 105,
            command_json=_json_text(["/tmp/sword2", "-i", "/tmp/in", "-o", "/tmp/out", "-j", "1"] + (["--use-factorized-ranker"] if selector == "factorized" else [])),
            normalized_command_json=_json_text(normalized),
            cwd_role="repository_root",
            input_structure_sha256=structures[entry_id],
            raw_summary_role=summary_role,
            raw_summary_sha256=summary_hash,
            stdout_sha256=stdout_hash,
            stderr_sha256=stderr_hash,
            selector_status_sha256=status_hash,
            requested_selector=selector,
            selector_used=selector,
            fallback=False,
            error_code="",
            selector_warning_code="",
            excluded_candidate_count=0 if selector == "legacy" else 1,
            binary_sha256="b" * 64,
            model_manifest_sha256=model_hash,
            runtime_manifest_sha256=runtime_hash,
        )

    def competitor_row(root: Path, entry_id: str, tool: str) -> LockedRunRow:
        raw_role, raw_hash = add_raw(root, f"raw/{tool}/{entry_id}.tsv", b"prediction\n")
        _, stdout_hash = add_raw(root, f"raw/{tool}/{entry_id}.stdout", b"")
        _, stderr_hash = add_raw(root, f"raw/{tool}/{entry_id}.stderr", b"")
        pdb_id, chain_id = metadata_by_id[entry_id]
        normalized_command = (
            [
                "<locked_artifact:merizo_python>",
                "predict.py",
                "-i",
                f"<input_structure:{entry_id}>",
                "-d",
                "cuda",
                "--return_indices",
                "--output_headers",
                "--pdb_chain",
                chain_id,
            ]
            if tool == "merizo"
            else [
                "<locked_artifact:chainsaw_python>",
                "get_predictions.py",
                "--structure_directory",
                "<results>/raw/chainsaw/_batch_stage",
                "--output",
                "<results>/raw/chainsaw/_batch.tsv",
            ]
        )
        return LockedRunRow(
            dataset="cath663", entry_id=entry_id, pdb_id=pdb_id, chain_id=chain_id,
            tool=tool, selector_variant="", locked_role="legacy-accuracy",
            pair_order="", pair_position="", reused_output=False, returncode=0,
            runtime_s=2.0, peak_rss_kb=200,
            command_json=_json_text(
                [f"/tmp/{tool}", *[f"arg{index}" for index in range(len(normalized_command) - 1)]]
            ),
            normalized_command_json=_json_text(normalized_command),
            cwd_role=f"locked_artifact:{tool}_source_tree",
            input_structure_sha256=structures[entry_id], raw_summary_role=raw_role,
            raw_summary_sha256=raw_hash, stdout_sha256=stdout_hash,
            stderr_sha256=stderr_hash, selector_status_sha256="",
            requested_selector="", selector_used="", fallback="", error_code="",
            selector_warning_code="", excluded_candidate_count="", binary_sha256="",
            model_manifest_sha256="", runtime_manifest_sha256="",
        )

    result: dict[str, Path] = {
        "dataset_metadata": metadata,
        "chainsaw_expected_success_ids": chainsaw_ids,
        "model_manifest": model,
        "runtime_manifest": runtime,
        "standalone_baseline": standalone,
    }
    for short, role in (
        ("legacy", "legacy-accuracy"),
        ("factorized", "factorized-accuracy"),
        ("resource", "paired-sword-resources"),
    ):
        root = tmp_path / short
        root.mkdir()
        runs: list[LockedRunRow] = []
        failures: list[FailureRow] = []
        order = None
        if short == "resource":
            assignments = counterbalanced_pair_orders(ids)
            assignment_rows = [[entry_id, first] for entry_id, first in assignments.items()]
            for entry_id, first in assignments.items():
                second = "factorized" if first == "legacy" else "legacy"
                runs.extend(
                    [
                        sword_row(root, entry_id, first, role, f"{first}_first", 1),
                        sword_row(root, entry_id, second, role, f"{first}_first", 2),
                    ]
                )
            order = {
                "seed": 37,
                "assignments": assignment_rows,
                "assignments_sha256": hashlib.sha256(canonical_json_bytes(assignment_rows)).hexdigest(),
                "legacy_first_count": sum(value == "legacy" for value in assignments.values()),
                "factorized_first_count": sum(value == "factorized" for value in assignments.values()),
            }
        else:
            selector = "legacy" if short == "legacy" else "factorized"
            runs.extend(sword_row(root, entry_id, selector, role) for entry_id in ids)
            if short == "legacy":
                runs.extend(competitor_row(root, entry_id, "merizo") for entry_id in ids)
                runs.append(competitor_row(root, "a", "chainsaw"))
                failures.append(
                    FailureRow(
                        dataset="cath663", entry_id="b", pdb_id="2bbb", chain_id="B",
                        tool="chainsaw", stage="expected_failure",
                        message="chainsaw_expected_failure", returncode=0,
                        command="[]", stderr=None,
                    )
                )
        runs_path = root / "runs.csv"
        failures_path = root / "failures.csv"
        write_locked_runs_csv(runs, runs_path)
        write_failures_csv(failures, failures_path)
        scores_path = root / "scores.csv"
        score_rows: list[dict[str, object]] = []
        if short != "resource":
            selector_tool_rows = ids
            for entry_id in selector_tool_rows:
                pdb_id, chain_id = metadata_by_id[entry_id]
                score_rows.append(
                    {"dataset": "cath663", "entry_id": entry_id, "pdb_id": pdb_id,
                     "chain_id": chain_id, "tool": "sword2-rust", "variant": "optimal",
                     "partition": "Optimal partition", "ndo": 0.9,
                     "n_pred_domains": 2, "n_true_domains": 2, "d_count_acc": 1.0,
                     "boundary_f1_10": 0.7}
                )
            if short == "legacy":
                for tool, tool_ids, ndo in (("merizo", ids, 0.2), ("chainsaw", ["a"], 0.2)):
                    for entry_id in tool_ids:
                        pdb_id, chain_id = metadata_by_id[entry_id]
                        score_rows.append(
                            {"dataset": "cath663", "entry_id": entry_id, "pdb_id": pdb_id,
                             "chain_id": chain_id, "tool": tool, "variant": "default",
                             "partition": tool, "ndo": ndo, "n_pred_domains": 2,
                             "n_true_domains": 2, "d_count_acc": 1.0,
                             "boundary_f1_10": 0.7}
                        )
            _write_score_rows(scores_path, score_rows)
            result[f"{short}_scores"] = scores_path

        raw_hashes = {
            path.relative_to(root).as_posix(): sha256_file(path)
            for path in sorted((root / "raw").rglob("*"))
            if path.is_file()
        }
        success = (
            {"sword2-rust:factorized": set(ids), "sword2-rust:legacy": set(ids)}
            if short == "resource"
            else {"sword2-rust": set(ids)}
        )
        if short == "legacy":
            success.update({"merizo": set(ids), "chainsaw": {"a"}})
        failure_sets = {"chainsaw": {"b"}} if short == "legacy" else {}
        artifacts = {}
        if short == "legacy":
            artifacts = {
                name: {"role": name, "kind": "file", "byte_count": 1,
                       "file_count": 1, "sha256": hashlib.sha256(name.encode()).hexdigest()}
                for name in {
                    "merizo_python", "merizo_source_tree", "merizo_model_tree",
                    "chainsaw_python", "chainsaw_source_tree", "chainsaw_model_tree",
                }
            }
        selector_counts = {
            "legacy": sum(row.selector_variant == "legacy" for row in runs),
            "factorized": sum(row.selector_variant == "factorized" for row in runs),
        }
        manifest = {
            **common,
            "locked_role": role,
            "normalized_argv": [
                "<repo>/benchmark/run_benchmark.py",
                "--dataset", "cath663",
                "--cache-dir", "<cache>",
                "--results-dir", "<results>",
                "--tools", "sword2-rust,merizo,chainsaw" if short == "legacy" else "sword2-rust",
                "--sword2-threads", "1",
                *(
                    ["--sword2-extra-args=--use-factorized-ranker"]
                    if short == "factorized" else []
                ),
                "--no-download", "--strict",
                "--locked-factorized-manifest", "<runtime_manifest>",
                "--locked-role", role,
                *(
                    [
                        "--chainsaw-expected-success-ids", "<chainsaw_expected_success_ids>",
                        *[
                            token
                            for artifact_role in sorted(
                                {
                                    "merizo_python", "merizo_source_tree", "merizo_model_tree",
                                    "chainsaw_python", "chainsaw_source_tree", "chainsaw_model_tree",
                                }
                            )
                            for token in (
                                "--locked-artifact",
                                f"{artifact_role}=<locked_artifact:{artifact_role}>",
                            )
                        ],
                    ]
                    if short == "legacy" else []
                ),
            ],
            "raw_evidence_sha256s": raw_hashes,
            "raw_evidence_tree_sha256": _mapping_hash(raw_hashes),
            "tools": ["sword2-rust", "merizo", "chainsaw"] if short == "legacy" else ["sword2-rust"],
            "sword2_extra_args": ["--use-factorized-ranker"] if short == "factorized" else [],
            "external_artifacts": artifacts,
            "expected_chainsaw": (
                {"file_sha256": sha256_file(chainsaw_ids), "success_count": 1,
                 "success_id_set_sha256": canonical_id_set_hash(["a"]), "failure_count": 1,
                 "failure_id_set_sha256": canonical_id_set_hash(["b"])}
                if short == "legacy" else None
            ),
            "order_design": order,
            "row_counts": {"scores": len(score_rows), "runs": len(runs), "failures": len(failures)},
            "success_id_counts": {name: len(values) for name, values in sorted(success.items())},
            "success_id_set_sha256s": {name: canonical_id_set_hash(values) for name, values in sorted(success.items())},
            "failure_id_counts": {name: len(values) for name, values in sorted(failure_sets.items())},
            "failure_id_set_sha256s": {name: canonical_id_set_hash(values) for name, values in sorted(failure_sets.items())},
            "selector_counts": selector_counts,
            "fallback_count": 0,
            "excluded_candidate_count": sum(
                int(row.excluded_candidate_count)
                for row in runs
                if type(row.excluded_candidate_count) is int
            ),
            "scores_sha256": sha256_file(scores_path) if short != "resource" else None,
            "runs_sha256": sha256_file(runs_path),
            "failures_sha256": sha256_file(failures_path),
        }
        assert set(manifest) == LOCKED_MANIFEST_KEYS
        manifest_path = root / "benchmark_manifest.json"
        manifest_path.write_bytes(canonical_json_bytes(manifest))
        result[f"{short}_manifest"] = manifest_path
        result[f"{short}_runs"] = runs_path
        result[f"{short}_failures"] = failures_path
    return result


def _evidence_cli_args(paths: dict[str, Path]) -> list[str]:
    result: list[str] = []
    for name in acceptance_module.EVIDENCE_ARGUMENT_NAMES:
        result.extend(["--" + name.replace("_", "-"), str(paths[name])])
    return result


def _accuracy_frames(ndo: float = 0.9):
    ids = ["a", "b", "c", "d"]
    factorized = pd.DataFrame(
        {
            "entry_id": ids,
            "ndo": [ndo] * 4,
            "d_count_acc": [1.0, 1.0, 1.0, 0.0],
            "n_pred_domains": [2, 2, 3, 2],
            "n_true_domains": [2, 2, 3, 3],
            "boundary_f1_10": [0.7] * 4,
            "continuity_cohort": ["contiguous", "contiguous", "discontinuous", "discontinuous"],
        }
    ).set_index("entry_id")
    standalone = pd.DataFrame(
        {"entry_id": ids, "ndo": [0.9, 0.9, 0.9, 0.9]}
    ).set_index("entry_id")
    merizo = pd.DataFrame(
        {"entry_id": ids, "ndo": [0.4, 0.4, 0.4, 0.4]}
    ).set_index("entry_id")
    chainsaw = pd.DataFrame(
        {"entry_id": ["a", "c"], "ndo": [0.3, 0.3]}
    ).set_index("entry_id")
    return factorized, standalone, {"merizo": merizo, "chainsaw": chainsaw}


def test_accuracy_gates_use_exact_signs_cohorts_and_thresholds():
    factorized, standalone, competitors = _accuracy_frames()
    result = evaluate_accuracy_gates(factorized, standalone, competitors)
    assert list(result["gates"]) == list(GATE_ORDER[:7])
    assert result["gates"]["overall_ndo"]["passed"] is True
    assert result["gates"]["merizo_ndo_ci_low"]["value"] > 0.0
    assert result["gates"]["chainsaw_ndo_ci_low"]["denominator"] == 2
    assert result["gates"]["domain_count_accuracy"]["value"] == pytest.approx(0.75)
    assert result["gates"]["domain_count_accuracy"]["passed"] is True
    assert result["gates"]["contiguous_standalone_ndo_delta"]["value"] == pytest.approx(0.0)
    assert result["gates"]["discontinuous_standalone_ndo_delta"]["value"] == pytest.approx(0.0)

    factorized, standalone, competitors = _accuracy_frames(ndo=0.8389)
    equality = evaluate_accuracy_gates(factorized, standalone, competitors)
    assert equality["gates"]["overall_ndo"]["passed"] is False


def test_resource_gates_use_median_runtime_ratio_and_ratio_of_rss_maxima():
    legacy = pd.DataFrame(
        {"entry_id": ["a", "b", "c"], "runtime_s": [1.0, 10.0, 10.0], "peak_rss_kb": [100, 100, 100]}
    ).set_index("entry_id")
    factorized = pd.DataFrame(
        {"entry_id": ["a", "b", "c"], "runtime_s": [10.0, 10.0, 10.0], "peak_rss_kb": [110, 90, 90]}
    ).set_index("entry_id")
    result = evaluate_resource_gates(legacy, factorized)
    assert list(result["gates"]) == list(GATE_ORDER[7:])
    assert result["gates"]["runtime_median_ratio"]["value"] == pytest.approx(1.0)
    assert result["gates"]["runtime_median_ratio"]["passed"] is True
    assert result["gates"]["rss_maxima_ratio"]["value"] == pytest.approx(1.1)
    assert result["gates"]["rss_maxima_ratio"]["passed"] is True


def test_gate_threshold_inclusivity_and_finite_inputs_are_fail_closed():
    factorized, standalone, competitors = _accuracy_frames(ndo=0.9)
    factorized["boundary_f1_10"] = 0.620
    factorized.loc["d", "d_count_acc"] = 1.0
    factorized.loc["d", "n_pred_domains"] = 3
    factorized["ndo"] = 0.005
    standalone["ndo"] = 0.01
    competitors["merizo"]["ndo"] = factorized["ndo"]
    competitors["chainsaw"]["ndo"] = factorized.loc[["a", "c"], "ndo"]
    result = evaluate_accuracy_gates(factorized, standalone, competitors)
    assert result["gates"]["domain_count_accuracy"]["passed"] is True
    assert result["gates"]["boundary_f1_10"]["passed"] is True
    assert result["gates"]["contiguous_standalone_ndo_delta"]["passed"] is True
    assert result["gates"]["discontinuous_standalone_ndo_delta"]["passed"] is True
    assert result["gates"]["merizo_ndo_ci_low"]["value"] == pytest.approx(0.0)
    assert result["gates"]["merizo_ndo_ci_low"]["passed"] is False
    assert result["gates"]["chainsaw_ndo_ci_low"]["passed"] is False

    legacy = pd.DataFrame(
        {"entry_id": ["a", "b", "c"], "runtime_s": [1.0, 1.0, 1.0], "peak_rss_kb": [100, 90, 80]}
    ).set_index("entry_id")
    factorized_resources = pd.DataFrame(
        {"entry_id": ["a", "b", "c"], "runtime_s": [1.15, 1.15, 1.15], "peak_rss_kb": [110, 99, 88]}
    ).set_index("entry_id")
    resources = evaluate_resource_gates(legacy, factorized_resources)
    assert resources["gates"]["runtime_median_ratio"]["value"] == pytest.approx(1.15)
    assert resources["gates"]["runtime_median_ratio"]["passed"] is True
    assert resources["gates"]["rss_maxima_ratio"]["passed"] is True

    factorized.loc["a", "ndo"] = float("nan")
    with pytest.raises(ValueError, match="nonfinite"):
        evaluate_accuracy_gates(factorized, standalone, competitors)


def test_acceptance_reports_are_canonical_deterministic_and_absent_only(tmp_path: Path):
    evidence = {
        name: hashlib.sha256(name.encode()).hexdigest()
        for name in acceptance_module.EVIDENCE_ARGUMENT_NAMES
    }
    rules = {
        "overall_ndo": (0.9, 0.8389, ">", "overall", "arithmetic_mean"),
        "merizo_ndo_ci_low": (0.1, 0.0, ">", "overall", "paired_chain_bootstrap_ci_low"),
        "chainsaw_ndo_ci_low": (0.1, 0.0, ">", "chainsaw", "paired_chain_bootstrap_ci_low"),
        "domain_count_accuracy": (0.8, 0.745, ">=", "overall", "arithmetic_mean"),
        "boundary_f1_10": (0.7, 0.620, ">=", "overall", "arithmetic_mean"),
        "contiguous_standalone_ndo_delta": (0.0, -0.005, ">=", "contiguous", "paired_cohort_mean_delta"),
        "discontinuous_standalone_ndo_delta": (0.0, -0.005, ">=", "discontinuous", "paired_cohort_mean_delta"),
        "runtime_median_ratio": (1.0, 1.15, "<=", "resources", "median_of_paired_ratios"),
        "rss_maxima_ratio": (1.0, 1.10, "<=", "resources", "ratio_of_selector_maxima"),
    }
    denominators = {
        "overall": 4,
        "chainsaw": 2,
        "contiguous": 2,
        "discontinuous": 2,
        "resources": 4,
    }
    gates = {
        name: {
            "name": name,
            "measured": True,
            "value": value,
            "threshold": threshold,
            "comparison": comparison,
            "denominator": denominators[denominator_role],
            "method": method,
            "passed": True,
            "source_sha256s": {"factorized_scores": evidence["factorized_scores"]},
        }
        for name, (value, threshold, comparison, denominator_role, method) in rules.items()
    }
    results = {
        "schema_version": 1,
        "gate_order": list(GATE_ORDER),
        "bootstrap": BOOTSTRAP,
        "model_manifest_sha256": "a" * 64,
        "runtime_manifest_sha256": "b" * 64,
        "coverage_attestation_sha256": "c" * 64,
        "evidence_sha256s": evidence,
        "denominators": denominators,
        "gates": gates,
        "diagnostics": {},
        "all_gates_measured": True,
        "all_gates_pass": True,
    }
    json_path = tmp_path / "a" / "acceptance.json"
    markdown_path = tmp_path / "a" / "acceptance.md"
    write_acceptance_report(json_path, markdown_path, results)
    assert json_path.read_bytes().endswith(b"\n")
    assert json.loads(json_path.read_bytes()) == results
    assert markdown_path.read_bytes().endswith(b"\n")

    with pytest.raises(FileExistsError):
        write_acceptance_report(json_path, tmp_path / "other.md", results)


def test_synthetic_coverage_precedes_metrics_and_evaluation_is_deterministic(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    paths = _synthetic_evidence(tmp_path, monkeypatch)
    coverage = tmp_path / "coverage.json"
    assert main(["coverage", *_evidence_cli_args(paths), "--coverage-out", str(coverage)]) == 0
    coverage_payload = json.loads(coverage.read_bytes())
    assert coverage_payload["valid"] is True
    assert coverage_payload["dataset_id_count"] == 2
    assert not ({"ndo", "boundary_f1_10", "runtime_median_ratio"} & set(coverage_payload))

    first_json = tmp_path / "first.json"
    first_markdown = tmp_path / "first.md"
    evaluate_args = [
        "evaluate",
        *_evidence_cli_args(paths),
        "--coverage-attestation", str(coverage),
        "--json-out", str(first_json),
        "--markdown-out", str(first_markdown),
        "--bootstrap-replicates", str(BOOTSTRAP["replicates"]),
        "--seed", str(BOOTSTRAP["seed"]),
    ]
    assert main(evaluate_args) == 0
    payload = json.loads(first_json.read_bytes())
    assert payload["gate_order"] == list(GATE_ORDER)
    assert set(payload["gates"]) == set(GATE_ORDER)
    assert payload["all_gates_measured"] is True
    assert payload["all_gates_pass"] is True
    assert payload["denominators"] == {
        "overall": 2,
        "chainsaw": 1,
        "contiguous": 1,
        "discontinuous": 1,
        "resources": 2,
    }

    second_json = tmp_path / "second.json"
    second_markdown = tmp_path / "second.md"
    evaluate_args[evaluate_args.index(str(first_json))] = str(second_json)
    evaluate_args[evaluate_args.index(str(first_markdown))] = str(second_markdown)
    assert main(evaluate_args) == 0
    assert second_json.read_bytes() == first_json.read_bytes()
    assert second_markdown.read_bytes() == first_markdown.read_bytes()


def test_coverage_rejects_duplicate_rank_one_and_tampered_raw_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    paths = _synthetic_evidence(tmp_path, monkeypatch)
    factorized_scores = paths["factorized_scores"]
    rows = list(csv.DictReader(factorized_scores.open(newline="")))
    _write_score_rows(factorized_scores, [*rows, rows[0]])
    manifest_path = paths["factorized_manifest"]
    manifest = json.loads(manifest_path.read_bytes())
    manifest["scores_sha256"] = sha256_file(factorized_scores)
    manifest["row_counts"]["scores"] += 1
    manifest_path.write_bytes(canonical_json_bytes(manifest))
    with pytest.raises(acceptance_module.InvalidEvidence, match="duplicated or incomplete"):
        main(
            [
                "coverage",
                *_evidence_cli_args(paths),
                "--coverage-out",
                str(tmp_path / "invalid.json"),
            ]
        )
    assert not (tmp_path / "invalid.json").exists()

    paths = _synthetic_evidence(tmp_path / "other", monkeypatch)
    raw = next((paths["resource_manifest"].parent / "raw").rglob("summary.json"))
    raw.write_text("tampered\n")
    with pytest.raises(acceptance_module.InvalidEvidence, match="raw evidence"):
        main(
            [
                "coverage",
                *_evidence_cli_args(paths),
                "--coverage-out",
                str(tmp_path / "tampered.json"),
            ]
        )


def test_evaluation_revalidates_coverage_before_loading_metrics(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    paths = _synthetic_evidence(tmp_path, monkeypatch)
    coverage = tmp_path / "coverage.json"
    assert main(["coverage", *_evidence_cli_args(paths), "--coverage-out", str(coverage)]) == 0
    paths["factorized_scores"].write_text("metric bytes changed\n")
    called = False

    def forbidden_read_csv(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("metrics were read before coverage validation")

    monkeypatch.setattr(acceptance_module.pd, "read_csv", forbidden_read_csv)
    with pytest.raises(acceptance_module.InvalidEvidence):
        main(
            [
                "evaluate",
                *_evidence_cli_args(paths),
                "--coverage-attestation", str(coverage),
                "--json-out", str(tmp_path / "acceptance.json"),
                "--markdown-out", str(tmp_path / "acceptance.md"),
                "--bootstrap-replicates", "10000",
                "--seed", "37",
            ]
        )
    assert called is False
    assert not (tmp_path / "acceptance.json").exists()


def test_promotion_rederives_committed_bytes_and_keeps_nonpromotable_cache_opt_in(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
):
    repo = tmp_path / "repo"
    paths = _synthetic_evidence(repo / "evidence", monkeypatch)
    coverage = repo / "coverage.json"
    acceptance = repo / "acceptance.json"
    markdown = repo / "acceptance.md"
    assert main(["coverage", *_evidence_cli_args(paths), "--coverage-out", str(coverage)]) == 0
    assert (
        main(
            [
                "evaluate",
                *_evidence_cli_args(paths),
                "--coverage-attestation", str(coverage),
                "--json-out", str(acceptance),
                "--markdown-out", str(markdown),
                "--bootstrap-replicates", "10000",
                "--seed", "37",
            ]
        )
        == 0
    )
    receipt = repo / "benchmark/models/factorized_ranker_v1_execution_receipt.json"
    report = repo / "benchmark/REPORT.md"
    receipt.parent.mkdir(parents=True)
    receipt.write_bytes(canonical_json_bytes({"receipt": "synthetic"}))
    report.write_text("synthetic report\n")
    subprocess.run(["git", "init", "-q", repo], check=True)
    subprocess.run(["git", "-C", repo, "config", "user.email", "test@example.com"], check=True)
    subprocess.run(["git", "-C", repo, "config", "user.name", "Test"], check=True)
    subprocess.run(["git", "-C", repo, "add", "."], check=True)
    subprocess.run(["git", "-C", repo, "commit", "-qm", "synthetic locked evidence"], check=True)
    commit = subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"], text=True).strip()
    monkeypatch.chdir(repo)
    exit_code = main(
        [
            "validate-promotion",
            *_evidence_cli_args(paths),
            "--coverage-attestation", str(coverage),
            "--acceptance", str(acceptance),
            "--markdown", str(markdown),
            "--task17-commit", commit,
        ]
    )
    assert exit_code == 1
    assert capsys.readouterr().out == "KEEP_OPT_IN cache_context_not_promotable\n"

    payload = json.loads(acceptance.read_bytes())
    payload["all_gates_pass"] = False
    acceptance.write_bytes(canonical_json_bytes(payload))
    with pytest.raises(acceptance_module.InvalidEvidence):
        main(
            [
                "validate-promotion",
                *_evidence_cli_args(paths),
                "--coverage-attestation", str(coverage),
                "--acceptance", str(acceptance),
                "--markdown", str(markdown),
                "--task17-commit", commit,
            ]
        )
