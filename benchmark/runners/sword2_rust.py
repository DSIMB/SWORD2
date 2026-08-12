from __future__ import annotations

import json
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from benchmark.numbering import (
    StructureNumbering,
    map_sword2_chopping,
    read_sword2_residue_mapping,
)
from benchmark.runners.base import PartitionPrediction, ToolRunResult, run_timed, run_timed_locked


_ALT_RE = re.compile(r"Alternative partition (\d+)")


def _partition_sort_key(name: str) -> tuple[int, int]:
    if name == "Optimal partition":
        return (0, 0)
    match = _ALT_RE.fullmatch(name)
    if match:
        return (1, int(match.group(1)))
    return (2, 0)


def _domains_to_chopping(partition: dict) -> str:
    domains = partition.get("Domains", {})
    chunks: list[str] = []
    for domain_name in sorted(domains, key=lambda name: int(re.findall(r"\d+", name)[0])):
        pu_ranges = domains[domain_name].get("PUs", {})
        chunks.append("_".join(pu_ranges.keys()))
    return ",".join(chunks)


def parse_summary_partitions(
    summary: dict,
    numbering: StructureNumbering | None,
    residue_mapping: dict[int, object] | None = None,
    tool: str = "sword2-rust",
    raw_path: Path | None = None,
) -> list[PartitionPrediction]:
    if "Optimal partition" not in summary:
        raise KeyError("summary JSON does not contain 'Optimal partition'")

    partitions: list[PartitionPrediction] = []
    for name in sorted(
        [key for key, value in summary.items() if isinstance(value, dict) and "Domains" in value],
        key=_partition_sort_key,
    ):
        raw_chopping = _domains_to_chopping(summary[name])
        common_chopping = map_sword2_chopping(raw_chopping, numbering, residue_mapping)
        variant = "optimal" if name == "Optimal partition" else "alternative"
        partitions.append(
            PartitionPrediction(
                tool=tool,
                name=name,
                variant=variant,
                chopping=common_chopping,
                n_domains=int(summary[name].get("Nb. domains", len(common_chopping.split(",")))),
                raw_path=raw_path,
            )
        )
    return partitions


def load_summary_partitions(
    summary_path: Path,
    numbering: StructureNumbering | None,
    residue_mapping_path: Path | None = None,
    chain_id: str | None = None,
    tool: str = "sword2-rust",
) -> list[PartitionPrediction]:
    residue_mapping = (
        read_sword2_residue_mapping(residue_mapping_path, chain_id=chain_id)
        if residue_mapping_path is not None and residue_mapping_path.exists()
        else None
    )
    return parse_summary_partitions(
        json.loads(summary_path.read_text()),
        numbering=numbering,
        residue_mapping=residue_mapping,
        tool=tool,
        raw_path=summary_path,
    )


@dataclass(frozen=True)
class Sword2RustRunner:
    repo_dir: Path = Path("/home/chili/cretin/PROJECTS/SWORD2")
    binary: Path = Path("/home/chili/cretin/PROJECTS/SWORD2/target/release/sword2")
    threads: int | None = None
    experiments: str | None = None
    extra_args: str | None = None
    fresh_only: bool = False
    locked_env: dict[str, str] | None = None
    selector_status_path: Path | None = None

    def run(self, structure_file: Path, output_dir: Path) -> ToolRunResult:
        structure_file = structure_file.resolve()
        output_dir = output_dir.resolve()
        if output_dir.exists():
            if self.fresh_only:
                raise FileExistsError(f"locked SWORD output already exists: {output_dir}")
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        command = [str(self.binary.resolve()), "-i", str(structure_file), "-o", str(output_dir)]
        if self.threads:
            command.extend(["-j", str(self.threads)])
        if self.extra_args:
            command.extend(self.extra_args.split())
        environment: dict[str, str] = dict(self.locked_env or {})
        if self.selector_status_path is not None:
            status_path = self.selector_status_path.absolute()
            if self.fresh_only and status_path.parent != output_dir:
                raise ValueError("locked selector status must be inside its fresh process directory")
            if status_path.exists():
                raise FileExistsError(f"selector status already exists: {status_path}")
            environment["SWORD2_SELECTOR_STATUS"] = os.fspath(status_path)
        if self.fresh_only:
            if self.selector_status_path is None:
                raise ValueError("locked SWORD execution requires a selector-status path")
            return run_timed_locked(command, cwd=self.repo_dir, env=environment)
        if self.experiments:
            return run_timed(command, cwd=self.repo_dir, env={"SWORD2_EXPERIMENTS": self.experiments})
        if environment:
            return run_timed(command, cwd=self.repo_dir, env=environment)
        return run_timed(command, cwd=self.repo_dir)
