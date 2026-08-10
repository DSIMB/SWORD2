from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from benchmark.numbering import StructureNumbering
from benchmark.runners.base import PartitionPrediction, ToolRunResult, run_timed
from benchmark.runners.sword2_rust import load_summary_partitions


DEFAULT_SWORD2_ORIG_DIR = Path("/home/chili/cretin/PROJECTS/SWORD2_ORIGINAL")
DEFAULT_SWORD2_CONDA_ENV = Path("/opt/apps/pkgs/anaconda3/2023.09-0/kvhhznl/envs/sword2")


def load_original_summary_partitions(
    summary_path: Path,
    numbering: StructureNumbering | None,
    residue_mapping_path: Path | None = None,
    chain_id: str | None = None,
) -> list[PartitionPrediction]:
    return load_summary_partitions(
        summary_path=summary_path,
        numbering=numbering,
        residue_mapping_path=residue_mapping_path,
        chain_id=chain_id,
        tool="sword2-original",
    )


@dataclass(frozen=True)
class Sword2OriginalRunner:
    repo_dir: Path = DEFAULT_SWORD2_ORIG_DIR
    conda_env: Path = DEFAULT_SWORD2_CONDA_ENV

    def run(
        self,
        pdb_id: str,
        output_dir: Path,
        structure_file: Path | None = None,
        chain_id: str | None = None,
    ) -> ToolRunResult:
        output_dir = output_dir.resolve()
        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        command = [
            "conda",
            "run",
            "-p",
            str(self.conda_env),
            "./SWORD2.py",
            "-o",
            str(output_dir),
        ]
        if structure_file is not None:
            command.extend(["-i", str(structure_file.resolve())])
            if chain_id:
                command.extend(["-c", chain_id])
        else:
            command.extend(["-p", pdb_id])
        return run_timed(command, cwd=self.repo_dir)
