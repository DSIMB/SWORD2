from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from pathlib import Path

from benchmark.numbering import (
    StructureNumbering,
    _split_segment,
    map_author_chopping,
    map_one_based_chopping,
    split_domains,
)
from benchmark.runners.base import ParsedPrediction, ToolRunResult, ensure_cuda_available, run_timed


DEFAULT_CHAINSAW_DIR = Path("/home/chili/cretin/PROJECTS/chainsaw")


def _write_tsv_rows_by_key(text: str, key_column: str, outputs_by_key: dict[str, Path]) -> None:
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return
    header = lines[0]
    columns = header.split("\t")
    try:
        key_index = columns.index(key_column)
    except ValueError as exc:
        raise ValueError(f"TSV header has no {key_column!r} column: {header!r}") from exc

    rows_by_key: dict[str, list[str]] = {key: [] for key in outputs_by_key}
    for line in lines[1:]:
        if line == header:
            continue
        values = line.split("\t")
        if len(values) <= key_index:
            continue
        key = values[key_index]
        if key in rows_by_key:
            rows_by_key[key].append(line)

    for key, output_path in outputs_by_key.items():
        rows = rows_by_key.get(key, [])
        if not rows:
            continue
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join([header, *rows]) + "\n")


def parse_chainsaw_tsv(text: str) -> ParsedPrediction:
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        raise ValueError("Chainsaw TSV contains no prediction rows")
    row = rows[0]
    chopping = row.get("chopping") or row.get("pred_chopping")
    if not chopping:
        raise ValueError(f"Chainsaw TSV row has no chopping column: {row!r}")
    return ParsedPrediction(chopping=chopping, n_domains=chopping.count("|") + 1, raw=row)


def _chopping_within_length(chopping: str, n_residues: int) -> bool:
    for domain in split_domains(chopping):
        for segment in domain:
            start, end = _split_segment(segment)
            if int(start) < 0 or int(end) >= n_residues:
                return False
    return True


def chainsaw_to_common_chopping(
    chopping: str,
    numbering: StructureNumbering | None = None,
    chain_id: str | None = None,
) -> str:
    """Chainsaw output is documented as 1-based sequential; convert to 0-based."""
    sequential = map_one_based_chopping(chopping)
    if numbering is None or _chopping_within_length(sequential, numbering.n_residues):
        return sequential
    return map_author_chopping(chopping, numbering=numbering, chain_id=chain_id)


@dataclass(frozen=True)
class ChainsawRunner:
    chainsaw_dir: Path = DEFAULT_CHAINSAW_DIR
    python: Path = DEFAULT_CHAINSAW_DIR / "chswEnv/bin/python"
    require_cuda: bool = False

    def _ensure_cuda_if_required(self) -> None:
        if self.require_cuda:
            ensure_cuda_available(self.python, self.chainsaw_dir)

    def run_batch(
        self,
        inputs: list[tuple[Path, Path]],
        stage_dir: Path,
        combined_output: Path,
    ) -> ToolRunResult:
        if not inputs:
            raise ValueError("Chainsaw batch requires at least one input")

        stage_dir = stage_dir.resolve()
        combined_output = combined_output.resolve()
        self._ensure_cuda_if_required()
        if stage_dir.exists():
            shutil.rmtree(stage_dir)
        stage_dir.mkdir(parents=True)
        combined_output.parent.mkdir(parents=True, exist_ok=True)
        if combined_output.exists():
            combined_output.unlink()

        resolved_inputs = [(structure_file.resolve(), output_tsv.resolve()) for structure_file, output_tsv in inputs]
        outputs_by_chain_id: dict[str, Path] = {}
        for structure_file, output_tsv in resolved_inputs:
            output_tsv.parent.mkdir(parents=True, exist_ok=True)
            if output_tsv.exists():
                output_tsv.unlink()
            chain_id = output_tsv.stem
            staged_path = stage_dir / f"{chain_id}{structure_file.suffix}"
            shutil.copy2(structure_file, staged_path)
            outputs_by_chain_id[chain_id] = output_tsv

        command = [
            str(self.python.absolute()),
            "get_predictions.py",
            "--structure_directory",
            str(stage_dir),
            "--output",
            str(combined_output),
        ]
        result = run_timed(command, cwd=self.chainsaw_dir)
        if combined_output.exists():
            _write_tsv_rows_by_key(combined_output.read_text(), "chain_id", outputs_by_chain_id)
        missing_outputs = [output_tsv for _, output_tsv in resolved_inputs if not output_tsv.exists()]
        if result.returncode == 0 and missing_outputs:
            missing = ", ".join(str(path) for path in missing_outputs)
            raise RuntimeError(f"Chainsaw completed but produced no TSV output for: {missing}")
        return result

    def run(self, structure_file: Path, output_tsv: Path) -> ToolRunResult:
        structure_file = structure_file.resolve()
        output_tsv = output_tsv.resolve()
        self._ensure_cuda_if_required()
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        if output_tsv.exists():
            output_tsv.unlink()
        command = [
            str(self.python.absolute()),
            "get_predictions.py",
            "--structure_file",
            str(structure_file),
            "--output",
            str(output_tsv),
        ]
        result = run_timed(command, cwd=self.chainsaw_dir)
        if result.returncode == 0 and not output_tsv.exists():
            raise RuntimeError(
                "Chainsaw completed but produced no TSV output; "
                f"check that the input path is visible to Chainsaw: {structure_file}"
            )
        return result
