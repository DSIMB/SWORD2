from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from benchmark.numbering import StructureNumbering, map_author_chopping
from benchmark.runners.base import ParsedPrediction, ToolRunResult, ensure_cuda_available, run_timed


DEFAULT_MERIZO_DIR = Path("/home/chili/cretin/PROJECTS/Merizo")


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
        key = Path(values[key_index]).name
        if key in rows_by_key:
            rows_by_key[key].append(line)

    for key, output_path in outputs_by_key.items():
        rows = rows_by_key.get(key, [])
        if not rows:
            continue
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join([header, *rows]) + "\n")


def parse_merizo_tsv(text: str) -> ParsedPrediction:
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        raise ValueError("Merizo TSV contains no prediction rows")
    row = rows[0]
    chopping = row.get("result") or row.get("chopping")
    if not chopping:
        raise ValueError(f"Merizo TSV row has no result/chopping column: {row!r}")
    n_domains = int(row["ndom"]) if row.get("ndom") else None
    return ParsedPrediction(chopping=chopping, n_domains=n_domains, raw=row)


def merizo_to_common_chopping(
    chopping: str,
    numbering: StructureNumbering,
    chain_id: str | None = None,
) -> str:
    return map_author_chopping(chopping, numbering=numbering, chain_id=chain_id)


@dataclass(frozen=True)
class MerizoRunner:
    merizo_dir: Path = DEFAULT_MERIZO_DIR
    python: Path = DEFAULT_MERIZO_DIR / "merizo/bin/python"
    device: str = "cpu"
    require_cuda: bool = True

    def _ensure_device_available(self) -> None:
        if self.device == "cuda" and self.require_cuda:
            ensure_cuda_available(self.python, self.merizo_dir)

    def run_batch(self, inputs: list[tuple[Path, Path, str]]) -> ToolRunResult:
        if not inputs:
            raise ValueError("Merizo batch requires at least one input")

        resolved_inputs = [(structure_file.resolve(), output_tsv.resolve(), chain_id) for structure_file, output_tsv, chain_id in inputs]
        chain_ids = {chain_id for _, _, chain_id in resolved_inputs}
        if len(chain_ids) != 1:
            raise ValueError("Merizo batch can only contain one chain id because predict.py accepts one --pdb_chain")

        self._ensure_device_available()

        for _, output_tsv, _ in resolved_inputs:
            output_tsv.parent.mkdir(parents=True, exist_ok=True)
            if output_tsv.exists():
                output_tsv.unlink()

        chain_id = next(iter(chain_ids))
        command = [
            str(self.python.absolute()),
            "predict.py",
            "-i",
            *[str(structure_file) for structure_file, _, _ in resolved_inputs],
            "-d",
            self.device,
            "--return_indices",
            "--output_headers",
            "--pdb_chain",
            chain_id,
        ]
        result = run_timed(command, cwd=self.merizo_dir)
        if result.stdout.strip():
            outputs_by_input = {
                structure_file.name: output_tsv for structure_file, output_tsv, _ in resolved_inputs
            }
            _write_tsv_rows_by_key(result.stdout, "input", outputs_by_input)
        missing_outputs = [output_tsv for _, output_tsv, _ in resolved_inputs if not output_tsv.exists()]
        if result.returncode == 0 and missing_outputs:
            missing = ", ".join(str(path) for path in missing_outputs)
            raise RuntimeError(f"Merizo completed but produced no TSV output for: {missing}")
        return result

    def run(self, structure_file: Path, output_tsv: Path, chain_id: str = "A") -> ToolRunResult:
        structure_file = structure_file.resolve()
        output_tsv = output_tsv.resolve()
        self._ensure_device_available()
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        if output_tsv.exists():
            output_tsv.unlink()
        command = [
            str(self.python.absolute()),
            "predict.py",
            "-i",
            str(structure_file),
            "-d",
            self.device,
            "--return_indices",
            "--output_headers",
            "--pdb_chain",
            chain_id,
        ]
        result = run_timed(command, cwd=self.merizo_dir)
        if result.stdout.strip():
            output_tsv.write_text(result.stdout)
        elif result.returncode == 0 and not output_tsv.exists():
            raise RuntimeError(
                "Merizo completed but produced no TSV output; "
                f"check that the input path is visible to Merizo: {structure_file}"
            )
        return result
