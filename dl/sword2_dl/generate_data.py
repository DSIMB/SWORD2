"""
Data generation pipeline: run SWORD2 on SwissProt proteins.

Steps:
1. Download SwissProt FASTA from UniProt
2. Parse sequences
3. Run SWORD2 on each protein (parallelized)
4. Parse SWORD2 JSON output
5. Convert to training format
6. Split into train/val/test
"""

import json
import logging
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
from tqdm import tqdm

logger = logging.getLogger(__name__)

RESULTS_SUBDIR = "results"
FAILURES_LOG = "failures.jsonl"


def download_swissprot(output_path: str) -> str:
    """Download SwissProt FASTA from UniProt."""
    import urllib.request

    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gz_path = output_path + ".gz"

    if Path(output_path).exists():
        logger.info(f"SwissProt FASTA already exists at {output_path}")
        return output_path

    logger.info(f"Downloading SwissProt from {url}...")
    urllib.request.urlretrieve(url, gz_path)

    logger.info("Decompressing...")
    import gzip

    with gzip.open(gz_path, "rt") as f_in, open(output_path, "w") as f_out:
        f_out.write(f_in.read())

    os.remove(gz_path)
    logger.info(f"SwissProt FASTA saved to {output_path}")
    return output_path


def parse_fasta(fasta_path: str) -> list[dict]:
    """Parse FASTA file into list of {id, sequence, description}."""
    sequences = []
    current_id = None
    current_desc = ""
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences.append(
                        {
                            "id": current_id,
                            "description": current_desc,
                            "sequence": "".join(current_seq),
                        }
                    )

                # Parse header: >sp|P12345|PROT_HUMAN Description OS=...
                header = line[1:]
                parts = header.split("|")
                if len(parts) >= 3:
                    current_id = parts[1]  # UniProt accession
                    current_desc = parts[2].split(" OS=")[0].strip()
                else:
                    current_id = header.split()[0]
                    current_desc = header

                current_seq = []
            elif line:
                current_seq.append(line)

    if current_id is not None:
        sequences.append(
            {
                "id": current_id,
                "description": current_desc,
                "sequence": "".join(current_seq),
            }
        )

    return sequences


def save_result_atomic(result: dict, results_dir: str) -> None:
    """Write a protein result JSON atomically via temp file + os.replace()."""
    os.makedirs(results_dir, exist_ok=True)
    protein_id = result["id"]
    final_path = os.path.join(results_dir, f"{protein_id}.json")
    tmp_path = os.path.join(results_dir, f".{protein_id}.json.tmp")
    with open(tmp_path, "w") as f:
        json.dump(result, f)
    os.replace(tmp_path, final_path)


def log_failure(protein_id: str, reason: str, output_dir: str) -> None:
    """Append a failure entry to the JSONL log."""
    log_path = os.path.join(output_dir, FAILURES_LOG)
    entry = json.dumps({
        "id": protein_id,
        "reason": reason,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    })
    with open(log_path, "a") as f:
        f.write(entry + "\n")


def get_completed_ids(output_dir: str) -> set[str]:
    """Return IDs of proteins that have a completed result JSON."""
    results_dir = os.path.join(output_dir, RESULTS_SUBDIR)
    if not os.path.isdir(results_dir):
        return set()
    return {
        p.stem for p in Path(results_dir).glob("*.json")
        if not p.name.startswith(".")
    }


def cleanup_incomplete_runs(output_dir: str, protein_ids: set[str]) -> int:
    """Remove intermediate dirs for proteins that lack a result JSON.

    These are from interrupted SWORD2 runs. Returns count of cleaned dirs.
    """
    completed = get_completed_ids(output_dir)
    cleaned = 0
    for pid in protein_ids:
        intermediate_dir = os.path.join(output_dir, pid)
        if os.path.isdir(intermediate_dir) and pid not in completed:
            shutil.rmtree(intermediate_dir, ignore_errors=True)
            cleaned += 1
    return cleaned


def run_sword2_on_sequence(
    protein_id: str,
    sequence: str,
    sword2_binary: str,
    base_dir: str,
    output_dir: str,
    timeout: int = 120,
    pdb_dir: Optional[str] = None,
) -> Optional[dict]:
    """Run SWORD2 on a single protein sequence.

    Uses a local PDB file if pdb_dir is provided and the file exists,
    otherwise fetches from AlphaFold via UniProt ID.

    Returns parsed result dict or None on failure.
    """
    result_dir = os.path.join(output_dir, protein_id)

    results_dir = os.path.join(output_dir, RESULTS_SUBDIR)

    try:
        # Prefer local PDB file over network fetch
        if pdb_dir:
            local_pdb = os.path.join(pdb_dir, f"AF-{protein_id}-F1-model_v4.pdb")
            input_args = ["-i", local_pdb] if os.path.exists(local_pdb) else ["-u", protein_id]
        else:
            input_args = ["-u", protein_id]

        cmd = [
            sword2_binary,
            *input_args,
            "-o", result_dir,
            "-e",  # disable energies for speed
            "-l",  # disable plots for speed
            "-q",  # quiet
            "--base-dir", base_dir,
        ]

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        if proc.returncode != 0:
            reason = f"nonzero exit {proc.returncode}: {proc.stderr[:200]}"
            logger.debug(f"SWORD2 failed for {protein_id}: {reason}")
            log_failure(protein_id, reason, output_dir)
            return None

        # Find and parse the JSON output
        result = parse_sword2_output(protein_id, sequence, result_dir)
        if result is None:
            log_failure(protein_id, "no parseable partitionings in SWORD2 output", output_dir)
            return None

        save_result_atomic(result, results_dir)
        return result

    except subprocess.TimeoutExpired:
        logger.debug(f"SWORD2 timed out for {protein_id}")
        log_failure(protein_id, f"timeout after {timeout}s", output_dir)
        return None
    except Exception as e:
        logger.debug(f"Error processing {protein_id}: {e}")
        log_failure(protein_id, str(e), output_dir)
        return None


def parse_contact_matrix(result_dir: str) -> Optional[list[list[int]]]:
    """Parse SWORD2 contact matrix into a list of contacting residue pairs.

    Reads the contact_matrix.mat file produced by SWORD2's peeling step.
    Returns a list of [i, j] pairs (0-indexed) where C-alpha distance < 8A.
    """
    mat_files = list(Path(result_dir).rglob("contact_matrix.mat"))
    if not mat_files:
        return None

    contacts = []
    try:
        with open(mat_files[0]) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    i = int(parts[0]) - 1  # convert to 0-indexed
                    j = int(parts[1]) - 1
                    val = float(parts[2])
                    if val > 0 and i != j:
                        contacts.append([i, j])
    except (ValueError, IndexError):
        return None

    return contacts if contacts else None


def parse_sword2_output(
    protein_id: str, sequence: str, result_dir: str
) -> Optional[dict]:
    """Parse SWORD2 JSON output into training format."""
    # Find the summary JSON (SWORD2 names it SWORD2_summary.json)
    json_files = list(Path(result_dir).rglob("SWORD2_summary.json"))
    if not json_files:
        return None

    summary_file = json_files[0]
    with open(summary_file) as f:
        summary = json.load(f)

    partitionings = []

    for key, value in summary.items():
        if key == "Ambiguity index":
            continue
        if not isinstance(value, dict) or "Domains" not in value:
            continue

        # Parse quality stars
        quality_str = value.get("Quality", "n/a")
        quality = quality_str.count("*") if quality_str != "n/a" else 0

        domains = []
        for domain_key, domain_data in value.get("Domains", {}).items():
            segments = []
            for pu_range in domain_data.get("PUs", {}):
                # Parse range like "1-100" or "1-50;151-200"
                for seg_str in pu_range.split(";"):
                    seg_str = seg_str.strip()
                    match = re.match(r"(\d+)-(\d+)", seg_str)
                    if match:
                        start = int(match.group(1)) - 1  # convert to 0-indexed
                        end = int(match.group(2)) - 1
                        segments.append([start, end])

            if segments:
                domains.append({"segments": segments})

        if domains:
            partitionings.append(
                {
                    "num_domains": value.get("Nb. domains", len(domains)),
                    "quality": quality,
                    "domains": domains,
                }
            )

    if not partitionings:
        return None

    result = {
        "id": protein_id,
        "sequence": sequence,
        "partitionings": partitionings,
    }

    # Parse contact map from SWORD2's contact matrix output
    contacts = parse_contact_matrix(result_dir)
    if contacts:
        result["contact_map"] = contacts

    return result


def process_batch(
    proteins: list[dict],
    sword2_binary: str,
    base_dir: str,
    output_dir: str,
    num_workers: int = 4,
    timeout: int = 120,
    pdb_dir: Optional[str] = None,
) -> list[dict]:
    """Process a batch of proteins in parallel with resume support."""
    # Resume: skip already-completed proteins
    completed = get_completed_ids(output_dir)
    remaining = [p for p in proteins if p["id"] not in completed]

    # Clean up incomplete runs (interrupted mid-SWORD2)
    remaining_ids = {p["id"] for p in remaining}
    n_incomplete = cleanup_incomplete_runs(output_dir, remaining_ids)

    logger.info(
        f"Resuming: {len(completed)} done, {n_incomplete} incomplete (will retry), "
        f"{len(remaining)} remaining"
    )

    if not remaining:
        logger.info("All proteins already processed")
        return []

    results = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {}
        for prot in remaining:
            future = executor.submit(
                run_sword2_on_sequence,
                prot["id"],
                prot["sequence"],
                sword2_binary,
                base_dir,
                output_dir,
                timeout,
                pdb_dir,
            )
            futures[future] = prot["id"]

        for future in tqdm(
            as_completed(futures), total=len(futures), desc="Running SWORD2"
        ):
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
            except Exception as e:
                logger.warning(f"Worker error for {futures[future]}: {e}")

    return results


def split_and_save(
    results_dir: str,
    output_dir: str,
    train_ratio: float = 0.9,
    val_ratio: float = 0.05,
) -> None:
    """Load per-protein results from disk, split into train/val/test, and save."""
    # Load all results from per-protein JSONs (sorted for deterministic ordering)
    results = []
    for json_file in sorted(Path(results_dir).glob("*.json")):
        with open(json_file) as f:
            results.append(json.load(f))

    if not results:
        logger.warning("No results to split")
        return

    logger.info(f"Loaded {len(results)} results for splitting")

    np.random.seed(42)
    indices = np.random.permutation(len(results))

    n_train = int(len(results) * train_ratio)
    n_val = int(len(results) * val_ratio)

    splits = {
        "train": indices[:n_train],
        "val": indices[n_train : n_train + n_val],
        "test": indices[n_train + n_val :],
    }

    for split_name, split_indices in splits.items():
        split_dir = Path(output_dir) / split_name
        split_dir.mkdir(parents=True, exist_ok=True)

        split_data = [results[i] for i in split_indices]

        # Save manifest
        manifest_path = Path(output_dir) / f"{split_name}.json"
        with open(manifest_path, "w") as f:
            json.dump(split_data, f)

        # Also save individual files for flexible loading
        for item in split_data:
            item_path = split_dir / f"{item['id']}.json"
            with open(item_path, "w") as f:
                json.dump(item, f)

        logger.info(f"{split_name}: {len(split_data)} proteins")


def main():
    """Main data generation pipeline."""
    import argparse

    parser = argparse.ArgumentParser(description="Generate SWORD3 training data")
    parser.add_argument(
        "--swissprot-fasta",
        default="data/swissprot.fasta",
        help="Path to SwissProt FASTA (downloaded if missing)",
    )
    parser.add_argument(
        "--sword2-binary",
        default="./target/release/sword2",
        help="Path to SWORD2 binary",
    )
    parser.add_argument(
        "--base-dir",
        default=".",
        help="SWORD2 base directory (for finding bin/)",
    )
    parser.add_argument(
        "--output-dir",
        default="data/sword2_results",
        help="Output directory for SWORD2 results",
    )
    parser.add_argument(
        "--processed-dir",
        default="data/processed",
        help="Output directory for processed training data",
    )
    parser.add_argument(
        "--max-proteins",
        type=int,
        default=None,
        help="Maximum number of proteins to process (None = all)",
    )
    parser.add_argument(
        "--min-seq-len",
        type=int,
        default=15,
        help="Minimum sequence length",
    )
    parser.add_argument(
        "--max-seq-len",
        type=int,
        default=2048,
        help="Maximum sequence length",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=40,
        help="Number of parallel workers",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Timeout per protein (seconds)",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download SwissProt if not present",
    )
    parser.add_argument(
        "--pdb-dir",
        default=None,
        help="Directory containing local AlphaFold PDB files (AF-{id}-F1-model_v4.pdb). Uses local files instead of downloading.",
    )
    parser.add_argument(
        "--skip-split",
        action="store_true",
        help="Skip the final train/val/test split (just run SWORD2 processing)",
    )
    parser.add_argument(
        "--rerun-failures",
        action="store_true",
        help="Re-process proteins that previously failed",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Step 1: Download SwissProt
    if args.download or not Path(args.swissprot_fasta).exists():
        os.makedirs(os.path.dirname(args.swissprot_fasta) or ".", exist_ok=True)
        download_swissprot(args.swissprot_fasta)

    # Step 2: Parse sequences
    logger.info("Parsing SwissProt FASTA...")
    proteins = parse_fasta(args.swissprot_fasta)
    logger.info(f"Found {len(proteins)} proteins")

    # Filter by length
    proteins = [
        p
        for p in proteins
        if args.min_seq_len <= len(p["sequence"]) <= args.max_seq_len
    ]
    logger.info(f"After length filter: {len(proteins)} proteins")

    if args.max_proteins:
        proteins = proteins[: args.max_proteins]
        logger.info(f"Using first {len(proteins)} proteins")

    # Step 3: Set up output directories
    results_dir = os.path.join(args.output_dir, RESULTS_SUBDIR)
    os.makedirs(results_dir, exist_ok=True)

    # Clean up stale temp files from interrupted atomic writes
    for tmp in Path(results_dir).glob(".*.json.tmp"):
        tmp.unlink()

    # Clear failures log if re-running failures
    if args.rerun_failures:
        failures_path = os.path.join(args.output_dir, FAILURES_LOG)
        if os.path.exists(failures_path):
            logger.info("Re-running failures: clearing failures log")
            os.remove(failures_path)

    # Step 4: Run SWORD2 (with resume support)
    results = process_batch(
        proteins,
        sword2_binary=args.sword2_binary,
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        num_workers=args.num_workers,
        timeout=args.timeout,
        pdb_dir=args.pdb_dir,
    )

    total_on_disk = len(list(Path(results_dir).glob("*.json")))
    logger.info(
        f"This run: {len(results)} new results. "
        f"Total on disk: {total_on_disk}/{len(proteins)} proteins"
    )

    # Step 5: Split and save
    if args.skip_split:
        logger.info("Skipping split (--skip-split)")
    else:
        os.makedirs(args.processed_dir, exist_ok=True)
        split_and_save(results_dir, args.processed_dir)

    logger.info("Data generation complete!")


if __name__ == "__main__":
    main()
