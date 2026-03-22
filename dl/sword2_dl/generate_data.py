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
import tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

import numpy as np
from tqdm import tqdm

logger = logging.getLogger(__name__)


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


def run_sword2_on_sequence(
    protein_id: str,
    sequence: str,
    sword2_binary: str,
    base_dir: str,
    output_dir: str,
    timeout: int = 120,
) -> Optional[dict]:
    """Run SWORD2 on a single protein sequence.

    Creates a temporary PDB-like file from AlphaFold or runs with UniProt ID.

    Returns parsed result dict or None on failure.
    """
    result_dir = os.path.join(output_dir, protein_id)

    try:
        # Run SWORD2 with UniProt ID (fetches from AlphaFold)
        cmd = [
            sword2_binary,
            "-u", protein_id,
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
            logger.debug(f"SWORD2 failed for {protein_id}: {proc.stderr[:200]}")
            return None

        # Find and parse the JSON output
        return parse_sword2_output(protein_id, sequence, result_dir)

    except subprocess.TimeoutExpired:
        logger.debug(f"SWORD2 timed out for {protein_id}")
        return None
    except Exception as e:
        logger.debug(f"Error processing {protein_id}: {e}")
        return None
    finally:
        # Clean up intermediate files to save disk space
        intermediate = os.path.join(result_dir, "intermediate")
        if os.path.exists(intermediate):
            shutil.rmtree(intermediate, ignore_errors=True)


def parse_sword2_output(
    protein_id: str, sequence: str, result_dir: str
) -> Optional[dict]:
    """Parse SWORD2 JSON output into training format."""
    # Find the summary JSON
    json_files = list(Path(result_dir).rglob("summary.json"))
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

    return {
        "id": protein_id,
        "sequence": sequence,
        "partitionings": partitionings,
    }


def process_batch(
    proteins: list[dict],
    sword2_binary: str,
    base_dir: str,
    output_dir: str,
    num_workers: int = 4,
    timeout: int = 120,
) -> list[dict]:
    """Process a batch of proteins in parallel."""
    results = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {}
        for prot in proteins:
            future = executor.submit(
                run_sword2_on_sequence,
                prot["id"],
                prot["sequence"],
                sword2_binary,
                base_dir,
                output_dir,
                timeout,
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
    results: list[dict],
    output_dir: str,
    train_ratio: float = 0.9,
    val_ratio: float = 0.05,
) -> None:
    """Split results into train/val/test and save."""
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

    parser = argparse.ArgumentParser(description="Generate SWORD2-DL training data")
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
        default=30,
        help="Minimum sequence length",
    )
    parser.add_argument(
        "--max-seq-len",
        type=int,
        default=1500,
        help="Maximum sequence length",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=8,
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

    # Step 3: Run SWORD2
    os.makedirs(args.output_dir, exist_ok=True)
    results = process_batch(
        proteins,
        sword2_binary=args.sword2_binary,
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        num_workers=args.num_workers,
        timeout=args.timeout,
    )
    logger.info(f"Successfully processed {len(results)}/{len(proteins)} proteins")

    # Step 4: Split and save
    os.makedirs(args.processed_dir, exist_ok=True)
    split_and_save(results, args.processed_dir)

    logger.info("Data generation complete!")


if __name__ == "__main__":
    main()
