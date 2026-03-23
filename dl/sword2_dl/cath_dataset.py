"""
CATH domain annotation parser and dataset generator.

Downloads CATH domain annotations and converts them to the same training
format used by SWORD2-DL, enabling training on gold-standard domain labels.

CATH domain list format (cath-domain-list-v4_4_0.txt):
    Each line: DOMAIN_ID CLASS ARCH TOPOL HOMOL SEG_COUNT SEGMENTS
    e.g.: 1oaiA00 2 60 40 10 1 2-126
          1cukA01 3 40 50 720 2 1-100 150-220

Segments are space-separated start-end ranges (1-indexed, inclusive).
"""

import json
import logging
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def download_cath_domain_list(output_path: str, version: str = "4_4_0") -> str:
    """Download CATH domain list file."""
    import urllib.request

    url = f"http://download.cathdb.info/cath/releases/all-releases/v{version}/cath-classification-data/cath-domain-list-v{version}.txt"

    if Path(output_path).exists():
        logger.info(f"CATH domain list already exists at {output_path}")
        return output_path

    logger.info(f"Downloading CATH domain list from {url}...")
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    urllib.request.urlretrieve(url, output_path)
    logger.info(f"Saved to {output_path}")
    return output_path


def download_cath_domain_boundaries(output_path: str, version: str = "4_4_0") -> str:
    """Download CATH domain boundaries file with segment info."""
    import urllib.request

    url = f"http://download.cathdb.info/cath/releases/all-releases/v{version}/cath-classification-data/cath-domain-boundaries-v{version}.txt"

    if Path(output_path).exists():
        logger.info(f"CATH boundaries file already exists at {output_path}")
        return output_path

    logger.info(f"Downloading CATH domain boundaries from {url}...")
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    urllib.request.urlretrieve(url, output_path)
    logger.info(f"Saved to {output_path}")
    return output_path


def parse_cath_domain_list(filepath: str) -> dict[str, list[dict]]:
    """Parse CATH domain list into per-chain domain annotations.

    Returns:
        Dict mapping chain_id (e.g. "1oaiA") to list of domain dicts,
        each with 'domain_id', 'cath_class', and 'segments'.
    """
    chains = defaultdict(list)

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 7:
                continue

            domain_id = parts[0]  # e.g. "1oaiA00"
            chain_id = domain_id[:5]  # e.g. "1oaiA"
            pdb_id = domain_id[:4]  # e.g. "1oai"
            chain_letter = domain_id[4]

            cath_class = tuple(int(x) for x in parts[1:5])

            n_segments = int(parts[5])
            segments = []

            seg_parts = parts[6:]
            for seg_str in seg_parts[:n_segments]:
                match = re.match(r"(-?\d+)-(-?\d+)", seg_str)
                if match:
                    start = int(match.group(1)) - 1  # convert to 0-indexed
                    end = int(match.group(2)) - 1
                    if start >= 0 and end >= start:
                        segments.append([start, end])

            if segments:
                chains[chain_id].append({
                    "domain_id": domain_id,
                    "cath_class": cath_class,
                    "segments": segments,
                })

    logger.info(f"Parsed {sum(len(v) for v in chains.values())} domains "
                f"across {len(chains)} chains from CATH")
    return dict(chains)


def parse_cath_domain_boundaries(filepath: str) -> dict[str, list[dict]]:
    """Parse CATH domain boundaries file (more detailed segment info).

    The boundaries file has format:
        DOMAIN_ID  D{domain_num}  FRAGMENT  start  stop  length
        with multiple segments per domain.

    Returns:
        Dict mapping chain_id to list of domain dicts.
    """
    chains = defaultdict(lambda: defaultdict(list))

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 6:
                continue

            domain_id = parts[0]
            chain_id = domain_id[:5]

            try:
                start = int(parts[3]) - 1  # 0-indexed
                end = int(parts[4]) - 1
                if start >= 0 and end >= start:
                    chains[chain_id][domain_id].append([start, end])
            except (ValueError, IndexError):
                continue

    # Convert to standard format
    result = {}
    for chain_id, domains in chains.items():
        result[chain_id] = [
            {"domain_id": did, "segments": segs}
            for did, segs in domains.items()
            if segs
        ]

    logger.info(f"Parsed {sum(len(v) for v in result.values())} domains "
                f"across {len(result)} chains from CATH boundaries")
    return result


def cath_to_training_format(
    chain_domains: list[dict],
    sequence: str,
    protein_id: str,
) -> Optional[dict]:
    """Convert CATH domain annotations for a chain to SWORD2-DL training format.

    CATH provides a single domain assignment (one partitioning).
    We create a single partitioning entry.

    Args:
        chain_domains: List of domain dicts with 'segments' field.
        sequence: Protein sequence.
        protein_id: Protein identifier.

    Returns:
        Training sample dict or None if invalid.
    """
    seq_len = len(sequence)

    # Validate segments are within sequence bounds
    domains = []
    for domain in chain_domains:
        valid_segments = []
        for start, end in domain["segments"]:
            start = max(0, start)
            end = min(seq_len - 1, end)
            if end >= start:
                valid_segments.append([start, end])
        if valid_segments:
            domains.append({"segments": valid_segments})

    if not domains:
        return None

    partitioning = {
        "num_domains": len(domains),
        "quality": 5,  # CATH annotations are gold-standard
        "domains": domains,
    }

    return {
        "id": protein_id,
        "sequence": sequence,
        "partitionings": [partitioning],
    }


def generate_cath_dataset(
    cath_file: str,
    sequence_source: str,
    output_dir: str,
    train_ratio: float = 0.9,
    val_ratio: float = 0.05,
    min_seq_len: int = 30,
    max_seq_len: int = 1500,
) -> None:
    """Generate SWORD2-DL training data from CATH annotations.

    Args:
        cath_file: Path to CATH domain list file.
        sequence_source: Path to FASTA file or directory with sequence files.
        output_dir: Where to save processed training data.
        train_ratio: Fraction for training set.
        val_ratio: Fraction for validation set.
        min_seq_len: Minimum sequence length.
        max_seq_len: Maximum sequence length.
    """
    import numpy as np

    # Parse CATH
    cath_chains = parse_cath_domain_list(cath_file)

    # Load sequences
    sequences = _load_sequences(sequence_source)
    logger.info(f"Loaded {len(sequences)} sequences")

    # Match CATH annotations with sequences
    results = []
    matched = 0
    for chain_id, domains in cath_chains.items():
        pdb_id = chain_id[:4].lower()
        chain_letter = chain_id[4]
        lookup_id = f"{pdb_id}_{chain_letter}"

        seq = sequences.get(lookup_id) or sequences.get(chain_id) or sequences.get(pdb_id)
        if seq is None:
            continue

        if not (min_seq_len <= len(seq) <= max_seq_len):
            continue

        sample = cath_to_training_format(domains, seq, chain_id)
        if sample is not None:
            results.append(sample)
            matched += 1

    logger.info(f"Matched {matched}/{len(cath_chains)} CATH chains with sequences")

    if not results:
        logger.error("No results to save!")
        return

    # Split and save
    np.random.seed(42)
    indices = np.random.permutation(len(results))
    n_train = int(len(results) * train_ratio)
    n_val = int(len(results) * val_ratio)

    splits = {
        "train": indices[:n_train],
        "val": indices[n_train:n_train + n_val],
        "test": indices[n_train + n_val:],
    }

    os.makedirs(output_dir, exist_ok=True)
    for split_name, split_indices in splits.items():
        split_data = [results[i] for i in split_indices]
        manifest_path = Path(output_dir) / f"{split_name}.json"
        with open(manifest_path, "w") as f:
            json.dump(split_data, f)
        logger.info(f"{split_name}: {len(split_data)} proteins")

    logger.info("CATH dataset generation complete!")


def _load_sequences(source: str) -> dict[str, str]:
    """Load sequences from FASTA file or directory."""
    from .generate_data import parse_fasta

    path = Path(source)
    sequences = {}

    if path.is_file() and path.suffix in (".fasta", ".fa", ".faa"):
        for entry in parse_fasta(str(path)):
            sequences[entry["id"]] = entry["sequence"]
    elif path.is_dir():
        for fasta_file in path.glob("*.fasta"):
            for entry in parse_fasta(str(fasta_file)):
                sequences[entry["id"]] = entry["sequence"]
    else:
        logger.warning(f"Cannot load sequences from {source}")

    return sequences


def main():
    """CLI entry point for CATH dataset generation."""
    import argparse

    parser = argparse.ArgumentParser(description="Generate CATH training data for SWORD2-DL")
    parser.add_argument("--cath-file", required=True, help="CATH domain list file")
    parser.add_argument("--sequences", required=True, help="FASTA file or directory with sequences")
    parser.add_argument("--output-dir", default="data/processed_cath", help="Output directory")
    parser.add_argument("--download", action="store_true", help="Download CATH file if not present")
    parser.add_argument("--min-seq-len", type=int, default=30)
    parser.add_argument("--max-seq-len", type=int, default=1500)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    if args.download and not Path(args.cath_file).exists():
        download_cath_domain_list(args.cath_file)

    generate_cath_dataset(
        cath_file=args.cath_file,
        sequence_source=args.sequences,
        output_dir=args.output_dir,
        min_seq_len=args.min_seq_len,
        max_seq_len=args.max_seq_len,
    )


if __name__ == "__main__":
    main()
