"""
Dataset and DataLoader for SWORD2-DL training.

Reads pre-processed SWORD2 results and converts them to training tensors.
Each sample contains:
- Protein sequence (tokenized by ESM-2)
- Multiple ground truth partitionings (as domain segment lists)
"""

import json
import logging
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

logger = logging.getLogger(__name__)


class Sword2Dataset(Dataset):
    """Dataset of protein sequences with SWORD2 domain partitionings.

    Each sample is stored as a JSON file with structure:
    {
        "id": "P12345",
        "sequence": "MVLSPADKTN...",
        "partitionings": [
            {
                "num_domains": 2,
                "quality": 5,
                "domains": [
                    {"segments": [[0, 99], [200, 250]]},  // discontinuous domain
                    {"segments": [[100, 199]]}             // continuous domain
                ]
            },
            ...  // alternative partitionings
        ]
    }
    """

    def __init__(
        self,
        data_dir: str | Path,
        esm_alphabet=None,
        max_seq_len: int = 1024,
        min_seq_len: int = 30,
        split: str = "train",
    ):
        self.data_dir = Path(data_dir)
        self.max_seq_len = max_seq_len
        self.min_seq_len = min_seq_len
        self.esm_alphabet = esm_alphabet

        # Load manifest
        manifest_file = self.data_dir / f"{split}.json"
        if manifest_file.exists():
            with open(manifest_file) as f:
                self.samples = json.load(f)
        else:
            # Scan directory for individual files
            self.samples = self._scan_directory(split)

        logger.info(f"Loaded {len(self.samples)} samples for {split} split")

    def _scan_directory(self, split: str) -> list[dict]:
        """Scan directory for sample files."""
        split_dir = self.data_dir / split
        if not split_dir.exists():
            logger.warning(f"Split directory {split_dir} not found")
            return []

        samples = []
        for f in sorted(split_dir.glob("*.json")):
            try:
                with open(f) as fh:
                    data = json.load(fh)
                seq_len = len(data["sequence"])
                if self.min_seq_len <= seq_len <= self.max_seq_len:
                    data["_path"] = str(f)
                    samples.append(data)
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"Skipping {f}: {e}")

        return samples

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> dict:
        sample = self.samples[idx]

        sequence = sample["sequence"]
        seq_len = len(sequence)

        # Tokenize with ESM-2 alphabet
        if self.esm_alphabet is not None:
            batch_converter = self.esm_alphabet.get_batch_converter()
            _, _, tokens = batch_converter([(sample["id"], sequence)])
            tokens = tokens.squeeze(0)  # (L+2,) includes BOS/EOS
        else:
            # Dummy tokens for testing without ESM
            tokens = torch.zeros(seq_len + 2, dtype=torch.long)

        # Parse partitionings
        partitionings = []
        for part in sample.get("partitionings", []):
            domains = []
            for domain in part["domains"]:
                segments = [tuple(seg) for seg in domain["segments"]]
                domains.append(segments)
            partitionings.append(domains)

        return {
            "id": sample["id"],
            "tokens": tokens,
            "seq_len": seq_len,
            "partitionings": partitionings,
        }


def collate_fn(batch: list[dict]) -> dict:
    """Custom collation: pad tokens to max length in batch.

    Returns:
        Dictionary with:
            - tokens: (B, max_len+2) padded token tensor
            - mask: (B, max_len) boolean mask
            - targets: list of target dicts for loss computation
    """
    max_len = max(item["seq_len"] for item in batch)
    B = len(batch)

    # Pad tokens
    tokens = torch.zeros(B, max_len + 2, dtype=torch.long)  # +2 for BOS/EOS
    mask = torch.zeros(B, max_len, dtype=torch.bool)

    targets = []

    for i, item in enumerate(batch):
        L = item["seq_len"]
        tok = item["tokens"]
        # Copy tokens: BOS + sequence + EOS + padding
        tokens[i, : len(tok)] = tok
        # Set padding token to 1 (ESM-2 padding token index)
        tokens[i, len(tok) :] = 1
        mask[i, :L] = True

        targets.append(
            {
                "id": item["id"],
                "seq_len": L,
                "partitionings": item["partitionings"],
            }
        )

    return {
        "tokens": tokens,
        "mask": mask,
        "targets": targets,
    }


def create_dataloaders(
    data_dir: str | Path,
    esm_alphabet=None,
    max_seq_len: int = 1024,
    min_seq_len: int = 30,
    batch_size: int = 4,
    num_workers: int = 4,
    pin_memory: bool = True,
) -> dict[str, DataLoader]:
    """Create train/val/test dataloaders."""
    loaders = {}

    for split in ["train", "val", "test"]:
        dataset = Sword2Dataset(
            data_dir=data_dir,
            esm_alphabet=esm_alphabet,
            max_seq_len=max_seq_len,
            min_seq_len=min_seq_len,
            split=split,
        )

        if len(dataset) == 0:
            logger.warning(f"Empty dataset for {split} split")
            continue

        loaders[split] = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=(split == "train"),
            num_workers=num_workers,
            pin_memory=pin_memory,
            collate_fn=collate_fn,
            drop_last=(split == "train"),
        )

    return loaders
