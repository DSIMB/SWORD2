"""
Dataset and DataLoader for SWORD2-DL training.

Loads pre-computed PLM embeddings from safetensors files and pairs them with
SWORD2 ground truth partitionings. Supports multiple embedding sources that
get concatenated per residue.

Expected directory layout:
    data/
    ├── processed/
    │   ├── train.json          # manifest with sequences + partitionings
    │   ├── val.json
    │   └── test.json
    └── embeddings/
        ├── esm2_650M/                  # one directory per PLM
        │   ├── P12345.safetensors      # key="P12345", value=(L, 1280)
        │   ├── P67890.safetensors
        │   └── ...
        ├── ankh2_large/
        │   ├── P12345.safetensors      # key="P12345", value=(L, 1536)
        │   └── ...
        └── esmc_600M/
            ├── P12345.safetensors      # key="P12345", value=(L, D)
            └── ...
"""

import json
import logging
from pathlib import Path

import torch
from safetensors.torch import load_file as load_safetensors
from torch.utils.data import Dataset, DataLoader

from .config import EmbeddingSource

logger = logging.getLogger(__name__)


class Sword2Dataset(Dataset):
    """Dataset of pre-computed PLM embeddings with SWORD2 domain partitionings.

    Each sample in the manifest JSON has structure:
    {
        "id": "P12345",
        "sequence": "MVLSPADKTN...",
        "partitionings": [
            {
                "num_domains": 2,
                "quality": 5,
                "domains": [
                    {"segments": [[0, 99], [200, 250]]},
                    {"segments": [[100, 199]]}
                ]
            },
            ...
        ],
        "contact_map": [[0,5], [0,6], ...]  // optional: list of contacting residue pairs
    }
    """

    def __init__(
        self,
        data_dir: str | Path,
        embedding_sources: list[EmbeddingSource],
        max_seq_len: int = 1024,
        min_seq_len: int = 30,
        split: str = "train",
    ):
        self.data_dir = Path(data_dir)
        self.embedding_sources = embedding_sources
        self.max_seq_len = max_seq_len
        self.min_seq_len = min_seq_len

        # Load manifest
        manifest_file = self.data_dir / f"{split}.json"
        if manifest_file.exists():
            with open(manifest_file) as f:
                all_samples = json.load(f)
        else:
            all_samples = self._scan_directory(split)

        # Filter by sequence length and embedding availability
        self.samples = []
        for s in all_samples:
            seq_len = len(s["sequence"])
            if not (self.min_seq_len <= seq_len <= self.max_seq_len):
                continue
            if not self._embeddings_exist(s["id"]):
                continue
            self.samples.append(s)

        logger.info(
            f"Loaded {len(self.samples)}/{len(all_samples)} samples for {split} "
            f"(filtered by length [{min_seq_len}, {max_seq_len}] and embedding availability)"
        )

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
                samples.append(data)
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"Skipping {f}: {e}")
        return samples

    def _embeddings_exist(self, protein_id: str) -> bool:
        """Check if safetensors embedding files exist for all sources."""
        for source in self.embedding_sources:
            path = Path(source.path) / f"{protein_id}.safetensors"
            if not path.exists():
                return False
        return True

    def _load_embedding(self, protein_id: str, source: EmbeddingSource) -> torch.Tensor:
        """Load a single PLM embedding from a safetensors file.

        Each safetensors file contains one tensor keyed by the protein ID.
        """
        path = Path(source.path) / f"{protein_id}.safetensors"
        tensors = load_safetensors(path)

        if protein_id in tensors:
            emb = tensors[protein_id]
        else:
            # Fallback: use the first (and presumably only) key
            emb = next(iter(tensors.values()))

        emb = emb.float()
        if emb.dim() == 3:
            emb = emb.squeeze(0)  # remove batch dim if present
        return emb

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> dict:
        sample = self.samples[idx]
        protein_id = sample["id"]
        seq_len = len(sample["sequence"])

        # Load and concatenate embeddings from all PLM sources
        embeddings = []
        for source in self.embedding_sources:
            emb = self._load_embedding(protein_id, source)
            # Truncate to sequence length if needed
            if emb.shape[0] > seq_len:
                emb = emb[:seq_len]
            elif emb.shape[0] < seq_len:
                pad = torch.zeros(seq_len - emb.shape[0], emb.shape[1])
                emb = torch.cat([emb, pad], dim=0)
            embeddings.append(emb)

        # Concatenate along feature dimension: (L, D1+D2+...+Dn)
        embeddings = torch.cat(embeddings, dim=-1)

        # Parse partitionings
        partitionings = []
        for part in sample.get("partitionings", []):
            domains = []
            for domain in part["domains"]:
                segments = [tuple(seg) for seg in domain["segments"]]
                domains.append(segments)
            partitionings.append(domains)

        # Parse contact map (list of [i, j] pairs → sparse representation)
        contact_pairs = sample.get("contact_map", None)

        result = {
            "id": protein_id,
            "embeddings": embeddings,  # (L, D_total)
            "seq_len": seq_len,
            "partitionings": partitionings,
        }

        if contact_pairs is not None:
            result["contact_pairs"] = contact_pairs

        return result


def collate_fn(batch: list[dict]) -> dict:
    """Custom collation: pad embeddings to max length in batch.

    Returns:
        Dictionary with:
            - embeddings: (B, max_len, D) padded embedding tensor
            - mask: (B, max_len) boolean mask
            - targets: list of target dicts for loss computation
    """
    max_len = max(item["seq_len"] for item in batch)
    B = len(batch)
    embed_dim = batch[0]["embeddings"].shape[-1]

    embeddings = torch.zeros(B, max_len, embed_dim)
    mask = torch.zeros(B, max_len, dtype=torch.bool)

    targets = []

    for i, item in enumerate(batch):
        L = item["seq_len"]
        embeddings[i, :L] = item["embeddings"]
        mask[i, :L] = True

        target = {
            "id": item["id"],
            "seq_len": L,
            "partitionings": item["partitionings"],
        }

        # Build contact map matrix if contact pairs are available
        if "contact_pairs" in item:
            contact_map = torch.zeros(max_len, max_len)
            for ci, cj in item["contact_pairs"]:
                if ci < L and cj < L:
                    contact_map[ci, cj] = 1.0
                    contact_map[cj, ci] = 1.0
            target["contact_map"] = contact_map

        targets.append(target)

    return {
        "embeddings": embeddings,
        "mask": mask,
        "targets": targets,
    }


def create_dataloaders(
    data_dir: str | Path,
    embedding_sources: list[EmbeddingSource],
    max_seq_len: int = 1024,
    min_seq_len: int = 30,
    batch_size: int = 8,
    num_workers: int = 4,
    pin_memory: bool = True,
) -> dict[str, DataLoader]:
    """Create train/val/test dataloaders."""
    loaders = {}

    for split in ["train", "val", "test"]:
        dataset = Sword2Dataset(
            data_dir=data_dir,
            embedding_sources=embedding_sources,
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
