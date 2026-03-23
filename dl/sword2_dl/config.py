"""Configuration dataclasses for SWORD2-DL."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


@dataclass
class EmbeddingSource:
    """A single PLM embedding source.

    Each source is a directory of safetensors files (one per protein).
    Each safetensors file contains a single tensor keyed by the protein ID,
    with shape (L, embed_dim).
    """

    name: str = "esm2_650M"
    embed_dim: int = 1280
    path: str = "data/embeddings/esm2_650M"  # directory of .safetensors files


@dataclass
class ModelConfig:
    # Embedding inputs (list of PLM sources to concatenate)
    embedding_sources: list[dict] = field(
        default_factory=lambda: [
            {"name": "esm2_650M", "embed_dim": 1280, "path": "data/embeddings/esm2_650M"},
        ]
    )

    # Pair module
    pair_dim: int = 128
    pair_num_blocks: int = 8
    pair_kernel_size: int = 3
    pair_dilation_cycle: list[int] = field(default_factory=lambda: [1, 2, 4, 8])
    pair_dropout: float = 0.1

    # Single (per-residue) refinement
    single_dim: int = 256
    single_num_layers: int = 4
    single_num_heads: int = 8
    single_dropout: float = 0.1

    # Multi-partitioning heads
    num_partitioning_slots: int = 10
    head_hidden_dim: int = 64

    # Auxiliary predictions
    predict_num_domains: bool = True
    predict_boundaries: bool = True
    predict_contact_map: bool = True
    max_num_domains: int = 20

    @property
    def total_embed_dim(self) -> int:
        """Sum of all embedding source dimensions."""
        return sum(s["embed_dim"] for s in self.embedding_sources)

    def get_embedding_sources(self) -> list[EmbeddingSource]:
        """Parse embedding source dicts into EmbeddingSource objects."""
        return [EmbeddingSource(**s) for s in self.embedding_sources]


@dataclass
class DataConfig:
    # SWORD2 data generation
    sword2_binary: str = "./target/release/sword2"
    swissprot_fasta: str = "data/swissprot.fasta"
    sword2_output_dir: str = "data/sword2_results"
    processed_data_dir: str = "data/processed"

    # Dataset
    max_seq_len: int = 1024
    min_seq_len: int = 30
    min_num_partitionings: int = 1
    train_split: float = 0.9
    val_split: float = 0.05
    # test_split = 1 - train_split - val_split

    # Dataloader
    batch_size: int = 8
    num_workers: int = 4
    pin_memory: bool = True


@dataclass
class TrainConfig:
    # Optimizer
    learning_rate: float = 3e-4
    weight_decay: float = 0.01
    warmup_steps: int = 1000
    max_steps: int = 200_000
    gradient_clip: float = 1.0

    # Loss weights
    co_membership_weight: float = 1.0
    confidence_weight: float = 0.1
    num_domains_weight: float = 0.1
    boundary_weight: float = 0.5
    contact_map_weight: float = 0.5
    dice_weight: float = 0.5  # weight of Dice loss vs BCE in co-membership loss
    boundary_smooth_width: int = 2  # label smoothing width for boundary targets

    # EMA
    use_ema: bool = True
    ema_decay: float = 0.999

    # Training
    fp16: bool = True
    gradient_accumulation: int = 4
    eval_every: int = 1000
    save_every: int = 5000
    log_every: int = 50

    # Checkpointing
    output_dir: str = "checkpoints"
    resume_from: Optional[str] = None

    # W&B
    wandb_project: str = "sword2-dl"
    wandb_run_name: Optional[str] = None


@dataclass
class Config:
    model: ModelConfig = field(default_factory=ModelConfig)
    data: DataConfig = field(default_factory=DataConfig)
    train: TrainConfig = field(default_factory=TrainConfig)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Config":
        with open(path) as f:
            raw = yaml.safe_load(f)
        config = cls()
        for section_name, section_data in (raw or {}).items():
            if hasattr(config, section_name) and isinstance(section_data, dict):
                section = getattr(config, section_name)
                for k, v in section_data.items():
                    if hasattr(section, k):
                        setattr(section, k, v)
        return config

    def to_yaml(self, path: str | Path) -> None:
        from dataclasses import asdict

        with open(path, "w") as f:
            yaml.dump(asdict(self), f, default_flow_style=False, sort_keys=False)
