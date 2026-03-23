"""
Inference script for DomainPartitionNet.

Takes pre-computed embeddings from safetensors files and predicts alternative
domain partitionings, plus a contact map for visualization.
"""

import json
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import torch
from safetensors.torch import load_file as load_safetensors

from .config import Config, ModelConfig, EmbeddingSource
from .model import DomainPartitionNet
from .postprocess import predict_partitionings

logger = logging.getLogger(__name__)


def load_model(
    checkpoint_path: str,
    config: ModelConfig | None = None,
    device: torch.device | None = None,
) -> DomainPartitionNet:
    """Load trained model from checkpoint."""
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if config is None:
        config = ModelConfig()

    model = DomainPartitionNet(config)

    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint["model"])
    model = model.to(device)
    model.eval()

    return model


def load_embeddings(
    protein_id: str,
    embedding_sources: list[EmbeddingSource],
) -> torch.Tensor:
    """Load and concatenate pre-computed embeddings from safetensors files."""
    embeddings = []
    for source in embedding_sources:
        filename = source.filename_template.format(id=protein_id)
        path = Path(source.path) / filename
        tensors = load_safetensors(str(path))

        if protein_id in tensors:
            emb = tensors[protein_id]
        else:
            emb = next(iter(tensors.values()))

        emb = emb.float()
        if emb.dim() == 3:
            emb = emb.squeeze(0)
        embeddings.append(emb)

    return torch.cat(embeddings, dim=-1)  # (L, D_total)


def predict_from_embeddings(
    model: DomainPartitionNet,
    embeddings: torch.Tensor,
    device: torch.device | None = None,
    min_confidence: float = 0.1,
    min_domain_size: int = 20,
) -> tuple[list[dict], Optional[np.ndarray]]:
    """Predict domain partitionings and contact map from pre-computed embeddings.

    Args:
        model: Trained DomainPartitionNet.
        embeddings: (L, D) concatenated PLM embeddings.
        device: Target device.
        min_confidence: Minimum confidence to report.
        min_domain_size: Minimum domain size in residues.

    Returns:
        Tuple of (partitionings list, contact_map array or None).
    """
    if device is None:
        device = next(model.parameters()).device

    seq_len = embeddings.shape[0]
    embeddings = embeddings.unsqueeze(0).to(device)  # (1, L, D)
    mask = torch.ones(1, seq_len, dtype=torch.bool, device=device)

    with torch.no_grad(), torch.cuda.amp.autocast(dtype=torch.float16):
        outputs = model(embeddings, mask)

    # Extract single-sample outputs
    single_outputs = {
        "co_membership": outputs["co_membership"][0],
        "confidence": outputs["confidence"][0],
    }
    if "num_domains_logits" in outputs:
        single_outputs["num_domains_logits"] = outputs["num_domains_logits"][0]
    if "boundary_logits" in outputs:
        single_outputs["boundary_logits"] = outputs["boundary_logits"][0]

    partitionings = predict_partitionings(
        single_outputs,
        seq_len=seq_len,
        min_confidence=min_confidence,
        min_domain_size=min_domain_size,
    )

    # Extract contact map
    contact_map = None
    if "contact_map_logits" in outputs:
        contact_map = (
            torch.sigmoid(outputs["contact_map_logits"][0, :seq_len, :seq_len])
            .cpu()
            .numpy()
        )

    return partitionings, contact_map


def format_output(
    protein_id: str,
    seq_len: int,
    partitionings: list[dict],
    contact_map: Optional[np.ndarray] = None,
) -> dict:
    """Format predictions in SWORD2-compatible JSON output."""
    output = {
        "id": protein_id,
        "sequence_length": seq_len,
        "method": "SWORD3",
    }

    for i, part in enumerate(partitionings):
        key = "Optimal partition" if i == 0 else f"Alternative partition {i}"
        partition_data = {
            "Partition": key,
            "Confidence": f"{part['confidence']:.3f}",
            "Nb. domains": part["num_domains"],
            "Domains": {},
        }

        for d_idx, domain_segments in enumerate(part["domains"]):
            domain_key = f"Domain {d_idx + 1}"
            segment_strs = [f"{s + 1}-{e + 1}" for s, e in domain_segments]
            partition_data["Domains"][domain_key] = {
                "segments": segment_strs,
                "residue_count": sum(e - s + 1 for s, e in domain_segments),
            }

        output[key] = partition_data

    if contact_map is not None:
        output["has_contact_map"] = True

    return output


def format_text_output(
    protein_id: str,
    seq_len: int,
    partitionings: list[dict],
    contact_map: Optional[np.ndarray] = None,
) -> str:
    """Format predictions as human-readable text (SWORD2-style)."""
    lines = [
        "SWORD3 Domain Prediction",
        f"Protein: {protein_id}",
        f"Sequence length: {seq_len}",
        f"Number of alternative partitionings: {len(partitionings)}",
        f"Contact map: {'yes' if contact_map is not None else 'no'}",
        "",
    ]

    for i, part in enumerate(partitionings):
        header = "Optimal partition" if i == 0 else f"Alternative partition {i}"
        lines.append("-" * 40)
        lines.append(header)
        lines.append(f"  Confidence: {part['confidence']:.3f}")
        lines.append(f"  Nb. domains: {part['num_domains']}")

        for d_idx, domain_segments in enumerate(part["domains"]):
            segment_strs = [f"{s + 1}-{e + 1}" for s, e in domain_segments]
            delineation = ";".join(segment_strs)
            n_residues = sum(e - s + 1 for s, e in domain_segments)
            lines.append(f"  Domain {d_idx + 1}: {delineation} ({n_residues} residues)")

        lines.append("")

    return "\n".join(lines)


def main():
    """Entry point for prediction."""
    import argparse

    parser = argparse.ArgumentParser(description="Predict protein domains with SWORD3")
    parser.add_argument(
        "--checkpoint", type=str, required=True,
        help="Path to model checkpoint",
    )
    parser.add_argument(
        "--protein-id", type=str, default=None,
        help="Protein ID (used to look up pre-computed embeddings)",
    )
    parser.add_argument(
        "--protein-ids", type=str, default=None,
        help="File with one protein ID per line",
    )
    parser.add_argument(
        "--embedding-file", type=str, default=None,
        help="Direct path to a single .safetensors embedding file",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output file path (default: stdout)",
    )
    parser.add_argument(
        "--contact-map-dir", type=str, default=None,
        help="Directory to save predicted contact maps as .npy files",
    )
    parser.add_argument(
        "--format", choices=["json", "text"], default="text",
        help="Output format",
    )
    parser.add_argument(
        "--config", type=str, default=None,
        help="Path to config YAML",
    )
    parser.add_argument(
        "--min-confidence", type=float, default=0.1,
        help="Minimum confidence threshold",
    )
    parser.add_argument(
        "--min-domain-size", type=int, default=20,
        help="Minimum domain size in residues",
    )
    parser.add_argument(
        "--device", type=str, default=None,
        help="Device (cuda/cpu)",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    if args.protein_id is None and args.protein_ids is None and args.embedding_file is None:
        parser.error("Must provide --protein-id, --protein-ids, or --embedding-file")

    # Load config
    config = Config.from_yaml(args.config) if args.config else Config()
    embedding_sources = config.model.get_embedding_sources()

    # Device
    if args.device:
        device = torch.device(args.device)
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Load model
    logger.info(f"Loading model from {args.checkpoint}")
    model = load_model(args.checkpoint, config.model, device)

    # Create contact map output dir
    if args.contact_map_dir:
        Path(args.contact_map_dir).mkdir(parents=True, exist_ok=True)

    # Collect protein IDs
    if args.embedding_file:
        protein_ids = [("direct", args.embedding_file)]
    elif args.protein_ids:
        with open(args.protein_ids) as f:
            protein_ids = [(pid.strip(), None) for pid in f if pid.strip()]
    else:
        protein_ids = [(args.protein_id, None)]

    logger.info(f"Predicting domains for {len(protein_ids)} protein(s)")

    all_results = []
    for pid, emb_file in protein_ids:
        # Load embeddings
        if emb_file:
            tensors = load_safetensors(emb_file)
            emb = next(iter(tensors.values())).float()
            if emb.dim() == 3:
                emb = emb.squeeze(0)
        else:
            emb = load_embeddings(pid, embedding_sources)

        seq_len = emb.shape[0]

        partitionings, contact_map = predict_from_embeddings(
            model, emb, device=device,
            min_confidence=args.min_confidence,
            min_domain_size=args.min_domain_size,
        )

        # Save contact map
        if contact_map is not None and args.contact_map_dir:
            np.save(Path(args.contact_map_dir) / f"{pid}.npy", contact_map)

        if args.format == "json":
            all_results.append(format_output(pid, seq_len, partitionings, contact_map))
        else:
            all_results.append(format_text_output(pid, seq_len, partitionings, contact_map))

    # Output
    output_file = open(args.output, "w") if args.output else sys.stdout

    if args.format == "json":
        json.dump(
            all_results if len(all_results) > 1 else all_results[0],
            output_file, indent=2,
        )
    else:
        for text in all_results:
            output_file.write(text + "\n")

    if args.output:
        output_file.close()
        logger.info(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
