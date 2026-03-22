"""
Inference script for DomainPartitionNet.

Takes a protein sequence (or FASTA file) and predicts alternative domain
partitionings. Outputs results in SWORD2-compatible format.
"""

import json
import logging
import sys
from pathlib import Path

import torch

from .config import Config, ModelConfig
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
    model.load_esm()

    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint["model"])
    model = model.to(device)
    model.eval()

    return model


def predict_sequence(
    model: DomainPartitionNet,
    sequence: str,
    device: torch.device | None = None,
    min_confidence: float = 0.1,
    min_domain_size: int = 20,
) -> list[dict]:
    """Predict domain partitionings for a single sequence.

    Args:
        model: Trained DomainPartitionNet.
        sequence: Amino acid sequence string.
        device: Target device.
        min_confidence: Minimum confidence to report.
        min_domain_size: Minimum domain size in residues.

    Returns:
        List of partitioning dicts sorted by confidence.
    """
    if device is None:
        device = next(model.parameters()).device

    # Tokenize
    batch_converter = model.esm_alphabet.get_batch_converter()
    _, _, tokens = batch_converter([("query", sequence)])
    tokens = tokens.to(device)

    seq_len = len(sequence)
    mask = torch.ones(1, seq_len, dtype=torch.bool, device=device)

    # Forward pass
    with torch.no_grad(), torch.cuda.amp.autocast(dtype=torch.float16):
        outputs = model(tokens, mask)

    # Extract single-sample outputs
    single_outputs = {
        "co_membership": outputs["co_membership"][0],  # (K, L, L)
        "confidence": outputs["confidence"][0],  # (K,)
    }
    if "num_domains_logits" in outputs:
        single_outputs["num_domains_logits"] = outputs["num_domains_logits"][0]
    if "boundary_logits" in outputs:
        single_outputs["boundary_logits"] = outputs["boundary_logits"][0]

    # Post-process
    partitionings = predict_partitionings(
        single_outputs,
        seq_len=seq_len,
        min_confidence=min_confidence,
        min_domain_size=min_domain_size,
    )

    return partitionings


def format_output(
    protein_id: str,
    sequence: str,
    partitionings: list[dict],
) -> dict:
    """Format predictions in SWORD2-compatible JSON output."""
    output = {
        "id": protein_id,
        "sequence_length": len(sequence),
        "method": "SWORD2-DL",
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
            # Format segments as "start-end" strings (1-indexed for compatibility)
            segment_strs = [
                f"{s + 1}-{e + 1}" for s, e in domain_segments
            ]
            partition_data["Domains"][domain_key] = {
                "segments": segment_strs,
                "residue_count": sum(e - s + 1 for s, e in domain_segments),
            }

        output[key] = partition_data

    return output


def format_text_output(
    protein_id: str,
    sequence: str,
    partitionings: list[dict],
) -> str:
    """Format predictions as human-readable text (SWORD2-style)."""
    lines = [
        f"SWORD2-DL Domain Prediction",
        f"Protein: {protein_id}",
        f"Sequence length: {len(sequence)}",
        f"Number of alternative partitionings: {len(partitionings)}",
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


def parse_fasta(fasta_path: str) -> list[tuple[str, str]]:
    """Parse FASTA file into list of (id, sequence) tuples."""
    sequences = []
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_id is not None:
        sequences.append((current_id, "".join(current_seq)))

    return sequences


def main():
    """Entry point for prediction."""
    import argparse

    parser = argparse.ArgumentParser(description="Predict protein domains with SWORD2-DL")
    parser.add_argument(
        "--checkpoint",
        type=str,
        required=True,
        help="Path to model checkpoint",
    )
    parser.add_argument(
        "--sequence",
        type=str,
        default=None,
        help="Single amino acid sequence",
    )
    parser.add_argument(
        "--fasta",
        type=str,
        default=None,
        help="Path to FASTA file",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file path (default: stdout)",
    )
    parser.add_argument(
        "--format",
        choices=["json", "text"],
        default="text",
        help="Output format",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to config YAML",
    )
    parser.add_argument(
        "--min-confidence",
        type=float,
        default=0.1,
        help="Minimum confidence threshold",
    )
    parser.add_argument(
        "--min-domain-size",
        type=int,
        default=20,
        help="Minimum domain size in residues",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Device (cuda/cpu)",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    if args.sequence is None and args.fasta is None:
        parser.error("Must provide --sequence or --fasta")

    # Load config
    if args.config:
        config = Config.from_yaml(args.config)
    else:
        config = Config()

    # Device
    if args.device:
        device = torch.device(args.device)
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Load model
    logger.info(f"Loading model from {args.checkpoint}")
    model = load_model(args.checkpoint, config.model, device)

    # Get sequences
    if args.sequence:
        sequences = [("query", args.sequence)]
    else:
        sequences = parse_fasta(args.fasta)

    logger.info(f"Predicting domains for {len(sequences)} sequence(s)")

    # Predict
    all_results = []
    for prot_id, seq in sequences:
        partitionings = predict_sequence(
            model,
            seq,
            device=device,
            min_confidence=args.min_confidence,
            min_domain_size=args.min_domain_size,
        )

        if args.format == "json":
            result = format_output(prot_id, seq, partitionings)
            all_results.append(result)
        else:
            text = format_text_output(prot_id, seq, partitionings)
            all_results.append(text)

    # Output
    output_file = open(args.output, "w") if args.output else sys.stdout

    if args.format == "json":
        json.dump(all_results if len(all_results) > 1 else all_results[0], output_file, indent=2)
    else:
        for text in all_results:
            output_file.write(text + "\n")

    if args.output:
        output_file.close()
        logger.info(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
