"""
End-to-end prediction pipeline: sequence → domain partitionings.

Single entry point that handles embedding computation, model inference,
and post-processing. Users only need a FASTA file or raw sequence.

Usage:
    # From sequence
    python -m sword2_dl.pipeline --sequence "MVLSPADKTN..." --checkpoint best.pt

    # From FASTA
    python -m sword2_dl.pipeline --fasta proteins.fasta --checkpoint best.pt

    # With pre-computed embeddings (faster)
    python -m sword2_dl.pipeline --embeddings P12345.safetensors --checkpoint best.pt
"""

import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import torch
from torch.cuda.amp import autocast

from .config import Config
from .model import DomainPartitionNet
from .postprocess import predict_partitionings

logger = logging.getLogger(__name__)


def compute_embeddings_esm2(
    sequences: dict[str, str],
    model_name: str = "esm2_t36_3B_UR50D",
    device: torch.device | None = None,
    batch_size: int = 4,
) -> dict[str, torch.Tensor]:
    """Compute ESM-2 embeddings for sequences on-the-fly.

    Args:
        sequences: Dict of protein_id -> sequence string.
        model_name: ESM-2 model to use.
        device: Compute device.
        batch_size: Batch size for embedding computation.

    Returns:
        Dict of protein_id -> (L, D) embedding tensor.
    """
    try:
        import esm
    except ImportError:
        logger.error(
            "ESM-2 not installed. Install with: pip install fair-esm\n"
            "Or provide pre-computed embeddings with --embeddings"
        )
        sys.exit(1)

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    logger.info(f"Loading ESM-2 model: {model_name}")
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    embeddings = {}
    items = list(sequences.items())

    for i in range(0, len(items), batch_size):
        batch_items = items[i:i + batch_size]
        data = [(pid, seq) for pid, seq in batch_items]

        _, _, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33])

        for j, (pid, seq) in enumerate(batch_items):
            # Remove BOS/EOS tokens
            emb = results["representations"][33][j, 1:len(seq) + 1].cpu()
            embeddings[pid] = emb

    return embeddings


def load_model(
    checkpoint_path: str,
    config: Config | None = None,
    device: torch.device | None = None,
    use_ema: bool = True,
) -> tuple[DomainPartitionNet, Config, torch.device]:
    """Load a trained model from checkpoint.

    Args:
        checkpoint_path: Path to checkpoint file.
        config: Config to use. If None, tries to load from checkpoint.
        device: Compute device. Auto-detects if None.
        use_ema: Whether to use EMA weights if available.

    Returns:
        Tuple of (model, config, device).
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if config is None:
        config = Config()

    model = DomainPartitionNet(config.model).to(device)

    checkpoint = torch.load(checkpoint_path, map_location=device)
    if use_ema and "ema" in checkpoint:
        for name, param in model.named_parameters():
            if name in checkpoint["ema"]["shadow"]:
                param.data.copy_(checkpoint["ema"]["shadow"][name])
        logger.info("Loaded EMA weights")
    else:
        model.load_state_dict(checkpoint["model"])

    model.eval()
    return model, config, device


def predict_from_embeddings(
    model: DomainPartitionNet,
    embeddings: torch.Tensor,
    device: torch.device,
    min_confidence: float = 0.1,
    min_domain_size: int = 20,
) -> dict:
    """Run prediction on pre-computed embeddings.

    Args:
        model: Trained model.
        embeddings: (L, D) embedding tensor.
        device: Compute device.
        min_confidence: Minimum confidence for output partitionings.
        min_domain_size: Minimum domain size in residues.

    Returns:
        Dict with 'partitionings', 'contact_map', and metadata.
    """
    seq_len = embeddings.shape[0]
    emb = embeddings.unsqueeze(0).to(device)
    mask = torch.ones(1, seq_len, dtype=torch.bool, device=device)

    with torch.no_grad(), autocast(dtype=torch.float16):
        outputs = model(emb, mask)

    # Unbatch
    single_outputs = {k: v[0] for k, v in outputs.items()}

    # Post-process partitionings
    partitionings = predict_partitionings(
        single_outputs, seq_len,
        min_confidence=min_confidence,
        min_domain_size=min_domain_size,
    )

    result = {"partitionings": []}
    for part in partitionings:
        result["partitionings"].append({
            "domains": [
                [{"start": s, "end": e} for s, e in segments]
                for segments in part["domains"]
            ],
            "num_domains": part["num_domains"],
            "confidence": float(part["confidence"]),
        })

    # Contact map
    if "contact_map_logits" in single_outputs:
        contact_probs = torch.sigmoid(
            single_outputs["contact_map_logits"][:seq_len, :seq_len]
        ).cpu().numpy()
        result["contact_map"] = contact_probs

    return result


def format_output(
    protein_id: str,
    result: dict,
    output_format: str = "text",
) -> str:
    """Format prediction results for output.

    Args:
        protein_id: Protein identifier.
        result: Prediction result dict.
        output_format: "text", "json", or "cath".
    """
    if output_format == "json":
        output = {
            "id": protein_id,
            "partitionings": result["partitionings"],
        }
        return json.dumps(output, indent=2)

    elif output_format == "cath":
        # CATH-style format: one line per domain with segments
        lines = []
        for i, part in enumerate(result["partitionings"]):
            lines.append(f"# Partitioning {i + 1} (confidence: {part['confidence']:.3f})")
            for j, domain in enumerate(part["domains"]):
                segs = ", ".join(f"{s['start']+1}-{s['end']+1}" for s in domain)
                lines.append(f"  Domain {j + 1}: {segs}")
        return "\n".join(lines)

    else:  # text
        lines = [f"Protein: {protein_id}"]
        lines.append(f"Number of alternative partitionings: {len(result['partitionings'])}")
        lines.append("")
        for i, part in enumerate(result["partitionings"]):
            lines.append(
                f"Partitioning {i + 1}: {part['num_domains']} domains "
                f"(confidence: {part['confidence']:.3f})"
            )
            for j, domain in enumerate(part["domains"]):
                segs = " + ".join(f"{s['start']+1}-{s['end']+1}" for s in domain)
                lines.append(f"  Domain {j + 1}: {segs}")
            lines.append("")
        return "\n".join(lines)


def parse_fasta(fasta_path: str) -> dict[str, str]:
    """Parse FASTA file into dict of id -> sequence."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                header = line[1:]
                parts = header.split("|")
                current_id = parts[1] if len(parts) >= 3 else header.split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences


def main():
    """CLI entry point for end-to-end prediction."""
    import argparse

    parser = argparse.ArgumentParser(
        description="SWORD3: Predict protein domain partitionings from sequence"
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--sequence", type=str, help="Raw protein sequence")
    input_group.add_argument("--fasta", type=str, help="FASTA file with sequences")
    input_group.add_argument(
        "--embeddings", type=str,
        help="Pre-computed embeddings (safetensors file)"
    )

    parser.add_argument("--checkpoint", required=True, help="Model checkpoint path")
    parser.add_argument("--config", default=None, help="Config YAML file")
    parser.add_argument("--output", default=None, help="Output file/directory")
    parser.add_argument(
        "--format", choices=["text", "json", "cath"], default="text",
        help="Output format"
    )
    parser.add_argument("--min-confidence", type=float, default=0.1)
    parser.add_argument("--min-domain-size", type=int, default=20)
    parser.add_argument("--save-contact-map", action="store_true",
                        help="Save contact map as .npy")
    parser.add_argument("--no-ema", action="store_true",
                        help="Don't use EMA weights even if available")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Load model
    config = Config.from_yaml(args.config) if args.config else Config()
    model, config, device = load_model(
        args.checkpoint, config, use_ema=not args.no_ema
    )

    embedding_sources = config.model.get_embedding_sources()

    # Determine input sequences
    if args.sequence:
        sequences = {"query": args.sequence}
    elif args.fasta:
        sequences = parse_fasta(args.fasta)
    elif args.embeddings:
        # Load pre-computed embeddings directly
        from safetensors.torch import load_file as load_safetensors
        tensors = load_safetensors(args.embeddings)
        for pid, emb in tensors.items():
            result = predict_from_embeddings(
                model, emb.float(), device,
                min_confidence=args.min_confidence,
                min_domain_size=args.min_domain_size,
            )
            output_text = format_output(pid, result, args.format)
            if args.output:
                os.makedirs(args.output, exist_ok=True)
                out_path = os.path.join(args.output, f"{pid}.{args.format}")
                with open(out_path, "w") as f:
                    f.write(output_text)
                if args.save_contact_map and "contact_map" in result:
                    np.save(os.path.join(args.output, f"{pid}_contacts.npy"),
                            result["contact_map"])
            else:
                print(output_text)
        return

    # Compute embeddings on-the-fly for sequence/fasta input
    logger.info(f"Computing embeddings for {len(sequences)} sequences...")
    all_embeddings = compute_embeddings_esm2(sequences, device=device)

    # Run predictions
    for pid, seq in sequences.items():
        if pid not in all_embeddings:
            logger.warning(f"No embeddings for {pid}, skipping")
            continue

        emb = all_embeddings[pid]

        # If using multi-PLM config with only ESM-2 computed,
        # pad missing dims with zeros
        expected_dim = config.model.total_embed_dim
        if emb.shape[-1] < expected_dim:
            pad = torch.zeros(emb.shape[0], expected_dim - emb.shape[-1])
            emb = torch.cat([emb, pad], dim=-1)

        result = predict_from_embeddings(
            model, emb, device,
            min_confidence=args.min_confidence,
            min_domain_size=args.min_domain_size,
        )

        output_text = format_output(pid, result, args.format)

        if args.output:
            os.makedirs(args.output, exist_ok=True)
            ext = "json" if args.format == "json" else "txt"
            out_path = os.path.join(args.output, f"{pid}.{ext}")
            with open(out_path, "w") as f:
                f.write(output_text)
            if args.save_contact_map and "contact_map" in result:
                np.save(os.path.join(args.output, f"{pid}_contacts.npy"),
                        result["contact_map"])
            logger.info(f"Saved {out_path}")
        else:
            print(output_text)


if __name__ == "__main__":
    main()
