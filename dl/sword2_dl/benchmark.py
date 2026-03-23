"""
Benchmark evaluation pipeline for SWORD3.

Evaluates on standard benchmarks (CATH-663, CASP targets) and compares
with published results from ChainSaw, Merizo, UniDoc, etc.

Usage:
    python -m sword2_dl.benchmark \
        --checkpoint checkpoints/checkpoint_best.pt \
        --config configs/default.yaml \
        --cath-file data/cath-domain-list.txt \
        --sequences data/cath_sequences.fasta \
        --output benchmark_results/
"""

import json
import logging
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch.cuda.amp import autocast

from .config import Config
from .metrics import (
    normalized_domain_overlap,
    domain_overlap_score,
    _partitioning_to_domains,
    _partitioning_to_assignment,
)
from .model import DomainPartitionNet
from .postprocess import predict_partitionings

logger = logging.getLogger(__name__)


# Published baseline results for comparison
PUBLISHED_RESULTS = {
    "ChainSaw": {"ndo": 0.78, "source": "Lau et al. 2024"},
    "Merizo": {"ndo": 0.72, "source": "Sheridan et al. 2024"},
    "UniDoc": {"ndo": 0.65, "source": "Wang et al. 2023"},
    "FUpred": {"ndo": 0.61, "source": "Zheng et al. 2021"},
}


@dataclass
class BenchmarkResult:
    """Result for a single benchmark protein."""
    protein_id: str
    seq_len: int
    gt_num_domains: int
    pred_num_domains: int
    ndo: float
    domain_overlap: float
    num_domains_correct: bool
    has_discontinuous: bool
    cath_class: str | None = None


def load_benchmark_targets(
    cath_file: str,
    sequence_source: str,
    test_set: str = "all",
) -> list[dict]:
    """Load benchmark targets from CATH domain annotations.

    Args:
        cath_file: Path to CATH domain list file.
        sequence_source: Path to FASTA with sequences.
        test_set: Which test set to use ("all", "cath663", or path to ID list).

    Returns:
        List of target dicts with 'id', 'sequence', 'partitionings', etc.
    """
    from .cath_dataset import parse_cath_domain_list, _load_sequences

    cath_chains = parse_cath_domain_list(cath_file)
    sequences = _load_sequences(sequence_source)

    # Filter to test set if specified
    if test_set != "all" and Path(test_set).exists():
        with open(test_set) as f:
            test_ids = set(line.strip() for line in f if line.strip())
        cath_chains = {k: v for k, v in cath_chains.items() if k in test_ids}

    targets = []
    for chain_id, domains in cath_chains.items():
        pdb_id = chain_id[:4].lower()
        chain_letter = chain_id[4]
        lookup_id = f"{pdb_id}_{chain_letter}"

        seq = sequences.get(lookup_id) or sequences.get(chain_id) or sequences.get(pdb_id)
        if seq is None:
            continue

        # Build partitioning in training format
        partitioning = []
        has_discontinuous = False
        cath_class = None
        for d in domains:
            segments = [(s, e) for s, e in d["segments"] if s < len(seq) and e < len(seq)]
            if segments:
                partitioning.append(segments)
                if len(segments) > 1:
                    has_discontinuous = True
            if "cath_class" in d and cath_class is None:
                c = d["cath_class"]
                cath_class = _cath_class_name(c[0]) if isinstance(c, tuple) else None

        if not partitioning:
            continue

        targets.append({
            "id": chain_id,
            "sequence": seq,
            "partitionings": [partitioning],
            "has_discontinuous": has_discontinuous,
            "cath_class": cath_class,
        })

    logger.info(f"Loaded {len(targets)} benchmark targets")
    return targets


def _cath_class_name(class_id: int) -> str:
    """Map CATH class number to name."""
    names = {1: "alpha", 2: "beta", 3: "alpha-beta", 4: "few-SS"}
    return names.get(class_id, f"class-{class_id}")


@torch.no_grad()
def run_benchmark(
    model: DomainPartitionNet,
    targets: list[dict],
    embedding_sources: list,
    device: torch.device,
    max_seq_len: int = 1024,
) -> list[BenchmarkResult]:
    """Run model on benchmark targets and compute metrics.

    Args:
        model: Trained model.
        targets: List of target dicts.
        embedding_sources: Embedding source configs.
        device: Compute device.
        max_seq_len: Maximum sequence length to process.

    Returns:
        List of BenchmarkResult for each target.
    """
    from safetensors.torch import load_file as load_safetensors

    model.eval()
    results = []

    for target in targets:
        seq_len = len(target["sequence"])
        if seq_len > max_seq_len:
            continue

        protein_id = target["id"]

        # Load embeddings
        embeddings = []
        skip = False
        for source in embedding_sources:
            filename = source.filename_template.format(id=protein_id)
            path = Path(source.path) / filename
            if not path.exists():
                skip = True
                break
            tensors = load_safetensors(str(path))
            emb = tensors.get(protein_id, next(iter(tensors.values()))).float()
            if emb.dim() == 3:
                emb = emb.squeeze(0)
            if emb.shape[0] > seq_len:
                emb = emb[:seq_len]
            elif emb.shape[0] < seq_len:
                pad = torch.zeros(seq_len - emb.shape[0], emb.shape[1])
                emb = torch.cat([emb, pad], dim=0)
            embeddings.append(emb)

        if skip:
            continue

        embeddings = torch.cat(embeddings, dim=-1).unsqueeze(0).to(device)
        mask = torch.ones(1, seq_len, dtype=torch.bool, device=device)

        with autocast(dtype=torch.float16):
            outputs = model(embeddings, mask)

        # Unbatch outputs
        single_outputs = {
            k: v[0] for k, v in outputs.items()
        }

        # Post-process
        partitionings = predict_partitionings(single_outputs, seq_len)

        if not partitionings:
            # Predict single domain as fallback
            pred_domains = [[(0, seq_len - 1)]]
        else:
            pred_domains = partitionings[0]["domains"]

        gt_domains = _partitioning_to_domains(target["partitionings"][0])

        # Compute metrics
        ndo = normalized_domain_overlap(pred_domains, gt_domains, seq_len)
        overlap = domain_overlap_score(pred_domains, gt_domains, seq_len)

        result = BenchmarkResult(
            protein_id=protein_id,
            seq_len=seq_len,
            gt_num_domains=len(gt_domains),
            pred_num_domains=len(pred_domains),
            ndo=ndo,
            domain_overlap=overlap,
            num_domains_correct=(len(pred_domains) == len(gt_domains)),
            has_discontinuous=target.get("has_discontinuous", False),
            cath_class=target.get("cath_class"),
        )
        results.append(result)

    return results


def analyze_results(results: list[BenchmarkResult]) -> dict:
    """Analyze benchmark results with breakdowns by category."""
    if not results:
        return {"error": "No results to analyze"}

    analysis = {}

    # Overall metrics
    ndos = [r.ndo for r in results]
    overlaps = [r.domain_overlap for r in results]
    analysis["overall"] = {
        "n_proteins": len(results),
        "mean_ndo": float(np.mean(ndos)),
        "median_ndo": float(np.median(ndos)),
        "std_ndo": float(np.std(ndos)),
        "mean_overlap": float(np.mean(overlaps)),
        "num_domains_accuracy": float(np.mean([r.num_domains_correct for r in results])),
    }

    # Breakdown by number of GT domains
    by_ndom = defaultdict(list)
    for r in results:
        key = str(r.gt_num_domains) if r.gt_num_domains <= 4 else "5+"
        by_ndom[key].append(r)
    analysis["by_num_domains"] = {
        k: {
            "n": len(v),
            "mean_ndo": float(np.mean([r.ndo for r in v])),
            "num_domains_accuracy": float(np.mean([r.num_domains_correct for r in v])),
        }
        for k, v in sorted(by_ndom.items())
    }

    # Breakdown by sequence length
    len_bins = [(0, 200), (200, 400), (400, 800), (800, float("inf"))]
    analysis["by_length"] = {}
    for lo, hi in len_bins:
        key = f"{lo}-{int(hi)}" if hi != float("inf") else f"{lo}+"
        subset = [r for r in results if lo <= r.seq_len < hi]
        if subset:
            analysis["by_length"][key] = {
                "n": len(subset),
                "mean_ndo": float(np.mean([r.ndo for r in subset])),
            }

    # Breakdown by continuous vs discontinuous
    cont = [r for r in results if not r.has_discontinuous]
    disc = [r for r in results if r.has_discontinuous]
    analysis["by_topology"] = {}
    if cont:
        analysis["by_topology"]["continuous"] = {
            "n": len(cont),
            "mean_ndo": float(np.mean([r.ndo for r in cont])),
        }
    if disc:
        analysis["by_topology"]["discontinuous"] = {
            "n": len(disc),
            "mean_ndo": float(np.mean([r.ndo for r in disc])),
        }

    # Breakdown by CATH class
    by_class = defaultdict(list)
    for r in results:
        if r.cath_class:
            by_class[r.cath_class].append(r)
    if by_class:
        analysis["by_cath_class"] = {
            k: {
                "n": len(v),
                "mean_ndo": float(np.mean([r.ndo for r in v])),
            }
            for k, v in sorted(by_class.items())
        }

    # Comparison with published results
    analysis["comparison"] = {}
    our_ndo = analysis["overall"]["mean_ndo"]
    for method, data in PUBLISHED_RESULTS.items():
        analysis["comparison"][method] = {
            "published_ndo": data["ndo"],
            "our_ndo": our_ndo,
            "delta": our_ndo - data["ndo"],
            "source": data["source"],
        }

    return analysis


def print_benchmark_report(analysis: dict) -> str:
    """Format benchmark analysis as a human-readable report."""
    lines = []
    lines.append("=" * 70)
    lines.append("SWORD3 Benchmark Report")
    lines.append("=" * 70)

    ov = analysis.get("overall", {})
    lines.append(f"\nOverall ({ov.get('n_proteins', 0)} proteins):")
    lines.append(f"  NDO:                {ov.get('mean_ndo', 0):.3f} +/- {ov.get('std_ndo', 0):.3f}")
    lines.append(f"  Domain overlap:     {ov.get('mean_overlap', 0):.3f}")
    lines.append(f"  Num domains acc:    {ov.get('num_domains_accuracy', 0):.1%}")

    lines.append("\nBy number of domains:")
    for k, v in analysis.get("by_num_domains", {}).items():
        lines.append(f"  {k:>3s} domains (n={v['n']:>4d}): NDO={v['mean_ndo']:.3f}, "
                     f"count_acc={v['num_domains_accuracy']:.1%}")

    lines.append("\nBy sequence length:")
    for k, v in analysis.get("by_length", {}).items():
        lines.append(f"  {k:>8s} (n={v['n']:>4d}): NDO={v['mean_ndo']:.3f}")

    if "by_topology" in analysis:
        lines.append("\nBy topology:")
        for k, v in analysis["by_topology"].items():
            lines.append(f"  {k:>15s} (n={v['n']:>4d}): NDO={v['mean_ndo']:.3f}")

    if "by_cath_class" in analysis:
        lines.append("\nBy CATH class:")
        for k, v in analysis["by_cath_class"].items():
            lines.append(f"  {k:>12s} (n={v['n']:>4d}): NDO={v['mean_ndo']:.3f}")

    lines.append("\nComparison with published methods:")
    lines.append(f"  {'Method':<15s} {'Published NDO':>15s} {'Our NDO':>10s} {'Delta':>10s}")
    lines.append("  " + "-" * 55)
    for method, data in analysis.get("comparison", {}).items():
        delta_str = f"{data['delta']:+.3f}"
        lines.append(f"  {method:<15s} {data['published_ndo']:>15.3f} "
                     f"{data['our_ndo']:>10.3f} {delta_str:>10s}")

    lines.append("\n" + "=" * 70)

    report = "\n".join(lines)
    return report


def main():
    """CLI entry point for benchmarking."""
    import argparse

    parser = argparse.ArgumentParser(description="Benchmark SWORD3")
    parser.add_argument("--checkpoint", required=True, help="Model checkpoint")
    parser.add_argument("--config", default=None, help="Config YAML")
    parser.add_argument("--cath-file", required=True, help="CATH domain list")
    parser.add_argument("--sequences", required=True, help="FASTA file with sequences")
    parser.add_argument("--test-set", default="all", help="Test set filter (ID list file)")
    parser.add_argument("--output", default="benchmark_results", help="Output directory")
    parser.add_argument("--use-ema", action="store_true", help="Use EMA weights from checkpoint")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Load config and model
    config = Config.from_yaml(args.config) if args.config else Config()
    model = DomainPartitionNet(config.model).to(device)

    checkpoint = torch.load(args.checkpoint, map_location=device)
    if args.use_ema and "ema" in checkpoint:
        # Load EMA weights
        for name, param in model.named_parameters():
            if name in checkpoint["ema"]["shadow"]:
                param.data.copy_(checkpoint["ema"]["shadow"][name])
        logger.info("Loaded EMA weights")
    else:
        model.load_state_dict(checkpoint["model"])
        logger.info("Loaded model weights")

    # Load benchmark targets
    targets = load_benchmark_targets(args.cath_file, args.sequences, args.test_set)

    # Run benchmark
    embedding_sources = config.model.get_embedding_sources()
    results = run_benchmark(model, targets, embedding_sources, device)

    # Analyze
    analysis = analyze_results(results)
    report = print_benchmark_report(analysis)

    # Save
    os.makedirs(args.output, exist_ok=True)
    with open(os.path.join(args.output, "analysis.json"), "w") as f:
        json.dump(analysis, f, indent=2)
    with open(os.path.join(args.output, "report.txt"), "w") as f:
        f.write(report)
    with open(os.path.join(args.output, "per_protein.json"), "w") as f:
        json.dump(
            [{"id": r.protein_id, "ndo": r.ndo, "overlap": r.domain_overlap,
              "gt_domains": r.gt_num_domains, "pred_domains": r.pred_num_domains}
             for r in results],
            f, indent=2,
        )

    print(report)
    logger.info(f"Results saved to {args.output}")


if __name__ == "__main__":
    main()
