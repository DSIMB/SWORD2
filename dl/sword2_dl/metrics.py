"""
Evaluation metrics for domain partitioning quality.

Standard metrics:
- NDO (Normalized Domain Overlap): CATH standard metric
- Boundary F1 with tolerance: domain boundary accuracy
- Number of domains accuracy
- V-measure / NMI: clustering quality
- Co-membership accuracy
- Domain overlap score (Hungarian-matched IoU)
"""

import numpy as np
import torch
import torch.nn.functional as F
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import normalized_mutual_info_score, v_measure_score

from .losses import build_co_membership_matrix, build_boundary_targets


@torch.no_grad()
def compute_metrics(
    outputs: dict[str, torch.Tensor],
    targets: list[dict],
) -> dict[str, float]:
    """Compute evaluation metrics for a batch.

    Args:
        outputs: Model outputs (co_membership, confidence, etc.)
        targets: List of target dicts with partitionings.

    Returns:
        Dictionary of metric names -> values (averaged over batch).
    """
    B = outputs["co_membership"].shape[0]
    device = outputs["co_membership"].device

    metrics = {
        "co_membership_acc": 0.0,
        "boundary_f1": 0.0,
        "boundary_f1_tol10": 0.0,
        "num_domains_acc": 0.0,
        "num_domains_acc_tol1": 0.0,
        "ndo": 0.0,
        "v_measure": 0.0,
        "nmi": 0.0,
    }
    n_valid = 0

    for b in range(B):
        seq_len = targets[b]["seq_len"]
        partitionings = targets[b]["partitionings"]

        if not partitionings:
            continue

        n_valid += 1

        # Use the first (optimal) partitioning as primary target
        gt_part = partitionings[0]
        gt_co_mem = build_co_membership_matrix(gt_part, seq_len, device)

        # Find best matching predicted slot by confidence
        pred_co_mem = outputs["co_membership"][b, :, :seq_len, :seq_len]  # (K, L, L)
        pred_conf = outputs["confidence"][b]  # (K,)

        best_k = pred_conf.argmax().item()
        pred_probs = torch.sigmoid(pred_co_mem[best_k])  # (L, L)

        # Co-membership accuracy (upper triangle)
        pred_binary = (pred_probs > 0.5).float()
        correct = (pred_binary == gt_co_mem).float()
        mask = torch.triu(torch.ones(seq_len, seq_len, device=device), diagonal=1)
        acc = (correct * mask).sum() / mask.sum()
        metrics["co_membership_acc"] += acc.item()

        # Boundary F1 (strict and tolerant)
        if "boundary_logits" in outputs:
            gt_boundaries = build_boundary_targets(gt_part, seq_len, device)
            pred_boundaries = (
                torch.sigmoid(outputs["boundary_logits"][b, best_k, :seq_len]) > 0.5
            ).float()

            f1_strict = _boundary_f1(pred_boundaries, gt_boundaries, tolerance=0)
            f1_tol10 = _boundary_f1(pred_boundaries, gt_boundaries, tolerance=10)
            metrics["boundary_f1"] += f1_strict
            metrics["boundary_f1_tol10"] += f1_tol10

        # Number of domains accuracy
        if "num_domains_logits" in outputs:
            pred_nd = outputs["num_domains_logits"][b, best_k].argmax().item() + 1
            gt_nd = len(gt_part)
            metrics["num_domains_acc"] += float(pred_nd == gt_nd)
            metrics["num_domains_acc_tol1"] += float(abs(pred_nd - gt_nd) <= 1)

        # Extract domain assignments for NDO/V-measure/NMI
        from .postprocess import extract_partitioning
        pred_part = extract_partitioning(
            pred_probs.cpu().numpy(),
            confidence=pred_conf[best_k].item(),
        )
        pred_assignment = pred_part["assignment"]
        gt_assignment = _partitioning_to_assignment(gt_part, seq_len)

        # NDO
        pred_domains = pred_part["domains"]
        gt_domains = _partitioning_to_domains(gt_part)
        ndo = normalized_domain_overlap(pred_domains, gt_domains, seq_len)
        metrics["ndo"] += ndo

        # V-measure and NMI
        metrics["v_measure"] += v_measure_score(gt_assignment, pred_assignment)
        metrics["nmi"] += normalized_mutual_info_score(gt_assignment, pred_assignment)

    if n_valid > 0:
        for k in metrics:
            metrics[k] /= n_valid

    return metrics


def _boundary_f1(
    pred: torch.Tensor,
    gt: torch.Tensor,
    tolerance: int = 0,
) -> float:
    """Compute boundary F1 with optional tolerance.

    With tolerance > 0, a predicted boundary is considered a true positive
    if it falls within +/- tolerance residues of any ground truth boundary.
    """
    pred_positions = pred.nonzero(as_tuple=True)[0].cpu().numpy()
    gt_positions = gt.nonzero(as_tuple=True)[0].cpu().numpy()

    if len(gt_positions) == 0:
        return 1.0 if len(pred_positions) == 0 else 0.0
    if len(pred_positions) == 0:
        return 0.0

    # Count true positives with tolerance
    tp = 0
    matched_gt = set()
    for p in pred_positions:
        for gi, g in enumerate(gt_positions):
            if gi not in matched_gt and abs(int(p) - int(g)) <= tolerance:
                tp += 1
                matched_gt.add(gi)
                break

    precision = tp / len(pred_positions) if len(pred_positions) > 0 else 0
    recall = tp / len(gt_positions) if len(gt_positions) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall + 1e-8) if (precision + recall) > 0 else 0
    return f1


def _partitioning_to_assignment(
    partitioning: list[list[tuple[int, int]]],
    seq_len: int,
) -> np.ndarray:
    """Convert domain partitioning to per-residue assignment array."""
    assignment = np.zeros(seq_len, dtype=int)
    for domain_idx, segments in enumerate(partitioning):
        for start, end in segments:
            start = max(0, start)
            end = min(seq_len - 1, end)
            assignment[start:end + 1] = domain_idx + 1
    return assignment


def _partitioning_to_domains(
    partitioning: list[list[tuple[int, int]]],
) -> list[list[tuple[int, int]]]:
    """Convert partitioning format to list of domains (list of segments)."""
    return [[(s, e) for s, e in segments] for segments in partitioning]


def normalized_domain_overlap(
    pred_domains: list[list[tuple[int, int]]],
    gt_domains: list[list[tuple[int, int]]],
    seq_len: int,
) -> float:
    """Compute Normalized Domain Overlap (NDO) — the CATH standard metric.

    For each ground truth domain, finds the best-matching predicted domain,
    computes overlap ratio. Returns average over all ground truth domains.

    This is the primary metric used to compare domain prediction methods
    (e.g., ChainSaw reports 78% NDO on CATH benchmarks).

    Args:
        pred_domains: List of predicted domains, each is list of (start, end) segments.
        gt_domains: List of ground truth domains.
        seq_len: Sequence length.

    Returns:
        NDO score in [0, 1].
    """
    if not gt_domains:
        return 1.0 if not pred_domains else 0.0
    if not pred_domains:
        return 0.0

    def segments_to_set(segments):
        residues = set()
        for start, end in segments:
            residues.update(range(start, end + 1))
        return residues

    pred_sets = [segments_to_set(d) for d in pred_domains]
    gt_sets = [segments_to_set(d) for d in gt_domains]

    # For each GT domain, find best overlap with any predicted domain
    total_overlap = 0.0
    for gt_set in gt_sets:
        if not gt_set:
            continue
        best_overlap = 0.0
        for pred_set in pred_sets:
            intersection = len(gt_set & pred_set)
            # NDO uses overlap ratio: |intersection| / |gt_domain|
            overlap = intersection / len(gt_set) if gt_set else 0
            best_overlap = max(best_overlap, overlap)
        total_overlap += best_overlap

    ndo = total_overlap / len(gt_domains)

    # Penalty for predicting wrong number of domains
    n_pred = len(pred_domains)
    n_gt = len(gt_domains)
    if n_pred > n_gt:
        # Penalize over-segmentation
        ndo *= n_gt / n_pred

    return ndo


def domain_overlap_score(
    pred_domains: list[list[tuple[int, int]]],
    gt_domains: list[list[tuple[int, int]]],
    seq_len: int,
) -> float:
    """Compute domain overlap score using Hungarian matching.

    Each domain is represented as a list of (start, end) segments.
    Computes IoU between predicted and ground truth domains,
    finds optimal matching, and returns average IoU.

    Args:
        pred_domains: Predicted domains.
        gt_domains: Ground truth domains.
        seq_len: Sequence length.

    Returns:
        Average IoU of best-matched domains.
    """
    n_pred = len(pred_domains)
    n_gt = len(gt_domains)

    if n_pred == 0 or n_gt == 0:
        return 0.0

    # Convert to residue sets
    def segments_to_set(segments):
        residues = set()
        for start, end in segments:
            residues.update(range(start, end + 1))
        return residues

    pred_sets = [segments_to_set(d) for d in pred_domains]
    gt_sets = [segments_to_set(d) for d in gt_domains]

    # Build IoU matrix
    n = max(n_pred, n_gt)
    iou_matrix = np.zeros((n, n))

    for i in range(n_pred):
        for j in range(n_gt):
            intersection = len(pred_sets[i] & gt_sets[j])
            union = len(pred_sets[i] | gt_sets[j])
            iou_matrix[i, j] = intersection / max(1, union)

    # Hungarian matching (maximize IoU = minimize -IoU)
    row_ind, col_ind = linear_sum_assignment(-iou_matrix)

    total_iou = iou_matrix[row_ind, col_ind].sum()
    return total_iou / max(n_pred, n_gt)
