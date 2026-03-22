"""
Evaluation metrics for domain partitioning quality.

Metrics:
- Domain overlap score (best Hungarian-matched IoU)
- Boundary F1 (precision/recall of domain boundaries)
- Number of domains accuracy
- Co-membership matrix AUC
"""

import numpy as np
import torch
import torch.nn.functional as F
from scipy.optimize import linear_sum_assignment

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
        "co_membership_auc": 0.0,
        "boundary_f1": 0.0,
        "num_domains_acc": 0.0,
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

        # Find best matching predicted slot
        pred_co_mem = outputs["co_membership"][b, :, :seq_len, :seq_len]  # (K, L, L)
        pred_conf = outputs["confidence"][b]  # (K,)

        best_k = pred_conf.argmax().item()
        pred_probs = torch.sigmoid(pred_co_mem[best_k])  # (L, L)

        # Co-membership AUC (approximate via thresholded accuracy)
        pred_binary = (pred_probs > 0.5).float()
        correct = (pred_binary == gt_co_mem).float()
        # Only count upper triangle (symmetric matrix)
        mask = torch.triu(torch.ones(seq_len, seq_len, device=device), diagonal=1)
        auc_approx = (correct * mask).sum() / mask.sum()
        metrics["co_membership_auc"] += auc_approx.item()

        # Boundary F1
        if "boundary_logits" in outputs:
            gt_boundaries = build_boundary_targets(gt_part, seq_len, device)
            pred_boundaries = (
                torch.sigmoid(outputs["boundary_logits"][b, best_k, :seq_len]) > 0.5
            ).float()

            tp = (pred_boundaries * gt_boundaries).sum()
            fp = (pred_boundaries * (1 - gt_boundaries)).sum()
            fn = ((1 - pred_boundaries) * gt_boundaries).sum()

            precision = tp / (tp + fp + 1e-8)
            recall = tp / (tp + fn + 1e-8)
            f1 = 2 * precision * recall / (precision + recall + 1e-8)
            metrics["boundary_f1"] += f1.item()

        # Number of domains accuracy
        if "num_domains_logits" in outputs:
            pred_nd = outputs["num_domains_logits"][b, best_k].argmax().item() + 1
            gt_nd = len(gt_part)
            metrics["num_domains_acc"] += float(pred_nd == gt_nd)

    if n_valid > 0:
        for k in metrics:
            metrics[k] /= n_valid

    return metrics


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
