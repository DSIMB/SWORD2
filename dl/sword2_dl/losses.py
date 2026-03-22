"""
Loss functions for DomainPartitionNet.

Uses Hungarian matching to assign predicted partitioning slots to ground truth
partitionings, then computes losses on matched pairs.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.optimize import linear_sum_assignment

from .config import TrainConfig


def build_co_membership_matrix(
    domain_assignments: list[list[list[tuple[int, int]]]],
    seq_len: int,
    device: torch.device,
) -> torch.Tensor:
    """Build binary co-membership matrix from domain assignments.

    Args:
        domain_assignments: List of domains, each domain is a list of (start, end)
            segments (0-indexed, inclusive).
        seq_len: Sequence length.
        device: Target device.

    Returns:
        (seq_len, seq_len) binary matrix where 1 = same domain.
    """
    # Assign each residue to a domain
    assignment = torch.zeros(seq_len, dtype=torch.long, device=device)
    for domain_idx, segments in enumerate(domain_assignments):
        for start, end in segments:
            assignment[start : end + 1] = domain_idx + 1  # 1-indexed domains

    # Build co-membership: same non-zero domain
    co_mem = (
        (assignment.unsqueeze(0) == assignment.unsqueeze(1))
        & (assignment.unsqueeze(0) > 0)
        & (assignment.unsqueeze(1) > 0)
    ).float()
    return co_mem


def build_boundary_targets(
    domain_assignments: list[list[list[tuple[int, int]]]],
    seq_len: int,
    device: torch.device,
) -> torch.Tensor:
    """Build per-residue boundary target.

    A residue is a boundary if it is at the start or end of any domain segment.

    Returns:
        (seq_len,) binary vector.
    """
    boundaries = torch.zeros(seq_len, device=device)
    for segments in domain_assignments:
        for start, end in segments:
            if 0 <= start < seq_len:
                boundaries[start] = 1.0
            if 0 <= end < seq_len:
                boundaries[end] = 1.0
    return boundaries


@torch.no_grad()
def hungarian_match(
    pred_co_mem: torch.Tensor,
    gt_co_mems: list[torch.Tensor],
    pred_conf: torch.Tensor,
) -> list[tuple[int, int]]:
    """Match predicted slots to ground truth partitionings via Hungarian algorithm.

    Args:
        pred_co_mem: (K, L, L) predicted co-membership logits
        gt_co_mems: list of M ground truth (L, L) binary co-membership matrices
        pred_conf: (K,) predicted confidence scores

    Returns:
        List of (pred_idx, gt_idx) matched pairs.
    """
    K = pred_co_mem.shape[0]
    M = len(gt_co_mems)

    # Compute cost matrix: BCE between each pred-gt pair
    cost = torch.zeros(K, M, device=pred_co_mem.device)
    pred_probs = torch.sigmoid(pred_co_mem)

    for j, gt in enumerate(gt_co_mems):
        # Binary cross-entropy cost
        bce = F.binary_cross_entropy(
            pred_probs,
            gt.unsqueeze(0).expand(K, -1, -1),
            reduction="none",
        ).mean(dim=(-2, -1))  # (K,)
        cost[:, j] = bce

    # Solve assignment
    cost_np = cost.cpu().numpy()
    row_ind, col_ind = linear_sum_assignment(cost_np)

    return list(zip(row_ind.tolist(), col_ind.tolist()))


class PartitioningLoss(nn.Module):
    """Combined loss for multi-partitioning domain prediction.

    Components:
    1. Co-membership BCE: matched prediction vs ground truth
    2. Confidence calibration: matched slots should have high confidence
    3. Number of domains: cross-entropy on domain count
    4. Boundary prediction: focal loss on boundary positions
    """

    def __init__(self, config: TrainConfig):
        super().__init__()
        self.config = config

    def focal_loss(
        self,
        logits: torch.Tensor,
        targets: torch.Tensor,
        alpha: float = 0.25,
        gamma: float = 2.0,
    ) -> torch.Tensor:
        """Focal loss for imbalanced boundary prediction."""
        probs = torch.sigmoid(logits)
        bce = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
        p_t = probs * targets + (1 - probs) * (1 - targets)
        focal_weight = (1 - p_t) ** gamma
        alpha_t = alpha * targets + (1 - alpha) * (1 - targets)
        return (alpha_t * focal_weight * bce).mean()

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        targets: list[dict],
    ) -> dict[str, torch.Tensor]:
        """Compute total loss.

        Args:
            outputs: Model outputs dict with keys:
                - co_membership: (B, K, L, L)
                - confidence: (B, K)
                - num_domains_logits: (B, K, max_domains) optional
                - boundary_logits: (B, K, L) optional
            targets: List of B target dicts, each with:
                - partitionings: list of domain assignments
                  (each is list of domains, each domain is list of (start,end) segments)
                - seq_len: int

        Returns:
            Dict of loss components and total loss.
        """
        B = outputs["co_membership"].shape[0]
        K = outputs["co_membership"].shape[1]
        device = outputs["co_membership"].device

        total_co_mem_loss = torch.tensor(0.0, device=device)
        total_conf_loss = torch.tensor(0.0, device=device)
        total_num_dom_loss = torch.tensor(0.0, device=device)
        total_boundary_loss = torch.tensor(0.0, device=device)
        num_matches = 0

        for b in range(B):
            seq_len = targets[b]["seq_len"]
            partitionings = targets[b]["partitionings"]

            if not partitionings:
                continue

            # Build ground truth co-membership matrices
            gt_co_mems = []
            gt_num_domains = []
            gt_boundaries = []
            for part in partitionings:
                gt_co_mems.append(
                    build_co_membership_matrix(part, seq_len, device)
                )
                gt_num_domains.append(len(part))
                gt_boundaries.append(
                    build_boundary_targets(part, seq_len, device)
                )

            # Crop predictions to actual sequence length
            pred_co_mem = outputs["co_membership"][b, :, :seq_len, :seq_len]  # (K, L, L)
            pred_conf = outputs["confidence"][b]  # (K,)

            # Hungarian matching
            matches = hungarian_match(pred_co_mem, gt_co_mems, pred_conf)

            # Compute losses for matched pairs
            matched_indices = set()
            for pred_idx, gt_idx in matches:
                matched_indices.add(pred_idx)

                # Co-membership BCE
                co_mem_loss = F.binary_cross_entropy_with_logits(
                    pred_co_mem[pred_idx], gt_co_mems[gt_idx]
                )
                total_co_mem_loss = total_co_mem_loss + co_mem_loss

                # Confidence: matched slots should be confident
                total_conf_loss = total_conf_loss + F.binary_cross_entropy(
                    pred_conf[pred_idx].unsqueeze(0),
                    torch.ones(1, device=device),
                )

                # Number of domains
                if "num_domains_logits" in outputs:
                    nd_logits = outputs["num_domains_logits"][b, pred_idx]
                    nd_target = torch.tensor(
                        min(gt_num_domains[gt_idx], nd_logits.shape[0]) - 1,
                        device=device,
                    )
                    total_num_dom_loss = total_num_dom_loss + F.cross_entropy(
                        nd_logits.unsqueeze(0), nd_target.unsqueeze(0)
                    )

                # Boundary prediction
                if "boundary_logits" in outputs:
                    bd_logits = outputs["boundary_logits"][b, pred_idx, :seq_len]
                    total_boundary_loss = total_boundary_loss + self.focal_loss(
                        bd_logits, gt_boundaries[gt_idx]
                    )

                num_matches += 1

            # Unmatched slots: should have low confidence
            for k in range(K):
                if k not in matched_indices:
                    total_conf_loss = total_conf_loss + F.binary_cross_entropy(
                        pred_conf[k].unsqueeze(0),
                        torch.zeros(1, device=device),
                    )

        # Average over matches
        if num_matches > 0:
            total_co_mem_loss = total_co_mem_loss / num_matches
            total_conf_loss = total_conf_loss / (B * K)
            total_num_dom_loss = total_num_dom_loss / num_matches
            total_boundary_loss = total_boundary_loss / num_matches

        # Weighted sum
        cfg = self.config
        total = (
            cfg.co_membership_weight * total_co_mem_loss
            + cfg.confidence_weight * total_conf_loss
            + cfg.num_domains_weight * total_num_dom_loss
            + cfg.boundary_weight * total_boundary_loss
        )

        return {
            "loss": total,
            "co_membership_loss": total_co_mem_loss,
            "confidence_loss": total_conf_loss,
            "num_domains_loss": total_num_dom_loss,
            "boundary_loss": total_boundary_loss,
        }
