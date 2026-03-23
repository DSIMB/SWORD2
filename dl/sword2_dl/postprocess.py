"""
Post-processing: extract domain partitionings from predicted co-membership matrices.

Uses spectral clustering on each co-membership matrix to assign residues to domains.
Handles discontinuous domains naturally since clustering is based on pairwise similarity.
"""

import numpy as np
import torch
from scipy.ndimage import label as connected_components
from sklearn.cluster import SpectralClustering


def extract_partitioning(
    co_membership: np.ndarray,
    confidence: float,
    num_domains_hint: int | None = None,
    min_domain_size: int = 20,
    max_domains: int = 20,
) -> dict:
    """Extract domain partitioning from a co-membership probability matrix.

    Args:
        co_membership: (L, L) symmetric probability matrix in [0, 1].
        confidence: Scalar confidence score for this partitioning.
        num_domains_hint: If provided, use this as number of clusters.
        min_domain_size: Minimum residues per domain.
        max_domains: Maximum number of domains.

    Returns:
        Dictionary with:
            - domains: list of domains, each is a list of (start, end) segments
            - num_domains: number of domains
            - confidence: confidence score
            - assignment: per-residue domain assignment array
    """
    L = co_membership.shape[0]

    # Determine number of domains
    if num_domains_hint is not None:
        n_clusters = max(1, min(num_domains_hint, max_domains))
    else:
        n_clusters = estimate_num_domains(co_membership, max_domains)

    if n_clusters <= 1:
        # Single domain
        return {
            "domains": [[(0, L - 1)]],
            "num_domains": 1,
            "confidence": confidence,
            "assignment": np.zeros(L, dtype=int),
        }

    # Spectral clustering
    # Ensure matrix is valid for clustering
    affinity = np.clip(co_membership, 0, 1)

    # Make sure affinity is symmetric (don't override diagonal —
    # let the model's predicted self-similarity be used directly)
    affinity = (affinity + affinity.T) / 2

    try:
        sc = SpectralClustering(
            n_clusters=n_clusters,
            affinity="precomputed",
            assign_labels="kmeans",
            random_state=42,
            n_init=10,
        )
        labels = sc.fit_predict(affinity)
    except Exception:
        # Fallback: threshold-based assignment
        labels = threshold_assignment(co_membership, n_clusters)

    # Merge small domains into nearest neighbor
    labels = merge_small_domains(labels, co_membership, min_domain_size)

    # Extract segments per domain
    domains = labels_to_segments(labels)

    return {
        "domains": domains,
        "num_domains": len(domains),
        "confidence": confidence,
        "assignment": labels,
    }


def estimate_num_domains(co_membership: np.ndarray, max_domains: int) -> int:
    """Estimate number of domains from eigenvalue gap of co-membership matrix.

    Uses the eigengap heuristic: the number of clusters is determined by
    the largest gap between consecutive eigenvalues of the Laplacian.
    """
    L = co_membership.shape[0]
    if L < 2:
        return 1

    affinity = np.clip(co_membership, 0, 1)
    affinity = (affinity + affinity.T) / 2

    # Compute normalized Laplacian eigenvalues
    degree = affinity.sum(axis=1)
    degree[degree == 0] = 1  # avoid division by zero
    d_inv_sqrt = 1.0 / np.sqrt(degree)
    laplacian = np.eye(L) - (d_inv_sqrt[:, None] * affinity * d_inv_sqrt[None, :])

    # Only need first few eigenvalues
    k = min(max_domains + 1, L)
    try:
        from scipy.sparse.linalg import eigsh
        eigenvalues = eigsh(laplacian, k=k, which="SM", return_eigenvectors=False)
        eigenvalues = np.sort(eigenvalues)
    except Exception:
        eigenvalues = np.sort(np.linalg.eigvalsh(laplacian))[:k]

    # Find largest eigengap
    gaps = np.diff(eigenvalues)
    if len(gaps) == 0:
        return 1

    # Skip first eigenvalue (always ~0)
    n_clusters = np.argmax(gaps[1:]) + 2 if len(gaps) > 1 else 1
    return max(1, min(n_clusters, max_domains))


def threshold_assignment(co_membership: np.ndarray, n_clusters: int) -> np.ndarray:
    """Fallback domain assignment via thresholding."""
    L = co_membership.shape[0]
    labels = np.zeros(L, dtype=int)

    # Use hierarchical thresholding
    threshold = 0.5
    binary = (co_membership > threshold).astype(int)
    components, n_found = connected_components(binary)

    if n_found >= n_clusters:
        labels = components - 1  # 0-indexed
    else:
        # Simple uniform split
        chunk_size = max(1, L // n_clusters)
        for i in range(L):
            labels[i] = min(i // chunk_size, n_clusters - 1)

    return labels


def merge_small_domains(
    labels: np.ndarray,
    co_membership: np.ndarray,
    min_size: int,
) -> np.ndarray:
    """Merge domains smaller than min_size into their most similar neighbor."""
    unique_labels = np.unique(labels)
    label_sizes = {l: (labels == l).sum() for l in unique_labels}

    labels = labels.copy()
    changed = True
    while changed:
        changed = False
        for lbl in list(label_sizes.keys()):
            if label_sizes.get(lbl, 0) < min_size and len(label_sizes) > 1:
                mask = labels == lbl
                if not mask.any():
                    continue

                # Find most similar other domain
                best_target = None
                best_sim = -1
                for other_lbl in label_sizes:
                    if other_lbl == lbl:
                        continue
                    other_mask = labels == other_lbl
                    sim = co_membership[np.ix_(mask, other_mask)].mean()
                    if sim > best_sim:
                        best_sim = sim
                        best_target = other_lbl

                if best_target is not None:
                    labels[mask] = best_target
                    label_sizes[best_target] += label_sizes[lbl]
                    del label_sizes[lbl]
                    changed = True

    # Re-index labels to be contiguous
    unique = np.unique(labels)
    mapping = {old: new for new, old in enumerate(unique)}
    return np.array([mapping[l] for l in labels])


def labels_to_segments(labels: np.ndarray) -> list[list[tuple[int, int]]]:
    """Convert per-residue labels to list of domain segments.

    Returns:
        List of domains, each domain is a list of (start, end) segments (inclusive).
        Handles discontinuous domains naturally.
    """
    n_domains = labels.max() + 1
    domains = []

    for d in range(n_domains):
        mask = labels == d
        if not mask.any():
            continue

        # Find contiguous segments
        segments = []
        in_segment = False
        start = 0

        for i in range(len(mask)):
            if mask[i] and not in_segment:
                start = i
                in_segment = True
            elif not mask[i] and in_segment:
                segments.append((start, i - 1))
                in_segment = False

        if in_segment:
            segments.append((start, len(mask) - 1))

        if segments:
            domains.append(segments)

    return domains


def deduplicate_partitionings(
    partitionings: list[dict],
    similarity_threshold: float = 0.9,
) -> list[dict]:
    """Remove near-duplicate partitionings.

    Two partitionings are considered duplicates if their domain assignments
    overlap by more than similarity_threshold.
    """
    if len(partitionings) <= 1:
        return partitionings

    kept = [partitionings[0]]

    for part in partitionings[1:]:
        is_duplicate = False
        for existing in kept:
            # Compare assignments
            a1 = part["assignment"]
            a2 = existing["assignment"]
            if len(a1) != len(a2):
                continue

            # Compute overlap (accounting for label permutation)
            overlap = compute_assignment_overlap(a1, a2)
            if overlap > similarity_threshold:
                is_duplicate = True
                break

        if not is_duplicate:
            kept.append(part)

    return kept


def compute_assignment_overlap(a1: np.ndarray, a2: np.ndarray) -> float:
    """Compute best overlap between two assignments (accounts for label permutation)."""
    from scipy.optimize import linear_sum_assignment

    labels1 = np.unique(a1)
    labels2 = np.unique(a2)

    n1, n2 = len(labels1), len(labels2)
    n = max(n1, n2)

    # Build overlap matrix
    overlap = np.zeros((n, n))
    for i, l1 in enumerate(labels1):
        for j, l2 in enumerate(labels2):
            overlap[i, j] = ((a1 == l1) & (a2 == l2)).sum()

    # Hungarian matching to maximize overlap
    row_ind, col_ind = linear_sum_assignment(-overlap)
    total_overlap = overlap[row_ind, col_ind].sum()

    return total_overlap / len(a1)


@torch.no_grad()
def predict_partitionings(
    model_outputs: dict[str, torch.Tensor],
    seq_len: int,
    min_confidence: float = 0.1,
    min_domain_size: int = 20,
    max_domains: int = 20,
    similarity_threshold: float = 0.9,
) -> list[dict]:
    """Full post-processing pipeline: model outputs -> domain partitionings.

    Args:
        model_outputs: Model forward pass outputs for a single sample.
        seq_len: Actual sequence length (without padding).
        min_confidence: Minimum confidence to keep a partitioning.
        min_domain_size: Minimum residues per domain.
        max_domains: Maximum domains per partitioning.
        similarity_threshold: Threshold for deduplication.

    Returns:
        List of partitioning dicts, sorted by confidence (descending).
    """
    co_mem = torch.sigmoid(model_outputs["co_membership"][:, :seq_len, :seq_len])
    confidences = model_outputs["confidence"]
    K = co_mem.shape[0]

    # Get number of domains hints if available
    num_domains_hints = None
    if "num_domains_logits" in model_outputs:
        num_domains_hints = (
            model_outputs["num_domains_logits"].argmax(dim=-1) + 1
        )  # (K,)

    partitionings = []
    for k in range(K):
        conf = confidences[k].item()
        if conf < min_confidence:
            continue

        hint = num_domains_hints[k].item() if num_domains_hints is not None else None

        part = extract_partitioning(
            co_membership=co_mem[k].cpu().numpy(),
            confidence=conf,
            num_domains_hint=hint,
            min_domain_size=min_domain_size,
            max_domains=max_domains,
        )
        partitionings.append(part)

    # Sort by confidence
    partitionings.sort(key=lambda x: x["confidence"], reverse=True)

    # Deduplicate
    partitionings = deduplicate_partitionings(partitionings, similarity_threshold)

    return partitionings
