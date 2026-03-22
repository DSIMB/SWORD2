"""
DomainPartitionNet: Lightweight model for protein domain partitioning.

Architecture (no PLM backbone — uses pre-computed frozen embeddings):
1. Pre-computed PLM embeddings (one or more, concatenated) as input
2. Single representation projection + transformer refinement
3. Pair representation via outer sum + dilated residual 2D convolutions
4. K partitioning heads predict co-membership matrices + confidence scores
5. Shared contact map head predicts residue-residue contacts (for visualization)
6. Auxiliary heads predict boundary probabilities and number of domains

The model outputs K alternative domain partitionings, each represented as
a symmetric L x L co-membership matrix where entry (i,j) indicates the
probability that residues i and j belong to the same domain, plus a shared
contact map prediction.

By consuming frozen embeddings instead of running a PLM backbone:
- All trainable parameters are in the lightweight head (~5-15M vs 650M+)
- Multiple PLM embeddings can be concatenated for richer features
- Training requires only 8-16 GB GPU memory
- Embedding computation is a one-time cost, amortized over all experiments
"""

from typing import Optional

import torch
import torch.nn as nn

from .config import ModelConfig


class RelativePositionEncoding(nn.Module):
    """Relative position encoding for the pair representation."""

    def __init__(self, pair_dim: int, max_relative_pos: int = 32):
        super().__init__()
        self.max_relative_pos = max_relative_pos
        num_embeddings = 2 * max_relative_pos + 1
        self.embedding = nn.Embedding(num_embeddings, pair_dim)

    def forward(self, seq_len: int, device: torch.device) -> torch.Tensor:
        pos = torch.arange(seq_len, device=device)
        rel_pos = pos.unsqueeze(0) - pos.unsqueeze(1)  # (L, L)
        rel_pos = rel_pos.clamp(-self.max_relative_pos, self.max_relative_pos)
        rel_pos = rel_pos + self.max_relative_pos  # shift to [0, 2*max]
        return self.embedding(rel_pos)  # (L, L, pair_dim)


class DilatedResBlock2D(nn.Module):
    """Dilated residual convolution block for pair representation."""

    def __init__(self, channels: int, kernel_size: int, dilation: int, dropout: float):
        super().__init__()
        padding = dilation * (kernel_size - 1) // 2
        self.net = nn.Sequential(
            nn.BatchNorm2d(channels),
            nn.GELU(),
            nn.Conv2d(channels, channels, kernel_size, padding=padding, dilation=dilation),
            nn.BatchNorm2d(channels),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv2d(channels, channels, kernel_size, padding=padding, dilation=dilation),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class PairModule(nn.Module):
    """Builds and refines pair representation from single embeddings."""

    def __init__(self, config: ModelConfig):
        super().__init__()
        self.proj_left = nn.Linear(config.single_dim, config.pair_dim)
        self.proj_right = nn.Linear(config.single_dim, config.pair_dim)
        self.rel_pos = RelativePositionEncoding(config.pair_dim)
        self.input_norm = nn.LayerNorm(config.pair_dim)

        # Dilated residual blocks
        blocks = []
        for i in range(config.pair_num_blocks):
            dilation = config.pair_dilation_cycle[i % len(config.pair_dilation_cycle)]
            blocks.append(
                DilatedResBlock2D(
                    config.pair_dim,
                    config.pair_kernel_size,
                    dilation,
                    config.pair_dropout,
                )
            )
        self.blocks = nn.ModuleList(blocks)

        # Row and column attention for global context
        self.row_attn = nn.MultiheadAttention(
            config.pair_dim, num_heads=4, dropout=config.pair_dropout, batch_first=True
        )
        self.col_attn = nn.MultiheadAttention(
            config.pair_dim, num_heads=4, dropout=config.pair_dropout, batch_first=True
        )
        self.row_norm = nn.LayerNorm(config.pair_dim)
        self.col_norm = nn.LayerNorm(config.pair_dim)

    def forward(
        self, single: torch.Tensor, mask: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        """
        Args:
            single: (B, L, single_dim)
            mask: (B, L) boolean, True = valid
        Returns:
            pair: (B, pair_dim, L, L)
        """
        B, L, _ = single.shape

        left = self.proj_left(single)   # (B, L, pair_dim)
        right = self.proj_right(single) # (B, L, pair_dim)

        # Outer sum
        pair = left.unsqueeze(2) + right.unsqueeze(1)  # (B, L, L, pair_dim)

        # Add relative position encoding
        pair = pair + self.rel_pos(L, single.device)

        pair = self.input_norm(pair)

        # Row attention: treat each row as a sequence
        pair_flat = pair.reshape(B * L, L, -1)
        pair_flat = pair_flat + self.row_attn(
            self.row_norm(pair_flat),
            self.row_norm(pair_flat),
            pair_flat,
        )[0]
        pair = pair_flat.reshape(B, L, L, -1)

        # Column attention: treat each column as a sequence
        pair_t = pair.transpose(1, 2).reshape(B * L, L, -1)
        pair_t = pair_t + self.col_attn(
            self.col_norm(pair_t),
            self.col_norm(pair_t),
            pair_t,
        )[0]
        pair = pair_t.reshape(B, L, L, -1).transpose(1, 2)

        # To channels-first for conv blocks: (B, L, L, C) -> (B, C, L, L)
        pair = pair.permute(0, 3, 1, 2).contiguous()

        # Dilated residual blocks
        for block in self.blocks:
            pair = block(pair)

        # Symmetrize
        pair = (pair + pair.transpose(-2, -1)) / 2

        return pair


class PartitioningHead(nn.Module):
    """Single partitioning prediction head.

    Outputs:
    - Co-membership matrix: L x L symmetric probability matrix
    - Confidence score: scalar in [0, 1]
    - Number of domains: classification over [1, max_domains]
    - Boundary probabilities: per-residue probability of being a domain boundary
    """

    def __init__(self, config: ModelConfig):
        super().__init__()
        self.predict_num_domains = config.predict_num_domains
        self.predict_boundaries = config.predict_boundaries

        # Co-membership prediction from pair representation
        self.co_membership = nn.Sequential(
            nn.Conv2d(config.pair_dim, config.head_hidden_dim, 1),
            nn.GELU(),
            nn.Conv2d(config.head_hidden_dim, 1, 1),
        )

        # Confidence from global pooling of pair representation
        self.confidence = nn.Sequential(
            nn.AdaptiveAvgPool2d(1),
            nn.Flatten(),
            nn.Linear(config.pair_dim, config.head_hidden_dim),
            nn.GELU(),
            nn.Linear(config.head_hidden_dim, 1),
            nn.Sigmoid(),
        )

        if self.predict_num_domains:
            self.num_domains_head = nn.Sequential(
                nn.AdaptiveAvgPool2d(1),
                nn.Flatten(),
                nn.Linear(config.pair_dim, config.head_hidden_dim),
                nn.GELU(),
                nn.Linear(config.head_hidden_dim, config.max_num_domains),
            )

        if self.predict_boundaries:
            # Boundary prediction from single + pair (diagonal) features
            self.boundary_head = nn.Sequential(
                nn.Linear(config.single_dim + config.pair_dim, config.head_hidden_dim),
                nn.GELU(),
                nn.Linear(config.head_hidden_dim, 1),
            )

    def forward(
        self,
        pair: torch.Tensor,
        single: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        """
        Args:
            pair: (B, pair_dim, L, L)
            single: (B, L, single_dim)
        Returns:
            dict with keys: co_membership, confidence, num_domains, boundaries
        """
        out = {}

        # Co-membership matrix
        co_mem = self.co_membership(pair).squeeze(1)  # (B, L, L)
        co_mem = (co_mem + co_mem.transpose(-2, -1)) / 2  # symmetrize
        out["co_membership"] = co_mem

        # Confidence
        out["confidence"] = self.confidence(pair).squeeze(-1)  # (B,)

        if self.predict_num_domains:
            out["num_domains_logits"] = self.num_domains_head(pair)  # (B, max_domains)

        if self.predict_boundaries:
            # Extract diagonal of pair representation
            B, C, L, _ = pair.shape
            diag = pair.diagonal(dim1=-2, dim2=-1)  # (B, C, L)
            diag = diag.permute(0, 2, 1)  # (B, L, C)
            boundary_input = torch.cat([single, diag], dim=-1)
            out["boundary_logits"] = self.boundary_head(boundary_input).squeeze(-1)  # (B, L)

        return out


class ContactMapHead(nn.Module):
    """Shared contact map prediction head.

    Predicts a symmetric L x L binary contact map from the pair representation.
    This is shared across all partitioning heads (contacts are a structural
    property independent of domain assignment).

    Useful for:
    - Visualization: users can view predicted contacts overlaid with domains
    - Auxiliary training signal: contact prediction regularizes the pair module
    - Structural validation: predicted contacts should be consistent with domains
    """

    def __init__(self, config: ModelConfig):
        super().__init__()
        self.contact_head = nn.Sequential(
            nn.Conv2d(config.pair_dim, config.head_hidden_dim, 1),
            nn.GELU(),
            nn.Conv2d(config.head_hidden_dim, config.head_hidden_dim, 3, padding=1),
            nn.GELU(),
            nn.Conv2d(config.head_hidden_dim, 1, 1),
        )

    def forward(self, pair: torch.Tensor) -> torch.Tensor:
        """
        Args:
            pair: (B, pair_dim, L, L)
        Returns:
            contact_logits: (B, L, L) symmetric contact map logits
        """
        contact = self.contact_head(pair).squeeze(1)  # (B, L, L)
        contact = (contact + contact.transpose(-2, -1)) / 2  # symmetrize
        return contact


class DomainPartitionNet(nn.Module):
    """Lightweight model for protein domain partitioning from pre-computed embeddings.

    Takes concatenated PLM embeddings as input (no backbone needed).
    All ~5-15M trainable parameters are in the projection, transformer,
    pair module, and partitioning heads.
    """

    def __init__(self, config: ModelConfig):
        super().__init__()
        self.config = config

        total_embed_dim = config.total_embed_dim

        # Project concatenated PLM embeddings to single representation
        self.single_proj = nn.Sequential(
            nn.LayerNorm(total_embed_dim),
            nn.Linear(total_embed_dim, config.single_dim),
            nn.GELU(),
            nn.Linear(config.single_dim, config.single_dim),
        )

        # Single representation refinement with transformer
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=config.single_dim,
            nhead=config.single_num_heads,
            dim_feedforward=config.single_dim * 4,
            dropout=config.single_dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.single_transformer = nn.TransformerEncoder(
            encoder_layer, num_layers=config.single_num_layers
        )

        # Pair module
        self.pair_module = PairModule(config)

        # K partitioning heads
        self.heads = nn.ModuleList(
            [PartitioningHead(config) for _ in range(config.num_partitioning_slots)]
        )

        # Shared contact map head
        if config.predict_contact_map:
            self.contact_map_head = ContactMapHead(config)

    def forward(
        self,
        embeddings: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
    ) -> dict[str, torch.Tensor]:
        """Forward pass.

        Args:
            embeddings: (B, L, D) concatenated pre-computed PLM embeddings
            mask: (B, L) boolean mask, True = valid residue

        Returns:
            Dictionary with:
                - co_membership: (B, K, L, L) co-membership logits
                - confidence: (B, K) confidence scores per slot
                - contact_map_logits: (B, L, L) contact map logits (if enabled)
                - num_domains_logits: (B, K, max_domains) if enabled
                - boundary_logits: (B, K, L) if enabled
        """
        # 1. Project to single representation
        single = self.single_proj(embeddings)  # (B, L, single_dim)

        # 2. Refine single representation
        if mask is not None:
            src_key_padding_mask = ~mask  # True = padding
        else:
            src_key_padding_mask = None
        single = self.single_transformer(
            single, src_key_padding_mask=src_key_padding_mask
        )  # (B, L, single_dim)

        # 3. Build pair representation
        pair = self.pair_module(single, mask)  # (B, pair_dim, L, L)

        # 4. Apply K partitioning heads
        all_co_mem = []
        all_conf = []
        all_num_dom = []
        all_bound = []

        for head in self.heads:
            head_out = head(pair, single)
            all_co_mem.append(head_out["co_membership"])
            all_conf.append(head_out["confidence"])
            if "num_domains_logits" in head_out:
                all_num_dom.append(head_out["num_domains_logits"])
            if "boundary_logits" in head_out:
                all_bound.append(head_out["boundary_logits"])

        outputs = {
            "co_membership": torch.stack(all_co_mem, dim=1),  # (B, K, L, L)
            "confidence": torch.stack(all_conf, dim=1),  # (B, K)
        }
        if all_num_dom:
            outputs["num_domains_logits"] = torch.stack(all_num_dom, dim=1)
        if all_bound:
            outputs["boundary_logits"] = torch.stack(all_bound, dim=1)

        # 5. Shared contact map prediction
        if hasattr(self, "contact_map_head"):
            outputs["contact_map_logits"] = self.contact_map_head(pair)  # (B, L, L)

        return outputs
