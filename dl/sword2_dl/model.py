"""
DomainPartitionNet: Lightweight model for protein domain partitioning.

Architecture:
1. Pre-computed PLM embeddings (one or more, concatenated) as input
2. Sinusoidal absolute position encoding + single representation transformer
3. Pair representation via outer sum + triangle updates + dilated 2D convolutions
4. Slot attention produces K diverse partitioning representations
5. Shared contact map head predicts residue-residue contacts (for visualization)
6. Multi-scale pair processing with downsampling/upsampling for large proteins
7. Auxiliary heads predict boundary probabilities and number of domains

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

import math
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

from .config import ModelConfig


class SinusoidalPositionEncoding(nn.Module):
    """Sinusoidal absolute position encoding for the single representation."""

    def __init__(self, dim: int, max_len: int = 4096):
        super().__init__()
        pe = torch.zeros(max_len, dim)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, dim, 2).float() * (-math.log(10000.0) / dim))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer("pe", pe.unsqueeze(0))  # (1, max_len, dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Add positional encoding. x: (B, L, D)"""
        return x + self.pe[:, :x.shape[1]]


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


class TriangleMultiplicationOutgoing(nn.Module):
    """Triangle multiplication (outgoing) from AlphaFold2.

    Updates pair (i,j) by aggregating over all intermediate positions k:
    pair(i,j) += sum_k gate(i,k) * proj(k,j)

    This propagates pairwise consistency through the pair representation.
    """

    def __init__(self, pair_dim: int, hidden_dim: int | None = None):
        super().__init__()
        hidden_dim = hidden_dim or pair_dim
        self.norm = nn.LayerNorm(pair_dim)
        self.proj_left = nn.Linear(pair_dim, hidden_dim)
        self.proj_right = nn.Linear(pair_dim, hidden_dim)
        self.gate_left = nn.Linear(pair_dim, hidden_dim)
        self.gate_right = nn.Linear(pair_dim, hidden_dim)
        self.out_proj = nn.Linear(hidden_dim, pair_dim)
        self.out_gate = nn.Linear(pair_dim, pair_dim)

    def forward(self, pair: torch.Tensor) -> torch.Tensor:
        """pair: (B, L, L, C) -> (B, L, L, C)"""
        z = self.norm(pair)
        left = self.proj_left(z) * torch.sigmoid(self.gate_left(z))  # (B, L, L, H)
        right = self.proj_right(z) * torch.sigmoid(self.gate_right(z))
        # Triangle: (B, i, k, H) @ (B, k, j, H)^T -> (B, i, j, H)
        out = torch.einsum("bikh,bjkh->bijh", left, right)
        out = self.out_proj(out) * torch.sigmoid(self.out_gate(z))
        return pair + out


class TriangleMultiplicationIncoming(nn.Module):
    """Triangle multiplication (incoming) from AlphaFold2.

    Updates pair (i,j) by aggregating: pair(i,j) += sum_k proj(k,i) * gate(k,j)
    """

    def __init__(self, pair_dim: int, hidden_dim: int | None = None):
        super().__init__()
        hidden_dim = hidden_dim or pair_dim
        self.norm = nn.LayerNorm(pair_dim)
        self.proj_left = nn.Linear(pair_dim, hidden_dim)
        self.proj_right = nn.Linear(pair_dim, hidden_dim)
        self.gate_left = nn.Linear(pair_dim, hidden_dim)
        self.gate_right = nn.Linear(pair_dim, hidden_dim)
        self.out_proj = nn.Linear(hidden_dim, pair_dim)
        self.out_gate = nn.Linear(pair_dim, pair_dim)

    def forward(self, pair: torch.Tensor) -> torch.Tensor:
        """pair: (B, L, L, C) -> (B, L, L, C)"""
        z = self.norm(pair)
        left = self.proj_left(z) * torch.sigmoid(self.gate_left(z))
        right = self.proj_right(z) * torch.sigmoid(self.gate_right(z))
        # Triangle: (B, k, i, H) @ (B, k, j, H)^T -> (B, i, j, H)
        out = torch.einsum("bkih,bkjh->bijh", left, right)
        out = self.out_proj(out) * torch.sigmoid(self.out_gate(z))
        return pair + out


class TriangleBlock(nn.Module):
    """A single triangle update block: outgoing + incoming + transition."""

    def __init__(self, pair_dim: int, dropout: float = 0.1):
        super().__init__()
        self.tri_out = TriangleMultiplicationOutgoing(pair_dim)
        self.tri_in = TriangleMultiplicationIncoming(pair_dim)
        self.transition = nn.Sequential(
            nn.LayerNorm(pair_dim),
            nn.Linear(pair_dim, pair_dim * 4),
            nn.GELU(),
            nn.Linear(pair_dim * 4, pair_dim),
            nn.Dropout(dropout),
        )

    def forward(self, pair: torch.Tensor) -> torch.Tensor:
        pair = self.tri_out(pair)
        pair = self.tri_in(pair)
        pair = pair + self.transition(pair)
        return pair


class SlotAttention(nn.Module):
    """Slot attention for diverse partitioning prediction (Locatello et al. 2020).

    K learnable slot vectors compete for pair features via iterative
    cross-attention, naturally producing diverse representations.
    """

    def __init__(
        self,
        pair_dim: int,
        slot_dim: int,
        num_slots: int,
        num_iters: int = 3,
        hidden_dim: int = 128,
    ):
        super().__init__()
        self.num_slots = num_slots
        self.num_iters = num_iters
        self.slot_dim = slot_dim

        # Learnable slot initializations
        self.slots_mu = nn.Parameter(torch.randn(1, num_slots, slot_dim) * (slot_dim ** -0.5))

        self.norm_input = nn.LayerNorm(pair_dim)
        self.norm_slots = nn.LayerNorm(slot_dim)

        # Project pair features to slot space
        self.to_k = nn.Linear(pair_dim, slot_dim)
        self.to_v = nn.Linear(pair_dim, slot_dim)
        self.to_q = nn.Linear(slot_dim, slot_dim)

        # GRU update for slots
        self.gru = nn.GRUCell(slot_dim, slot_dim)
        self.mlp = nn.Sequential(
            nn.LayerNorm(slot_dim),
            nn.Linear(slot_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, slot_dim),
        )

    def forward(self, pair: torch.Tensor) -> torch.Tensor:
        """
        Args:
            pair: (B, C, L, L) pair representation in channels-first format
        Returns:
            slots: (B, K, slot_dim) — K slot representations
        """
        B, C, L, _ = pair.shape

        # Flatten pair to (B, L*L, C), take upper triangle for efficiency
        pair_flat = pair.permute(0, 2, 3, 1).reshape(B, L * L, C)
        pair_flat = self.norm_input(pair_flat)

        k = self.to_k(pair_flat)  # (B, L*L, slot_dim)
        v = self.to_v(pair_flat)  # (B, L*L, slot_dim)

        # Initialize slots
        slots = self.slots_mu.expand(B, -1, -1).clone()  # (B, K, slot_dim)

        for _ in range(self.num_iters):
            slots_prev = slots
            slots = self.norm_slots(slots)

            q = self.to_q(slots)  # (B, K, slot_dim)

            # Attention: softmax over slots (competition)
            scale = self.slot_dim ** -0.5
            attn = torch.einsum("bkd,bnd->bkn", q, k) * scale  # (B, K, N)
            attn = attn.softmax(dim=1)  # normalize over slots (competition)

            # Weighted sum of values
            attn_norm = attn / (attn.sum(dim=-1, keepdim=True) + 1e-8)
            updates = torch.einsum("bkn,bnd->bkd", attn_norm, v)

            # GRU update
            slots = self.gru(
                updates.reshape(B * self.num_slots, self.slot_dim),
                slots_prev.reshape(B * self.num_slots, self.slot_dim),
            ).reshape(B, self.num_slots, self.slot_dim)

            slots = slots + self.mlp(slots)

        return slots


class MultiScalePairProcessor(nn.Module):
    """U-Net style multi-scale processing for pair representations.

    Processes the pair matrix at multiple resolutions for better
    receptive field coverage on large proteins.
    """

    def __init__(self, pair_dim: int, num_scales: int = 2, dropout: float = 0.1):
        super().__init__()
        self.num_scales = num_scales

        # Downsampling path
        self.down_convs = nn.ModuleList()
        self.down_blocks = nn.ModuleList()
        for _ in range(num_scales):
            self.down_convs.append(
                nn.Conv2d(pair_dim, pair_dim, 3, stride=2, padding=1)
            )
            self.down_blocks.append(nn.Sequential(
                nn.BatchNorm2d(pair_dim),
                nn.GELU(),
                nn.Conv2d(pair_dim, pair_dim, 3, padding=1),
                nn.BatchNorm2d(pair_dim),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Conv2d(pair_dim, pair_dim, 3, padding=1),
            ))

        # Upsampling path
        self.up_convs = nn.ModuleList()
        self.up_blocks = nn.ModuleList()
        for _ in range(num_scales):
            self.up_convs.append(
                nn.ConvTranspose2d(pair_dim, pair_dim, 2, stride=2)
            )
            self.up_blocks.append(nn.Sequential(
                nn.BatchNorm2d(pair_dim * 2),  # concat with skip connection
                nn.GELU(),
                nn.Conv2d(pair_dim * 2, pair_dim, 1),
                nn.BatchNorm2d(pair_dim),
                nn.GELU(),
                nn.Conv2d(pair_dim, pair_dim, 3, padding=1),
            ))

    def forward(self, pair: torch.Tensor) -> torch.Tensor:
        """pair: (B, C, L, L) -> (B, C, L, L)"""
        skips = []

        # Downsampling
        x = pair
        for down_conv, down_block in zip(self.down_convs, self.down_blocks):
            skips.append(x)
            x = down_conv(x)
            x = x + down_block(x)

        # Upsampling with skip connections
        for up_conv, up_block, skip in zip(
            self.up_convs, self.up_blocks, reversed(skips)
        ):
            x = up_conv(x)
            # Handle size mismatch from stride-2 downsampling of odd sizes
            if x.shape[-2:] != skip.shape[-2:]:
                x = F.pad(x, [0, skip.shape[-1] - x.shape[-1],
                              0, skip.shape[-2] - x.shape[-2]])
            x = up_block(torch.cat([x, skip], dim=1))

        return pair + x  # residual


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
    """Builds and refines pair representation from single embeddings.

    Pipeline: outer sum → relative pos → row/col attention → triangle updates
    → dilated conv blocks → multi-scale processing → symmetrize.
    """

    def __init__(self, config: ModelConfig):
        super().__init__()
        self.proj_left = nn.Linear(config.single_dim, config.pair_dim)
        self.proj_right = nn.Linear(config.single_dim, config.pair_dim)
        self.rel_pos = RelativePositionEncoding(config.pair_dim)
        self.input_norm = nn.LayerNorm(config.pair_dim)

        # Row and column attention for global context
        self.row_attn = nn.MultiheadAttention(
            config.pair_dim, num_heads=4, dropout=config.pair_dropout, batch_first=True
        )
        self.col_attn = nn.MultiheadAttention(
            config.pair_dim, num_heads=4, dropout=config.pair_dropout, batch_first=True
        )
        self.row_norm = nn.LayerNorm(config.pair_dim)
        self.col_norm = nn.LayerNorm(config.pair_dim)

        # Triangle update blocks (AlphaFold2-inspired)
        num_triangle = getattr(config, "num_triangle_blocks", 2)
        self.triangle_blocks = nn.ModuleList(
            [TriangleBlock(config.pair_dim, config.pair_dropout)
             for _ in range(num_triangle)]
        )

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

        # Multi-scale processing
        use_multiscale = getattr(config, "use_multiscale", True)
        if use_multiscale:
            num_scales = getattr(config, "num_scales", 2)
            self.multiscale = MultiScalePairProcessor(
                config.pair_dim, num_scales, config.pair_dropout
            )
        else:
            self.multiscale = None

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

        # Triangle updates (channels-last format: B, L, L, C)
        for tri_block in self.triangle_blocks:
            pair = tri_block(pair)

        # To channels-first for conv blocks: (B, L, L, C) -> (B, C, L, L)
        pair = pair.permute(0, 3, 1, 2).contiguous()

        # Dilated residual blocks
        for block in self.blocks:
            pair = block(pair)

        # Multi-scale processing (for large-protein receptive field)
        if self.multiscale is not None:
            pair = self.multiscale(pair)

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


class SlotPartitioningHead(nn.Module):
    """Slot-conditioned partitioning head.

    Takes a slot vector and the shared pair representation, produces
    a slot-specific co-membership matrix via FiLM conditioning.
    """

    def __init__(self, config: ModelConfig, slot_dim: int):
        super().__init__()
        self.predict_num_domains = config.predict_num_domains
        self.predict_boundaries = config.predict_boundaries

        # FiLM conditioning: slot → scale + bias for pair features
        self.film = nn.Linear(slot_dim, config.pair_dim * 2)

        # Per-slot refinement convolutions
        self.refine = nn.Sequential(
            nn.Conv2d(config.pair_dim, config.head_hidden_dim, 3, padding=1),
            nn.GELU(),
            nn.Conv2d(config.head_hidden_dim, config.head_hidden_dim, 3, padding=1),
            nn.GELU(),
        )

        # Co-membership prediction
        self.co_membership = nn.Conv2d(config.head_hidden_dim, 1, 1)

        # Confidence from slot vector directly
        self.confidence = nn.Sequential(
            nn.Linear(slot_dim, config.head_hidden_dim),
            nn.GELU(),
            nn.Linear(config.head_hidden_dim, 1),
            nn.Sigmoid(),
        )

        if self.predict_num_domains:
            self.num_domains_head = nn.Sequential(
                nn.Linear(slot_dim, config.head_hidden_dim),
                nn.GELU(),
                nn.Linear(config.head_hidden_dim, config.max_num_domains),
            )

        if self.predict_boundaries:
            self.boundary_head = nn.Sequential(
                nn.Linear(config.single_dim + config.pair_dim + slot_dim, config.head_hidden_dim),
                nn.GELU(),
                nn.Linear(config.head_hidden_dim, 1),
            )

    def forward(
        self,
        pair: torch.Tensor,
        single: torch.Tensor,
        slot: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        """
        Args:
            pair: (B, pair_dim, L, L) shared pair representation
            single: (B, L, single_dim)
            slot: (B, slot_dim) slot vector for this partitioning
        """
        out = {}

        # FiLM conditioning: modulate pair features with slot
        film_params = self.film(slot)  # (B, pair_dim * 2)
        gamma, beta = film_params.chunk(2, dim=-1)  # each (B, pair_dim)
        conditioned_pair = pair * gamma[:, :, None, None] + beta[:, :, None, None]

        # Refine
        refined = self.refine(conditioned_pair)

        # Co-membership
        co_mem = self.co_membership(refined).squeeze(1)
        co_mem = (co_mem + co_mem.transpose(-2, -1)) / 2
        out["co_membership"] = co_mem

        # Confidence from slot
        out["confidence"] = self.confidence(slot).squeeze(-1)

        if self.predict_num_domains:
            out["num_domains_logits"] = self.num_domains_head(slot)

        if self.predict_boundaries:
            B, C, L, _ = pair.shape
            diag = pair.diagonal(dim1=-2, dim2=-1).permute(0, 2, 1)  # (B, L, C)
            slot_expanded = slot.unsqueeze(1).expand(-1, L, -1)  # (B, L, slot_dim)
            boundary_input = torch.cat([single, diag, slot_expanded], dim=-1)
            out["boundary_logits"] = self.boundary_head(boundary_input).squeeze(-1)

        return out


class DomainPartitionNet(nn.Module):
    """Lightweight model for protein domain partitioning from pre-computed embeddings.

    Takes concatenated PLM embeddings as input.
    Uses slot attention for diverse partitioning prediction and triangle
    updates for pairwise consistency.
    """

    def __init__(self, config: ModelConfig):
        super().__init__()
        self.config = config
        self.use_slot_attention = getattr(config, "use_slot_attention", True)
        K = config.num_partitioning_slots

        total_embed_dim = config.total_embed_dim

        # Absolute position encoding
        self.pos_encoding = SinusoidalPositionEncoding(config.single_dim)

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

        # Pair module (now with triangle updates + multi-scale)
        self.pair_module = PairModule(config)

        # Slot attention + slot-conditioned heads
        slot_dim = getattr(config, "slot_dim", config.pair_dim)
        if self.use_slot_attention:
            self.slot_attention = SlotAttention(
                pair_dim=config.pair_dim,
                slot_dim=slot_dim,
                num_slots=K,
                num_iters=getattr(config, "slot_iters", 3),
                hidden_dim=config.head_hidden_dim,
            )
            self.slot_head = SlotPartitioningHead(config, slot_dim)
        else:
            # Fallback: K independent heads (original architecture)
            self.heads = nn.ModuleList(
                [PartitioningHead(config) for _ in range(K)]
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

        # 2. Add absolute position encoding
        single = self.pos_encoding(single)

        # 3. Refine single representation
        if mask is not None:
            src_key_padding_mask = ~mask  # True = padding
        else:
            src_key_padding_mask = None
        single = self.single_transformer(
            single, src_key_padding_mask=src_key_padding_mask
        )  # (B, L, single_dim)

        # 4. Build pair representation (with triangle updates + multi-scale)
        pair = self.pair_module(single, mask)  # (B, pair_dim, L, L)

        # 5. Generate K partitioning predictions
        if self.use_slot_attention:
            # Slot attention: K slots compete for pair features
            slots = self.slot_attention(pair)  # (B, K, slot_dim)

            all_co_mem = []
            all_conf = []
            all_num_dom = []
            all_bound = []

            for k in range(self.config.num_partitioning_slots):
                head_out = self.slot_head(pair, single, slots[:, k])
                all_co_mem.append(head_out["co_membership"])
                all_conf.append(head_out["confidence"])
                if "num_domains_logits" in head_out:
                    all_num_dom.append(head_out["num_domains_logits"])
                if "boundary_logits" in head_out:
                    all_bound.append(head_out["boundary_logits"])
        else:
            # Fallback: K independent heads
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

        # 6. Shared contact map prediction
        if hasattr(self, "contact_map_head"):
            outputs["contact_map_logits"] = self.contact_map_head(pair)  # (B, L, L)

        return outputs
