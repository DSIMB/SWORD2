"""
Training loop for DomainPartitionNet.

Features:
- Mixed precision training (fp16)
- Gradient accumulation
- Cosine LR schedule with linear warmup
- W&B logging
- Periodic evaluation and checkpointing
"""

import logging
import math
import os
from pathlib import Path

import torch
import torch.nn as nn
from torch.cuda.amp import GradScaler, autocast

from .config import Config
from .dataset import create_dataloaders
from .losses import PartitioningLoss
from .metrics import compute_metrics
from .model import DomainPartitionNet

logger = logging.getLogger(__name__)


def get_cosine_schedule_with_warmup(
    optimizer: torch.optim.Optimizer,
    warmup_steps: int,
    total_steps: int,
    min_lr_ratio: float = 0.01,
):
    """Cosine annealing with linear warmup."""

    def lr_lambda(step):
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        progress = (step - warmup_steps) / max(1, total_steps - warmup_steps)
        return min_lr_ratio + (1 - min_lr_ratio) * 0.5 * (1 + math.cos(math.pi * progress))

    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)


@torch.no_grad()
def evaluate(
    model: DomainPartitionNet,
    dataloader,
    criterion: PartitioningLoss,
    device: torch.device,
    max_batches: int = 50,
) -> dict[str, float]:
    """Run evaluation on validation set."""
    model.eval()
    total_losses = {}
    all_metrics = []
    n_batches = 0

    for batch in dataloader:
        if n_batches >= max_batches:
            break

        embeddings = batch["embeddings"].to(device)
        mask = batch["mask"].to(device)
        targets = batch["targets"]

        with autocast(dtype=torch.float16):
            outputs = model(embeddings, mask)
            losses = criterion(outputs, targets)

        for k, v in losses.items():
            total_losses[k] = total_losses.get(k, 0.0) + v.item()

        batch_metrics = compute_metrics(outputs, targets)
        all_metrics.append(batch_metrics)
        n_batches += 1

    avg_losses = {k: v / max(1, n_batches) for k, v in total_losses.items()}

    if all_metrics:
        avg_metrics = {}
        for key in all_metrics[0]:
            avg_metrics[key] = sum(m[key] for m in all_metrics) / len(all_metrics)
        avg_losses.update(avg_metrics)

    model.train()
    return avg_losses


def train(config: Config) -> None:
    """Main training function."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Training on {device}")

    # Create model
    model = DomainPartitionNet(config.model)
    model = model.to(device)

    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    logger.info(f"Parameters: {total_params:,} total, {trainable_params:,} trainable")

    # Create dataloaders
    embedding_sources = config.model.get_embedding_sources()
    loaders = create_dataloaders(
        data_dir=config.data.processed_data_dir,
        embedding_sources=embedding_sources,
        max_seq_len=config.data.max_seq_len,
        min_seq_len=config.data.min_seq_len,
        batch_size=config.data.batch_size,
        num_workers=config.data.num_workers,
        pin_memory=config.data.pin_memory,
    )

    if "train" not in loaders:
        raise RuntimeError("No training data found!")

    train_loader = loaders["train"]
    val_loader = loaders.get("val")

    # Loss, optimizer, scheduler
    criterion = PartitioningLoss(config.train)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.train.learning_rate,
        weight_decay=config.train.weight_decay,
    )
    scheduler = get_cosine_schedule_with_warmup(
        optimizer, config.train.warmup_steps, config.train.max_steps
    )

    scaler = GradScaler() if config.train.fp16 else None

    # W&B
    use_wandb = False
    try:
        import wandb

        wandb.init(
            project=config.train.wandb_project,
            name=config.train.wandb_run_name,
            config={
                "model": config.model.__dict__,
                "data": config.data.__dict__,
                "train": config.train.__dict__,
            },
        )
        use_wandb = True
    except ImportError:
        logger.info("wandb not available, logging to console only")

    # Resume from checkpoint
    global_step = 0
    if config.train.resume_from:
        checkpoint = torch.load(config.train.resume_from, map_location=device)
        model.load_state_dict(checkpoint["model"])
        optimizer.load_state_dict(checkpoint["optimizer"])
        scheduler.load_state_dict(checkpoint["scheduler"])
        global_step = checkpoint["global_step"]
        if scaler and "scaler" in checkpoint:
            scaler.load_state_dict(checkpoint["scaler"])
        logger.info(f"Resumed from step {global_step}")

    # Training loop
    os.makedirs(config.train.output_dir, exist_ok=True)
    model.train()
    optimizer.zero_grad()

    epoch = 0
    running_loss = 0.0
    best_val_loss = float("inf")

    while global_step < config.train.max_steps:
        epoch += 1

        for batch in train_loader:
            if global_step >= config.train.max_steps:
                break

            embeddings = batch["embeddings"].to(device)
            mask = batch["mask"].to(device)
            targets = batch["targets"]

            # Forward
            if config.train.fp16:
                with autocast(dtype=torch.float16):
                    outputs = model(embeddings, mask)
                    losses = criterion(outputs, targets)
                    loss = losses["loss"] / config.train.gradient_accumulation

                scaler.scale(loss).backward()
            else:
                outputs = model(embeddings, mask)
                losses = criterion(outputs, targets)
                loss = losses["loss"] / config.train.gradient_accumulation
                loss.backward()

            running_loss += losses["loss"].item()

            # Gradient accumulation step
            if (global_step + 1) % config.train.gradient_accumulation == 0:
                if config.train.fp16:
                    scaler.unscale_(optimizer)

                nn.utils.clip_grad_norm_(
                    model.parameters(), config.train.gradient_clip
                )

                if config.train.fp16:
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    optimizer.step()

                scheduler.step()
                optimizer.zero_grad()

            global_step += 1

            # Logging
            if global_step % config.train.log_every == 0:
                avg_loss = running_loss / config.train.log_every
                lr = scheduler.get_last_lr()[0]
                log_dict = {
                    "train/loss": avg_loss,
                    "train/lr": lr,
                    "train/step": global_step,
                    "train/epoch": epoch,
                }
                for k, v in losses.items():
                    if k != "loss":
                        log_dict[f"train/{k}"] = v.item()

                logger.info(
                    f"Step {global_step} | loss={avg_loss:.4f} | lr={lr:.2e}"
                )
                if use_wandb:
                    wandb.log(log_dict, step=global_step)
                running_loss = 0.0

            # Evaluation
            if val_loader and global_step % config.train.eval_every == 0:
                val_metrics = evaluate(model, val_loader, criterion, device)
                logger.info(
                    f"Step {global_step} | val_loss={val_metrics['loss']:.4f}"
                )
                if use_wandb:
                    wandb.log(
                        {f"val/{k}": v for k, v in val_metrics.items()},
                        step=global_step,
                    )

                if val_metrics["loss"] < best_val_loss:
                    best_val_loss = val_metrics["loss"]
                    save_checkpoint(
                        model, optimizer, scheduler, scaler,
                        global_step, config.train.output_dir, "best",
                    )

            # Checkpointing
            if global_step % config.train.save_every == 0:
                save_checkpoint(
                    model, optimizer, scheduler, scaler,
                    global_step, config.train.output_dir, f"step_{global_step}",
                )

    # Final save
    save_checkpoint(
        model, optimizer, scheduler, scaler, global_step, config.train.output_dir, "final"
    )
    logger.info(f"Training complete at step {global_step}")

    if use_wandb:
        wandb.finish()


def save_checkpoint(model, optimizer, scheduler, scaler, global_step, output_dir, name):
    """Save training checkpoint."""
    path = os.path.join(output_dir, f"checkpoint_{name}.pt")
    state = {
        "model": model.state_dict(),
        "optimizer": optimizer.state_dict(),
        "scheduler": scheduler.state_dict(),
        "global_step": global_step,
    }
    if scaler:
        state["scaler"] = scaler.state_dict()
    torch.save(state, path)
    logger.info(f"Saved checkpoint to {path}")


def main():
    """Entry point for training."""
    import argparse

    parser = argparse.ArgumentParser(description="Train DomainPartitionNet")
    parser.add_argument("--config", type=str, default=None, help="Path to config YAML file")
    parser.add_argument("--resume", type=str, default=None, help="Checkpoint to resume from")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    if args.config:
        config = Config.from_yaml(args.config)
    else:
        config = Config()

    if args.resume:
        config.train.resume_from = args.resume

    train(config)


if __name__ == "__main__":
    main()
