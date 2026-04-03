"""
Step 2: Train Sparse Autoencoders on scGPT activations.

Trains SAEs at configurable layers, dictionary sizes, and sparsity levels.

Usage:
    python -m src.train_saes \
        --activations-dir outputs/activations \
        --output-dir outputs/sae_models \
        --layers 0 6 11 \
        --dict-sizes 2048 \
        --lambdas 1.0 3.0 10.0
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

from src.model import SparseAutoencoder


def train_sae(
    activations: np.ndarray,
    d_hidden: int,
    lam: float,
    epochs: int = 40,
    lr: float = 3e-4,
    batch_size: int = 512,
) -> SparseAutoencoder:
    """Train a single SAE on activation data."""
    dataset = TensorDataset(torch.tensor(activations, dtype=torch.float32))
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)

    sae = SparseAutoencoder(activations.shape[1], d_hidden)
    optimizer = torch.optim.Adam(sae.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    for epoch in range(epochs):
        for (batch,) in loader:
            loss, mse, l1, _ = sae.loss(batch, lam)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(sae.parameters(), 1.0)
            optimizer.step()
            sae.normalize_decoder()
        scheduler.step()

        if epoch % 10 == 0:
            sae.eval()
            with torch.no_grad():
                sample = torch.tensor(
                    activations[: min(1000, len(activations))], dtype=torch.float32
                )
                _, f = sae(sample)
                avg_l0 = (f > 0).float().sum(dim=1).mean().item()
                n_alive = ((f > 0).float().mean(dim=0) > 0.01).sum().item()
            sae.train()
            print(f"    Epoch {epoch}: L0={avg_l0:.0f}, alive={n_alive}/{d_hidden}")

    return sae


def main(args):
    acts_dir = Path(args.activations_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    configs = [
        (layer, d, lam)
        for layer in args.layers
        for d in args.dict_sizes
        for lam in args.lambdas
    ]

    results = {}
    for layer_idx, d_hidden, lam in configs:
        name = f"layer{layer_idx}_d{d_hidden}_lam{lam}"
        model_path = out_dir / f"{name}.pt"
        if model_path.exists() and not args.overwrite:
            print(f"  {name}: exists, skipping (use --overwrite to retrain)")
            continue

        token_path = acts_dir / f"layer_{layer_idx}_token_acts.npy"
        if not token_path.exists():
            print(f"  {name}: no activations at {token_path}, skipping")
            continue

        acts = np.load(token_path)
        acts_mean = acts.mean(axis=0, keepdims=True)
        acts_std = acts.std(axis=0, keepdims=True)
        acts_std[acts_std < 1e-8] = 1.0
        acts_normed = (acts - acts_mean) / acts_std

        print(f"\nTraining {name} on {len(acts_normed)} tokens...")
        sae = train_sae(acts_normed, d_hidden, lam, epochs=args.epochs)

        torch.save(
            {
                "state_dict": sae.state_dict(),
                "d_input": acts.shape[1],
                "d_hidden": d_hidden,
                "lambda": lam,
                "layer": layer_idx,
                "acts_mean": acts_mean,
                "acts_std": acts_std,
            },
            model_path,
        )

        # Compute cell-level features
        cell_path = acts_dir / f"layer_{layer_idx}_cell_acts.npy"
        if cell_path.exists():
            cell_acts = np.load(cell_path)
            cell_normed = (cell_acts - acts_mean) / acts_std
            sae.eval()
            with torch.no_grad():
                cell_features = sae.encode(
                    torch.tensor(cell_normed, dtype=torch.float32)
                ).numpy()
            np.save(out_dir / f"{name}_cell_features.npy", cell_features)

            active = cell_features > 0
            avg_l0 = active.sum(axis=1).mean()
            alive = (active.mean(axis=0) > 0.01).sum()
            results[name] = {
                "layer": layer_idx,
                "d_hidden": d_hidden,
                "lambda": lam,
                "avg_L0": float(avg_l0),
                "alive": int(alive),
                "dead": int(d_hidden - alive),
            }
            print(f"  Cell-level: L0={avg_l0:.0f}, alive={alive}/{d_hidden}")

    with open(out_dir / "training_summary.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nDone. Models saved to {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train Sparse Autoencoders")
    parser.add_argument("--activations-dir", default="outputs/activations")
    parser.add_argument("--output-dir", default="outputs/sae_models")
    parser.add_argument("--layers", nargs="+", type=int, default=[0, 6, 11])
    parser.add_argument("--dict-sizes", nargs="+", type=int, default=[2048])
    parser.add_argument("--lambdas", nargs="+", type=float, default=[1.0, 3.0, 10.0])
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--overwrite", action="store_true")
    main(parser.parse_args())
