"""
Retrain SAEs with proper sparsity (target L0 ≈ 10-50).
The previous λ values (0.01-0.1) gave L0=1200-2400 — essentially not sparse at all.
We need λ ≈ 1-100 to achieve real sparsity.
"""
import json
import numpy as np
import torch
import torch.nn as nn
from pathlib import Path
from torch.utils.data import DataLoader, TensorDataset
import time

import platform
if platform.system() == "Windows":
    ACT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\activations")
    OUT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\sae_models")
else:
    ACT_DIR = Path("/mnt/d/openclaw/biodyn-nmi-paper/brain-sae-paper/experiments/activations")
    OUT_DIR = Path("/mnt/d/openclaw/biodyn-nmi-paper/brain-sae-paper/experiments/sae_models")
OUT_DIR.mkdir(parents=True, exist_ok=True)

D_MODEL = 512
DICT_SIZES = [2048, 4096]
# Much higher λ values for actual sparsity
LAMBDAS = [1.0, 3.0, 10.0, 30.0]
LAYERS_TO_TRAIN = [0, 6, 11]
EPOCHS = 50
BATCH_SIZE = 512
MAX_TRAIN_SAMPLES = 50000
LR = 3e-4
DEVICE = "cpu"


class SparseAutoencoder(nn.Module):
    """Anthropic-style ReLU SAE."""
    def __init__(self, d_input, d_hidden):
        super().__init__()
        self.d_input = d_input
        self.d_hidden = d_hidden

        self.W_enc = nn.Linear(d_input, d_hidden, bias=True)
        self.W_dec = nn.Linear(d_hidden, d_input, bias=True)

        # Initialize decoder columns to unit norm
        with torch.no_grad():
            self.W_dec.weight.data = nn.functional.normalize(self.W_dec.weight.data, dim=0)

    def encode(self, x):
        x_centered = x - self.W_dec.bias
        return torch.relu(self.W_enc(x_centered))

    def decode(self, f):
        return self.W_dec(f)

    def forward(self, x):
        f = self.encode(x)
        x_hat = self.decode(f)
        return x_hat, f

    def loss(self, x, lam):
        x_hat, f = self.forward(x)
        mse = (x - x_hat).pow(2).mean()
        l1 = f.abs().mean()
        return mse + lam * l1, mse, l1, f


def compute_r_squared(sae, data, batch_size=1024):
    """Compute R² reconstruction fidelity."""
    sae.eval()
    ss_res, ss_tot = 0.0, 0.0
    data_mean = data.mean(axis=0)
    with torch.no_grad():
        for start in range(0, len(data), batch_size):
            batch = torch.tensor(data[start:start+batch_size], dtype=torch.float32)
            x_hat, _ = sae(batch)
            ss_res += (batch - x_hat).pow(2).sum().item()
            ss_tot += (batch - torch.tensor(data_mean, dtype=torch.float32)).pow(2).sum().item()
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def train_sae(activations, d_hidden, lam, epochs=EPOCHS, lr=LR):
    """Train a single SAE."""
    dataset = TensorDataset(torch.tensor(activations, dtype=torch.float32))
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)

    sae = SparseAutoencoder(D_MODEL, d_hidden).to(DEVICE)
    optimizer = torch.optim.Adam(sae.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    history = []
    for epoch in range(epochs):
        total_loss, total_mse, total_l1, n = 0, 0, 0, 0
        sae.train()
        for (batch,) in loader:
            batch = batch.to(DEVICE)
            loss, mse, l1, _ = sae.loss(batch, lam)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(sae.parameters(), 1.0)
            optimizer.step()

            # Normalize decoder weights
            with torch.no_grad():
                sae.W_dec.weight.data = nn.functional.normalize(sae.W_dec.weight.data, dim=0)

            total_loss += loss.item() * batch.shape[0]
            total_mse += mse.item() * batch.shape[0]
            total_l1 += l1.item() * batch.shape[0]
            n += batch.shape[0]

        scheduler.step()
        avg_loss = total_loss / n
        avg_mse = total_mse / n
        avg_l1 = total_l1 / n
        history.append({"epoch": epoch, "loss": avg_loss, "mse": avg_mse, "l1": avg_l1})

        if epoch % 10 == 0 or epoch == epochs - 1:
            sae.eval()
            with torch.no_grad():
                sample = torch.tensor(activations[:min(2000, len(activations))], dtype=torch.float32)
                _, f = sae(sample)
                alive = (f > 0).float().mean(dim=0)
                n_alive = (alive > 0.01).sum().item()
                avg_l0 = (f > 0).float().sum(dim=1).mean().item()
                # Dead features (never activate on >1% of data)
                dead = d_hidden - n_alive
            sae.train()
            print(f"  Epoch {epoch:3d}: loss={avg_loss:.6f} mse={avg_mse:.6f} l1={avg_l1:.6f} "
                  f"alive={n_alive}/{d_hidden} dead={dead} avg_L0={avg_l0:.1f}")

    return sae, history


def main():
    results = {}
    start_time = time.time()

    for layer_idx in LAYERS_TO_TRAIN:
        token_path = ACT_DIR / f"layer_{layer_idx:02d}_token_activations.npy"
        cell_path = ACT_DIR / f"layer_{layer_idx:02d}_cell_activations.npy"

        if token_path.exists():
            acts = np.load(token_path)
            print(f"\nLayer {layer_idx}: loaded token activations {acts.shape}")
        elif cell_path.exists():
            acts = np.load(cell_path)
            print(f"\nLayer {layer_idx}: loaded cell activations {acts.shape}")
        else:
            print(f"\nLayer {layer_idx}: no activations found, skipping")
            continue

        if len(acts) > MAX_TRAIN_SAMPLES:
            rng = np.random.RandomState(42)
            idx = rng.choice(len(acts), MAX_TRAIN_SAMPLES, replace=False)
            acts = acts[idx]
            print(f"  Subsampled to {len(acts)} tokens")

        # Normalize
        acts_mean = acts.mean(axis=0, keepdims=True)
        acts_std = acts.std(axis=0, keepdims=True)
        acts_std[acts_std < 1e-8] = 1.0
        acts_normed = (acts - acts_mean) / acts_std

        # Save normalization stats (overwrite with potentially updated ones)
        np.save(OUT_DIR / f"layer_{layer_idx:02d}_mean.npy", acts_mean)
        np.save(OUT_DIR / f"layer_{layer_idx:02d}_std.npy", acts_std)

        for d_hidden in DICT_SIZES:
            for lam in LAMBDAS:
                name = f"layer{layer_idx:02d}_d{d_hidden}_lam{lam}"
                print(f"\n--- Training {name} ---")
                elapsed = time.time() - start_time
                print(f"  (elapsed: {elapsed/60:.1f} min)")

                sae, history = train_sae(acts_normed, d_hidden, lam)

                # Compute final stats
                sae.eval()
                with torch.no_grad():
                    all_acts_tensor = torch.tensor(acts_normed, dtype=torch.float32)
                    features_list = []
                    for start in range(0, len(acts_normed), 512):
                        batch = all_acts_tensor[start:start+512]
                        _, f = sae(batch)
                        features_list.append(f.numpy())
                    features = np.concatenate(features_list, axis=0)

                alive_frac = ((features > 0).mean(axis=0) > 0.01).mean()
                avg_l0 = (features > 0).sum(axis=1).mean()
                final_mse = history[-1]["mse"]
                r2 = compute_r_squared(sae, acts_normed)

                # Save model
                torch.save({
                    "state_dict": sae.state_dict(),
                    "d_input": D_MODEL,
                    "d_hidden": d_hidden,
                    "lambda": lam,
                    "layer": layer_idx,
                    "history": history,
                }, OUT_DIR / f"{name}.pt")

                # Save features
                np.save(OUT_DIR / f"{name}_features.npy", features)

                # Also compute and save cell-level features
                cell_acts = np.load(cell_path)
                cell_acts_normed = (cell_acts - acts_mean) / acts_std
                with torch.no_grad():
                    cell_feats_list = []
                    for start in range(0, len(cell_acts_normed), 512):
                        batch = torch.tensor(cell_acts_normed[start:start+512], dtype=torch.float32)
                        _, f = sae(batch)
                        cell_feats_list.append(f.numpy())
                    cell_features = np.concatenate(cell_feats_list, axis=0)
                np.save(OUT_DIR / f"{name}_cell_features.npy", cell_features)

                results[name] = {
                    "layer": layer_idx,
                    "d_hidden": d_hidden,
                    "lambda": lam,
                    "final_mse": final_mse,
                    "R2": float(r2),
                    "alive_fraction": float(alive_frac),
                    "avg_L0": float(avg_l0),
                    "n_dead_features": int((1 - alive_frac) * d_hidden),
                    "n_training_samples": len(acts_normed),
                }
                print(f"  Final: MSE={final_mse:.6f}, R²={r2:.4f}, alive={alive_frac:.2%}, L0={avg_l0:.1f}")

                # Early skip: if L0 is already in good range, don't need higher λ
                # (but train all for completeness)

    # Save results
    with open(OUT_DIR / "training_results_sparse.json", "w") as f:
        json.dump(results, f, indent=2)

    elapsed = time.time() - start_time
    print(f"\n=== Sparse SAE Training Complete ({elapsed/60:.1f} min) ===")
    print(f"\nSummary:")
    for name, r in sorted(results.items()):
        print(f"  {name}: R²={r['R2']:.4f} MSE={r['final_mse']:.6f} L0={r['avg_L0']:.1f} "
              f"alive={r['alive_fraction']:.2%} dead={r['n_dead_features']}")

    # Identify best models (R² > 0.8 and L0 in 10-50 range)
    print(f"\n--- Best models (L0 in [5, 100] range) ---")
    for name, r in sorted(results.items(), key=lambda x: abs(x[1]['avg_L0'] - 30)):
        if 5 <= r['avg_L0'] <= 100:
            print(f"  ★ {name}: R²={r['R2']:.4f} L0={r['avg_L0']:.1f}")


if __name__ == "__main__":
    main()
