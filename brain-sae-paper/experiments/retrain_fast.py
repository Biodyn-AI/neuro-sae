"""
Fast targeted SAE retraining — focus on best configs only.
Train d2048 on layers 0, 6, 11 with high λ to achieve L0 ≈ 10-50.
Uses binary search on λ to find good sparsity.
"""
import sys
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

D_MODEL = 512
EPOCHS = 40
BATCH_SIZE = 512
MAX_TRAIN_SAMPLES = 50000
LR = 3e-4
DEVICE = "cpu"

# Flush print
def log(msg):
    print(msg, flush=True)


class SparseAutoencoder(nn.Module):
    def __init__(self, d_input, d_hidden):
        super().__init__()
        self.d_input = d_input
        self.d_hidden = d_hidden
        self.W_enc = nn.Linear(d_input, d_hidden, bias=True)
        self.W_dec = nn.Linear(d_hidden, d_input, bias=True)
        with torch.no_grad():
            self.W_dec.weight.data = nn.functional.normalize(self.W_dec.weight.data, dim=0)

    def encode(self, x):
        return torch.relu(self.W_enc(x - self.W_dec.bias))

    def decode(self, f):
        return self.W_dec(f)

    def forward(self, x):
        f = self.encode(x)
        return self.decode(f), f

    def loss(self, x, lam):
        x_hat, f = self.forward(x)
        mse = (x - x_hat).pow(2).mean()
        l1 = f.abs().mean()
        return mse + lam * l1, mse, l1, f


def train_sae(activations, d_hidden, lam, epochs=EPOCHS):
    dataset = TensorDataset(torch.tensor(activations, dtype=torch.float32))
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)

    sae = SparseAutoencoder(D_MODEL, d_hidden).to(DEVICE)
    optimizer = torch.optim.Adam(sae.parameters(), lr=LR)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    for epoch in range(epochs):
        sae.train()
        total_mse, n = 0, 0
        for (batch,) in loader:
            loss, mse, l1, _ = sae.loss(batch, lam)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(sae.parameters(), 1.0)
            optimizer.step()
            with torch.no_grad():
                sae.W_dec.weight.data = nn.functional.normalize(sae.W_dec.weight.data, dim=0)
            total_mse += mse.item() * batch.shape[0]
            n += batch.shape[0]
        scheduler.step()

        if epoch % 10 == 0 or epoch == epochs - 1:
            sae.eval()
            with torch.no_grad():
                sample = torch.tensor(activations[:2000], dtype=torch.float32)
                _, f = sae(sample)
                avg_l0 = (f > 0).float().sum(dim=1).mean().item()
                alive = ((f > 0).float().mean(dim=0) > 0.01).sum().item()
            log(f"    epoch {epoch:3d}: mse={total_mse/n:.6f} L0={avg_l0:.1f} alive={alive}/{d_hidden}")
    return sae


def get_stats(sae, data):
    sae.eval()
    with torch.no_grad():
        feats_list = []
        ss_res, ss_tot = 0.0, 0.0
        data_mean = data.mean(axis=0)
        for start in range(0, len(data), 1024):
            batch = torch.tensor(data[start:start+1024], dtype=torch.float32)
            x_hat, f = sae(batch)
            feats_list.append(f.numpy())
            ss_res += (batch - x_hat).pow(2).sum().item()
            ss_tot += (batch - torch.tensor(data_mean)).pow(2).sum().item()
        features = np.concatenate(feats_list)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        alive_frac = ((features > 0).mean(axis=0) > 0.01).mean()
        avg_l0 = (features > 0).sum(axis=1).mean()
        dead = int((1 - alive_frac) * sae.d_hidden)
    return features, float(r2), float(alive_frac), float(avg_l0), dead


def main():
    results = {}
    t0 = time.time()
    
    configs = [
        # (layer, d_hidden, lambdas_to_try)
        (0, 2048, [1.0, 3.0, 10.0]),
        (6, 2048, [1.0, 3.0, 10.0]),
        (11, 2048, [1.0, 3.0, 10.0]),
    ]

    for layer_idx, d_hidden, lambdas in configs:
        token_path = ACT_DIR / f"layer_{layer_idx:02d}_token_activations.npy"
        cell_path = ACT_DIR / f"layer_{layer_idx:02d}_cell_activations.npy"

        log(f"\n=== Layer {layer_idx} ===")
        log(f"Loading activations...")
        
        if token_path.exists():
            acts = np.load(token_path)
            log(f"Loaded token activations {acts.shape}")
        else:
            log(f"No token activations, skipping")
            continue

        if len(acts) > MAX_TRAIN_SAMPLES:
            rng = np.random.RandomState(42)
            idx = rng.choice(len(acts), MAX_TRAIN_SAMPLES, replace=False)
            acts = acts[idx]
            log(f"Subsampled to {len(acts)}")

        acts_mean = acts.mean(axis=0, keepdims=True)
        acts_std = acts.std(axis=0, keepdims=True)
        acts_std[acts_std < 1e-8] = 1.0
        acts_normed = (acts - acts_mean) / acts_std
        
        np.save(OUT_DIR / f"layer_{layer_idx:02d}_mean.npy", acts_mean)
        np.save(OUT_DIR / f"layer_{layer_idx:02d}_std.npy", acts_std)

        for lam in lambdas:
            name = f"layer{layer_idx:02d}_d{d_hidden}_lam{lam}"
            elapsed = (time.time() - t0) / 60
            log(f"\n--- {name} ({elapsed:.1f} min elapsed) ---")
            
            sae = train_sae(acts_normed, d_hidden, lam)
            features, r2, alive_frac, avg_l0, dead = get_stats(sae, acts_normed)

            # Save model
            torch.save({
                "state_dict": sae.state_dict(),
                "d_input": D_MODEL,
                "d_hidden": d_hidden,
                "lambda": lam,
                "layer": layer_idx,
            }, OUT_DIR / f"{name}.pt")

            # Save features
            np.save(OUT_DIR / f"{name}_features.npy", features)

            # Cell-level features
            cell_acts = np.load(cell_path)
            cell_normed = (cell_acts - acts_mean) / acts_std
            sae.eval()
            with torch.no_grad():
                cf_list = []
                for s in range(0, len(cell_normed), 512):
                    b = torch.tensor(cell_normed[s:s+512], dtype=torch.float32)
                    _, f = sae(b)
                    cf_list.append(f.numpy())
                cell_features = np.concatenate(cf_list)
            np.save(OUT_DIR / f"{name}_cell_features.npy", cell_features)

            results[name] = {
                "layer": layer_idx, "d_hidden": d_hidden, "lambda": lam,
                "R2": r2, "alive_fraction": alive_frac,
                "avg_L0": avg_l0, "n_dead": dead,
            }
            log(f"  RESULT: R²={r2:.4f} L0={avg_l0:.1f} alive={alive_frac:.2%} dead={dead}")

    # Save
    with open(OUT_DIR / "training_results_sparse.json", "w") as f:
        json.dump(results, f, indent=2)
    
    elapsed = (time.time() - t0) / 60
    log(f"\n=== DONE ({elapsed:.1f} min) ===")
    for name, r in sorted(results.items()):
        flag = "★" if 5 <= r["avg_L0"] <= 100 else " "
        log(f"  {flag} {name}: R²={r['R2']:.4f} L0={r['avg_L0']:.1f} alive={r['alive_fraction']:.2%}")


if __name__ == "__main__":
    main()
