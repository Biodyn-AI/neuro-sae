"""
Step 2: Train Sparse Autoencoders on scGPT hidden state activations.
Anthropic-style ReLU SAE with L1 sparsity.
"""
import json
import numpy as np
import torch
import torch.nn as nn
from pathlib import Path
from torch.utils.data import DataLoader, TensorDataset

ACT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\activations")
OUT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\sae_models")
OUT_DIR.mkdir(parents=True, exist_ok=True)

D_MODEL = 512
DICT_SIZES = [2048, 4096]  # 4x and 8x
LAMBDAS = [0.01, 0.03, 0.1]
LAYERS_TO_TRAIN = [0, 6, 11]  # early, middle, late
EPOCHS = 30
BATCH_SIZE = 512
MAX_TRAIN_SAMPLES = 50000  # subsample tokens for speed on CPU
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
        # x - b_dec, then encode
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
        for (batch,) in loader:
            batch = batch.to(DEVICE)
            loss, mse, l1, _ = sae.loss(batch, lam)
            optimizer.zero_grad()
            loss.backward()
            # Clip gradients
            torch.nn.utils.clip_grad_norm_(sae.parameters(), 1.0)
            optimizer.step()
            
            # Normalize decoder weights to unit norm
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
            # Compute sparsity stats
            sae.eval()
            with torch.no_grad():
                sample = torch.tensor(activations[:min(1000, len(activations))], dtype=torch.float32)
                _, f = sae(sample)
                alive = (f > 0).float().mean(dim=0)
                n_alive = (alive > 0.01).sum().item()
                avg_l0 = (f > 0).float().sum(dim=1).mean().item()
            sae.train()
            print(f"  Epoch {epoch:3d}: loss={avg_loss:.6f} mse={avg_mse:.6f} l1={avg_l1:.6f} "
                  f"alive={n_alive}/{d_hidden} avg_L0={avg_l0:.1f}")
    
    return sae, history


def main():
    results = {}
    
    for layer_idx in LAYERS_TO_TRAIN:
        # Prefer token-level activations for more training data
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
        
        # Subsample for CPU training speed
        if len(acts) > MAX_TRAIN_SAMPLES:
            rng = np.random.RandomState(42)
            idx = rng.choice(len(acts), MAX_TRAIN_SAMPLES, replace=False)
            acts = acts[idx]
            print(f"  Subsampled to {len(acts)} tokens")
        
        # Normalize activations (zero mean, unit variance per dimension)
        acts_mean = acts.mean(axis=0, keepdims=True)
        acts_std = acts.std(axis=0, keepdims=True)
        acts_std[acts_std < 1e-8] = 1.0
        acts_normed = (acts - acts_mean) / acts_std
        
        # Save normalization stats
        np.save(OUT_DIR / f"layer_{layer_idx:02d}_mean.npy", acts_mean)
        np.save(OUT_DIR / f"layer_{layer_idx:02d}_std.npy", acts_std)
        
        for d_hidden in DICT_SIZES:
            for lam in LAMBDAS:
                name = f"layer{layer_idx:02d}_d{d_hidden}_lam{lam}"
                print(f"\n--- Training {name} ---")
                
                sae, history = train_sae(acts_normed, d_hidden, lam)
                
                # Save model
                torch.save({
                    "state_dict": sae.state_dict(),
                    "d_input": D_MODEL,
                    "d_hidden": d_hidden,
                    "lambda": lam,
                    "layer": layer_idx,
                    "history": history,
                }, OUT_DIR / f"{name}.pt")
                
                # Compute and save feature activations on all data
                sae.eval()
                with torch.no_grad():
                    all_acts_tensor = torch.tensor(acts_normed, dtype=torch.float32)
                    # Process in batches
                    features_list = []
                    for start in range(0, len(acts_normed), 512):
                        batch = all_acts_tensor[start:start+512]
                        _, f = sae(batch)
                        features_list.append(f.numpy())
                    features = np.concatenate(features_list, axis=0)
                
                np.save(OUT_DIR / f"{name}_features.npy", features)
                
                # Summary stats
                alive_frac = ((features > 0).mean(axis=0) > 0.01).mean()
                avg_l0 = (features > 0).sum(axis=1).mean()
                final_mse = history[-1]["mse"]
                
                results[name] = {
                    "layer": layer_idx,
                    "d_hidden": d_hidden,
                    "lambda": lam,
                    "final_mse": final_mse,
                    "alive_fraction": float(alive_frac),
                    "avg_L0": float(avg_l0),
                    "n_training_samples": len(acts_normed),
                }
                print(f"  Final: MSE={final_mse:.6f}, alive={alive_frac:.2%}, L0={avg_l0:.1f}")
    
    # Save results summary
    with open(OUT_DIR / "training_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("\n=== Step 2 Complete ===")
    print(f"Trained {len(results)} SAEs, saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
