"""Finish training the last model (layer11 lam10.0) and save all results."""
import json, numpy as np, torch, torch.nn as nn, platform, time
from pathlib import Path
from torch.utils.data import DataLoader, TensorDataset

if platform.system() == "Windows":
    ACT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\activations")
    OUT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\sae_models")
else:
    ACT_DIR = Path("/mnt/d/openclaw/biodyn-nmi-paper/brain-sae-paper/experiments/activations")
    OUT_DIR = Path("/mnt/d/openclaw/biodyn-nmi-paper/brain-sae-paper/experiments/sae_models")

D_MODEL = 512

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
        return (x - x_hat).pow(2).mean() + lam * f.abs().mean(), (x - x_hat).pow(2).mean(), f

# Train layer11 lam10.0
print("Loading layer 11 activations...", flush=True)
acts = np.load(ACT_DIR / "layer_11_token_activations.npy")
rng = np.random.RandomState(42)
idx = rng.choice(len(acts), 50000, replace=False)
acts = acts[idx]
acts_mean = acts.mean(axis=0, keepdims=True)
acts_std = acts.std(axis=0, keepdims=True)
acts_std[acts_std < 1e-8] = 1.0
acts_normed = (acts - acts_mean) / acts_std

print("Training layer11_d2048_lam10.0...", flush=True)
d_hidden = 2048
lam = 10.0
dataset = TensorDataset(torch.tensor(acts_normed, dtype=torch.float32))
loader = DataLoader(dataset, batch_size=512, shuffle=True, drop_last=True)
sae = SparseAutoencoder(D_MODEL, d_hidden)
optimizer = torch.optim.Adam(sae.parameters(), lr=3e-4)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=40)

for epoch in range(40):
    sae.train()
    for (batch,) in loader:
        loss, mse, f = sae.loss(batch, lam)
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(sae.parameters(), 1.0)
        optimizer.step()
        with torch.no_grad():
            sae.W_dec.weight.data = nn.functional.normalize(sae.W_dec.weight.data, dim=0)
    scheduler.step()
    if epoch % 10 == 0 or epoch == 39:
        sae.eval()
        with torch.no_grad():
            _, f = sae(torch.tensor(acts_normed[:2000], dtype=torch.float32))
            avg_l0 = (f > 0).float().sum(dim=1).mean().item()
            alive = ((f > 0).float().mean(dim=0) > 0.01).sum().item()
        print(f"  epoch {epoch}: L0={avg_l0:.1f} alive={alive}/2048", flush=True)

# Save
sae.eval()
name = "layer11_d2048_lam10.0"
torch.save({"state_dict": sae.state_dict(), "d_input": D_MODEL, "d_hidden": d_hidden, "lambda": lam, "layer": 11}, OUT_DIR / f"{name}.pt")

# Features on training data
with torch.no_grad():
    fl = []
    for s in range(0, len(acts_normed), 512):
        _, f = sae(torch.tensor(acts_normed[s:s+512], dtype=torch.float32))
        fl.append(f.numpy())
    features = np.concatenate(fl)
np.save(OUT_DIR / f"{name}_features.npy", features)

# Cell-level features
cell_acts = np.load(ACT_DIR / "layer_11_cell_activations.npy")
cell_normed = (cell_acts - acts_mean) / acts_std
with torch.no_grad():
    cl = []
    for s in range(0, len(cell_normed), 512):
        _, f = sae(torch.tensor(cell_normed[s:s+512], dtype=torch.float32))
        cl.append(f.numpy())
    cell_features = np.concatenate(cl)
np.save(OUT_DIR / f"{name}_cell_features.npy", cell_features)

# R2
ss_res, ss_tot = 0.0, 0.0
dm = acts_normed.mean(axis=0)
with torch.no_grad():
    for s in range(0, len(acts_normed), 512):
        b = torch.tensor(acts_normed[s:s+512], dtype=torch.float32)
        xh, _ = sae(b)
        ss_res += (b - xh).pow(2).sum().item()
        ss_tot += (b - torch.tensor(dm)).pow(2).sum().item()
r2 = 1 - ss_res / ss_tot

alive_frac = ((features > 0).mean(axis=0) > 0.01).mean()
avg_l0 = (features > 0).sum(axis=1).mean()
print(f"\nFinal: R²={r2:.4f} L0={avg_l0:.1f} alive={alive_frac:.2%}", flush=True)

# Now compile ALL results
results = {}
all_data = [
    ("layer00_d2048_lam1.0", 0, 2048, 1.0),
    ("layer00_d2048_lam3.0", 0, 2048, 3.0),
    ("layer00_d2048_lam10.0", 0, 2048, 10.0),
    ("layer06_d2048_lam1.0", 6, 2048, 1.0),
    ("layer06_d2048_lam3.0", 6, 2048, 3.0),
    ("layer06_d2048_lam10.0", 6, 2048, 10.0),
    ("layer11_d2048_lam1.0", 11, 2048, 1.0),
    ("layer11_d2048_lam3.0", 11, 2048, 3.0),
    ("layer11_d2048_lam10.0", 11, 2048, 10.0),
]

for model_name, layer, dh, l in all_data:
    feat_path = OUT_DIR / f"{model_name}_features.npy"
    if not feat_path.exists():
        continue
    feats = np.load(feat_path)
    af = ((feats > 0).mean(axis=0) > 0.01).mean()
    al0 = (feats > 0).sum(axis=1).mean()
    dead = int((1 - af) * dh)
    
    # Quick R2 from model
    pt_path = OUT_DIR / f"{model_name}.pt"
    ckpt = torch.load(pt_path, map_location="cpu")
    m = SparseAutoencoder(D_MODEL, dh)
    m.load_state_dict(ckpt["state_dict"])
    m.eval()
    
    cell_path = ACT_DIR / f"layer_{layer:02d}_cell_activations.npy"
    ca = np.load(cell_path)
    mn = np.load(OUT_DIR / f"layer_{layer:02d}_mean.npy")
    st = np.load(OUT_DIR / f"layer_{layer:02d}_std.npy")
    ca_n = (ca - mn) / st
    
    ssr, sst = 0.0, 0.0
    camean = ca_n.mean(axis=0)
    with torch.no_grad():
        for s in range(0, len(ca_n), 512):
            b = torch.tensor(ca_n[s:s+512], dtype=torch.float32)
            xh, _ = m(b)
            ssr += (b - xh).pow(2).sum().item()
            sst += (b - torch.tensor(camean)).pow(2).sum().item()
    cell_r2 = 1 - ssr / sst
    
    results[model_name] = {
        "layer": layer, "d_hidden": dh, "lambda": l,
        "R2_token": float(r2) if model_name == "layer11_d2048_lam10.0" else None,
        "R2_cell": float(cell_r2),
        "alive_fraction": float(af),
        "avg_L0": float(al0),
        "n_dead": dead,
    }
    print(f"{model_name}: R²_cell={cell_r2:.4f} L0={al0:.1f} alive={af:.2%} dead={dead}", flush=True)

with open(OUT_DIR / "training_results_sparse.json", "w") as f:
    json.dump(results, f, indent=2)

print("\nDone!", flush=True)
