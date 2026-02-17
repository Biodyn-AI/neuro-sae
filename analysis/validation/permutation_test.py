"""Permutation null - use normal approximation for null + 1000 perms for validation."""
import json, numpy as np, pandas as pd
from scipy import stats
np.random.seed(42)
print("Start", flush=True)

with open(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\outputs\head_layer_checkpoint_replication\gene_list_1200.json") as f:
    gene_list = json.load(f)
gene_to_idx = {g: i for i, g in enumerate(gene_list)}
n_genes = len(gene_list)

trrust = pd.read_csv(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\external\networks\trrust_human.tsv", 
                      sep='\t', header=None, names=['tf','target','mode','pmid'])
trrust_edges = set()
for _, row in trrust.iterrows():
    if row['tf'] in gene_to_idx and row['target'] in gene_to_idx:
        i, j = gene_to_idx[row['tf']], gene_to_idx[row['target']]
        if i != j: trrust_edges.add((i, j))
n_pos = len(trrust_edges)
print(f"TRRUST edges: {n_pos}", flush=True)

# The JSON reports per-layer AUROCs. The `all_mean` is likely mean of per-layer AUROCs.
# Let me replicate per-layer AUROC computation.
print("Memory-mapping attention...", flush=True)
attn = np.lib.format.open_memmap(
    r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\outputs\head_layer_checkpoint_replication\brain\attention_scores_head_layer.npy",
    mode='r')

mask = ~np.eye(n_genes, dtype=bool)
n_total = n_genes * (n_genes - 1)
n_neg = n_total - n_pos

# Compute per-layer AUROC (mean across heads per layer)
print("Computing per-layer AUROCs...", flush=True)
layer_aurocs = []
layer_ranks_list = []  # store for permutation test

for L in range(12):
    layer_mean = np.zeros((n_genes, n_genes), dtype=np.float64)
    for H in range(8):
        layer_mean += attn[L, H].astype(np.float64)
    layer_mean /= 8
    
    scores = layer_mean[mask]
    ranks = scores.argsort().argsort().astype(np.float64) + 1
    
    pos_flat = []
    for (i, j) in trrust_edges:
        flat = i * (n_genes - 1) + (j if j < i else j - 1)
        pos_flat.append(flat)
    pos_flat = np.array(pos_flat)
    
    rank_sum = ranks[pos_flat].sum()
    U = rank_sum - n_pos * (n_pos + 1) / 2
    auroc = U / (n_pos * n_neg)
    layer_aurocs.append(auroc)
    print(f"  L{L}: AUROC = {auroc:.6f}", flush=True)
    
    if L == 0:
        # Save ranks for L0 permutation test (representative)
        saved_ranks = ranks.copy()
        saved_pos_flat = pos_flat.copy()

print(f"\nMean of per-layer AUROCs: {np.mean(layer_aurocs):.6f}", flush=True)

# The JSON all_mean for brain trrust is 0.7175 - let's see if per-layer AUROCs match JSON
# JSON L0 brain trrust = 0.7205
# If our L0 matches, great. If not, different scoring method.

# Now do permutation test on L0 (representative layer with highest AUROC)
# Use 1000 perms - should be fast enough with pre-computed ranks
best_layer = int(np.argmax(layer_aurocs))
print(f"\nBest layer: L{best_layer} (AUROC={layer_aurocs[best_layer]:.4f})", flush=True)

# Recompute ranks for best layer
layer_mean = np.zeros((n_genes, n_genes), dtype=np.float64)
for H in range(8):
    layer_mean += attn[best_layer, H].astype(np.float64)
layer_mean /= 8
scores = layer_mean[mask]
ranks = scores.argsort().argsort().astype(np.float64) + 1
observed_auroc = layer_aurocs[best_layer]

# Analytical null: under random labeling, AUROC ~ N(0.5, sigma^2)
# sigma^2 = (n_total + 1) / (12 * n_pos * n_neg)
sigma = np.sqrt((n_total + 1) / (12 * n_pos * n_neg))
z_analytical = (observed_auroc - 0.5) / sigma
p_analytical = 1 - stats.norm.cdf(z_analytical)
print(f"\n=== ANALYTICAL NULL ===")
print(f"Observed AUROC (L{best_layer}): {observed_auroc:.4f}")
print(f"Null: N(0.5, {sigma:.6f})")
print(f"Z-score: {z_analytical:.2f}")
print(f"P-value (one-sided): {p_analytical:.2e}", flush=True)

# Also for all-layer mean AUROC
mean_auroc = np.mean(layer_aurocs)
z_mean = (mean_auroc - 0.5) / sigma
p_mean = 1 - stats.norm.cdf(z_mean)
print(f"\nMean AUROC across layers: {mean_auroc:.4f}")
print(f"Z-score: {z_mean:.2f}, P-value: {p_mean:.2e}", flush=True)

# Quick empirical permutation (1000 perms) for validation
print(f"\nRunning 1000 empirical permutations on L{best_layer}...", flush=True)
n_perms = 1000
null_aurocs = np.zeros(n_perms)
for p in range(n_perms):
    rand_idx = np.random.choice(n_total, size=n_pos, replace=False)
    rank_sum = ranks[rand_idx].sum()
    null_aurocs[p] = (rank_sum - n_pos*(n_pos+1)/2) / (n_pos * n_neg)
    if (p+1) % 200 == 0:
        print(f"  {p+1}/{n_perms}", flush=True)

emp_p = (null_aurocs >= observed_auroc).sum() / n_perms
print(f"\n=== EMPIRICAL PERMUTATION (L{best_layer}) ===")
print(f"Observed: {observed_auroc:.4f}")
print(f"Null: {null_aurocs.mean():.4f} +/- {null_aurocs.std():.4f}")
print(f"Null range: [{null_aurocs.min():.4f}, {null_aurocs.max():.4f}]")
print(f"Empirical P: {emp_p}")
if emp_p == 0:
    print(f"P < {1/n_perms}")
print("\nDONE", flush=True)
