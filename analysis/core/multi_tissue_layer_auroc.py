"""
Multi-tissue layer-stratified AUROC analysis for NMI paper.
"""
import numpy as np
import json
import os
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

BASE = r'D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp'
HLC = f'{BASE}\\outputs\\head_layer_checkpoint_replication'
TRRUST = f'{BASE}\\external\\networks\\trrust_human.tsv'

# Load gene list
with open(f'{HLC}\\gene_list_1200.json') as f:
    genes = json.load(f)
gene2idx = {g: i for i, g in enumerate(genes)}
n_genes = len(genes)

# Load TRRUST and map
trrust_mapped = set()
with open(TRRUST) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2 and parts[0] in gene2idx and parts[1] in gene2idx:
            trrust_mapped.add((gene2idx[parts[0]], gene2idx[parts[1]]))

labels = np.zeros((n_genes, n_genes), dtype=np.int8)
for i, j in trrust_mapped:
    labels[i, j] = 1

mask = ~np.eye(n_genes, dtype=bool)
y_true = labels[mask]
print(f"TRRUST mapped: {len(trrust_mapped)} edges, pos={y_true.sum()}, total={y_true.size}")

def layer_aurocs(scores_path, counts_path, name):
    scores = np.load(scores_path)  # (12, 8, 1200, 1200)
    counts = np.load(counts_path)  # (1200, 1200)
    counts_safe = np.maximum(counts, 1)
    n_layers = scores.shape[0]
    
    print(f"\n=== {name} ===")
    results = {}
    for L in range(n_layers):
        s = scores[L].mean(axis=0) / counts_safe
        s = np.nan_to_num(s, nan=0.0, posinf=0.0, neginf=0.0)
        a = roc_auc_score(y_true, s[mask])
        
        # Best head
        best_h = max(roc_auc_score(y_true, np.nan_to_num(scores[L, h] / counts_safe, nan=0.0, posinf=0.0, neginf=0.0)[mask]) for h in range(scores.shape[1]))
        
        results[L] = {'auroc': a, 'best_head': best_h}
        print(f"  L{L:2d}: AUROC={a:.4f}  best_head={best_h:.4f}")
    
    return results

tissues = {}
for name, subdir in [('Brain', 'brain'), ('Kidney', 'kidney'), ('Whole Human', 'whole_human')]:
    sp = f'{HLC}\\{subdir}\\attention_scores_head_layer.npy'
    cp = f'{HLC}\\{subdir}\\attention_counts_head_layer.npy'
    if os.path.exists(sp):
        tissues[name] = (sp, cp)

all_results = {}
for name, (sp, cp) in tissues.items():
    all_results[name] = layer_aurocs(sp, cp, name)

# Comparison table
print("\n\n=== LAYER-STRATIFIED AUROC COMPARISON ===")
header = f"{'Layer':<8}" + "".join(f"{t:<20}" for t in all_results)
print(header)
print("-" * len(header))
for L in range(12):
    row = f"L{L:<7}"
    for t in all_results:
        row += f"{all_results[t][L]['auroc']:<20.4f}"
    print(row)

# Depth-dependent pattern
print("\n=== DEPTH-DEPENDENT PATTERN ===")
from scipy.stats import spearmanr
for t, res in all_results.items():
    aurocs = [res[L]['auroc'] for L in range(12)]
    rho, p = spearmanr(range(12), aurocs)
    best = np.argmax(aurocs)
    print(f"{t}: rho={rho:.3f} (p={p:.4f}), best=L{best} ({aurocs[best]:.4f}), early={np.mean(aurocs[:4]):.4f}, late={np.mean(aurocs[8:]):.4f}")

# Bootstrap CI on brain best layer
print("\n=== BOOTSTRAP CIs (brain best layer) ===")
brain_scores = np.load(tissues['Brain'][0])
brain_counts = np.maximum(np.load(tissues['Brain'][1]), 1)
brain_aurocs = [all_results['Brain'][L]['auroc'] for L in range(12)]
best_L = np.argmax(brain_aurocs)
best_flat = np.nan_to_num((brain_scores[best_L].mean(axis=0) / brain_counts)[mask], nan=0.0, posinf=0.0, neginf=0.0)

rng = np.random.RandomState(42)
n = len(y_true)
n_boot = 2000
boot_attn = np.array([roc_auc_score(y_true[idx := rng.choice(n, n, replace=True)], best_flat[idx]) for _ in range(n_boot)])
print(f"Attention L{best_L}: {np.mean(boot_attn):.4f} [{np.percentile(boot_attn, 2.5):.4f}, {np.percentile(boot_attn, 97.5):.4f}]")

# Correlation baseline using brain processed h5ad
print("\nComputing correlation baselines...")
import scanpy as sc

def corr_baseline(h5_path, name):
    adata = sc.read_h5ad(h5_path)
    shared = [g for g in genes if g in adata.var_names]
    shared_idx = np.array([gene2idx[g] for g in shared])
    adata_sub = adata[:, shared]
    if adata_sub.n_obs > 500:
        adata_sub = adata_sub[rng.choice(adata_sub.n_obs, 500, replace=False)]
    X = adata_sub.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    corr = np.abs(np.corrcoef(X.T))
    corr = np.nan_to_num(corr, nan=0.0)
    
    # Map to full matrix using advanced indexing
    full_corr = np.zeros((n_genes, n_genes))
    ix = np.ix_(shared_idx, shared_idx)
    full_corr[ix] = corr
    
    flat = full_corr[mask]
    auroc = roc_auc_score(y_true, flat)
    print(f"{name} correlation AUROC: {auroc:.4f} ({len(shared)} shared genes, {adata_sub.n_obs} cells)")
    
    # Bootstrap
    boot = np.array([roc_auc_score(y_true[idx := rng.choice(n, n, replace=True)], flat[idx]) for _ in range(n_boot)])
    print(f"  CI: [{np.percentile(boot, 2.5):.4f}, {np.percentile(boot, 97.5):.4f}]")
    return auroc, boot

brain_corr, boot_corr_brain = corr_baseline(f'{BASE}\\outputs\\tabula_sapiens_processed.h5ad', 'Brain')

# Difference CI
diff = boot_attn - boot_corr_brain
print(f"\nAUROC difference (attn L{best_L} - corr): {np.mean(diff):.4f} [{np.percentile(diff, 2.5):.4f}, {np.percentile(diff, 97.5):.4f}]")
print(f"P(attn <= corr): {np.mean(diff <= 0):.4f}")

# Kidney correlation
kidney_h5 = f'{BASE}\\data\\raw\\tabula_sapiens_kidney.h5ad'
if os.path.exists(kidney_h5):
    corr_baseline(kidney_h5, 'Kidney')

# Save all results
output = {
    'layer_aurocs': {t: {str(L): r for L, r in res.items()} for t, res in all_results.items()},
    'n_trrust_mapped': len(trrust_mapped),
    'n_genes': n_genes,
    'brain_best_layer': int(best_L),
    'brain_best_auroc': float(brain_aurocs[best_L]),
    'brain_corr_auroc': float(brain_corr),
}
with open(f'{HLC}\\multi_tissue_results.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nSaved to {HLC}\\multi_tissue_results.json")
