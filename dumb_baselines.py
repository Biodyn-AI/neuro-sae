"""
Dumb baselines: Can simple gene-level statistics explain TRRUST/DoRothEA AUROC?
Baselines:
  1. Detection rate product
  2. Mean expression product  
  3. HVG variance product
Compare against attention-based AUROC on the same edges/negatives.
"""
import json, numpy as np, pandas as pd, scanpy as sc
from sklearn.metrics import roc_auc_score
from pathlib import Path

BASE = Path(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp")
ATLAS_DIR = BASE / "outputs" / "head_layer_checkpoint_replication"

# Load gene list
with open(ATLAS_DIR / "gene_list_1200.json") as f:
    genes = json.load(f)
gene_to_idx = {g: i for i, g in enumerate(genes)}
n_genes = len(genes)

# Load h5ad for gene statistics
print("Loading h5ad...")
adata = sc.read_h5ad(ATLAS_DIR / "immune_shared_hvg1200_processed.h5ad")
print(f"  adata: {adata.shape}")

# Align gene order to our 1200 list
# Check if adata.var_names matches
adata_genes = list(adata.var_names)
print(f"  adata genes: {len(adata_genes)}, first 5: {adata_genes[:5]}")

# Build gene stats for our 1200 genes
from scipy.sparse import issparse
X = adata.X
if issparse(X):
    X_dense = None  # work with sparse
    detection_rate = np.array((X > 0).mean(axis=0)).flatten()
    mean_expr = np.array(X.mean(axis=0)).flatten()
    # Variance: E[X^2] - E[X]^2
    X_sq = X.copy()
    X_sq.data **= 2
    variance = np.array(X_sq.mean(axis=0)).flatten() - mean_expr**2
else:
    detection_rate = (X > 0).mean(axis=0)
    mean_expr = X.mean(axis=0)
    variance = X.var(axis=0)

# Map to our gene indices
adata_gene_to_pos = {g: i for i, g in enumerate(adata_genes)}
det_rate_vec = np.zeros(n_genes)
mean_expr_vec = np.zeros(n_genes)
var_vec = np.zeros(n_genes)
mapped = 0
for g, idx in gene_to_idx.items():
    if g in adata_gene_to_pos:
        pos = adata_gene_to_pos[g]
        det_rate_vec[idx] = detection_rate[pos]
        mean_expr_vec[idx] = mean_expr[pos]
        var_vec[idx] = variance[pos]
        mapped += 1
print(f"Mapped {mapped}/{n_genes} genes from adata")
del adata

# Load TRRUST and DoRothEA
df = pd.read_csv(BASE / "external" / "networks" / "trrust_human.tsv", sep="\t", header=None,
                 names=["source","target","mode","pmid"])
trrust_edges = set()
for _, row in df.iterrows():
    s, t = row["source"], row["target"]
    if s in gene_to_idx and t in gene_to_idx:
        trrust_edges.add((gene_to_idx[s], gene_to_idx[t]))

df2 = pd.read_csv(BASE / "external" / "networks" / "dorothea_chipseq_human.tsv", sep="\t")
doro_edges = set()
for _, row in df2.iterrows():
    s, t = str(row["source"]), str(row["target"])
    if s in gene_to_idx and t in gene_to_idx:
        doro_edges.add((gene_to_idx[s], gene_to_idx[t]))

print(f"TRRUST edges: {len(trrust_edges)}, DoRothEA edges: {len(doro_edges)}")

# Also try FULL TRRUST (all genes, not just 1200 HVGs)
all_trrust_genes = set()
for _, row in df.iterrows():
    all_trrust_genes.add(row["source"])
    all_trrust_genes.add(row["target"])
trrust_in_vocab = sum(1 for g in all_trrust_genes if g in gene_to_idx)
print(f"TRRUST unique genes: {len(all_trrust_genes)}, in vocab: {trrust_in_vocab}")

# Evaluation function
def eval_all(pos_edges, ref_name, attn_tissues, seed=42):
    pos_list = list(pos_edges)
    rng = np.random.default_rng(seed)
    n_neg = min(len(pos_list) * 10, 50000)
    neg = set()
    while len(neg) < n_neg:
        i, j = rng.integers(0, n_genes, size=2)
        if i != j and (int(i),int(j)) not in pos_edges and (int(i),int(j)) not in neg:
            neg.add((int(i),int(j)))
    all_pairs = pos_list + list(neg)
    labels = np.array([1]*len(pos_list) + [0]*len(neg))
    si = np.array([p[0] for p in all_pairs])
    ti = np.array([p[1] for p in all_pairs])
    
    print(f"\n--- {ref_name} ({len(pos_list)} pos, {len(neg)} neg) ---")
    
    # Dumb baselines
    for bname, vec in [("det_rate_product", det_rate_vec), 
                        ("mean_expr_product", mean_expr_vec),
                        ("variance_product", var_vec)]:
        scores = vec[si] * vec[ti]
        auc = roc_auc_score(labels, scores) if np.std(scores) > 0 else 0.5
        print(f"  {bname}: AUROC = {auc:.4f}")
    
    # Attention baselines per tissue
    for tissue, attn in attn_tissues.items():
        # Best layer and pooled
        sc_all = attn.mean(axis=(0,1))[si, ti].astype(float)
        auc_pool = roc_auc_score(labels, sc_all) if np.std(sc_all) > 0 else 0.5
        best_auc = 0
        best_l = -1
        for layer in range(attn.shape[0]):
            sc_l = attn[layer].mean(axis=0)[si, ti].astype(float)
            a = roc_auc_score(labels, sc_l) if np.std(sc_l) > 0 else 0.5
            if a > best_auc:
                best_auc = a
                best_l = layer
        print(f"  attention_{tissue}_pooled: AUROC = {auc_pool:.4f}")
        print(f"  attention_{tissue}_best_L{best_l}: AUROC = {best_auc:.4f}")

# Load attention for one tissue to compare
print("\nLoading attention matrices...")
attn_tissues = {}
for tissue in ["brain", "kidney", "whole_human"]:
    p = ATLAS_DIR / tissue / "attention_scores_head_layer.npy"
    attn_tissues[tissue] = np.load(p)
    print(f"  {tissue}: {attn_tissues[tissue].shape}")

for ref_name, ref_edges in [("TRRUST", trrust_edges), ("DoRothEA", doro_edges)]:
    eval_all(ref_edges, ref_name, attn_tissues)

# Cleanup
del attn_tissues
print("\n=== DONE ===")
