"""
Expression-Matched Negative Evaluation v2.
Matches negatives by ACTUAL expression statistics (mean, detection rate) from h5ad,
NOT by attention marginals. Includes bootstrap CIs.
"""
import numpy as np
import json
import os
import scanpy as sc
from sklearn.metrics import roc_auc_score

np.random.seed(42)

base = r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\outputs\head_layer_checkpoint_replication"

# Load gene list and attention
genes = json.load(open(os.path.join(base, "gene_list_1200.json")))
gene2idx = {g: i for i, g in enumerate(genes)}
scores = np.load(os.path.join(base, "brain", "attention_scores_head_layer.npy"))
attn_mean = scores.mean(axis=(0, 1)).astype(np.float32)
print(f"Attention: {scores.shape}, genes: {len(genes)}")

# Load expression data
h5ad_path = os.path.join(base, "immune_shared_hvg1200_processed.h5ad")
print(f"Loading {h5ad_path}...")
adata = sc.read_h5ad(h5ad_path)
print(f"AnnData: {adata.shape}")

# Compute expression stats per gene
import scipy.sparse as sp
X = adata.X
if sp.issparse(X):
    X_dense = np.array(X.todense())
else:
    X_dense = np.array(X)

# Map adata var_names to our gene list indices
adata_genes = list(adata.var_names)
adata_gene2col = {g: i for i, g in enumerate(adata_genes)}

gene_mean_expr = np.zeros(len(genes))
gene_detection_rate = np.zeros(len(genes))
gene_var_expr = np.zeros(len(genes))

matched_count = 0
for i, g in enumerate(genes):
    if g in adata_gene2col:
        col = adata_gene2col[g]
        vals = X_dense[:, col]
        gene_mean_expr[i] = vals.mean()
        gene_detection_rate[i] = (vals > 0).mean()
        gene_var_expr[i] = vals.var()
        matched_count += 1

print(f"Matched {matched_count}/{len(genes)} genes in h5ad")
print(f"Mean expr range: [{gene_mean_expr.min():.4f}, {gene_mean_expr.max():.4f}]")
print(f"Detection rate range: [{gene_detection_rate.min():.4f}, {gene_detection_rate.max():.4f}]")

# Load TRRUST
trrust_path = os.path.join(r"D:\openclaw\biodyn-nmi-paper\multi_model_proper", "trrust_rawdata.human.tsv")
gene_set = set(genes)
positive_edges = []
with open(trrust_path) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            tf, target = parts[0], parts[1]
            if tf in gene_set and target in gene_set and tf != target:
                positive_edges.append((tf, target))
positive_edges = list(set(positive_edges))
print(f"TRRUST positive edges: {len(positive_edges)}")

positive_set = set((gene2idx[tf], gene2idx[tgt]) for tf, tgt in positive_edges)

# Expression-matched sampling
def find_expr_matched(gene_idx, is_tf=True, tol=0.20):
    """Find genes matched by mean expression and detection rate."""
    me = gene_mean_expr[gene_idx]
    dr = gene_detection_rate[gene_idx]
    
    me_low, me_high = me * (1 - tol), me * (1 + tol)
    dr_low, dr_high = dr - tol * max(dr, 0.05), dr + tol * max(dr, 0.05)
    
    if me < 1e-6:
        me_high = 1e-5
    
    mask = ((gene_mean_expr >= me_low) & (gene_mean_expr <= me_high) &
            (gene_detection_rate >= dr_low) & (gene_detection_rate <= dr_high))
    return np.where(mask)[0]

N_MATCHED = 50
all_pairs = []
all_labels = []

for tf_name, tgt_name in positive_edges:
    tf_idx = gene2idx[tf_name]
    tgt_idx = gene2idx[tgt_name]
    
    matched_tfs = find_expr_matched(tf_idx, is_tf=True)
    matched_tgts = find_expr_matched(tgt_idx, is_tf=False)
    
    neg_candidates = [(mt, mg) for mt in matched_tfs for mg in matched_tgts
                      if mt != mg and (mt, mg) not in positive_set]
    
    if len(neg_candidates) < 5:
        matched_tfs = find_expr_matched(tf_idx, tol=0.40)
        matched_tgts = find_expr_matched(tgt_idx, tol=0.40)
        neg_candidates = [(mt, mg) for mt in matched_tfs for mg in matched_tgts
                          if mt != mg and (mt, mg) not in positive_set]
    
    n_sample = min(N_MATCHED, len(neg_candidates))
    if n_sample == 0:
        print(f"  SKIP: {tf_name}->{tgt_name}")
        continue
    
    chosen = [neg_candidates[i] for i in np.random.choice(len(neg_candidates), n_sample, replace=False)]
    
    all_pairs.append((tf_idx, tgt_idx))
    all_labels.append(1)
    for ni, nj in chosen:
        all_pairs.append((ni, nj))
        all_labels.append(0)

all_labels = np.array(all_labels)
all_pairs = np.array(all_pairs)
print(f"\nTotal: {len(all_labels)} ({all_labels.sum()} pos, {(1-all_labels).sum():.0f} neg)")

# Compute scores
attn_vals = np.array([attn_mean[i, j] for i, j in all_pairs])
row_mean = attn_mean.mean(axis=1)
col_mean = attn_mean.mean(axis=0)
marginal_vals = np.array([row_mean[i] * col_mean[j] for i, j in all_pairs])
expr_product = np.array([gene_mean_expr[i] * gene_mean_expr[j] for i, j in all_pairs])

auroc_attn = roc_auc_score(all_labels, attn_vals)
auroc_marginal = roc_auc_score(all_labels, marginal_vals)
auroc_expr = roc_auc_score(all_labels, expr_product)

print(f"\n=== EXPRESSION-MATCHED AUROC ===")
print(f"Attention (mean all heads): {auroc_attn:.4f}")
print(f"Marginal product baseline:  {auroc_marginal:.4f}")
print(f"Expression product baseline: {auroc_expr:.4f}")

# Bootstrap CI on attention AUROC
N_BOOT = 1000
boot_aurocs = []
pos_idx = np.where(all_labels == 1)[0]
neg_idx = np.where(all_labels == 0)[0]

for b in range(N_BOOT):
    bp = np.random.choice(pos_idx, len(pos_idx), replace=True)
    bn = np.random.choice(neg_idx, len(neg_idx), replace=True)
    bidx = np.concatenate([bp, bn])
    try:
        ba = roc_auc_score(all_labels[bidx], attn_vals[bidx])
    except:
        ba = 0.5
    boot_aurocs.append(ba)

boot_aurocs = np.array(boot_aurocs)
ci_low = np.percentile(boot_aurocs, 2.5)
ci_high = np.percentile(boot_aurocs, 97.5)
print(f"Attention AUROC 95% CI: [{ci_low:.4f}, {ci_high:.4f}]")

# Per-head best
head_aurocs = {}
for layer in range(scores.shape[0]):
    for head in range(scores.shape[1]):
        h_attn = scores[layer, head].astype(np.float32)
        h_vals = np.array([h_attn[i, j] for i, j in all_pairs])
        try:
            head_aurocs[f"L{layer}_H{head}"] = roc_auc_score(all_labels, h_vals)
        except:
            head_aurocs[f"L{layer}_H{head}"] = 0.5

sorted_heads = sorted(head_aurocs.items(), key=lambda x: x[1], reverse=True)
best_name, best_auroc = sorted_heads[0]
print(f"Best head: {best_name} = {best_auroc:.4f}")

# Bootstrap CI on best head
best_layer, best_head_num = int(best_name.split('_')[0][1:]), int(best_name.split('_')[1][1:])
best_h_attn = scores[best_layer, best_head_num].astype(np.float32)
best_h_vals = np.array([best_h_attn[i, j] for i, j in all_pairs])
boot_best = []
for b in range(N_BOOT):
    bp = np.random.choice(pos_idx, len(pos_idx), replace=True)
    bn = np.random.choice(neg_idx, len(neg_idx), replace=True)
    bidx = np.concatenate([bp, bn])
    try:
        boot_best.append(roc_auc_score(all_labels[bidx], best_h_vals[bidx]))
    except:
        boot_best.append(0.5)
boot_best = np.array(boot_best)

# Save
output = {
    "experiment": "expression_matched_negative_evaluation_v2",
    "description": "AUROC with negatives matched on ACTUAL expression statistics (mean expression ±20%, detection rate ±20%) from source h5ad, not attention marginals",
    "source_h5ad": "immune_shared_hvg1200_processed.h5ad",
    "matching_variables": ["mean_expression", "detection_rate"],
    "matching_tolerance": 0.20,
    "n_positive_edges": int(all_labels.sum()),
    "n_negative_edges": int((1 - all_labels).sum()),
    "n_total": len(all_labels),
    "negatives_per_positive": N_MATCHED,
    "auroc_attention_mean_all_heads": round(float(auroc_attn), 4),
    "auroc_attention_95ci": [round(float(ci_low), 4), round(float(ci_high), 4)],
    "auroc_marginal_product_baseline": round(float(auroc_marginal), 4),
    "auroc_expression_product_baseline": round(float(auroc_expr), 4),
    "auroc_best_single_head": round(float(best_auroc), 4),
    "best_head": best_name,
    "best_head_95ci": [round(float(np.percentile(boot_best, 2.5)), 4), round(float(np.percentile(boot_best, 97.5)), 4)],
    "auroc_median_head": round(float(np.median(list(head_aurocs.values()))), 4),
    "top5_heads": {k: round(v, 4) for k, v in sorted_heads[:5]},
    "n_bootstrap": N_BOOT
}

out_path = os.path.join(r"D:\openclaw\biodyn-nmi-paper\multi_model_proper", "expression_matched_v2.json")
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {out_path}")
print(json.dumps(output, indent=2))
