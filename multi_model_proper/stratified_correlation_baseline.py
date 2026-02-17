"""Stratified correlation baseline: within-cell-type Spearman correlations vs TRRUST."""

import json, os, sys, warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')

# 1. Load data
print("Loading DLPFC data...")
adata = sc.read_h5ad(r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad")
print(f"Full dataset: {adata.shape}")

# Convert Ensembl to gene symbols
if adata.var_names[0].startswith('ENSG') and 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].values
    adata.var_names_make_unique()

ct_col = 'cell_type'
print(f"Cell types: {dict(adata.obs[ct_col].value_counts())}")

# Subsample to 497 cells
np.random.seed(42)
sc.pp.subsample(adata, n_obs=497, random_state=42)
print(f"Subsampled to {adata.shape[0]} cells")

# 2. Load TRRUST
trrust_path = r"D:\openclaw\biodyn-nmi-paper\multi_model_proper\trrust_rawdata.human.tsv"
if not os.path.exists(trrust_path):
    import urllib.request
    urllib.request.urlretrieve("https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv", trrust_path)

trrust = pd.read_csv(trrust_path, sep='\t', header=None, names=['TF','Target','Type','PMID'])
trrust = trrust.drop_duplicates(subset=['TF','Target'])

gene_names = set(adata.var_names)
trrust_f = trrust[trrust['TF'].isin(gene_names) & trrust['Target'].isin(gene_names)]
trrust_edges = set(zip(trrust_f['TF'], trrust_f['Target']))
print(f"TRRUST edges in data: {len(trrust_edges)}")

# Use top 1000 HVGs (matching paper setup)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1000)
adata = adata[:, adata.var['highly_variable']].copy()
print(f"After HVG selection: {adata.shape}")

gene_names = set(adata.var_names)
trrust_f = trrust[trrust['TF'].isin(gene_names) & trrust['Target'].isin(gene_names)]
trrust_edges = set(zip(trrust_f['TF'], trrust_f['Target']))
tfs = sorted(trrust_f['TF'].unique())
targets = sorted(set(trrust_f['Target'].unique()))
print(f"TRRUST edges after HVG: {len(trrust_edges)}, TFs: {len(tfs)}, Targets: {len(targets)}")

# 3. Get expression
if hasattr(adata.X, 'toarray'):
    X = adata.X.toarray()
else:
    X = np.array(adata.X)

gene_list = list(adata.var_names)
gene_to_idx = {g: i for i, g in enumerate(gene_list)}
cell_types = adata.obs[ct_col].values
unique_cts = [ct for ct in np.unique(cell_types) if (cell_types == ct).sum() >= 5]
print(f"Cell types with >=5 cells: {len(unique_cts)}")

# 4. For AUROC: score ALL possible TF-target pairs (TFs x all genes)
# This matches how attention AUROC is computed
all_target_genes = sorted(gene_names)
print(f"\nScoring {len(tfs)} TFs x {len(all_target_genes)} genes = {len(tfs)*len(all_target_genes)} pairs...")

labels = []
scores_pooled = []
scores_strat_max = []
scores_strat_mean = []

for i, tf in enumerate(tfs):
    if i % 20 == 0:
        print(f"  TF {i+1}/{len(tfs)}: {tf}")
    tf_idx = gene_to_idx[tf]
    tf_expr = X[:, tf_idx]
    
    for tgt in all_target_genes:
        if tf == tgt:
            continue
        tgt_idx = gene_to_idx[tgt]
        tgt_expr = X[:, tgt_idx]
        
        labels.append(1 if (tf, tgt) in trrust_edges else 0)
        
        # Pooled
        r, _ = spearmanr(tf_expr, tgt_expr)
        scores_pooled.append(abs(r) if not np.isnan(r) else 0)
        
        # Stratified
        ct_corrs = []
        ct_sizes = []
        for ct in unique_cts:
            mask = cell_types == ct
            n_ct = mask.sum()
            r_ct, _ = spearmanr(tf_expr[mask], tgt_expr[mask])
            if not np.isnan(r_ct):
                ct_corrs.append(abs(r_ct))
                ct_sizes.append(n_ct)
        
        if ct_corrs:
            scores_strat_max.append(max(ct_corrs))
            total = sum(ct_sizes)
            scores_strat_mean.append(sum(c*s/total for c,s in zip(ct_corrs, ct_sizes)))
        else:
            scores_strat_max.append(0)
            scores_strat_mean.append(0)

labels = np.array(labels)
print(f"\nPairs: {len(labels)}, Positive: {labels.sum()}, Negative: {(1-labels).sum()}")

auroc_pool = roc_auc_score(labels, np.array(scores_pooled))
auroc_max = roc_auc_score(labels, np.array(scores_strat_max))
auroc_mean = roc_auc_score(labels, np.array(scores_strat_mean))

print(f"\n{'='*50}")
print(f"Pooled Spearman AUROC:           {auroc_pool:.4f}")
print(f"Stratified (CSSI-max) AUROC:     {auroc_max:.4f}")
print(f"Stratified (CSSI-mean) AUROC:    {auroc_mean:.4f}")
print(f"L13-L14 attention AUROC:         0.694")
print(f"Gap (attention - strat_max):     {0.694 - auroc_max:.4f}")
print(f"{'='*50}")

results = {
    "n_cells": int(adata.shape[0]),
    "n_genes": int(adata.shape[1]),
    "n_cell_types": len(unique_cts),
    "cell_types": {str(ct): int((cell_types==ct).sum()) for ct in unique_cts},
    "n_trrust_edges_evaluated": int(labels.sum()),
    "n_total_pairs": int(len(labels)),
    "auroc_pooled_spearman": round(float(auroc_pool), 4),
    "auroc_stratified_max": round(float(auroc_max), 4),
    "auroc_stratified_mean": round(float(auroc_mean), 4),
    "auroc_L13_L14_attention": 0.694,
}

with open(r"D:\openclaw\biodyn-nmi-paper\multi_model_proper\stratified_correlation_baseline.json", 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved.")
