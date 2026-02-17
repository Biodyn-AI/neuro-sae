import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from pathlib import Path

BASE = Path(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp")
ATLAS_DIR = BASE / "outputs" / "head_layer_checkpoint_replication"

with open(ATLAS_DIR / "gene_list_1200.json") as f:
    genes = json.load(f)
gene_to_idx = {g: i for i, g in enumerate(genes)}
n_genes = len(genes)

# TRRUST
df = pd.read_csv(BASE / "external" / "networks" / "trrust_human.tsv", sep="\t", header=None,
                 names=["source","target","mode","pmid"])
trrust_edges = set()
for _, row in df.iterrows():
    s, t = row["source"], row["target"]
    if s in gene_to_idx and t in gene_to_idx:
        trrust_edges.add((gene_to_idx[s], gene_to_idx[t]))
print(f"TRRUST: {len(trrust_edges)} edges in vocabulary")

# DoRothEA
df2 = pd.read_csv(BASE / "external" / "networks" / "dorothea_chipseq_human.tsv", sep="\t")
doro_edges = set()
for _, row in df2.iterrows():
    s, t = str(row["source"]), str(row["target"])
    if s in gene_to_idx and t in gene_to_idx:
        doro_edges.add((gene_to_idx[s], gene_to_idx[t]))
print(f"DoRothEA: {len(doro_edges)} edges in vocabulary")

def eval_auroc(attn, pos_edges, seed=42):
    pos_list = list(pos_edges)
    rng = np.random.default_rng(seed)
    n_neg = min(len(pos_list) * 10, 50000)
    neg = set()
    while len(neg) < n_neg:
        i, j = rng.integers(0, n_genes, size=2)
        if i != j and (int(i),int(j)) not in pos_edges and (int(i),int(j)) not in neg:
            neg.add((int(i),int(j)))
    all_p = pos_list + list(neg)
    labels = np.array([1]*len(pos_list) + [0]*len(neg))
    si = np.array([p[0] for p in all_p])
    ti = np.array([p[1] for p in all_p])
    
    results = {}
    for layer in range(attn.shape[0]):
        sc = attn[layer].mean(axis=0)[si, ti].astype(float)
        results["L%d" % layer] = roc_auc_score(labels, sc) if np.std(sc) > 0 else 0.5
    sc_all = attn.mean(axis=(0,1))[si, ti].astype(float)
    results["all_mean"] = roc_auc_score(labels, sc_all) if np.std(sc_all) > 0 else 0.5
    return results

all_results = {}
for tissue in ["brain", "kidney", "whole_human"]:
    p = ATLAS_DIR / tissue / "attention_scores_head_layer.npy"
    print("\n=== %s ===" % tissue)
    attn = np.load(p)
    print("  Shape: %s" % str(attn.shape))
    
    for ref_name, ref_edges in [("trrust", trrust_edges), ("dorothea", doro_edges)]:
        res = eval_auroc(attn, ref_edges)
        layer_keys = ["L%d" % i for i in range(12)]
        best_l = max(layer_keys, key=lambda k: res[k])
        am = res["all_mean"]
        bv = res[best_l]
        print("  %s: best=%s (%.4f), all_mean=%.4f" % (ref_name, best_l, bv, am))
        for k in layer_keys:
            print("    %s: %.4f" % (k, res[k]))
        all_results.setdefault(tissue, {})[ref_name] = res
    del attn

with open(r"D:\openclaw\biodyn-nmi-paper\multi_tissue_layer_auroc.json", "w") as f:
    json.dump(all_results, f, indent=2)
print("\nSaved results.")
