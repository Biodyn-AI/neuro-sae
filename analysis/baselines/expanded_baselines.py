"""
Expanded evaluation: dumb baselines vs attention on multiple edge sets.
Including full DoRothEA (all confidence levels) = 483 edges.
Also try downloading OmniPath TF-target interactions.
"""
import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from pathlib import Path
import urllib.request

BASE = Path(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp")
ATLAS_DIR = BASE / "outputs" / "head_layer_checkpoint_replication"

with open(ATLAS_DIR / "gene_list_1200.json") as f:
    genes = json.load(f)
gene_to_idx = {g: i for i, g in enumerate(genes)}
n_genes = len(genes)
gene_set = set(genes)

# Gene stats from h5ad
import scanpy as sc
from scipy.sparse import issparse
adata = sc.read_h5ad(ATLAS_DIR / "immune_shared_hvg1200_processed.h5ad")
X = adata.X
if issparse(X):
    detection_rate = np.array((X > 0).mean(axis=0)).flatten()
    mean_expr = np.array(X.mean(axis=0)).flatten()
    X_sq = X.copy(); X_sq.data **= 2
    variance = np.array(X_sq.mean(axis=0)).flatten() - mean_expr**2
else:
    detection_rate = (X > 0).mean(axis=0); mean_expr = X.mean(axis=0); variance = X.var(axis=0)
adata_genes = list(adata.var_names)
adata_g2p = {g: i for i, g in enumerate(adata_genes)}
det_vec = np.zeros(n_genes); mean_vec = np.zeros(n_genes); var_vec = np.zeros(n_genes)
for g, idx in gene_to_idx.items():
    if g in adata_g2p:
        p = adata_g2p[g]
        det_vec[idx] = detection_rate[p]; mean_vec[idx] = mean_expr[p]; var_vec[idx] = variance[p]
del adata

# Edge sets
def load_trrust():
    df = pd.read_csv(BASE / "external" / "networks" / "trrust_human.tsv", sep="\t", header=None,
                     names=["source","target","mode","pmid"])
    edges = set()
    for _, r in df.iterrows():
        if r["source"] in gene_to_idx and r["target"] in gene_to_idx:
            edges.add((gene_to_idx[r["source"]], gene_to_idx[r["target"]]))
    return edges

def load_dorothea_full():
    df = pd.read_csv(BASE / "external" / "networks" / "dorothea_human.tsv", sep="\t")
    edges = {}  # by confidence
    for _, r in df.iterrows():
        s, t = str(r["source"]), str(r["target"])
        if s in gene_to_idx and t in gene_to_idx:
            c = str(r["confidence"])[0]  # A, B, C, D
            edges.setdefault(c, set()).add((gene_to_idx[s], gene_to_idx[t]))
    return edges

def try_omnipath():
    """Try to download OmniPath TF-target interactions."""
    url = "https://omnipathdb.org/interactions?datasets=tfregulons&fields=sources,references&genesymbols=yes&format=tsv"
    out = BASE / "external" / "networks" / "omnipath_tfreg.tsv"
    if not out.exists():
        print("Downloading OmniPath TF regulons...")
        try:
            urllib.request.urlretrieve(url, out)
            print(f"  Downloaded to {out}")
        except Exception as e:
            print(f"  Failed: {e}")
            return set()
    df = pd.read_csv(out, sep="\t")
    print(f"OmniPath TF regulons: {len(df)} rows, columns: {list(df.columns)[:8]}")
    edges = set()
    src_col = "source_genesymbol" if "source_genesymbol" in df.columns else "source"
    tgt_col = "target_genesymbol" if "target_genesymbol" in df.columns else "target"
    for _, r in df.iterrows():
        s, t = str(r[src_col]), str(r[tgt_col])
        if s in gene_to_idx and t in gene_to_idx:
            edges.add((gene_to_idx[s], gene_to_idx[t]))
    return edges

trrust = load_trrust()
doro_by_conf = load_dorothea_full()
doro_AB = doro_by_conf.get("A", set()) | doro_by_conf.get("B", set())
doro_ABC = doro_AB | doro_by_conf.get("C", set())
doro_all = doro_ABC | doro_by_conf.get("D", set())
omni = try_omnipath()

edge_sets = {
    "TRRUST (28)": trrust,
    "DoRothEA A+B (27)": doro_AB,
    "DoRothEA A+B+C (65)": doro_ABC,
    "DoRothEA all (483)": doro_all,
}
if omni:
    edge_sets[f"OmniPath TF-reg ({len(omni)})"] = omni

# Load one attention tissue for comparison
print("\nLoading brain attention...")
attn = np.load(ATLAS_DIR / "brain" / "attention_scores_head_layer.npy")
attn_pooled = attn.mean(axis=(0, 1))  # (1200, 1200)
del attn

def evaluate(pos_edges, name, seed=42):
    pos_list = list(pos_edges)
    if len(pos_list) < 5:
        print(f"\n{name}: too few edges ({len(pos_list)}), skipping")
        return
    rng = np.random.default_rng(seed)
    n_neg = min(len(pos_list) * 10, 50000)
    neg = set()
    while len(neg) < n_neg:
        i, j = rng.integers(0, n_genes, size=2)
        if i != j and (int(i),int(j)) not in pos_edges:
            neg.add((int(i),int(j)))
    all_p = pos_list + list(neg)
    labels = np.array([1]*len(pos_list) + [0]*len(neg))
    si = np.array([p[0] for p in all_p])
    ti = np.array([p[1] for p in all_p])
    
    print(f"\n{name}:")
    results = {}
    for bname, vec in [("det_rate", det_vec), ("mean_expr", mean_vec), ("variance", var_vec)]:
        scores = vec[si] * vec[ti]
        auc = roc_auc_score(labels, scores) if np.std(scores) > 0 else 0.5
        results[bname] = auc
        print(f"  {bname}_product: {auc:.4f}")
    
    sc_attn = attn_pooled[si, ti].astype(float)
    auc_attn = roc_auc_score(labels, sc_attn) if np.std(sc_attn) > 0 else 0.5
    results["attention_brain_pooled"] = auc_attn
    print(f"  attention_brain_pooled: {auc_attn:.4f}")
    return results

all_results = {}
for name, edges in edge_sets.items():
    r = evaluate(edges, name)
    if r:
        all_results[name] = r

with open(r"D:\openclaw\biodyn-nmi-paper\dumb_baseline_results.json", "w") as f:
    json.dump(all_results, f, indent=2)
print("\n=== RESULTS SAVED ===")
