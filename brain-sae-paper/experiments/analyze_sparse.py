"""
Analyze sparse SAE features: cell-type AUROC, gene correlations, enrichment tests.
Designed to work with the retrained sparse SAEs (λ >= 1).
"""
import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy import stats
import platform

if platform.system() == "Windows":
    BASE = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments")
else:
    BASE = Path("/mnt/d/openclaw/biodyn-nmi-paper/brain-sae-paper/experiments")

ACT_DIR = BASE / "activations"
SAE_DIR = BASE / "sae_models"
OUT_DIR = BASE / "analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Immune cell marker gene sets (matching actual data)
IMMUNE_MARKERS = {
    "B_cell": ["CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "BANK1", "CD22", "BLK",
               "BLNK", "IGHM", "IGHD", "TCL1A", "FCER2"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "CD2", "CD7", "LCK", "ZAP70", "ITK",
               "TRAC", "TRBC1", "TRBC2", "LEF1", "TCF7"],
    "CD4_T": ["CD4", "IL7R", "CCR7", "FOXP3", "IL2RA", "CTLA4"],
    "CD8_T": ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "NKG7", "IFNG"],
    "macrophage": ["CD68", "CD14", "CD163", "MRC1", "MARCO", "MSR1", "CSF1R",
                   "ITGAM", "FCGR1A", "ADGRE1"],
    "monocyte": ["CD14", "FCGR3A", "S100A9", "S100A8", "LYZ", "VCAN"],
    "neutrophil": ["S100A8", "S100A9", "S100A12", "FCGR3B", "CSF3R", "CXCR2",
                   "MMP9", "ELANE", "MPO"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "KLRB1", "NCAM1", "KLRC1", "NCR1",
                "FCGR3A", "PRF1"],
    "plasma_cell": ["SDC1", "IGHG1", "IGHG2", "MZB1", "JCHAIN", "XBP1",
                    "PRDM1", "IRF4"],
    "erythrocyte": ["HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1", "GYPA"],
}

FUNCTIONAL_SETS = {
    "cytokine_signaling": ["IL1B", "IL6", "TNF", "IFNG", "IL10", "IL2", "IL4",
                           "IL17A", "TGFB1", "CSF2"],
    "antigen_presentation": ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
                             "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
                             "B2M", "TAP1", "TAP2", "CIITA"],
    "TCR_signaling": ["CD3D", "CD3E", "ZAP70", "LCK", "LAT", "SLP76",
                      "PLCG1", "ITK", "FYN", "NFATC1"],
    "BCR_signaling": ["CD79A", "CD79B", "SYK", "BTK", "BLNK", "PLCG2",
                      "PIK3CD", "VAV1", "LYN"],
    "innate_immunity": ["TLR2", "TLR4", "TLR7", "TLR8", "MYD88", "IRAK1",
                        "NFKB1", "NLRP3", "CASP1", "IL1B"],
}


def load_metadata():
    with open(ACT_DIR / "metadata.json") as f:
        return json.load(f)


def compute_auroc(y_true, y_score):
    """Simple AUROC computation without sklearn."""
    pos = y_score[y_true]
    neg = y_score[~y_true]
    if len(pos) == 0 or len(neg) == 0:
        return 0.5
    # Mann-Whitney U
    n_pos, n_neg = len(pos), len(neg)
    u_stat = 0
    for p in pos:
        u_stat += (neg < p).sum() + 0.5 * (neg == p).sum()
    return u_stat / (n_pos * n_neg)


def fisher_exact_enrichment(feature_genes, target_genes, all_genes):
    """Fisher's exact test for enrichment."""
    fs = set(feature_genes)
    ts = set(target_genes)
    ag = set(all_genes)
    a = len(fs & ts)
    b = len(fs - ts)
    c = len(ts - fs)
    d = len(ag - fs - ts)
    if a == 0:
        return 1.0, 0.0
    odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    return pval, odds


def benjamini_hochberg(pvals, alpha=0.05):
    """BH FDR correction."""
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = np.array(pvals)[sorted_idx]
    adjusted = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if i == n - 1:
            adjusted[sorted_idx[i]] = sorted_pvals[i]
        else:
            adjusted[sorted_idx[i]] = min(
                sorted_pvals[i] * n / (i + 1),
                adjusted[sorted_idx[i + 1]]
            )
    return adjusted


def analyze_model(name, metadata):
    """Full analysis of one SAE model."""
    print(f"\nAnalyzing {name}...", flush=True)
    
    # Load cell-level features
    cell_feat_path = SAE_DIR / f"{name}_cell_features.npy"
    if not cell_feat_path.exists():
        print(f"  No cell features found, skipping", flush=True)
        return None
    
    cell_features = np.load(cell_feat_path)
    n_cells, n_features = cell_features.shape
    print(f"  Cell features: {cell_features.shape}", flush=True)
    
    # Load cell types
    cell_types = np.load(ACT_DIR / "cell_types.npy", allow_pickle=True)
    n_use = min(len(cell_types), n_cells)
    cell_types = cell_types[:n_use]
    cell_features = cell_features[:n_use]
    
    # Basic stats
    active_mask = cell_features > 0
    feature_freq = active_mask.mean(axis=0)
    alive_mask = feature_freq > 0.01
    n_alive = alive_mask.sum()
    avg_l0_cell = active_mask.sum(axis=1).mean()
    
    result = {
        "name": name,
        "n_cells": int(n_use),
        "n_features": n_features,
        "n_alive": int(n_alive),
        "avg_L0_cell": float(avg_l0_cell),
    }
    
    # ---- Cell-type AUROC ----
    print(f"  Computing cell-type AUROC...", flush=True)
    unique_types = sorted(set(cell_types))
    ct_auroc = {}
    
    for ct in unique_types:
        ct_mask = cell_types == ct
        if ct_mask.sum() < 5:
            continue
        
        best_auroc = 0
        best_feature = -1
        auroc_list = []
        
        for fi in range(n_features):
            if not alive_mask[fi]:
                continue
            auroc = compute_auroc(ct_mask, cell_features[:, fi])
            auroc_list.append((fi, auroc))
            if auroc > best_auroc:
                best_auroc = auroc
                best_feature = fi
        
        # Also count features above thresholds
        n_above_08 = sum(1 for _, a in auroc_list if a > 0.8)
        n_above_07 = sum(1 for _, a in auroc_list if a > 0.7)
        
        ct_auroc[ct] = {
            "best_auroc": float(best_auroc),
            "best_feature": int(best_feature),
            "n_above_0.8": n_above_08,
            "n_above_0.7": n_above_07,
            "n_cells": int(ct_mask.sum()),
        }
    
    result["cell_type_auroc"] = ct_auroc
    
    # Summary stats
    max_aurocs = [v["best_auroc"] for v in ct_auroc.values()]
    result["mean_best_auroc"] = float(np.mean(max_aurocs)) if max_aurocs else 0
    result["n_types_above_0.8"] = sum(1 for v in ct_auroc.values() if v["best_auroc"] > 0.8)
    
    # ---- Gene correlations (for top features) ----
    print(f"  Computing gene correlations...", flush=True)
    gene_ids_per_cell = metadata.get("gene_ids_per_cell", [])
    id_to_gene = metadata.get("id_to_gene", {})
    all_gene_names = metadata.get("gene_names", [])
    
    # Build gene presence matrix (only for genes in our marker sets)
    marker_genes = set()
    for gs in list(IMMUNE_MARKERS.values()) + list(FUNCTIONAL_SETS.values()):
        marker_genes.update(gs)
    
    # Map gene names to IDs in vocab
    gene_name_to_id = {v: int(k) for k, v in id_to_gene.items()}
    marker_gene_ids = {g: gene_name_to_id[g] for g in marker_genes if g in gene_name_to_id}
    
    # Build gene presence matrix for all genes in each cell
    # For efficiency, compute correlation of features with presence of marker genes
    gene_feature_corr = {}
    
    if gene_ids_per_cell and len(gene_ids_per_cell) >= n_use:
        # For each marker gene, compute presence vector
        for gene_name, gene_id in marker_gene_ids.items():
            presence = np.zeros(n_use, dtype=np.float32)
            for ci in range(n_use):
                if ci < len(gene_ids_per_cell) and gene_id in gene_ids_per_cell[ci]:
                    presence[ci] = 1.0
            
            if presence.sum() < 3 or presence.std() < 1e-8:
                continue
            
            # Correlation with each alive feature
            for fi in range(n_features):
                if not alive_mask[fi] or cell_features[:, fi].std() < 1e-8:
                    continue
                corr = np.corrcoef(cell_features[:, fi], presence)[0, 1]
                if abs(corr) > 0.1:
                    if fi not in gene_feature_corr:
                        gene_feature_corr[fi] = []
                    gene_feature_corr[fi].append({
                        "gene": gene_name,
                        "correlation": float(corr),
                    })
    
    # Sort correlations per feature
    for fi in gene_feature_corr:
        gene_feature_corr[fi].sort(key=lambda x: -abs(x["correlation"]))
    
    result["gene_feature_correlations"] = {
        str(k): v for k, v in gene_feature_corr.items()
    }
    result["n_features_with_gene_corr"] = len(gene_feature_corr)
    
    # ---- Gene set enrichment ----
    print(f"  Running enrichment tests...", flush=True)
    enrichment_results = {}
    all_pvals = []
    all_pval_ids = []
    
    for fi in gene_feature_corr:
        top_genes = [g["gene"] for g in gene_feature_corr[fi] if g["correlation"] > 0.1]
        if len(top_genes) < 2:
            continue
        
        for set_name, gene_set in {**IMMUNE_MARKERS, **FUNCTIONAL_SETS}.items():
            pval, odds = fisher_exact_enrichment(
                top_genes, gene_set,
                list(marker_gene_ids.keys())
            )
            all_pvals.append(pval)
            all_pval_ids.append((fi, set_name, pval, odds))
    
    # FDR correction
    if all_pvals:
        fdr = benjamini_hochberg(all_pvals)
        for idx, (fi, set_name, pval, odds) in enumerate(all_pval_ids):
            if fdr[idx] < 0.05:
                if str(fi) not in enrichment_results:
                    enrichment_results[str(fi)] = {"enrichments": {}, "genes": []}
                    if fi in gene_feature_corr:
                        enrichment_results[str(fi)]["genes"] = [
                            g["gene"] for g in gene_feature_corr[fi][:10]
                        ]
                enrichment_results[str(fi)]["enrichments"][set_name] = {
                    "pval": float(pval),
                    "fdr": float(fdr[idx]),
                    "odds_ratio": float(odds),
                }
    
    result["enrichment_results"] = enrichment_results
    result["n_enriched_features"] = len(enrichment_results)
    
    # Summary: which gene sets are most represented
    set_counts = defaultdict(int)
    for fi_data in enrichment_results.values():
        for sn in fi_data["enrichments"]:
            set_counts[sn] += 1
    result["gene_set_summary"] = dict(sorted(set_counts.items(), key=lambda x: -x[1]))
    
    # ---- PCA comparison ----
    print(f"  Computing PCA comparison...", flush=True)
    # Load raw cell activations for PCA
    parts = name.split("_")
    layer_idx = int(parts[0].replace("layer", ""))
    cell_acts = np.load(ACT_DIR / f"layer_{layer_idx:02d}_cell_activations.npy")
    cell_acts = cell_acts[:n_use]
    
    # PCA via SVD
    cell_acts_centered = cell_acts - cell_acts.mean(axis=0)
    U, S, Vt = np.linalg.svd(cell_acts_centered, full_matrices=False)
    pca_scores = U[:, :50] * S[:50]  # top 50 PCs
    
    pca_auroc = {}
    for ct in unique_types:
        ct_mask = cell_types == ct
        if ct_mask.sum() < 5:
            continue
        best_auroc = 0
        for pc in range(50):
            auroc = compute_auroc(ct_mask, pca_scores[:, pc])
            auroc_rev = compute_auroc(ct_mask, -pca_scores[:, pc])
            auroc = max(auroc, auroc_rev)
            best_auroc = max(best_auroc, auroc)
        pca_auroc[ct] = float(best_auroc)
    
    result["pca_comparison"] = {
        ct: {
            "sae_auroc": ct_auroc.get(ct, {}).get("best_auroc", 0),
            "pca_auroc": pca_auroc.get(ct, 0),
        }
        for ct in unique_types
        if ct in ct_auroc and ct in pca_auroc
    }
    
    # Mean advantage
    advantages = [
        v["sae_auroc"] - v["pca_auroc"]
        for v in result["pca_comparison"].values()
    ]
    result["mean_sae_advantage"] = float(np.mean(advantages)) if advantages else 0
    
    print(f"  Done. Alive={n_alive}, enriched={len(enrichment_results)}, "
          f"mean_best_auroc={result['mean_best_auroc']:.3f}", flush=True)
    
    return result


def main():
    metadata = load_metadata()
    
    # Find sparse SAE models (λ >= 1)
    results_file = SAE_DIR / "training_results_sparse.json"
    if results_file.exists():
        with open(results_file) as f:
            training_results = json.load(f)
    else:
        training_results = {}
    
    # Also check for cell features files directly
    all_results = {}
    
    for cell_feat_file in sorted(SAE_DIR.glob("*_cell_features.npy")):
        name = cell_feat_file.stem.replace("_cell_features", "")
        # Only analyze λ >= 1 models
        if "lam1.0" in name or "lam3.0" in name or "lam10.0" in name:
            result = analyze_model(name, metadata)
            if result:
                # Merge training stats
                if name in training_results:
                    result["training"] = training_results[name]
                all_results[name] = result
    
    # Save
    with open(OUT_DIR / "sparse_analysis_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    # Print summary
    print(f"\n{'='*60}", flush=True)
    print(f"ANALYSIS SUMMARY", flush=True)
    print(f"{'='*60}", flush=True)
    
    for name, r in sorted(all_results.items()):
        tr = r.get("training", {})
        print(f"\n{name}:", flush=True)
        print(f"  L0={tr.get('avg_L0', r.get('avg_L0_cell', '?')):.1f} "
              f"R²={tr.get('R2', '?')} "
              f"alive={r['n_alive']}/{r['n_features']}", flush=True)
        print(f"  Cell-type AUROC: mean_best={r['mean_best_auroc']:.3f}, "
              f"types>0.8: {r['n_types_above_0.8']}", flush=True)
        print(f"  Enriched features: {r['n_enriched_features']}", flush=True)
        if r.get("gene_set_summary"):
            for gs, cnt in list(r["gene_set_summary"].items())[:5]:
                print(f"    {gs}: {cnt} features", flush=True)
        if r.get("pca_comparison"):
            print(f"  SAE vs PCA mean advantage: {r['mean_sae_advantage']:.3f}", flush=True)
    
    print(f"\n=== Analysis Complete ===", flush=True)


if __name__ == "__main__":
    main()
