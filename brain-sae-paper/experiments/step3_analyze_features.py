"""
Step 3: Analyze SAE features — gene enrichment, cell type mapping, GO terms.
"""
import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy import stats

ACT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\activations")
SAE_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\sae_models")
OUT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\analysis")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- Known neuroscience gene sets ----------
CELL_TYPE_MARKERS = {
    "excitatory_neurons": [
        "SLC17A7", "SATB2", "CUX2", "RORB", "TLE4", "FEZF2", "THEMIS",
        "SLC17A6", "CAMK2A", "GRIN1", "GRIA1", "GRIA2", "NRGN", "VGLUT1",
    ],
    "inhibitory_neurons": [
        "GAD1", "GAD2", "SLC32A1", "SST", "PVALB", "VIP", "LAMP5",
        "ADARB2", "LHX6", "RELN", "CALB1", "CALB2", "NPY", "CCK",
    ],
    "astrocytes": [
        "GFAP", "AQP4", "SLC1A2", "SLC1A3", "ALDH1L1", "S100B", "GJA1",
        "NDRG2", "SOX9", "ID3", "FABP7", "GLUL",
    ],
    "oligodendrocytes": [
        "MBP", "MOG", "OLIG1", "OLIG2", "PLP1", "MAG", "CNP",
        "CLDN11", "SOX10", "MOBP", "ERMN",
    ],
    "microglia": [
        "CX3CR1", "P2RY12", "TMEM119", "CSF1R", "AIF1", "ITGAM",
        "CD68", "TREM2", "TYROBP", "C1QA", "C1QB", "C1QC",
    ],
    "OPCs": [
        "PDGFRA", "CSPG4", "VCAN", "GPR17", "OLIG2",
    ],
    "endothelial": [
        "CLDN5", "FLT1", "PECAM1", "VWF", "CDH5",
    ],
}

NEURO_GWAS_GENES = {
    "intelligence_EA": [
        "FOXP2", "BDNF", "NRXN1", "SHANK3", "DISC1", "DYRK1A", "MECP2",
        "FMR1", "CNTNAP2", "NRXN3", "NLGN4X", "DLG4", "GRIN2A", "GRIN2B",
        "CACNA1C", "TCF4", "ZNF804A", "ANK3", "RBFOX1", "SYNGAP1",
        "PTEN", "TSC1", "TSC2", "UBE3A", "AUTS2", "RERE",
    ],
    "synaptic_transmission": [
        "SYN1", "SYN2", "SYP", "SNAP25", "STX1A", "VAMP2", "SYT1",
        "CPLX1", "CPLX2", "STXBP1", "NSF", "RIMS1",
    ],
    "myelination": [
        "MBP", "PLP1", "MOG", "MAG", "CNP", "OLIG1", "OLIG2",
        "NKX2-2", "SOX10", "MYRF",
    ],
    "neuroinflammation": [
        "TNF", "IL1B", "IL6", "TGFB1", "CCL2", "CXCL10",
        "NFKB1", "STAT3", "TLR4", "NLRP3",
    ],
}


def compute_feature_gene_associations(feature_acts, gene_ids_per_cell, id_to_gene, n_top=20):
    """For cell-level features, find which genes are associated with each feature."""
    n_features = feature_acts.shape[1]
    n_cells = feature_acts.shape[0]
    
    # Build gene expression matrix: which genes are present in which cells
    all_genes = set()
    for ids in gene_ids_per_cell:
        for gid in ids:
            all_genes.add(gid)
    
    gene_list = sorted(all_genes)
    gene_idx_map = {g: i for i, g in enumerate(gene_list)}
    gene_presence = np.zeros((n_cells, len(gene_list)), dtype=np.float32)
    
    for ci, ids in enumerate(gene_ids_per_cell):
        if ci >= n_cells:
            break
        for gid in ids:
            gene_presence[ci, gene_idx_map[gid]] = 1.0
    
    # For each feature, compute correlation with each gene's presence
    results = {}
    for fi in range(n_features):
        f_vals = feature_acts[:, fi]
        if f_vals.std() < 1e-8:
            continue
        
        # Correlation with gene presence
        corrs = np.array([
            np.corrcoef(f_vals, gene_presence[:, gi])[0, 1]
            if gene_presence[:, gi].std() > 1e-8 else 0.0
            for gi in range(len(gene_list))
        ])
        
        # Top positive and negative
        top_pos = np.argsort(corrs)[-n_top:][::-1]
        top_neg = np.argsort(corrs)[:n_top]
        
        results[fi] = {
            "top_positive_genes": [
                {"gene": id_to_gene.get(str(gene_list[gi]), f"id_{gene_list[gi]}"),
                 "correlation": float(corrs[gi])}
                for gi in top_pos if abs(corrs[gi]) > 0.05
            ],
            "top_negative_genes": [
                {"gene": id_to_gene.get(str(gene_list[gi]), f"id_{gene_list[gi]}"),
                 "correlation": float(corrs[gi])}
                for gi in top_neg if abs(corrs[gi]) > 0.05
            ],
        }
    
    return results


def enrichment_test(feature_genes, gene_set, all_genes_in_data):
    """Fisher's exact test for enrichment of gene_set in feature_genes."""
    feature_set = set(feature_genes)
    target_set = set(gene_set)
    all_set = set(all_genes_in_data)
    
    a = len(feature_set & target_set)  # in both
    b = len(feature_set - target_set)  # in feature only
    c = len(target_set - feature_set)  # in target only
    d = len(all_set - feature_set - target_set)  # in neither
    
    if a == 0:
        return 1.0, 0.0
    
    odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    return pval, odds


def analyze_sae_model(name, feature_path, layer_idx, metadata):
    """Full analysis of one SAE model."""
    print(f"\nAnalyzing {name}...")
    features = np.load(feature_path)
    n_samples, n_features = features.shape
    
    # Use cell-level activations for gene associations
    cell_acts_path = ACT_DIR / f"layer_{layer_idx:02d}_cell_activations.npy"
    
    # Determine if features are token-level or cell-level
    n_cells = metadata["n_cells"]
    
    # Load cell-level features by averaging token features per cell
    # For now use the feature activations directly
    # If token-level, we need to aggregate to cell level
    gene_ids_per_cell = metadata.get("gene_ids_per_cell", [])
    id_to_gene = metadata.get("id_to_gene", {})
    all_gene_names = metadata.get("gene_names", [])
    
    result = {
        "name": name,
        "layer": layer_idx,
        "n_features": n_features,
        "n_samples": n_samples,
    }
    
    # Feature activity statistics
    active_mask = features > 0
    feature_freq = active_mask.mean(axis=0)  # how often each feature fires
    alive_features = (feature_freq > 0.01).sum()
    dead_features = n_features - alive_features
    avg_l0 = active_mask.sum(axis=1).mean()
    
    result["alive_features"] = int(alive_features)
    result["dead_features"] = int(dead_features)
    result["avg_L0"] = float(avg_l0)
    result["feature_freq_stats"] = {
        "mean": float(feature_freq.mean()),
        "median": float(np.median(feature_freq)),
        "max": float(feature_freq.max()),
        "min": float(feature_freq[feature_freq > 0].min()) if (feature_freq > 0).any() else 0,
    }
    
    # For cell-level analysis, aggregate token features to cell level if needed
    # We'll use the feature activations on cell-level activations
    cell_features_path = SAE_DIR / f"{name}_cell_features.npy"
    
    if n_samples != n_cells and cell_acts_path.exists():
        print(f"  Recomputing cell-level features...")
        import torch
        sae_ckpt = torch.load(SAE_DIR / f"{name}.pt", map_location="cpu")
        from step2_train_saes import SparseAutoencoder
        sae = SparseAutoencoder(sae_ckpt["d_input"], sae_ckpt["d_hidden"])
        sae.load_state_dict(sae_ckpt["state_dict"])
        sae.eval()
        
        cell_acts = np.load(cell_acts_path)
        mean = np.load(SAE_DIR / f"layer_{layer_idx:02d}_mean.npy")
        std = np.load(SAE_DIR / f"layer_{layer_idx:02d}_std.npy")
        cell_acts_normed = (cell_acts - mean) / std
        
        with torch.no_grad():
            cell_features = sae.encode(torch.tensor(cell_acts_normed, dtype=torch.float32)).numpy()
        np.save(cell_features_path, cell_features)
    elif n_samples == n_cells:
        cell_features = features
    else:
        cell_features = features[:n_cells] if n_samples >= n_cells else features
    
    # Gene-feature associations (using cell-level)
    if gene_ids_per_cell and len(gene_ids_per_cell) >= min(n_cells, cell_features.shape[0]):
        print(f"  Computing gene-feature associations...")
        n_use = min(cell_features.shape[0], len(gene_ids_per_cell))
        gene_assoc = compute_feature_gene_associations(
            cell_features[:n_use], gene_ids_per_cell[:n_use], id_to_gene, n_top=15
        )
        
        # For each alive feature, test enrichment against known gene sets
        print(f"  Running enrichment tests...")
        enrichment_results = {}
        
        for fi in sorted(gene_assoc.keys()):
            if feature_freq[fi] < 0.01:
                continue
            
            top_genes = [g["gene"] for g in gene_assoc[fi]["top_positive_genes"]]
            
            fi_enrichments = {}
            for set_name, gene_set in {**CELL_TYPE_MARKERS, **NEURO_GWAS_GENES}.items():
                pval, odds = enrichment_test(top_genes, gene_set, all_gene_names)
                if pval < 0.05:
                    fi_enrichments[set_name] = {"pval": pval, "odds_ratio": odds}
            
            if fi_enrichments:
                enrichment_results[str(fi)] = {
                    "top_genes": top_genes[:10],
                    "enrichments": fi_enrichments,
                }
        
        result["n_enriched_features"] = len(enrichment_results)
        result["enrichment_results"] = enrichment_results
        
        # Top-level summary: which gene sets are most represented?
        set_counts = defaultdict(int)
        for fi_data in enrichment_results.values():
            for set_name in fi_data["enrichments"]:
                set_counts[set_name] += 1
        result["gene_set_summary"] = dict(sorted(set_counts.items(), key=lambda x: -x[1]))
    
    # Cell type mapping (if cell types available)
    cell_types_path = ACT_DIR / "cell_types.npy"
    if cell_types_path.exists():
        print(f"  Mapping features to cell types...")
        cell_types = np.load(cell_types_path, allow_pickle=True)
        n_use = min(len(cell_types), cell_features.shape[0])
        cell_types = cell_types[:n_use]
        cell_feats = cell_features[:n_use]
        
        unique_types = list(set(cell_types))
        type_feature_map = {}
        
        for ct in unique_types:
            ct_mask = cell_types == ct
            if ct_mask.sum() < 5:
                continue
            
            # For each feature, test if it's more active in this cell type
            ct_scores = []
            for fi in range(n_features):
                if feature_freq[fi] < 0.01:
                    continue
                in_ct = cell_feats[ct_mask, fi]
                out_ct = cell_feats[~ct_mask, fi]
                if in_ct.std() < 1e-8 and out_ct.std() < 1e-8:
                    continue
                try:
                    t_stat, pval = stats.ttest_ind(in_ct, out_ct, alternative='greater')
                    if pval < 0.01:
                        ct_scores.append({"feature": fi, "t_stat": float(t_stat), "pval": float(pval)})
                except:
                    pass
            
            ct_scores.sort(key=lambda x: -x["t_stat"])
            type_feature_map[ct] = ct_scores[:20]
        
        result["cell_type_features"] = type_feature_map
    
    return result


def main():
    # Load metadata
    with open(ACT_DIR / "metadata.json") as f:
        metadata = json.load(f)
    
    # Find all SAE models
    all_results = {}
    
    for pt_file in sorted(SAE_DIR.glob("layer*_d*_lam*.pt")):
        name = pt_file.stem
        feature_path = SAE_DIR / f"{name}_features.npy"
        if not feature_path.exists():
            print(f"Skipping {name}: no features file")
            continue
        
        # Parse name
        parts = name.split("_")
        layer_idx = int(parts[0].replace("layer", ""))
        
        result = analyze_sae_model(name, feature_path, layer_idx, metadata)
        all_results[name] = result
    
    # Save all results
    with open(OUT_DIR / "all_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    # Generate summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    summary_lines = []
    for name, r in sorted(all_results.items()):
        line = (f"{name}: alive={r['alive_features']}/{r['n_features']}, "
                f"L0={r['avg_L0']:.1f}, enriched={r.get('n_enriched_features', 0)}")
        print(line)
        summary_lines.append(line)
        
        if r.get("gene_set_summary"):
            for gs, count in list(r["gene_set_summary"].items())[:5]:
                print(f"  {gs}: {count} features")
    
    with open(OUT_DIR / "summary.txt", "w") as f:
        f.write("\n".join(summary_lines))
    
    print(f"\n=== Step 3 Complete ===")
    print(f"Results saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
