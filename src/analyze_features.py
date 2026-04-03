"""
Step 3: Analyze SAE features.

Computes cell-type AUROC, gene set enrichment, PCA comparison, bootstrap CIs,
feature ablation, hyperparameter sensitivity, and layer characterization.

Usage:
    python -m src.analyze_features \
        --activations-dir outputs/activations \
        --sae-dir outputs/sae_models \
        --output-dir results
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

from src.model import SparseAutoencoder


def load_sae(model_path: Path) -> tuple:
    """Load a trained SAE and its normalization stats."""
    ckpt = torch.load(model_path, map_location="cpu", weights_only=False)
    sae = SparseAutoencoder(ckpt["d_input"], ckpt["d_hidden"])
    sae.load_state_dict(ckpt["state_dict"])
    sae.eval()
    return sae, ckpt["acts_mean"], ckpt["acts_std"]


def best_auroc_per_type(features, cell_types):
    """Compute best single-feature AUROC for each cell type."""
    unique_types = sorted(set(ct for ct in cell_types if np.sum(cell_types == ct) >= 5))
    results = {}
    for ct in unique_types:
        labels = (cell_types == ct).astype(int)
        best = 0
        best_j = 0
        for j in range(features.shape[1]):
            if features[:, j].std() < 1e-8:
                continue
            try:
                auroc = roc_auc_score(labels, features[:, j])
                auroc = max(auroc, 1 - auroc)
                if auroc > best:
                    best = auroc
                    best_j = j
            except ValueError:
                pass
        results[ct] = {"auroc": best, "feature_idx": best_j}
    return results


def pca_analysis(cell_acts, cell_types):
    """R1.M4: PCA comparison metrics."""
    print("\n--- PCA comparison ---")
    pca = PCA(n_components=50)
    pca_features = pca.fit_transform(cell_acts)
    unique_types = sorted(set(ct for ct in cell_types if np.sum(cell_types == ct) >= 5))

    aurocs = {}
    for ct in unique_types:
        labels = (cell_types == ct).astype(int)
        best = 0
        for j in range(50):
            try:
                a = roc_auc_score(labels, pca_features[:, j])
                best = max(best, a, 1 - a)
            except ValueError:
                pass
        aurocs[ct] = best

    result = {
        "mean_best_auroc": float(np.mean(list(aurocs.values()))),
        "n_types_above_0.8": sum(1 for v in aurocs.values() if v > 0.8),
        "per_type_auroc": {k: float(v) for k, v in aurocs.items()},
    }
    print(f"  PCA: mean AUROC={result['mean_best_auroc']:.3f}, "
          f"types>0.8={result['n_types_above_0.8']}")
    return result


def bootstrap_ci(features, cell_types, n_bootstrap=1000):
    """R1.M2: Bootstrap confidence intervals for AUROC."""
    print("\n--- Bootstrap CIs ---")
    auroc_data = best_auroc_per_type(features, cell_types)
    results = {}
    for ct, info in auroc_data.items():
        labels = (cell_types == ct).astype(int)
        n_pos = int(labels.sum())
        j = info["feature_idx"]
        boot_aurocs = []
        for _ in range(n_bootstrap):
            idx = np.random.choice(len(labels), len(labels), replace=True)
            if len(set(labels[idx])) < 2:
                continue
            try:
                a = roc_auc_score(labels[idx], features[idx, j])
                boot_aurocs.append(max(a, 1 - a))
            except ValueError:
                pass
        ci_low = np.percentile(boot_aurocs, 2.5) if boot_aurocs else info["auroc"]
        ci_high = np.percentile(boot_aurocs, 97.5) if boot_aurocs else info["auroc"]
        results[ct] = {
            "n": n_pos,
            "best_auroc": info["auroc"],
            "ci_low": float(ci_low),
            "ci_high": float(ci_high),
            "ci_width": float(ci_high - ci_low),
            "well_powered": bool(n_pos >= 100),
        }
        tag = "well-powered" if n_pos >= 100 else "exploratory"
        print(f"  {ct}: n={n_pos}, AUROC={info['auroc']:.3f} "
              f"[{ci_low:.3f}, {ci_high:.3f}] ({tag})")
    return results


def ablation_experiment(sae, cell_normed, features, cell_types):
    """R1.M3: Feature ablation to test causal specificity."""
    print("\n--- Ablation experiment ---")
    results = {}
    for target_ct in ["Macrophage", "B cell"]:
        labels = (cell_types == target_ct).astype(int)
        if labels.sum() < 5:
            continue
        auroc_data = best_auroc_per_type(features, cell_types)
        if target_ct not in auroc_data:
            continue
        best_j = auroc_data[target_ct]["feature_idx"]
        best_auroc = auroc_data[target_ct]["auroc"]

        with torch.no_grad():
            x = torch.tensor(cell_normed, dtype=torch.float32)
            z = sae.encode(x)
            x_recon = sae.decode(z)
            normal_mse = (x - x_recon).pow(2).mean(dim=1).numpy()

            z_ablated = z.clone()
            z_ablated[:, best_j] = 0
            x_ablated = sae.decode(z_ablated)
            ablated_mse = (x - x_ablated).pow(2).mean(dim=1).numpy()

        target_mask = labels == 1
        delta_target = ablated_mse[target_mask].mean() - normal_mse[target_mask].mean()
        delta_other = ablated_mse[~target_mask].mean() - normal_mse[~target_mask].mean()

        # AUROC after ablation
        ablated_features = features.copy()
        ablated_features[:, best_j] = 0
        ablated_auroc_data = best_auroc_per_type(ablated_features, cell_types)
        ablated_auroc = ablated_auroc_data.get(target_ct, {}).get("auroc", 0)

        results[target_ct] = {
            "feature_idx": int(best_j),
            "original_auroc": float(best_auroc),
            "ablated_best_auroc": float(ablated_auroc),
            "delta_mse_target": float(delta_target),
            "delta_mse_other": float(delta_other),
            "selectivity_ratio": float(delta_target / max(delta_other, 1e-10)),
        }
        print(f"  {target_ct}: feature={best_j}, "
              f"selectivity={delta_target / max(delta_other, 1e-10):.1f}x")
    return results


def hyperparameter_sensitivity(sae_dir, acts_dir, cell_types):
    """R2.2: Expansion factor sensitivity."""
    print("\n--- Hyperparameter sensitivity ---")
    cell_acts = np.load(acts_dir / "layer_11_cell_acts.npy")
    results = {}
    for d_hidden in [1024, 2048, 4096]:
        name = f"layer11_d{d_hidden}_lam3.0"
        feat_path = sae_dir / f"{name}_cell_features.npy"
        model_path = sae_dir / f"{name}.pt"
        if not feat_path.exists():
            continue
        features = np.load(feat_path)
        sae, acts_mean, acts_std = load_sae(model_path)
        cell_normed = (cell_acts - acts_mean) / acts_std

        active = features > 0
        with torch.no_grad():
            x = torch.tensor(cell_normed, dtype=torch.float32)
            x_hat, _ = sae(x)
            mse = (x - x_hat).pow(2).mean().item()
        r2 = 1 - mse / np.var(cell_normed)
        aurocs = best_auroc_per_type(features, cell_types)
        mean_auroc = np.mean([v["auroc"] for v in aurocs.values()])

        results[f"expansion_{d_hidden // 512}x"] = {
            "d_hidden": d_hidden,
            "expansion": d_hidden / 512,
            "avg_L0": float(active.sum(axis=1).mean()),
            "alive": int((active.mean(axis=0) > 0.01).sum()),
            "dead": int(d_hidden - (active.mean(axis=0) > 0.01).sum()),
            "R2": float(r2),
            "mean_best_auroc": float(mean_auroc),
        }
        print(f"  {d_hidden//512}x: R2={r2:.3f}, AUROC={mean_auroc:.3f}")
    return results


def layer_characterization(sae_dir, acts_dir, cell_types):
    """R2.3: Per-layer feature analysis."""
    print("\n--- Layer characterization ---")
    results = {}
    for layer_idx in [0, 6, 11]:
        # Try lambda=3 first, then lambda=1
        for lam in [3.0, 1.0]:
            feat_path = sae_dir / f"layer{layer_idx}_d2048_lam{lam}_cell_features.npy"
            if feat_path.exists():
                break
        if not feat_path.exists():
            continue
        features = np.load(feat_path)
        aurocs = best_auroc_per_type(features, cell_types)
        resolved = [ct for ct, v in aurocs.items() if v["auroc"] > 0.8]
        unresolved = [ct for ct, v in aurocs.items() if v["auroc"] <= 0.8]
        mean_auroc = np.mean([v["auroc"] for v in aurocs.values()])

        results[f"layer_{layer_idx}"] = {
            "per_type_auroc": {k: v["auroc"] for k, v in aurocs.items()},
            "resolved_types": resolved,
            "unresolved_types": unresolved,
            "mean_auroc": float(mean_auroc),
        }
        print(f"  Layer {layer_idx}: mean AUROC={mean_auroc:.3f}, "
              f"resolved={len(resolved)}, unresolved={len(unresolved)}")
    return results


def activation_statistics(acts_dir):
    """R1.M1: scGPT activation distribution statistics."""
    print("\n--- Activation statistics (scGPT) ---")
    scgpt_stats = {}
    for layer_idx in [0, 6, 11]:
        token_path = acts_dir / f"layer_{layer_idx}_token_acts.npy"
        if not token_path.exists():
            continue
        acts = np.load(token_path)
        norms = np.linalg.norm(acts, axis=1)
        flat = acts.flatten()
        var_per_dim = np.var(acts, axis=0)
        pr = float(var_per_dim.sum() ** 2 / (var_per_dim**2).sum()) if (var_per_dim**2).sum() > 0 else 0
        scgpt_stats[f"layer_{layer_idx}"] = {
            "mean_norm": float(norms.mean()),
            "kurtosis": float(stats.kurtosis(flat)),
            "fraction_near_zero": float(np.mean(np.abs(flat) < 0.1)),
            "participation_ratio": pr,
        }
        print(f"  scGPT L{layer_idx}: norm={norms.mean():.1f}, PR={pr:.1f}")

    # GPT-2 comparison
    gpt2_stats = {}
    try:
        from transformers import GPT2Model, GPT2Tokenizer

        print("  Loading GPT-2...")
        tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
        gpt2 = GPT2Model.from_pretrained("gpt2", output_hidden_states=True)
        gpt2.eval()
        texts = [
            "The immune system protects against pathogens through innate and adaptive responses.",
            "Single-cell sequencing reveals cellular heterogeneity in complex tissues.",
            "Transformer models learn contextual representations of sequential data.",
            "Gene expression programs define cellular identity and function.",
            "Deep learning has transformed biological data analysis.",
        ] * 20
        all_hidden = {i: [] for i in range(13)}
        for text in texts:
            inputs = tokenizer(text, return_tensors="pt", truncation=True, max_length=128)
            with torch.no_grad():
                outputs = gpt2(**inputs)
            for i, h in enumerate(outputs.hidden_states):
                all_hidden[i].append(h.squeeze(0).numpy())
        for layer_idx in [0, 6, 11]:
            acts = np.concatenate(all_hidden[layer_idx], axis=0)
            if len(acts) > 50000:
                acts = acts[np.random.choice(len(acts), 50000, replace=False)]
            norms = np.linalg.norm(acts, axis=1)
            flat = acts.flatten()
            var_per_dim = np.var(acts, axis=0)
            pr = float(var_per_dim.sum() ** 2 / (var_per_dim**2).sum()) if (var_per_dim**2).sum() > 0 else 0
            gpt2_stats[f"layer_{layer_idx}"] = {
                "mean_norm": float(norms.mean()),
                "kurtosis": float(stats.kurtosis(flat)),
                "fraction_near_zero": float(np.mean(np.abs(flat) < 0.1)),
                "participation_ratio": pr,
            }
            print(f"  GPT-2 L{layer_idx}: norm={norms.mean():.1f}, PR={pr:.1f}")
    except ImportError:
        print("  transformers not installed, skipping GPT-2 comparison")

    return {"scgpt": scgpt_stats, "gpt2": gpt2_stats}


def main(args):
    acts_dir = Path(args.activations_dir)
    sae_dir = Path(args.sae_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cell_types = np.load(acts_dir / "cell_types.npy", allow_pickle=True)
    cell_acts = np.load(acts_dir / "layer_11_cell_acts.npy")

    # Load main SAE (layer 11, lambda=3)
    main_model_path = sae_dir / "layer11_d2048_lam3.0.pt"
    main_feat_path = sae_dir / "layer11_d2048_lam3.0_cell_features.npy"

    if main_model_path.exists() and main_feat_path.exists():
        sae, acts_mean, acts_std = load_sae(main_model_path)
        features = np.load(main_feat_path)
        cell_normed = (cell_acts - acts_mean) / acts_std

        # PCA comparison
        pca_res = pca_analysis(cell_acts, cell_types)
        with open(out_dir / "pca_analysis.json", "w") as f:
            json.dump(pca_res, f, indent=2)

        # Bootstrap CIs
        boot_res = bootstrap_ci(features, cell_types)
        with open(out_dir / "bootstrap_ci.json", "w") as f:
            json.dump(boot_res, f, indent=2)

        # Ablation
        abl_res = ablation_experiment(sae, cell_normed, features, cell_types)
        with open(out_dir / "ablation_results.json", "w") as f:
            json.dump(abl_res, f, indent=2)

    # Hyperparameter sensitivity
    sens_res = hyperparameter_sensitivity(sae_dir, acts_dir, cell_types)
    with open(out_dir / "hyperparameter_sensitivity.json", "w") as f:
        json.dump(sens_res, f, indent=2)

    # Layer characterization
    layer_res = layer_characterization(sae_dir, acts_dir, cell_types)
    with open(out_dir / "layer_characterization.json", "w") as f:
        json.dump(layer_res, f, indent=2)

    # Activation statistics
    act_res = activation_statistics(acts_dir)
    with open(out_dir / "activation_statistics.json", "w") as f:
        json.dump(act_res, f, indent=2)

    print(f"\nAll results saved to {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze SAE features")
    parser.add_argument("--activations-dir", default="outputs/activations")
    parser.add_argument("--sae-dir", default="outputs/sae_models")
    parser.add_argument("--output-dir", default="results")
    main(parser.parse_args())
