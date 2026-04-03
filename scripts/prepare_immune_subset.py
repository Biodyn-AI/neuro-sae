"""
Download Tabula Sapiens blood + bone marrow data from CZ CELLxGENE and
create the 1,000-cell immune subset used in the paper.

Requires: pip install requests scanpy anndata
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import requests
import scanpy as sc

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "blood": {
        "id": "b225ee37-5e06-4e49-9c25-c3d7b5008dab",
        "file": "tabula_sapiens_blood.h5ad",
    },
    "bone_marrow": {
        "id": "c7f0c3ea-2083-4d87-a8e0-7f69626aa40d",
        "file": "tabula_sapiens_bone_marrow.h5ad",
    },
}

# Cell-type mapping from CZ CELLxGENE ontology names to paper labels
TYPE_MAP = {
    "B cell": "B cell",
    "CD4-positive, alpha-beta T cell": "CD4+ T cell",
    "naive thymus-derived CD4-positive, alpha-beta T cell": "Naive CD4+ T cell",
    "CD8-positive, alpha-beta T cell": "CD8+ T cell",
    "macrophage": "Macrophage",
    "neutrophil": "Neutrophil",
    "plasma cell": "Plasma cell",
    "monocyte": "Monocyte",
    "classical monocyte": "Classical monocyte",
    "natural killer cell": "NK cell",
    "erythrocyte": "Erythrocyte",
    "non-classical monocyte": "Non-classical monocyte",
    "intermediate monocyte": "Intermediate monocyte",
    "mature NK T cell": "NK T cell",
    "hematopoietic precursor cell": "Hematopoietic precursor",
    "hematopoietic stem cell": "Hematopoietic stem cell",
    "basophil": "Basophil",
    "regulatory T cell": "Regulatory T cell",
    "plasmacytoid dendritic cell": "Plasmacytoid DC",
    "myeloid dendritic cell": "Myeloid DC",
    "granulocyte": "Granulocyte",
    "erythroid progenitor cell": "Erythroid progenitor",
    "common myeloid progenitor": "Common myeloid progenitor",
    "T cell": "T cell",
    "platelet": "Platelet",
}

# Target composition matching the paper's Table 2
TARGET_COUNTS = {
    "B cell": 189,
    "CD4+ T cell": 180,
    "Macrophage": 127,
    "Neutrophil": 116,
    "CD8+ T cell": 112,
    "Plasma cell": 44,
    "Monocyte": 42,
    "NK cell": 28,
    "Erythrocyte": 30,
    "Classical monocyte": 25,
    "Naive CD4+ T cell": 20,
}


def get_download_url(dataset_id: str, asset_type: str = "H5AD") -> str:
    """Get CZ CELLxGENE download URL for a dataset."""
    collection_url = (
        "https://api.cellxgene.cziscience.com/dp/v1/collections/"
        "e5f58829-1a66-40b5-a624-9046778e74f5"
    )
    resp = requests.get(collection_url)
    resp.raise_for_status()
    for ds in resp.json().get("datasets", []):
        ds_id = ds.get("dataset_id", ds.get("id", ""))
        if ds_id == dataset_id:
            for asset in ds.get("dataset_assets", []):
                if asset["filetype"] == asset_type:
                    dl = requests.get(
                        f"https://api.cellxgene.cziscience.com/dp/v1/"
                        f"datasets/{ds_id}/asset/{asset['id']}"
                    )
                    dl.raise_for_status()
                    return dl.json()["url"]
    raise ValueError(f"Dataset {dataset_id} not found")


def download_dataset(name: str, info: dict) -> Path:
    """Download a single tissue dataset if not already present."""
    path = DATA_DIR / info["file"]
    if path.exists():
        print(f"  {name}: already downloaded")
        return path
    url = get_download_url(info["id"])
    print(f"  {name}: downloading from CZ CELLxGENE...")
    import subprocess
    subprocess.run(["curl", "-L", "-o", str(path), url], check=True)
    return path


def main():
    print("Downloading Tabula Sapiens tissue datasets...")
    paths = {}
    for name, info in DATASETS.items():
        paths[name] = download_dataset(name, info)

    print("\nBuilding 1,000-cell immune subset...")
    datasets = []
    for name, path in paths.items():
        adata = sc.read_h5ad(path)
        adata.var_names = pd.Index(adata.var["feature_name"].astype(str).values)
        adata.var_names_make_unique()
        datasets.append(adata)

    combined = ad.concat(datasets, join="outer")
    combined.obs["cell_type_mapped"] = (
        combined.obs["cell_type"].map(TYPE_MAP).fillna(combined.obs["cell_type"])
    )

    rng = np.random.RandomState(42)
    indices = []
    for ct, n in TARGET_COUNTS.items():
        available = np.where(combined.obs["cell_type_mapped"] == ct)[0]
        indices.extend(rng.choice(available, min(n, len(available)), replace=False))

    remaining = set(combined.obs["cell_type_mapped"].unique()) - set(TARGET_COUNTS)
    other = []
    for ct in sorted(remaining):
        available = np.where(combined.obs["cell_type_mapped"] == ct)[0]
        other.extend(rng.choice(available, min(15, len(available)), replace=False))
    rng.shuffle(other)
    indices.extend(other[:87])

    subset = combined[indices].copy()
    out_path = DATA_DIR / "tabula_sapiens_immune_1000.h5ad"
    subset.write_h5ad(out_path)
    print(f"Saved {subset.n_obs} cells to {out_path}")
    print(f"Cell types: {subset.obs['cell_type_mapped'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
