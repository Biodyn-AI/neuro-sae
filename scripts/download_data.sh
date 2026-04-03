#!/usr/bin/env bash
# Download all external data and model checkpoints.
#
# Usage:
#   bash scripts/download_data.sh
#
# Requirements: git, pip install gdown cellxgene-census scanpy anndata

set -euo pipefail
cd "$(dirname "$0")/.."

echo "=== 1. Cloning scGPT repository ==="
if [ -d "external/scGPT" ]; then
    echo "  external/scGPT already exists, skipping."
else
    mkdir -p external
    git clone https://github.com/bowang-lab/scGPT.git external/scGPT
fi

echo ""
echo "=== 2. Downloading scGPT brain checkpoint ==="
if [ -f "external/scGPT_checkpoints/brain/best_model.pt" ]; then
    echo "  Checkpoint already exists, skipping."
else
    mkdir -p external/scGPT_checkpoints/brain
    pip install -q gdown
    gdown --folder "https://drive.google.com/drive/folders/1vf1ijfQSk7rGdDGpBntR5bi5g6gNt-Gx" \
        -O external/scGPT_checkpoints/brain/
    # gdown puts files in a subfolder; move them up if needed
    if [ -d "external/scGPT_checkpoints/brain/scGPT_brain" ]; then
        mv external/scGPT_checkpoints/brain/scGPT_brain/* external/scGPT_checkpoints/brain/
        rmdir external/scGPT_checkpoints/brain/scGPT_brain
    fi
fi

echo ""
echo "=== 3. Downloading Tabula Sapiens immune data ==="
if [ -f "data/tabula_sapiens_immune_1000.h5ad" ]; then
    echo "  data/tabula_sapiens_immune_1000.h5ad already exists, skipping."
else
    mkdir -p data
    python scripts/prepare_immune_subset.py
fi

echo ""
echo "=== Done ==="
echo "Files:"
echo "  external/scGPT/                          -- scGPT source code"
echo "  external/scGPT_checkpoints/brain/         -- model checkpoint"
echo "  data/tabula_sapiens_immune_1000.h5ad      -- 1000-cell immune subset"
