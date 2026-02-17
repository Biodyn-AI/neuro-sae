# Geneformer GRN Inference — Proper Multi-Model Validation Report

## Overview

We replicated the scGPT attention→GRN pipeline for **Geneformer V1-10M** using the same methodology:
1. Load real DLPFC brain scRNA-seq data (11,202 cells × 32,344 genes)
2. Tokenize via Geneformer's rank-value encoding (19,905/32,344 genes mapped)
3. Run inference, extract attention weights from all 6 layers × 4 heads
4. Average attention across layers/heads, keep top-20 attended targets per gene per cell
5. Aggregate across cells to produce gene-gene edge scores
6. Evaluate against TRRUST (8,427 edges) and DoRothEA (276,562 edges) reference GRNs

## Key Results

| Cells | Total Edges | Time (s) | TRRUST AUROC | DoRothEA AUROC |
|-------|------------|----------|--------------|----------------|
| 200   | 1,564,127  | 8.3      | 0.444        | 0.473          |
| 500   | 3,269,295  | 17.6     | 0.549        | 0.486          |
| 1000  | 5,382,836  | 36.9     | 0.522        | 0.486          |

**All AUROC values are near random (0.5).** Top-K precision/recall/F1 are essentially zero across all thresholds (500, 1000, 5000 edges).

## Comparison with scGPT

From our scGPT analysis (NMI paper):
- **scGPT** attention-based GRN: AUROC ~0.50-0.52 against TRRUST (also near-random)
- **Geneformer**: AUROC 0.44-0.55 against TRRUST (also near-random)

**Both models show the same pattern**: raw attention weights from foundation models do not encode biologically meaningful gene regulatory relationships. This is consistent with the NMI paper's central thesis — that attention weights in these models lack the interpretable structure that mechanistic interpretability claims would require.

## Methodological Details

### Tokenization
- Geneformer uses rank-value encoding: genes ranked by expression normalized to corpus median
- Standard BERT architecture (BertForMaskedLM), 10.3M parameters
- Input: top 512 genes per cell (truncated from median ~2048 for GPU memory)
- 19,905/32,344 DLPFC genes had valid Geneformer token mappings

### Attention Extraction
- All 6 layers × 4 heads averaged to produce per-cell attention matrices
- For each gene position, top-20 most-attended target positions retained
- Edge scores averaged across cells

### Evaluation
- **TRRUST**: 8,427 curated TF→target regulatory edges (human)
- **DoRothEA**: 276,562 edges including chip-seq validated interactions
- Edges filtered to genes present in both predicted and reference sets
- Metrics: AUROC, AUPRC, top-K precision/recall/F1

### Reference gene overlap
- With TRRUST: 27,859–94,896 filtered edges (39–115 true positives), ~0.12-0.14% positive rate
- With DoRothEA: 1.5M–5.2M filtered edges (2,140–7,483 true positives), ~0.14% positive rate

## Interpretation

The near-random performance confirms a fundamental issue: **attention in masked language models reflects statistical co-occurrence patterns in the training corpus, not causal regulatory relationships**. Both scGPT and Geneformer show this pattern despite different:
- Training data (scGPT: 33M cells; Geneformer: 30M cells)
- Input representations (scGPT: binned expression; Geneformer: rank-value)
- Architectures (scGPT: GPT-style decoder; Geneformer: BERT encoder)

This convergent evidence strengthens the paper's argument that attention-based mechanistic interpretability in single-cell foundation models is currently unreliable for GRN inference.

## Files

- `geneformer_grn_pipeline.py` — Full pipeline script
- `geneformer_grn_results.json` — Raw numerical results
- `geneformer_edges_{200,500,1000}cells.csv` — Top 50K predicted edges per cell count
- `stdout3.log` — Full execution log
