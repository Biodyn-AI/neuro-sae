# neuro-sae

This repository is a research workspace for mechanistic-interpretability experiments on single-cell foundation models (mainly scGPT and Geneformer), plus manuscript assets and generated reports.

It is not a packaged Python library. Most workflows are script-driven, and several scripts rely on local absolute paths from the original development environment.

## Scope

The code and artifacts here cover several related tracks:

- Attention-based gene regulatory network (GRN) inference and evaluation
- Scaling behavior characterization across cell counts
- Baseline method comparisons (correlation, MI, GENIE3/GRNBoost2)
- CSSI (Cell-State Stratified Interpretability) cross-dataset and cross-tissue validation
- CRISPRi perturbation validation pipelines
- Synthetic validation studies with known ground truth
- Paper drafts, sections, and generated figures/reports

## Repository Layout

- `analysis/`: validation and baseline analysis scripts (including dependency smoke tests)
- `baseline_comparison/`: GRN baseline scripts and summary outputs
- `brain-sae-paper/`: SAE-focused manuscript and experiment scripts
- `cssi_real/`: CSSI validation scripts and outputs
- `experiments/`: scaling and perturbation validation experiments
- `figures/`: generated figure assets
- `multi_model/`, `multi_model_extended/`, `multi_model_proper/`: multi-model/Geneformer analysis scripts and reports
- `novel_method/`: CSSI method experiments and figures
- `paper_docs/root_paper/`: unified NMI manuscript source
- `reports/`: generated markdown reports summarizing major analyses
- `results/`: legacy result snapshots
- `synthetic_validation/`: synthetic data validation pipeline

## Environment Setup

Python 3.10+ is recommended.

Install core dependencies:

```bash
python -m pip install \
  numpy pandas scipy scikit-learn matplotlib seaborn \
  scanpy anndata h5py torch transformers statsmodels \
  networkx requests
```

Optional dependencies used by specific scripts:

```bash
python -m pip install arboreto geneformer
```

Quick dependency check:

```bash
python analysis/tests/test_deps.py
```

## Quick Start (Self-Contained / Low-Friction)

These scripts can run without external model checkpoints:

### 1) Synthetic validation pipeline

```bash
python synthetic_validation/synthetic_validation.py
python synthetic_validation/reproduce_analysis.py --analysis scaling
python synthetic_validation/reproduce_analysis.py --analysis sweep
```

### 2) Perturbation validation demo (simulated)

```bash
python experiments/perturbation_validation/demo_perturbation_validation.py
```

Note: this demo is intentionally simulated and writes outputs under a hardcoded path in the script.

## Data-Dependent Pipelines

Most core pipelines require local datasets/checkpoints and often hardcoded paths. Before running, inspect and update path constants at the top of each script.

### Scaling characterization (Geneformer attention)

```bash
python experiments/scaling_characterization/scaling_characterization.py
```

Expected resources include:

- Tabula Sapiens immune `.h5ad`
- Geneformer model + token dictionaries
- TRRUST reference network

### CSSI validation

```bash
python cssi_real/run_cssi_validation.py --both
# or
python cssi_real/cssi_cross_dataset_validation.py
python cssi_real/cssi_cross_tissue_validation.py --source-tissue immune --target-tissue kidney
```

### Baseline comparison

```bash
python baseline_comparison/fast_baseline.py
python baseline_comparison/baseline_comparison.py
```

### Saturation analysis

```bash
python analysis/validation/saturation_analysis.py \
  --data_path /path/to/raw_data_dir \
  --output_dir /path/to/output_dir
```

## Manuscripts

Two main manuscript roots are present:

- `brain-sae-paper/main.tex`
- `paper_docs/root_paper/main.tex`

Compile from the manuscript directory, for example:

```bash
cd paper_docs/root_paper
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Note: the top-level `Makefile` currently expects `main.tex` at repo root and may not match active manuscript locations.

## Key Existing Reports

- `reports/SCALING_CHARACTERIZATION_REPORT.md`
- `reports/SATURATION_ANALYSIS_REPORT.md`
- `reports/PERTURBATION_VALIDATION_REPORT.md`
- `reports/CSSI_CROSS_VALIDATION_REPORT.md`
- `multi_model_proper/REPORT.md`

These files summarize already-run analyses and are a good starting point for understanding current findings.

## Practical Notes

- Several scripts assume WSL + CUDA setup.
- Many scripts use absolute Windows/WSL paths (for example `/mnt/d/...` or `D:\...`).
- This repository contains exploratory and report-generation code; not all scripts are productionized or parameterized.

If you want to make this repo fully reproducible on a new machine, the first step is to centralize path/config handling (for example, via CLI args or a shared YAML config) and remove hardcoded filesystem assumptions.
