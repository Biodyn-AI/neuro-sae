# Scaling Characterization Report: Geneformer Attention-Based GRN Recovery

## Experiment Overview

This experiment addresses the reviewer concern that "only four cell-count points for the core scaling claim" was insufficient to characterize the relationship between cell count and GRN recovery performance. We tested Geneformer's attention-based gene regulatory network (GRN) recovery at **9 cell count points** spanning 25 to 1,000 cells, each with 2 independent repeats and 50-iteration bootstrap confidence intervals.

### Configuration

| Parameter | Value |
|---|---|
| Dataset | Tabula Sapiens immune subset (20,000 cells x 34,293 genes post-filter) |
| Cell counts tested | 25, 50, 100, 150, 200, 300, 500, 750, 1000 |
| Attention layer | Layer 13 (pooled across all 18 heads) |
| Gene vocabulary | Top 1,000 most frequent genes |
| Ground truth | TRRUST (8,330 mapped edges) |
| Positive edges evaluated | 161 per condition |
| Repeats per cell count | 2 (stratified sampling by cell type) |
| Bootstrap iterations | 50 per repeat |
| Sampling | Stratified by cell type, proportional to dataset distribution |

---

## Results Table

| Cell Count | Mean AUROC | Std Dev | 95% CI Lower | 95% CI Upper | CI Width |
|---:|---:|---:|---:|---:|---:|
| 25 | 0.5359 | 0.0091 | 0.5144 | 0.5648 | 0.0504 |
| 50 | 0.5419 | 0.0062 | 0.5172 | 0.5813 | 0.0641 |
| 100 | 0.5445 | 0.0011 | 0.5182 | 0.5942 | 0.0760 |
| 150 | 0.5557 | 0.0026 | 0.5209 | 0.6053 | 0.0844 |
| 200 | 0.5924 | 0.0085 | 0.5504 | 0.6336 | 0.0832 |
| 300 | 0.5936 | 0.0214 | 0.5564 | 0.6398 | 0.0834 |
| 500 | 0.6072 | 0.0046 | 0.5593 | 0.6600 | 0.1007 |
| 750 | 0.6106 | 0.0066 | 0.5635 | 0.6551 | 0.0916 |
| 1000 | 0.5960 | 0.0028 | 0.5391 | 0.6389 | 0.0998 |

---

## Curve Fitting Analysis

Three scaling models were fit to the 9 data points. Power-law fitting did not converge (maxfev exceeded).

### Model Comparison

| Model | R² | Equation |
|---|---:|---|
| **Exponential saturation** | **0.8998** | `AUROC = 0.0888 * (1 - exp(-0.00554 * cells)) + 0.5185` |
| Logarithmic | 0.8313 | `AUROC = 0.0224 * ln(cells) + 0.4572` |
| Linear | 0.5910 | `AUROC = 0.0000683 * cells + 0.5520` |

**Best-fit model: Exponential saturation (R² = 0.90)**

The exponential saturation model implies:
- **Baseline AUROC** (zero-cell extrapolation): ~0.518
- **Asymptotic ceiling**: ~0.607 (baseline + amplitude = 0.518 + 0.089)
- **Characteristic scale**: ~181 cells (1/rate = 1/0.00554), meaning ~63% of the maximum improvement is achieved by ~181 cells

### Scaling Pattern Assessment

The script's pattern analysis classified the scaling as **logarithmic** based on the decreasing rate of marginal improvement (change trend = -0.0036).

Consecutive AUROC differences:

| Interval | Delta AUROC | Relative Change |
|---|---:|---:|
| 25 -> 50 | +0.0061 | +1.1% |
| 50 -> 100 | +0.0025 | +0.5% |
| 100 -> 150 | +0.0112 | +2.1% |
| 150 -> 200 | +0.0367 | +6.6% |
| 200 -> 300 | +0.0012 | +0.2% |
| 300 -> 500 | +0.0136 | +2.3% |
| 500 -> 750 | +0.0034 | +0.6% |
| 750 -> 1000 | **-0.0146** | **-2.4%** |

---

## Key Findings

### 1. Diminishing returns with clear saturation

AUROC rises from 0.536 at 25 cells to a peak of ~0.611 at 750 cells, then **declines slightly** to 0.596 at 1,000 cells. The total dynamic range is modest (~0.075 AUROC units), and most of the gain occurs by 200 cells. The exponential saturation model captures this pattern well (R² = 0.90).

### 2. Inflection around 150-200 cells

The largest single jump occurs between 150 and 200 cells (+0.037 AUROC, +6.6% relative). Below 150 cells, performance hovers near chance (0.536-0.556). Above 200 cells, gains are incremental. This suggests ~200 cells is a critical threshold for the attention-based GRN recovery method.

### 3. Performance plateau and possible overfitting at high cell counts

Performance peaks at 750 cells (0.611) then drops at 1,000 cells (0.596). This reversal suggests either:
- Noise dilution from heterogeneous cell types at high counts
- Attention signal averaging that washes out cell-type-specific regulatory patterns
- Stochastic variation (the 300-cell point also shows high inter-repeat variance: std = 0.021)

### 4. Confidence intervals overlap substantially

All 95% bootstrap CIs overlap across conditions, indicating that while the trend is consistent, the per-condition differences are not individually statistically significant. The overall trend from 25 to 750 cells is clear, but adjacent points cannot be reliably distinguished.

### 5. Not a simple linear scaling

The linear model is a poor fit (R² = 0.59). The relationship is sublinear, better described by logarithmic (R² = 0.83) or exponential saturation (R² = 0.90) curves. This confirms that **doubling cell count does not proportionally improve GRN recovery** -- there are strongly diminishing returns.

---

## Implications for the NMI Paper

1. **The 9-point scaling curve replaces the original 4-point curve**, providing substantially more granular evidence for the scaling claim.
2. **The curve shape is sublinear / saturating**, not linear. The exponential saturation model (R² = 0.90) is the best descriptor, with a characteristic scale of ~181 cells.
3. **Performance ceiling is modest** (~0.61 AUROC), suggesting that scaling cell count alone is insufficient for strong GRN recovery from attention weights.
4. **200 cells is a practical minimum** for meaningful above-chance performance; below this threshold, AUROC stays in the 0.53-0.56 range.
5. **The 750 -> 1000 cell decline** warrants discussion: it may reflect the increasing dominance of common cell types in larger samples, diluting regulatory signals from rarer populations.

---

## Artifacts

- Full results JSON: `experiments/scaling_characterization/scaling_characterization_results.json`
- Experiment log: `experiments/scaling_characterization/scaling_characterization.log`
- Script: `experiments/scaling_characterization/scaling_characterization.py`
