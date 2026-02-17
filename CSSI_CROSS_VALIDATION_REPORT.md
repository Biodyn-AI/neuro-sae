# CSSI Cross-Dataset Validation Report (scGPT Brain)

## Experiment
- **Model:** scGPT brain checkpoint (12 layers, 8 heads, d=512)
- **Dataset:** DLPFC 11k brain scRNA-seq
- **Cells:** 497 total -> Half A (248), Half B (249)
- **Genes (Half A):** 1455
- **Genes (Half B):** 1458
- **Ground truth:** TRRUST (8427 directed edges)
- **Metric:** AUROC (attention score vs TRRUST edge label)

## Per-Layer AUROC

| Layer | Half A | Half B |
|-------|--------|--------|
|  0 | 0.6208 | 0.6467 |
|  1 | 0.5966 | 0.6043 |
|  2 | 0.4608 | 0.4465 |
|  3 | 0.5382 | 0.5086 |
|  4 | 0.6100 | 0.5664 |
|  5 | 0.5434 | 0.4888 |
|  6 | 0.5315 | 0.4723 |
|  7 | 0.4761 | 0.4250 |
|  8 | 0.5022 | 0.4383 |
|  9 | 0.5582 | 0.5462 |
| 10 | 0.6079 | 0.6440 |
| 11 | 0.5941 | 0.6205 |

## Cross-Validation

- **Best layers on Half A:** [0, 4, 10]
- **Best layers on Half B:** [0, 10, 11]
- **Layer overlap:** 2/3

### Layers selected on A, evaluated on B
| Layer | AUROC (A, train) | AUROC (B, test) |
|-------|------------------|-----------------|
| 0 | 0.6208 | 0.6467412830075172 |
| 4 | 0.6100 | 0.5663731779809281 |
| 10 | 0.6079 | 0.643999786385038 |

### Layers selected on B, evaluated on A
| Layer | AUROC (B, train) | AUROC (A, test) |
|-------|------------------|-----------------|
| 0 | 0.6467 | 0.6207683698046115 |
| 10 | 0.6440 | 0.6079084804427164 |
| 11 | 0.6205 | 0.5940514095809466 |

### Comparison to Random Baseline
- Mean AUROC of best-A layers on B: **0.6190**
- Random 3-layer baseline on B: **0.5392** +/- 0.0397
- Improvement: **0.0799** (2.0 sigma)

## Conclusion
Cross-dataset validation **PASSED**. Layers identified on one half of the data
generalize to the held-out half, with 2/3 layer overlap
and 0.0799 AUROC improvement over random layer selection
(2.0 sigma).

*Generated 2026-02-14 20:19:14 | Runtime: 840s*