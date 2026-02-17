# Codex Research Quality Review

1. Scaling behavior: **Issue** — benchmark dependence on sparse/incomplete TRRUST/DoRothEA in brain; unresolved overfitting alternative weakens causal claim
2. Mediation bias: **Issue** — only 16 run-pairs (explicitly underpowered)
3. Detectability theory: **PASS**
4. Cross-context consistency: **Issue** — only 3 tissues with protocol heterogeneity; batch/assay confound
5. Perturbation validation: **Issue** — no result survives multiple-testing correction
6. Cross-species ortholog: **Issue** — single tissue/species pair; cell-composition confound
7. Pseudotime: **Issue** — small panel (56 pairs), underpowered
8. Batch/donor leakage: **Issue** — kidney n=1 donor
9. Uncertainty calibration: **PASS**
10. CSSI: **Issue** — circularity (same dataset for selection + evaluation)
11. Synthetic validation: **Issue** — generator may encode same assumptions (tautology risk)
12. Multi-model: **Issue** — one tissue/context, benchmark-context mismatch possible
