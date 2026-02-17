# Adversarial Review — Round 2

**Paper:** "A Comprehensive Framework for Mechanistic Interpretability of Single-Cell Foundation Models"
**Target:** Nature Machine Intelligence
**Reviewer stance:** Hostile but constructive

---

## MAJOR CONCERNS (potential grounds for rejection)

### 1. Paper length vastly exceeds NMI limits
NMI Articles allow ~3,000–5,000 words of main text (~6–8 typeset pages). This manuscript is **16 pages in two-column format** with ~8,000+ words. This is approximately 2× the allowed length. **The paper must either be restructured as a longer-format Review/Perspective or drastically cut.** Recommendation: move Methods entirely to Supplementary, condense Results subsections, and trim the Discussion.

### 2. CSSI validated only on synthetic data — no real-data validation
The paper's main constructive contribution (CSSI) is validated exclusively on synthetic experiments with oracle cell-state labels. There is **zero demonstration on real scGPT or Geneformer attention weights**. A reviewer would ask: "Does CSSI actually improve GRN recovery from real foundation model attention?" Without this, the central constructive claim is unsubstantiated on biological data. This is the single most important gap.

### 3. No external/independent validation
All analyses use Tabula Sapiens or the same few Perturb-seq datasets (Dixit, Adamson, Shifrut). No independent cohort or dataset validates the findings. For a framework paper claiming generality, at least one fully independent dataset (e.g., Replogle 2022 genome-scale Perturb-seq, or a non-Tabula Sapiens atlas) should be used.

### 4. Geneformer comparison is superficial
The multi-model section (Section 3.11) tests only 200 vs. 500 cells and reports "stable performance" and "average similarity of 0.979" without showing the same detailed metrics (F1, AUROC, TRRUST/DoRothEA benchmarks, seed robustness) used for scGPT. The comparison is asymmetric and unconvincing. A hostile reviewer would call this a strawman: scGPT gets exhaustive stress-testing while Geneformer gets a paragraph. Either provide equivalent depth or frame more carefully.

### 5. Extremely low absolute F1 scores undermine practical significance
The paper acknowledges F1 scores are very low (order 10⁻⁵ to 10⁻⁴) and argues the *direction* of change matters. While the defense ("needle in a haystack") has merit, a reviewer could argue: if absolute performance is at noise floor, measuring *trends* in noise is meaningless. The paper should show that the 200-cell F1 is statistically significantly above a random baseline (not just enrichment ratio) with proper confidence intervals.

### 6. Single-author paper making very broad claims
Twelve analyses spanning scaling, mediation theory, detectability, cross-species transfer, pseudotime, batch effects, calibration, synthetic validation, and multi-model comparison — all by a single author. Reviewers will question depth vs. breadth. Some analyses feel like separate papers stapled together rather than a cohesive framework. Consider whether the narrative could be tightened.

---

## MINOR CONCERNS (revision required)

### 7. Inconsistency fixed: "nine" vs. "twelve" analyses
~~The Introduction stated "nine complementary analyses" but enumerated twelve.~~ **Fixed in this round** — changed to "twelve."

### 8. Citation error: Scanpy pipeline cited as Stuart et al. (Seurat)
~~The Methods cited `stuart2019comprehensive` (Seurat/integration) for "Scanpy-based quality control."~~ **Fixed in this round** — now cites `wolf2018scanpy`.

### 9. BibTeX warning: niculescumizil2005predicting
~~Entry typed as `@article` but has `booktitle` field.~~ **Fixed in this round** — changed to `@inproceedings`.

### 10. Missing discussion of SERGIO
The paper cites SERGIO (Dibaeinia 2020) for ground-truth simulation but uses its own custom synthetic generator instead. Why not use the established SERGIO simulator? This choice should be justified — SERGIO is the standard benchmark.

### 11. Correlation-based edges ≠ attention-based edges in ortholog/batch sections
Sections 3.6 (ortholog), 3.7 (pseudotime), and 3.8 (batch) use Spearman/Pearson correlation-based edge scores, NOT attention-derived edges. Yet the paper frames itself as evaluating "mechanistic interpretability of foundation models." These sections evaluate correlation-based GRN inference generally, not foundation model interpretability specifically. This disconnect should be acknowledged more clearly.

### 12. Calibration ground truth has 45.8% positive rate
The paper acknowledges this in Limitations but it's a serious concern. Real GRN sparsity is ~1-5%. A 45.8% positive rate means the calibration analysis operates in a fundamentally different regime. Calibration curves and ECE values may not transfer to realistic sparsity levels.

### 13. No figure numbering scheme or supplementary organization
18 figures in the main text is excessive for NMI. Many should move to Extended Data or Supplementary. NMI typically allows 6–8 main figures.

### 14. Data/Code availability statements are vague
"Will be made available upon publication" is insufficient for NMI. Reviewers expect a GitHub link (even if private during review) or at minimum a detailed description of what will be released.

### 15. No ethics/competing interests statement
NMI requires explicit competing interests and ethics declarations.

---

## FACTUAL ERRORS OR INCONSISTENCIES

### 16. Enumeration count now consistent (fixed)
Introduction now correctly says "twelve complementary analyses."

### 17. Abstract mentions "twelve critical challenges" — the list in Introduction has 12 items. Consistent. ✓

### 18. Table 1 (robustness) shows 3 tiers × 2 metrics — no confidence intervals
The seed robustness table reports point estimates without uncertainty. Bootstrap CIs should be added for consistency with the rest of the paper.

---

## STATISTICAL CONCERNS

### 19. Multiple testing across the entire framework
The paper performs dozens of statistical tests across 12 analyses. There is no framework-level multiple testing correction. While each section internally applies FDR, the overall false discovery rate across the paper is uncontrolled.

### 20. Synthetic validation confirms predictions — but could be circular
The synthetic data is generated to match the theoretical model (steady-state GRN dynamics, attention = tanh(true + noise)). Of course the theory works when data matches the theory. A stronger test would use SERGIO or a generative model NOT designed to match the theoretical framework.

### 21. Cohen's d = 1.58 but p = 0.068 (pseudotime)
A large effect size with non-significant p-value suggests low power. This should be discussed — the sample of 56 pairs may be underpowered.

---

## WRITING QUALITY

### 22. Generally strong
The writing is clear, well-structured, and appropriately technical. Section headers with Evidence/Inference/Decision Implication structure is effective.

### 23. Abstract is long (~350 words)
NMI abstracts should be ≤150 words. Current abstract is ~2.3× the limit.

### 24. Some redundancy between Results and Discussion
The Discussion largely restates Results findings. Could be condensed.

---

## FIGURE/TABLE ISSUES

### 25. All 18 figures present in figures/ directory ✓
### 26. All figures referenced in text ✓
### 27. Several figures from earlier work (fig4–fig7, fig_ext_*, fig_matched_*, fig_snapshot_*, fig_simulation_*, fig_theory_*) are in figures/ but NOT referenced in main.tex
These orphan figures should either be referenced (as Extended Data) or removed from the directory.

---

## REFERENCE ISSUES

### 28. All \citep keys present in references.bib ✓
### 29. Several bib entries unused in text
`szklarczyk2021string`, `aibar2017scenic`, `huynh2010inferring`, `adebayo2018sanity`, `doshivelez2017rigorous`, `sundararajan2017axiomatic`, `pearl2009causality`, `robins1992identifiability`, `vanderweele2015explanation`, `spirtes2000causation`, `lopez2018deep`, `lotfollahi2023predicting`, `roohani2024predicting`, `bunne2023learning`, `norman2019exploring`, `geiger2024causal`, `bills2023language`, `zheng2024benchmarking` — these are in .bib but not cited. BibTeX will ignore them, but the .bib should be cleaned.

---

## CHANGES MADE IN THIS ROUND

1. **Fixed** "nine complementary analyses" → "twelve complementary analyses" in Introduction
2. **Fixed** Scanpy citation from `stuart2019comprehensive` → `wolf2018scanpy`
3. **Fixed** BibTeX entry type for `niculescumizil2005predicting` (`@article` → `@inproceedings`)

---

## PRIORITY ACTIONS FOR NEXT ROUND

1. **CRITICAL:** Validate CSSI on real scGPT/Geneformer attention data, not just synthetic
2. **CRITICAL:** Cut paper to NMI length (~5,000 words main text, ≤8 main figures, ≤150-word abstract) or reformat as Review/Perspective
3. **HIGH:** Deepen Geneformer comparison with equivalent metrics to scGPT
4. **HIGH:** Add independent external validation dataset
5. **HIGH:** Justify custom synthetic generator vs. SERGIO
6. **MEDIUM:** Acknowledge correlation-vs-attention disconnect in ortholog/batch sections
7. **MEDIUM:** Add competing interests and ethics statements
8. **MEDIUM:** Provide concrete code/data availability (GitHub link)
9. **LOW:** Clean unused .bib entries
10. **LOW:** Add CIs to Table 1
