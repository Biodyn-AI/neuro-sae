# Adversarial Review — Round 4

**Paper:** "A Comprehensive Framework for Mechanistic Interpretability of Single-Cell Foundation Models"
**Target:** Nature Machine Intelligence
**Reviewer stance:** Hostile but constructive

---

## STATUS OF PREVIOUS ISSUES

### Round 2 Issues — Resolved
- ✅ **SERGIO justification** — Now explained in Methods (Section 2.10): custom generator needed for attention matrix modeling
- ✅ **Cohen's d / p-value discrepancy** — Now explained as low power from 56-pair sample
- ✅ **Multiple testing across framework** — Addressed in new Discussion paragraph
- ✅ **Ethics/competing interests** — Present
- ✅ **Code availability** — GitHub link provided with reviewer access offer
- ✅ **Correlation-vs-attention disconnect** — Acknowledged in Limitations
- ✅ **External validation discussion** — New subsection in Discussion (Section 4.4)
- ✅ **CSSI real-data validation** — PBMC-like simulation with 22 known edges added
- ✅ **Multi-model expansion** — scVI and C2S-Pythia added (Table 5)
- ✅ **Geneformer comparison deepened** — Table 4 with systematic property comparison

### Round 2 Issues — Partially Addressed
- ⚠️ **Paper length** — Still ~18 pages in two-column. Exceeds NMI Article limits but acceptable for Review/Perspective format. *Instruction says to ignore length concerns.*
- ⚠️ **Abstract length** — Still ~350 words vs NMI's 150-word limit. Same caveat.
- ⚠️ **No CSSI on real attention weights** — PBMC validation uses Spearman correlation as proxy, not actual scGPT/Geneformer attention. Acknowledged in External Validation subsection as future work.
- ⚠️ **Unused bib entries** — Still present but harmless (BibTeX ignores them)

---

## NEW ISSUES IDENTIFIED IN ROUND 4

### MAJOR

#### 1. Multi-model comparison depth remains asymmetric
scVI and C2S-Pythia get 1 paragraph and a table row each. The "metric change" values (-1.0% for scVI, +3.0% for C2S-Pythia) lack specificity: metric of *what*? Cell range for C2S-Pythia is only 50-200 cells (not comparable to the 200-1000 range used for scGPT). The Geneformer comparison (200-500) is also a narrower range than scGPT's 200-1000. While the directional finding is clear (scGPT degrades, others don't), the comparison is not fully apples-to-apples.

**Severity:** Medium. The qualitative conclusion holds but quantitative comparison is weakened.

#### 2. CSSI real-data validation is still simulation-based
The "real-data validation" in CSSI (Section 3.10) uses simulated data with PBMC-like structure, not actual PBMC scRNA-seq data. The term "real-data" is misleading. The improvement (1.05-1.16x) is modest compared to synthetic (1.85x). This should be framed more carefully — it's a "realistic simulation" not "real data."

**Severity:** Medium. Already improved from prior round but terminology could mislead.

#### 3. Table 1 (robustness) still lacks confidence intervals
Round 2 flagged this. Point estimates without uncertainty remain.

**Severity:** Low-medium.

### MINOR

#### 4. Introduction lists 9 assumptions but 12 analyses
The Introduction lists nine "largely untested assumptions" (First through Ninth) but twelve analyses. The mapping between assumptions and analyses is not explicit. Three analyses (CSSI, synthetic validation, multi-model) are constructive contributions rather than assumption tests, which is fine but should be stated.

#### 5. C2S-Pythia citation missing
C2S-Pythia is mentioned but not cited. No bib entry visible for it.

#### 6. scVI citation
The paper cites `lopez2018deep` for scVI, which is correct. ✓

#### 7. Multi-model Methods section is brief
Section 2.12 (Multi-Model Validation) mentions scVI and C2S-Pythia were tested but doesn't describe the specific metrics or evaluation protocol used for them.

#### 8. Extended model comparison table (Table 5) uses "Δ Metric" without defining the metric
What metric changed by -1.0% for scVI? Reconstruction loss? Latent similarity? This should be specified.

---

## FACTUAL CHECKS

- ✅ Abstract: "twelve critical challenges" matches 12 enumerated analyses
- ✅ All 18 figure references resolve to files in figures/
- ✅ All \citep keys in text appear in references.bib (spot-checked: lopez2018deep, replogle2022mapping, dibaeinia2020sergio)
- ✅ Bootstrap CI [0.008, 0.091] and p=0.030 from CSSI real-data match script output
- ✅ CSSI-max TP@25=22 vs pooled TP@25=19 matches script output
- ✅ Equation numbering: 3 equations, all referenced
- ✅ Table numbering: 5 tables referenced (robustness, topk_ortholog, cssi_scaling, cssi_realdata, model_comparison) + extended_model_comparison = 6 tables

---

## STATISTICAL CHECKS

- ✅ Sign test p=0.002 for 9/9 unanimous: correct (2^{-9} ≈ 0.002)
- ✅ Wilcoxon p=0.002 consistent with 9 paired observations
- ✅ Bootstrap 100 resamples: p=0.030 with 97% win rate → 3/100 resamples with diff ≤ 0 → p=0.03. ✓
- ✅ Permutation null ρ=0.011±0.008, z=92.6 for ortholog transfer: plausible given n=25,876

---

## WRITING QUALITY

- Strong overall. Evidence/Inference/Decision Implication structure is effective and consistent.
- The new multi-model paragraph integrates smoothly.
- The multiple testing discussion in Discussion is well-placed and well-reasoned.
- Minor: "scVI~\citep{lopez2018deep}" appears in both Methods and Results — consistent. ✓

---

## OVERALL ASSESSMENT

The paper has improved substantially from rounds 1-2. Key additions:
1. CSSI real-data validation (modest but significant improvement)
2. Extended multi-model comparison (scVI, C2S-Pythia)  
3. Multiple testing justification
4. Cohen's d power discussion
5. External validation roadmap

**Remaining weaknesses** are primarily: (a) the gap between CSSI's synthetic success (1.85x) and realistic-simulation success (1.05-1.16x) raises questions about practical utility; (b) no analysis on actual foundation model attention weights (everything uses correlation-based edge scores or synthetic attention); (c) multi-model comparison remains shallower for non-scGPT models.

**These are addressable in revision** and do not undermine the paper's core contributions: the systematic identification of failure modes, the theoretical framework, and the CSSI methodology.

**Recommendation:** Accept with minor revision (assuming format compliance with NMI).

---

## PRIORITY ACTIONS FOR NEXT ROUND

1. **MEDIUM:** Clarify "Δ Metric" in Table 5 — specify what metric is being measured for each model
2. **MEDIUM:** Add C2S-Pythia citation or note it's unpublished
3. **LOW:** Add bootstrap CIs to Table 1 (robustness)
4. **LOW:** Rename CSSI "real-data validation" to "realistic-simulation validation" for accuracy
5. **LOW:** Clean unused bib entries
