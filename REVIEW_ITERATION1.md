# Adversarial Review — Iteration 1
**Paper:** A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models  
**Reviewer model:** Claude Opus 4  
**Date:** 2026-02-15  
**Overall Rating: 5.5 / 10**

---

## Executive Summary

The paper attempts something ambitious: a 12-analysis framework evaluating attention-based GRN inference from single-cell foundation models. The breadth is impressive. However, the paper suffers from a fundamental tension: **it makes strong negative claims about attention-based GRN inference while simultaneously claiming CSSI "fixes" the problem—but the fix is validated primarily on custom synthetic data designed to confirm the theory, and on real data where CSSI adds essentially zero improvement at the best-performing layers.** The paper is honest about many limitations (commendably so), but this honesty sometimes reads as preemptive defense rather than resolution.

---

## Critical Issues (Ranked by Impact)

### 1. **CSSI's Real-Data Performance Contradicts the Narrative** ⚠️ CRITICAL
**Impact: Very High**

The central constructive claim is that CSSI mitigates scaling failure. But Table 5 (real scGPT attention) reveals:
- At the best layers (L13-L14), CSSI provides **Δ = 0.000 to −0.001**—literally no improvement or slight degradation.
- The largest real-data improvement is Δ = +0.033 at Layer 10, which is modest.
- The 1.85× improvement figure comes entirely from synthetic data with oracle cell-state labels.

**The paper acknowledges this** ("CSSI's primary contribution is diagnostic") but continues to frame CSSI as a "constructive resolution to scaling failure." This framing mismatch is the paper's biggest credibility problem. A diagnostic tool that identifies good layers is useful but is not the same as "mitigating scaling failure."

**Fix:** Reframe CSSI honestly as a diagnostic/layer-selection framework, not a performance-enhancing method. The 1.85× number should not appear in the abstract without the caveat "synthetic only, oracle labels."

### 2. **Synthetic Validation Is Circular by Design** ⚠️ CRITICAL  
**Impact: Very High**

The synthetic generator creates attention matrices as `A_attention = tanh(A_true + noise)`. This *encodes the paper's theory* into the ground truth. The paper admits this ("confirms internal consistency of our theoretical framework rather than providing independent external validation") but still uses synthetic results as a primary evidence pillar.

The SERGIO justification is unconvincing: the claim that SERGIO "cannot model attention matrix formation" is true but misses the point. The issue is that the custom generator *assumes the relationship the paper is trying to prove*. Any method that pools attention will fail on data generated to have state-specific attention—this is a tautology, not validation.

**Fix:** Either (a) use a real foundation model's attention on data with known ground-truth perturbation effects (Replogle 2022 exists!), or (b) clearly label all synthetic results as "internal consistency checks" rather than "validation."

### 3. **The 0.694 AUROC Claim Needs Much Stronger Controls** ⚠️ HIGH
**Impact: High**

The headline real-data result—Layer 13 achieving AUROC 0.694 on 497 brain cells against TRRUST—is the paper's strongest positive finding. But:

- **N = 497 cells, single tissue, single dataset.** This is one experiment.
- **Layer selection on the same data used for reporting.** The paper acknowledges circularity but the cross-validation (248/249 split) produced *different* top layers (L17, L16, L13 vs. L13, L14). This is not "robust"—it's unstable.
- **No comparison to simple baselines at the layer level.** Does computing Spearman correlations within cell types at these same 497 cells achieve similar AUROC? The baseline comparison (Section 3.2) was done with pooled 500 cells on different tissue. A fair comparison would run GENIE3/GRNBoost2 on the same 497 brain cells, stratified by cell type, against the same TRRUST edges.
- **The scGPT cross-validation section contradicts itself:** it says top layers from Half A were [0, 4, 10] while the primary analysis found L13-L14. Layers 0 and 4 being "top performers" in one split while L13-L14 dominate in the full data suggests extreme instability.

**Fix:** Run the layer-stratified analysis on at least 2 additional tissues/datasets. Add stratified baselines (correlation within cell types) to demonstrate that attention adds value beyond what simple methods achieve with the same stratification.

### 4. **Claims–Evidence Mismatch Throughout** ⚠️ HIGH
**Impact: High**

Several results are presented with language stronger than the evidence supports:

| Claim | Evidence | Problem |
|-------|----------|---------|
| "Scaling failure" | Non-monotonic curve peaking at 750 cells, all CIs overlap | Not "failure"—it's a plateau with noisy decline at one point |
| "Mediation is non-additive" | 10/16 run-pairs, 3 tissues | N=16 explicitly acknowledged as underpowered |
| "Perturbation validation" | Nothing survives FDR correction | Honest, but then why present it as a framework component? |
| "Pseudotime fails for ~79%" | p=0.124 after FDR | A null result, not a "failure" finding |
| "CSSI mitigates scaling failure" | Only on synthetic data with oracle labels | See Issue #1 |

**Fix:** Systematically downgrade language. "Scaling plateau with possible decline" not "failure." "Preliminary evidence of non-additivity" not "systematic bias." Frame perturbation and pseudotime sections as establishing that current validation approaches are underpowered, not that they "fail."

### 5. **Reference Database Saturation Analysis Is Not Convincing** ⚠️ MEDIUM
**Impact: Medium**

The saturation analysis (claiming scaling failure is *worse* with complete references) uses synthetic data with `σ = 0.05 + 0.002 × cell_count`—noise that *increases linearly with cell count by construction*. Of course performance degrades with more cells when you add more noise! This is not evidence that scaling failure is "genuine"—it's evidence that the simulation assumes scaling failure.

**Fix:** Remove or reframe. The noise model must be justified from real data characteristics, not chosen to produce the desired result.

### 6. **Scope vs. Depth Trade-off Hurts the Paper** ⚠️ MEDIUM
**Impact: Medium**

12 analyses × limited depth = many weak findings rather than a few strong ones. The cross-species, pseudotime, batch leakage, and calibration sections are each individually publishable topics that get insufficient treatment here. The result is a paper where:
- 3-4 analyses produce null results after correction
- 2-3 analyses are acknowledged as underpowered
- The strongest positive result (L13-L14 AUROC) is from a single dataset

**Fix:** Consider focusing the paper on 4-5 analyses with deeper treatment: scaling behavior, CSSI with real-data validation on multiple tissues, multi-model comparison, and baseline comparison. Move the rest to supplementary.

### 7. **Multiple Testing Correction Is Applied Inconsistently** ⚠️ MEDIUM
**Impact: Medium**

The paper applies BH-FDR across all 47 tests framework-wide, which is admirably conservative. But:
- The scaling analysis reports "adjusted p = 0.011" but the scaling curve has overlapping CIs at every adjacent point. What exactly is being tested?
- The CSSI Wilcoxon test (p = 2.5 × 10⁻¹¹) is on synthetic data and shouldn't be in the same correction family as real-data tests.
- Mixing synthetic-data tests with real-data tests in the same FDR correction family is methodologically questionable.

**Fix:** Separate the correction into real-data and synthetic-data families. Be explicit about what null hypothesis each p-value tests.

### 8. **The "Attention Encodes Co-expression Not Regulation" Finding Is Underexploited** ⚠️ MEDIUM
**Impact: Medium**

The finding that attention correlates with co-expression (ρ = 0.31–0.42) but not regulation (ρ = −0.01–0.02) is actually the most important mechanistic insight in the paper. But it gets 3 sentences in Section 3.10. This should be a central result with deeper analysis: Which layers show the strongest co-expression correlation? Does CSSI change this pattern? Does the co-expression signal in L13-L14 look different from early layers?

**Fix:** Expand this analysis significantly. It could unify the scaling, CSSI, and multi-model findings.

---

## Moderate Issues

### 9. **Geneformer Evaluation on Wrong Tissue**
Geneformer is evaluated on DLPFC brain tissue, which the paper itself acknowledges has poor reference database coverage. The paper then compares TRRUST AUROC between scGPT (evaluated on Tabula Sapiens) and Geneformer (DLPFC)—different tissues, different datasets. Table 4 presenting them side-by-side implies direct comparability that doesn't exist.

### 10. **The Abstract Is Misleading**
The abstract states "CSSI cell-state stratification provides... modest improvements in intermediate layers" but leads with "improving GRN recovery by up to 1.85×"—a synthetic-only result. The casual reader will remember 1.85× and miss "synthetic only, oracle labels."

### 11. **Missing Critical Baseline: Stratified Correlation**
The paper never tests whether simple Spearman correlation *within cell types* (i.e., the same stratification CSSI uses but without attention) achieves comparable AUROC to CSSI on attention matrices. If stratified correlation matches CSSI-on-attention, then attention adds nothing—the improvement comes entirely from stratification, applicable to any method.

### 12. **Power Analysis for Perturbation Section Is Post-Hoc**
The paper computes post-hoc power analysis showing n ≈ 300 targets needed. Post-hoc power analysis is widely criticized (Hoenig & Heisey 2001) and doesn't add information beyond the observed p-values and effect sizes.

---

## Minor Issues

13. **Two-column format with 0.8in margins** makes the paper very dense and hard to read.
14. **"9 assumptions" in intro** is a list that could be a paragraph. The enumeration suggests systematic testing of each, but assumptions 6-9 get uneven coverage.
15. **The Discussion** re-summarizes results rather than synthesizing. It reads as "Section X showed Y" repeated 12 times.
16. **Code/data availability** says "DOI: to be assigned" and "[repository URL]"—placeholder text.

---

## What the Paper Does Well

- **Honest about limitations:** The paper flags its own weaknesses (underpowered mediation, circular synthetic validation, post-hoc power). This is rare and commendable.
- **Breadth of evaluation:** Testing across scGPT, Geneformer, scVI, C2S-Pythia is valuable.
- **Practical recommendations:** The 10-point recommendation list is genuinely useful for practitioners.
- **Layer-stratification insight:** The finding that L13-L14 encode regulatory signal while early layers don't is interesting and novel, even if undervalidated.
- **Baseline comparison:** Showing GENIE3/GRNBoost2 also fail on brain tissue is an important context.

---

## Actionable Fixes Ranked by Impact

| Rank | Fix | Effort | Impact |
|------|-----|--------|--------|
| 1 | Validate L13-L14 finding on ≥2 additional tissues with stratified correlation baselines | High | Critical |
| 2 | Reframe CSSI as diagnostic/layer-selection tool, not performance enhancer | Low | High |
| 3 | Remove or relabel synthetic results as "internal consistency" not "validation" | Low | High |
| 4 | Expand attention-encodes-coexpression analysis as central mechanistic finding | Medium | High |
| 5 | Downgrade language throughout to match evidence strength | Low | High |
| 6 | Add stratified-correlation baseline (same cell-type splits, no attention) | Medium | High |
| 7 | Run real-data CSSI on Replogle 2022 perturbation atlas | High | High |
| 8 | Fix noise model in saturation analysis or remove section | Low | Medium |
| 9 | Separate FDR correction into real vs. synthetic families | Low | Medium |
| 10 | Consolidate to 5-6 analyses with deeper treatment, move rest to supplement | Medium | Medium |
| 11 | Fix Geneformer/scGPT comparison to use matched tissue | High | Medium |
| 12 | Remove post-hoc power analysis from perturbation section | Low | Low |

---

## Rating Justification: 5.5/10

**Above average for ambition and honesty; below average for evidence strength.**

The paper identifies real problems (attention doesn't encode regulation, scaling is non-trivial, reference databases are limited) but overclaims its solutions. CSSI is framed as a fix but provides no improvement on the best real-data results. The strongest finding (layer-specific regulatory signal) rests on a single dataset with unstable cross-validation. The synthetic validation is circular. Multiple analyses produce null results after correction.

For a top venue (JMLR), I would recommend **major revision**: focus on fewer claims with stronger evidence, validate on additional datasets, and honestly reframe CSSI as a diagnostic tool rather than a scaling-failure solution.
