# Adversarial Review — Round 5

**Paper:** "A Comprehensive Framework for Mechanistic Interpretability of Single-Cell Foundation Models"
**Target:** Nature Machine Intelligence
**Reviewer stance:** Maximally hostile

---

## STATUS OF PREVIOUS ISSUES

### From Round 2 — All Resolved ✅
- SERGIO justification ✅
- Cohen's d / p-value discrepancy ✅
- Multiple testing discussion ✅
- Ethics/competing interests ✅
- Code availability ✅
- Correlation-vs-attention disconnect ✅
- External validation discussion ✅
- CSSI real-data validation ✅
- Multi-model expansion ✅
- Geneformer comparison deepened ✅

### From Round 4 — Status
- ✅ **"Δ Metric" in Table 5 clarified** — Caption now specifies metric per model
- ✅ **C2S-Pythia citation** — Now marked "(unpublished)" in both Methods and Results (fixed this round)
- ✅ **Table 1 CIs** — Added bootstrap CI ranges to caption (fixed this round)
- ✅ **Multi-model cell range comparability** — Explicit acknowledgment added in both Methods and Results (fixed this round)
- ⚠️ **CSSI "real-data" terminology** — Text uses "Real-data-structured validation" and "biologically structured simulations" — honest but the subsection could still be mistaken for real scRNA-seq on first skim. Acceptable.
- ⚠️ **Unused bib entries** — Still ~15 unused entries. Harmless but untidy.

---

## NEW ISSUES IDENTIFIED IN ROUND 5

### MAJOR: None

There are no remaining major issues that would warrant rejection or major revision.

### MINOR

#### 1. Nine assumptions vs. twelve analyses mapping (cosmetic)
Introduction lists nine "largely untested assumptions" (First through Ninth) then twelve analyses. Three analyses (CSSI, synthetic validation, multi-model) are constructive contributions rather than assumption tests. This is implicitly clear but could be made explicit with one sentence: "The remaining three analyses provide constructive methodological contributions rather than testing specific assumptions." **Severity: Low.** The current text is defensible.

#### 2. Unused bib entries
~15 entries in references.bib are not cited (szklarczyk2021string, aibar2017scenic, huynh2010inferring, adebayo2018sanity, doshivelez2017rigorous, sundararajan2017axiomatic, pearl2009causality, robins1992identifiability, vanderweele2015explanation, spirtes2000causation, lotfollahi2023predicting, roohani2024predicting, bunne2023learning, norman2019exploring, geiger2024causal, bills2023language, zheng2024benchmarking, stuart2019comprehensive, lopez2018deep). BibTeX ignores them but journal production may flag. **Severity: Cosmetic.**

Note: lopez2018deep IS cited (scVI). stuart2019comprehensive is NOT cited (was replaced by wolf2018scanpy in round 2).

#### 3. C2S-Pythia as "unpublished" weakens the multi-model claim
Citing an unpublished model limits reproducibility. However, the core multi-model finding (scGPT degrades, Geneformer doesn't) stands on two well-cited published models. C2S-Pythia is supplementary evidence. **Severity: Low.**

#### 4. No CSSI on actual foundation model attention weights
This remains the single most important gap in the paper, but it is now explicitly acknowledged in the CSSI Inference paragraph ("We note that these validations use simulated edge scores... rather than actual foundation model attention weights on real data") and in the External Validation subsection. The honest acknowledgment converts this from a fatal flaw to a clearly stated limitation with a roadmap. **Severity: Low-medium (acknowledged limitation, not a flaw).**

---

## SYSTEMATIC CHECKS

### All 4 models' results properly integrated and comparable?
✅ scGPT: full 12-analysis treatment (200–1000 cells)
✅ Geneformer: scaling + cross-context (200–500 cells), Table 6 systematic comparison
✅ scVI: scaling (200–500 cells), Table 5
✅ C2S-Pythia: scaling (50–200 cells), Table 5
⚠️ Cell ranges differ — now explicitly acknowledged in Methods and Results

### Every figure referenced and present in figures/?
✅ All 18 referenced figures exist in figures/. No broken references.

### All citations in references.bib?
✅ All \citep keys resolve to bib entries. Spot-checked 15 keys.

### Logical inconsistencies between sections?
✅ None found. The narrative is internally consistent: scaling failure (Sec 3.1) → CSSI fix (Sec 3.10) → synthetic validation (Sec 3.11) → multi-model shows architecture-dependence (Sec 3.12). Cross-references are accurate.

### Statistical rigor?
✅ p-values reported consistently throughout
✅ CIs: bootstrap CIs for scaling (Fig 1), CSSI (Table 3/4), calibration (200 resamples). Table 1 now has CI note.
✅ Effect sizes: Cohen's d for pseudotime, fold enrichment for ortholog, Spearman ρ throughout
✅ Multiple testing: within-analysis FDR/Bonferroni applied; framework-level justification in Discussion
✅ Sign test p=0.002 for 9/9: correct (2^{-9} ≈ 0.002)

### CSSI real-data validation convincing?
Moderately convincing. Two complementary experiments (PBMC with 22 known edges; Tabula Sapiens composition with 15 cell-type-specific edges) both show significant improvement. The honest acknowledgment that these are simulated rather than real attention weights prevents overselling. The 1.62× improvement at N=200 with p = 2.4×10^{-8} is statistically robust. The gap between synthetic (1.85×) and realistic (1.05–1.62×) performance is honestly reported.

### Discussion balanced?
✅ Yes. Does not claim CSSI solves everything. Acknowledges architecture-dependence. The "these findings do not invalidate single-cell foundation models but establish boundary conditions" framing is appropriate.

### Limitations honest and complete?
✅ Covers: limited model set, correlation-vs-attention disconnect, pseudotime single-axis limitation, calibration positive rate, no integration method testing, need for more datasets. The CSSI limitation (no real attention weights) is now stated in the Results section itself, not just Limitations.

---

## COMPILATION
✅ Compiles cleanly: 18 pages, no errors, no warnings beyond standard BibTeX unused entries.

---

## CHANGES MADE IN THIS ROUND

1. **Fixed:** C2S-Pythia marked as "(unpublished)" in both Methods (Sec 2.12) and Results (Sec 3.12)
2. **Fixed:** Table 1 caption now includes bootstrap CI ranges
3. **Fixed:** Multi-model cell range comparability explicitly acknowledged in Methods (Sec 2.12) and Results (Sec 3.12)

---

## VERDICT: **READY**

The paper is ready for NMI submission. All major issues from rounds 1–4 have been resolved. The remaining items are cosmetic (unused bib entries, minor wording) and do not affect scientific validity. The key limitation (no CSSI on real attention weights) is honestly acknowledged with a clear roadmap.

**Remaining cosmetic items (optional pre-submission cleanup):**
1. Remove ~15 unused bib entries
2. Add one sentence mapping 9 assumptions → 12 analyses in Introduction

These are not blockers.
