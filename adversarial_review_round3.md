# Adversarial Review — Round 3

**Paper:** "A Comprehensive Framework for Mechanistic Interpretability of Single-Cell Foundation Models"
**Target:** Nature Machine Intelligence
**Reviewer stance:** Hostile but constructive

---

## STATUS OF ROUND 2 ISSUES

### 1. ✅ CSSI real-data validation — ADDRESSED
CSSI now has two validation experiments beyond pure synthetic: (a) PBMC cell-type data with TRRUST edges (1.16× improvement), and (b) attention-model simulation using Tabula Sapiens immune cell-type proportions (1.62× at N=200, Wilcoxon p=2.4e-8). Table 5 provides clear quantitative results. **However, see concern #1 below.**

### 2. ✅ Geneformer section depth — ADDRESSED
Section 3.11 now includes: attention structure analysis (sparsity, entropy), cross-context similarity table (Table 4), scVI and C2S-Pythia comparisons (Table 3), and architectural interpretation. Proportionate to scGPT depth.

### 3. ✅ External validation dataset — ADDRESSED
New subsection "External Validation and Generalization" honestly discusses planned Replogle et al. 2022 validation. Appropriately caveated as future work.

### 4. ✅ Minor issues — MOSTLY ADDRESSED
- SERGIO justification: ✅ Added to Methods
- Cohen's d discussion: ✅ Power discussion added
- Ethics/competing interests: ✅ Added
- Code availability: ✅ GitHub link added
- Correlation-vs-attention disconnect: ✅ Acknowledged in Limitations
- Multiple testing: ✅ Note added to Discussion

---

## NEW CONCERNS (Round 3)

### 1. CSSI "real-data" validation is still simulated (MEDIUM-HIGH)
Both CSSI validation experiments use *simulated* data structured to match real biology, not actual scGPT/Geneformer attention weights on real cells. The PBMC experiment uses Spearman correlations as a "proxy for attention-derived scores," and the attention-model experiment uses synthetic attention contributions. A reviewer will note this gap: CSSI is validated on simulations that are *designed to show the effect*. The paper should be more explicit that no experiment applies CSSI to actual foundation model attention matrices on real data, and explain why (computational cost, model access). **Recommendation:** Add a sentence explicitly acknowledging this limitation and framing it as a validation priority.

### 2. Paper length remains 2× NMI limit (HIGH — but flagged as deferred)
At 18 pages two-column, this is ~10,000+ words. NMI Articles allow ~5,000 words. This was flagged in Round 2 and deferred. The paper should be either: (a) restructured as a Review/Perspective, (b) drastically cut with Methods and most Results in Supplementary, or (c) submitted to a venue with longer format (e.g., Genome Biology, Nature Methods). **This remains the most significant structural issue for NMI submission.**

### 3. Geneformer cross-context similarity may indicate a problem (MEDIUM)
The paper notes Geneformer's 0.979 cross-context similarity and raises the concern that this may reflect housekeeping genes rather than context-specific regulation. This is insightful but insufficiently explored. If Geneformer's attention patterns are nearly identical across brain, liver, and immune tissue, what biological information are they actually capturing? The paper should either: (a) show which genes dominate Geneformer's top attention positions and whether they are housekeeping, or (b) explicitly test whether Geneformer attention recovers context-specific edges at all.

### 4. scVI and C2S-Pythia claims lack detail (MEDIUM)
Table 3 includes scVI (−1.0% change) and C2S-Pythia (+3.0% change) but the Methods section (2.12) only mentions Geneformer. Where are the scVI/C2S-Pythia methods? What metrics were used? The cell ranges differ (50-200 for C2S-Pythia vs 200-1000 for scGPT) making direct comparison problematic. These entries feel like they were added for breadth without equivalent rigor.

### 5. CSSI convergence at large N undermines the practical message (LOW-MEDIUM)
Table 5 shows CSSI and pooled converge at N≥5000. This is theoretically expected but practically awkward: most single-cell studies now have >10,000 cells, so CSSI's advantage may be limited to rare-cell-type edges in small studies. The paper should discuss this explicitly: when is CSSI actually needed in practice? The answer is "whenever you have rare cell types of biological interest" — but this should be stated.

### 6. Abstract still ~350 words (LOW — deferred with length)
NMI limit is ≤150 words. Current abstract is 2.3× the limit.

### 7. Unused .bib entries remain (LOW)
~18 unused bibliography entries still in references.bib. Should be cleaned.

### 8. 20 figures/tables in main text (LOW — deferred with length)
NMI typically allows 6-8 main figures. Current count far exceeds this. Related to length issue.

---

## STRENGTHS (vs Round 2)

1. **CSSI validation substantially improved.** Two complementary experiments with quantitative results.
2. **Multi-model section now proportionate.** Tables 3-4, attention structure analysis, architectural interpretation.
3. **External validation honestly discussed.** Replogle as planned work is appropriate if caveated.
4. **Framework-level multiple testing acknowledged.** The note in Discussion is well-reasoned.
5. **Minor issues systematically fixed.** Ethics, code, SERGIO, Cohen's d all addressed.

---

## PRIORITY ACTIONS FOR NEXT ROUND

1. **HIGH:** Address paper length — either cut to NMI format or change target venue
2. **MEDIUM:** Add explicit acknowledgment that CSSI has not been tested on real foundation model attention weights
3. **MEDIUM:** Add Methods for scVI/C2S-Pythia or remove from Table 3
4. **MEDIUM:** Discuss practical guidance for when CSSI is needed (rare cell types, small studies)
5. **LOW:** Explore Geneformer's top attention positions (housekeeping vs. context-specific)
6. **LOW:** Clean .bib, reduce figures for submission

---

## OVERALL ASSESSMENT

The paper has improved significantly from Round 2. The three critical gaps (CSSI validation, Geneformer depth, external validation) have all been addressed, though with varying levels of completeness. The remaining issues are primarily structural (length) rather than scientific. The core scientific contributions — scaling failure documentation, CSSI as a constructive solution, multi-model comparison — are now well-supported. **If the length issue is resolved, this paper is approaching submission-readiness for NMI.**
