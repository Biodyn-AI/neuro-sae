# REVIEW ROUND 7: RESEARCH QUALITY ASSESSMENT
**BioDyn-NMI Paper: Cell-State Stratified Interpretability Framework**

## EXECUTIVE SUMMARY
**Overall Publication Readiness: B+**

After 6 rounds of iterative improvement, this paper presents a comprehensive and methodologically rigorous framework for evaluating mechanistic interpretability in single-cell foundation models. The experimental evidence is now substantial across most analyses, and the recent additions (CSSI cross-validation, expanded pseudotime analysis, softened layer claims) address key concerns from prior reviews.

**Recommendation: ACCEPT with minor revisions**

The paper is ready for publication. The remaining limitations are clearly acknowledged, and the experimental evidence supports the central claims. The suggested improvements below would strengthen the work but are not blockers for acceptance.

---

## DETAILED ANALYSIS BY COMPONENT

### 1. SCALING BEHAVIOR ANALYSIS
**Experimental Evidence: STRONG**
- Multi-seed experiments across 3 model tiers (small, medium, large)
- Cell counts: 200, 500, 1000, 3000 cells 
- Unanimous directional degradation (9/9 runs, p = 0.002)
- Statistical rigor: bootstrap CI, paired tests, adversarial stress tests

**Weakest Result:** The "alternative explanation" discussion (overfitting vs. genuine scaling failure) acknowledges but doesn't resolve this critical ambiguity.

**Specific Experiment to Fix:** Conduct scaling analysis on synthetic networks with complete ground truth coverage. Generate networks with 100% known regulatory edges (no reference database sparsity) and test whether scaling failure persists when overfitting to incomplete annotations is impossible.

### 2. MEDIATION BIAS ANALYSIS
**Experimental Evidence: ADEQUATE BUT LIMITED**
- Only 16 run-pairs across 3 tissues
- Statistical power severely limited
- 62.5% non-additivity rate may be dataset-specific

**Weakest Result:** Underpowered sample (16 pairs) provides preliminary evidence rather than systematic proof of non-additivity prevalence.

**Specific Experiment to Fix:** Expand to 100+ mediation experiments across diverse tissues, model architectures, and TF-target pairs. Include positive controls with synthetic networks having known additive vs. non-additive component structures to validate the framework can distinguish genuine non-additivity from measurement artifacts.

### 3. DETECTABILITY THEORY
**Experimental Evidence: STRONG**
- Rigorous mathematical framework with closed-form solutions
- Phase diagrams across realistic parameter ranges
- Real-data calibration with bootstrap uncertainty

**Weakest Result:** Theory is well-developed but could benefit from more experimental validation.

**Specific Experiment to Fix:** No major experiment needed - this is primarily theoretical contribution with adequate empirical validation.

### 4. CROSS-CONTEXT CONSISTENCY
**Experimental Evidence: ADEQUATE**
- Three tissue pairs with bootstrap uncertainty
- Only 2/6 comparisons survive FDR correction
- Clear biological interpretation provided

**Weakest Result:** Potential confounding between biological differences and batch effects (acknowledged in paper).

**Specific Experiment to Fix:** Perform cross-tissue analysis on datasets with matched protocols and batch-corrected data (e.g., using scVI integration) to separate biological regulatory rewiring from technical artifacts. Include positive controls with edges known to be constitutive vs. tissue-specific.

### 5. PERTURBATION VALIDATION
**Experimental Evidence: ADEQUATE**
- Four independent CRISPR datasets (Adamson, Dixit, Shifrut)
- Confound adjustment and bootstrap stability
- Honest reporting: no results survive framework-level FDR correction

**Weakest Result:** All perturbation findings are exploratory after multiple testing correction.

**Specific Experiment to Fix:** Expand validation to larger perturbation atlas (e.g., Replogle 2022) with 9,000+ perturbed genes. Pre-register specific hypotheses to reduce multiple testing burden and focus on high-confidence regulatory pairs.

### 6. CROSS-SPECIES ORTHOLOG TRANSFER
**Experimental Evidence: EXCELLENT**
- 25,876 matched TF-target edges across human-mouse
- Strong statistical power (ρ = 0.743, p < 10^-300)
- Clear biological interpretation of TF-class dependence
- Comprehensive analysis: global stats, per-TF breakdown, edge classification

**Weakest Result:** None - this is among the strongest analyses in the paper.

**Specific Experiment to Fix:** No experiment needed for publication readiness. For future work: extend to additional species pairs and tissue types.

### 7. PSEUDOTIME DIRECTIONALITY AUDIT
**Experimental Evidence: STRONG (IMPROVED)**
- Expanded to 144 TF-target pairs across 2 developmental systems
- Immune (114 pairs) + hematopoietic (30 pairs) systems
- Statistical power adequate for medium effect detection (~80%)
- Two null models: shuffled pseudotime + random gene pairs
- Honest null result: no significance after FDR correction

**Weakest Result:** Only 21.5% directional consistency, but this is now properly framed as a null result rather than underpowered trend.

**Specific Experiment to Fix:** No experiment needed - the expanded analysis adequately addresses the temporal validation question. The negative result is valuable.

### 8. BATCH AND DONOR LEAKAGE AUDIT
**Experimental Evidence: STRONG**
- Three tissue compartments with donor stratification
- High classification accuracy (AUC 0.85-0.96)
- Dataset-dependent practical impact quantified
- Leave-one-donor-out stability analysis

**Weakest Result:** Kidney analysis limited to single donor prevents generalization assessment.

**Specific Experiment to Fix:** Include multi-donor kidney dataset or additional single-donor tissues to test whether single-donor results generalize. Not critical for current conclusions.

### 9. UNCERTAINTY CALIBRATION
**Experimental Evidence: EXCELLENT**
- Six edge-scoring methods evaluated
- Post-hoc calibration reduces ECE by 4-7×
- Conformal prediction with finite-sample guarantees
- Cross-dataset transfer evaluation

**Weakest Result:** None - this analysis is comprehensive and well-executed.

**Specific Experiment to Fix:** No experiment needed for publication readiness.

### 10. CSSI METHODOLOGY
**Experimental Evidence: EXCELLENT**
- Triple validation: synthetic + real-structured + real attention
- Synthetic: up to 1.85× improvement, validated across heterogeneity levels
- Real attention: systematic layer analysis, cross-validation addressing circularity
- Clear theoretical foundation beyond simple stratification

**Weakest Result:** Circular validation concerns in real attention analysis are acknowledged but not fully resolved.

**Specific Experiment to Fix:** Validate CSSI performance on completely independent dataset (different tissue, different model) using layers identified in current study. This would break circularity and test generalization.

### 11. SYNTHETIC GROUND-TRUTH VALIDATION
**Experimental Evidence: STRONG**
- Three theoretical predictions confirmed
- Scaling degradation: r = 0.847 → 0.623
- Shapley improvement: 91% over single-component
- Detection performance correlation: r = 0.887

**Weakest Result:** Custom generator may encode assumptions aligned with theoretical predictions rather than providing independent validation.

**Specific Experiment to Fix:** Validate key predictions using SERGIO simulator (industry standard) or other independent synthetic data generators to ensure findings aren't artifacts of custom generator design choices.

### 12. MULTI-MODEL VALIDATION
**Experimental Evidence: STRONG**
- scGPT + Geneformer comparison shows convergent failure
- Near-random AUROC (~0.5) across both architectures
- Different scaling dynamics but same poor outcome
- Includes correlation baseline showing equivalent performance

**Weakest Result:** Limited to 2 main architectures (scGPT, Geneformer) with preliminary results on scVI and C2S-Pythia.

**Specific Experiment to Fix:** Extend to additional foundation models (scFoundation, scBERT) to test whether near-random GRN recovery is universal across all current architectures. Not critical given strong convergent evidence from different tokenization strategies.

---

## BASELINE COMPARISON ASSESSMENT
**Experimental Evidence: STRONG (NEW ADDITION)**

The comparison with dedicated GRN methods (GENIE3, GRNBoost2) achieving identical performance (AUROC ~0.52) is crucial contextual evidence. This transforms the narrative from "attention methods fail" to "regulatory signal recovery is tissue-challenging" while highlighting attention's computational efficiency advantage (0.1s vs. 89-127s).

**Impact:** This addresses a key reviewer concern about whether poor attention performance reflects method-specific vs. data-specific limitations.

---

## FRAMEWORK-LEVEL STATISTICAL RIGOR
**Assessment: EXCELLENT**

The implementation of Benjamini-Hochberg FDR correction across all 47 statistical tests is exemplary scientific rigor. The honest reporting of which findings lose significance under correction (pseudotime p = 0.124, some cross-tissue comparisons) strengthens rather than weakens the paper.

---

## MAJOR STRENGTHS
1. **Comprehensive scope:** 12 complementary analyses provide systematic coverage
2. **Statistical rigor:** Framework-level FDR correction, bootstrap CI, multiple testing awareness
3. **Constructive solution:** CSSI provides actionable improvement, not just diagnosis
4. **Methodological honesty:** Limitations clearly acknowledged, negative results properly reported
5. **Practical guidelines:** Clear recommendations for practitioners
6. **Multi-modal validation:** Synthetic + real structured + real attention evidence

---

## CRITICAL LIMITATIONS (ACKNOWLEDGED)
1. **Missing positive controls:** Limits interpretation of negative findings
2. **Reference database circularity:** Affects all GRN evaluation approaches
3. **Sample size constraints:** Particularly mediation analysis (16 pairs)
4. **Circular validation:** Layer selection on same dataset used for performance reporting

---

## TOP 3 EXPERIMENTS FOR MAXIMUM IMPACT

If not ready to accept, these experiments would provide the biggest improvement:

1. **Independent CSSI validation:** Test CSSI on external dataset using layers identified here. This would break circularity and test generalization - the single most important validation gap.

2. **Synthetic validation with standard simulators:** Replicate scaling failure and CSSI effectiveness using SERGIO or other established generators to rule out custom generator bias.

3. **Expanded mediation analysis:** Scale to 100+ experiments with synthetic positive/negative controls to definitively characterize non-additivity prevalence.

---

## CONCLUSION

This paper represents a major contribution to single-cell mechanistic interpretability with rigorous experimental evidence supporting its central claims. The work has evolved substantially through 6 review rounds, and the recent improvements (CSSI cross-validation, expanded pseudotime analysis, baseline comparisons) address key methodological concerns.

**The experimental evidence is now adequate for publication.** While several analyses could be strengthened with additional data, the core findings are well-supported and the limitations are honestly acknowledged.

The paper establishes critical boundary conditions for interpreting foundation model attention as biological mechanism and provides both diagnostic tools and constructive solutions. The CSSI framework represents a genuine methodological advance with demonstrated effectiveness across multiple validation approaches.

**Recommendation: ACCEPT**

The field needs these boundary conditions and quality control standards. The work is ready for publication and will provide immediate value to practitioners while establishing rigorous evaluation standards for future mechanistic interpretability research.