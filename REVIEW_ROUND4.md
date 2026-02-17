# ROUND 4 ADVERSARIAL REVIEW: NMI Paper
## "A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models"

**Date**: February 14, 2026  
**Reviewer**: Subagent Round 4 Adversarial Review

---

## EXECUTIVE SUMMARY

After three comprehensive rounds of review and revision, this paper has substantially improved in statistical rigor, methodological transparency, and narrative coherence. The authors have successfully addressed most major criticisms from previous reviews, implementing framework-level FDR correction, adding proper baseline comparisons, and providing theoretical grounding for CSSI.

However, **several critical issues remain** that prevent acceptance at a top-tier venue like JMLR. Most significantly, the paper now suffers from **internal inconsistencies introduced during the revision process**, **persistent overclaims despite attempts at softening**, and **missing external validation** that a JMLR reviewer would demand.

**RECOMMENDATION: MAJOR REVISION REQUIRED**

---

## PART 1: INCONSISTENCIES INTRODUCED BY PREVIOUS FIXES

### **1.1 Contradictory Scaling Claims**

**CRITICAL INCONSISTENCY**: The paper now contains contradictory claims about scaling behavior that were introduced during attempts to reconcile previous reviews:

- **Lines 36-37** (Abstract): Claims "inverse scaling for interpretability in scGPT"
- **Lines 643-645** (Multi-model section): States "Geneformer's rank-based tokenization produces stable but equally uninformative attention patterns"
- **Lines 79** (Methods): Lists scaling points as "200, 500, 1000, 3000" 
- **Lines 715** (Limitations): States "only 3 points (200/1000/3000)"

**Problem**: The methods and limitations sections contradict each other on basic experimental design. This appears to be an artifact of multiple revisions where different sections were updated independently.

### **1.2 CSSI Effectiveness Contradictions**

**CRITICAL INCONSISTENCY**: The paper simultaneously claims CSSI works and doesn't work on real attention data:

- **Lines 535, 565, 567** (CSSI section): Claims "AUROC 0.694-0.706 in later layers of scGPT"
- **Lines 553-555, 563** (CSSI results): Shows "~zero or slightly negative CSSI delta" for strongest layers
- **Line 722** (Conclusions): Still claims CSSI "eliminates scaling failure" broadly

**Problem**: These are materially contradictory claims about the same experiments. The revision process appears to have left both positive and negative findings in the text without resolution.

### **1.3 Mixed Edge Score Methodologies**

**INTRODUCED CONFUSION**: Previous reviews correctly identified that the paper mixed attention-based and correlation-based analyses. The authors attempted to fix this but created new confusion:

- **Section 6.7** (Cross-species): Now explicitly states it uses "Spearman correlation-based edge scores" but is still framed as validating "attention-based GRN inference"
- **Section 6.8** (Pseudotime): Uses "expression correlations" but conclusions are applied to attention methods
- **Lines 676** (Discussion): States "CSSI-enhanced interpretability pipelines that stratify cells before any mechanistic analysis" without distinguishing correlation vs. attention

**Problem**: The attempted fix created a hybrid methodology that is neither purely attention-based nor purely correlation-based, making interpretation unclear.

---

## PART 2: STILL-OVERSTATED CLAIMS AFTER SOFTENING

### **2.1 "Fundamental" and "Architecture-Independent" Claims**

**OVERSTATED**: Despite attempts to soften language, the paper still makes broad architectural generalizations from limited evidence:

- **Line 36**: "architecture-independent near-random failure"
- **Line 645**: "attention-based regulatory inference fails across architectures"

**Reality check**: Evidence exists for only scGPT + Geneformer deeply, with scVI/C2S-Pythia explicitly acknowledged as "preliminary" (Line 641). This does not support "architecture-independent" or "fundamental" claims.

### **2.2 CSSI as "Constructive Solution" 

**OVERSTATED**: The paper continues to present CSSI as a major methodological breakthrough:

- **Line 722**: "CSSI eliminates scaling failure" 
- **Lines 136-157**: Extensive theoretical justification claiming CSSI is "far beyond naive partitioning"

**Reality check**: CSSI is cell-type stratification before analysis - a standard preprocessing step that any competent practitioner would consider. The "1.85× improvement" is on already near-random baselines (AUROC ~0.5), making the absolute improvement marginal.

### **2.3 "Comprehensive Framework" Claim**

**OVERSTATED**: The title and abstract claim a "comprehensive" and "systematic" framework:

- **Title**: "Systematic Framework for Mechanistic Interpretability"
- **Line 45**: "comprehensive framework comprising twelve complementary analyses"

**Reality check**: The framework has significant gaps (see Part 3) and several analyses use correlation-based rather than attention-based methods, limiting its relevance to foundation model interpretability specifically.

---

## PART 3: HELD-OUT VALIDATION INTEGRATION PROBLEMS

### **3.1 Missing True External Validation**

**CRITICAL GAP**: While the paper mentions plans for external validation, no true held-out validation exists:

- **Lines 703-715** (External Validation section): Discusses "planning validation on the Replogle et al. (2022)" but provides no actual results
- The paper relies on the same datasets (Tabula Sapiens, Dixit/Adamson/Shifrut) used throughout development

**JMLR expectation**: A methods paper of this scope requires demonstration on completely independent datasets, not just promises of future validation.

### **3.2 Layer/Head Selection Circularity Unresolved**

**CRITICAL GAP**: Previous reviews identified circular validation in CSSI layer selection. The current text acknowledges this but doesn't fix it:

- **Lines 535-567**: Still reports results on pre-selected "best layers" (L13-L14) identified from the same benchmark used for evaluation
- No nested cross-validation or truly held-out layer selection provided

**JMLR expectation**: Proper out-of-sample validation of the layer selection strategy, not just acknowledgment of the limitation.

---

## PART 4: STATISTICAL TEST REGISTRY COMPLETENESS

### **4.1 Registry vs. Text Discrepancies**

**PROBLEMS IDENTIFIED**: The statistical test registry (appendix_statistical_tests.tex) has several issues:

1. **Missing tests**: The registry lists 47 tests but several analyses in the main text are not registered (e.g., bootstrap CIs reported as ranges without specific tests)
2. **Inconsistent correction**: Some results report raw p-values in text but corrected values in registry
3. **Effect size inconsistencies**: Registry shows different effect sizes than main text for some analyses

**Specific example**: Pseudotime analysis reports "adjusted p = 0.124" in discussion but registry shows different values.

### **4.2 Double-Penalization Continuation**

**UNRESOLVED**: Previous Review Round 2 identified "double-penalization" in FDR correction (within-analysis + framework-level). The current version still shows evidence of this:

- **Lines 109-110**: Claims "framework-wide BH correction over 47 tests"
- **Individual sections**: Still apply within-analysis corrections (e.g., "No individual pair reached significance after FDR correction")

**Problem**: This conservative approach may be over-correcting, leading to inflated Type II error rates.

---

## PART 5: WRITING QUALITY ISSUES FROM EDITS

### **5.1 Repetitive Methodological Blocks**

**QUALITY ISSUE**: The paper contains repetitive "important methodological note" blocks that break flow:

- **Lines 136-157**: Extended CSSI justification repeats points made elsewhere
- **Lines 183-196**: Synthetic validation justification rehashes the same anti-SERGIO arguments from methods
- **Lines 676-690**: CSSI discussion repeats claims from earlier sections

**Impact**: These repetitive blocks make the paper feel bloated and hurt readability.

### **5.2 Inconsistent Terminology**

**QUALITY ISSUE**: Previous reviews noted terminology issues that remain:

- "Edge scores" sometimes means attention weights, sometimes correlations (inconsistent throughout)
- "GRN inference" conflates causal discovery with edge ranking
- "Validation" includes both reference comparison and experimental perturbation

**Impact**: Readers cannot be certain what specific method is being discussed in any given analysis.

### **5.3 Broken Flow from Piecemeal Edits**

**QUALITY ISSUE**: The revision process has created awkward transitions:

- **Lines 617-645** (Multi-model section): Abruptly introduces Geneformer without proper motivation
- **Lines 676-690** (Discussion): Jumps between unrelated points without clear logical progression
- **Lines 703-715** (External validation): Feels like an afterthought rather than integrated discussion

---

## PART 6: MISSING EXPERIMENTS FOR JMLR ACCEPTANCE

### **6.1 True Cross-Model Layer Analysis**

**MISSING**: The paper claims layer-dependent effects (L13-L14 for scGPT) but doesn't test whether this generalizes to Geneformer:

- **Required**: Layer-by-layer analysis of Geneformer attention patterns
- **Required**: Cross-model consistency of which layers contain regulatory signal
- **Current gap**: Geneformer analysis is limited to overall AUROC, not layer-specific patterns

### **6.2 Attention-Based Cross-Species Analysis**

**MISSING**: The cross-species analysis (Section 6.7) uses correlation-based edges, not attention-derived edges:

- **Required**: Actual attention-weight-based cross-species comparison
- **Required**: Testing whether attention-derived human→mouse transfer follows the same biological patterns as correlation-based transfer
- **Current gap**: All cross-species findings may not apply to attention-based methods

### **6.3 Nested CSSI Validation**

**MISSING**: Previous reviews identified circular CSSI validation that remains unfixed:

- **Required**: Split data → identify best layers on training set → apply CSSI on held-out test set
- **Required**: Cross-dataset validation (train layer selection on one tissue, test on another)
- **Current gap**: All CSSI results use the same data for layer selection and performance evaluation

### **6.4 Scaling Density Analysis**

**MISSING**: The scaling analysis uses only 3-4 cell count points, insufficient for characterizing scaling curves:

- **Required**: At least 8-10 cell count points from 200-5000 range
- **Required**: Identification of potential scaling regime transitions
- **Current gap**: Cannot distinguish monotonic scaling failure from more complex scaling relationships

### **6.5 Proper Negative Controls**

**MISSING**: The paper lacks proper negative controls to validate that observed patterns aren't statistical artifacts:

- **Required**: Apply CSSI to random attention matrices and show it doesn't improve performance
- **Required**: Test whether "good" layers (L13-L14) are consistently good across independent datasets
- **Current gap**: Cannot rule out that observed patterns are dataset-specific noise

---

## PART 7: SPECIFIC LINE-LEVEL ISSUES

### **7.1 Reference Inconsistencies**

- **Line 533**: "TRRUST/TRRust" - inconsistent spelling within same paragraph
- **Lines 79 vs 715**: Scaling points inconsistency mentioned above
- **Lines 553-555 vs 565-567**: Contradictory CSSI effectiveness claims

### **7.2 Statistical Reporting Issues**

- **Lines 617-620**: Reports AUROC values without confidence intervals that are mentioned elsewhere
- **Table 6**: Mentioned in text but confidence intervals not shown in actual table
- **Various sections**: Mix of raw and adjusted p-values without clear indication

### **7.3 Methods-Results Mismatches**

- **Methods Section 2.2**: Describes 4 scaling points, Results only use 3 in some analyses
- **Methods Section 2.11**: Claims CSSI has "rigorous theoretical grounding" but theory is just mathematical formalization of stratification
- **Methods vs Results**: Various parameter mismatches (HVG counts, bootstrap iterations, etc.)

---

## RECOMMENDATIONS FOR REVISION

### **ESSENTIAL FIXES (Required for Acceptance)**

1. **Resolve internal contradictions**: Carefully reconcile all contradictory claims about CSSI effectiveness, scaling behavior, and experimental parameters

2. **True external validation**: Provide actual results (not just plans) on independent datasets like Replogle et al. (2022)

3. **Fix CSSI circular validation**: Implement proper nested cross-validation for layer/head selection

4. **Attention-only analyses**: Either convert correlation-based analyses to attention-based or clearly scope claims to correlation methods only

5. **Expand model comparison**: Provide layer-level analysis for Geneformer and test cross-model generalization of layer effects

### **STRONGLY RECOMMENDED**

6. **Scale up scaling analysis**: Test 8-10 cell count points to properly characterize scaling relationships

7. **Add proper negative controls**: Test CSSI on random attention matrices, validate layer consistency across datasets

8. **Streamline writing**: Remove repetitive blocks, fix terminology consistency, improve logical flow

9. **Complete statistical registry**: Ensure all tests in appendix match main text reporting

### **MINOR BUT IMPORTANT**

10. **Fix citation inconsistencies**: TRRUST spelling, figure references, etc.
11. **Methods-results alignment**: Ensure all parameters match between sections
12. **Confidence interval display**: Show CIs in figures, not just text

---

## CONCLUSION

This paper has made substantial improvements through three rounds of revision and addresses an important problem in single-cell foundation model interpretability. The statistical corrections, baseline comparisons, and theoretical development of CSSI represent genuine contributions.

However, **the revision process has introduced new critical problems**: internal contradictions, unresolved circular validation, and missing external validation that would be expected for JMLR acceptance. The paper feels like it has been over-edited, with different sections providing conflicting information about the same experiments.

Most critically, the **lack of true external validation** is a deal-breaker for a methods paper of this scope. Claims about "architecture-independent" failures and CSSI effectiveness **must be demonstrated on independent held-out datasets**, not just planned.

The paper is much closer to acceptance than in earlier rounds, but requires one more careful revision to resolve the inconsistencies introduced during previous fixes and provide the external validation that JMLR reviewers will demand.

**RECOMMENDATION: MAJOR REVISION REQUIRED**

**Key message to authors**: You've addressed most of the original scientific issues successfully. The remaining problems are primarily about internal consistency and external validation - both fixable with careful attention to detail and additional computation rather than fundamental redesign.
