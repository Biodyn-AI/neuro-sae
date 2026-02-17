# ROUND 5 FINAL ADVERSARIAL REVIEW: NMI Paper
## "A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models: Attention-Based Gene Regulatory Network Inference and Cell-State Stratified Analysis"

**Date**: February 14, 2026  
**Final Adversarial Review - Round 5**

---

## EXECUTIVE SUMMARY: SUBSTANTIAL PROGRESS BUT CRITICAL GAPS REMAIN

After four comprehensive rounds of review and revision, this paper has undergone remarkable improvement in statistical rigor, experimental design, and theoretical grounding. The authors have successfully addressed most of the fundamental flaws identified in earlier rounds, implementing proper framework-level FDR correction, adding comprehensive baseline comparisons, and providing substantial theoretical justification for CSSI.

However, **several critical issues persist that prevent confident acceptance at a top-tier venue**. Most significantly, the paper suffers from **unresolved internal contradictions introduced during the revision process**, **methodological inconsistencies that blur the core contributions**, and **missing external validation** that limits confidence in generalizability.

**RECOMMENDATION: MINOR REVISION REQUIRED**

The paper is now very close to publication quality. The remaining issues are primarily about consistency, clarity, and completing the validation framework rather than fundamental scientific flaws.

---

## PART 1: REMARKABLE PROGRESS FROM PREVIOUS ROUNDS

### **1.1 Successfully Addressed Major Issues**

**✅ Statistical Rigor**: The framework-level Benjamini-Hochberg FDR correction across all 47 statistical tests (Section 2.3) is exemplary. The honest reporting of which findings survive correction demonstrates scientific integrity.

**✅ Baseline Integration**: The comprehensive baseline comparison (Section 5.2) that shows all methods achieve near-random performance on brain tissue fundamentally reframes the narrative from "attention methods fail" to "tissue-specific challenges exist." This is excellent scientific reasoning.

**✅ Multi-Model Validation**: The Geneformer analysis (Section 6.12) provides crucial evidence that near-random GRN recovery is not scGPT-specific but occurs across architecturally distinct models. This substantially strengthens the generalizability claims.

**✅ CSSI Theoretical Grounding**: While still fundamentally cell-type stratification, the mathematical formalization and systematic validation across synthetic, biologically-structured, and real attention data provides solid methodological foundation.

**✅ Narrative Coherence**: The transition from "attention fails universally" to "unstratified approaches fail to account for cellular heterogeneity" is well-supported and logically consistent with the evidence.

---

## PART 2: REMAINING CRITICAL ISSUES (MINOR BUT IMPORTANT)

### **2.1 Internal Contradictions from Revision Artifacts**

**ISSUE**: The paper contains contradictory statements about experimental parameters that appear to be artifacts of multiple revision cycles:

- **Lines 79** (Methods): Lists scaling points as "200, 500, 1000, 3000"
- **Lines 715** (Limitations): States "only 3 points (200/1000/3000)"

**Evidence**: The methods section clearly describes four scaling points, but limitations section claims only three were tested.

**Impact**: This creates confusion about the actual experimental design and appears to be a simple editing oversight.

**Fix**: Reconcile these statements - either clarify that different analyses used different point sets or ensure consistency.

### **2.2 CSSI Effectiveness Presentation Inconsistencies**

**ISSUE**: The paper presents seemingly contradictory information about CSSI performance on real attention data:

- **Claims strong performance**: "AUROC 0.694-0.706 in later layers of scGPT" (lines 535, 565, 567)
- **Shows minimal CSSI improvement**: "~zero or slightly negative CSSI delta" for strongest layers (Table 7, lines 553-555)
- **Broad conclusions**: "CSSI eliminates scaling failure" (line 722)

**Analysis**: These aren't actually contradictory - they reflect that CSSI's value is diagnostic (identifying good layers) rather than uniformly improving all layers. However, the presentation obscures this important distinction.

**Fix**: Clarify that CSSI's primary contribution is identifying which architectural components contain extractable regulatory signal, with improvements concentrated in intermediate layers.

### **2.3 Mixed Methodology Clarity**

**ISSUE**: Several key analyses use correlation-based rather than attention-based edge scores, but this is not consistently highlighted:

- **Cross-species analysis** (Section 6.7): Uses "Spearman correlation-based edge scores"  
- **Pseudotime analysis** (Section 6.8): Uses "expression correlations"
- **Batch leakage** (Section 6.9): Uses "Pearson correlations"

**Analysis**: The authors do acknowledge this limitation in methodological notes, but the broader implications for interpretation are not sufficiently emphasized.

**Fix**: Add a clear summary table or section distinguishing which analyses apply directly to attention-based interpretability vs. general GRN inference methodology.

---

## PART 3: MINOR CONSISTENCY AND PRESENTATION ISSUES

### **3.1 Reference and Spelling Inconsistencies**

**Minor but noticeable issues**:
- **Line 533**: "TRRUST/TRRust" - inconsistent spelling within same paragraph
- Various minor parameter mismatches between methods and results sections
- Inconsistent confidence interval reporting (sometimes in text, sometimes in figures)

**Fix**: Careful proofreading pass to ensure consistency throughout.

### **3.2 Statistical Reporting Clarity**

**ISSUE**: The paper mixes raw and adjusted p-values in places without clear indication:

- Some sections report framework-level adjusted p-values
- Others report within-analysis corrections
- Not always clear which correction level is being used

**Fix**: Consistently indicate whether p-values are raw, within-analysis corrected, or framework-level corrected.

---

## PART 4: ASSESSMENT OF CLAIMED CONTRIBUTIONS

### **4.1 CSSI as Methodological Contribution**

**Assessment**: While CSSI is fundamentally cell-type stratification, the systematic validation and mathematical framework do constitute a meaningful contribution:

- **Synthetic validation**: Demonstrates up to 1.85× improvement with known ground truth
- **Real attention matrices**: Shows diagnostic value in identifying regulatory-rich layers
- **Theoretical framework**: Provides mathematical basis for when stratification helps

**Verdict**: Legitimate methodological contribution, though not revolutionary.

### **4.2 Comprehensive Framework Assessment**

**Assessment**: The twelve complementary analyses do constitute a comprehensive evaluation framework:

- Covers scaling, bias, detectability, cross-context validation, perturbation consistency
- Includes proper statistical controls and multiple testing correction
- Integrates constructive solution (CSSI) with diagnostic findings

**Verdict**: The "comprehensive framework" claim is well-supported.

### **4.3 Architecture-Independent Claims**

**Assessment**: Evidence exists for scGPT and Geneformer deeply, with preliminary evidence for scVI and C2S-Pythia:

- Both major models show near-random GRN recovery (AUROC ≈ 0.5)
- Different scaling dynamics but similar outcomes
- Limited evidence for other architectures

**Verdict**: "Architecture-independent" is overstated but "consistent across tested architectures" is defensible.

---

## PART 5: EXTERNAL VALIDATION AND GENERALIZABILITY

### **5.1 Missing External Validation**

**REMAINING GAP**: The paper mentions plans for validation on Replogle et al. (2022) but provides no actual results (lines 703-715).

**Assessment**: While this limits immediate confidence in generalizability, the consistency across multiple datasets within the current framework (Tabula Sapiens, multiple Perturb-seq studies, cross-species data) provides reasonable evidence of robustness.

**Recommendation**: This should be addressed in future work rather than blocking publication.

### **5.2 Layer Selection Circularity**

**REMAINING CONCERN**: The CSSI evaluation still uses the same data for layer identification and performance reporting.

**Assessment**: The authors acknowledge this limitation and provide held-out validation that shows similar patterns despite ranking sensitivity. The biological interpretation (later layers encode regulatory signals) is plausible and consistent across analyses.

**Recommendation**: The acknowledgment of this limitation is sufficient for current publication.

---

## PART 6: BIOLOGICAL AND METHODOLOGICAL SOUNDNESS

### **6.1 Biological Interpretability**

**STRENGTH**: The paper's biological interpretations are sound:
- Cell-state-specific regulation is well-established
- Layer-dependent signal concentration makes biological sense
- Cross-species conservation patterns align with known regulatory biology

### **6.2 Statistical Rigor**

**STRENGTH**: The statistical approach is now exemplary:
- Proper multiple testing correction
- Honest reporting of negative results
- Comprehensive confidence intervals and effect size reporting
- Framework-level registration of all statistical tests

### **6.3 Methodological Transparency**

**STRENGTH**: The methods are detailed and reproducible:
- Complete code availability promises
- Detailed parameter specifications
- Honest discussion of limitations

---

## PART 7: COMPARISON TO PREVIOUS REVIEW ROUNDS

### **What Has Been Successfully Fixed**:
1. **Major statistical flaws** (Round 1) - Completely resolved
2. **Circular synthetic validation concerns** (Round 1) - Well-justified
3. **Narrative contradictions** (Rounds 1-2) - Largely resolved
4. **Cherry-picked baseline issues** (Round 1) - Transformed into strength
5. **Methodological inconsistency problems** (Round 2) - Mostly resolved
6. **Over-claimed contributions** (Rounds 2-3) - Significantly softened

### **What Remains to be Polished**:
1. **Internal parameter inconsistencies** - Minor editing issues
2. **CSSI presentation clarity** - Need clearer framing of diagnostic vs. improvement value
3. **Mixed methodology emphasis** - Need clearer delineation of attention vs. correlation analyses

---

## RECOMMENDATIONS FOR FINAL REVISION

### **ESSENTIAL FIXES (Required)**

1. **Resolve parameter inconsistencies**: Ensure methods and limitations sections agree on experimental design details (scaling points, sample sizes, etc.)

2. **Clarify CSSI contribution**: Explicitly state that CSSI's primary value is diagnostic (identifying regulatory-rich layers) with improvements concentrated in intermediate layers, not uniform improvement across all layers

3. **Methodology delineation**: Add a summary table or clear section distinguishing analyses that directly test attention-based interpretability vs. those testing general GRN inference principles

### **STRONGLY RECOMMENDED**

4. **Consistent statistical reporting**: Clearly indicate correction level for all p-values (raw, within-analysis, framework-level)

5. **Proofreading pass**: Fix spelling inconsistencies (TRRUST/TRRust) and ensure parameter consistency throughout

6. **External validation roadmap**: If Replogle validation isn't ready, provide clear timeline and specific hypotheses to be tested

### **MINOR IMPROVEMENTS**

7. **Figure consistency**: Ensure all confidence intervals are shown in figures rather than just mentioned in text
8. **Reference formatting**: Check citation consistency throughout
9. **Terminology glossary**: Consider adding a brief glossary of key terms to improve accessibility

---

## CONCLUSION: A STRONG PAPER NEEDING FINAL POLISH

This paper represents a substantial contribution to single-cell foundation model interpretability. After four rounds of rigorous revision, it has evolved from a deeply flawed initial submission to a methodologically sound, statistically rigorous, and scientifically valuable framework.

**The core scientific contributions are solid:**
- Systematic documentation of attention-based GRN inference limitations
- Cell-State Stratified Interpretability as a constructive methodological improvement
- Comprehensive evaluation framework with proper statistical controls
- Multi-model validation showing generalizability of findings

**The remaining issues are primarily cosmetic:**
- Internal consistency problems from multiple revision cycles
- Presentation clarity around CSSI's diagnostic vs. improvement functions
- Minor methodological delineation needs

**The paper is very close to acceptance**. The scientific soundness is now strong, the statistical rigor is exemplary, and the contributions are meaningful. The remaining fixes are primarily editorial rather than requiring new experiments or fundamental restructuring.

**FINAL VERDICT: MINOR REVISION REQUIRED**

The paper should be accepted after addressing the internal consistency issues and clarifying the presentation of CSSI's contribution. This represents genuine progress toward establishing quality-control standards for mechanistic interpretability in computational biology.

---

## CONFIDENCE ASSESSMENT

**High confidence in recommendation**: The paper has been thoroughly reviewed across four rounds by multiple perspectives. The remaining issues are well-identified and addressable through careful editing rather than fundamental changes.

**Publication readiness**: 90% - Very close to acceptance pending minor revisions.

**Scientific contribution value**: High - Establishes important boundary conditions and methodological improvements for an emerging field.

**Reproducibility**: High - Detailed methods, promised code availability, comprehensive parameter specifications.

**Impact potential**: High - Likely to influence best practices for single-cell foundation model interpretability and establish evaluation standards for the field.