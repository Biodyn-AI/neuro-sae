# ROUND 8 ADVERSARIAL RESEARCH REVIEW
## NMI Paper: Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models

**Reviewer Position**: HOSTILE - Attempting to REJECT this paper
**Date**: February 14, 2026
**Review Target**: 12+ analyses spanning scaling, bias, detectability, cross-context validation, perturbation consistency, ortholog transfer, pseudotime, batch leakage, calibration, CSSI, synthetic validation, and multi-model comparison

---

## EXECUTIVE SUMMARY

While this paper attempts a comprehensive framework for evaluating mechanistic interpretability in single-cell foundation models, it suffers from fundamental methodological flaws, circular reasoning, inadequate statistical power, and overclaims that undermine its scientific validity. The authors present what appears to be a rigorous multi-faceted analysis but fail to address critical confounds, lack proper positive controls, and make causal claims unsupported by their correlational evidence.

**RECOMMENDATION: REJECT**

---

## DETAILED ANALYSIS BY SECTION

### 1. SCALING BEHAVIOR ANALYSIS - Evidence Quality: C

**Current Evidence**: Monotonic degradation in F1 scores from 200→3000 cells across model tiers, with statistical significance (p=0.002).

**Critical Weakness**: **Overfitting Alternative Not Adequately Addressed**
The paper mentions but does not rigorously test the alternative explanation that superior performance at smaller cell counts reflects overfitting to sparse reference databases rather than genuine regulatory signal recovery. This is a fatal flaw.

**Alternative Explanation**: The "inverse scaling" could be **improved generalization** rather than signal loss. As cell counts increase, models may learn more robust representations that don't conform to potentially incomplete TRRUST/DoRothEA annotations. The authors acknowledge this but dismiss it without proper testing.

**New Experiment Needed**: **Reference Database Saturation Analysis**
- Generate synthetic networks with KNOWN ground truth at different sparsity levels (0.001%, 0.01%, 0.1% true edges)
- Test whether scaling failure persists when ground truth coverage is complete
- Expected output: If overfitting explanation is correct, synthetic complete networks should show improved scaling
- **Script needed**: `experiments/synthetic_saturation_test.py` - generate complete ground truth networks at varying sparsity, test scGPT scaling behavior

**Specific Implementation**: Create 5 synthetic networks with identical statistical properties but different reference database coverage (10%, 25%, 50%, 75%, 100% of true edges annotated). If scaling failure only occurs with incomplete coverage, this suggests overfitting rather than genuine signal loss.

---

### 2. MEDIATION BIAS ANALYSIS - Evidence Quality: D

**Current Evidence**: 62.5% non-additivity rate (10/16 run-pairs), but severely underpowered.

**Critical Weakness**: **Sample Size Catastrophically Small**
With only 16 run-pairs, this analysis has ~20% power to detect even large effects. The authors acknowledge this but still draw strong conclusions.

**Alternative Explanation**: Observed non-additivity could be **measurement noise** rather than biological interaction structure. With such small samples, random measurement error could easily produce apparent non-additivity.

**New Experiment Needed**: **Large-Scale Mediation Archive Construction**
- Expand to 200+ run-pairs across 10+ tissues and 5+ model architectures  
- Include synthetic positive controls with known additive vs non-additive structure
- Expected output: Definitive characterization of non-additivity prevalence with adequate statistical power
- **Script needed**: `experiments/mediation_expansion.py` - automated pipeline for large-scale mediation analysis across Tabula Sapiens
- **Data needed**: Tabula Sapiens full atlas, multiple scGPT checkpoints

**Power Calculation**: To detect medium effect size (d=0.5) at 80% power requires n≥128 pairs. Current n=16 provides <30% power.

---

### 3. DETECTABILITY THEORY - Evidence Quality: B

**Current Evidence**: Mathematical framework with phase diagrams showing intervention vs attention signal advantages.

**Critical Weakness**: **Theoretical Framework Lacks Experimental Validation**
The detectability equations assume specific noise models and effect sizes that may not reflect biological reality. No direct comparison with real foundation model components.

**Alternative Explanation**: The theoretical advantages may **disappear under biological constraints** not captured in the simplified mathematical framework (e.g., network topology constraints, gene expression constraints).

**New Experiment Needed**: **Real Foundation Model Calibration**
- Extract actual effect sizes and noise parameters from scGPT attention matrices across 20+ tissues
- Compare empirical detectability with theoretical predictions  
- Test whether theoretical sample complexity estimates hold for real biological networks
- Expected output: Calibrated detectability framework with biologically realistic parameters
- **Script needed**: `experiments/detectability_calibration.py` - extract noise parameters from real scGPT runs, validate theoretical predictions

---

### 4. CROSS-CONTEXT CONSISTENCY - Evidence Quality: C

**Current Evidence**: Variable cross-tissue correlations (-0.44 to 0.71), only 2/6 comparisons survive FDR correction.

**Critical Weakness**: **Batch Effects vs Biology Confound Unresolved**
The authors acknowledge but don't adequately address that technical batch effects between tissues could explain low consistency rather than genuine biological differences.

**Alternative Explanation**: **Technical artifacts** (different protocols, processing pipelines, sequencing depths) create spurious tissue-specific patterns that manifest as biological regulatory differences.

**New Experiment Needed**: **Batch-Controlled Cross-Tissue Analysis**
- Process samples from same donors across multiple tissues using identical protocols
- Apply computational batch correction (Harmony, scVI) to existing data and re-evaluate consistency
- Include technical replicates to separate biological vs technical variance
- Expected output: Quantification of batch vs biological contributions to cross-tissue inconsistency
- **Script needed**: `experiments/batch_controlled_cross_tissue.py` - implement Harmony/scVI correction, re-run consistency analysis
- **Data needed**: Tabula Sapiens with batch metadata, technical replicates where available

---

### 5. PERTURBATION VALIDATION - Evidence Quality: D

**Current Evidence**: No results survive framework-level multiple testing correction. All findings are exploratory.

**Critical Weakness**: **Multiple Testing Overcorrection**
Applying FDR correction across all 47 tests in the entire paper dilutes power for individual analyses. This is statistically conservative to the point of being uninformative.

**Alternative Explanation**: The correction may be **inappropriately stringent** - perturbation validation should be treated as an independent hypothesis with its own error rate.

**New Experiment Needed**: **Focused Perturbation Study**
- Design standalone perturbation validation study with appropriate sample size
- Pre-register specific hypotheses with dedicated statistical analysis plan
- Include positive controls (known regulatory relationships) and negative controls (non-regulatory gene pairs)
- Expected output: Definitive assessment of attention-perturbation alignment with appropriate statistical power
- **Script needed**: `experiments/focused_perturb_validation.py` - power analysis, sample size calculation, standalone validation framework

**Implementation**: Target 20 well-characterized TF-target pairs with both attention scores and perturbation outcomes. Power analysis suggests n≥50 perturbations needed for 80% power.

---

### 6. CROSS-SPECIES ORTHOLOG TRANSFER - Evidence Quality: B+

**Current Evidence**: Strong global conservation (ρ=0.743) but high per-TF variability.

**Critical Weakness**: **Single Tissue, Single Species Pair Limits Generalizability**
Only human-mouse lung tested. Unknown whether findings generalize to other tissues/species.

**Alternative Explanation**: **Lung-specific or mammalian-specific patterns** may not represent general principles of regulatory conservation.

**New Experiment Needed**: **Multi-Species, Multi-Tissue Conservation Atlas**
- Extend to 5+ tissue types (immune, liver, kidney, brain, heart) 
- Include non-mammalian species (zebrafish, Drosophila) where available
- Test evolutionary distance effects on conservation patterns
- Expected output: General principles of regulatory conservation across phylogenetic distances
- **Script needed**: `experiments/multi_species_conservation.py` - automated ortholog mapping, conservation analysis across species/tissues
- **Data needed**: Cross-species single-cell atlases (human, mouse, zebrafish), ortholog databases

---

### 7. PSEUDOTIME DIRECTIONALITY AUDIT - Evidence Quality: C

**Current Evidence**: Only 21.5% directional consistency (31/144 pairs), but statistical significance lost after FDR correction (p=0.124).

**Critical Weakness**: **Single Pseudotime Algorithm Tested**
Only diffusion pseudotime used. Alternative algorithms (Monocle3, PAGA, RNA velocity) might show different results.

**Alternative Explanation**: **Methodological limitation** of pseudotime rather than biological invalidity of regulatory relationships.

**New Experiment Needed**: **Multi-Algorithm Temporal Validation**
- Test 5+ pseudotime methods (DPT, Monocle3, PAGA, Slingshot, RNA velocity)
- Include time-course perturbation data as positive control
- Test synthetic trajectories with known temporal ordering
- Expected output: Method-specific vs universal temporal validation limitations
- **Script needed**: `experiments/multi_pseudotime_validation.py` - implement multiple pseudotime algorithms, synthetic trajectory testing
- **Data needed**: Time-course perturbation datasets, trajectory simulation tools

---

### 8. BATCH AND DONOR LEAKAGE AUDIT - Evidence Quality: B

**Current Evidence**: High leakage in imbalanced datasets (54.6% flagged edges), stable in balanced datasets (10.1% flagged).

**Critical Weakness**: **Practical Impact Assessment Missing**
High classifier accuracy doesn't necessarily mean biological conclusions are invalid if leakage is biologically meaningful (donor-specific biology).

**Alternative Explanation**: **Genuine biological variation** between donors (genetic differences, environmental exposures, disease states) rather than pure technical artifacts.

**New Experiment Needed**: **Biological vs Technical Leakage Decomposition**
- Include genetic variants (eQTLs) to separate genetic from technical donor effects
- Test identical twins or clonal populations to isolate technical variation
- Validate key regulatory findings in independent donor cohorts
- Expected output: Quantification of biological vs technical contributions to donor leakage
- **Script needed**: `experiments/bio_tech_leakage_decomp.py` - eQTL integration, genetic variant analysis
- **Data needed**: Genotype data for Tabula Sapiens donors, population genetics databases

---

### 9. UNCERTAINTY CALIBRATION - Evidence Quality: B

**Current Evidence**: Severe miscalibration (ECE 0.27-0.47) improved 4-7× by post-hoc calibration.

**Critical Weakness**: **Calibration Doesn't Transfer Between Datasets**
K562-trained calibrators fail on T-cell data (ECE 0.32-0.42), limiting practical utility.

**Alternative Explanation**: **Fundamental unsuitability** of edge scores as probabilities, even with calibration.

**New Experiment Needed**: **Universal Calibration Framework**
- Develop meta-calibration approaches that generalize across cell types
- Test whether feature-based calibration (using cell-type, tissue metadata) improves transfer
- Compare with Bayesian approaches that model uncertainty explicitly
- Expected output: Transferable calibration framework for regulatory edge scores
- **Script needed**: `experiments/universal_calibration.py` - meta-learning calibration, feature-based approaches
- **Data needed**: Multiple cell types with perturbation ground truth for training

---

### 10. CELL-STATE STRATIFIED INTERPRETABILITY (CSSI) - Evidence Quality: B-

**Current Evidence**: Up to 1.85× improvement in synthetic validation, modest real-data gains.

**Critical Weaknesses**: 
1. **Circular Validation Problem**: Layer selection performed on full dataset, then performance reported on same data
2. **Limited Real-Data Improvement**: Best layers already perform well (AUROC 0.69-0.71), CSSI adds minimal value
3. **Implementation Complexity vs Benefit Ratio**: Requires cell-type annotation, clustering, may be impractical

**Alternative Explanations**: 
- **Confound control** rather than signal enhancement - CSSI might just be sophisticated batch correction
- **Cherry-picked examples** - synthetic validation may be optimized for CSSI success

**New Experiments Needed**:

**A. Independent Validation Framework**
- Hold out 50% of data for layer selection, evaluate CSSI on remaining 50%
- Cross-validate across different tissues for layer transferability  
- Expected output: Unbiased estimate of CSSI improvement
- **Script needed**: `experiments/cssi_holdout_validation.py` - proper train/test split for layer selection

**B. CSSI vs Alternative Controls**  
- Compare CSSI against sophisticated batch correction (Harmony, scVI)
- Test whether simple k-means clustering achieves similar benefits to biological annotation
- Include negative controls (random cell groupings)
- Expected output: Specific value of biological stratification vs general clustering
- **Script needed**: `experiments/cssi_control_comparison.py` - systematic comparison with alternative approaches

**Implementation Priority**: This is the paper's main constructive contribution but needs rigorous validation to avoid overclaims.

---

### 11. SYNTHETIC GROUND-TRUTH VALIDATION - Evidence Quality: C

**Current Evidence**: Validates theoretical predictions on custom-generated data.

**Critical Weakness**: **Custom Generator Bias**
Authors acknowledge the generator encodes assumptions aligned with their predictions, limiting independent validation value.

**Alternative Explanation**: **Circular validation** - synthetic results confirm theoretical assumptions rather than biological reality.

**New Experiment Needed**: **Independent Synthetic Validation**
- Use established generators (SERGIO, scDesign3) as external validation
- Compare predictions across multiple simulation frameworks
- Include adversarial examples designed to break theoretical predictions
- Expected output: Theory validation on independently generated data
- **Script needed**: `experiments/independent_synthetic_validation.py` - SERGIO integration, scDesign3 comparison
- **Data needed**: SERGIO framework, scDesign3 installation

---

### 12. MULTI-MODEL GENERALIZATION - Evidence Quality: B+

**Current Evidence**: Both scGPT and Geneformer show near-random GRN recovery (AUROC ≈0.5).

**Critical Weakness**: **Limited Architecture Coverage**
Only 2 foundation models tested (+ brief scVI, C2S-Pythia experiments). Conclusions about "foundation models generally" are overstated.

**Alternative Explanation**: **Tissue-specific failure** rather than architecture-general failure (both tested on challenging brain tissue).

**New Experiment Needed**: **Comprehensive Multi-Model, Multi-Tissue Analysis**
- Test 5+ architectures (scFoundation, scBERT, Geneformer variants, scGPT variants)
- Include 10+ tissue types with varying regulatory complexity
- Test whether any architecture×tissue combination succeeds
- Expected output: Architecture-tissue performance matrix identifying optimal combinations
- **Script needed**: `experiments/comprehensive_multimodel.py` - standardized evaluation across models/tissues
- **Data needed**: Multiple foundation model checkpoints, tissue-matched ground truth

---

## TOP 5 CRITICAL IMPROVEMENTS RANKED BY IMPACT/EFFORT

### 1. **Reference Database Saturation Analysis** (Impact: HIGH, Effort: MEDIUM)
- **Why Critical**: Directly addresses the core scaling failure alternative explanation
- **Implementation**: Synthetic networks with varying ground truth coverage
- **Expected Timeline**: 2-3 weeks
- **Scripts needed**: `experiments/synthetic_saturation_test.py`
- **Impact**: Could completely change interpretation of scaling results

### 2. **CSSI Holdout Validation** (Impact: HIGH, Effort: LOW)  
- **Why Critical**: Fixes circular validation in main constructive contribution
- **Implementation**: Proper train/test split for layer selection
- **Expected Timeline**: 1 week
- **Scripts needed**: `experiments/cssi_holdout_validation.py`
- **Impact**: Essential for validating paper's primary methodological claim

### 3. **Focused Perturbation Validation Study** (Impact: HIGH, Effort: MEDIUM)
- **Why Critical**: Current perturbation results are statistically uninformative  
- **Implementation**: Standalone study with proper power analysis
- **Expected Timeline**: 3-4 weeks (requires new data collection/processing)
- **Scripts needed**: `experiments/focused_perturb_validation.py`
- **Impact**: Could provide definitive validation or invalidation of attention-perturbation alignment

### 4. **Batch-Controlled Cross-Tissue Analysis** (Impact: MEDIUM, Effort: LOW)
- **Why Critical**: Resolves major confound in cross-tissue consistency findings
- **Implementation**: Apply Harmony/scVI correction to existing data
- **Expected Timeline**: 1-2 weeks  
- **Scripts needed**: `experiments/batch_controlled_cross_tissue.py`
- **Impact**: Clarifies biological vs technical contributions to inconsistency

### 5. **Large-Scale Mediation Archive** (Impact: MEDIUM, Effort: HIGH)
- **Why Critical**: Current mediation analysis is severely underpowered
- **Implementation**: Expand to 200+ run-pairs across tissues
- **Expected Timeline**: 4-6 weeks
- **Scripts needed**: `experiments/mediation_expansion.py`  
- **Impact**: Provides definitive characterization of non-additivity prevalence

---

## ADDITIONAL METHODOLOGY CONCERNS

### Statistical Power Issues
- Multiple analyses underpowered (mediation n=16, perturbation corrected to non-significance)
- Framework-level FDR correction may be overly conservative
- Post-hoc power analysis needed for all major claims

### Missing Positive Controls
- No synthetic networks with known attention-regulatory relationships
- No time-course data with established temporal ordering
- No validated regulatory relationships as positive controls for each tissue

### Circular Reasoning Patterns
- Reference databases derived from correlational methods, validated with correlational approaches
- CSSI layer selection on full dataset, performance reported on same data
- Theoretical framework validated on synthetic data generated to match theory

### Reproducibility Concerns  
- Many analyses rely on specific random seeds or data splits
- Cross-validation results show sensitivity to data composition  
- Limited documentation of hyperparameter sensitivity

---

## VERDICT

While this paper attempts a comprehensive evaluation framework, it suffers from fundamental methodological flaws that undermine its conclusions:

1. **Core scaling failure results** may reflect overfitting rather than genuine signal loss
2. **Primary constructive contribution (CSSI)** uses circular validation  
3. **Statistical power** inadequate for many key claims
4. **Missing positive controls** limit interpretation of negative findings
5. **Alternative explanations** inadequately addressed throughout

The paper requires major revision addressing the top 5 critical improvements before it can support its strong claims about mechanistic interpretability limitations and solutions.

**RECOMMENDED DECISION: MAJOR REVISION or REJECT**

*Note: This adversarial review intentionally adopts a hostile stance to identify maximum weaknesses. Many findings may still be valuable despite these limitations, but the methodological concerns require addressing for publication in a high-impact venue.*