# Adversarial Review (Codex R11) — Research Quality Focus

**Manuscript:** *Attention in Single-Cell Foundation Models Encodes Co-Expression, Not Regulation: A Diagnostic Framework with Cross-Tissue Validation*

## Summary of what the paper claims
The paper argues (i) attention-derived gene–gene edges from scGPT/Geneformer track **co-expression** rather than **regulation** (TRRUST/DoRothEA), explaining AUROC≈0.5; (ii) a layer-/state-stratified diagnostic (CSSI) can recover regulatory signal (scGPT-18L late layers AUROC≈0.70 vs pooled ≈0.54); and (iii) this “signal recovery” replicates across tissues (scGPT-12L brain/kidney/whole-human AUROC≈0.72 TRRUST; ≈0.84 DoRothEA) with permutation p<1e-4.

The topic is timely (interpretability claims around attention) and the paper is unusually explicit about caveats in several places. However, in its current form the empirical evidence is **not strong enough** to support the headline framing at a top computational biology venue.

## Overall evaluation (top-venue bar)
### Would it be accepted at NMI / Nature Methods / Genome Biology?
**Unlikely in current form.** The central empirical “cross-tissue replication” and some “core” conclusions rely on *extremely small* labeled edge sets (27–28 TRRUST edges) and evaluation designs that are vulnerable to **circularity**, **selection effects**, and **benchmark mismatch**. There are promising ideas (layer diagnostics; caution about attention-as-explanation), but the experimental validation is not at “Methods / GB” rigor yet.

**Score (1–10): 5.5/10**
- Strong: clear problem statement, multi-pronged auditing mindset, admits limitations.
- Weak: key quantitative claims hinge on underpowered / low-coverage benchmarks; several analyses are not properly controlled; statistics sometimes give an illusion of rigor (bootstrap CIs, FDR across 47 tests) while the dominant uncertainty is *ground-truth incompleteness and evaluation set construction*.

---

## 1) Are experimental results properly powered and statistically sound?

### Major power/statistics issues
1. **Cross-tissue “replication” uses 27–28 TRRUST edges** (and 123 DoRothEA). With AUROC, such small positive sets can produce seemingly stable means while being dominated by a handful of TFs/targets. The reported tightness across tissues (0.716±0.002) is not meaningful when the edge labels are tiny and likely biased.
   - The paper acknowledges the limitation, but still treats it as “strongest evidence” that signal is real.
   - Bootstrap CIs over the *edge set* are reported (±0.03–0.05), but the main narrative still leans on point estimates.

2. **Multiple testing correction applied “framework-wide” (47 tests)** is not a substitute for *proper experimental design*. It may even be counterproductive: you combine heterogeneous hypotheses (scaling, pseudotime, batch leakage, etc.) then apply BH globally, which makes some biologically central validations (perturbation alignment) “non-significant” by construction.
   - For a methods/interpretability paper, what matters is preregistered primary endpoints (few) + properly controlled secondary analyses, not a giant omnibus FDR across everything.

3. **Scaling analysis statistical support is weak as written.** You report 2 repeats per cell-count and 50-iteration bootstrap CIs; then you state non-monotonicity (750→1000 decline) even though “all 95% bootstrap confidence intervals overlap” and adjacent comparisons are not significant.
   - Curve fitting (R²=0.90 exponential saturation) on 9 points with small replicates is descriptive, not inferential.
   - Also, the evaluated positives are “161 per condition on average” (unclear why it varies); that itself implies unstable AUROC estimates.

4. **Mediation non-additivity is explicitly underpowered (N=16 run-pairs)** yet highlighted as a key corrected-significant finding in the Discussion (“p=0.003”). With N=16, p-values are extremely sensitive to distributional assumptions and to how run-pairs were constructed.

5. **Pseudotime directionality audit**: you mix counts inconsistently (e.g., “19 of 75 pairs analyzed” but earlier “30 hematopoietic pairs”; totals don’t reconcile cleanly) and interpret a failure rate as an inherent limitation of pseudotime rather than (a) choice of TF/targets measured at mRNA not activity, (b) pseudotime binning/lag choices, (c) lineage topology mismatch.

### Statistical reporting gaps
- No clear statement of **unit of replication** for key AUROC values (edges? cell resamples? donors? independent datasets?). Many “tight” uncertainties appear to be within-dataset resampling, not across biological replicates.
- AUROC is emphasized even when positives are tiny; AUROC can be inflated by class imbalance and is insensitive to precision at the top of the ranking.

**Bottom line:** Many analyses are numerically careful (bootstraps, permutation tests), but the dominant limitation is *benchmark coverage* and *independence*, not lack of p-values.

---

## 2) Are there circular validations or data leakage?

### High-risk circularity / leakage points
1. **Layer selection circularity in scGPT-18L is real and central.** The paper states that layer identification and performance reporting use the same dataset and calls this “exploratory,” but the headline CSSI/scGPT-18L AUROC≈0.70 is still used as a core contribution.
   - The split-half attempt (248/249) is a minimal mitigation but still within one dataset and does not protect against dataset-specific biases.

2. **Reference database circularity** is acknowledged but under-addressed. TRRUST/DoRothEA edges are literature-derived and enriched for well-studied TFs; DoRothEA includes ChIP-seq binding which is *not regulation* and can be correlated with expression/co-expression. Since your core claim is attention≈co-expression, validation against partially co-expression-derived resources risks self-fulfilling “above chance” AUROC.

3. **Cross-tissue replication with scGPT-12L almost certainly shares a common selection bottleneck.** The fact that every layer in every tissue gives ~0.72 (TRRUST) and ~0.84 (DoRothEA) is suspiciously uniform—suggesting:
   - either the evaluation set is dominated by constitutive edges that are easy to rank by simple expression frequency/variance;
   - or the scoring method effectively collapses to a per-gene marginal statistic (e.g., gene frequency / token rank effects) that is invariant across layers.

4. **Donor/batch leakage audit is done on correlation-based edges**, not on the attention-derived edges that the paper is about. This is framed as “boundary conditions,” but the paper then uses it to make recommendations about mechanistic interpretability generally.

5. **Synthetic generator encodes the desired mapping**: you explicitly generate synthetic attention matrices as `tanh(A_true + structured noise + expression-bias noise)` and then show CSSI and Shapley win. That is not an independent validation; it is a positive control at best. Currently it functions as a *circular* support for the theoretical story.

---

## 3) Do claims match evidence strength?

### Overstated / mismatched claims
1. **Title and Abstract overreach** given benchmark limitations.
   - “Attention encodes co-expression, not regulation” is too categorical. Your evidence is primarily correlation to TRRUST edges and AUROC≈0.5 on heterogeneous tissue. That supports: *pooled attention edges are dominated by co-expression and do not recover curated regulatory edges in this setup.* It does **not** establish that attention cannot encode regulatory relationships in principle or in other tasks/finetunes.

2. **“Cross-tissue replication confirms signal is real”** is not justified with 27–28 TRRUST edges. At best: “across tissues, performance on a small, shared, constitutive edge set is consistent.” That is a very different claim.

3. **GRN failure attribution**: you argue near-random AUROC is “intrinsic to regulatory signal recovery from heterogeneous tissue,” but the baseline comparison uses (i) only 500 cells, (ii) only top 500 genes, (iii) generic methods (GENIE3/GRNBoost2) without context-specific tuning, and (iv) TRRUST/DoRothEA which may have poor brain coverage. This does not establish an intrinsic limit; it establishes a benchmark mismatch + design choice.

4. **Scaling plateau / decline**: you frame heterogeneity-driven attention dilution as explanation. But you do not directly test heterogeneity as a mediator in the real scaling experiment (e.g., fixed composition vs increasing heterogeneity with cell count). The synthetic experiments do, but they’re constructed.

5. **“Later layers encode regulation”**: For scGPT-18L, AUROC improves in deep layers on one dataset with circular selection. For scGPT-12L, all layers are equal. A more accurate claim is: “layer dependence is architecture- and dataset-dependent; layer stratification is necessary as a diagnostic.”

### Under-claimed / unclear contributions
- The paper’s strongest contribution might actually be the **negative result / cautionary audit framework** rather than CSSI performance. But the current framing tries to sell a “recovery” story that isn’t fully supported.

---

## 4) Alternative explanations not addressed (or not ruled out)

1. **Token frequency / gene prevalence confounds in transformers.** Attention patterns can strongly reflect token frequency, embedding norms, and positional biases. If “edge scores” correlate with co-occurrence, they may also correlate with *marginal gene frequency* and HVG selection artifacts. This is not audited.

2. **Evaluation set bias toward highly expressed TFs/targets.** With only 27 TRRUST edges, AUROC could be driven by a simple heuristic: edges where both genes are consistently detected / high variance. You mention this as a concern but don’t test it.
   - A simple baseline like `score(TF,target)=mean(expr(TF))*mean(expr(target))` or detection rate product might match your AUROC; that would demolish the “regulatory signal” claim.

3. **Cell-type composition differences drive both co-expression and “regulation” labels.** If TRRUST edges are enriched for cell-identity programs (e.g., STAT3/MYC/TP53), AUROC may reflect cell-type mixture rather than regulatory activity.

4. **DoRothEA AUROC≈0.84 could be “binding accessibility / housekeeping”** rather than regulation. You note ChIP-seq vs functional regulation, but don’t do a stratified analysis (e.g., DoRothEA confidence levels, tissue relevance, stimulatory vs inhibitory edges).

5. **Attention extraction details** (normalization across heads/layers, symmetrization, directionality) could dominate results. The paper does not give enough detail for a reviewer to judge whether the edge score is a biologically meaningful directed measure or just an undirected similarity.

---

## 5) Actionable fixes (ranked by impact)

### Tier 1 (highest impact; needed for top-venue acceptance)
1. **Replace the scGPT-12L “cross-tissue replication” benchmark with a non-toy evaluation set.**
   - 27–28 TRRUST edges is not publishable as primary evidence.
   - Options: (a) expand gene set / mapping to obtain hundreds–thousands of edges per tissue; (b) evaluate at the TF level (per-TF AUROC) with enough targets; (c) use additional resources (OmniPath TF→target, ChEA, RegNetwork) and show consistent patterns.
   - Report performance with **confidence intervals across independent donors/datasets**, not just edge bootstraps.

2. **Pre-register a primary endpoint and do proper held-out validation for layer selection.**
   - Split by donors (or by dataset/tissue) for selecting “best layer” and test on held-out donors/datasets.
   - Alternatively: choose layer based on an unsupervised criterion (e.g., minimal correlation with co-expression) and then evaluate on TRRUST to avoid label leakage.

3. **Add “dumb” baselines that test the constitutive-expression confound.**
   - Examples: detection-rate product, mean expression product, HVG rank-based score, PCA loading co-membership.
   - If these baselines match AUROC≈0.72/0.84, your “regulatory signal” interpretation collapses; you must then reframe as “benchmark bias.”

4. **Directly test whether CSSI improves *generalization* rather than within-dataset fit.**
   - Evaluate CSSI-selected layers/strata on an external dataset (different atlas or cohort) with similar cell types.
   - If unavailable, do leave-one-cell-type-out or leave-one-donor-out layer selection and report degradation.

5. **Clarify and standardize the evaluation universe.**
   - Define the candidate edge set (all TF–target pairs among HVGs? all gene–gene pairs?) and how negatives are sampled.
   - Report AUROC/AUPRC with a fixed universe; avoid “161 evaluated per condition on average” ambiguity.

### Tier 2 (strong improvements; address key skepticism)
6. **Quantify how much of attention-edge variance is explained by marginal gene statistics.**
   - Regress edge scores on TF frequency, target frequency, detection rates, and co-expression; report partial correlations.
   - Show that “regulatory layers” retain signal beyond these confounds.

7. **Demonstrate directedness (or explicitly drop direction claims).**
   - If attention-derived edges are undirected, stop calling it regulation unless you validate direction with perturbations/time.
   - If directed: show that TF→target differs from target→TF and improves against signed/causal benchmarks.

8. **Strengthen perturbation validation with properly powered datasets.**
   - You mention Replogle (>9k targets). This is exactly what is needed.
   - Make one perturbation atlas a central validation, not a future-work note.

9. **Make the synthetic generator a true positive control rather than a supporting pillar.**
   - Add alternative synthetic generators not aligned with your assumptions, including SERGIO for expression and then pass those expressions through the real model to generate attention weights (even if noisy). The current justification for not using SERGIO reads like convenience.

### Tier 3 (presentation/rigor polishing)
10. **Rework the global “47-test BH correction” strategy.**
   - Identify 2–3 primary hypotheses and correct within-family; keep exploratory analyses clearly labeled.

11. **Report effect sizes with uncertainty at the right level.**
   - For cross-tissue claims, CIs across donors/datasets matter more than edge bootstrap CIs.

12. **Tighten internal consistency in the pseudotime section** (pair counts, denominators, and interpretation of TF activity vs mRNA).

---

## Specific high-risk statements to tone down (or support)
- “Cross-tissue replication confirms signal is real” → should become conditional on stronger benchmarks.
- “Attention encodes co-expression, not regulation” → should be softened to “pooled attention edges are dominated by co-expression under standard extraction; regulatory signal is not recoverable without stratification and may not correspond to causal regulation.”
- Any numerical thresholds (e.g., “200 cells threshold”, “1.85× improvement”) → currently too context-specific; require external validation.

---

## Final recommendation
**Reject at top-tier in current form; potentially accept at a workshop/interpretability venue after major revisions.** The core idea (diagnostic auditing of attention-based GRN claims) is valuable, but the paper must (1) eliminate toy-sized evaluation sets, (2) implement genuinely held-out validation for layer selection, and (3) rule out trivial baselines explaining the reported AUROCs.
