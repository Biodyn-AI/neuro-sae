# Codex Adversarial Review R19 (main.tex)

**Overall rating (1–10): 8.7 / 10**  
(Previous: 8.5. This revision is meaningfully clearer and more internally consistent; remaining weaknesses are now mostly about study design/power and external validation, not prose.)

## What improved since last round
1. **“Metric map” reconciliation paragraph works**: the cascade (pooled ≈0.5 → per-layer ≈0.72 → gene-level nulls ≥0.75 → expression-matched ≈0.65) now reads as a *single coherent story* rather than a set of apparently contradictory metrics.
2. **Leave-one-TF-out limitation is correctly acknowledged**: you preempt the obvious robustness critique instead of leaving reviewers to infer it.
3. **Statistical tests table inclusion** is a nontrivial credibility win: the framework-level “47 tests” claim no longer feels hand-wavy.

## Summary judgment
The paper is now close to the ceiling for *text-only* improvements. The central argument is crisp and reads like a diagnostic framework rather than a loose collection of results. The main residual risks are not rhetorical—they are **(i) evaluation set size/power**, **(ii) dependence on curated benchmarks with known biases**, and **(iii) boundary between “diagnostic” and “biological” claims**, especially when numbers look strong (AUROC ~0.72/0.84).

If you changed nothing except text, I doubt you can push acceptance odds much higher than this, because the remaining objections will be “do more experiments / stronger validation.” That said, there are still a few *surgical* text edits that could reduce reviewer attack surface.

## Major remaining vulnerabilities (adversarial)

### 1) Small positive edge set in cross-tissue scGPT-12L analysis is still the weak link
You do flag this, but a skeptical reviewer can still argue:
- With **27–28 TRRUST edges across 14 TFs**, AUROC estimates are unstable and TF-level dependence violates i.i.d.
- The dramatic invariance across tissues/layers could be a **quantization/rounding artifact** (values mostly 0.71–0.72 / 0.84–0.85) rather than a deep “degeneracy” phenomenon.

**Text-only mitigation:** explicitly state (in the cross-tissue section) that the near-constant per-layer AUROC values could be driven by (a) small edge set, (b) averaging/rounding to two decimals, (c) the particular negative universe, and that the *interpretation* rests primarily on the null-model comparison + expression-matched evaluation, not on the flatness alone.

### 2) Risk of over-strong language (“degenerate evaluation design”)
The critique is valid, but “degenerate” can read as rhetorically loaded and may antagonize.

**Text-only mitigation:** replace with something like **“benchmark-confounded evaluation regime”** or **“evaluation dominated by gene-level confounds”** and reserve “degenerate” for a sentence that is explicitly framed as an interpretation, not a fact.

### 3) Multiple-testing correction messaging: could invite pedantic objections
You state “47 tests” and that some key findings survive BH. But throughout the paper you still report many **raw p-values** (and some q’s). A hostile reviewer may nitpick consistency:
- In abstract you cite `p_raw < 1e-4` (fine), but the paper also claims framework-level BH is applied.
- The reader may wonder: *which* results are raw vs corrected in each section?

**Text-only mitigation:** Add a one-sentence reminder in Results/Discussion that **raw p-values are used in descriptive correlational statements, while inferential claims are referenced with q** (or explicitly label q where it matters). Right now it’s *mostly* clear but not uniformly.

### 4) Pooled vs per-layer comparisons still mix models/datasets/edge universes
You acknowledge differences, but the narrative sometimes reads like a single system:
- Pooled AUROC ~0.5: cross-model scGPT+Geneformer table at 200/500/1000 (immune?)
- Per-layer AUROC ~0.72: scGPT-12L cross-tissue with 27 edges
- scGPT-18L per-layer ~0.70: brain tissue 8,330 TRRUST edges

A reviewer can argue you are comparing incomparable evaluation regimes.

**Text-only mitigation:** In the “Reconciling” paragraph, add a compact bullet list stating explicitly that these four AUROC values differ in **(model, tissue, gene vocabulary size, positive edge count, negative universe)** and therefore are not directly numerically comparable—only conceptually.

### 5) Expression-matched evaluation: reviewers may ask whether matching is “too loose”
Matching within ±20% on mean expression and detection rate is reasonable, but reviewers could ask:
- sensitivity to tolerance (±10%, ±5%)
- how many matched negatives actually found per positive (sometimes <<50?)

**Text-only mitigation:** even without new analysis, add: “we report the median number of matched negatives per positive in Appendix/footnote” *if that number is already available elsewhere*. If not available, add as a stated limitation.

## Minor issues / nitpicks worth fixing (text-only)
1. **TRRUST edge count phrasing inconsistency:** “27–28 mapped edges” vs “27-edge evaluation set” appears multiple times. Pick one and explain why it varies (mapping differences across tissues?).
2. **Grammar/precision:** “pseudotime validation is underpowered” vs earlier “fails for 79%” — you do nuance it, but still could be read as contradictory. Consider rephrasing “fails” → “does not recover expected directionality for most pairs under this operationalization.”
3. **CSSI on real data:** you call the findings exploratory but also present many decimals and variants; could add one line that CSSI is mainly a *layer-finding tool* unless validated out-of-sample (you already do, but reinforce).
4. **Baseline comparison conclusion:** “issue is tissue-specific rather than attention-specific” is plausible, but could be softened: benchmarks may be mismatched / context incomplete; avoid sounding like a definitive diagnosis.

## Is this the ceiling for text-only improvements?
**Almost.** The manuscript is now well-structured and defensible in writing. Further large gains likely require **new evidence**, not rewording:
- larger/less biased regulatory ground truth or perturbation-scale validation
- stronger positive controls
- sensitivity analyses for expression-matching and small-edge-set instability

You can still squeeze a small boost (maybe +0.1–0.2) by implementing the mitigation edits above, but the next substantive jump in perceived quality needs experiments.

## Bottom line
- **Current state:** strong, coherent diagnostic paper with unusually good null-model framing; very close to a “final” narrative.
- **Remaining reviewer pressure points:** small edge sets, benchmark bias, comparability across regimes, and wording that could be read as overconfident.
