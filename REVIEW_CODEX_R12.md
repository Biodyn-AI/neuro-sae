# Adversarial Review (Codex R12) — Round 2

**Manuscript:** *Attention in Single-Cell Foundation Models Reflects Co-Expression, Not Regulation: A Diagnostic Framework with Null-Model Benchmarking*

## Overall score (1–10)
**6.5/10** (up from **5.5/10** in R11)

The manuscript improved materially: the title/abstract are less hype-y, the new “dumb baseline” table is a strong self-skeptical control, and the narrative now more clearly acknowledges that the apparent AUROC “signal” is likely dominated by gene-level prominence/expression statistics.

However, the paper is still **not at top-venue (Nature Methods / Genome Biology) rigor** because several “core” conclusions still lean on underpowered / biased benchmarks, and some remaining text *still* implies “recovered regulation” despite the null-model results directly undermining that interpretation.

---

## What improved since R11 (real progress)

1. **Null-model benchmarking (new Table \ref{tab:dumb_baselines}) is a major upgrade.**
   - Showing that **variance/mean-expression/detection-rate products beat attention AUROC** is exactly the kind of control that prevents over-interpretation.
   - This also resolves the prior “tension” (pooled attention correlates ~0 with TRRUST but per-layer AUROC ~0.72): the AUROC can be achieved without pairwise structure.

2. **Claims are toned down in several places** (more caveats, explicit “exploratory” notes).

3. **Expanded edge sets** (OmniPath mention, full DoRothEA “all” at 483 edges) is directionally correct.

Net: the paper now reads more like a **diagnostic audit + benchmark critique**, which is its strongest lane.

---

## Major remaining issues (still blocking a strong acceptance)

### 1) The paper still contains an internal contradiction: it *both* concedes benchmark confounding *and* claims “genuine regulatory recovery”
Even after adding null models, the Discussion states:
> “…proper layer and cell-state stratification recovers genuine regulatory signal.”

But Finding 3 + Table \ref{tab:dumb_baselines} indicates the opposite: **the “regulatory signal” in TRRUST/DoRothEA AUROC is explainable by gene-level statistics**.

**Action:** You need to pick a coherent central claim.
- If the null-model result is true (and it looks plausible), then the right conclusion is:
  - *Attention AUROC on these benchmarks is not evidence of regulation; it’s evidence of benchmark enrichment for high-expression/high-variance genes.*
  - CSSI is a **diagnostic** for where *some* benchmark-discriminative signal resides, not proof of regulatory structure.

Right now, the “Core Contribution 2” section still reads like a recovery story (AUROC 0.69–0.71) rather than “this AUROC may be non-mechanistic.”

### 2) Cross-tissue “replication” is still built on a toy TRRUST set (27–28 edges)
You now acknowledge this limitation, but the abstract still foregrounds per-layer AUROC ~0.72 TRRUST and permutation p<1e-4.

With **27–28 positives**, AUROC is extremely sensitive to:
- which TFs are present,
- whether those TFs/targets are “constitutive / prominent genes,”
- and the negative sampling universe.

The fact that the table shows **every layer ~0.72 and every tissue ~0.72** looks less like biology and more like a degenerate evaluation design (or a marginal-statistic proxy that’s invariant to layer).

**Action:** Either:
- (A) **Reframe** this as a *benchmark pathology demonstration* (recommended), OR
- (B) Replace it with a benchmark with hundreds–thousands of positives per tissue, or per-TF evaluation with adequate targets.

### 3) “Dumb baseline beats attention” is a strong result, but you need to tighten the experimental logic around it
Table \ref{tab:dumb_baselines} is persuasive, yet several critical clarifications are missing for a reviewer:
- **What is the negative set?** All non-edges among HVGs? TF×target universe? Randomly sampled?
- Are negatives **expression-matched** or unconstrained? (Unconstrained negatives essentially guarantee that expression-driven methods win.)
- How much of attention’s AUROC is explained by **marginal gene statistics** (e.g., TF frequency + target frequency) vs any pairwise term?

**Action (high impact):** Add an “expression-matched evaluation”:
- Sample negatives matched by (mean, variance, detection rate) of TF and target.
- Or compute AUROC within strata of gene-expression quantiles.

If attention still performs above matched baselines under matched evaluation, that would be the first real evidence of learned pairwise structure.

### 4) Several supporting analyses are still not directly about attention, yet are used to generalize about “mechanistic interpretability of foundation models”
Batch leakage, ortholog transfer, pseudotime, etc. are explicitly correlation-based. That’s fine as *boundary conditions*, but the paper frequently slides from “correlation-based edges behave like X” to broad implications about attention-based mechanistic claims.

**Action:** Keep them, but tighten the rhetoric:
- clearly separate “attention-specific findings” from “general edge-score caveats.”
- Avoid implying these audits validate/falsify attention mechanisms.

### 5) The global “47-test BH correction” is still a methodological smell
Applying BH across a heterogeneous grab-bag of tests isn’t wrong mathematically, but it reads like “statistics theater” rather than good experimental design.

**Action:** Predefine 2–4 primary endpoints, with family-wise correction within those families; treat the rest as exploratory.

### 6) Some internal inconsistencies persist
Examples:
- **Pseudotime counts** still don’t reconcile cleanly (“30 hematopoietic pairs” vs “19 of 75 pairs analyzed”).
- scGPT-12L description says “27–28 TRRUST edges” for 1,200 HVGs: that mapping rate is so low that it needs explanation (TF list? directed mapping? filtering?).

These things erode trust.

---

## What I think the paper *is* now (and how to position it)
The best and most defensible “main story” after the new baselines is:

1. **Pooled attention edges strongly reflect co-expression / co-occurrence** (robust across models).
2. Apparent “GRN AUROC” can be misleading; **standard curated benchmarks are confounded by gene prominence**.
3. CSSI/layer stratification is a useful *diagnostic*, but **AUROC improvements do not imply regulation unless gene-level confounds are controlled**.

If you position the paper as **a cautionary methods/benchmarking paper** (“why attention-based GRNs are easy to fool; how to diagnose confounds”), it becomes much more publishable.

---

## Concrete, high-impact revisions to target next

### Tier 1 (would most increase credibility)
1. **Add expression-matched benchmarking** (or stratified AUROC by expression quantiles). This is the single most important next experiment.
2. **Define the evaluation universe and negative sampling** precisely and consistently.
3. **Rewrite abstract/discussion to remove any implication that AUROC ≈ regulation** unless supported under confound-controlled evaluation.
4. **Strengthen held-out validation** for layer selection (donor split or dataset split).

### Tier 2
5. Quantify variance explained:
   - regress attention edge score on marginal stats + co-expression; report partial correlations.
6. Report per-TF performance (AUROC per TF with adequate targets), not just global AUROC.

---

## Bottom line
Yes—the paper improved substantially and is more intellectually honest. The new null-model comparison is a real step toward a solid mechanistic/benchmark critique.

But to clear a high bar, you must **stop treating benchmark AUROC as evidence of regulatory signal** and implement **confound-controlled evaluation** (expression-matched negatives or stratified analysis). Without that, the central “recovery” narrative remains unsupported by the paper’s own strongest new result.
