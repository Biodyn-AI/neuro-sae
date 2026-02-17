# Codex Adversarial Review (Round 17)

**Paper:** `main.tex` (updated)

## Overall rating: **8.3 / 10**
This is now a cohesive, high-signal diagnostic paper with a clear through-line (co-expression dominance → layer localization via CSSI → null-model benchmarking revealing gene-level confounds, with an expression-matched “escape hatch” that preserves modest pairwise structure). Compared to earlier rounds, the narrative is substantially tighter and the evaluation protocol consolidation helps a lot.

What still prevents **9/10** is not “missing experiments” so much as **remaining credibility/clarity gaps in evaluation design and claims calibration**: a few places still read like strong conclusions are supported by **very small or heterogeneous evaluation sets**, and some analyses mix **attention-derived** and **correlation-derived** edges in a way that can confuse what is being claimed about *foundation-model attention specifically*.

---

## What’s strongest (and now basically in good shape)

1. **Core thesis is coherent and defensible.** The co-expression vs regulation dissociation is clearly stated, repeated consistently, and reconciled with the per-layer AUROC paradox via null models.
2. **CSSI repositioned correctly.** Calling it diagnostic (signal localization / layer selection) rather than “proof of regulation” is the right move.
3. **Null-model benchmarking is the paper’s differentiator.** The “gene-level products match/exceed AUROC” point is persuasive and important.
4. **Methods §2.10 helps a lot.** Centralizing candidate edge universe / negatives / directedness / attention score definition reduces reader uncertainty.

---

## Main blockers to 9/10 (actionable)

### 1) **TRRUST/DoRothEA evaluation set sizes and comparability are still fragile / potentially misleading**

- The scGPT-12L cross-tissue TRRUST panel is **27–28 positive edges** across **14 TFs**, yet the text repeatedly discusses AUROC ≈ 0.72 as if it is a stable, generalizable “uniform” signal. You do acknowledge the limitation, but it still drives major narrative weight.
- This interacts badly with the paper’s own key argument: **gene-level prominence confounding** becomes *easier* to exploit when the positive set is dominated by a handful of famous TFs.

**What to do (minimal):**
- Add a **prominent “small-n edge panel” warning box** right at first mention of Table  \ref{tab:cross_tissue_layers} (not later).
- Report **AUPRC** alongside AUROC (at least for the small-edge settings) because AUROC can look “good” even when the positive rate is tiny or when separation is driven by TF frequency effects.
- Report **per-TF grouped evaluation** (e.g., leave-one-TF-out AUROC, or macro-averaged per-TF AUROC). This would directly test whether performance is simply “did we rank STAT3/TP53/MYC edges highly.”

**Why this blocks 9/10:** A skeptical reviewer will argue: *“Your ‘uniform cross-tissue AUROC’ is an artifact of evaluating ~27 curated edges dominated by famous TFs; the invariance may not mean what you say it means.”* You already basically agree—make it mathematically harder to dismiss.

---

### 2) **Title/abstract headline is still a bit over-absolute relative to your own “expression-matched AUROC 0.646” result**

- Title: “Attention … Reflects Co-Expression, Not Regulation.”
- But later: expression-matched negatives show attention retains AUROC 0.646 vs 0.522 baselines → “direct evidence” of pairwise structure beyond gene-level confounds.

This is not a contradiction, but it invites a reviewer to say you are overselling the negative claim. A more defensible framing is: **“mostly co-expression; regulatory signal is modest and confounded under standard benchmarks.”**

**What to do:**
- Slightly soften the title or add a subtitle clause acknowledging residual pairwise signal (even a single phrase).
- In the abstract’s first sentence-level claim, add a qualifier like **“in standard pooled/benchmark settings”**.

**Why this blocks 9/10:** It’s “reviewer psychology”: absolutist titles draw fire, especially when the paper itself contains counterevidence (even if modest).

---

### 3) **Mixing attention-derived and correlation-derived edges in supporting analyses muddies the scope**

You do add methodological notes (“correlation-based to establish boundary conditions”), which helps. Still, as written, a reader can come away thinking the paper *tests attention* in sections where it doesn’t.

**What to do:**
- Add an explicit, repeated labeling convention: e.g., each supporting subsection header includes **[ATTN]** vs **[CORR]**.
- In Discussion, explicitly list which supporting analyses are **not attention-specific** and state what would need to change to make them attention-specific.

**Why this blocks 9/10:** Clarity of claims is a top-tier acceptance criterion. Right now it’s good, but still not “unassailable.”

---

### 4) **Multiple testing correction claims feel under-documented / potentially inconsistent**

- You state: “All 47 tests … BH FDR … Key findings survive correction: scaling failure (p=0.011), mediation non-additivity (p=0.003), cross-species conservation (p=0.001).”
- But in the main text, many p-values are presented as raw/huge; it’s not always clear which are corrected and which are not.

**What to do:**
- Ensure every reported p-value is annotated as **q** (corrected) vs **p** (raw), consistently.
- Consider a single table (or appendix pointer) listing the 47 tests with raw p and BH q, and referencing it in-line.

**Why this blocks 9/10:** Reviewers often distrust “we corrected globally” unless they can audit it quickly.

---

### 5) **A few quantitative claims still read too “clean” given heterogeneity and underpowered settings**

Examples:
- “critical threshold around 200 cells” is plausible but presented with a precision that’s hard to justify when adjacent CIs overlap.
- Some synthetic results (CSSI ratios, Shapley +91%) are very strong and may read as “generator aligned to theory” (you acknowledge this once, but the paper relies on them heavily).

**What to do:**
- Rephrase “200 cells is a practical threshold” into “**order-of-magnitude guidance**; in our setting ~200.”
- For synthetic sections, add a one-liner specifying which generator knobs are “pro-theory” and which are neutral.

**Why this blocks 9/10:** It’s about epistemic humility matching the demonstrated uncertainty.

---

## Minor issues / polish (won’t block 9/10 but worth fixing)

- **Terminology:** You use “co-expression” and “co-occurrence” as Spearman on raw counts; some readers will object to raw-count correlations. Consider a brief sensitivity note (log1p, Pearson, distance correlation) or justify raw counts explicitly.
- **Candidate universe definition:** Methods say TF status defined by TRRUST/DoRothEA membership; ensure it’s clear whether TF list differs per benchmark and how overlaps are handled.
- **Directedness:** Good that you state directed AUROC; consider reminding readers that many baselines (correlation) are symmetric and how you assign direction.

---

## Bottom line
This is **very close** to a 9/10. The remaining work is mainly:

1) Make the evaluation panels **harder to dismiss** (per-TF grouping, AUPRC, clearer small-n caveat),
2) Align the **headline claims** (title/abstract) with the nuanced conclusion that **residual pairwise signal exists**,
3) Improve auditability of **global multiple testing correction**,
4) Make the boundary between **attention-specific** and **general GRN-inference boundary-condition** analyses impossible to misread.

If those are tightened, the paper becomes not just compelling but *reviewer-resistant*.
