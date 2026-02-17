# Hostile Review of `main.tex` (Codex, Round 2)

## Overall score: **3/10**

This reads like a paper trying to turn a set of mostly negative, underpowered, and partially circular analyses into a "systematic framework" by sheer volume. The strongest claims are repeatedly qualified by the manuscript's own limitations, and the one purported constructive contribution (CSSI) is validated mainly in settings that are either synthetic or vulnerable to selection/circularity.

## (1) Underpowered results

- **Mediation bias is severely underpowered**: only **16 run-pairs** (`main.tex:263`, `main.tex:281`). You still promote "62.5% non-additivity" as a central finding (`main.tex:699`). That is not robust prevalence estimation; it is pilot-scale.
- **Perturbation validation collapses after correction**: "No perturbation validation tests survived framework-level multiple testing correction" (`main.tex:336`). Yet the paper keeps leveraging these outcomes narratively.
- **Pseudotime result is null after correction**: adjusted p-value reported as non-significant (`main.tex:407`, `main.tex:416`, `main.tex:699`). This is framed as insight, but statistically it is an inconclusive negative.
- **Scaling characterization is coarse**: only four cell-count points for the core scaling claim (`main.tex:736`), insufficient to justify confident statements about curve shape/mechanism.

## (2) Circular validation

- You explicitly acknowledge **reference database circularity** (`main.tex:734`) and then continue to ground key conclusions in TRRUST/DoRothEA recovery throughout.
- The paper admits CSSI real-data validation remains vulnerable to **same-dataset layer selection/reporting circularity** (`main.tex:544`, `main.tex:550`, `main.tex:578`), but still presents layer-specific AUROC gains as strong mechanistic support.
- Synthetic validation is admitted to be **framework-consistent by construction** (`main.tex:607`), i.e., not independent validation.

## (3) Alternative explanations not adequately addressed

You list alternatives but do not resolve them experimentally:

- **Scaling failure vs reference mismatch/overfitting** (`main.tex:231`) is "addressed" by a synthetic saturation analysis (`main.tex:233-237`) that still depends on your proxy assumptions, not independent biological ground truth.
- **Cross-tissue inconsistency vs batch/protocol confounding** is acknowledged (`main.tex:328`), and batch leakage is high (`main.tex:430`), but no decisive deconfounded design is shown.
- **Near-random GRN recovery** could reflect benchmark mismatch in brain context (`main.tex:251`, `main.tex:619`), not necessarily method failure; this remains unresolved.

## (4) Single weakest result

**Condition-specific perturbation validation** (`main.tex:333-340`) is the weakest result.

Why this is weakest:
- It is the only section where you directly admit complete failure after multiplicity control (`main.tex:336`).
- It is central to causal credibility, yet ends as exploratory noise.
- Without perturbation support, many "mechanistic" interpretations reduce to correlation-pattern storytelling.

## (5) One experiment to improve the weakest result

Run a **pre-registered, out-of-domain perturbation validation on Replogle-scale CRISPRi data**, with strict train/calibration/test separation by perturbation target and cell context.

Minimal design:
- Train any mapping/hyperparameters on existing datasets (Dixit/Adamson).
- Freeze the full pipeline (including layer/head/CSSI choices) before touching Replogle.
- Evaluate on held-out genes not seen in tuning.
- Primary endpoint: signed effect concordance + rank correlation, with BH correction for a pre-specified small set of tests.
- Report negative controls (matched random genes, permuted targets) and effect sizes, not just p-values.

If this succeeds, the paper gains causal credibility. If it fails, the manuscript should be reframed as a cautionary benchmark paper, not a constructive mechanistic framework.
