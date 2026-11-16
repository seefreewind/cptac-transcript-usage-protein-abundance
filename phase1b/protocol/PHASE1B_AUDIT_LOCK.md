# Phase 1B audit lock

Locked on 2026-11-14 before the Phase 1B calculations. This file fixes the audit rules. It does not replace the Phase 1A model, the 0.02 threshold, or the co-primary status of Literal A and Literal B.

Phase 1A files are read only. A reproducibility failure stops the audit and sets TECHNICAL_HOLD. Protein models are not refit.

## Reproducibility

Cancer-level means, medians, standard deviations, and quantiles already stored in Phase 1A are recomputed from `primary_literal_A_gene_fits.tsv.gz` and `primary_literal_B_gene_fits.tsv.gz`. Integer fit counts must match exactly. Floating comparisons use absolute tolerance 1e-7 for means, medians, and standard deviations, and 1e-6 for quantiles. The gene-level partial R² values were stored at 10 significant figures, so this tolerance accepts that rounding.

Jackknife standard errors are taken from the Phase 1A cancer summary because the gene file does not store case membership. Each stored interval must equal mean ± t_{n-1} × jackknife SE within 1e-8. REML, Hartung–Knapp intervals, Cochran Q, I², and prediction intervals are recomputed from those means and standard errors. The REML objective and the t_{k-2} prediction interval follow the Phase 1A estimator. A second one-dimensional grid over τ² must agree with the optimizer within 1e-6.

A mismatch beyond these tolerances is TECHNICAL_HOLD. The audit then stops.

## What I² is allowed to mean

I² is reported as a description of between-cohort dispersion relative to the jackknife standard errors. It is not treated as evidence of biological heterogeneity by itself. The audit separates that dispersion into within-gene uncertainty, between-gene dispersion, jackknife uncertainty of the cancer mean, and between-cancer τ².

High I² is called robust only when the recomputed I² stays at or above 0.75 for both inverses in all three of: the ten-cancer analysis, the analysis excluding OV, and the paired-support analysis, and each of those prediction intervals still includes values below 0.02.

## Gene-distribution summaries

The unweighted arithmetic mean remains the primary estimand. The following are descriptive:

- percentiles 10, 25, 50, 75, 90, 95, 99
- fractions at 0.01, 0.02, 0.05, 0.10
- fraction with q < 0.05, and with both q < 0.05 and partial R² ≥ 0.02
- top 1%, 5%, and 10% shares of the sum of nonnegative partial R², using k = max(1, ceil(p × n))
- Gini coefficient of max(partial R², 0)
- effective gene count 1 / Σ s_i², where s_i is the share of that nonnegative sum
- 10% and 20% trimmed means
- 10% winsorized mean, labeled as a sensitivity summary only

No gene identifier is written into a ranked table.

Distribution class, applied separately to each cancer and inverse:

- BROAD_WEAK when median / mean > 0.70 and the top 10% share is below 0.35
- CONCENTRATED when the top 1% share exceeds 0.25 or the top 10% share exceeds 0.60
- MIXED_DISTRIBUTION otherwise

## Sample size and power

Within-cancer N tertiles are computed on successful fits. A tertile contrast is interpreted only when the within-cancer 90th minus 10th percentile of N exceeds 2. Otherwise the N distribution is recorded as concentrated and the tertile gap is not used as evidence.

The planning detectable partial R² uses the frozen anchor λ = 7.848880, f² = λ / N, partial R² = f² / (1 + f²). Exact power at partial R² of 0.02 and 0.10 uses a noncentral F with 1 and N − 5 degrees of freedom and noncentrality f² × N, where f² = partial R² / (1 − partial R²). N = 71 is not treated as a power guarantee for a partial R² of 0.02.

The independence standard error sd / sqrt(n_genes) is reported only as an anti-conservative approximation. The jackknife standard error remains the cancer-level uncertainty.

## Measurement definitions

Close to zero means principal proportion < 0.05. Close to one means principal proportion > 0.95. Entropy is the base-2 Shannon entropy of the isoform composition, with 0 log 0 taken as 0. It is undefined when the linear total is 0.

Primary usage summaries give each successful gene equal weight. A case-pooled mean is stored as a labeled extra column.

Residual variance of the principal proportion is the residual mean square after linear projection on an intercept, log2(G), WES purity, and CNV. A gene is called nearly determined by M0 when that residual variance is below 1e-4.

Correlations are within-gene Pearson and Spearman correlations, then summarized by the median across successful genes. Stored CNV VIF and condition numbers are summarized from the Phase 1A fits. CNV is not removed from the primary model.

Protein values are never imputed and never replaced by zero. Missing protein is not read as functional loss. Phosphoproteomics is not entered into any calculation.

## Meta-analysis

The primary cross-cancer estimator remains REML with the Hartung–Knapp interval. Fixed-effect inverse-variance summaries, the median of the ten cancer means, and a 20% trimmed mean of those ten means are descriptive. Leave-one-cancer-out estimates include prediction intervals. No cancer is dropped because it weakens the pooled mean.

## Meta-regression qualification

There are 10 independent cohorts. That number is below the lock’s minimum of 15 independent cohorts for a qualified meta-regression. The audit therefore records META_REGRESSION_NOT_QUALIFIED and does not fit a moderator model.

The fraction of genes with partial R² ≥ 0.02 is not an eligible moderator, because it is a functional of the same distribution as the cancer-level mean. Eligible descriptive correlates, each exploratory and each with n = 10, are median N, the fraction of candidate genes with N ≥ 71, median principal proportion, median residual variance of the principal proportion, protein-cell missingness, and median CNV VIF. Pairwise |Spearman| > 0.70 is reported as confounding. Benjamini–Hochberg adjustment is applied within the transformation-specific family of those eligible correlates. None of them becomes the primary explanation.

## Explanation class

Applied only after the calculations, using the same rule for the evidence that is actually in hand.

1. STRONG_EVIDENCE_OF_MEASUREMENT_OR_POWER_CONTRIBUTION when, for both inverses, either the cross-cancer range of high-N-tertile means is at most half the unrestricted range, or |Spearman| between the cancer mean and median N is at least 0.80 with a same-signed within-cancer tertile gap in at least 8 cancers and a median absolute tertile gap of at least 0.005. Tertile evidence is ignored where the N distribution is concentrated.

2. EVIDENCE_OF_CANCER_SPECIFIC_EFFECT_HETEROGENEITY when, for both inverses, the high-N-tertile range remains at least 80% of the unrestricted range where tertiles are defined, every eligible |Spearman| is below 0.30, the Gini range across cancers is below 0.10, and the unrestricted range exceeds 4 times the median jackknife standard error.

3. MIXED_MEASUREMENT_AND_CANCER_SPECIFIC_HETEROGENEITY when, for both inverses, the high-N-tertile range remains at least half the unrestricted range and at least one eligible |Spearman| is at least 0.70.

4. INSUFFICIENT_EVIDENCE_TO_EXPLAIN_HETEROGENEITY in every remaining case, including the case where meta-regression is unqualified and the patterns above do not fire.

## Phase 2 readiness

This label does not authorize Phase 2.

- STOP_CURRENT_DIRECTION if reproducibility fails, or if the Literal A all-zero sensitivity moves any cancer mean by more than 0.01 and across 0.02.
- PHASE2_REQUIRES_REDESIGN if the top 1% share exceeds 0.50 in at least 6 cancer–inverse results, or if Literal A and Literal B disagree on which cancers have a jackknife interval entirely above 0.02.
- QUALIFIED_FOR_PHASE2 only if reproducibility passes, both prediction-interval lower bounds exceed 0.02, and the explanation class is STRONG_EVIDENCE_OF_MEASUREMENT_OR_POWER_CONTRIBUTION. A qualified pan-cancer claim is what that label would have required. The Phase 1A prediction intervals already sit on both sides of 0.02, so this label is expected to remain closed unless the reproduction itself differs.
- QUALIFIED_FOR_NARROW_PHASE2 if reproducibility passes, both pooled Hartung–Knapp lower bounds exceed 0.02, A and B agree on the cancer intervals, the all-zero sensitivity does not move the means by more than 0.002, and the redesign conditions do not fire.

Narrow means a later, separately authorized analysis would keep the same endpoint inside each cancer. It does not authorize gene ranking, pathway analysis, isoform switching, or a mechanism claim.

## Prohibitions

No DGE, DTU, isoform-switch discovery, gene ranking, enrichment, survival analysis, machine learning, phosphoproteomic rescue, threshold change, inverse selection, or pooled patient-level regression. Phase 2 is not started.

PHASE 2 AUTHORIZATION remains FALSE at the end of this audit.
