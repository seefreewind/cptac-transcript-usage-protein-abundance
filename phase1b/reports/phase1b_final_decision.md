# Phase 1B final decision

## A. Is the Phase 1A result reproducible?

PASS

Gene-level counts, means, medians, standard deviations, and the stored percentiles were recomputed from the Phase 1A gene files. Jackknife intervals match mean ± t_(n−1) times the stored jackknife standard error. REML, Hartung–Knapp intervals, I², Cochran Q, and prediction intervals were recomputed from the cancer means and those standard errors.

## B. Is high heterogeneity robust?

YES

I² stays high for both inverses in the ten-cancer analysis and in the analysis excluding OV. The prediction intervals continue to include values below 0.02. High I² here is the ratio of a small between-cancer variance to still smaller jackknife standard errors. The absolute gap between the lowest and highest cancer means is about one percentage point of partial R².

## C. Is the heterogeneity explained?

MIXED_MEASUREMENT_AND_CANCER_SPECIFIC_HETEROGENEITY

Meta-regression was not qualified with 10 cohorts, and no moderator model was fit. The locked mixed class fired because the cross-cancer range remained among higher-N genes and because protein-cell missingness was the only eligible correlate with absolute Spearman correlation at least 0.70. That correlation is about −0.72 for Literal A and −0.75 for Literal B. Its Benjamini–Hochberg q values are 0.11 and 0.08. The same moderator is correlated with the fraction of candidate genes at N ≥ 71 (Spearman −0.92) and with the median principal proportion (Spearman about 0.73 to 0.76). Within a cancer, a gene’s protein observation fraction and its partial R² have Spearman correlations of about 0.01 to 0.07. The median gene is fully observed in every cohort. This is an ecological association on 10 cohorts, not a causal assignment.

The five candidate explanations, each left as a description of the audit rather than as a selected cause:

1. Cancer means differ by about 0.01 partial R², which is several jackknife standard errors, and the same four cohorts fall at or below 0.02 under both inverses.
2. Median N differs by cohort, while power at a gene-level partial R² of 0.02 remains low throughout. The cohort with the smallest N is not the cohort with the smallest mean.
3. Principal-proportion location, entropy, and residual variance differ across cohorts; the locked correlation threshold for calling that the source of the mean differences was not met.
4. Cohort-level protein-cell missingness is the correlate that crossed the locked absolute Spearman screen of 0.70. It is confounded with eligibility and with the principal-proportion level. The within-cancer association between coverage and partial R² is about 0.01 to 0.07. Missing protein values were left missing.
5. The gene-level distribution is right-skewed in every cohort. The top 1% share and the Gini coefficient describe a mixed distribution, not a mean carried by a handful of genes and not a flat weak effect.

## D. Is the project ready for Phase 2?

QUALIFIED_FOR_NARROW_PHASE2

This label does not start Phase 2. A narrow qualification would keep the same gene-level endpoint inside each cancer. It does not authorize gene ranking, pathway analysis, isoform switching, or a mechanism claim.

## E. Scientific interpretation

Transcript usage provides a small average incremental explanatory contribution to gene-level protein abundance beyond reconstructed gene-level RNA expression, tumor purity and CNV, but the magnitude varies substantially across cancer cohorts.

The measurements are case-level. The outcome is gene-level protein abundance. Literal A and Literal B remain co-primary. Partial R² is incremental statistical explanation. The audit does not identify isoform-specific protein abundance, protein functional loss, a causal effect of transcript switching, clinical utility, or a therapeutic target.

## Audit status

PHASE 1B AUDIT = COMPLETED

PHASE 2 AUTHORIZATION = FALSE
