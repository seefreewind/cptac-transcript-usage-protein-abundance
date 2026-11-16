# Phase 1A lock

Locked before any protein regression. Phase 1A authorization is limited to the frozen gene-level models. Phase 1B authorization is false.

The authoritative records, in order, are the Phase 0.6 handover, the statistical-model freeze, the transformation-resolution report, the data-readiness table, the usage and reconstruction audits, the Phase 0.5 analysis plan, the cohort-integrity report, and the two SHA256 ledgers. Phase 0, Phase 0.5, and Phase 0.6 files are not edited.

## Co-primary transformations

Literal A is X = 2^Y. A stored zero contributes 1. An all-zero gene–case stays in the Literal A model, with G equal to the isoform-row count and P equal to 1 over that count, and is flagged `stored_all_zero`.

Literal B is X = expm1(Y ln 2), which is 2^Y − 1. A stored zero contributes 0. A gene–case with G = 0 is excluded. No pseudocount is added.

Neither transformation is chosen from a protein association. Candidate C is a sensitivity only.

## Model

Within each cancer, for each eligible gene and each co-primary transformation, M0 and M1 use the same cases.

M0: protein ~ log2(G) + WES purity + CNV.

M1: protein ~ log2(G) + P + WES purity + CNV.

The outcome stays on the stored protein scale. CNV stays on the stored log2-ratio scale. There is an intercept. Cancers are not pooled. No batch or plex term is added.

A gene is primary-eligible only when it has at least two isoform rows, the frozen principal transcript is present, protein, log2(G), P, WES purity, and CNV are numeric, N is at least 71, both fits are full rank, the protein outcome has nonzero variance, and P has nonzero residual variation after the M0 predictors. Rank-deficient fits are labeled `RANK_DEFICIENT` and are not given partial R² = 0. A predictor P with no residual variation is labeled `NO_RESIDUAL_P_VARIATION`.

The principal transcript is the unique GENCODE v34 MANE Select transcript, otherwise the unique APPRIS principal-1 transcript. Matrix row order is not a selection rule. GENCODE v50 is not used.

## Numerical choices fixed here

Continuous predictors are centered and scaled to sample standard deviation with divisor n (ddof = 0) inside each fit. The intercept is not scaled. The protein outcome is not scaled. This does not change fitted values or partial R². The reported coefficient of P is the coefficient after that scaling.

Partial R² = (RSS_M0 − RSS_M1) / RSS_M0. The signed partial correlation is the Pearson correlation of the M0 residuals of protein and of P. Their squares must agree within 1e-8. A negative RSS gap smaller in magnitude than 1e-8 times max(RSS_M0, 1) may be clipped to zero and is logged. A larger negative gap stops the affected analysis.

Design rank uses an SVD tolerance of max(n, p) × machine epsilon × the largest singular value. The CNV variance-inflation factor is computed in M1. The Pearson residual for the outlier sensitivity is the M0 residual divided by sqrt(RSS_M0 / residual degrees of freedom), without a leverage adjustment.

The cancer-level mean is the unweighted mean of partial R² over successful fits. Its interval is the delete-one-case jackknife specified in the Phase 1A directive, because the model freeze calls for a case-resampled interval and does not define the resampling design. The point estimate is the full-data mean. The raw interval is compared with 0.02. A copy clipped to [0, 1] is display-only.

The cross-cancer summary uses REML for τ² and a Hartung–Knapp interval. The model freeze requires a random-effects summary weighted by the within-cancer sampling variance and does not name the estimator. The jackknife standard error is that sampling variance. REML and Hartung–Knapp are the Phase 1A operational choices.

## Sensitivities fixed here

S1 drops Literal A all-zero gene–cases and refits. S2 replaces log2(G) with stored RSEM on the primary cases whenever RSEM is numeric. S3 is Candidate C. S4 is not run: the frozen plan gives the CLR pseudocount as half the smallest positive linear value in the gene, and it does not define the one-degree-of-freedom contrast that would replace P. S5 replaces WES purity with ESTIMATE_ESTIMATEScore on the primary cases. S6 refits without CNV only when the M1 CNV VIF exceeds 10, on the same cases. S7 drops every M0 case with absolute Pearson residual above 4, once, and refits when N remains at least 71. S8 repeats the cross-cancer summary without OV and keeps the OV-specific analysis.

## Decision rule fixed here

The practical threshold is partial R² = 0.02, applied to the Hartung–Knapp interval of the cancer-level means.

- `TRANSFORMATION_DEPENDENT` when one transformation has a pooled lower bound above 0.02 and the other does not.
- `ROBUST_PRACTICAL_INCREMENT` when both do, the paired-support analysis also does, no leave-one-cancer-out analysis pulls a lower bound to 0.02 or below, and neither prediction interval extends below 0.02.
- `CANCER_SPECIFIC_OR_HETEROGENEOUS` when both pooled lower bounds are above 0.02 and at least one of those consistency conditions fails.
- `SMALL_OR_UNCERTAIN_INCREMENT` when both pooled intervals include a value at or below 0.02.
- `TECHNICAL_NO_RESULT` when the primary fits are not completed.

P values do not replace the threshold. Genes are not ranked.
