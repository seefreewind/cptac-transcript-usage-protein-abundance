# Phase 1A final decision

## A. Execution status

EXECUTION_PASS

## B. Scientific result

CANCER_SPECIFIC_OR_HETEROGENEOUS

The classification uses the frozen 0.02 threshold on the Hartung–Knapp interval, not on a P value. For all ten cancers, the REML pooled mean partial R² is 0.02481 for Literal A (95% CI 0.02222 to 0.02740) and 0.02396 for Literal B (95% CI 0.02162 to 0.02631). Both lower confidence bounds are above 0.02. The paired-support analysis also stays above 0.02, and every leave-one-cancer-out confidence interval stays above 0.02.

The 95% prediction intervals do not. Literal A spans 0.01642 to 0.03320 and Literal B spans 0.01642 to 0.03151. I² is 0.918 for A and 0.910 for B. Cancer-level jackknife intervals for CCRCC, COAD, LUAD and PDAC include values at or below 0.02 under both inverses. BRCA, GBM, HNSCC, LSCC, OV and UCEC have intervals entirely above 0.02 under both inverses. The average increment is therefore not a claim about every cancer.

Gene-level medians are 0.0080 to 0.0121. The share of successful genes with partial R² at least 0.02 is 0.289 to 0.391. The cancer-level mean is the unweighted mean of a right-skewed gene distribution.

## C. Evidence boundary

The 1,027 tumors are matched at case level. The protein outcome is the stored gene-level abundance. Literal A and Literal B are co-primary transformations, and neither was selected from the protein fits. This analysis does not identify an isoform-specific protein abundance and does not demonstrate protein functional loss. Partial R² is incremental statistical explanation. M0 and M1 use RNA-derived variables from the same isoform matrix. Under Literal A, a gene–case whose stored isoform values are all zero remains in the model with equal transcript weights; that behavior is tested by sensitivity S1.

## D. Next-phase authorization

PHASE 1B AUTHORIZATION = FALSE
