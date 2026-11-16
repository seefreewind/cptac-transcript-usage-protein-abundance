# Phase 1B covariate structure

Correlations are within-gene correlations on the primary complete cases, summarized by the median across successful genes. CNV VIF and condition numbers come from the stored Phase 1A fits. CNV stays in the primary model.

Residual variance is the mean square of the principal proportion after projection on intercept, log2(G), WES purity, and CNV. A value below 1e-4 marks a proportion that is nearly determined by M0.

RSEM correlations use the primary cases that also have numeric RSEM. RSEM is not treated as the same scale as the reconstructed total.

