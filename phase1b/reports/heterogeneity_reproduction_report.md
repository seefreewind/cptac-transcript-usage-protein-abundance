# Phase 1B heterogeneity reproduction

Recomputed from the Phase 1A gene-level files at 2026-11-14T11:21:37Z.

Integer counts and stored moments match the Phase 1A cancer summary.
Failed rows: 0.

I² compares between-cancer dispersion with the jackknife standard errors. Those standard errors are about 0.001, while the cancer means span roughly 0.019 to 0.029. A high I² can therefore come from precise cancer means rather than from a large biological gap. The between-gene standard deviation within a cancer is about 0.03 to 0.05, much larger than the jackknife standard error of the mean.

The independence approximation sd / sqrt(n_genes) is smaller than the jackknife standard error because genes share cases. The jackknife standard error is the uncertainty used for the cancer mean.

