# Phase 1B power and precision

At N = 71 the frozen planning formula gives a detectable partial R² of 0.100 for 80% power at two-sided alpha 0.05. The observed cancer means are near 0.02. Exact power for a gene-level partial R² of 0.02 at the observed median N is reported in `results/power_sensitivity.csv` and is low.

The cancer-level mean is an average over thousands of genes, so its jackknife standard error is near 0.001 even though a single gene is not powered for a partial R² of 0.02. Between-gene standard deviations are about 0.03 to 0.05. Between-cancer τ² is on the order of 1e-5.

A cancer with fewer eligible genes is not assigned a weaker biological effect from that count alone. OV remains in every cancer-specific display.

Sample-size rows: 20. Power rows: 20.

