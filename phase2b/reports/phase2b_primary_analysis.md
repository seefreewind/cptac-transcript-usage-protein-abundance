# Phase 2B primary analysis

The outcome is the frozen Phase 1A partial R². S2 is a predicted structural difference. S1 is a predicted sequence difference without a supported Pfam, SignalP, or TMHMM difference. The coefficient is an adjusted difference in means. It is an association.

Literal A: adjusted S2 minus S1 = 0.004644, cluster-robust interval 0.003269 to 0.006019, genes 6887, rows 50164. Adjustment set: principal_length, n_pfam, pfam_coverage, n_alt_protein_coding, n_complete, protein_fraction, median_log2_g.
Literal B: adjusted S2 minus S1 = 0.004697, cluster-robust interval 0.003346 to 0.006048, genes 6859, rows 49690. Adjustment set: principal_length, n_pfam, pfam_coverage, n_alt_protein_coding, n_complete, protein_fraction, median_log2_g.

The practical magnitude locked before fitting is 0.005 partial R². Both adjusted point estimates are below that magnitude. Both cluster-robust intervals lie above zero and extend above 0.005. Common-support status: qualified.

Literal A uses 3,339 S2 genes and 3,548 S1 genes. Literal B uses 3,327 S2 genes and 3,532 S1 genes. These are the genes with at least one successful complete fit. The category universe remains 4,829 S2 genes and 4,779 S1 genes.

The unadjusted cancer-indicator estimates are 0.004899 for Literal A and 0.004798 for Literal B. The gene-equal-weight differences are 0.005430 and 0.005010. The gene-level cluster-bootstrap intervals are 0.003304 to 0.006020 and 0.003413 to 0.005953. The within-cancer stratified bootstrap intervals are 0.003904 to 0.005419 and 0.003957 to 0.005488. One thousand draws used seed 20261114.

On the gene-mean scale, Literal A S2 has mean 0.0289, median 0.0189, and 47.1% of genes at or above 0.02. S1 has mean 0.0235, median 0.0170, and 41.0% at or above 0.02. The arithmetic mean remains the location used for the primary coefficient. The standardized gene-equal difference is 0.178 under A and 0.177 under B.

Every leave-one-cancer-out estimate stays between 0.004179 and 0.004921. Removing 20 genes with a FASTA header conflict or more than one UniProt row changes the estimates to 0.004655 and 0.004709. No gene crossed the one-step influence screen of 2/sqrt(G).
