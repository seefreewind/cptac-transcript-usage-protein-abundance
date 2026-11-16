# Gene-reconstruction consistency

Audit date: 2026-11-14. For each gene and each shared RNA case, G is the sum of the inverted isoform rows. The comparison is between log2(G) and the stored RSEM value in `*_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt`. The RSEM matrix is an upper-quartile normalization of coding genes scaled to 1,500, then log2. The isoform matrix is an FPKM quantification. Equality of the two numbers is not assumed. A higher correlation is not a reason to prefer either inverse.

Every isoform-matrix gene is present in the RSEM matrix. The pair counts below equal the gene-sample counts in the usage report. Full stratum counts are in `phase0_6/manifests/cancer_specific_quantification_audit.csv`.

## Overall agreement and the scale offset

Residual means log2(G) minus the stored RSEM value.

| Cancer | Pairs A | Spearman A | Pearson A | Median residual A | Residual p05 A | Residual p95 A | Pairs B | Spearman B | Pearson B | Median residual B | Residual p05 B | Residual p95 B |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCA | 7,340,949 | 0.916 | 0.902 | 0 | −6.625 | 1.000 | 3,377,147 | 0.925 | 0.922 | −6.001 | −7.837 | −3.404 |
| CCRCC | 6,248,907 | 0.926 | 0.877 | −1.495 | −6.932 | 0 | 3,832,917 | 0.899 | 0.887 | −5.996 | −7.973 | −2.998 |
| COAD | 6,430,914 | 0.913 | 0.900 | 0 | −6.630 | 1.000 | 2,905,000 | 0.927 | 0.925 | −6.102 | −7.963 | −3.531 |
| GBM | 6,006,231 | 0.930 | 0.880 | −1.769 | −6.962 | 0 | 3,894,818 | 0.899 | 0.887 | −5.936 | −7.926 | −2.910 |
| HNSCC | 6,552,252 | 0.928 | 0.876 | −1.530 | −6.972 | 0 | 4,085,818 | 0.899 | 0.890 | −6.064 | −8.027 | −3.056 |
| LSCC | 6,552,252 | 0.931 | 0.877 | −1.598 | −6.944 | 0 | 4,118,878 | 0.900 | 0.889 | −6.002 | −7.983 | −2.978 |
| LUAD | 6,673,590 | 0.930 | 0.878 | −1.592 | −6.952 | 0 | 4,172,134 | 0.900 | 0.888 | −6.007 | −7.982 | −3.007 |
| OV | 4,974,858 | 0.916 | 0.902 | 0 | −6.602 | 1.000 | 2,185,158 | 0.928 | 0.924 | −6.055 | −7.860 | −3.479 |
| PDAC | 6,370,245 | 0.931 | 0.877 | −1.607 | −7.017 | 0 | 3,979,422 | 0.901 | 0.890 | −6.071 | −8.035 | −3.094 |
| UCEC | 5,763,555 | 0.929 | 0.881 | −1.508 | −6.866 | 0 | 3,554,876 | 0.898 | 0.886 | −5.916 | −7.887 | −2.897 |

Spearman correlation of log2(G) with RSEM is 0.913–0.931 under A and 0.898–0.928 under B. Pearson correlation is 0.876–0.902 under A and 0.886–0.925 under B. The rank agreement says the reconstructed gene sum and the RSEM column order genes and samples similarly. The residual location says they are on different scales. Under B the median residual is about −6 in every cancer: log2(2^Y − 1) sits well below a log2 RSEM value that was scaled to an upper quartile of 1,500. Under A the median residual is 0 in BRCA, COAD, and OV, and about −1.5 to −1.8 in the other seven. Those three cancers are the ones in which more than half of gene-samples have RSEM equal to 0, and a single-row gene with stored isoform 0 has log2(G_A) = log2(1) = 0.

## The shared zero

When every isoform row of a gene-sample is stored 0, the RSEM value is also 0 in 0.9942–0.9991 of those gene-samples. The lowest share is GBM (0.9942) and the highest is OV (0.9991). The two matrices share a detection floor. That agreement is a property of the zero code. It does not identify which inverse produced the zero.

## Expression level

RSEM quartiles are cut on the pairs that enter each inverse. Where the stored RSEM zero mass is larger than 25 percent, the lowest A quartile is that zero mass, RSEM has no variance inside it, and Spearman correlation is undefined. Where the zero mass is also larger than 50 percent (BRCA, COAD, OV), the second A quartile is empty. In the remaining cancers the second A quartile, which is the band just above zero, has Spearman correlation 0.277–0.330. The upper A quartile has Spearman correlation 0.751–0.831.

Literal B drops the all-zero gene-samples, so its four quartiles are populated. Spearman correlation rises from 0.327–0.435 in the lowest quartile to 0.660–0.686 in the highest.

A second cut uses the largest stored isoform value in the gene-sample. At or below 1, Spearman correlation is 0.788–0.853 under A and 0.708–0.738 under B. Above 1, it is 0.797–0.828 under A and 0.797–0.838 under B. The two inverses are closest where the stored values are high, because 2^Y and 2^Y − 1 differ by 1 and that constant is small relative to a large linear value. They differ most where zeros and small positive values share a gene.

## Transcript-row count

Spearman correlation under A declines as the number of isoform rows increases: 0.968–0.994 for one row, 0.962–0.979 for two, 0.917–0.930 for three to five, 0.897–0.914 for six to ten, and 0.828–0.870 for eleven or more. For one row, log2(G_A) equals the stored isoform value, including the shared zero, so this stratum mostly restates the agreement between one isoform row and the gene-level RSEM column.

Under B the one-row stratum excludes stored zeros, and its Spearman correlation is 0.817–0.852. In the eleven-or-more stratum, B’s Spearman correlation is 0.926–0.948. That stratum contains only gene-samples with at least one positive isoform row. The A figure for the same row-count bin still includes all-zero gene-samples, whose log2(G_A) is the log2 of the row count.

## What the audit supports

The isoform rows can be summed, after either literal inverse, into a gene-level quantity that tracks RSEM by rank across all ten cancers. The tracking is weaker just above the zero floor and stronger at high RSEM. The numerical level of log2(G) is not the RSEM level. The primary gene covariate in the frozen model remains log2(G) from the isoform sum. The RSEM column remains the sensitivity covariate. Both statements are applied twice, once under each literal inverse.
