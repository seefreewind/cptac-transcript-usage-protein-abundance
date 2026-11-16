# Annotation opportunity

Longer proteins and proteins with more Pfam accessions have more chances to be placed in S2. The primary model adjusts for length, domain count, and the other covariates that survived the locked missingness and VIF rules. The measured-isoform-row count was removed because its variance-inflation factor exceeded 10. For these genes it is almost the same number as the count of alternative protein-coding transcripts, and that count stays in the model.

Length standardized difference, S2 minus S1: 0.165. Pfam-count standardized difference: 0.454. Overlapping window retains 2267 S2 genes and 2460 S1 genes.

The length quintiles were cut on the pooled S1 and S2 lengths, without using partial R². Under Literal A the gene-equal differences from the shortest to the longest quintile are 0.000457, 0.005982, 0.005392, 0.004979, and 0.009128. Under Literal B they are 0.001326, 0.004929, 0.005754, 0.004411, and 0.007274. The shortest proteins show a difference near zero. The longest proteins show a larger difference. S2 is also more common in the longest quintile. The linear length adjustment and the overlapping window are the prespecified protections against that pattern. The window estimates are 0.005756 under A and 0.005731 under B. The primary coefficient remains the full-sample adjusted estimate.

The measured-isoform-row count had variance-inflation factors of 7,732 under A and 7,671 under B. It was removed by the locked rule. The alternative protein-coding transcript count, which carries almost the same information, was retained. Its retained VIF is about 1.07.

A SwissProt accession that differs from the TrEMBL accession was not treated as an ambiguous protein sequence. The ambiguous-mapping sensitivity removes genes with a FASTA header conflict or more than one UniProt row.
