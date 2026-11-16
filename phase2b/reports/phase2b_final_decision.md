# Phase 2B final decision

## A. Structural category reproducibility

PASS

The Phase 2A counts reproduced: 10,795 candidate genes, 39,490 measured protein-coding transcripts, 10,742 genes with a principal SwissProt or TrEMBL accession, 9,608 genes with a different predicted sequence, 4,829 genes with a supported structural difference, and 10,228 genes with a Pfam accession on the principal protein.

## B. Common annotation-support adequacy

ADEQUATE

## C. Primary S2-versus-S1 association

Literal A: adjusted S2 minus S1 = 0.004644, cluster-robust interval 0.003269 to 0.006019, genes 6887, rows 50164. Adjustment set: principal_length, n_pfam, pfam_coverage, n_alt_protein_coding, n_complete, protein_fraction, median_log2_g.
Literal B: adjusted S2 minus S1 = 0.004697, cluster-robust interval 0.003346 to 0.006048, genes 6859, rows 49690. Adjustment set: principal_length, n_pfam, pfam_coverage, n_alt_protein_coding, n_complete, protein_fraction, median_log2_g.

## D. Cancer-specific interaction

Every cancer-specific adjusted estimate is positive. Under both inverses the smallest is COAD, 0.001991 and 0.001803, and those intervals include zero. The largest is LSCC, 0.008059 and 0.007571, and those intervals lie above 0.005. GBM, HNSCC, and UCEC also have estimates above 0.005 with intervals above zero. The locked cancer-specific label requires, in both inverses, one cancer at or below zero whose interval stays below 0.005. No cancer meets that low side. The magnitude varies. The sign does not.

## E. Literal A/B consistency

The two adjusted estimates have the same sign: True.

## F. Robustness to annotation opportunity

Length standardized difference 0.165. Pfam-count standardized difference 0.454. Window S2 2267, S1 2460. The window estimates remain positive. The shortest length quintile is near zero and the longest is larger, so a single adjusted mean does not describe every protein length. That gradient was not used to change the category rule or the 0.005 magnitude.

## G. Independent evidence limitations

Level 3 evidence remains zero. The association uses gene-level protein abundance and predicted annotation.

## H. Scientific novelty

PARTIAL

## Overall

NO_PRACTICAL_STRUCTURAL_ASSOCIATION

Genes with predicted principal-versus-alternative protein structural differences were compared with genes whose predicted sequences differ without a supported structural annotation. The comparison is an association with transcript-usage incremental explanatory values for gene-level protein abundance. The locked practical magnitude is 0.005. The adjusted estimates are 0.004644 under Literal A and 0.004697 under Literal B.

The analysis does not establish that an alternative isoform is expressed as a protein, that a domain is experimentally absent, or that a structural difference causes the partial R².

PHASE 2B EXECUTION = COMPLETED

PHASE 2C AUTHORIZATION = FALSE
