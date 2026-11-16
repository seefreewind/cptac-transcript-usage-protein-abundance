# Structural difference feasibility

Four evidence layers stay separate.

A. Transcript structure. 10532 of 10795 candidate genes have a matrix alternative with an exon-composition difference, a CDS-length difference, an alternative CDS boundary, or an NMD biotype.

B. Predicted protein sequence. 9608 candidate genes have at least one protein-coding matrix alternative whose GENCODE v34 translation string differs from the principal translation. Genes with a protein-coding alternative in the matrix: 10789.

C. Predicted functional-region annotation. 4829 candidate genes have a Level 2 pair: the sequences differ and a supported Pfam, signal-peptide, or transmembrane difference is present.

D. Experimentally supported protein difference. The count is zero. No isoform-specific peptide sequence is in the frozen inputs.

The sequence method is exact string equality, `phase2a_exact_v1`. Unequal strings have a blank alignment identity and a length-difference lower bound. That lower bound is not a residue alignment.

NMD translations are kept at Level 0. A predicted ORF change is not counted as a stable protein product.

The comparison uses transcripts present in the isoform matrix. For the 10,795 candidate genes, GENCODE v34 also annotates 10947 NMD transcripts in 4676 genes and 17466 protein-coding transcripts that are absent from the matrix. Those unmeasured transcripts were not compared.
