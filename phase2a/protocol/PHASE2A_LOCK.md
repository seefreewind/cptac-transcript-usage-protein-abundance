# Phase 2A audit lock

Locked on 2026-11-14 before the annotation counts were calculated.

This phase asks whether the frozen GENCODE v34 transcript universe can be joined to protein sequences and to functional-region coordinates. It does not refit the Phase 1A protein models, and it does not rank genes.

## Sources

- Transcript coordinates and biotypes: `phase0/raw/annotation/gencode.v34.annotation.gtf.gz`, SHA256 `941926bde63e86f4eebbd29a8f747937288750da6b7dda8b21a920980c38b6e1`.
- Principal rule: the Phase 1A rule. Unique MANE Select, otherwise unique APPRIS principal-1 when the gene is absent from the MANE tag set. A non-unique MANE tag blocks the APPRIS fallback.
- Measured transcripts: the isoform-matrix transcript identifiers. The repaired BRCA matrix is the reference set. The other nine matrices must contain the same identifiers.
- Protein sequences: GENCODE v34 `gencode.v34.pc_translations.fa.gz`, retrieved from the GENCODE release 34 FTP on 2026-11-14. Last-Modified on that server: Wed, 29 Apr 2020 15:45:23 GMT.
- Reviewed accessions: the local GENCODE v34 SwissProt metadata file.
- Unreviewed accessions: GENCODE v34 `gencode.v34.metadata.TrEMBL.gz` from the same FTP and the same Last-Modified time.
- Domain coordinates: Ensembl release 100 `protein_feature`, `translation`, and `analysis` tables. GENCODE v34 is the Ensembl 100 gene build, so this is the matching coordinate release. Last-Modified observed for `protein_feature.txt.gz`: Wed, 20 May 2020 10:53:36 GMT.
- UniProt release 2026_04 is the count-only response already stored in Phase 0. It is recorded as a current resource and is not joined to domains.

Version suffixes stay attached to the original identifiers. A database join may also store the stable identifier and the integer version. Those are additional columns.

## Sequence comparison

The sequence method is exact string equality on the GENCODE v34 translation. Identity is 1 when the strings are equal. When they differ, the audit records the two lengths and the absolute length difference as a lower bound on residues that cannot occupy the same position. It does not report a residue-alignment identity, and it does not call that length difference an alignment.

A predicted sequence difference requires both transcripts to be `protein_coding`, both translations to be present, and the strings to differ. An NMD biotype can carry a translation string in the FASTA. That string stays at evidence Level 0. The audit does not treat it as a detected truncated protein.

## Domain comparison

Pfam accessions come from Ensembl 100 features whose analysis `logic_name` is `pfam`. The accession is the hit name before a dot. Signal peptides are `signalp` features. Transmembrane spans are `tmhmm` features. A coordinate is supported when `1 <= start <= end <= sequence length` on the GENCODE translation of that ENSP.

A principal-versus-alternative domain call is made only when both translations exist, each ENSP has one Ensembl translation row, the Ensembl version equals the GENCODE version, and every Pfam, SignalP, and TMHMM feature on those rows has supported coordinates. Otherwise the pair is `mapping_uncertain` or `no_annotation_available`.

The pair classes are mutually exclusive:

- `both_only`: the Pfam accession sets are equal and nonempty.
- `neither`: both proteins have a valid map and both Pfam sets are empty.
- `principal_only`: at least one accession is limited to the principal protein, and the alternative has no accession of its own.
- `alternative_only`: the reverse.
- `mixed`: each protein has an accession the other lacks.
- `mapping_uncertain`.
- `no_annotation_available`.

A Pfam set difference under the valid-map rule is a predicted domain-composition difference. The report name for that event is `predicted_domain_composition_difference`. The audit does not relabel it as protein dysfunction. Absence of a feature file is `no_annotation_available`, not a domain loss.

Signal-peptide presence and transmembrane-span counts are separate annotation comparisons. They can support evidence Level 2 when the map is valid and the protein sequences differ. They are not catalytic-residue evidence. Catalytic residues are counted only when an analysis logic name is `catalytic_residue` or `active_site`.

## Evidence levels

- Level 0: transcript annotation, including noncoding and NMD biotypes, and any pair whose translation is missing.
- Level 1: two protein-coding translations with different sequences and no supported functional-region difference.
- Level 2: a Level 1 sequence difference plus a supported Pfam, signal-peptide, or transmembrane difference.
- Level 3: an isoform-specific peptide or other protein identification that distinguishes the isoform. Gene-level CPTAC abundance is not Level 3. The stored SEPEP matrices have no peptide sequences and no protein accessions, so this audit assigns no Level 3 events.

## Phase 1A linkage

The result universe is every gene–cancer–transformation row with `fit_status = OK`. The candidate universe is the frozen 10,795-gene set, counted once. The intersection is the set of genes with OK fits under both literal inverses. Partial R² is not read and is not a filter. The intersection count must match `n_paired` in `phase1a/results/common_support_summary.csv`. The OK count must match `n_success` in `phase1a/results/eligibility_by_cancer_transform.csv`.

## Decision thresholds

These thresholds were set before the coverage fractions were known.

A. Transcript annotation completeness

- PASS: the GTF hash matches, the ten matrices share one transcript set, and every historical count matches: 108,777 rows, 60,669 genes, 42,424 single-row genes, 18,245 multi-row genes, 0 unmapped identifiers, 0 duplicate identifiers, 14,620 unique MANE genes, 2,819 unique APPRIS principal-1 genes, 1,001 tagged genes without a unique principal, 17,439 principals, 37 principals absent from the matrix, 10,795 candidate genes. The frozen principal table matches gene by gene.
- FAIL: any of those checks fail.

B. Transcript-to-protein mapping

- PASS: at least 99% of protein-coding transcripts in the candidate measured universe have one protein identifier, at least 99% of those identifiers have a GENCODE v34 translation, and no transcript is attached to two different translation strings.
- PARTIAL: the sequence coverage is at least 80% and below 99%, and the mapping classes can still be counted.
- FAIL: sequence coverage is below 80% or the translation file is absent.

UniProt accession coverage is reported beside the sequence result. Incomplete SwissProt coverage does not, by itself, move a PASS sequence map to PARTIAL.

C. Protein domain annotation

- PASS: at least 80% of candidate principal proteins with a translation have one version-matched Ensembl 100 row, and at least 95% of the Pfam features used on those version-matched rows have supported coordinates.
- PARTIAL: the version-matched fraction is at least 40% and below 80%, and uncertain pairs are counted separately.
- FAIL: the feature tables are absent or the version-matched fraction is below 40%.

D. Phase 1A linkage

- PASS: the OK counts and the paired counts match the frozen Phase 1A tables, and the eight annotation counts are present for both inverses, for the intersection, and for the candidate universe.
- FAIL: a count cannot be joined or a frozen count disagrees.

E. Independent validation feasibility

- PASS: a public cohort has cancer and control groups, patient-level identifiers, transcript-level quantification, an annotation compatible with GENCODE v34, at least 50 cases, and no sample overlap with the 1,027 CPTAC cases.
- PARTIAL: at least one public accession is documented with cancer type, platform, a sample count, and an explicit independence statement, while the PASS conditions are unmet.
- FAIL: no accession can be named.

F. Novelty

- PASS is not available when the audited literature already contains pan-cancer isoform-switch or CPTAC proteogenomic studies. Those studies make the qualification PARTIAL when their estimand is different from the frozen partial-R² estimand.
- FAIL: a paper reports the same estimand, the same principal rule, and the same covariates.
- Unknown is reserved for a question the search did not cover.

## Overall label

- `STOP_CURRENT_DIRECTION` when A or D is FAIL.
- `PHASE2B_REQUIRES_REDESIGN` when B is FAIL.
- `QUALIFIED_FOR_PHASE2B` when A, B, C, and D are PASS. This label means a transcript-to-protein structural annotation study can be specified. It does not qualify a functional-loss study.
- `QUALIFIED_FOR_NARROW_PHASE2B` when A and D are PASS, B is PASS or PARTIAL, and C is PARTIAL. The narrow object is predicted sequence difference, with domain coordinates left unresolved.
- Any other combination is `PHASE2B_REQUIRES_REDESIGN`.

Phase 2B authorization stays false under every label.

## Prohibitions

No differential expression, no differential transcript usage, no isoform-switch discovery, no gene ranking, no pathway enrichment, no survival analysis, no machine-learning model, no protein-functional-loss claim, and no use of gene-level protein abundance as isoform-specific evidence.
