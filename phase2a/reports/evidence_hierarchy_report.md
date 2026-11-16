# Evidence hierarchy

Level 0 is transcript annotation. Level 1 is a predicted protein-sequence difference. Level 2 is a predicted functional-region difference on supported coordinates. Level 3 is independent isoform-specific protein evidence.

- LEVEL_0: 4861 pair events. Transcript annotation, a missing translation, an NMD biotype, or an identical pair of translations.
- LEVEL_1: 14811 pair events. Two protein-coding GENCODE v34 translations differ, and no supported functional-region difference is called.
- LEVEL_2: 9050 pair events. A Level 1 sequence difference plus a supported Pfam, signal-peptide, or transmembrane difference on version-matched coordinates.
- LEVEL_3: 0 pair events. An isoform-specific peptide or other protein identification distinguishes the isoform. Gene-level abundance does not qualify.

Of the Level 0 events, 4834 are identical protein-coding translations and 27 are matrix alternatives without a paired protein-coding translation.

The current SEPEP files are sample-level quantitative groupings. Six accessible files were assessed in Phase 0.5, four returned HTTP 403, and none of the accessible files contain peptide sequences or protein accessions. This audit assigns no Level 3 events and does not promote gene-level CPTAC abundance to isoform evidence.
