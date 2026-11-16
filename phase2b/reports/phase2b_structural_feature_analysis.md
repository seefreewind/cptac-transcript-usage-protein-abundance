# Structural feature analysis

Pfam, SignalP, and TMHMM indicators can overlap. Their coefficients are secondary and are not independent discoveries.

Domain truncation and domain-boundary change are not estimated. A residue alignment was not qualified. The feature summary counts an accession present on only one of the two compared sequences.

- pfam_diff A: genes 3626, adjusted difference 0.005175574673915991, interval 0.0035447162086513877 to 0.006806433139180595.
- signal_diff A: genes 810, adjusted difference 0.00047766030576042406, interval -0.0021875098387381165 to 0.0031428304502589643.
- tm_diff A: genes 1603, adjusted difference 0.0012223276310379158, interval -0.0005350849949098135 to 0.0029797402569856453.
- pfam_diff B: genes 3626, adjusted difference 0.005405418295081263, interval 0.0038007907800195317 to 0.007010045810142994.
- signal_diff B: genes 810, adjusted difference -0.0002067632335526585, interval -0.0027309522222044706 to 0.0023174257550991535.
- tm_diff B: genes 1603, adjusted difference 0.001295899083783178, interval -0.00043131581641067 to 0.0030231139839770257.

The Pfam indicator is the secondary feature whose interval lies above zero under both inverses. The signal-peptide and transmembrane intervals include zero. A gene can carry more than one feature, so these three coefficients are not separate tests of independent biology.

Among S2 genes, 3,548 have a Pfam accession present on the principal sequence and absent from the compared alternative, and 291 have an accession in the opposite direction. Signal peptide presence is principal-only for 663 genes and alternative-only for 147. A predicted transmembrane span is principal-only for 732 genes and alternative-only for 166. These are annotation counts. They are not experimental domain loss, catalytic-site loss, or demonstrated localization changes.
