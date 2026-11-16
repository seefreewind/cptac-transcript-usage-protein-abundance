# Phase 1A primary results

Scientific classification: CANCER_SPECIFIC_OR_HETEROGENEOUS

Cancer-level means are unweighted across successful genes. Intervals are delete-one-case jackknife intervals on the raw partial-R² scale.

| Cancer | Transform | Successful genes | Mean partial R² | 95% interval | Proportion ≥ 0.02 |
| --- | --- | ---: | ---: | --- | ---: |
| BRCA | A | 5988 | 0.02678968485 | 0.0239790817, 0.029600288 | 0.3648964596 |
| CCRCC | A | 5415 | 0.02105973612 | 0.01947129008, 0.02264818216 | 0.2969529086 |
| COAD | A | 4010 | 0.02039011042 | 0.0179595204, 0.02282070043 | 0.3077306733 |
| GBM | A | 6216 | 0.02813654751 | 0.02560718609, 0.03066590894 | 0.3687258687 |
| HNSCC | A | 5710 | 0.02551038232 | 0.02375186519, 0.02726889945 | 0.3406304729 |
| LSCC | A | 6577 | 0.02807236922 | 0.02611839053, 0.03002634792 | 0.3697734529 |
| LUAD | A | 6263 | 0.0223158786 | 0.01925972323, 0.02537203397 | 0.3151844164 |
| OV | A | 4317 | 0.02709222564 | 0.02426902422, 0.02991542705 | 0.3699328237 |
| PDAC | A | 5282 | 0.01973622736 | 0.0179664341, 0.02150602063 | 0.2904202953 |
| UCEC | A | 5549 | 0.02924008814 | 0.0272035081, 0.03127666819 | 0.3907010272 |
| BRCA | B | 5906 | 0.02628626261 | 0.02322173248, 0.02935079275 | 0.3603115476 |
| CCRCC | B | 5375 | 0.02079917613 | 0.01929208323, 0.02230626903 | 0.2917209302 |
| COAD | B | 3949 | 0.02016812512 | 0.01768995638, 0.02264629385 | 0.3069131426 |
| GBM | B | 6171 | 0.02643134322 | 0.02421554325, 0.02864714318 | 0.3527791282 |
| HNSCC | B | 5659 | 0.02423849734 | 0.0226061458, 0.02587084887 | 0.3304470755 |
| LSCC | B | 6528 | 0.02663410155 | 0.02473048117, 0.02853772194 | 0.3544730392 |
| LUAD | B | 6224 | 0.02149922433 | 0.01875490084, 0.02424354783 | 0.307840617 |
| OV | B | 4230 | 0.02617590985 | 0.02342760485, 0.02892421484 | 0.368321513 |
| PDAC | B | 5246 | 0.01929293663 | 0.01758769001, 0.02099818325 | 0.2893633244 |
| UCEC | B | 5498 | 0.02836682291 | 0.02657527672, 0.0301583691 | 0.3861404147 |

Cross-cancer REML / Hartung–Knapp summaries are in `results/random_effects_summary.csv`. OV remains in its own rows above.

All-ten pooled means: Literal A 0.02481 (95% CI 0.02222 to 0.02740; 95% prediction interval 0.01642 to 0.03320) and Literal B 0.02396 (95% CI 0.02162 to 0.02631; 95% prediction interval 0.01642 to 0.03151). The paired-support pooled intervals also lie above 0.02. Removing OV leaves the pooled intervals above 0.02. The prediction intervals include values below 0.02, which is why the frozen rule classifies the result as CANCER_SPECIFIC_OR_HETEROGENEOUS.

Literal B pre-fit counts on the analysis inputs match the Phase 0.5 eligibility table in every cancer: 5933, 5385, 3971, 6189, 5678, 6548, 6239, 4250, 5260 and 5512. Rank-deficient fits and fits with no residual variation in the principal-transcript proportion are excluded from these means and are counted in `results/cancer_level_summary.csv`.
