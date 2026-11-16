# Figure and table legends

Figure 1. Study design and measurement levels. A, ten CPTAC cohorts with 1,027 case-level matched tumor observations and individual cohort counts. B, published gene RNA, transcript RNA and gene-level protein layers. C, GENCODE v34 principal assignment uses unique MANE Select and a unique APPRIS principal-1 fallback only when MANE is absent, yielding 10,795 candidates. D, within-cancer nested models estimate gene-level partial R², with cancer-level jackknife uncertainty followed by REML synthesis, Hartung–Knapp confidence intervals and prediction intervals. Gene-specific complete-case counts are smaller than cohort intersections. Case-level matching; aliquot-level correspondence was not established.

Figure 2. Nested models and primary pooled incremental explanation. A, G is the sum of reconstructed linear transcript values and P is the principal-transcript value divided by G. The outcome is gene-level protein abundance. M1 adds P to M0 with log₂(G), WES purity and CNV. Partial R² is (RSS_M0 − RSS_M1)/RSS_M0. B, pooled cancer-level means and 95% Hartung–Knapp CIs from ten cohorts. Blue denotes Literal A, red Literal B, and the dashed line marks 0.02. Two co-primary measurement representations use the same underlying tumors. Prediction intervals appear in Figure 4.

Figure 3. Primary gene-level partial R² distributions. A, main-range ECDFs from 0 to 0.10. Full-range ECDFs from 0 to 1 retain the upper tail in Supplementary Figure S2. B, frozen cancer means (circles) and medians (open squares), with the 0.02 reference line. Rows separate Literal A and B. The same cancer colors and order are used throughout. Every successful Phase 1A gene fit contributes, with counts in Table 2. Viewing limits in A leave the full distribution in Supplementary Figure S2. No uncertainty intervals are attached to descriptive ECDFs or medians.

Figure 4. Primary cancer-specific heterogeneity. Panels A and B show Literal A and B. Circles are cancer-level mean partial R² with 95% delete-one-case jackknife CIs. n counts successful gene fits. Diamonds and thick bars give REML pooled means and 95% Hartung–Knapp CIs. Thin grey lines give 95% prediction intervals. The dashed vertical line marks 0.02. The confidence interval describes uncertainty in the pooled mean; the prediction interval describes variation among cancer-level means under the random-effects model. All ten cancers are retained in fixed order. I² is reported as a percentage.

Figure 5. Primary transformation robustness. A, paired gene–case-support cancer means for Literal A and B with an identity line and cancer colors defined in Figure 3. B, all-ten, paired-support and excluding-OV REML means with 95% Hartung–Knapp CIs. C, all ten leave-one-cancer-out analyses for each transformation. B and C use identical horizontal scales and a dashed 0.02 threshold. Blue denotes A and red B. All-ten and paired-support syntheses use ten cancer means, exclusion analyses use nine. Two co-primary measurement representations retain their shared biological observations.

Supplementary Figure S1. Secondary predicted annotation comparison. S2 and S1 definitions appear in A. B shows adjusted S2–S1 differences of 0.004644 (95% CI 0.003269–0.006019) for Literal A and 0.004697 (0.003346–0.006048) for Literal B, with gene-cluster-robust intervals and the 0.005 practical threshold. Both point estimates were below the prespecified practical magnitude of 0.005, while both confidence intervals extended above it. Product-level function requires direct experimental evidence.

Supplementary Figure S2. Full-range ECDFs use all successful gene fits under A and B, colored by cancer as in Figure 3. The display spans zero to one and retains the full upper tail.

Supplementary Figure S3. Complete stored measurement-sensitivity means. Each dot is a stored cancer-level mean, except S8, whose stored row is a pooled excluding-OV mean. Blue denotes Literal A, red B and grey triangles Candidate C. Candidate C is a sensitivity representation. S4 CLR has no plotted estimate because it is recorded as NOT_RUN_UNDERDEFINED. Missing CI fields remain unestimated. The dashed reference is 0.02.

Table 1. Counts describe case-level intersections. CNV header coverage is distinct from gene-specific complete-case coverage. Protein missingness is the fraction of matrix cells, with no patient-level missingness interpretation.

Table 2. Primary within-cancer summaries use every successful gene fit under Literal A or B and at least 71 complete cases per fit. Intervals are 95% delete-one-case jackknife CIs.

Table 3. Primary synthesis of ten cancer-level means. The confidence interval describes uncertainty in the pooled mean; the prediction interval describes variation among cancer-level means under the random-effects model. Practical threshold 0.02.

Table 4. Primary sensitivity summaries retain stored REML means, 95% Hartung–Knapp CIs, 95% PIs and I². All-ten and paired-support analyses use ten cancers, excluding OV uses nine. Both inverses are co-primary measurement representations. Leave-one-cancer-out results are in Supplementary Table S2.

## Supplementary table legends

Supplementary Table S1. Complete frozen primary cancer-level summaries and gene-fit archives for Literal A and B. Rows identify cancer, transformation and successful gene fits. Cancer means and medians retain stored precision. Delete-one-case jackknife uncertainty accompanies cancer means. The separate heterogeneity file retains REML, Hartung–Knapp, prediction-interval and Cochran Q fields. The practical threshold is 0.02. Gene-level archives retain eligibility and model outputs in their original schema.

Supplementary Table S2. Every stored leave-one-cancer-out synthesis for both co-primary transformations. Each row pools nine cancer-level means. Stored 95% Hartung–Knapp confidence limits are retained. Fields absent from the frozen output remain blank. The practical threshold is 0.02.

Supplementary Table S3. Full prespecified measurement-sensitivity record and transformation comparison. Rows identify the cancer, transformation and sensitivity setting. Cancer-level means and pooled summaries retain their original units and status fields. Candidate C is a sensitivity representation. CLR remains NOT_RUN_UNDERDEFINED. Available uncertainty fields are preserved. The practical threshold is 0.02.

Supplementary Table S4. Original gene eligibility and fit-exclusion records for Literal A and B. Rows preserve the frozen cancer, gene or exclusion-reason level. These are model-support counts, and no inferential intervals apply. The primary complete-case minimum was 71.

Supplementary Table S5. Frozen GENCODE v34 principal-transcript annotation and transcript-structure coverage records. Counts describe transcript or gene annotation coverage according to each source field. These descriptive counts have no confidence interval. Both transformations share the reference annotations.

Supplementary Table S6. Frozen structural-annotation coverage and category counts. S2 denotes different predicted sequences with a supported Pfam, signal-peptide or transmembrane annotation difference. S1 denotes different predicted sequences without a supported difference in these features. Counts are descriptive gene-level annotation summaries shared by both transformations.

Supplementary Table S7. Full frozen Phase 2B structural regression outputs under Literal A and B. The secondary adjusted contrast compares S2 with S1 and retains stored gene-cluster-robust standard errors and 95% confidence intervals. The prespecified practical magnitude is 0.005. Original covariates, support counts and status fields are retained.

Supplementary Table S8. Stored null-expected and null-excess partial R² summaries for each cancer and co-primary transformation. Values are copied from the frozen cancer-level summaries. These descriptive fields have no additional confidence interval. The primary inference continues to use raw partial R² and the 0.02 practical threshold.
