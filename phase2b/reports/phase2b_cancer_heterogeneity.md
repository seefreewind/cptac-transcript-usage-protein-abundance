# Cancer-specific structural contrasts

Phase 1B left the partial-R² heterogeneity as mixed measurement and cancer-specific structure. These contrasts ask whether the S2-versus-S1 difference itself varies by cancer. They do not assign a cause to that variation. Inside one cancer, the protein-observation fraction is a rescaling of the complete-case count, so the rank check keeps only one of those two columns.

| Cancer | Transform | S2 genes | S1 genes | Adjusted difference | Interval |
|---|---|---:|---:|---:|---|
| BRCA | A | 2532 | 2892 | 0.003445 | 0.001156 to 0.005735 |
| CCRCC | A | 2298 | 2626 | 0.003086 | 0.000911 to 0.005261 |
| COAD | A | 1680 | 1967 | 0.001991 | -0.000388 to 0.004370 |
| GBM | A | 2685 | 2938 | 0.005565 | 0.002922 to 0.008208 |
| HNSCC | A | 2408 | 2768 | 0.006052 | 0.003429 to 0.008676 |
| LSCC | A | 2831 | 3117 | 0.008059 | 0.005535 to 0.010582 |
| LUAD | A | 2686 | 2996 | 0.004623 | 0.002593 to 0.006654 |
| OV | A | 1780 | 2109 | 0.003570 | 0.000814 to 0.006327 |
| PDAC | A | 2232 | 2576 | 0.003042 | 0.000992 to 0.005093 |
| UCEC | A | 2325 | 2718 | 0.005354 | 0.002634 to 0.008075 |
| BRCA | B | 2497 | 2857 | 0.003864 | 0.001577 to 0.006150 |
| CCRCC | B | 2283 | 2605 | 0.003012 | 0.000811 to 0.005214 |
| COAD | B | 1649 | 1945 | 0.001803 | -0.000595 to 0.004202 |
| GBM | B | 2669 | 2913 | 0.005974 | 0.003472 to 0.008476 |
| HNSCC | B | 2389 | 2742 | 0.005917 | 0.003393 to 0.008441 |
| LSCC | B | 2807 | 3100 | 0.007571 | 0.005125 to 0.010017 |
| LUAD | B | 2670 | 2979 | 0.004826 | 0.002827 to 0.006825 |
| OV | B | 1740 | 2070 | 0.004235 | 0.001533 to 0.006936 |
| PDAC | B | 2215 | 2562 | 0.002959 | 0.000989 to 0.004929 |
| UCEC | B | 2305 | 2693 | 0.005174 | 0.002533 to 0.007815 |

Leave-one-cancer-out estimates:

- Omit BRCA A: 0.004804 (0.003393 to 0.006215).
- Omit CCRCC A: 0.004807 (0.003391 to 0.006224).
- Omit COAD A: 0.004856 (0.003445 to 0.006267).
- Omit GBM A: 0.004530 (0.003138 to 0.005922).
- Omit HNSCC A: 0.004481 (0.003117 to 0.005845).
- Omit LSCC A: 0.004179 (0.002819 to 0.005538).
- Omit LUAD A: 0.004652 (0.003250 to 0.006053).
- Omit OV A: 0.004740 (0.003340 to 0.006140).
- Omit PDAC A: 0.004812 (0.003400 to 0.006225).
- Omit UCEC A: 0.004569 (0.003208 to 0.005929).
- Omit BRCA B: 0.004813 (0.003428 to 0.006198).
- Omit CCRCC B: 0.004876 (0.003489 to 0.006264).
- Omit COAD B: 0.004921 (0.003540 to 0.006302).
- Omit GBM B: 0.004538 (0.003167 to 0.005910).
- Omit HNSCC B: 0.004555 (0.003213 to 0.005897).
- Omit LSCC B: 0.004304 (0.002967 to 0.005641).
- Omit LUAD B: 0.004685 (0.003307 to 0.006062).
- Omit OV B: 0.004737 (0.003365 to 0.006110).
- Omit PDAC B: 0.004884 (0.003496 to 0.006272).
- Omit UCEC B: 0.004649 (0.003313 to 0.005985).
