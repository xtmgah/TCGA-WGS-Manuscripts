# Recomputed statistics for main Figures 1–5

The build calculates the tractable main-figure statistics from packaged measured or processed values, checks them against the immutable saved results, then passes the newly calculated tables to the renderers. Saved tables supply labels, test specifications and row order, not numerical answers. A mismatch stops the build before figure drawing. Work is confined to the standalone output workspace; the packaged reference inputs are never overwritten.

`functions/R/recompute.R` fits the survival and RNA-state models. `recompute_extended.R` handles chronology, metric tests, TP53 regression, displayed Cox contrasts, feature effects and timing permutations. `ac2_statistics.R` handles AC2 and EGFR tests. `outputs/verification/extended_statistics.json` records methods, correction families, inputs and renderer destinations for 15 regenerated tables or selected table rows. `statistical_verification.json` records column-level numeric comparisons; `ac2_verification.json` records the 30 AC2 table checks. Counts of regression checks are not counts of biological hypotheses.

## Test families and data boundaries

| Figure/panel | Independent calculation | Cohort and correction family |
|---|---|---|
| 1g | Segment-length-weighted mean gain times; range permutation test | 729 specimens: 241 ASTRO, 95 OLIGO, 393 GBM. 6,665 gain segments; the duplicated `category=all` atlas rows are excluded. |
| 1h–j | 9 two-sided Wilcoxon rank-sum comparisons | BH separately for the 3 subtype pairs within each metric. The packaged 1x chronology is retained. MRCA/latency counts: 147/107/251; diagnosis-age counts: 262/148/397 (ASTRO/OLIGO/GBM). |
| 2c–d | KM/log-rank; adjusted Efron Cox model; displayed HR orientation | KM: 252 patients, 59 deaths. Cox: 228 patients, 50 deaths. Group and sex contrasts are reversed for display, with reciprocal HR and swapped reciprocal CI. No multiplicity correction. |
| 2f | OLS mutant-copy dosage on total CN centered at 2, group and interaction | 226 unique primary patients: 163/63. Uses **uncapped** copy numbers; the display cap of 5 does not enter the model. Four terms, ordinary standard errors, 95% t intervals. |
| 2g–h | Wilcoxon and Welch tests, medians, means, quartiles and group differences | 254 primary patients; metric-specific finite observations. Age: 167/86; MRCA and latency: 95/47. BH separately across **all 3 metrics**, including nondisplayed latency, for each test type. Packaged 2.5x chronology is retained. |
| 3a–c, e–f; complete saved test inventory | 2 Fisher and 26 Wilcoxon tests, effects, intervals and summaries; separate TERT-promoter Fisher test | 263 specimens: 172/91, with finite observations per continuous metric. BH across the 28-test inventory; the separate TERT-promoter test is outside that family. Panels retain their original nominal or adjusted labels. |
| 3d | Weighted gain-time difference and permutation test | 241 specimens: 161/80. Uses the ASTRO events already packaged for Figure 1, joined by specimen to the Figure 3 group roster. |
| 3h | 12 pooled-SD standardized mean differences and Wilcoxon tests | Recorded identity/log10(x+1) transformations; BH across the 12 features for display. Historical nominal significance columns are also independently regenerated. |
| 4d–e | KM/log-rank; Efron Cox model; displayed HR orientation | 382 patients, 277 deaths. Two group contrasts are reciprocated for display; the other covariates retain the fitted orientation. |
| 4g | Weighted gain-time means and range permutation test | 392 specimens: 157/132/103 (DN1/DN2/DN3). Uses specimen-level weighted sums in external input 169. |
| 4h–k; complete six-metric inventory | 18 Wilcoxon comparisons across MRCA age, diagnosis age, latency, MATH, PGA and ploidy | BH separately for the 3 group pairs within each metric. Includes nondisplayed MRCA and age tests. Finite observations and saved group assignments are retained. |
| 5a–g and associated AC2 tables | Architecture/count summaries; Fisher, Kruskal–Wallis and Wilcoxon tests; EGFR CN, expression and fusion tests | Preserves the existing AC2 filters and patient denominators. EGFR expression uses 208 curated values. Fusion tests use the 3×3 class-count table. BH across 3 pairs separately for each EGFR metric. Other AC2 correction families are unchanged in `ac2_statistics.R`. |
| 5i | Multinomial likelihood-ratio test, full versus reduced model | 205 patients: 85/66/54. Full model: dominant RNA state ~ sequencing center + standardized Battenberg purity + group; reduced model omits group. No multiplicity correction. Only the displayed row is retained in the minimized RNA reference table; the 79 unused upstream rows were removed during data review. Composition counts are checked directly against the same metadata. |

Wilcoxon tests use `exact=FALSE` with R's continuity correction. Quartiles use R's default type 7. The Figure 3 binary effects are sample odds ratios with log-Wald intervals, adding 0.5 to all four cells if any cell is zero; the P value comes from Fisher's exact test. These are the manuscript's conventions, including the distinction between this odds ratio and Fisher's conditional maximum-likelihood estimate.

The two EGFR measurement/class-count tables remain checked copies: this package does not rerun expression normalization or fusion calling. Their **tests** are now recomputed. All 28 other AC2 tables are regenerated. DESeq2 volcano statistics and GSEA enrichment statistics remain frozen with their upstream fits, as do ordering/mixture fits and molecular/chronological timing estimates.

## Timing permutations

All use 5,000 iterations with explicit Mersenne-Twister / Inversion / Rejection RNG settings. Permuting whole specimens preserves their segment weights. The sufficient statistics are `sum(time × segment_length)`, `sum(segment_length)` and segment count; no sequencing reads or upstream timing-model fits are needed.

| Panel | Seed and sampling scheme | Observed statistic | Exceedances | Original P convention |
|---|---|---:|---:|---|
| 1g | 101; shuffle subtype labels across specimens | 0.3339714671332766 | 0 | (0 + 1) / 5001 = 0.0001999600079984003 |
| 3d | 11; select 161 specimen IDs for group 1 each iteration | −0.1061437916083917 | 130 | 130 / 5000 = 0.026; two-sided absolute difference |
| 4g | 11; shuffle DN labels across specimens | 0.4432993653367281 | 0 | 0 / 5000 = 0 |

The Figure 4 zero records no exceedances in 5,000 trials; it does not establish a zero underlying probability. The original correction conventions are retained for reproduction rather than silently replacing them. Changing them would be a manuscript analysis revision.

## Validation contract

P and adjusted P values use relative tolerance `1e-8` with **zero absolute tolerance**. This prevents, for example, a P value near `1e-140` from being accepted as zero. Other numeric fields use absolute tolerance `1e-10` plus relative tolerance `1e-8`. Missing-value locations, unique comparison keys, complete schemas, table sizes and categorical annotations must match. Expected tables remain protected by the input manifest's SHA-256 checksums.

Run the focused failure tests with:

```sh
Rscript --vanilla tests/test_statistics.R
```

The 12 tests cover tiny-probability corruption, exact-zero references, missingness/nonfinite values, duplicated or missing comparisons, key-based alignment, BH family boundaries, the permutation plus-one convention and duplicate-specimen rejection. Synthetic test artifacts go under ignored `outputs/test-statistics/`.

Then run the normal complete build and verification:

```sh
Rscript Rscripts/reproduce.R --figures 1-5 --verify --render-reports
```

Passing regression checks establishes independent recalculation from the packaged inputs and agreement with the original analysis. It does not independently establish the validity of upstream classifications, measurements, statistical assumptions or expensive frozen fits. File-level Git/request routing and minimization are documented in DATA_DISTRIBUTION.md. Study-specific access permissions, environment portability and the scientific release license are separate from statistical agreement.
