# Methodology and manuscript review

These items require an author decision before treating the code package as a fully reconciled version of the final manuscript. They distinguish what the recovered source computes from what the manuscript says. Reproducing a historical plot does not resolve a conflicting statistical description. The source formulas and selection rules are retained; this document does not authorize changing the manuscript or rerunning a different methodology.

Historical source names and line ranges refer to the audited analysis scripts recorded in the provenance table. DOCX paragraph numbers refer to the extracted paragraph inventory of `PTC_Evolution_Manuscript_MainText.SJC.docx`; they are evidence locators, not manuscript section numbers.

## 1. Multiple-testing method: Holm versus Benjamini–Hochberg

**Observed code.** Bare `p.adjust(p.value)` calls use R's default Holm correction, regardless of a column being named `FDR`.

| Figure/source | Computation recovered | Manuscript description | Review action |
| --- | --- | --- | --- |
| Fig2d, `TP_bb_heatmap_freq.R` 439–511 | Fisher tests; Holm-adjusted values determine labels, but the source y axis uses raw P | Fig2 legend paragraph 143 says −log10(FDR) and labels FDR<0.05 | Confirm the final plotting/statistics version; reconcile both the method and the y-axis quantity. |
| Related MATH analyses, `TP_Comparision.R` 2474–2513 | Wilcoxon tests; Holm within Study2; raw-P y axis in cohort plots | Methods paragraphs68/85 say BH; combined Fig3b legend says FDR | Recover exact combined-panel code and distinguish raw-P display from adjusted-P decision rules. |
| Fig4c/d, `TP_Comparision.R` 2260–2336 | Fisher tests and a Holm-valued `FDR` column; plotted raw P | Methods paragraph 77 says BH | Specify actual adjustment for the reported inferential claims; do not infer BH from the variable name. |
| Fig5f, `TP_Comparision.R` 3495–3553 | Regression effects with Holm-adjusted P on the y axis | Methods paragraph 81 and legend paragraph 146 say BH/FDR | Historical effects and adjusted values reproduce exactly. Decide whether to describe Holm or deliberately reanalyze with BH. The cleaned plot identifies its axis as Holm-adjusted P. |

Fig5a differs: its source explicitly applies BH within predictor/effect type (`TP_Comparision.R` 3189–3193). Extended Data Fig4a/b also use explicit BH-adjusted pairwise Wilcoxon tests. They should not be changed to Holm for uniformity.

## 2. Fig4 upper threshold is a fixed Bonferroni-style cutoff

**Observed code:** the upper line is `-log10(0.05/12)` at `TP_Comparision.R` 2295 and 2336. This is a fixed raw-P cutoff of approximately 0.004167 for 12 tests. The y axis is −log10(raw P).

**Manuscript:** Fig4 legend paragraph 145 calls the upper line “FDR <0.05.”

**Review action:** decide whether the upper line should be described as the fixed 0.05/12 cutoff, or whether a different intended final figure/statistical rule exists. An FDR procedure does not generally correspond to that fixed line. Confirm the intended family of 12 tests independently of the number of displayed categories.

## 3. MATH input is CCF in the traced plots, while the manuscript describes VAF

**Observed code:** `TP_Comparision.R` 2366–2387 computes both `MATH_CCF` and `MATH_VAF`. The direct related export at 2398 and driver analyses at 2468–2534 use `MATH_CCF = 100 × mad(CCF)/median(CCF)`, including R's default MAD scaling factor 1.4826.

**Manuscript:** Results paragraph21 and Methods paragraph67 define the score using somatic variant-allele fractions (VAF). CCF adjusts for tumor purity/copy number and is not the same input variable as VAF. Exact final ridgeline and combined-cohort Fig3a–c generators remain partially identified; the related files alone cannot prove that the final panels use the same object.

**Review action:** locate the final Fig3a–c code and verify the intended score definition, MAD convention and inputs. Then reconcile Methods, labels and interpretation. Do not silently rename CCF-derived dispersion as VAF-derived MATH or substitute a new VAF calculation.

## 4. Fig2c contains segment observations, not one independent point per tumor

**Observed code:** `TP_bb_heatmap_freq.R` 519–538 filters qualifying subclonal-loss segments without deduplicating tumors. The current input snapshot yields 69 qualifying segment rows from 63 unique tumors.

**Manuscript:** Fig2c legend paragraph 143 says each point represents one tumor.

**Review action:** confirm whether segment-level sampling was intended. If so, state the observation unit and the 69/63 distinction. If one point per tumor was intended, an explicit selection or aggregation rule is needed and must be reviewed before changing the analysis. The reproducibility package preserves the original rows.

## 5. Fig5a model denominator and exposure coding need explicit description

**Observed code:** `TP_Comparision.R` 3046–3058 creates the alteration matrix from tumors with at least one curated alteration. The all-cohort `right_join` is commented out. The model therefore does not automatically include every cohort tumor as an alteration-negative control. The cleaned implementation in [Fig5.R](../Rscripts/Fig5.R) preserves this denominator; independent original-expression replay produces identical coefficients and P values.

The logistic block caps dose at 1000, fills missing dose with 0 and missing age at exposure with 100 (3180–3183). Latency blocks first fill missing values in matched clinical rows with 0 and only later give unmatched subjects the age-at-exposure default 100 (3453–3459 and 3508–3514). The manuscript's reference-category description in paragraph85 does not fully describe these distinctions. Original and cleaned logistic execution also emit fitted-probabilities 0/1 warnings for two models; this may reflect sparse/extreme fitted strata.

**Review action:** document the actual modeling denominator and missing-data/reference coding; inspect sparse model strata before interpreting individual associations. No logistic estimator or exclusion rule was changed during code cleaning.

## 6. Fig6 highlighted tumors are the BRAF–TERT intersection

**Observed code:** `Segmented_breakpoint_analysis.R` 10–11 selects TERT mutation-driver tumors and then restricts to BRAF mutation-driver tumors within that set. The highlighted points are therefore BRAF-plus-TERT co-mutants. The reproduced plot highlights 36 TCGA and 1 Chornobyl tumors.

**Manuscript:** Fig6 legend paragraph 147 describes all TERT-promoter-mutant tumors as highlighted; Results paragraph 35 similarly discusses TERT-mutant points.

**Review action:** confirm the intended highlighted population, then align the legend and interpretation with the executed selection. The [Fig6.R](../Rscripts/Fig6.R) code retains the intersection.

The user confirmed that panels 6a and 6b are the left TCGA and right Chornobyl regions of the same imported `tp_age_clock_mutation_acc.pdf`, with panel letters added in Keynote. This is a verified mapping and requires no correction.

## 7. Fig6b slope label and drawn line use different adjustment sets

**Observed code:** the Chornobyl slope label comes from `lm(value ~ age_at_diagnosis + Tumor_Purity + Sex)`; the displayed `geom_smooth(method="lm")` fits value against age alone. These produce slopes 4.581181 and 4.632625 in the current inputs. Both round to 4.6, masking the distinction. TCGA predictions are explicitly evaluated at mean purity and modal sex.

**Manuscript:** Methods paragraph 83 specifies the adjusted Chornobyl model, and Fig6b describes a fitted regression line with uncertainty.

**Review action:** decide whether the Chornobyl line/band should display adjusted predictions or the caption should distinguish the descriptive unadjusted line from the adjusted coefficient. The historical rendering mismatch is preserved and documented; it is not silently corrected.

## 8. Extended Data Fig2 and Supplementary Fig3 reuse the same plots

**Observed evidence:** Keynote slides 9 and 14 reuse identical media IDs 12560 and 12567, corresponding to the same two cohort-specific MATH PDFs. Both sets of manuscript identifiers remain in the provenance table.

**Manuscript:** paragraph 23 cites Extended Data Fig2 for minimum detected CCF, although its legend paragraph 150 describes MATH and minimum CCF is in Supplementary Fig2. Extended Data Fig2b describes points/boxes; Supplementary Fig3b describes violins/boxes for the same source visual.

**Review action:** decide whether the repeated figures are intentional, fix the minimum-CCF cross-reference, and harmonize the descriptions of the duplicated panels. Nothing has been dropped from the package solely because it is duplicated.

## 9. Age-bin endpoints and APOBEC thresholds must remain explicit

**Observed code:** age groups are <20, [20,30), [30,40), [40,50), [50,60), [60,70], and >70 (`TP_Comparision.R` 30–41). A tumor diagnosed at exactly 70 is in the 60–70 group. TERT frequencies use all records with Gene==TERT, whereas BRAF age frequencies require Type==Mutation_Driver. All tumors within each cohort/age group provide the frequency denominator.

APOBEC positivity in Extended Data Fig4c requires SBS2+SBS13 strictly greater than 50 assigned mutations **and** strictly greater than 0.05 of the total (1640–1655). It is not an inclusive ≥50 or ≥5% rule.

**Review action:** retain these endpoint/operator definitions in input documentation and any final methods detail. Current reproduced APOBEC percentages agree with paragraph 33 and the legend (9.1%, 11.9%, 21.8%, 39.0%).

## 10. Chronological interpretation depends on a partially traced upstream chain

**Observed code:** the principal Fig5/Extended Data Fig4 plots use `MRCAdata` filtered to acceleration==`1x`. The stored intermediate was generated by `MutationTime_Analysis_thyroid.R` 1400–1470. Reproducing its downstream plots does not rerun the mutation-cluster/rate/simulation chain. Exact code for the imported Fig5g MRCA-driver regression was not located, despite recovering its original PDF and nearby latency-model preparation.

**Manuscript:** paragraph 81 states cross-scenario consistency, while paragraphs 32/41 interpret a roughly 24-year later MRCA association. Those claims require their own source/model evidence; a recovered image does not supply it.

**Review action:** provide the final Fig5g generator and the saved scenario-comparison analysis. Preserve the distinction between inferred population ancestry/latency and the acquisition time of TERT; the manuscript already states that TERT mutation timing is not directly inferred.

## 11. Additional analysis-state and missing-data decisions

- **Fig2a clone/subclone states:** the recovered script assigns both plot objects from the same current state, while the alternate CN-state choice is commented (`TP_bb_heatmap_freq.R` 55–56, 331–333). Final state-specific execution must be recovered before claiming the distinct displayed profiles reproduce.
- **Extended Data Fig1 missingness:** the burden preparation replaces missing values by0 (`TP_Comparision.R` 623–651). This may conflate unavailable measurements with genuine zero burden; preserve it as source behavior and verify the input-level meaning before interpretation.
- **Fig3e accepted cluster count:** related subclone-count code does not explicitly apply the manuscript's at-least 50-mutation condition. Recover the exact driver-stratified generator before treating its counts as validated.
- **Fig2d self-association filtering:** source removes Gene==chr22q but can retain the curated chr22q_Clone/chr22q_Subclone categories. Review whether such rows belong in the reported event comparisons.

## Changes made during reproducibility packaging

Code cleaning made paths portable, isolated input loading, declared packages, removed debugging/clipboard actions, set Roboto Condensed and legible output dimensions, and fixed namespace collisions. The Fig5f axis now names its actually computed Holm-adjusted P value; the numbers are unchanged. A recorded seed stabilizes segmented bootstrap restarts, and the historical fitted optimum reproduces. No alternative statistical model, tumor selection rule, genome build or data interpretation was substituted to resolve the review items above.
