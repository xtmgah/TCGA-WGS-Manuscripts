# First runnable version: scientific scope

All 49 current lettered panels in main Figures 1–5 are drawn from packaged numeric/categorical inputs. All five pages are composed at the manuscript dimensions. Reference artwork is used only by verification. The panel-level visual review is complete (see `VALIDATION.md`), with small quantified rendering residuals. This does not certify pixel-identical reproduction or public release readiness.

| Figure | Recreated in this package | Frozen upstream results |
|---|---|---|
| 1 | Ordering/density drawing, subtype composition, gain-timing atlas/curves and permutation test, chronology distributions and 9 pairwise tests, BIC display | Ordering samples/KDE, gain-timing and chronology estimates, mixture fits |
| 2 | Ordering, KM/log-rank and Cox fitting, CIC/onco plot aggregation, TP53 dosage regression, chronology tests/BH, plots, volcano | Molecular calls, dosage/timing estimates, adjusted DESeq2 results |
| 3 | Full 28-test inventory, TERT-promoter Fisher test, 12-feature standardized effects and Wilcoxon/BH tests, timing permutation, CT/CN matrix rendering | Genomic measurements, CT calls and timing model estimates |
| 4 | Ordering, KM/log-rank and Cox fitting, displayed Cox contrasts, gain-timing permutation, all 18 metric comparisons and plots | Ordering/timing estimates |
| 5 | Current AC2 architecture, oncogene burden/prevalence/co-occurrence and copy-number aggregation; Fisher/Kruskal-Wallis/Wilcoxon/BH tests; EGFR expression/fusion tests; RNA state composition and adjusted multinomial LRT | AC2 2.0.0 feature classification, normalized expression/fusion calls, GSEA results |

The numerical regression checks use the existing project's outputs as expected values. A passing check establishes agreement with that project, not an independent scientific validation of its models. The AC2 check includes 28 regenerated tables and two checked expression/fusion input tables; all EGFR tests are recomputed. See [STATISTICAL_RECOMPUTATION.md](STATISTICAL_RECOMPUTATION.md) for methods, cohorts, exact correction families and the permutation conventions.

Keep cohort denominators distinct:

| Analysis | Denominator |
|---|---|
| Pan-glioma ordering | 814 specimens; groups 406 / 251 / 157 |
| ASTRO ordering | 172 C17p / 91 CTR specimens |
| ASTRO survival | KM 252 patients, 59 deaths; adjusted Cox 228 patients, 50 deaths |
| Figure 3 CT/CN matrix | 88 primary CT-positive patients: 69 / 19; 704 cells, 2 incomplete CN calls |
| GBM ordering | 399 specimens: 160 / 134 / 105 |
| GBM survival | 382 patients, 277 deaths |
| Figure 5 architecture | 389 primary patients: 159 / 129 / 101; includes 399 GBM specimens assigned to those patients |
| Figure 5 RNA state | 205 patients: 85 / 66 / 54 |

Figure 5 preserves the current AC2 update's different filters: passing features in architecture and EGFR copy-number panels; candidate ecDNA including LowCN in oncogene panels b–d. Patient and specimen denominators are not interchangeable. Reclassification, RNA cohort choice and correction families are retained.

The page sizes are 180 × 242.833333 mm, 180 × 250 mm, 180 × 200 mm, 180 × 250 mm and 180 × 252 mm. Figure 3 has a–h and Figure 5 a–i. Older manifests or drafts must not override these values.

The data-minimization revision uses 105 inputs: 70 Git figure-support/style files and 35 author-request files. It removes unused inputs/columns while preserving the numerical checks and reviewed figure pixels. See [DATA_DISTRIBUTION.md](DATA_DISTRIBUTION.md); Box/Biowulf are author-held and have no public download link.
