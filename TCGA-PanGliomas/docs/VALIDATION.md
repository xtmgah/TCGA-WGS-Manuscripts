# Validation of the standalone reproduction

Validated locally on 2026-10-06. The complete five-figure build succeeds from **105 packaged inputs**. The statistical recomputation, panel-level visual review, data-classification/minimization and fresh macOS R/Python environment milestones are complete. The freshly installed environment's drawing/statistics phase took **101.21 seconds**, excluding HTML reports and initial dependency setup. A second platform is deferred at the author's request.

## Passed checks

- **A freshly installed macOS arm64 environment passes.** R 4.5.1, Python 3.12.13, Pandoc 3.12 and all R/Python dependencies were installed into a new prefix. The original project, primary package, old R framework, Homebrew and old bundled Python runtime were blocked during the isolated run; reference artwork was also blocked during drawing. All five pages match the previously reviewed pixels exactly at 150 and 300 dpi. All six HTML reports render, and the executable Figure 1 R Markdown entrypoint independently rebuilds with a new font cache and matches the reviewed pixels. See [ENVIRONMENT_VALIDATION.md](ENVIRONMENT_VALIDATION.md), `environment_checks.json` and `environment/validated-macos/`.
- **All 14 Python tests, 12 statistical R tests and two entrypoint-isolation checks pass in that fresh environment.** New tests verify Cairo's empty-space font normalization preserves pixels and text and never substitutes painted glyphs. The font-selection test prevents a competing installed font from changing Matplotlib's selected files. The explicit-Python selection test prevents accidental reuse of the old convenience dependency directory. No scientific tolerance or visual-review threshold was loosened.
- **141 numerical regression checks covering 1,064 numeric entries (including expected missing values) pass.** They include all chronology/metric pairwise tests and BH families, TP53 regression, full Figure 3 tests/effects, displayed Cox contrasts, three seeded timing permutations, survival fits and the displayed adjusted RNA-state multinomial test. Newly computed values replace the materialized renderer inputs only after verification. See [STATISTICAL_RECOMPUTATION.md](STATISTICAL_RECOMPUTATION.md).
- **All 30 AC2 tables match** the saved release. Of these, 28 are regenerated; two are checked expression/fusion measurement tables. EGFR expression and fusion tests are now independently recomputed. Candidate/passing filters and patient/specimen denominators are preserved.
- **All 12 focused statistical regression tests pass.** These exercise tiny-P corruption, exact zeros, missing values, duplicate/missing comparisons, BH family boundaries and permutation conventions. P/q comparisons have no absolute tolerance, so a tiny nonzero P cannot be silently accepted as zero.
- All five pages remain vector PDFs with exact manuscript page dimensions and 49 panel letters. **Every page retains the reviewed pixel fingerprint at both 150 and 300 dpi** after expanded recomputation. The combined PDF matches its standalone pages at both resolutions; numerical text remains unchanged.
- **All 169 former inputs have a recorded distribution decision:** 70 Git figure-support/style inputs, 35 author-request inputs, and 64 retired files. The cleanup removed 416 columns from retained tables. The request set is 0.763 MB compressed; no raw sequencing data were added. All **207 tracked original project file hashes** are unchanged. Only standalone copies were minimized, and the previous package is preserved privately. See [DATA_DISTRIBUTION.md](DATA_DISTRIBUTION.md).
- **All seven release-contract tests pass.** The exporter rejects unreviewed files, individual-record columns and participant identifiers in proposed Git data. Git-only input validation works without the private bundle. These checks enforce the chosen routes; they do not certify study-specific source-data permissions.
- The build runs under an OS sandbox denying reads/writes to the original project outside this folder and denying manuscript reference-PDF reads. Verification runs separately with access to the reference copies.
- Both PDF-text regression tests pass. Regenerated PDFs have no MuPDF syntax/rendering warnings. The original Figure 4 reference's empty ligature wrapper remains an upstream artifact and does not affect the compared labels.
- All five R Markdown reports and the data-format report render as self-contained HTML. No software was installed into the original project.
- The input corruption/absence and external-data relocation contract tests pass with the new 105-input manifest. The per-figure R Markdown rebuild smoke test was performed at the first runnable milestone; that unchanged interface retains its earlier evidence in `milestone_checks.json`.
- A fresh extraction of the minimized Git/code archive plus the private author-request archive rebuilt all five figures in 74.22 seconds, with reads of the original project, primary package and manuscript reference PDFs blocked during rendering. All three statistical audits and all five reviewed image fingerprints at 150 and 300 dpi match the primary run exactly. All runtime, dependency, report and input files match the tested export. Git-only validation also passes without the private bundle; a full build correctly reports the 35 missing inputs and author-request instructions. Results are recorded in `data_minimization_checks.json` and `milestone_checks.json`.


## Completed panel review

The review covered full pages, all 49 panel crops, exact pixel-difference maps at 150 and 300 dpi, and enlarged comparisons of residual regions. All **1,556 text spans** match by content, font, size and color. The greatest matched text-origin difference is 0.000153 pt (less than 0.000054 mm). There are no missing or extra numerical tokens, no page-edge text clipping, and no embedded raster artwork. Minimum extracted type size is 6 pt.

| Figure | Before: changed at 150 dpi | After: changed at 150 dpi | After: changed at 300 dpi | After: differences >16/255 at 300 dpi |
|---|---:|---:|---:|---:|
| 1 | 1.030% | 0.00459% | 0.00405% | 0.00000% |
| 2 | 5.115% | 0.00102% | 0.00107% | 0.00000% |
| 3 | 7.106% | 0.00000% | 0.00464% | 0.00000% |
| 4 | 3.604% | 0.00191% | 0.00115% | 0.00088% |
| 5 | 0.863% | 0.02262% | 0.01414% | 0.00043% |

The initial differences were resolved in the rendering code: native Cairo padding, title centering and wrapping, first-row Figure 3 margins, inset and matrix legends, risk-table headings, consistent black cohort labels, and the italic gene headings/fusion key in Figure 5. Points and jitter now agree at panel level. Labels that depend on counts or test results still come from the packaged data and calculations.

The remaining differences are small glyph/path-edge and rasterization residuals, including a few legend-edge pixels. They were inspected and accepted for this visual reproduction milestone. **Exact pixel identity is not claimed**: Figure 3 is exact at 150 dpi but has low-intensity differences at 300 dpi. Changed-pixel fractions are rendering diagnostics, not a measure of scientific accuracy.

`outputs/verification/panel_review.html` shows reference/reproduction/difference views for each panel. Full metrics, text checks and whole-page difference maps are in `outputs/verification/`. `tests/visual_review.json` records the reviewed image fingerprints at both resolutions and the reference PDF hashes. Verification reports `REVIEWED_MATCH` only when the new output has exactly those reviewed pixels; changed output requires renewed review. There is no blanket tolerance that silently approves a future altered plot.

`functions/styles/title_positions.json` contains static title typography/positions only. `functions/python/polish.py` translates or redraws presentation elements on the newly rendered pages. Neither file embeds reference artwork or substitutes saved figure images for data-driven plotting.

## Remaining release work

1. Decide the scientific code/data license, citation and release version.

Second-platform validation is deferred at the author's request. No container test is requested. The fresh-environment recipe and exact installed dependency inventory are recorded; a fully tested cross-platform lockfile remains optional future work.

The execution, visual-review, tractable statistical recomputation, data-routing/minimization and fresh macOS environment milestones are complete. Git is the only public distribution location; request-only data remain author-held on Box/Biowulf. Source/institutional conditions still govern any requested transfer. These checks do not independently validate frozen upstream models.
