## Driver-defined evolutionary routes in papillary thyroid cancer

This package traces the PTC manuscript's **13 figures and 46 panel-level entries**
from the manuscript and final Keynote to original plot filenames, historical R code,
required inputs and cleaned executable code. Its organization follows the repository's
[TCGA-UCEC](https://github.com/xtmgah/TCGA-WGS-Manuscripts/tree/main/TCGA-UCEC) implementation.

- **R scripts:** one module per main, Extended Data and Supplementary figure.
- **R Markdown and HTML:** figure-by-figure source, legends, provenance and limitations.
- **Data:** documented private input schemas and small GRCh38 reference resources.
- **Audit:** every panel has one of the four required reproduction states.

**Download the HTML pages and open them locally.** Public HTML is a self-contained
code/provenance report; no private patient-level data or plot previews are embedded.
The report renders without controlled inputs. Analysis execution and optional
private preview rendering are separate, explicit operations.

### Start from a manuscript panel

Open [figure_provenance.tsv](provenance/figure_provenance.tsv) and find its exact
`manuscript_id`, for example `Fig. 1a` or `Supplementary Fig. 2`. Each row gives:

**manuscript panel → Keynote slide/object → exact original filename → source code
location/hash → inputs → cleaned R module → output → validated status.**

Figure 6a/b share one imported composite PDF; the author confirmed that the two
panel letters were added manually in Keynote. Extended Data Figure 2 and
Supplementary Figure 3 reuse the same two imported PDFs. Figure 4e was native
Keynote artwork and now has an explicitly identified R recreation.

### Figures included

| Figure | Code | HTML |
| --- | --- | --- |
| Fig. 1: Driver landscape | [R](Rscripts/Fig1.R) / [Rmd](Rscripts/Fig1.Rmd) | [HTML](Rscripts/Fig1.html) |
| Fig. 2: Clonal/subclonal 22q loss | [R](Rscripts/Fig2.R) / [Rmd](Rscripts/Fig2.Rmd) | [HTML](Rscripts/Fig2.html) |
| Fig. 3: Intratumor heterogeneity | [R](Rscripts/Fig3.R) / [Rmd](Rscripts/Fig3.Rmd) | [HTML](Rscripts/Fig3.html) |
| Fig. 4: Subclonal architecture | [R](Rscripts/Fig4.R) / [Rmd](Rscripts/Fig4.Rmd) | [HTML](Rscripts/Fig4.html) |
| Fig. 5: Age, exposure and chronology | [R](Rscripts/Fig5.R) / [Rmd](Rscripts/Fig5.Rmd) | [HTML](Rscripts/Fig5.html) |
| Fig. 6: Clonal mutation clock | [R](Rscripts/Fig6.R) / [Rmd](Rscripts/Fig6.Rmd) | [HTML](Rscripts/Fig6.html) |
| Extended Data Fig. 1: Burdens | [R](Rscripts/ExtendedDataFig1.R) / [Rmd](Rscripts/ExtendedDataFig1.Rmd) | [HTML](Rscripts/ExtendedDataFig1.html) |
| Extended Data Fig. 2: MATH associations | [R](Rscripts/ExtendedDataFig2.R) / [Rmd](Rscripts/ExtendedDataFig2.Rmd) | [HTML](Rscripts/ExtendedDataFig2.html) |
| Extended Data Fig. 3: DPClust examples | [R](Rscripts/ExtendedDataFig3.R) / [Rmd](Rscripts/ExtendedDataFig3.Rmd) | [HTML](Rscripts/ExtendedDataFig3.html) |
| Extended Data Fig. 4: BRAF/TERT groups | [R](Rscripts/ExtendedDataFig4.R) / [Rmd](Rscripts/ExtendedDataFig4.Rmd) | [HTML](Rscripts/ExtendedDataFig4.html) |
| Supplementary Fig. 1: SCNA examples | [R](Rscripts/SupplementaryFig1.R) / [Rmd](Rscripts/SupplementaryFig1.Rmd) | [HTML](Rscripts/SupplementaryFig1.html) |
| Supplementary Fig. 2: Minimum CCF | [R](Rscripts/SupplementaryFig2.R) / [Rmd](Rscripts/SupplementaryFig2.Rmd) | [HTML](Rscripts/SupplementaryFig2.html) |
| Supplementary Fig. 3: Reused MATH panels | [R](Rscripts/SupplementaryFig3.R) / [Rmd](Rscripts/SupplementaryFig3.Rmd) | [HTML](Rscripts/SupplementaryFig3.html) |

### Software dependencies

Validation used R 4.5.1. Install the packages in
[software_versions.tsv](provenance/software_versions.tsv), plus Pandoc for HTML.
Roboto Condensed Regular, Bold, Italic and Bold Italic must be installed; the runtime
checks and registers them with showtext, and stops if they are unavailable.
PDF text is outlined through showtext; PNG exports use 300 dpi. Canvas sizes and
font sizes are explicit. No personal `.Rprofile` is sourced.

### Instructions for reproducibility

Run from `TCGA-PTC/`. Obtain the private inputs listed in
[data/README.md](data/README.md) and configure the input/output roots as described
in [config/README.md](config/README.md).

```sh
# Data-free public documentation build
Rscript Rscripts/render_all.R

# Run all figures, or pass a comma-separated subset as the third argument
Rscript Rscripts/run_all.R "$PTC_INPUT_DIR" "$PTC_OUTPUT_DIR"
Rscript Rscripts/run_all.R "$PTC_INPUT_DIR" "$PTC_OUTPUT_DIR" Fig4,Fig6

# Optional local review HTML with private previews; output must be separate
Rscript Rscripts/render_all.R "$PTC_REVIEW_DIR" "$PTC_OUTPUT_DIR"

# Independent inventory, relative-link and publication-safety checks
python3 tools/validate_package.py
```

Each figure runs in a fresh R process. Inspect `*_execution.tsv`, `*_run.log` and
`*_sessionInfo.txt` under your output directory. A successful process exit does
not mean that an unresolved panel was reproduced. Full execution objects and
plots remain local and are excluded by `.gitignore`.

### Validation, gaps and review

The audit records **30 reproduced, 9 partially identified and 7 unresolved entries**.
See [VALIDATION.md](VALIDATION.md) for execution, comparison and QA evidence.

- [Master provenance and status](provenance/figure_provenance.tsv)
- [All unresolved or non-reproduced panels, with reasons](provenance/unresolved_figures.tsv)
- [Per-panel validation](provenance/validation.tsv)
- [Completeness counts](provenance/completeness_summary.json)
- [Scientific and source-version issues requiring review](provenance/methodology_review.md)
- [Keynote map](provenance/keynote_panel_map.tsv) and [embedded/original hashes](provenance/embedded_vs_local_files.tsv)
- [Exact export-line evidence](provenance/source_code_evidence.tsv) and [historical source hashes](provenance/source_catalog.tsv)

`reproduced` means executed and compared with the corresponding embedded plot
for data/analysis/visual correspondence; styling may differ. Native artwork is
explicitly described as a recreation. Related historical variants are not promoted
to final-plot reproduction. Source scripts are named as logical archived filenames
in the provenance table; they are not dependencies located outside this package.

The package preserves recovered statistical methods. In particular, Holm-adjusted
variables labelled FDR, raw-P plot variants, source sample restrictions and model
conventions are documented instead of silently altered. Reproducing a historical
plot does not resolve those scientific/manuscript discrepancies.
