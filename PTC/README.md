# Driver-defined evolutionary routes in papillary thyroid cancer

Analysis code for whole-genome sequencing of papillary thyroid carcinomas from
TCGA-THCA and the Chornobyl exposed and unexposed cohorts.

The analyses cover somatic drivers, copy-number alterations, intratumor
heterogeneity, mutational signatures, and tumor evolutionary timing.
Figure-specific R scripts and R Markdown reports are provided below.

## Figures

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
| Supplementary Fig. 3: MATH associations | [R](Rscripts/SupplementaryFig3.R) / [Rmd](Rscripts/SupplementaryFig3.Rmd) | [HTML](Rscripts/SupplementaryFig3.html) |

Download the HTML files to view the analysis documentation locally. The reports
display the corresponding R code. Generating figures requires the study inputs
described in [data/README.md](data/README.md).

## Software

The code was run with R 4.5.1. Package versions are listed in
[config/software_versions.tsv](config/software_versions.tsv). Rendering HTML
requires rmarkdown, knitr, and Pandoc. Plotting requires Roboto Condensed with
Regular, Bold, Italic, and Bold Italic styles.

## Run the analyses

Run commands from the `PTC/` directory. Set the input, reference, and output
locations to directories outside the repository; see
[config/README.md](config/README.md).

```sh
export PTC_INPUT_DIR="/path/to/private/analysis_inputs"
export PTC_REFERENCE_DIR="/path/to/GRCh38/reference_files"
export PTC_OUTPUT_DIR="/path/to/local/ptc_figures"

# Run all available analyses
Rscript Rscripts/run_all.R "$PTC_INPUT_DIR" "$PTC_OUTPUT_DIR"

# Run selected figures
Rscript Rscripts/run_all.R "$PTC_INPUT_DIR" "$PTC_OUTPUT_DIR" Fig4,Fig6

# Render documentation without loading study data
Rscript Rscripts/render_all.R

# Optional reports with locally generated figures
Rscript Rscripts/render_all.R "/path/to/local/ptc_reports" "$PTC_OUTPUT_DIR"
```

Each figure runs in a separate R process. Outputs include PDF/PNG figures, execution
logs, and local R results. The figure pages describe the available analyses and
any additional inputs or implementations required for individual panels.

## Data

[data/wgs_data.RData](data/wgs_data.RData) contains a three-row cohort summary
(`wgs_data`: cohort name and tumor count) and cohort plotting colors
(`study_color2`). It contains no individual sample records or genomic measurements.

Private analysis inputs, genomic reference files, and generated results are kept
outside the repository. The [input specification](data/input_manifest.tsv) lists
expected filenames, object names, and column schemas without data values.
All genomic coordinates must use GRCh38. Oncoplot colors are defined in
[config/oncoplot_colors.csv](config/oncoplot_colors.csv).
