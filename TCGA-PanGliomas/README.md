# TCGA Pan-Glioma main figures

This repository contains the R and Python code used to reproduce **main Figures 1–5** of the TCGA Pan-Glioma project. It includes:

- R scripts for each main figure, supported by shared R and Python functions.
- Processed figure-support data and documentation of the required inputs.
- R Markdown notebooks and rendered HTML reports for each figure.

**To view a figure without running the code, download its HTML report and open it locally in a web browser.** Full regeneration requires additional inputs requested from the authors; these are not included in Git.

## Software dependencies

- **R 4.5.1** and **Python 3.12.13**, with packages listed in the [environment recipe](environment/clean-macos.yml) and [Python requirements](requirements.txt).
- **Pandoc 3.12** for rendering HTML reports.
- The supplied **Roboto Condensed fonts**, including the variable regular and italic faces registered with macOS Font Book for R/Cairo.

Follow the [tested installation instructions](docs/ENVIRONMENT_VALIDATION.md#repeat-the-tested-setup-on-macos-arm64) for the complete setup, including `hrbrthemes` and fonts. Validation used a freshly installed R/Python environment on **macOS with Apple silicon**; other platforms have not yet been tested.

**Hardware requirements:** The figure workflow runs on a standard desktop or laptop. The validated build took approximately two minutes, excluding dependency installation and HTML reports; runtime will vary by machine.

## Figures included

**Fig. 1: Inferred genomic event ordering across TCGA gliomas.**
Pooled and group-specific event orderings, glioma subtype composition, chromosome-gain timing, and inferred tumor chronology.

[R script](Rscripts/Fig1.R) · [R Markdown](Rscripts/Fig1.Rmd) · [HTML report](Rscripts/Fig1.html)

**Fig. 2: Inferred event ordering and survival in astrocytoma.**
Event-ordering groups, survival analyses, TP53/ATRX/CIC/TERT alterations, TP53 dosage, tumor chronology, and differential expression.

[R script](Rscripts/Fig2.R) · [R Markdown](Rscripts/Fig2.Rmd) · [HTML report](Rscripts/Fig2.html)

**Fig. 3: Selective structural remodeling in astrocytomas.**
Chromothripsis, structural-variant burden, telomere and TERT measurements, chromosome-gain timing, and genomic-feature comparisons between astrocytoma groups.

[R script](Rscripts/Fig3.R) · [R Markdown](Rscripts/Fig3.Rmd) · [HTML report](Rscripts/Fig3.html)

**Fig. 4: Inferred event ordering and molecular phenotypes in GBM.**
Event-ordering groups in glioblastoma (GBM), survival, chromosome-gain timing, latency, intratumor heterogeneity, genome alteration, and ploidy.

[R script](Rscripts/Fig4.R) · [R Markdown](Rscripts/Fig4.Rmd) · [HTML report](Rscripts/Fig4.html)

**Fig. 5: Amplicon architecture, EGFR amplification and transcriptional states in GBM.**
Amplicon classes, extrachromosomal DNA (ecDNA) oncogene burden and co-carriage, EGFR copy number/expression/fusions, gene-set enrichment, and transcriptional states.

[R script](Rscripts/Fig5.R) · [R Markdown](Rscripts/Fig5.Rmd) · [HTML report](Rscripts/Fig5.html)

See the [manuscript figure captions](reference/main_figure_captions.txt) for panel descriptions, cohorts, and statistical annotations.

## Instructions for reproducibility

Keep the repository folders together and run the commands below **from the folder containing this README**:

```text
TCGA-PanGliomas/
├── README.md
├── Rscripts/          # Figure scripts, notebooks, and HTML reports
├── functions/         # Shared code and fonts; loaded by the scripts
├── data/
│   ├── processed/     # Inputs included in Git
│   └── external/      # Author-request inputs; excluded from Git
├── environment/      # Dependency recipe and validated versions
├── docs/             # Setup, methods, and validation details
├── reference/        # Manuscript figures/captions for comparison
├── tests/            # Numerical and visual verification
└── outputs/          # Created when the figures are regenerated
```

### 1. Set up the software

Complete the [installation instructions](docs/ENVIRONMENT_VALIDATION.md#repeat-the-tested-setup-on-macos-arm64), including environment activation and fonts, before running the figure commands. They do not install dependencies automatically. Check the R packages, Pandoc, and fonts with:

```sh
Rscript Rscripts/check_environment.R
```

### 2. Obtain the required data

Git is the project's only public data location. It includes **70 aggregate figure-support and style files**. Full reproduction also requires **35 author-request files**, approximately 0.763 MB compressed. Contact the corresponding authors listed in the manuscript, following the [data-request instructions](data/external/README.md).

The authors retain these inputs privately on Box and/or Biowulf; neither provides public access for this project. Requests are handled under the applicable source and institutional permissions. No raw sequencing data are included or needed for this figure workflow.

After receiving the matching `panglioma-request-data.tar.gz` archive, extract it as follows (replace the example archive path):

```sh
mkdir -p data/external
tar -xzf /path/to/panglioma-request-data.tar.gz -C data/external
Rscript Rscripts/reproduce.R --check-inputs
```

The request-only data and generated outputs are excluded by `.gitignore`. **Do not upload the private archive to Git.** If the inputs are already installed and pass the check, proceed directly to the build.

### 3. Reproduce the figures

```sh
Rscript Rscripts/reproduce.R --figures 1-5 --verify --render-reports
```

The build produces:

- Individual PDFs and PNG previews in `outputs/figures/`.
- A combined five-page PDF at `outputs/figures/main_figures.pdf`.
- Calculated tables, logs, and verification results under `outputs/`.
- Updated HTML reports in `Rscripts/`.

To build one figure, run `Rscript Rscripts/Fig3.R`, for example. To build a selection, use `Rscript Rscripts/reproduce.R --figures 2,4 --verify --render-reports`. Single-figure builds still require the complete input set because the figures share calculations. Knitting a figure's `.Rmd` also rebuilds that figure by default.

## Data format

Download and open [data_format.html](Rscripts/data_format.html) for an overview of the inputs, or read the [data README](data/README.md). The [input manifest](data/manifest.tsv) lists every required file, its columns, checksums, and distribution route.

## WGS analysis and reproduction scope

This package starts from processed whole-genome sequencing (WGS) and expression-analysis results. It redraws all **49 main-figure panels** and recomputes the documented statistical comparisons. Upstream event-ordering, timing, DESeq2, and GSEA results are retained as saved inputs; the package does not rerun those analyses from raw sequencing data.

All five figures passed numerical checks and panel-level review against the manuscript, with small documented rendering differences. See the [validation record](docs/VALIDATION.md), [reproduction scope](docs/REPRODUCTION_SCOPE.md), and [statistical methods](docs/STATISTICAL_RECOMPUTATION.md) for the evidence and limits of reproduction. Data redistribution decisions are documented in the [distribution review](docs/DATA_DISTRIBUTION.md).
