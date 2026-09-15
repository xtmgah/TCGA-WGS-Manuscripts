# Analysis inputs

The included [wgs_data.RData](wgs_data.RData) is a small cohort-level summary:

| Object | Contents |
| --- | --- |
| `wgs_data` | Three rows containing cohort names and tumor counts: TCGA-THCA (455), Chornobyl-unexposed (70), and Chornobyl-exposed (283). |
| `study_color2` | A named color vector for those three cohorts. |

This file contains no sample identifiers, individual participant records,
mutation calls, copy-number measurements, or sequence data. It does not replace
the private inputs required to run the analyses.

[input_manifest.tsv](input_manifest.tsv) specifies the filenames, R object names,
and column schemas expected by the figure scripts. These are input descriptions
without participant-level values.

## Obtain and prepare inputs

Access to sequencing data is described in the manuscript's Data availability
statement. Use approved study data and the matching author-prepared analysis
intermediates. The figure scripts start from processed tables rather than raw
BAM/CRAM files. Obtain the required intermediates from the study authors under
the applicable data-access conditions.

Keep the inputs outside this repository. Set `PTC_INPUT_DIR` to their root directory
and preserve the filenames and subdirectories in the input specification. The
workflow reads the inputs and writes results to a separate output directory.

## Genome reference

The analyses use **GRCh38**. Provide the following reference annotations in an
external directory and set `PTC_REFERENCE_DIR` to that directory:

| File | Required contents |
| --- | --- |
| `hg38_primary.fai` | GRCh38 FASTA index with chromosome names and lengths; the first 24 records must be chr1 through chr22, chrX, and chrY in that order. |
| `hg38_cytoBand.txt.gz` | GRCh38 cytobands with chromosome, start, end, band name, and stain columns. |

Use matching chromosome names and primary chromosome lengths throughout the input
tables and reference annotations. Do not mix GRCh37 or CHM13 coordinates with GRCh38.
The plotting workflow does not align reads or convert genome coordinates.

## Input details

- `Genome_landscape_manual_final.RData` supplies the curated somatic-event tables.
- `DP_info_data.RData` supplies mutation and cluster cancer-cell fractions; the
  per-sample DPClust density plots require additional plotting inputs.
- `thyroid_evolution_analysis/Chronological_timing_short.RData` supplies timing
  estimates. The principal age/latency analyses select the `1x` acceleration scenario.
- Probe-level LogR/BAF measurements are needed for sample-specific copy-number plots;
  segment summaries cannot replace those measurements.
- Filenames, R object names, cohort labels, and column names must match the input
  specification. Individual scripts apply their analysis-specific filters.
