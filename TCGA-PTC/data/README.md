# Required data and reference build

No participant-level datasets are included. [input_manifest.tsv](input_manifest.tsv)
lists expected filenames, required R objects, schemas, source-file hashes and panels.
Hashes identify the locally audited snapshot; they do not grant access to it.

## Access and preparation

The manuscript's data-availability paragraph reports TCGA-THCA access under
**phs000178.v11.p8** and Chornobyl access under **phs001134**. Obtain the necessary
dbGaP authorization and follow the study's approved distribution instructions.
These accessions are transcribed from the supplied manuscript; this package has
not independently verified its claim that both downloads are available through GDC.
Raw sequencing data are not direct inputs to this figure package.

The package starts from author-prepared, QC-filtered RData snapshots. Request the
matching approved analysis intermediates from the authors, or regenerate them with
the original WGS pipeline and the upstream scripts listed in the master provenance.
Preserve the cohort identifiers, curation rules, filters, reference build and
snapshot versions. A downloaded raw BAM/CRAM alone cannot replace these objects.

Place files under the configured private input root, preserving the documented
subdirectories. Code checks file/object presence before plotting. Private inputs
are read-only; outputs are written to a separate configured directory.

## GRCh38 compatibility

The supplied manuscript (paragraph 48) specifies **GRCh38**. Historical R scripts
reference assembly38 FASTA indexes, hg38 gene annotations and hg38 exclusion ranges.
The small reference subset in `refs/` was checked against GRCh38 chromosome ends.
No alignment, liftover or coordinate conversion is performed here. Inputs derived
from GRCh37 or CHM13 must not be mixed with these resources.

Reference/resource hashes and provenance are listed in
[reference_manifest.tsv](reference_manifest.tsv). The [software versions](../provenance/software_versions.tsv)
record packages present during validation, rather than asserting a fully frozen
historical environment.

## Important input distinctions

- `DP_info_data.RData` contains assigned mutation/cluster CCFs; it does not preserve
  every original DPClust density-fit or sample-plot input.
- `Genome_landscape_manual_final.RData` is a curated event snapshot. Rebuilding it
  requires the original curation, not simply re-running a plotting script.
- `Chronological_timing_short.RData` contains multiple acceleration scenarios.
  Primary latency plots explicitly select `1x`; upstream chronological inference
  is documented separately from figure replay.
- Summary copy-number objects cannot replace probe-level LogR/BAF tracks required
  by the unresolved sample-specific SCNA panels.
- Source plotting inputs may contain extra or excluded samples; each module
  preserves the recovered source filtering and documents denominator differences.
