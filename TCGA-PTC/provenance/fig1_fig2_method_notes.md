# Fig. 1 and Fig. 2: reproduction scope and retained decisions

The cleaned modules are [Fig1.R](../Rscripts/Fig1.R) and [Fig2.R](../Rscripts/Fig2.R). They read author-supplied input snapshots and write only to the caller's output directory. The source analysis directory is never used as an output directory by the package. Historical input assembly and manual curation remain separate from regenerating plots from these snapshots.

## Execution

From the `TCGA-PTC` directory, source `functions/ptc_runtime.R` and the required figure module, then call `run_Fig1(input_dir, output_dir)` or `run_Fig2(input_dir, output_dir)`. Each panel returns its own result or failure, so one unavailable dependency does not prevent other panels from being attempted. An optional `panels` argument selects lower-case panel letters for a targeted rerun.

Both modules use Roboto Condensed. Plot dimensions are enlarged where needed to preserve readable text. The plotting helpers use a `ragg` temporary device for layout measurement, avoiding fallback font metrics from the default PostScript device. PDF text produced by `showtext` can be stored as vector outlines rather than selectable text.

## What can be reproduced

| Manuscript panel | Cleaned output basename | Scope |
|---|---|---|
| Fig. 1a | `Fig1a_related_age`, `Fig1a_related_depth`, `Fig1a_related_purity` | Related historical components; the final `study2_overview.pdf` assembly code was not recovered. |
| Fig. 1b | `Fig1b` | IntOGen driver frequencies from saved curated events and driver annotations; figure legend styling is reconstructed from known cohort/role labels. |
| Fig. 1c | `Fig1c_related_…` | Recovered cohort-specific 5-Mb breakpoint densities; the final combined layout and hotspot annotation code were not recovered. |
| Fig. 1d–f | `Fig1d`–`Fig1f` | Recovered Circos plot method with explicit cohort filters. |
| Fig. 1g–i | `Fig1g`–`Fig1i` | Recovered oncoprint method, exact saved curated events, shared gene order and explicit cohort filters. |
| Fig. 2a | `Fig2a_state_reconstruction` | Explicit clone/subclone state reconstruction; final generator remains incompletely recovered. |
| Fig. 2b | `Fig2b` | Historical chr22 copy-number heatmap. |
| Fig. 2c | `Fig2c` | Historical subclonal-fraction distribution; segment-level observations are retained. |
| Fig. 2d | `Fig2d_related_raw_p` | Historical raw-P variant, which differs from the final embedded adjusted-P plot. |

Every basename has a PDF and a 300-dpi PNG. Related outputs and reconstructions must not be represented as an exact recovery of missing final export code.

## Fig. 1 decisions

- **Cohort overview:** `Study2` defines TCGA-THCA, Chornobyl-unexposed and Chornobyl-exposed groups. Historical comparison blocks use `MEDIAN_COVERAGE`, `BB_Purity` and `age_at_diagnosis`. The historical depth export contains tumor and normal facets; this differs from the final overview's tumor-only component. The reusable comparison helper preserves pairwise Wilcoxon tests with BH adjustment and seed 1 for quasirandom point placement.
- **IntOGen:** the original mounted `drivers.tsv` is unavailable in the source snapshot, but `intogene_drivers.RData` contains its imported `intogene` object. Frequency calculations preserve distinct tumor–gene events, the covariance-snapshot inclusion filter, the overall cohort denominator and cohort-specific denominators. `wgs_covdata` is loaded explicitly because the historical script used this object without loading it. A vector comparison of the 57 cohort data points (19 genes × 3 cohorts) confirmed that the recovered curves correspond to the embedded final PDF after accounting for different canvas dimensions.
- **Breakpoint density:** both breakends are counted in 5,000,000-bp intervals, normalized by cohort sample count. The historical 0.5 upper y-axis limit is retained; points beyond that limit can trigger warnings. The missing combined-panel exporter and manually placed hotspot labels are not invented.
- **Circos:** the source contains three commented cohort filters and matching commented export filenames. The cleaned module selects these pairs explicitly. Links use chromosomes 1–22 and X, 1,000-bp intervals and the original SV-type colors. The ideogram contains the 24 primary hg38 chromosomes, including Y, consistent with the historical hg38 ideogram. Local cytobands replace an implicit external lookup.
- **Oncoprints:** the manually curated `data_top0` and saved `data_tmb0` objects are used directly. The original three `if(TRUE)` cohort blocks overwrite one another; the cleaned code iterates their intended cohort settings. The original gene-order logic, split tiles, alteration aggregation and sample ordering are retained. The historical helper may warn about tumor–gene combinations with multiple alteration records; this behavior was not silently reduced to one alteration. Only the three helper functions actually used by these panels are included.

## Fig. 2 decisions and unresolved differences

- Clonal versus subclonal states are assigned from the larger versus smaller `frac*_A` fraction, using the exact historical rules, including tied fractions and missing values. Relative copy number uses a baseline of 2 or 4 according to WGD status and retains the original ±4 caps and homozygous-deletion sentinel.
- The genome-wide frequency method retains 1-Mb bins, CNTools mean values, centromere treatment and the historical exclusion of chromosome 21 values above 5% in absolute frequency. These are documented analysis decisions, not newly introduced quality filters.
- The saved genome-wide script sets both `p_freq_clone` and `p_freq_subclone` to the same plot. A commented state switch demonstrates how the historical operator could generate the two inputs separately. The cleaned function exposes that switch as an argument, but the result is labelled a **state reconstruction** because the final execution sequence and Illustrator edits are not completely recovered.
- For the chr22 heatmap and deletion classifications, the longest chr22 segment defines the clonal and subclonal loss sets. No new arm-overlap threshold is introduced. The historical heatmap retains neutral samples as white rows. The available snapshots yield 100 clonal-loss and 63 subclonal-loss tumors.
- The subclone-fraction plot retains one observation per qualifying segment, without silent deduplication: the available snapshot contains **69 observations from 63 unique tumors**. The manuscript legend refers to individual tumors; this difference requires author review even though the historical plotted distribution is reproduced.
- The historical association code uses two-sided Fisher tests and `p.adjust()` with its default **Holm** method, despite naming the resulting column `FDR`. The cleaned code writes `method='holm'` explicitly. It retains the source's **raw-P y-axis** and labels based on adjusted significance. The final embedded PDF uses an adjusted-P axis and a narrower displayed event set. Its exact generating export code was not recovered, so the historical result remains a clearly named related output.

## Supporting references and helper provenance

[reference_checks.json](../data/refs/reference_checks.json) records the bundled reference validation and source hashes. The cytoband maximal endpoints for chromosomes 1–22, X and Y match the local GRCh38 reference FAI lengths. The primary FAI subset also retains chrM for the historical genomic-burden denominator used elsewhere in the package. These are small public reference resources; they contain no cohort observations.

[fig1_fig2_helper_sources.json](fig1_fig2_helper_sources.json) records source hashes and extraction ranges for the reused functions. Cleaning removes debugging calls and inactive unrelated code, makes selected paths and cohorts explicit, and increases typography. The unused data-frame p-value-override branch in the comparison helper had a misspelled variable name; that spelling is corrected. No p-value override is used by Fig. 1 or Fig. 2.
