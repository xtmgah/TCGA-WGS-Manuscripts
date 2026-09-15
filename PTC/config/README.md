# Configuration

Run commands from `PTC/`. Configure three directories outside the Git repository:

| Variable | Contents |
| --- | --- |
| `PTC_INPUT_DIR` | Private analysis inputs listed in [input_manifest.tsv](../data/input_manifest.tsv). |
| `PTC_REFERENCE_DIR` | GRCh38 reference files described in [data/README.md](../data/README.md). |
| `PTC_OUTPUT_DIR` | Locally generated figures, execution logs, and R results. |

Input and output directories must be separate. Keep optional HTML reports that
include generated plots outside the repository as well.

[software_versions.tsv](software_versions.tsv) lists the R package versions used
for these scripts. [oncoplot_colors.csv](oncoplot_colors.csv) defines plotting
colors and contains no participant measurements.
