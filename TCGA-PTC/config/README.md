# Input and output configuration

Run commands from `TCGA-PTC/`. Pass input/output directories to `Rscripts/run_all.R`,
or set `PTC_INPUT_DIR` and `PTC_OUTPUT_DIR` in your shell. The defaults are `data/`
and `Rscripts/Figures/`. Private inputs and generated figures are ignored by git.

`PTC_INPUT_DIR` is the root containing the documented RData files and the
`thyroid_evolution_analysis/` subdirectory. It is read-only to this workflow.
The code never submits cluster jobs or writes into the input directory.

All coordinate resources must be GRCh38. The supplied reference subset is
validated against GRCh38 chromosome lengths; no coordinate conversion occurs.
