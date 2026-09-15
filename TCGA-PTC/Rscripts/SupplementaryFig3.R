# The final Keynote reuses identical imported PDFs in ED2 and Supplementary Fig3.
# Shared implementation prevents duplicated code from drifting.
run_SupplementaryFig3 <- function(input_dir,output_dir) {
  run_ExtendedDataFig2(input_dir,output_dir,prefix='SupplementaryFig3',figure_label='Supplementary Fig. 3')
}
