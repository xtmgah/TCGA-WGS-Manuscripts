# Original Battenberg LogR/BAF/copy-number/CCF sample tracks.
# Segmented copy-number summary objects do not replace probe-level LogR/BAF.
run_SupplementaryFig1 <- function(input_dir,output_dir) {
  setNames(lapply(c('a','b','c'),function(p)ptc_blocked(paste0('Supplementary Fig. 1',p),
    'Exact SCNA PDF export code and per-sample LogR/BAF tracks absent. Requires the original Battenberg plotting workflow, version, probe tracks and segment/purity outputs.','unresolved')),
    paste0('Supplementary Fig. 1',c('a','b','c')))
}
