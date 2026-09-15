# Supplementary Figure 1: Battenberg LogR/BAF, copy-number and CCF tracks.
# Segment summaries cannot replace the required probe-level measurements.
run_SupplementaryFig1 <- function(input_dir,output_dir) {
  setNames(lapply(c('a','b','c'),function(p)ptc_blocked(paste0('Supplementary Fig. 1',p),
    'This release requires the Battenberg sample plotting implementation, software version, probe-level LogR/BAF tracks and segment/purity outputs.','unavailable')),
    paste0('Supplementary Fig. 1',c('a','b','c')))
}
