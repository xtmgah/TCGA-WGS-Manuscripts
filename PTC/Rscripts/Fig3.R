# Figure 3: CCF-based heterogeneity preparation.
# Plotting implementations for these panels are not included in this release.
prepare_Fig3_math <- function(input_dir) {
  e <- ptc_evolution_inputs(input_dir,math=TRUE)
  ptc_math(e$wgsdata,e$DP_info_data)
}
run_Fig3 <- function(input_dir,output_dir) {
  reasons <- c(
    a='The ridgeline plotting implementation is not included; CCF-based heterogeneity preparation is available.',
    b='The pooled alteration-association model and plot are not included.',
    c='The pooled alteration-group comparison plot is not included.',
    d='The NRPCC threshold calculation and plot are not included. Coverage and purity values alone do not define the calculation.',
    e='The driver-specific subclone-count plotting implementation is not included.'
  )
  setNames(lapply(names(reasons),function(p)ptc_blocked(paste0('Fig. 3',p),reasons[[p]],'unavailable')),paste0('Fig. 3',names(reasons)))
}
