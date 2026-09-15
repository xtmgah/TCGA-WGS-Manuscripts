# Extended Data Figure 3: sample-specific DPClust examples.
# Complete density-fit inputs and the sample plotting implementation are required.
run_ExtendedDataFig3 <- function(input_dir,output_dir) {
  setNames(lapply(c('a','b'),function(p)ptc_blocked(paste0('Extended Data Fig. 3',p),
    'This release requires the DPClust sample plotting implementation, software version and complete density-fit inputs.','unavailable')),
    paste0('Extended Data Fig. 3',c('a','b')))
}
