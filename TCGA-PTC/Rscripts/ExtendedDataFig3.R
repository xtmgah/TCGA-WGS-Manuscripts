# Original sample-specific DPClust exports were embedded in Keynote.
# No export call, plotting function or complete upstream plotting input was
# recovered in the supplied 27 scripts. The DP_info_data snapshot contains
# assigned clusters but does not establish the original density fit/renderer.
run_ExtendedDataFig3 <- function(input_dir,output_dir) {
  setNames(lapply(c('a','b'),function(p)ptc_blocked(paste0('Extended Data Fig. 3',p),
    'DPClust original export generator and complete per-sample plotting inputs are absent. Recover the historical pipeline/version and plotting invocation before replay.','unresolved')),
    paste0('Extended Data Fig. 3',c('a','b')))
}
