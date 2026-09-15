# Final Fig.3 has five recovered filenames, but none has a generating export
# call in the 27 supplied source scripts. Do not substitute older MATH plots.
# The shared ptc_math() preparation is traceable to TP_Comparision.R:2363–2387.
prepare_Fig3_math <- function(input_dir) {
  e <- ptc_evolution_inputs(input_dir,math=TRUE)
  ptc_math(e$wgsdata,e$DP_info_data)
}
run_Fig3 <- function(input_dir,output_dir) {
  reasons <- c(
    a='tp_tumor_hetrogenerity_math_ccf2.pdf: related MATH_CCF calculation found; final ridge/box layout generator absent.',
    b='tp_tumor_hetrogenerity_math_ccf_drivers_all.pdf: cohort-specific predecessor found; final pooled model and plotting block absent.',
    c='tp_tumor_hetrogenerity_math_ccf_drivers2_all.pdf: cohort-specific predecessor found; final pooled comparison generator absent.',
    d='NRPCC_threshold.pdf: exact output and NRPCC threshold calculation absent; coverage/purity objects alone do not define the calculation.',
    e='drivers_number_subclones.pdf: cluster counts available; driver-specific final plot/model generator absent.'
  )
  setNames(lapply(names(reasons),function(p)ptc_blocked(paste0('Fig. 3',p),reasons[[p]],if(p=='d')'unresolved' else 'source code partially identified')),paste0('Fig. 3',names(reasons)))
}
