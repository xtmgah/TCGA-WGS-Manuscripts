# Rscript Rscripts/run_all.R INPUT_DIRECTORY OUTPUT_DIRECTORY [Fig1,Fig2,...]
args <- commandArgs(TRUE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),'..'),mustWork=TRUE)
source(file.path(root,'functions','ptc_runtime.R'))
input <- if(length(args)>=1)args[1] else Sys.getenv('PTC_INPUT_DIR')
output <- if(length(args)>=2)args[2] else Sys.getenv('PTC_OUTPUT_DIR')
if(!nzchar(input) || !nzchar(output))
  stop('Provide input and output directories, or set PTC_INPUT_DIR and PTC_OUTPUT_DIR. Both must be outside the repository.')
input <- ptc_external_directory(input,'Input directory')
output <- ptc_external_directory(output,'Output directory',must_exist=FALSE)
reference <- Sys.getenv('PTC_REFERENCE_DIR')
if(nzchar(reference)) {
  reference <- ptc_external_directory(reference,'Reference directory')
  Sys.setenv(PTC_REFERENCE_DIR=reference)
}
for(path in c(input,if(nzchar(reference)) reference)) {
  if(identical(output,path) || startsWith(output,paste0(path,'/')) || startsWith(path,paste0(output,'/')))
    stop('Output directory must not overlap the input or reference directory.')
}
output <- ptc_external_directory(output,'Output directory',must_exist=FALSE,create=TRUE)
figures <- c(paste0('Fig',1:6),paste0('ExtendedDataFig',1:4),paste0('SupplementaryFig',1:3))
if(length(args)>=3) {
  selected <- strsplit(args[3],',',fixed=TRUE)[[1]]
  if(!all(selected %in% figures)) stop('Unknown figure selection: ',paste(setdiff(selected,figures),collapse=', '))
  figures <- selected
}
if(!requireNamespace('ragg',quietly=TRUE)) stop('Install the packages listed in config/software_versions.tsv.')
codes <- integer()
for(fig in figures) {
  log <- file.path(output,paste0(fig,'_run.log'))
  codes <- c(codes,system2(file.path(R.home('bin'),'Rscript'),
    args=c('--vanilla',shQuote(file.path(root,'tools','run_one.R')),shQuote(root),fig,shQuote(input),shQuote(output)),stdout=log,stderr=log))
  cat(fig,': process exit ',tail(codes,1),' — ',basename(log),'\n',sep='')
}
cat('Finished. Check *_execution.tsv for generated plots, unavailable panels and input errors.\n')
if(any(codes!=0)) quit(status=1)
