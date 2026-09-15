# One R process per figure keeps package and plotting state independent.
a <- commandArgs(TRUE)
if(length(a)!=4L) stop('Arguments: PTC_ROOT FIGURE INPUT_DIR OUTPUT_DIR')
root <- normalizePath(a[1],mustWork=TRUE)
figure <- a[2]
setwd(root)
source('functions/ptc_runtime.R')
input <- ptc_external_directory(a[3],'Input directory')
output <- ptc_external_directory(a[4],'Output directory',must_exist=FALSE)
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
for(f in setdiff(list.files('functions',pattern='[.]R$',full.names=TRUE),'functions/ptc_runtime.R')) source(f)
for(f in list.files('Rscripts',pattern='^(Fig[0-9]+|Fig4_schematic|ExtendedDataFig[0-9]+|SupplementaryFig[0-9]+)[.]R$',full.names=TRUE)) source(f)
result <- tryCatch(get(paste0('run_',figure))(input,output),error=function(e)list(run_error=conditionMessage(e)))
# Results and intermediate objects are saved only to the external output directory.
saveRDS(result,file.path(output,paste0(figure,'_execution.rds')))
extract_paths <- function(x) {
  if(is.character(x)) return(x[grepl('[.](pdf|png)$',x)])
  if(is.list(x)) return(unlist(lapply(x,extract_paths),use.names=FALSE))
  character()
}
rows <- lapply(names(result),function(n) {
  z <- result[[n]]
  if(!is.list(z)) return(data.frame(panel_id=n,execution_status='run_error',outputs='',notes=as.character(z)))
  s <- z$status
  if(is.null(s)) s <- z$reproducibility_status
  data.frame(panel_id=if(is.null(z$panel_id))n else paste(z$panel_id,collapse=';'),
    execution_status=if(is.null(s))'generated' else s,
    outputs=paste(unique(basename(extract_paths(z))),collapse=';'),
    notes=paste(c(z$reason,z$notes,z$warnings),collapse='; '))
})
write.table(do.call(rbind,rows),file.path(output,paste0(figure,'_execution.tsv')),sep='\t',quote=TRUE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,paste0(figure,'_sessionInfo.txt')))
cat(figure,'finished; inspect the per-panel output status.\n')
