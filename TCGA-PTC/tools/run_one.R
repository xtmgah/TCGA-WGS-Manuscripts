# Invoked by Rscripts/run_all.R. One R process per figure prevents state leakage.
a <- commandArgs(TRUE)
if(length(a)!=4L)stop('Arguments: PTC_ROOT FIGURE INPUT_DIR OUTPUT_DIR')
root<-normalizePath(a[1],mustWork=TRUE);figure<-a[2]
input<-normalizePath(a[3],mustWork=FALSE);output<-normalizePath(a[4],mustWork=FALSE)
setwd(root);dir.create(output,recursive=TRUE,showWarnings=FALSE)
source('functions/ptc_runtime.R')
for(f in setdiff(list.files('functions',pattern='[.]R$',full.names=TRUE),'functions/ptc_runtime.R'))source(f)
for(f in list.files('Rscripts',pattern='^(Fig[0-9]+|Fig4_schematic|ExtendedDataFig[0-9]+|SupplementaryFig[0-9]+)[.]R$',full.names=TRUE))source(f)
result<-tryCatch(get(paste0('run_',figure))(input,output),error=function(e)list(run_error=conditionMessage(e)))
# Full execution objects are private local validation artifacts, never committed.
saveRDS(result,file.path(output,paste0(figure,'_execution.rds')))
extract_paths<-function(x){
 if(is.character(x))return(x[grepl('[.](pdf|png)$',x)])
 if(is.list(x))return(unlist(lapply(x,extract_paths),use.names=FALSE))
 character()
}
rows<-lapply(names(result),function(n){
 z<-result[[n]]
 if(!is.list(z))return(data.frame(panel_id=n,execution_status='run_error',outputs='',notes=as.character(z)))
 s<-z$status
 if(is.null(s))s<-z$reproducibility_status
 data.frame(panel_id=if(is.null(z$panel_id))n else paste(z$panel_id,collapse=';'),
  execution_status=if(is.null(s))'generated_unverified' else s,
  outputs=paste(unique(basename(extract_paths(z))),collapse=';'),
  notes=paste(c(z$reason,z$notes,z$warnings),collapse='; '))
})
write.table(do.call(rbind,rows),file.path(output,paste0(figure,'_execution.tsv')),sep='\t',quote=TRUE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,paste0(figure,'_sessionInfo.txt')))
cat(figure,'finished; inspect per-panel execution status.\n')
