# Rscript Rscripts/run_all.R INPUT_DIRECTORY OUTPUT_DIRECTORY [Fig1,Fig2,...]
args<-commandArgs(TRUE)
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),'..'),mustWork=TRUE)
input<-if(length(args)>=1)args[1] else Sys.getenv('PTC_INPUT_DIR',file.path(root,'data'))
output<-if(length(args)>=2)args[2] else Sys.getenv('PTC_OUTPUT_DIR',file.path(root,'Rscripts','Figures'))
input<-normalizePath(input,mustWork=FALSE)
dir.create(output,recursive=TRUE,showWarnings=FALSE);output<-normalizePath(output)
figures<-c(paste0('Fig',1:6),paste0('ExtendedDataFig',1:4),paste0('SupplementaryFig',1:3))
if(length(args)>=3){selected<-strsplit(args[3],',',fixed=TRUE)[[1]];stopifnot(all(selected %in% figures));figures<-selected}
if(!requireNamespace('ragg',quietly=TRUE))stop('Install dependencies listed in provenance/software_versions.tsv.')
codes<-integer()
for(fig in figures){
 log<-file.path(output,paste0(fig,'_run.log'))
 codes<-c(codes,system2(file.path(R.home('bin'),'Rscript'),
   args=c('--vanilla',shQuote(file.path(root,'tools','run_one.R')),shQuote(root),fig,shQuote(input),shQuote(output)),stdout=log,stderr=log))
 cat(fig,': process exit ',tail(codes,1),' — ',basename(log),'\n',sep='')
}
cat('Completed. A process exit of zero does not mean every panel reproduced; inspect *_execution.tsv.\n')
if(any(codes!=0))quit(status=1)
