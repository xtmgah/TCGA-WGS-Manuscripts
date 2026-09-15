# Rscript Rscripts/render_all.R [OUTPUT_DIRECTORY] [PRIVATE_FIGURE_DIRECTORY]
# Default is a data-free documentation build next to the .Rmd files.
args<-commandArgs(TRUE)
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),'..'),mustWork=TRUE)
if(!requireNamespace('rmarkdown',quietly=TRUE))stop('Install rmarkdown and knitr.')
if(!rmarkdown::pandoc_available())stop('Pandoc is required. Install Pandoc or set RSTUDIO_PANDOC to its directory.')
output<-if(length(args))args[1] else file.path(root,'Rscripts')
dir.create(output,recursive=TRUE,showWarnings=FALSE);output<-normalizePath(output)
figure_dir<-if(length(args)>1)normalizePath(args[2],mustWork=TRUE) else file.path(root,'Rscripts','Figures')
if(length(args)>1 && output==normalizePath(file.path(root,'Rscripts')))
 stop('Private figure previews must be rendered to a separate output directory outside the public Rscripts directory.')
files<-list.files(file.path(root,'Rscripts'),pattern='[.]Rmd$',full.names=TRUE)
for(f in files){
 rmarkdown::render(f,output_dir=output,knit_root_dir=root,
  params=list(include_figures=length(args)>1,figure_dir=figure_dir),
  envir=new.env(parent=globalenv()),quiet=TRUE)
 cat('Rendered ',basename(f),'\n',sep='')
}
cat('Documentation render finished; private analysis execution is a separate step.\n')
