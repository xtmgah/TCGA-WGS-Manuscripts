# Rscript Rscripts/render_all.R [OUTPUT_DIRECTORY] [PRIVATE_FIGURE_DIRECTORY]
# With no arguments, build documentation without plot images beside the Rmd files.
args <- commandArgs(TRUE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),'..'),mustWork=TRUE)
source(file.path(root,'functions','ptc_runtime.R'))
include_figures <- length(args)>1
output <- if(length(args))args[1] else file.path(root,'Rscripts')
if(include_figures) {
  figure_dir <- ptc_external_directory(args[2],'Private figure directory')
  output <- ptc_external_directory(output,'Private HTML output directory',must_exist=FALSE,create=TRUE)
} else {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  output <- normalizePath(output,mustWork=TRUE)
  figure_dir <- ''
}
if(!requireNamespace('rmarkdown',quietly=TRUE)) stop('Install rmarkdown and knitr.')
if(!rmarkdown::pandoc_available()) stop('Pandoc is required. Install Pandoc or set RSTUDIO_PANDOC to its directory.')
files <- list.files(file.path(root,'Rscripts'),pattern='[.]Rmd$',full.names=TRUE)
for(f in files) {
  rmarkdown::render(f,output_dir=output,knit_root_dir=root,
    params=list(include_figures=include_figures,figure_dir=figure_dir),
    envir=new.env(parent=globalenv()),quiet=TRUE)
  cat('Rendered ',basename(f),'\n',sep='')
}
cat('HTML generation finished.\n')
