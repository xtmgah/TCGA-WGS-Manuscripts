args <- commandArgs(trailingOnly=TRUE)
script <- sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
repo <- normalizePath(file.path(dirname(script),'../..'))
out <- if(length(args))normalizePath(args[1])else file.path(repo,'outputs')
manifest <- jsonlite::fromJSON(file.path(out,'run_manifest.json'))
for(n in manifest$figures)rmarkdown::render(file.path(repo,'Rscripts',paste0('Fig',n,'.Rmd')),
 params=list(rebuild=FALSE,output_dir=out),envir=new.env(parent=globalenv()),quiet=TRUE)
rmarkdown::render(file.path(repo,'Rscripts/data_format.Rmd'),envir=new.env(parent=globalenv()),quiet=TRUE)
