#!/usr/bin/env Rscript
script <- sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
repo <- normalizePath(file.path(dirname(script),'..'))
source(file.path(repo,'functions/R/entrypoint.R'))
run_reproduction(repo,c('--figures','5',commandArgs(trailingOnly=TRUE)))
