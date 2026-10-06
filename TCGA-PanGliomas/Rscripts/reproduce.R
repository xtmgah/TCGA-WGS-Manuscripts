#!/usr/bin/env Rscript
args0 <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '',args0[grepl('^--file=',args0)][1])
repo <- normalizePath(file.path(dirname(script),'..'),mustWork=TRUE)
source(file.path(repo,'functions/R/entrypoint.R'))
run_reproduction(repo,commandArgs(trailingOnly=TRUE))
