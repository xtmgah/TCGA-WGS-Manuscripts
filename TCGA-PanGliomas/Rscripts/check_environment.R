#!/usr/bin/env Rscript
script <- sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
repo <- normalizePath(file.path(dirname(script),'..'))
d <- read.delim(file.path(repo,'environment/R-direct-packages.tsv'))
missing <- d$Package[!vapply(d$Package,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
if(!rmarkdown::pandoc_available())stop('Pandoc is required for HTML reports')
if(!any(systemfonts::system_fonts()$family=='Roboto Condensed'))stop('Install the bundled Roboto Condensed fonts for Cairo/R')
if(Sys.info()[['sysname']]=='Darwin') {
 selected <- systemfonts::match_fonts(rep('Roboto Condensed',3),italic=c(FALSE,FALSE,TRUE),weight=c('normal','bold','normal'))$path
 expected <- file.path(repo,'functions/fonts/system',c('RobotoCondensed-VariableFont_wght.ttf','RobotoCondensed-VariableFont_wght.ttf','RobotoCondensed-Italic-VariableFont_wght.ttf'))
 if(!identical(unname(tools::md5sum(selected)),unname(tools::md5sum(expected))))stop('Register the exact bundled variable Roboto Condensed faces in functions/fonts/system/ with macOS Font Book; another font version is selected.')
}
actual <- vapply(d$Package,function(p)as.character(packageVersion(p)),character(1))
changed <- actual!=gsub('-','.',d$Version,fixed=TRUE)
if(any(changed))message('Versions differ from validated snapshot: ',paste(d$Package[changed],collapse=', '))
message('R dependencies, Pandoc, and Roboto Condensed: PASS')
message('Run Rscripts/reproduce.R --check-inputs to validate the data bundle.')
