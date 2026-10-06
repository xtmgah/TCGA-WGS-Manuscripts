# Shared entry point for scripts and executable R Markdown reports.
run_reproduction <- function(repo,args=character()) {
 runtime_path<-file.path(repo,'.local/runtime.json')
 runtime<-if(file.exists(runtime_path))jsonlite::fromJSON(runtime_path)else list()
 python<-Sys.getenv('PANGLIOMA_PYTHON',unset=if(!is.null(runtime$python))runtime$python else 'python3')
 rscript<-Sys.getenv('PANGLIOMA_RSCRIPT',unset=if(!is.null(runtime$rscript))runtime$rscript else file.path(R.home('bin'),'Rscript'))
 deps<-file.path(repo,'.local/python_deps')
 old<-Sys.getenv('PYTHONPATH');on.exit(Sys.setenv(PYTHONPATH=old),add=TRUE)
 # An explicitly selected interpreter must use its own installed packages.
 if(!nzchar(Sys.getenv('PANGLIOMA_PYTHON')) && dir.exists(deps))Sys.setenv(PYTHONPATH=paste(c(deps,old[nzchar(old)]),collapse=.Platform$path.sep))
 Sys.setenv(PYTHONDONTWRITEBYTECODE='1')
 cmd<-c(file.path(repo,'functions/python/run.py'),args,'--rscript',rscript)
 result<-system2(python,vapply(cmd,shQuote,character(1)))
 if(result!=0)stop('Reproduction failed; inspect outputs/logs.',call.=FALSE)
 invisible(result)
}
