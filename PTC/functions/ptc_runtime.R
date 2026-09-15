# Shared input, output, reference and plotting support.
.ptc_runtime_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork=TRUE), error=function(e) NULL)
.ptc_package_root <- if(length(.ptc_runtime_file)==1L) dirname(dirname(.ptc_runtime_file)) else getwd()

ptc_resolve_path <- function(path) {
  if(length(path)!=1L || is.na(path) || !nzchar(path)) stop('Specify a nonempty directory path.')
  path <- path.expand(path)
  if(!grepl('^(/|[A-Za-z]:)',path)) path <- file.path(getwd(),path)
  tail_parts <- character()
  ancestor <- path
  while(!file.exists(ancestor) && !dir.exists(ancestor)) {
    parent <- dirname(ancestor)
    if(identical(parent,ancestor)) stop('Cannot resolve directory path: ',path)
    tail_parts <- c(basename(ancestor),tail_parts)
    ancestor <- parent
  }
  resolved <- normalizePath(ancestor,winslash='/',mustWork=TRUE)
  for(part in tail_parts) {
    if(part=='.') next
    resolved <- if(part=='..') dirname(resolved) else file.path(resolved,part)
    if(file.exists(resolved) || dir.exists(resolved))
      resolved <- normalizePath(resolved,winslash='/',mustWork=TRUE)
  }
  resolved
}
ptc_repository_root <- function() {
  path <- ptc_resolve_path(.ptc_package_root)
  repeat {
    if(file.exists(file.path(path,'.git')) || dir.exists(file.path(path,'.git'))) return(path)
    parent <- dirname(path)
    if(identical(path,parent)) return(ptc_resolve_path(.ptc_package_root))
    path <- parent
  }
}
ptc_external_directory <- function(path,label,must_exist=TRUE,create=FALSE) {
  resolved <- ptc_resolve_path(path)
  repository <- ptc_repository_root()
  check_location <- function(x) {
    if(identical(x,repository) || startsWith(x,paste0(repository,'/')))
      stop(label,' must be outside the Git repository: ',x)
  }
  check_location(resolved)
  if(must_exist && !dir.exists(resolved)) stop(label,' does not exist: ',resolved)
  if(create && !dir.exists(resolved)) {
    if(!dir.create(resolved,recursive=TRUE,showWarnings=FALSE)) stop('Cannot create ',label,': ',resolved)
    resolved <- ptc_resolve_path(resolved)
    check_location(resolved)
  }
  resolved
}
ptc_reference <- function(file,reference_dir=Sys.getenv('PTC_REFERENCE_DIR')) {
  if(!nzchar(reference_dir)) stop('Set PTC_REFERENCE_DIR to an external directory containing the required GRCh38 reference files.')
  reference_dir <- ptc_external_directory(reference_dir,'Reference directory')
  path <- file.path(reference_dir,file)
  if(!file.exists(path)) stop('Required GRCh38 reference is missing: ',file)
  path <- ptc_resolve_path(path)
  ptc_external_directory(dirname(path),'Reference file directory')
  path
}
ptc_font <- function() {
  if (!requireNamespace('systemfonts', quietly=TRUE) ||
      !requireNamespace('showtext', quietly=TRUE) ||
      !requireNamespace('sysfonts', quietly=TRUE)) stop('Install systemfonts, showtext and sysfonts.')
  f <- systemfonts::system_fonts()
  f <- f[f$family == 'Roboto Condensed', , drop=FALSE]
  if (!nrow(f)) stop('Roboto Condensed is required. Install the font before plotting.')
  choose <- function(style) {
    z <- f$path[f$style == style]
    if (!length(z)) stop('Missing Roboto Condensed font style: ', style)
    z[[1]]
  }
  sysfonts::font_add('Roboto Condensed', regular=choose('Regular'), bold=choose('Bold'), italic=choose('Italic'), bolditalic=choose('Bold Italic'))
  showtext::showtext_auto()
  showtext::showtext_opts(dpi=300)
  ggplot2::update_geom_defaults('text', list(family='Roboto Condensed', size=4.5))
  ggplot2::update_geom_defaults('label', list(family='Roboto Condensed', size=4.5))
  invisible(f$path)
}
ptc_theme <- function() {
  ggplot2::theme_minimal(base_family='Roboto Condensed', base_size=13) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=18,face='bold',margin=ggplot2::margin(b=8)),
      plot.subtitle=ggplot2::element_text(size=13,margin=ggplot2::margin(b=10)),
      axis.title=ggplot2::element_text(size=13),axis.text=ggplot2::element_text(size=11),
      legend.title=ggplot2::element_text(size=12,face='bold'),legend.text=ggplot2::element_text(size=11),
      strip.text=ggplot2::element_text(size=12,face='bold'),plot.caption=ggplot2::element_text(size=9),
      panel.grid.major=ggplot2::element_blank(),panel.grid.minor=ggplot2::element_blank(),
      panel.spacing=grid::unit(0.8,'lines'),plot.margin=ggplot2::margin(8,12,8,8))
}
ptc_load <- function(file, input_dir, envir=parent.frame(), objects=NULL) {
  input_dir <- ptc_external_directory(input_dir,'Input directory')
  p <- file.path(input_dir,file)
  if (!file.exists(p)) stop('Required private input is missing: ', file)
  p <- ptc_resolve_path(p)
  ptc_external_directory(dirname(p),'Input file directory')
  e <- new.env(parent=baseenv())
  available <- load(p,envir=e)
  wanted <- if(is.null(objects)) available else objects
  absent <- setdiff(wanted,available)
  if(length(absent)) stop('Missing object(s) in ', file, ': ',paste(absent,collapse=', '))
  for(n in wanted) assign(n,get(n,e),envir=envir)
  invisible(wanted)
}
ptc_save <- function(plot, output_dir, id, width=10, height=6) {
  if (!grepl('^[A-Za-z0-9_-]+$',id)) stop('Unsafe output identifier: ',id)
  output_dir <- ptc_external_directory(output_dir,'Output directory',must_exist=FALSE,create=TRUE)
  paths <- setNames(file.path(output_dir,paste0(id,c('.pdf','.png'))),c('pdf','png'))
  ggplot2::ggsave(paths[['pdf']],plot=plot,width=width,height=height,units='in',device=grDevices::cairo_pdf,bg='white',limitsize=FALSE)
  ggplot2::ggsave(paths[['png']],plot=plot,width=width,height=height,units='in',device=ragg::agg_png,dpi=300,bg='white',limitsize=FALSE)
  paths
}
ptc_capture <- function(id, expr) {
  warnings <- character()
  tryCatch({
    value <- withCallingHandlers(eval(substitute(expr),envir=parent.frame()),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    list(panel_id=id,status='generated',result=value,warnings=unique(warnings))
  },error=function(e) list(panel_id=id,status='unavailable',reason=conditionMessage(e),warnings=unique(warnings)))
}
ptc_blocked <- function(id, reason, status='unavailable') list(panel_id=id,status=status,reason=reason)
