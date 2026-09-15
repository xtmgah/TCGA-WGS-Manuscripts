# Common runtime. No source data are bundled; paths are passed explicitly.
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
  p <- file.path(input_dir,file)
  if (!file.exists(p)) stop('Required private input is missing: ', file)
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
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  paths <- setNames(file.path(output_dir,paste0(id,c('.pdf','.png'))),c('pdf','png'))
  ggplot2::ggsave(paths[['pdf']],plot=plot,width=width,height=height,units='in',device=grDevices::cairo_pdf,bg='white',limitsize=FALSE)
  ggplot2::ggsave(paths[['png']],plot=plot,width=width,height=height,units='in',device=ragg::agg_png,dpi=300,bg='white',limitsize=FALSE)
  paths
}
ptc_capture <- function(id, expr) {
  warnings <- character()
  tryCatch({
    value <- withCallingHandlers(eval(substitute(expr),envir=parent.frame()),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    list(panel_id=id,status='generated_unverified',result=value,warnings=unique(warnings))
  },error=function(e) list(panel_id=id,status='blocked',reason=conditionMessage(e),warnings=unique(warnings)))
}
ptc_blocked <- function(id, reason, status='source code partially identified') list(panel_id=id,status=status,reason=reason)
