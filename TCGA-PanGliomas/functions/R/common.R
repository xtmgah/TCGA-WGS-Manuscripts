# Display-only publication styling; no statistical fitting or tests.
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(jsonlite)})
REF_ROOT <- normalizePath(Sys.getenv('PANGLIOMA_WORKSPACE'))
REF_WORK <- REF_ROOT
REF_STAGE <- file.path(REF_ROOT,'stage')
REF_STYLE <- fromJSON(file.path(Sys.getenv('PANGLIOMA_PY_CODE'),'style.json'))
REF_FONT <- REF_STYLE$font_family
if (!any(systemfonts::system_fonts()$family==REF_FONT)) stop('Roboto Condensed is required')
ggplot2::update_geom_defaults('text',list(family=REF_FONT,size=6/ggplot2::.pt))
ggplot2::update_geom_defaults('label',list(family=REF_FONT,size=6/ggplot2::.pt))
ref_theme <- function(grid=FALSE) {
 theme_classic(base_family=REF_FONT,base_size=6.5)+theme(
  text=element_text(family=REF_FONT,color='#222222'),
  axis.text=element_text(size=6.5,color='#303030'),axis.title=element_text(size=6.5),
  axis.line=element_line(linewidth=.22,color='#505050'),axis.ticks=element_line(linewidth=.22,color='#505050'),
  axis.ticks.length=unit(1.2,'pt'),
  panel.grid.major.y=if(grid)element_line(linewidth=.18,color='#E8E8E8')else element_blank(),
  panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
  plot.title.position='plot',
  plot.title=element_text(size=8,hjust=0,face='plain',margin=margin(b=2)),
  plot.subtitle=element_text(size=6.5,hjust=0,margin=margin(b=2)),
  strip.background=element_blank(),strip.text=element_text(size=7.5),
  legend.title=element_text(size=6.5),legend.text=element_text(size=6.5),
  legend.key.height=unit(7,'pt'),legend.key.width=unit(10,'pt'),
  legend.margin=margin(0,0,0,0),legend.box.spacing=unit(1,'pt'),
  plot.margin=margin(1,2,1,2),plot.background=element_rect(fill='white',color=NA))
}
ref_save <- function(p,id,width_mm,height_mm) {
 ggsave(file.path(REF_STAGE,'panels',paste0(id,'.pdf')),p,width=width_mm,height=height_mm,units='mm',device=cairo_pdf,bg='white',limitsize=FALSE)
}
ref_p <- function(x,prefix='P') {
 vapply(x,function(v) {
  if(is.na(v))return(paste0(prefix,' = NA'))
  if(v<.001)return(paste0(prefix,' < 0.001'))
  paste0(prefix,' = ',formatC(v,format='f',digits=if(v<.01)3 else 2))
 },character(1))
}
