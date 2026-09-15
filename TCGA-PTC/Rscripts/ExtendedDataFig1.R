# Extended Data Fig.1 — genomic burdens. Source: TP_Comparision.R:621–712.
# Replay uses the saved thyroid_burden intermediate. Its upstream code divides
# TMB by canonical GRCh38 chr1–22,X,Y,M length, scales PGA by 100 and transforms
# SV/TE counts with log2(n+1). See input_manifest for exact upstream inputs.
run_ExtendedDataFig1 <- function(input_dir,output_dir) {
  suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(scales);library(hrbrthemes);library(cowplot)})
  ptc_font();theme_set(ptc_theme())
  ptc_load('thyroid_burden_tp.RData',input_dir,objects='thyroid_burden')
  ptc_load('wgsdata_tp.RData',input_dir,objects=c('wgsdata','study_color2'))
  dat <- thyroid_burden |> filter(Tumor_Barcode %in% wgsdata$Tumor_Barcode)
  list('Extended Data Fig. 1'=ptc_capture('Extended Data Fig. 1',{
    r <- plot_group_compare(tdata=dat,value=value,group='Study2',group_name='Group',facet='name',
      FDR_facet=TRUE,facet_ncol=4,ylab='Value',hide_ns=TRUE,palette=study_color2,
      pcol_sig='#BB0E3D',facet_scales='free_y',output_file=NULL,width=14,height=7,
      legend_position='bottom',y_pos_adjust=c(1.9,2.5,0.8,65,63,122.5624,4,5,6,3,2,1.5))
    p <- r$plot + ptc_theme()+theme(legend.position='bottom',axis.text.x=element_blank(),axis.ticks.x=element_blank())
    list(files=ptc_save(p,output_dir,'ExtendedDataFig1',14,7),statistics=r$stats)
  }))
}
