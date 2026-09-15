# Extended Data Figure 2: cohort-specific CCF-based MATH associations.
# The same analysis is also shown in Supplementary Figure 3.
run_ExtendedDataFig2 <- function(input_dir,output_dir,prefix='ExtendedDataFig2',figure_label='Extended Data Fig. 2') {
  suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2)})
  ptc_font();theme_set(ptc_theme())
  e <- ptc_evolution_inputs(input_dir,math=TRUE)
  math <- ptc_math(e$wgsdata,e$DP_info_data)
  dat <- ptc_event_matrix(e$wgsdata,e$data_top0) |>
    left_join(select(math,Tumor_Barcode,MATH_CCF,MATH_VAF),by='Tumor_Barcode')
  stat <- dat |> group_by(Key,Study2) |> group_modify(function(x,key){
    fit <- tryCatch(wilcox.test(MATH_CCF~value,data=x),error=function(e)NULL)
    if(is.null(fit)) return(tibble(FC=numeric(),p.value=numeric()))
    tibble(FC=log2(median(x$MATH_CCF[x$value=='Yes'])/median(x$MATH_CCF[x$value=='No'])),p.value=fit$p.value)
  }) |> ungroup() |> arrange(p.value) |> group_by(Study2) |>
    mutate(FDR=p.adjust(p.value,method='holm')) |> ungroup() |>
    separate(Key,c('Gene','Type'),sep='@')
  results <- list()
  results[[paste0(figure_label,'a')]] <- ptc_capture(paste0(figure_label,'a'),{
    files <- ptc_save(ptc_volcano(stat,'FC',facet=TRUE),output_dir,paste0(prefix,'a'),14,6)
    list(files=files,statistics=stat)
  })
  # Select alterations with nominal P < 0.05 in at least one cohort.
  selected <- stat |> filter(p.value<0.05) |> distinct(Gene,Type) |>
    mutate(Key=paste0(Gene,'@',Type)) |> arrange(Type)
  dd <- selected |> left_join(dat,by='Key') |> mutate(Key=forcats::fct_inorder(Key))
  results[[paste0(figure_label,'b')]] <- ptc_capture(paste0(figure_label,'b'),{
    p <- ggplot(dd,aes(value,MATH_CCF,fill=Study2))+
      geom_violin(trim=FALSE,width=0.7)+geom_boxplot(outlier.colour=NA,width=0.1,fill='white')+
      facet_grid(Study2~Key)+scale_fill_manual(values=e$study_color2)+
      labs(x='Alteration status',y='Tumor heterogeneity MATH score (CCF)',fill='Group')+ptc_theme()
    ptc_save(p,output_dir,paste0(prefix,'b'),max(14,3.2*nrow(selected)),10)
  })
  results
}
