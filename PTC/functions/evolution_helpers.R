# Shared mutation-cluster, alteration-matrix and heterogeneity calculations.
ptc_evolution_inputs <- function(input_dir, math=FALSE) {
  e <- new.env(parent=baseenv())
  ptc_load('wgsdata_tp.RData',input_dir,e,c('wgsdata','study_color2'))
  ptc_load('DP_info_data.RData',input_dir,e,if(math)c('DP_info_data','clone_data') else 'clone_data')
  ptc_load('Genome_landscape_manual_final.RData',input_dir,e,'data_top0')
  e
}
ptc_event_matrix <- function(wgsdata, data_top0) {
  # Include categories with more than 10 event rows and a pooled fusion category.
  top <- dplyr::bind_rows(data_top0,dplyr::distinct(dplyr::mutate(dplyr::filter(data_top0,Type=='Fusion'),Gene='All_Gene')))
  top <- dplyr::mutate(top,Key=paste0(Gene,'@',Type))
  keys <- dplyr::count(top,Key) |> dplyr::filter(n>10) |> dplyr::pull(Key)
  events <- top |> dplyr::filter(Key %in% keys) |> dplyr::distinct(Tumor_Barcode,Key) |> dplyr::mutate(value='Yes')
  tidyr::expand_grid(Tumor_Barcode=wgsdata$Tumor_Barcode,Key=keys) |>
    dplyr::left_join(events,by=c('Tumor_Barcode','Key')) |>
    dplyr::mutate(value=tidyr::replace_na(value,'No')) |>
    dplyr::left_join(dplyr::select(wgsdata,Tumor_Barcode,Study2),by='Tumor_Barcode')
}
ptc_math <- function(wgsdata, DP_info_data) {
  # R mad() retains its default normal-consistency multiplier (1.4826).
  DP_info_data |> dplyr::filter(Analysis_Barcode %in% wgsdata$Analysis_Barcode) |>
    dplyr::group_by(Tumor_Barcode,Analysis_Barcode) |>
    dplyr::summarise(MATH_CCF=100*stats::mad(CCF,na.rm=TRUE)/stats::median(CCF,na.rm=TRUE),
      MATH_VAF=100*stats::mad(VAF,na.rm=TRUE)/stats::median(VAF,na.rm=TRUE),.groups='drop') |>
    dplyr::left_join(dplyr::select(wgsdata,Tumor_Barcode,Analysis_Barcode,Study2),by=c('Tumor_Barcode','Analysis_Barcode'))
}
ptc_fisher <- function(data, group_cols) {
  data |> dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(x,key) {
      fit <- tryCatch(stats::fisher.test(factor(x$Larger_Subclone,levels=c('No','Yes')),factor(x$value,levels=c('No','Yes'))),error=function(e)NULL)
      if(is.null(fit)) return(tibble::tibble(estimate=numeric(),p.value=numeric()))
      tibble::tibble(estimate=unname(fit$estimate),p.value=fit$p.value)
    }) |> dplyr::ungroup()
}
ptc_volcano <- function(data, x, facet=FALSE, bonferroni=FALSE) {
  p <- ggplot2::ggplot(data,ggplot2::aes(x=.data[[x]],y=-log10(p.value),fill=Type))+
    ggplot2::geom_point(shape=21,size=3,stroke=0.2)+
    ggrepel::geom_text_repel(data=dplyr::filter(data,p.value<0.05),ggplot2::aes(label=Gene),size=4.3,family='Roboto Condensed',seed=1)+
    ggplot2::geom_vline(xintercept=0,linewidth=0.2)+
    ggplot2::geom_hline(yintercept=-log10(0.05),linetype=2,linewidth=0.3)+
    ggsci::scale_fill_d3()+ptc_theme()+
    ggplot2::labs(x=if(x=='log2OR') 'log2 odds ratio (Fisher exact test)' else 'log2 fold change',y='−log10(P value)',fill='Alteration')
  if(facet) p <- p+ggplot2::facet_wrap(~Study2)
  # The upper reference line uses the fixed raw-P threshold 0.05/12.
  if(bonferroni) p <- p+ggplot2::geom_hline(yintercept=-log10(0.05/12),linetype=2,colour='#BB0E3D',linewidth=0.3)
  p
}
