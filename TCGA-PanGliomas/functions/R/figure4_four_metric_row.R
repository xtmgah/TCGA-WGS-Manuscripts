source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'plot_helpers_34.R'))
cairo_pdf(tempfile(fileext='.pdf'),width=7,height=10)
cowplot::set_null_device(function(width,height)grDevices::cairo_pdf(tempfile(fileext='.pdf'),width=width,height=height))
D <- 'results/analysis/bcor_reclassification_2026-09-21/figures/figure4/derived-data/publication-panels'
mv <- display_read(file.path(D,'figure4_six_metric_values.tsv'))
mt <- display_read(file.path(D,'figure4_six_metric_pairwise_tests.tsv'))
spec <- data.frame(letter=letters[8:11],metric=c('Latency','MATH','PGA','Ploidy'),ylab=c('Years','Score','Genome altered (%)','Copies'))
box_audit <- list()
for(i in seq_len(nrow(spec))) {
 metric <- spec$metric[i]
 z <- mv |> filter(.data$metric==.env$metric) |> transmute(group=unname(dn_names[DN_Group]),value=value)
 a <- mt |> filter(.data$metric==.env$metric,q_value<.05) |> transmute(g1=unname(dn_names[group1]),g2=unname(dn_names[group2]),label=ref_p(q_value,'q'))
 p <- distribution(z,metric,spec$ylab[i],dn_cols,a,mean_diamond=TRUE,zero=metric!='Ploidy')+theme(plot.margin=margin(3,4,2,4))
 before <- ggplot_build(p)$data[[1]]
 p$layers[[1]] <- geom_boxplot(width=.30,outlier.shape=NA,linewidth=.21,color='#333333',fill='white',alpha=1)
 after <- ggplot_build(p)$data[[1]]
 keep <- c('ymin','lower','middle','upper','ymax','outliers')
 stopifnot(identical(before[keep],after[keep]))
 box_audit[[metric]] <- list(values=after[keep],width=unique(after$xmax-after$xmin),quantiles_and_outliers_unchanged=TRUE)
 ref_save_original_34(p,paste0('M4_',spec$letter[i]),41.5,29.5)
}
ref_export_semantics('M4_four_metric_row',list(no_models_or_tests_run=TRUE,metric_counts=mv |> filter(metric %in% spec$metric) |> count(metric,DN_Group),pairwise_tests=mt |> filter(metric %in% spec$metric),box_audit=box_audit))
