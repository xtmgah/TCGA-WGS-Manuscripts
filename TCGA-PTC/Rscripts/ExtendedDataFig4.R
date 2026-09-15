# Extended Data Fig4 a/b/c: MRCA age, diagnosis age and APOBEC by BRAF/TERT.
# Sources: TP_Comparision.R 3363-3370 + 3408-3437; 3007-3038; 1640-1689.
run_ExtendedDataFig4 <- function(input_dir, output_dir) {
  ptc_age_packages(c('rstatix','ggpubr','ggbeeswarm','ggsci'))
  x <- ptc_age_common(input_dir)
  results <- list()
  results[['Extended Data Fig. 4a']] <- ptc_age_result('Extended Data Fig. 4a', {
    dat <- ptc_age_mrca(input_dir,x$cohort) %>% filter(Study2=='TCGA-THCA') %>%
      ptc_age_mutant_group(x$alterations,'MRCA_age')
    result <- ptc_age_box(dat,'MRCA_age','Estimated MRCA age',-10)
    list(files=ptc_save(result$plot,output_dir,'ExtendedDataFig4a',10,6),statistics=result$tests,
      notes='TCGA-THCA, acceleration=1x; groups reordered by median MRCA age; BH pairwise Wilcoxon tests.')
  })
  results[['Extended Data Fig. 4b']] <- ptc_age_result('Extended Data Fig. 4b', {
    dat <- x$covariates %>% filter(Study=='TCGA-THCA') %>% ptc_age_mutant_group(x$alterations,'age_at_diagnosis')
    result <- ptc_age_box(dat,'age_at_diagnosis','Age at diagnosis',c(85,86,87,88,89,90))
    list(files=ptc_save(result$plot,output_dir,'ExtendedDataFig4b',10,6),statistics=result$tests,
      notes='TCGA-THCA groups reordered by diagnosis age; BH pairwise Wilcoxon tests.')
  })
  results[['Extended Data Fig. 4c']] <- ptc_age_result('Extended Data Fig. 4c', {
    activity <- ptc_age_load(input_dir,'thyroid_mutational_signatures_tp_final.RData','thyroid_activity_all')
    ratio <- ptc_age_load(input_dir,'thyroid_mutational_signatures_tp_final.RData','thyroid_activity_all_ratio')
    apobec <- activity %>% transmute(Tumor_Barcode,APOBEC=SBS2+SBS13) %>%
      left_join(ratio %>% transmute(Tumor_Barcode,APOBEC_ratio=SBS2+SBS13),by='Tumor_Barcode') %>%
      filter(APOBEC>50,APOBEC_ratio>.05) %>% pull(Tumor_Barcode)
    dat <- x$covariates %>% filter(Study=='TCGA-THCA') %>%
      ptc_age_mutant_group(x$alterations,'age_at_diagnosis') %>%
      mutate(APOBEC=if_else(Tumor_Barcode %in% apobec,'Present','Absent'))
    counts <- dat %>% count(Group,APOBEC) %>% group_by(Group) %>% mutate(Freq=n/sum(n)) %>%
      ungroup() %>% arrange(desc(APOBEC),desc(Freq)) %>% mutate(Group=fct_rev(fct_inorder(as.character(Group))))
    p <- ggplot(counts,aes(Group,n,fill=APOBEC)) + geom_col(position='fill',width=.75) +
      geom_text(aes(label=percent(Freq,accuracy=.1)),position=position_fill(vjust=.5),size=4.5,color='white',family='Roboto Condensed') +
      scale_y_continuous(labels=label_percent(),expand=expansion(add=0)) +
      labs(x=NULL,y='Frequency (%)',fill='APOBEC mutagenesis') + ggsci::scale_fill_jama() +
      ptc_theme() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
    list(files=ptc_save(p,output_dir,'ExtendedDataFig4c',10,6),statistics=counts,
      notes='APOBEC present requires SBS2+SBS13 strictly >50 mutations AND fraction strictly >0.05. Four BRAF/TERT groups within TCGA-THCA.')
  })
  results
}
