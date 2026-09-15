# Figure 4 — mutation-cluster architecture and driver associations.
# Source: TP_Comparision.R:2157–2352. Statistical definitions are preserved.
run_Fig4 <- function(input_dir, output_dir) {
  suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2)})
  ptc_font(); theme_set(ptc_theme())
  e <- ptc_evolution_inputs(input_dir)
  wgsdata <- e$wgsdata; clone_data <- e$clone_data
  trunk <- clone_data |> filter(Clone=='Yes') |> group_by(Tumor_Barcode) |>
    summarise(Clonal_Mutation=sum(nMutations,na.rm=TRUE),.groups='drop')
  cluster <- clone_data |> left_join(trunk,by='Tumor_Barcode') |>
    mutate(flag=case_when(Clone=='Yes'~'Clone',Clone=='No' & nMutations>Clonal_Mutation~'Large_Subclone',TRUE~'Small_Subclone'))
  rich <- unique(cluster$Tumor_Barcode[cluster$flag=='Large_Subclone'])
  dat <- wgsdata |> select(Tumor_Barcode,Study2) |> left_join(cluster,by='Tumor_Barcode') |>
    mutate(Large_Subclone=if_else(Tumor_Barcode %in% rich,'High Subclonal Mutation Burden','Low Subclonal Mutation Burden'))
  results <- list()
  # Fig. 4a: original kernel density and scale limits, no reclassification.
  results[['Fig. 4a']] <- ptc_capture('Fig. 4a',{
    p <- ggplot(dat,aes(CCF,colour=Study2))+geom_density(linewidth=0.8)+
      geom_vline(xintercept=c(0.2,0.265,1),linetype=2,linewidth=0.25)+
      scale_x_continuous(breaks=scales::pretty_breaks(n=6),limits=c(0,1.2))+
      scale_colour_manual(values=e$study_color2)+labs(x='CCF',y='Density',colour='Group')+ptc_theme()
    ptc_save(p,output_dir,'Fig4a',10,6)
  })
  # Fig. 4b: mutation fraction is not cellular fraction; CCF is shown separately.
  results[['Fig. 4b']] <- ptc_capture('Fig. 4b',{
    p <- ggplot(dat,aes(CCF,nMutations/DPClust_Mutations,size=nMutations/DPClust_Mutations,fill=flag))+
      geom_point(shape=21,colour='black',stroke=0.2)+facet_grid(Large_Subclone~Study2)+
      ggsci::scale_fill_npg()+scale_size_binned(breaks=scales::pretty_breaks())+
      scale_x_continuous(breaks=scales::pretty_breaks(n=6),limits=c(0,1.2))+
      labs(x='CCF',y='Proportion of mutations',fill='Clonality',size='Proportion of mutations')+ptc_theme()
    ptc_save(p,output_dir,'Fig4b',14,8)
  })
  events <- ptc_event_matrix(wgsdata,e$data_top0) |>
    mutate(Larger_Subclone=if_else(Tumor_Barcode %in% rich,'Yes','No'))
  for(by_cohort in c(FALSE,TRUE)) {
    id <- if(by_cohort)'Fig. 4d' else 'Fig. 4c'
    results[[id]] <- ptc_capture(id,{
      stat <- ptc_fisher(events,if(by_cohort)c('Key','Study2') else 'Key') |> arrange(p.value)
      if(by_cohort) stat <- group_by(stat,Study2)
      stat <- stat |> mutate(FDR=p.adjust(p.value,method='holm'),log2OR=log2(estimate)) |>
        ungroup() |> separate(Key,c('Gene','Type'),sep='@')
      files <- ptc_save(ptc_volcano(stat,'log2OR',by_cohort,TRUE),output_dir,if(by_cohort)'Fig4d' else 'Fig4c',if(by_cohort)14 else 10,6)
      list(files=files,statistics=stat)
    })
  }
  results[['Fig. 4e']] <- run_Fig4_schematic(input_dir,output_dir)
  results
}
