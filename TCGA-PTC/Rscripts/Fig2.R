# Fig. 2: chromosome 22q loss. Historical rules are retained explicitly.
# Fig2a uses a recovered manual state switch, not a verified final generator.
# Fig2d deliberately exports the recovered raw-P variant; final embedded PDF
# has an adjusted-P axis, whose exact generating export code is unavailable.
ptc_chr22_inputs <- function(input_dir) {
  ptc_load('wgsdata_tp.RData',input_dir,objects=c('wgsdata','study_color2'))
  ptc_load('ngspurity_data.RData',input_dir,objects=c('ngspurity','BBprofile'))
  cnv <- BBprofile %>% filter(Tumor_Barcode %in% wgsdata$Tumor_Barcode) %>%
    mutate(frac1_A=as.numeric(frac1_A),frac2_A=as.numeric(frac2_A),
           nMaj1_A=as.integer(nMaj1_A),nMaj2_A=as.integer(nMaj2_A),nMin1_A=as.integer(nMin1_A),nMin2_A=as.integer(nMin2_A)) %>%
    filter(chr %in% c(1:22,'X','Y')) %>% mutate(frac1_A=replace_na(frac1_A,0),frac2_A=replace_na(frac2_A,0)) %>%
    mutate(clone_frac=if_else(frac1_A>frac2_A,frac1_A,frac2_A),
           clone_nMaj=if_else(frac1_A>frac2_A,nMaj1_A,nMaj2_A),clone_nMin=if_else(frac1_A>frac2_A,nMin1_A,nMin2_A),
           subclone_frac=if_else(frac1_A<frac2_A,frac1_A,frac2_A),
           subclone_nMaj=if_else(frac1_A<frac2_A,nMaj1_A,nMaj2_A),subclone_nMin=if_else(frac1_A<frac2_A,nMin1_A,nMin2_A)) %>%
    select(Tumor_Barcode,chr:ntot,contains('clone'))
  state_data <- function(state) {
    cnv %>% filter(chr==22) %>% mutate(clone_total=if(state=='clone') clone_nMin+clone_nMaj else subclone_nMin+subclone_nMaj) %>%
      left_join(ngspurity %>% select(Tumor_Barcode,WGD_Status=MCN_WGD,Tumor_Purity=BB_Purity,Tumor_Ploidy=BB_Ploidy),by='Tumor_Barcode') %>%
      mutate(WGD_Status=replace_na(WGD_Status,'nWGD'),relative_copy=clone_total-if_else(WGD_Status=='WGD',4,2),
             relative_copy=pmin(4,pmax(-4,relative_copy)),relative_copy=if_else(clone_nMaj==0,-4,relative_copy),
             startpos=as.integer(startpos),endpos=as.integer(endpos))
  }
  clone <- state_data('clone');subclone <- state_data('subclone')
  combined <- bind_rows(clone %>% mutate(chr='Clone'),subclone %>% mutate(chr='Subclone'))
  # Exact historical classification uses the longest chr22 segment, including
  # its missing-to-neutral rule. It is not a newly introduced arm-overlap rule.
  states <- combined %>% group_by(Tumor_Barcode,chr) %>% mutate(relative_copy=replace_na(relative_copy,0)) %>%
    arrange(desc(abs(endpos-startpos))) %>% slice(1) %>% ungroup() %>%
    select(Tumor_Barcode,chr,relative_copy) %>% arrange(relative_copy) %>% pivot_wider(names_from=chr,values_from=relative_copy)
  sample_order <- states %>% left_join(wgsdata %>% select(Study2,Tumor_Barcode),by='Tumor_Barcode') %>%
    arrange(Study2,desc(Clone<0),desc(Subclone<0),Clone,Subclone) %>% mutate(Seq=seq_along(Tumor_Barcode)) %>% select(Study2,Tumor_Barcode,Seq)
  list(wgsdata=wgsdata,colors=study_color2,ngspurity=ngspurity,cnv=cnv,clone=clone,subclone=subclone,
       combined=combined,order=sample_order,clonal_loss=states$Tumor_Barcode[states$Clone<0],subclonal_loss=states$Tumor_Barcode[states$Subclone<0])
}

run_Fig2 <- function(input_dir,output_dir,panels=letters[1:4]) {
  suppressPackageStartupMessages({
    library(dplyr);library(tidyr);library(tibble);library(forcats);library(stringr)
    library(ggplot2);library(scales);library(cowplot);library(hrbrthemes)
    library(purrr);library(broom);library(rlang)
  })
  source('functions/group_comparison.R',local=TRUE)
  source('functions/cnv_helpers.R',local=TRUE)
  ptc_font()
  old_null_device <- cowplot::set_null_device('agg')
  on.exit(cowplot::set_null_device(old_null_device),add=TRUE)
  results <- list()
  if ('a' %in% panels) results[['Fig. 2a']] <- ptc_capture('Fig. 2a', {
    suppressPackageStartupMessages({library(CNTools);library(valr)})
    x <- ptc_chr22_inputs(input_dir)
    clone <- ptc_cnv_frequency(x$cnv,x$ngspurity,'data/refs/hg38_cytoBand.txt.gz','clone')
    subclone <- ptc_cnv_frequency(x$cnv,x$ngspurity,'data/refs/hg38_cytoBand.txt.gz','subclone')
    p <- cowplot::plot_grid(clone,subclone,ncol=1,align='v',labels=c('Clonal','Subclonal'),label_fontfamily='Roboto Condensed',label_size=14)
    list(paths=ptc_save(p,output_dir,'Fig2a_state_reconstruction',14,8),
         scope='Explicit reconstruction of historical clone/subclone state switch; original saved script assigns both rows the same plot. Requires author confirmation before final reproduction claim.')
  })
  if ('b' %in% panels) results[['Fig. 2b']] <- ptc_capture('Fig. 2b', {
    x <- ptc_chr22_inputs(input_dir)
    dat <- x$combined %>% left_join(x$order,by='Tumor_Barcode') %>%
      mutate(Study2=fct_recode(Study2,'Chornobyl-\nUnexposed'='Chornobyl-Unexposed','Chornobyl-\nExposed'='Chornobyl-Exposed'))
    p <- ggplot(dat,aes(fill=relative_copy)) + geom_rect(aes(xmin=startpos,xmax=endpos,ymin=Seq-.5,ymax=Seq+.5)) +
      facet_grid(fct_rev(Study2)~chr,scales='free_y',space='free',switch='both') +
      scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,na.value='white') +
      scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
      labs(x=NULL,y=NULL,fill='Relative copy number') + ptc_theme() +
      theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.spacing=grid::unit(.05,'lines'),
            strip.placement='outside',strip.text.y.left=element_text(angle=0),strip.background=element_rect(color='gray50',fill='white',linewidth=.3),legend.position='bottom')
    list(paths=ptc_save(p,output_dir,'Fig2b',8,10),clonal_loss_n=length(x$clonal_loss),subclonal_loss_n=length(x$subclonal_loss),
         scope='Exact historical chr22 classification and heatmap data; includes neutral samples as white rows.')
  })
  if ('c' %in% panels) results[['Fig. 2c']] <- ptc_capture('Fig. 2c', {
    x <- ptc_chr22_inputs(input_dir)
    dat <- x$subclone %>% filter(Tumor_Barcode %in% x$subclonal_loss,!is.na(subclone_frac),subclone_frac>0,relative_copy<0) %>%
      left_join(x$wgsdata %>% select(Tumor_Barcode,Study2),by='Tumor_Barcode')
    p <- plot_group_compare(dat,subclone_frac,group='Study2',group_name='Group',ylab='Subclone fraction',palette=x$colors,p_off=TRUE)$plot
    list(paths=ptc_save(p,output_dir,'Fig2c',9,6),n_segment_observations=nrow(dat),n_unique_tumors=n_distinct(dat$Tumor_Barcode),
         scope='Historical segment-level observations retained; compare counts with legend describing individual tumors.')
  })
  if ('d' %in% panels) results[['Fig. 2d']] <- ptc_capture('Fig. 2d', {
    x <- ptc_chr22_inputs(input_dir)
    ptc_load('Genome_landscape_manual_final.RData',input_dir,objects='data_top0')
    labels <- data_top0 %>% filter(Gene!='chr22q') %>% mutate(label=paste(Gene,Alteration,Type,sep='@')) %>% select(Tumor_Barcode,label) %>% mutate(value=TRUE)
    event_matrix <- x$wgsdata %>% select(Tumor_Barcode) %>% left_join(labels,by='Tumor_Barcode') %>% pivot_wider(names_from=label,values_from=value,values_fill=FALSE) %>% select(-any_of('NA')) %>% pivot_longer(-Tumor_Barcode)
    fit_state <- function(loss,label) {
      event_matrix %>% mutate(value2=Tumor_Barcode %in% loss) %>% group_by(name) %>%
        group_modify(~{fit <- tryCatch(fisher.test(.x$value,.x$value2),error=function(e)NULL);if(is.null(fit))tibble() else broom::tidy(fit)}) %>%
        ungroup() %>% arrange(p.value) %>% mutate(FDR=p.adjust(p.value,method='holm')) %>%
        separate(name,into=c('Gene','Alteration','Type'),sep='@') %>% mutate(Group=label)
    }
    dat <- bind_rows(fit_state(x$clonal_loss,'Clonal chr22q deletion'),fit_state(x$subclonal_loss,'Subclonal chr22q deletion'))
    p <- ggplot(dat,aes(log2(estimate),-log10(p.value),fill=Type)) + geom_point(shape=21,size=3,stroke=.2) +
      ggrepel::geom_text_repel(data=dat %>% filter(FDR<.05),aes(label=Gene),size=4.5,family='Roboto Condensed',seed=1) +
      geom_vline(xintercept=0,linewidth=.2) + facet_wrap(~Group) + ggsci::scale_fill_d3() +
      scale_x_continuous(breaks=pretty_breaks()) + scale_y_continuous(breaks=pretty_breaks()) +
      labs(x=expression(log[2]('Odds ratio')),y=expression(-log[10]('P value')),fill='Alteration') + ptc_theme()
    list(paths=ptc_save(p,output_dir,'Fig2d_related_raw_p',12,7),n_tests=nrow(dat),adjustment='Holm, matching stats::p.adjust default in historical source',
         scope='Recovered historical raw-P variant only. Final embedded PDF uses adjusted-P values; exact final export code was not recovered.')
  })
  results
}
