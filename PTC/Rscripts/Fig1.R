# Figure 1: cohort characteristics, driver mutations and structural variation.
# Run from PTC after loading the shared functions.
run_Fig1 <- function(input_dir, output_dir, panels=letters[1:9]) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(forcats); library(stringr)
    library(tibble); library(ggplot2); library(scales); library(cowplot)
    library(hrbrthemes); library(rlang); library(data.table)
  })
  source('functions/group_comparison.R', local = TRUE)
  source('functions/landscape_helpers.R', local = TRUE)
  ptc_font()
  old_null_device <- cowplot::set_null_device('agg')
  on.exit(cowplot::set_null_device(old_null_device),add=TRUE)
  results <- list()

  if ('a' %in% panels) results[['Fig. 1a']] <- ptc_capture('Fig. 1a', {
    ptc_load('wgsdata_tp.RData', input_dir, objects=c('wgsdata','study_color2'))
    ptc_load('wgs_clinical.RData', input_dir, objects='wgs_clinical')
    ptc_load('wgs_metrics.RData', input_dir, objects='wgs_metrics')
    ptc_load('ngspurity_data.RData', input_dir, objects='ngspurity')
    age <- wgsdata %>% left_join(wgs_clinical %>% select(Subject,value=age_at_diagnosis),by='Subject')
    age_plot <- plot_group_compare(age,value,group='Study2',group_name='Group',ylab='Age at diagnosis',palette=study_color2,hide_ns=FALSE)$plot
    depth <- bind_rows(
      wgsdata %>% select(Study2,Barcode=Tumor_Barcode,Type=Tumor_Type) %>% left_join(wgs_metrics %>% select(Barcode,value=MEDIAN_COVERAGE),by='Barcode'),
      wgsdata %>% select(Study2,Barcode=Normal_Barcode,Type=Normal_Type) %>% left_join(wgs_metrics %>% select(Barcode,value=MEDIAN_COVERAGE),by='Barcode') %>% mutate(Type='NB/NT'))
    depth_plot <- plot_group_compare(depth,value,group='Study2',facet='Type',facet_scales='free_y',ylab='Sequencing depth',hide_ns=FALSE,palette=study_color2)$plot
    purity <- wgsdata %>% select(Study2,Tumor_Barcode) %>% left_join(ngspurity %>% select(Tumor_Barcode,value=BB_Purity),by='Tumor_Barcode')
    purity_plot <- plot_group_compare(purity,value,group='Study2',ylab='Tumor purity',hide_ns=FALSE,palette=study_color2)$plot
    # Export age, sequencing-depth and purity components separately.
    list(age=ptc_save(age_plot,output_dir,'Fig1a_age',8,6),
         depth=ptc_save(depth_plot,output_dir,'Fig1a_depth',10,6),
         purity=ptc_save(purity_plot,output_dir,'Fig1a_purity',8,6),
         scope='Age, sequencing-depth and purity components; combined overview layout is not available.')
  })

  if ('b' %in% panels) results[['Fig. 1b']] <- ptc_capture('Fig. 1b', {
    ptc_load('wgsdata_tp.RData',input_dir,objects=c('wgsdata','study_color2'))
    ptc_load('wgs_covdata_tp.RData',input_dir,objects='wgs_covdata')
    ptc_load('Genome_landscape_manual_final.RData',input_dir,objects='data_top0')
    ptc_load('intogene_drivers.RData',input_dir,objects='intogene')
    events <- data_top0 %>% filter(Tumor_Barcode %in% wgs_covdata$Tumor_Barcode,Type=='Mutation_Driver',Gene %in% intogene$SYMBOL) %>% select(Tumor_Barcode,Gene) %>% distinct()
    overall <- events %>% count(Gene,name='nMutSample') %>% mutate(Frequency=nMutSample/length(wgsdata$Tumor_Barcode))
    cohort <- events %>% left_join(wgsdata %>% select(Tumor_Barcode,Study2),by='Tumor_Barcode') %>% count(Study2,Gene,name='nMutSample') %>% left_join(wgsdata %>% count(Study2,name='Sample_Size'),by='Study2') %>% mutate(Freq=nMutSample/Sample_Size)
    dat <- intogene %>% left_join(overall %>% select(SYMBOL=Gene,Frequency),by='SYMBOL') %>%
      left_join(cohort %>% select(Study2,SYMBOL=Gene,Freq) %>% pivot_wider(names_from=Study2,values_from=Freq,values_fill=0),by='SYMBOL') %>%
      mutate(SYMBOL=fct_reorder(SYMBOL,Frequency)) %>% arrange(SYMBOL) %>%
      mutate(ROLE_color=case_when(ROLE=='Act'~'#01665e',ROLE=='LoF'~'#C50084FF',TRUE~'#F0AB00FF'),
             SYMBOL_lab=fct_inorder(paste0("<span style='color:",ROLE_color,"'>",SYMBOL,'</span>')))
    cohort_lines <- dat %>% select(SYMBOL_lab,all_of(names(study_color2))) %>%
      pivot_longer(-SYMBOL_lab,names_to='Cohort',values_to='cohort_frequency')
    p <- ggplot(dat) + geom_col(aes(SYMBOL_lab,Frequency),fill='#cccccc',color='black',linewidth=.2) +
      geom_point(data=cohort_lines,aes(SYMBOL_lab,cohort_frequency,color=Cohort)) +
      geom_line(data=cohort_lines,aes(SYMBOL_lab,cohort_frequency,color=Cohort,group=Cohort)) +
      scale_color_manual(values=study_color2,breaks=names(study_color2))
    p <- p + scale_y_continuous(breaks=pretty_breaks(),labels=percent_format()) +
      labs(x=NULL,y='Mutation frequency',color='Cohort') + ptc_theme() +
      theme(axis.text.x=ggtext::element_markdown(angle=90,hjust=1,vjust=.5,face='bold.italic',size=11),
            axis.text.x.bottom=ggtext::element_markdown(angle=90,hjust=1,vjust=.5,face='bold.italic',size=11),legend.position='top')
    list(paths=ptc_save(p,output_dir,'Fig1b',10,6),n_cohort=nrow(wgsdata),n_driver_genes=nrow(dat),scope='Driver frequencies calculated from the supplied IntOGen and curated event tables.')
  })

  if ('c' %in% panels) results[['Fig. 1c']] <- ptc_capture('Fig. 1c', {
    ptc_load('wgsdata_tp.RData',input_dir,objects=c('wgsdata','study_color2'))
    sv <- readr::read_tsv(file.path(input_dir,'1170_manta_meerkat_union_window50_highquality.txt'),show_col_types=FALSE)
    fai <- readr::read_tsv(ptc_reference('hg38_primary.fai'),col_names=FALSE,show_col_types=FALSE) %>%
      transmute(chr=str_remove(X1,'chr'),len=X2) %>% slice(1:24) %>% mutate(start2=cumsum(c(0,head(len,-1)))+1,end2=cumsum(len))
    bins <- bind_rows(lapply(seq_len(nrow(fai)),function(i) {
      starts <- seq(1,fai$len[i],by=5000000)
      tibble(chr=fai$chr[i],start=starts,end=c(starts[-1]-1,fai$len[i]),start2=starts+fai$start2[i]-1,end2=c(starts[-1]-1,fai$len[i])+fai$start2[i]-1)
    }))
    outputs <- list()
    for (cohort_name in names(study_color2)) {
      samples <- wgsdata %>% filter(Study2==cohort_name) %>% pull(Analysis_Barcode)
      events <- sv %>% filter(SAMPLE %in% samples,!is.na(CHROM1),CHROM1 %in% fai$chr,CHROM2 %in% fai$chr)
      breaks <- bind_rows(events %>% transmute(chrom=as.character(CHROM1),start=POS1),events %>% transmute(chrom=as.character(CHROM2),start=POS2)) %>% mutate(end=start+1)
      hit <- valr::bed_intersect(bins %>% select(chrom=chr,start,end),breaks,suffix=c('','.y')) %>% count(chrom,start,end) %>% rename(chr=chrom)
      dat <- bins %>% left_join(hit,by=c('chr','start','end')) %>% mutate(n=replace_na(n,0L)/length(samples),pos=(start2+end2)/2)
      p <- ggplot(dat,aes(pos,n)) + geom_area(fill='gray30',linewidth=.5) + geom_line(linewidth=.1) +
        scale_x_continuous(breaks=(fai$start2+fai$end2)/2,labels=fai$chr,expand=c(0,0)) +
        scale_y_continuous(limits=c(0,.5),breaks=pretty_breaks(),expand=c(0,0)) +
        labs(x='Chromosome',y='Breakpoints per 5 Mb per sample',title=cohort_name) + ptc_theme()
      outputs[[cohort_name]] <- ptc_save(p,output_dir,paste0('Fig1c_',gsub('[^A-Za-z0-9]','_',cohort_name)),12,4)
    }
    list(paths=outputs,scope='Per-cohort breakpoint profiles; combined layout and kinase annotations are not included.')
  })

  for (panel in intersect(c('d','e','f'),panels)) {
    id <- paste0('Fig. 1',panel)
    cohort_name <- c(d='TCGA-THCA',e='Chornobyl-Unexposed',f='Chornobyl-Exposed')[[panel]]
    results[[id]] <- ptc_capture(id, {
      ptc_load('wgsdata_tp.RData',input_dir,objects='wgsdata')
      ptc_load('svdata.RData',input_dir,objects='svdata')
      sv_colors <- c(DEL='#007BBD',DUP='#BB0E3D',TRA='#01665e',h2hINV='#947100',t2tINV='#984ea3')
      dat <- svdata %>% mutate(chrA=CHROM1,chrB=CHROM2,posA=POS1,posB=POS2) %>%
        filter(chrA %in% c(1:22,'X'),chrB %in% c(1:22,'X')) %>%
        mutate(chrA=paste0('chr',chrA),chrB=paste0('chr',chrB)) %>%
        filter(Analysis_Barcode %in% wgsdata$Analysis_Barcode) %>%
        left_join(wgsdata %>% select(Analysis_Barcode,Study2),by='Analysis_Barcode') %>% filter(Study2==cohort_name)
      bed1 <- as.data.frame(dat %>% transmute(chr=chrA,start=posA,end=posA+1000))
      bed2 <- as.data.frame(dat %>% transmute(chr=chrB,start=posB,end=posB+1000))
      cytoband <- read.table(gzfile(ptc_reference('hg38_cytoBand.txt.gz')),sep='\t',header=FALSE,stringsAsFactors=FALSE)
      draw <- function() {
        graphics::par(family='Roboto Condensed')
        circlize::circos.clear()
        on.exit(circlize::circos.clear())
        circlize::circos.initializeWithIdeogram(cytoband=cytoband,chromosome.index=paste0('chr',c(1:22,'X','Y')))
        circlize::circos.genomicLink(bed1,bed2,col=scales::alpha(unname(sv_colors[dat$SVTYPE]),1),border=NA,directional=-1,arr.length=0,lwd=.1)
      }
      output_dir <- ptc_external_directory(output_dir,'Output directory',must_exist=FALSE,create=TRUE)
      paths <- setNames(file.path(output_dir,paste0('Fig1',panel,c('.pdf','.png'))),c('pdf','png'))
      grDevices::cairo_pdf(paths[['pdf']],width=8,height=8,family='Roboto Condensed')
      tryCatch(draw(),finally=grDevices::dev.off())
      ragg::agg_png(paths[['png']],width=8,height=8,units='in',res=300)
      tryCatch(draw(),finally=grDevices::dev.off())
      list(paths=paths,n_events=nrow(dat),scope='Cohort-specific structural-variant links with external GRCh38 cytobands.')
    })
  }

  for (panel in intersect(c('g','h','i'),panels)) {
    id <- paste0('Fig. 1',panel)
    cohort_name <- c(g='TCGA-THCA',h='Chornobyl-Exposed',i='Chornobyl-Unexposed')[[panel]]
    results[[id]] <- ptc_capture(id, {
      ptc_load('wgsdata_tp.RData',input_dir,objects='wgsdata')
      ptc_load('Genome_landscape_manual_final.RData',input_dir,objects=c('data_top0','data_tmb0'))
      colors <- readr::read_csv('config/oncoplot_colors.csv',show_col_types=FALSE)
      landscape_colors <- setNames(colors$Color,colors$Name)
      landscape_colors['Nonsense_Mutation'] <- '#FACE00';landscape_colors['Fusion'] <- '#BB0E3D'
      data_top0 <- data_top0 %>% filter(Tumor_Barcode %in% wgsdata$Tumor_Barcode)
      overall <- oncoplot(data_top0 %>% select(-Study),sample_level0=wgsdata$Tumor_Barcode,landscape_colors=landscape_colors,GeneSortOnly=TRUE,sample_name=FALSE,cell_height=1,p2_hidden=TRUE)
      fusion_genes <- data_top0 %>% count(Gene,Type) %>% pivot_wider(names_from=Type,values_from=n,values_fill=0) %>% filter(Fusion>Mutation_Driver) %>% pull(Gene)
      g <- overall$gene_level
      gene_order <- unique(c(g[g=='chr22q_Subclone'],g[g=='chr22q_Clone'],g[g=='TERT'],g[!g %in% c(fusion_genes,'chr22q_Subclone','chr22q_Clone','TERT')],g[g %in% fusion_genes]))
      samples <- wgsdata %>% filter(Study2==cohort_name) %>% pull(Tumor_Barcode) %>% unique()
      top <- data_top0 %>% filter(Tumor_Barcode %in% samples) %>% select(-Study)
      burden <- data_tmb0 %>% filter(Tumor_Barcode %in% samples) %>% select(-Study)
      top_plot <- oncoplot(top,sample_level0=samples,landscape_colors=landscape_colors,GeneSortOnly=TRUE,sample_name=FALSE,cell_height=1,p2_hidden=TRUE,gene_level=gene_order[gene_order %in% top$Gene])
      burden_plot <- oncoplot3(burden,landscape_colors=landscape_colors,sample_level0=samples,sample_level=top_plot$sample_level$Tumor_Barcode,tmar=.2,bmar=.05,height=6)
      p <- oncoplot_combined(burden_plot,top_plot)
      list(paths=ptc_save(p,output_dir,paste0('Fig1',panel),12,8),n_tumors=length(samples),scope='Curated events in a shared gene order, with an explicit cohort selection.')
    })
  }
  results
}
