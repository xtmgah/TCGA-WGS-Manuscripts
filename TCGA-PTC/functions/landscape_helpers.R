# Only the three oncoplot functions used by Fig1g-i are retained.
# Sources: oncoplot_revised.R:1-271; Oncoplots_functions.R:405-494,665-end.
# Scientific ordering, split tiles, summaries and composition are preserved.
# Cleaning removes disabled/debug comments and makes legend lookup explicit.
# Typography increased to >=11pt; Roboto Condensed inherited from original theme.
# Caller libraries: tidyverse, data.table, cowplot, hrbrthemes, scales.
oncoplot <- function(data,data_clone = NULL, Gene_Sig=NULL,landscape_colors=NULL,gene_level=NULL,sample_level=NULL,sample_level0=sample_level0,GeneSortOnly=FALSE,namemax=NULL,tmar=0,bmar=0,nbreaks=NULL,p2_axis_hidden=FALSE,p2_hidden=FALSE,cell_height=0.95,legend_adjust=FALSE,removeAlterbyColor=FALSE,sample_name=TRUE){

  data <- data %>% filter(Tumor_Barcode %in% sample_level0)

  if(removeAlterbyColor & !is.null(landscape_colors)){
    data <- data %>% filter(Alteration %in%  names(landscape_colors[!is.na(landscape_colors)]))
  }


  blankgene <- gene_level[!(gene_level %in% data$Gene)]
  if(length(blankgene)>0){
    n = length(blankgene)
    blankdata <- tibble(Subject=rep(data$Subject[1],n),Tumor_Barcode=rep(data$Tumor_Barcode[1],n),Gene=blankgene,Alteration=rep(NA_character_,n),Type=rep(data$Type[1],n))
    data <- bind_rows(data,blankdata)
  }

  samplesize <- length(sample_level0)
  altertype <- data$Type[1]

  data0 <- data %>%
    group_by(Subject,Tumor_Barcode,Gene,Type) %>%
    summarise(Alteration_info=paste0(Alteration,collapse = ',')) %>%
    ungroup() %>%
    mutate(Alteration=if_else(Alteration_info=="NA",NA_character_,if_else(str_detect(Alteration_info,","),'Multi_Hit',Alteration_info))) %>%
    mutate(Alteration_info=if_else(str_detect(Alteration_info,","),Alteration_info,'')) %>%
    select(Subject,Tumor_Barcode, Gene, Type,Alteration, Alteration_info )

  data <- data %>%
    group_by(Subject,Tumor_Barcode,Gene,Type) %>%
    summarise(Alteration_info=paste0(Alteration,collapse = '/')) %>%
    ungroup() %>%
    mutate(Alteration=if_else(Alteration_info=="NA",NA_character_,Alteration_info)) %>%
    select(Subject,Tumor_Barcode, Gene, Type,Alteration, Alteration_info )

  data <- as.data.table(data)
  data <- data[
    , .(Alteration = unlist(strsplit(Alteration, "/", fixed = TRUE))),
    by = .(Subject, Tumor_Barcode, Gene, Type, Alteration_info)
  ]
  data[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(Tumor_Barcode, Gene)]
  data[, height:=cell_height/.N, by=list(Tumor_Barcode, Gene)]

  data <- as_tibble(data)%>% mutate(shift=shift*cell_height) # %>% rename(Alteration=V1)


  Gene_Freq <- data0 %>% count(Gene,sort=T) %>% mutate(n=if_else(Gene %in% blankgene,0L,n)) %>% mutate(Freq=n/samplesize) %>% arrange(Freq) %>%  mutate(Freq=percent(Freq,accuracy = 0.1))
  if(is.null(gene_level)){
    gene_level <- Gene_Freq %>% pull(Gene)
  }else{
    Gene_Freq <- Gene_Freq %>% mutate(Gene=factor(Gene,levels = gene_level)) %>% arrange(Gene)
  }

  Gene_Freq_alte <- data0 %>% count(Gene,Alteration) %>%  mutate(Gene=factor(Gene,levels = gene_level)) %>% arrange(Gene)

  input_sample_level=TRUE

  if(is.null(sample_level)){
    input_sample_level=FALSE
    sampleorder <- data0 %>%
      mutate(Gene=factor(Gene,levels = rev(gene_level))) %>%
      select(Tumor_Barcode,Gene,Alteration) %>%
      pivot_wider(id_cols = Tumor_Barcode,names_from = 'Gene',values_from = 'Alteration') %>%
      arrange_at(vars(one_of(rev(gene_level)))) %>%
      pull(Tumor_Barcode) %>%
      as.character()

    sample_level <- c(sampleorder,sample_level0[!c(sample_level0 %in% sampleorder)])
  }

  databg <- crossing(Tumor_Barcode=sample_level,Gene=gene_level) %>%
    mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>%
    mutate(Gene=factor(Gene,levels = gene_level))

  data <- data %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>%
    mutate(Gene=factor(Gene,levels = gene_level))

  data0 <- data0 %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>%
    mutate(Gene=factor(Gene,levels = gene_level))

  if(!is.null(namemax)){
    gene_level2 <- str_pad(gene_level,namemax,pad = " ")
  }else{
    gene_level2 <- gene_level
  }



  Gene_Freq <- Gene_Freq %>% mutate(Freq=if_else(Freq=='100.0%','',Freq))

  p1 <- databg %>%
    ggplot(aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene))) +
    geom_tile(height=cell_height,fill="gray95")+
    theme_ipsum_rc(base_size = 13,axis_text_size = 11,grid = '',axis = '')+
    theme(axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(face = "italic"),panel.background = element_blank(),legend.position = 'top',legend.justification = "top",legend.direction = 'vertical')+
    scale_y_continuous(expand = c(0,0),breaks = 1:length(gene_level2),labels = gene_level2,sec.axis = dup_axis(breaks = 1:length(gene_level2),labels = Gene_Freq$Freq))+
    geom_tile(data = data,aes(factor(Tumor_Barcode,levels = sample_level),as.integer(Gene)+shift,fill=Alteration,height=height),size=0)+
    scale_fill_manual(values = landscape_colors,breaks = sort(unique(data$Alteration)))+
    geom_vline(xintercept = 1:samplesize-0.5,color="white",size=0.2)+
    panel_border(size = 0.3,color = 'gray70')+
    guides(fill = guide_legend(ncol = 1,title.position = "top",title=altertype))

  if(legend_adjust){
    p1 <- p1+guides(fill = guide_legend(nrow = 1,title.position = 'left',title.hjust = 0.5,title.vjust = 0.5))+theme(legend.position = "bottom")
  }

  if(sample_name){
    p1 <- p1 + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = 11))
  }else{
    p1 <- p1 + theme(axis.text.x = element_blank())
  }

  if(!is.null(data_clone)){
    data_clone <- data_clone %>% filter(Tumor_Barcode %in% sample_level0)
    data_clone <- data_clone %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>%
      mutate(Gene=factor(Gene,levels = gene_level))
    p1 <- p1+geom_point(data=data_clone,pch=16,col="gray20",size=0.7)
  }

  oncoplot_legend <- cowplot::get_legend(p1+theme(plot.margin=margin(b=-12,l=0,r=0,unit="cm"),legend.key.size = unit(0.18, "cm")))
  p1 <- p1+theme(legend.position="none")


  if(p2_hidden){
    p2 <- NULL
  }else{
    p2 <- Gene_Freq_alte %>%
      ggplot(aes(Gene,y=n,fill=fct_rev(Alteration)))+geom_bar(stat="identity",width = 0.5,size=0)+
      scale_x_discrete(expand = expand_scale(add=c(0.5,0.5)))+
      coord_flip()+
      scale_fill_manual(values = landscape_colors,breaks = sort(unique(Gene_Freq_alte$Alteration)))+
      theme(legend.position = 'none',panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.text.x=element_text(size = 11),axis.line.x = element_line(size = 0.2),axis.ticks.x=element_line(size=0.2))

    if(is.null(nbreaks)){
      freq_max <- Gene_Freq_alte %>% group_by(Gene) %>% summarise(t=sum(n)) %>% arrange(desc(t)) %>% slice(1) %>% pull(t)
      freq_max <- ceiling(freq_max/5)*5
      p2 <- p2 + scale_y_continuous(breaks=pretty_breaks(),expand = c(0,0),position = 'right',limits = c(0,freq_max))
    }else{
      p2 <- p2 + scale_y_continuous(breaks=nbreaks,expand = c(0,0),position = 'right',limits = c(0,max(nbreaks)))
    }


    if(p2_axis_hidden){
      p2 <- p2+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank())
    }


  }




  if(!is.null(Gene_Sig)){
    p3 <- Gene_Sig %>%
      ggplot(aes(Gene,y=pvalue))+geom_bar(stat="identity",width = 0.5,size=0,fill="gray60")+
      scale_x_discrete(expand = expand_scale(add=c(0.5,0.5)),position = 'top')+
      scale_y_reverse(breaks=pretty_breaks(n = 2),expand = c(0,0),position = 'right')+
      coord_flip()+
      theme(legend.position = 'none',panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.text.x=element_text(size = 11),axis.line.x = element_line(size = 0.2),axis.ticks.x=element_line(size=0.2))
  }

  if(!is.null(Gene_Sig)){
    p_com <- align_plots(
      p3+theme(plot.margin=margin(l=0.2,r=0,t=tmar,b=bmar,unit="cm")),
      p1+theme(plot.margin=margin(r=0,t=tmar,b=bmar,unit="cm")),
      p2+theme(plot.margin=margin(l=-0.1,r=0.1,t=tmar,b=bmar,unit="cm")),
      align = 'h',
      axis = 'tb'
    )

  } else {
    p_com <- align_plots(p1+theme(plot.margin=margin(r=0,l=0.2,t=tmar,b=bmar,unit="cm")),
                         p2+theme(plot.margin=margin(l=-0.1,r=0.1,t=tmar,b=bmar,unit="cm")),
                         align = 'h',
                         axis = 'tb')
  }

  data <- data0
  if(!input_sample_level){

    if(GeneSortOnly){
      rankinfo <- data %>%
        select(Gene) %>%
        unique() %>%
        arrange(desc(Gene)) %>%
        mutate(Order=seq_along(Gene))
    }else{
      rankinfo <- data %>%
        select(Gene,Alteration) %>%
        unique() %>%
        arrange(desc(Gene),Alteration) %>%
        mutate(Order=seq_along(Gene))
    }

    maxrnak <- max(rankinfo$Order)+1

    rankinfo <- tibble(Tumor_Barcode=sample_level) %>%
      left_join(
        data %>%
          select(Tumor_Barcode,Gene,Alteration) %>%
          left_join(rankinfo) %>%
          group_by(Tumor_Barcode) %>%
          arrange(Order) %>%
          slice(1) %>%
          ungroup()
      ) %>%
      mutate(Order=if_else(is.na(Order),as.integer(maxrnak),as.integer(Order)))

    rankinfo <- rankinfo %>% select(Tumor_Barcode,Order)
    colnames(rankinfo)[2] <- altertype
  } else{

    rankinfo <- tibble(Tumor_Barcode=sample_level) %>%
      mutate(Order=as.integer(seq_along(Tumor_Barcode)))
    colnames(rankinfo)[2] <- altertype
  }

  return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level))

}

oncoplot3 <- function(data,sample_level0,Gene_Sig=NULL, landscape_colors=NULL,gene_level=NULL,sample_level=NULL,bmar=0,tmar=0,height=10){
  data <- data %>% filter(Tumor_Barcode %in% sample_level0)

  samplesize <- length(sample_level0)
  altertype <- data$Type[1]

  input_sample_level=TRUE

  if(is.null(sample_level)){
    input_sample_level=FALSE
    sampleorder <- data %>%
      select(Tumor_Barcode,Gene,Alteration) %>%
      arrange(desc(Alteration)) %>%
      pull(Tumor_Barcode)

    sample_level <- c(sampleorder,sample_level0[!c(sample_level0 %in% sampleorder)])
  }


  if(is.null(gene_level)){
    gene_level <- data %>% group_by(Gene) %>% summarise(mean=mean(Alteration)) %>% arrange(mean) %>% pull(Gene)
  }


  data <- data %>% mutate(Tumor_Barcode=factor(Tumor_Barcode,levels = sample_level)) %>%
    mutate(Gene=factor(Gene,levels = gene_level))


  p1 <- data %>% ggplot(aes(Tumor_Barcode,Alteration,fill=Gene))+
    geom_bar(position="stack", stat="identity",size=0)+
    scale_fill_manual(values = landscape_colors, breaks = sort(unique(data$Gene)))+
    theme_ipsum_rc(base_size = 13,axis_text_size = 11,grid = 'Y',axis = 'XY',axis_col = 'black')+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.title.y.right = element_blank(),axis.text.x = element_blank(),panel.background = element_blank(),legend.position = 'top',legend.justification = "top",legend.direction = 'vertical',legend.key.height =unit(1,"line"),axis.ticks.y = element_line(color = 'black') )+
    scale_y_continuous(expand = c(0,0),breaks = pretty_breaks(n = 3),labels = comma,sec.axis = dup_axis())+
    labs(fill=altertype)

  oncoplot_legend <- cowplot::get_legend(p1+theme(plot.margin=margin(b=-12,l=0,r=0,unit="cm"),legend.key.size = unit(0.18, "cm")))
  p1 <- p1+theme(legend.position="none")

  p2 <- ggplot()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),panel.background = element_blank())

  if(!is.null(Gene_Sig)){
    p3 <- ggplot()+theme(panel.background = element_blank())+scale_y_discrete(expand = expand_scale(add=c(0.5,0.5)))+scale_x_continuous(expand = c(0,0))
  }

  if(!is.null(Gene_Sig)){
    p_com <- align_plots(
      p3+theme(plot.margin=margin(l=0.2,r=0,t=tmar,b=bmar,unit="cm")),
      p1+theme(plot.margin=margin(r=0,b=bmar,t=tmar,unit="cm")),
      p2+theme(plot.margin=margin(l=-0.1,r=0.1,t=tmar,b=bmar,unit="cm")),
      align = 'h',
      axis = 'tb'
    )
  } else {
    p_com <- align_plots(p1+theme(plot.margin=margin(r=0,l=0.2,b=bmar,t=tmar,unit="cm")),
                         p2+theme(plot.margin=margin(l=-0.1,r=0.1,b=bmar,t=tmar,unit="cm")),
                         align = 'h',
                         axis = 'tb')
  }


  if(!input_sample_level){
    rankinfo <- tibble(Tumor_Barcode=sample_level) %>%
      left_join(
        data %>%
          select(Tumor_Barcode,Alteration)
      ) %>%
      arrange(desc(Alteration)) %>%
      mutate(Order=as.integer(seq_along(Tumor_Barcode)))

    rankinfo <- rankinfo %>% select(Tumor_Barcode,Order)
    colnames(rankinfo)[2] <- altertype
  } else{

    rankinfo <- tibble(Tumor_Barcode=sample_level) %>%
      mutate(Order=as.integer(seq_along(Tumor_Barcode)))
    colnames(rankinfo)[2] <- altertype
  }

  return(list(oncoplot=p_com,oncoplot_legend=oncoplot_legend,sample_level=rankinfo,gene_level=gene_level,height=height))

}

oncoplot_combined <- function(...,oncolist=NULL){
  oncolist <- c(list(...), oncolist)
  num_plots <- length(oncolist)
  plotlist1 <- list()
  plotlist2 <- list()
  leglist <- list()
  sizelist <- vector(length = num_plots)
  for(i in 1:num_plots){
    plotlist1[[i]] <- (oncolist[[i]]$oncoplot)[[1]]
    plotlist2[[i]] <- (oncolist[[i]]$oncoplot)[[2]]
    leglist[[i]] <- oncolist[[i]]$oncoplot_legend
    if(is.null(oncolist[[i]]$height)){
      sizelist[i] <- length(oncolist[[i]]$gene_level)
    }else {
      sizelist[i] <- oncolist[[i]]$height
    }

  }

  oncoplot_com1 <- plot_grid(plotlist = plotlist1,
                             align = 'v',
                             axis = 'lr',
                             rel_widths = rep(1,num_plots),
                             rel_heights = sizelist,
                             ncol= 1)
  oncoplot_com2 <- plot_grid(plotlist = plotlist2,
                             align = 'v',
                             axis = 'lr',
                             rel_widths = rep(1,num_plots),
                             rel_heights = sizelist,
                             ncol= 1)



  oncoplot_com <- plot_grid(oncoplot_com1,
                            oncoplot_com2,
                            align = 'h',
                            axis = 'tb',
                            rel_widths = c(10,1),
                            rel_heights = c(1,1),
                            nrow = 1)
  oncoplot_legend <- plot_grid(plotlist = c(leglist,list(NULL)),
                               align = 'h',
                               rel_widths = c(rep(1,num_plots),1),
                               nrow = 1)
  oncoplot_final <- plot_grid(
    oncoplot_com,
    oncoplot_legend,rel_heights = c(1,0.2),
    align = 'v',axis = 'l',ncol = 1)
  return(oncoplot_com)
}
