set_wd()
libztw()

load('wgs_data.RData')
load('~/NIH-Work/R/annotable_ens90.RData')

# 7q11.1/EGFR-------------------------------------------------------------------------
#samplelist <- sherlock_data %>% left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group)) %>% filter(Gene=='7q11.1',Alteration == 'Amplification') %>% pull(Tumor_Barcode)

for(barcode in wgs_groups_info$Tumor_Barcode){
  SCNA_Chr(barcode = barcode,chrq = '7',gene = 'EGFR',TMP = 'EGFR/')
}

# MDM2 ------------------------------------------------------------------------
for(barcode in wgs_groups_info$Tumor_Barcode){
  SCNA_Chr(barcode = barcode,chrq = '12',gene = 'MDM2',TMP = 'MDM2/')
}

# MDM2 ecDNA ------------------------------------------------------------------------
for(barcode in tdata$Tumor_Barcode){
  SCNA_Chr2(barcode = barcode,chrq = '12',gene = 'MDM2',TMP = 'MDM2_ecDNA/')
}


# CDKN2A ------------------------------------------------------------------------
for(barcode in wgs_groups_info$Tumor_Barcode){
  SCNA_Chr(barcode = barcode,chrq = '9',gene = 'CDKN2A',TMP = 'CDKN2A/')
}



## functions
barcode = 'TCGA-2G-AAER-01A-11D-A734-36'
chrq = 'X'
gene = c('SPRY3')
TMP = 'SPRY3/'
startpos0 <- 150e06
endpos0 <- 157e06


SCNA_Chr_Pos(barcode = 'TCGA-2G-AAER-01A-11D-A734-36',chrq = 'X',gene = c('SPRY3'),startpos0=154e06,endpos0=156e06,TMP = 'SPRY3/')
SCNA_Chr_Pos(barcode = 'TCGA-2G-AAEZ-01A-11D-A734-36',chrq = 'X',gene = c('SPRY3'),startpos0=154e06,endpos0=156e06,TMP = 'SPRY3/')
SCNA_Chr_Pos(barcode = 'TCGA-2G-AAEZ-01A-11D-A734-36',chrq = 'X',gene = c('SPRY3'),startpos0=154e06,endpos0=156e06,TMP = 'SPRY3/')



barcode = 'TCGA-2G-AAF7-01A-11D-A734-36'
chrq = '12'
gene = c('KRAS')
TMP = 'KRAS/'
startpos0 <- 150e06
endpos0 <- 157e06

SCNA_Chr_Pos(barcode = 'TCGA-2G-AAF7-01A-11D-A734-36',chrq = '12',gene = c('KRAS','MDM2'),startpos0=24e06,endpos0=26e06,TMP = 'KRAS/')




SCNA_Chr_Pos <- function(barcode,chrq,startpos0,endpos0,gene,ofilename=NULL,TMP='./'){
  
  NGSpurity_folder='/Volumes/Zhang_Group/ZTW/TCGA-TGCT/NGSpurity/Result/'
  bbprofile <- 'Battenberg_0'
  
  if(is.null(ofilename)){
    ofilename <-  paste0(TMP,barcode,'_',chrq,'_',gene,'.pdf')
  }
  
  #if(file.exists(ofilename)){ return(NULL)}
  
  genedata <- grch38 %>% filter(symbol %in% gene) %>% select(symbol,Chromosome = chr,start,end)
  
  #bbprofile <- BBsolution4 %>% filter(Tumor_Barcode == barcode) %>% pull(Battenberg)
  #bbprofile <- 'Battenberg_0'
  
  
  
  BB_SUBCLONE <- read_delim(file = paste0(NGSpurity_folder,barcode,'/',bbprofile,'/',barcode,'_subclones.txt'))
  BB_SUBCLONE <- BB_SUBCLONE %>%
    mutate(frac1_A=if_else(is.na(frac1_A),0,frac1_A)) %>%
    mutate(frac2_A=if_else(is.na(frac2_A),0,frac2_A)) %>%
    mutate(
      clone_frac=if_else(frac1_A>frac2_A,frac1_A,frac2_A),
      clone_nMaj=if_else(frac1_A>frac2_A,nMaj1_A,nMaj2_A),
      clone_nMin=if_else(frac1_A>frac2_A,nMin1_A,nMin2_A)
    ) %>% 
    mutate(
      subclone_frac=if_else(frac1_A>frac2_A,frac2_A,frac1_A),
      subclone_nMaj=if_else(frac1_A>frac2_A,nMaj2_A,nMaj1_A),
      subclone_nMin=if_else(frac1_A>frac2_A,nMin2_A,nMin1_A)
    )
  
  
  BB_SUBCLONE2 <- BB_SUBCLONE %>%
    select(Chromosome=chr,startpos,endpos,BAF,LogR,contains('clone')) %>%
    filter(Chromosome == chrq) # startpos>=startpos0,endpos<=endpos0
  
  BB_SUBCLONE2 <- BB_SUBCLONE2 %>%
    mutate(startpos =if_else(startpos<startpos0,startpos0,startpos),endpos =if_else(endpos>endpos0,endpos0,endpos))
  
  datafile <- paste0(NGSpurity_folder,barcode,'/',bbprofile,'/NGSpurity/',barcode,'_mutant_LogRgc_BAF_nsample.txt')
  
  # if(!file.exists(datafile)){
  #   datafile <- paste0('/Volumes/data/NSLC2/NGSpurity2/Result_pub/',barcode,'/',bbprofile,'/NGSpurity/',barcode,'_mutant_LogRgc_BAF_nsample.txt')
  # }
  # 
  
  jdata <- read_delim(datafile,delim = '\t',col_names = F,col_types = cols('X1'='c'))
  colnames(jdata) <- c("Chromosome","Position","LogR","BAF")
  # jdata <-jdata %>%  mutate(Chromosome=factor(Chromosome,levels = chrs))
  jdata <- jdata %>% filter(Chromosome == chrq,Position >= startpos0, Position <= endpos0)
  
  chrs=c(1:22,"X")
  
  chrs_data <- bind_rows(
    jdata %>% group_by(Chromosome) %>% dplyr::arrange(Position) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=Inf,BAF=Inf) %>% dplyr::filter(Chromosome!=chrs[1]),
    jdata %>% group_by(Chromosome) %>% dplyr::arrange(desc(Position)) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=-Inf,BAF=-Inf) %>% dplyr::filter(Chromosome!=chrs[length(chrs)]),
  )
  
  chrs_data2 <- bind_rows(
    jdata %>% group_by(Chromosome) %>% arrange(Position) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=Inf,BAF=Inf) %>% dplyr::filter(Chromosome==chrs[1]),
    jdata %>% group_by(Chromosome) %>% arrange(desc(Position)) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=-Inf,BAF=-Inf) %>% dplyr::filter(Chromosome==chrs[length(chrs)]),
  )
  
  chrs_data3 <- left_join(
    jdata %>% group_by(Chromosome) %>% arrange(Position) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=NA_real_,BAF=NA_real_) %>% dplyr::rename(startpos=Position),
    jdata %>% group_by(Chromosome) %>% arrange(desc(Position)) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=NA_real_,BAF=NA_real_) %>% dplyr::rename(endpos=Position)
  )
  
  
  
  # logR
  p1 <- jdata %>% 
    ggplot(aes(Position,LogR))+#facet_grid(.~Chromosome,space = "free_x",scales = "free")+
    #geom_bin2d(bins = 200) +
    #scale_fill_viridis_c()+
    geom_point(pch=21,fill=ncicolpal[8],col='white',size=2,stroke=0.1)+
    labs(x="",y="LogR")+
    scale_x_continuous(expand = c(0,0),limits = c(startpos0,endpos0))+
    scale_y_continuous(expand = c(0,0),limits = c(floor(min(jdata$LogR)),ceiling(max(jdata$LogR))))+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = FALSE, grid = "Y",ticks = TRUE)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_text(hjust = 0.5),strip.background = element_blank())+
    geom_vline(data=chrs_data ,aes(xintercept = Position),col="#cccccc",size=0.5)+
    geom_vline(data=chrs_data2,aes(xintercept = Position),col="black",size=c(0.8))+
    geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)+
    geom_hline(yintercept = 0,col="#cccccc",size=0.5)+
    geom_vline(xintercept = c(genedata$start,genedata$end),linetype=2,col='#f03b20')+
    ggrepel::geom_text_repel(data = genedata,aes(x = (start+end)/2,y=max(jdata$LogR),label=symbol,col='#f03b20'),
                             force_pull   = 0,
                             nudge_y       = 0,
                             segment.size  = 0.2,
                             segment.color = "grey40",
                             direction     = "x",
                             angle=90,
                             size=3
    )
  
  ## add logr segement
  p1 <- p1+geom_segment(data=BB_SUBCLONE2,aes(x=startpos,xend=endpos,y=LogR, yend=LogR),col="#e7298a",size=1.2)+coord_cartesian(clip = 'off')
  
  
  # BAF plot
  p2 <- jdata %>%
    ggplot(aes(Position,BAF))+#facet_grid(.~Chromosome,space = "free_x",scales = "free")+
    # geom_bin2d(bins = 200) +
    # scale_fill_viridis_c()+
    geom_point(pch=21,fill=ncicolpal[8],col='white',size=2,stroke=0.1)+
    labs(x="",y="BAF")+
    scale_x_continuous(expand = c(0,0),limits = c(startpos0,endpos0))+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = FALSE, grid = "Y",ticks = TRUE)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_blank())+
    geom_vline(data=chrs_data ,aes(xintercept = Position),col="#cccccc",size=0.5)+
    geom_vline(data=chrs_data2,aes(xintercept = Position),col="black",size=c(0.8))+
    geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)+
    geom_hline(yintercept = 0,col="white",size=0.5)+
    geom_vline(xintercept = c(genedata$start,genedata$end),linetype=2,col='#f03b20')
  
  p2 <- p2+geom_segment(data=BB_SUBCLONE2,aes(x=startpos,xend=endpos,y=BAF, yend=BAF),col="#e7298a",size=1.2)+
    geom_segment(data=BB_SUBCLONE2 %>% mutate(BAF=1-BAF),aes(x=startpos,xend=endpos,y=BAF, yend=BAF),col="#e7298a",size=1.2)
  
  
  # Copy number
  nmax <- max(c(BB_SUBCLONE2$clone_nMaj+BB_SUBCLONE2$clone_nMin),c(BB_SUBCLONE2$clone_nMaj+BB_SUBCLONE2$clone_nMin),na.rm = TRUE)
  nmax_flag <- if_else(nmax>5,TRUE,FALSE)
  nmax <- if_else(nmax>5,5L,nmax)
  nmax <- if_else(nmax<3,3L,nmax)
  swidth=(0.1*nmax)/3
  
  swidth=0.1
  data_p3 <- BB_SUBCLONE2 %>%
    bind_rows(chrs_data3) %>%
    mutate(clone_nTot=clone_nMaj+clone_nMin,subclone_nTot=subclone_nMaj+subclone_nMin) %>%
    mutate(clone_nTot=if_else(clone_nTot>nmax,nmax,clone_nTot),subclone_nTot=if_else(subclone_nTot>nmax,nmax,subclone_nTot)) %>%
    mutate(clone_nTot_adjust=if_else(subclone_frac!=0 & (clone_nTot==subclone_nTot|clone_nTot==subclone_nMin),subclone_frac*swidth,0)) %>%
    mutate(clone_nMin_adjust=if_else(subclone_frac!=0 & (clone_nMin==subclone_nTot|clone_nMin==subclone_nMin),subclone_frac*swidth,0)) %>%
    mutate(subclone_nTot_adjust=if_else(subclone_frac!=0 & (subclone_nTot==clone_nTot|subclone_nTot==clone_nMin),clone_frac*swidth,0)) %>%
    mutate(subclone_nMin_adjust=if_else(subclone_frac!=0 & (subclone_nMin==clone_nTot|subclone_nMin==clone_nMin),clone_frac*swidth,0))
  
  data_p3_tmp <- data_p3 %>% mutate(ymin=clone_nMin-clone_frac*swidth+clone_nMin_adjust,ymin2=subclone_nMin-subclone_frac*swidth-subclone_nMin_adjust,ymax=clone_nTot+clone_frac*swidth+clone_nTot_adjust,ymax2=subclone_nTot+subclone_frac*swidth-subclone_nTot_adjust)
  nmin <- min(c(data_p3_tmp$ymin,data_p3_tmp$ymin2),na.rm = TRUE)
  nmax2 <- max(c(data_p3_tmp$ymax,data_p3_tmp$ymax2),na.rm = TRUE)+0.1
  nmin <- if_else(nmin>=0,0,nmin-0.1)
  
  
  # p3 <- data_p3 %>%
  #   ggplot()+facet_grid(.~Chromosome,space = "free_x",scales = "free")+
  #   geom_segment(aes(x=startpos,xend=endpos,y=clone_nMaj+clone_nMin, yend=clone_nMaj+clone_nMin),col="#FF7F0EFF",size=4)+
  #   geom_segment(aes(x=startpos,xend=endpos,y=clone_nMin, yend=clone_nMin),col="#006d2c",size=4)+
  #   labs(x="",y="Battenberg Copy Number\n")+
  #   scale_x_continuous(expand = c(0,0))+
  #   ylim(c(0,nmax))+
  #   theme_ipsum_rc(axis_title_just = "m",axis_title_face = "bold",axis_title_size = 12,axis = FALSE, grid = "Y",ticks = TRUE)+
  #   theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_blank())+
  #   geom_vline(data=chrs_data ,aes(xintercept = Position),col="#cccccc",size=0.5)+
  #   geom_vline(data=chrs_data2 ,aes(xintercept = Position),col="black",size=0.8)+
  #   geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)
  #pay attention to the order for the overlap of homozygous deletion
  
  p3 <- data_p3 %>%
    ggplot()+#facet_grid(.~Chromosome,space = "free_x",scales = "free")+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=subclone_nMin-subclone_frac*swidth-subclone_nMin_adjust, ymax=subclone_nMin+subclone_frac*swidth-subclone_nMin_adjust),fill="#2171b5",size=0)+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=subclone_nTot-subclone_frac*swidth-subclone_nTot_adjust, ymax=subclone_nTot+subclone_frac*swidth-subclone_nTot_adjust),fill="#c51b8a",size=0)+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=clone_nMin-clone_frac*swidth+clone_nMin_adjust, ymax=clone_nMin+clone_frac*swidth+clone_nMin_adjust),fill="#006d2c",size=0)+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=clone_nTot-clone_frac*swidth+clone_nTot_adjust, ymax=clone_nTot+clone_frac*swidth+clone_nTot_adjust),fill="#FF7F0EFF",size=0)+
    labs(x="",y="Copy Number")+
    scale_x_continuous(expand = c(0,0),limits = c(startpos0,endpos0))+
    scale_y_continuous(breaks = seq(0,nmax),labels = c(seq(0,nmax-1),paste0("\u2265",nmax)),limits = c(nmin,nmax2))+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = FALSE, grid = "Y",ticks = TRUE)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_blank())+
    geom_vline(data=chrs_data ,aes(xintercept = Position),col="#cccccc",size=0.5)+
    geom_vline(data=chrs_data2 ,aes(xintercept = Position),col="black",size=0.8)+
    geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)+
    geom_vline(xintercept = c(genedata$start,genedata$end),linetype=2,col='#f03b20')
  
  
  p4 <- BB_SUBCLONE2 %>%
    bind_rows(chrs_data3) %>%
    bind_rows(chrs_data3 %>% dplyr::slice(1) %>% dplyr::mutate(endpos=startpos,clone_frac=0)) %>%
    ggplot(aes(xmin=startpos/1e6,xmax=endpos/1e6,ymin=0, ymax=1,fill=clone_frac))+#facet_grid(.~Chromosome,space = "free_x",scales = "free")+
    geom_rect()+
    labs(x=paste0('Chromosome ',chrq,' (Mb)'),y="CCF")+
    scale_fill_material(palette = "blue",na.value = "transparent",breaks=c(0,1))+
    scale_x_continuous(expand = c(0,0),breaks = pretty_breaks(n = 7),limits = c(startpos0/1e6,endpos0/1e6))+
    scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = FALSE, grid = FALSE,ticks = TRUE)+
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_blank())+
    geom_vline(data=chrs_data ,aes(xintercept = Position/1e6),col="#cccccc",size=0.5)+
    geom_vline(data=chrs_data2 ,aes(xintercept = Position/1e6),col="black",size=0.8)+
    geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)+
    geom_vline(xintercept = c(genedata$start/1e6,genedata$end/1e6),linetype=2,col='#f03b20')
  
  p_com <- plot_grid(p1+theme(plot.margin = margin(b=-5,r=1,l=1)),p2+theme(plot.margin = margin(t=-5,b=-5,r=1,l=1)),p3+theme(plot.margin = margin(t=-5,r=1,l=1)),p4+theme(plot.margin = margin(t=-10,b=5,r=1,l=1)),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(1,1,0.7,0.2))
  
  ggsave(filename = ofilename,plot = p_com,width = 8,height = 6,device = cairo_pdf)
  
}







SCNA_Chr <- function(barcode,chrq,gene,ofilename=NULL,TMP='./'){
  
  if(is.null(ofilename)){
    ofilename <-  paste0(TMP,barcode,'_',chrq,'_',gene,'.pdf')
  }
  
  if(file.exists(ofilename)){ return(NULL)}
  
  genedata <- grch38 %>% filter(symbol==gene) %>% select(Chromosome = chr,start,end)
  
  bbprofile <- BBsolution4 %>% filter(Tumor_Barcode == barcode) %>% pull(Battenberg)
  
  BB_SUBCLONE2 <- BBprofile %>% filter(Tumor_Barcode == barcode) %>% 
    select(Chromosome=chr,startpos,endpos,BAF,LogR,contains('clone')) %>% 
    filter(Chromosome == chrq)
  
  datafile <- paste0('/Volumes/data/NSLC2/NGSpurity2/Result/',barcode,'/',bbprofile,'/NGSpurity/',barcode,'_mutant_LogRgc_BAF_nsample.txt')
  
  if(!file.exists(datafile)){
    datafile <- paste0('/Volumes/data/NSLC2/NGSpurity2/Result_pub/',barcode,'/',bbprofile,'/NGSpurity/',barcode,'_mutant_LogRgc_BAF_nsample.txt')
  }
  
  
  jdata <- read_delim(datafile,delim = '\t',col_names = F,col_types = cols('X1'='c'))
  colnames(jdata) <- c("Chromosome","Position","LogR","BAF")
  # jdata <-jdata %>%  mutate(Chromosome=factor(Chromosome,levels = chrs))
  jdata <- jdata %>% filter(Chromosome == chrq)
  
  chrs=c(1:22,"X")
  
  chrs_data <- bind_rows(
    jdata %>% group_by(Chromosome) %>% dplyr::arrange(Position) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=Inf,BAF=Inf) %>% dplyr::filter(Chromosome!=chrs[1]),
    jdata %>% group_by(Chromosome) %>% dplyr::arrange(desc(Position)) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=-Inf,BAF=-Inf) %>% dplyr::filter(Chromosome!=chrs[length(chrs)]),
  )
  
  chrs_data2 <- bind_rows(
    jdata %>% group_by(Chromosome) %>% arrange(Position) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=Inf,BAF=Inf) %>% dplyr::filter(Chromosome==chrs[1]),
    jdata %>% group_by(Chromosome) %>% arrange(desc(Position)) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=-Inf,BAF=-Inf) %>% dplyr::filter(Chromosome==chrs[length(chrs)]),
  )
  
  chrs_data3 <- left_join(
    jdata %>% group_by(Chromosome) %>% arrange(Position) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=NA_real_,BAF=NA_real_) %>% dplyr::rename(startpos=Position),
    jdata %>% group_by(Chromosome) %>% arrange(desc(Position)) %>% dplyr::slice(1) %>% ungroup() %>% mutate(LogR=NA_real_,BAF=NA_real_) %>% dplyr::rename(endpos=Position)
  )
  
  
  
  # logR
  p1 <- jdata %>% 
    ggplot(aes(Position,LogR))+#facet_grid(.~Chromosome,space = "free_x",scales = "free")+
    geom_bin2d(bins = 200) +
    scale_fill_viridis_c()+
    labs(x="",y="LogR")+
    scale_x_continuous(expand = c(0,0))+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = FALSE, grid = "Y",ticks = TRUE)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_text(hjust = 0.5),strip.background = element_blank())+
    geom_vline(data=chrs_data ,aes(xintercept = Position),col="#cccccc",size=0.5)+
    geom_vline(data=chrs_data2,aes(xintercept = Position),col="black",size=c(0.8,0.5))+
    geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)+
    geom_hline(yintercept = 0,col="#cccccc",size=0.5)+
    geom_vline(xintercept = c(genedata$start,genedata$end),linetype=2,col='#f03b20')
  
  ## add logr segement
  p1 <- p1+geom_segment(data=BB_SUBCLONE2,aes(x=startpos,xend=endpos,y=LogR, yend=LogR),col="#e7298a",size=1.2)
  
  
  # BAF plot
  p2 <- jdata %>%
    ggplot(aes(Position,BAF))+#facet_grid(.~Chromosome,space = "free_x",scales = "free")+
    geom_bin2d(bins = 200) +
    #geom_point()+
    scale_fill_viridis_c()+
    labs(x="",y="BAF")+
    scale_x_continuous(expand = c(0,0))+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = FALSE, grid = "Y",ticks = TRUE)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_blank())+
    geom_vline(data=chrs_data ,aes(xintercept = Position),col="#cccccc",size=0.5)+
    geom_vline(data=chrs_data2,aes(xintercept = Position),col="black",size=c(0.8,0.5))+
    geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)+
    geom_hline(yintercept = 0,col="white",size=0.5)+
    geom_vline(xintercept = c(genedata$start,genedata$end),linetype=2,col='#f03b20')
  
  p2 <- p2+geom_segment(data=BB_SUBCLONE2,aes(x=startpos,xend=endpos,y=BAF, yend=BAF),col="#e7298a",size=1.2)+
    geom_segment(data=BB_SUBCLONE2 %>% mutate(BAF=1-BAF),aes(x=startpos,xend=endpos,y=BAF, yend=BAF),col="#e7298a",size=1.2)
  
  
  # Copy number
  nmax <- max(c(BB_SUBCLONE2$clone_nMaj+BB_SUBCLONE2$clone_nMin),c(BB_SUBCLONE2$clone_nMaj+BB_SUBCLONE2$clone_nMin),na.rm = TRUE)
  nmax_flag <- if_else(nmax>5,TRUE,FALSE)
  nmax <- if_else(nmax>5,5L,nmax)
  nmax <- if_else(nmax<3,3L,nmax)
  swidth=(0.1*nmax)/3
  
  swidth=0.1
  data_p3 <- BB_SUBCLONE2 %>%
    bind_rows(chrs_data3) %>%
    mutate(clone_nTot=clone_nMaj+clone_nMin,subclone_nTot=subclone_nMaj+subclone_nMin) %>%
    mutate(clone_nTot=if_else(clone_nTot>nmax,nmax,clone_nTot),subclone_nTot=if_else(subclone_nTot>nmax,nmax,subclone_nTot)) %>%
    mutate(clone_nTot_adjust=if_else(subclone_frac!=0 & (clone_nTot==subclone_nTot|clone_nTot==subclone_nMin),subclone_frac*swidth,0)) %>%
    mutate(clone_nMin_adjust=if_else(subclone_frac!=0 & (clone_nMin==subclone_nTot|clone_nMin==subclone_nMin),subclone_frac*swidth,0)) %>%
    mutate(subclone_nTot_adjust=if_else(subclone_frac!=0 & (subclone_nTot==clone_nTot|subclone_nTot==clone_nMin),clone_frac*swidth,0)) %>%
    mutate(subclone_nMin_adjust=if_else(subclone_frac!=0 & (subclone_nMin==clone_nTot|subclone_nMin==clone_nMin),clone_frac*swidth,0))
  
  data_p3_tmp <- data_p3 %>% mutate(ymin=clone_nMin-clone_frac*swidth+clone_nMin_adjust,ymin2=subclone_nMin-subclone_frac*swidth-subclone_nMin_adjust,ymax=clone_nTot+clone_frac*swidth+clone_nTot_adjust,ymax2=subclone_nTot+subclone_frac*swidth-subclone_nTot_adjust)
  nmin <- min(c(data_p3_tmp$ymin,data_p3_tmp$ymin2),na.rm = TRUE)
  nmax2 <- max(c(data_p3_tmp$ymax,data_p3_tmp$ymax2),na.rm = TRUE)+0.1
  nmin <- if_else(nmin>=0,0,nmin-0.1)
  
  
  # p3 <- data_p3 %>%
  #   ggplot()+facet_grid(.~Chromosome,space = "free_x",scales = "free")+
  #   geom_segment(aes(x=startpos,xend=endpos,y=clone_nMaj+clone_nMin, yend=clone_nMaj+clone_nMin),col="#FF7F0EFF",size=4)+
  #   geom_segment(aes(x=startpos,xend=endpos,y=clone_nMin, yend=clone_nMin),col="#006d2c",size=4)+
  #   labs(x="",y="Battenberg Copy Number\n")+
  #   scale_x_continuous(expand = c(0,0))+
  #   ylim(c(0,nmax))+
  #   theme_ipsum_rc(axis_title_just = "m",axis_title_face = "bold",axis_title_size = 12,axis = FALSE, grid = "Y",ticks = TRUE)+
  #   theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_blank())+
  #   geom_vline(data=chrs_data ,aes(xintercept = Position),col="#cccccc",size=0.5)+
  #   geom_vline(data=chrs_data2 ,aes(xintercept = Position),col="black",size=0.8)+
  #   geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)
  #pay attention to the order for the overlap of homozygous deletion
  
  p3 <- data_p3 %>%
    ggplot()+#facet_grid(.~Chromosome,space = "free_x",scales = "free")+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=subclone_nMin-subclone_frac*swidth-subclone_nMin_adjust, ymax=subclone_nMin+subclone_frac*swidth-subclone_nMin_adjust),fill="#2171b5",size=0)+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=subclone_nTot-subclone_frac*swidth-subclone_nTot_adjust, ymax=subclone_nTot+subclone_frac*swidth-subclone_nTot_adjust),fill="#c51b8a",size=0)+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=clone_nMin-clone_frac*swidth+clone_nMin_adjust, ymax=clone_nMin+clone_frac*swidth+clone_nMin_adjust),fill="#006d2c",size=0)+
    geom_rect(aes(xmin=startpos,xmax=endpos,ymin=clone_nTot-clone_frac*swidth+clone_nTot_adjust, ymax=clone_nTot+clone_frac*swidth+clone_nTot_adjust),fill="#FF7F0EFF",size=0)+
    labs(x="",y="Copy Number")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(breaks = seq(0,nmax),labels = c(seq(0,nmax-1),paste0("\u2265",nmax)),limits = c(nmin,nmax2))+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = FALSE, grid = "Y",ticks = TRUE)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_blank())+
    geom_vline(data=chrs_data ,aes(xintercept = Position),col="#cccccc",size=0.5)+
    geom_vline(data=chrs_data2 ,aes(xintercept = Position),col="black",size=0.8)+
    geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)+
    geom_vline(xintercept = c(genedata$start,genedata$end),linetype=2,col='#f03b20')
  
  
  p4 <- BB_SUBCLONE2 %>%
    bind_rows(chrs_data3) %>%
    bind_rows(chrs_data3 %>% dplyr::slice(1) %>% dplyr::mutate(endpos=startpos,clone_frac=0)) %>%
    ggplot(aes(xmin=startpos/1e6,xmax=endpos/1e6,ymin=0, ymax=1,fill=clone_frac))+#facet_grid(.~Chromosome,space = "free_x",scales = "free")+
    geom_rect()+
    labs(x=paste0('Chromosome ',chrq,' (Mb)'),y="CCF")+
    scale_fill_material(palette = "blue",na.value = "transparent",breaks=c(0,1))+
    scale_x_continuous(expand = c(0,0),breaks = pretty_breaks(n = 7))+
    scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = FALSE, grid = FALSE,ticks = TRUE)+
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.spacing = unit(0,"cm"),legend.position = "none",strip.text.x = element_blank())+
    geom_vline(data=chrs_data ,aes(xintercept = Position/1e6),col="#cccccc",size=0.5)+
    geom_vline(data=chrs_data2 ,aes(xintercept = Position/1e6),col="black",size=0.8)+
    geom_hline(yintercept = c(-Inf,Inf),col="black",size=0.8)+
    geom_vline(xintercept = c(genedata$start/1e6,genedata$end/1e6),linetype=2,col='#f03b20')
  
  p_com <- plot_grid(p1+theme(plot.margin = margin(b=-5,r=1,l=1)),p2+theme(plot.margin = margin(t=-5,b=-5,r=1,l=1)),p3+theme(plot.margin = margin(t=-5,r=1,l=1)),p4+theme(plot.margin = margin(t=-10,b=5,r=1,l=1)),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(1,1,0.7,0.2))
  
  ggsave(filename = ofilename,plot = p_com,width = 8,height = 6,device = cairo_pdf)
  
}

