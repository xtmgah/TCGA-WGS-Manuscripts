set_wd()
library(tidyverse)
library(tidylog)
library(CNTools)
library(ggsci)
library(scales)
library(hrbrthemes)
library(cowplot)
library(valr)

cyto <- read_delim('cytoBand.txt.gz',delim = '\t',col_names = F)
colnames(cyto) <- c('Chromosome', 'Start' , 'End', 'Band', 'gieStain')
cyto <- cyto %>% mutate(Chromosome = str_remove(Chromosome,"chr"),Centro=if_else(gieStain=="acen",1,0), Arm=paste0(Chromosome,str_sub(Band,1,1))) %>% filter(Chromosome %in% c(1:22,"X","Y"))
hg38centro <- cyto %>% filter(Centro!=0,Chromosome %in% c(1:22)) %>% select(chrom=Chromosome,start=Start,end=End,Centro) 


load('wgs_data.RData')
load('ngspurity_data.RData')


sample_order <- wgs_data %>% select(Tumor_Barcode) %>% mutate(Seq=seq_along(Tumor_Barcode)) 

# read subclone data  -----------------------------------------------------

cnv <- BBprofile %>% mutate(frac1_A=as.numeric(frac1_A),frac2_A=as.numeric(frac2_A),nMaj1_A=as.integer(nMaj1_A),nMaj2_A=as.integer(nMaj2_A), nMin1_A=as.integer(nMin1_A),nMin2_A=as.integer(nMin2_A))
#cnv %>% filter(chr %in% c(1:22,"X","Y"),Tumor_Barcode %in% sherlock_samples_unique$Tumor_Barcode) %>% select(Tumor_Barcode,everything()) %>% write_delim('bb_subclone_232.txt',delim = '\t',col_names = T)
cnv <- cnv %>% filter(chr %in% c(1:22,"X","Y"),Tumor_Barcode %in% sample_order$Tumor_Barcode)
cnv <- cnv %>% 
  mutate(frac1_A=if_else(is.na(frac1_A),0,frac1_A)) %>% 
  mutate(frac2_A=if_else(is.na(frac2_A),0,frac2_A)) %>% 
  mutate(
    clone_frac=if_else(frac1_A>frac2_A,frac1_A,frac2_A),
    clone_nMaj=if_else(frac1_A>frac2_A,nMaj1_A,nMaj2_A),
    clone_nMin=if_else(frac1_A>frac2_A,nMin1_A,nMin2_A),
    subclone_frac=if_else(frac1_A<frac2_A,frac1_A,frac2_A),
    subclone_nMaj=if_else(frac1_A<frac2_A,nMaj1_A,nMaj2_A),
    subclone_nMin=if_else(frac1_A<frac2_A,nMin1_A,nMin2_A)
  ) %>% 
  select(Tumor_Barcode,chr:ntot,contains('clone'))


# Define the value --------------------------------------------------------
totalsamples <- cnv %>% count(Tumor_Barcode) %>% pull(Tumor_Barcode) %>% length()
allsamples <- cnv %>% count(Tumor_Barcode) %>% pull(Tumor_Barcode) %>% unique()


cnvdata <- cnv %>%
  filter(chr %in% c(1:22,"X")) %>% 
  mutate(clone_total=clone_nMin+clone_nMaj) %>% 
  left_join(
    ngspurity %>% select(Tumor_Barcode,WGD_Status=MCN_WGD,Tumor_Purity=BB_Purity,Tumor_Ploidy=BB_Ploidy)
  ) %>% 
  mutate(WGD_Status=if_else(is.na(WGD_Status),'nWGD',WGD_Status)) %>% 
  mutate(relative_copy=clone_total-if_else(WGD_Status=="WGD",2*if_else(chr=='X',1,2),if_else(chr=='X',1,2))) %>% 
  mutate(relative_copy=if_else(relative_copy>4,4,relative_copy)) %>% 
  mutate(relative_copy=if_else(relative_copy< -4,-4,relative_copy)) %>% 
  mutate(relative_copy=if_else(clone_nMaj==0,-4,relative_copy)) %>% 
  mutate(startpos=as.integer(startpos),endpos=as.integer(endpos)) 


# Dendrograms -------------------------------------------------------------
source('~/NIH-Work/EAGLE_NSLC/FirstPaper/Biowulf/ggdendro.R')
segdata <- cnvdata %>% 
  mutate(num.mark=1000) %>% 
  select(ID=Tumor_Barcode,chrom=chr,loc.start=startpos,loc.end=endpos,num.mark,seg.mean=relative_copy) %>%
  mutate(seg.mean=if_else(is.na(seg.mean),0,seg.mean)) %>% 
  as.data.frame()
segdata %>% filter_all(any_vars(is.na(.)))
cnseg <- CNSeg(segdata)
rdseg <- getRS(cnseg, by = "region", imput = FALSE, XY = FALSE, what = "mean")
reducedseg <- rs(rdseg)

msegdata <- as.matrix(reducedseg[,-(1:3)])

msegdata <- apply(msegdata, 2, as.numeric)

### add samples without any change
if(dim(msegdata)[2]<totalsamples){
  numbr <- dim(msegdata)[1]
  numbc <- totalsamples-dim(msegdata)[2]
  cnvtmp <- matrix(rep(0,numbr*numbc),nrow = numbr,ncol =numbc )
  extrasamples <- allsamples[!(allsamples %in% colnames(msegdata))]
  colnames(cnvtmp) <- extrasamples
  msegdata <- cbind(msegdata,cnvtmp)
  
  cnvdata <- bind_rows(
    cnvdata,
    tibble(Tumor_Barcode=rep(extrasamples,each=22),chr=rep(1:22,numbc)) %>% mutate(chr=as.character(chr))
  )
  #change cnv data 
  
  
}


hc <- hclust(dist(t(msegdata),method = 'euclidean'),method = 'ward.D')
hcdata <- dendro_data_k(hc, 3)

cols <- c("#a9a9a9", "#2ca02c", "#ff7f0e","#1f77b4")

p_dend <- plot_ggdendro(hcdata,
                        direction   = "lr",
                        scale.color = cols,
                        label.size  = 2.5,
                        branch.size = 0.5,
                        expand.y    = 0,nudge.label = 0)

p_dend <- p_dend+
  scale_x_continuous(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))+
  theme_void()

sample_order <- as_tibble(hcdata$labels) %>% select(Tumor_Barcode=label) %>% unique() %>% mutate(Seq=seq_along(Tumor_Barcode)) 

#chr19loss <- sample_order %>% tail(n=11) %>% mutate(chr19loss="Y") %>% select(-Seq)
#chr19loss <- cnvdata %>% filter(chr=="19",clone_nMaj==0) %>% mutate(size=endpos-startpos) %>% filter(size>30000000) %>% mutate(chr19loss="Y") %>% select(Tumor_Barcode,chr19loss)
#save(chr19loss,file='chr19loss.RData')

cnvclust <- as_tibble(hcdata$labels) %>% select(Tumor_Barcode=label,CNV_Clust=clust) %>% mutate(CNV_Clust=paste0('C',CNV_Clust))
#cnvclust <- cnvclust %>% mutate(SCNA_Group=if_else(CNV_Clust == 'C1', 'Forte',if_else(CNV_Clust == 'C2', 'Piano', 'Mezzo-forte'))) 
#save(cnvclust,file='cnvclust_wgs_all.RData')
save(cnvclust,file='cnvclust_wgs_all.RData')

## C3 and WGD
library(broom)
cnvclust %>% 
  mutate(Group1=if_else(CNV_Clust == 'C3','Yes','No')) %>% 
  left_join(
    ngspurity %>% select(Tumor_Barcode,MCN_WGD) %>% mutate(Group2=if_else(MCN_WGD=='WGD','Yes','No'))
  ) %>% 
  do(tidy(fisher.test(.$Group1,.$Group2)))
  


# Heatmap plot ------------------------------------------------------------

chrlevels <- c(seq(1:22),"X")
p_heatmap <- cnvdata %>% 
  mutate(chr=factor(chr,levels = chrlevels)) %>%
  mutate(startpos=as.integer(startpos),endpos=as.integer(endpos)) %>% 
  left_join(sample_order) %>% 
  arrange(Seq,chr) %>% 
  #filter(Seq %in% c(1:20)) %>% 
  ggplot(aes(fill=relative_copy))+
  #geom_segment(aes(x=startpos,xend=endpos,y=Tumor_Barcode,yend=Tumor_Barcode),size=3)+
  geom_rect(aes(xmin=startpos,xmax=endpos,ymin=Seq-0.5,ymax=Seq+0.5))+
  facet_grid(~chr,scales = 'free_x',space = 'free',switch="both")+
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",midpoint = 0)+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
  labs(fill="Relative Copy Number")+
  scale_y_continuous(breaks = sample_order$Seq,labels = sample_order$Tumor_Barcode,expand = expand_scale(add = c(0, 0)), sec.axis = dup_axis())+
  #facet_wrap(~chr,nrow = 1,scales = 'free_x',strip.position = "bottom")+
  labs(x="",y="")+
  #theme_ipsum_rc()+
  theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0, "lines"),axis.text.x = element_blank(),axis.ticks = element_blank(),axis.text.y = element_text(size = 6),strip.placement = "outside",strip.switch.pad.wrap = unit(0, "cm"),strip.text = element_text(face = "bold",size=10),strip.background = element_rect(color = "gray50",fill="white"),panel.background = element_rect(fill = 'white'),axis.line= element_line(colour = 'black'),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank() )+
  coord_cartesian(clip="off")

ggsave('heatmap_all_segements_raw.pdf',width = 14,height = 8,plot = p_heatmap,device = cairo_pdf)

#  coord_fixed(ratio = 1)
#panel_border(color = "black") 

#panel.background = element_rect(fill = 'green')
## the following code to change the color for strip
# g <- ggplot_gtable(ggplot_build(p_heatmap))
# stripr <- which(grepl('strip-b', g$layout$name))
# fills <- c(rep(c("gray50","gray80"),11))
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# p_heatmap <- ggdraw(g,)
# p_heatmap


p_heatmap_legend <- get_legend(
  # create some space to the left of the legend
  p_heatmap + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_heatmap <- p_heatmap+theme(legend.position = "none")



# Sample Annotation --------------------------------------------------------------
annodata <- sample_order %>% 
  left_join(
    ngspurity %>% select(Tumor_Barcode,WGD_Status=MCN_WGD,Tumor_Purity=BB_Purity,Tumor_Ploidy=BB_Ploidy) 
  ) %>% 
  left_join(wgs_data %>% select(Tumor_Barcode,Subtype) %>% unique())
#left_join(clinical_data)

# purity
annodata_1 <- annodata %>% select(Seq,Tumor_Purity) 
p_a1 <- annodata_1 %>% 
  ggplot(aes("anno1",Seq,fill=Tumor_Purity))+
  geom_tile()+
  scale_fill_viridis_c(option = "D")+
  labs(fill="Purity")+
  theme_void()+
  scale_y_discrete(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))

p_a1_legend <- get_legend(
  # create some space to the left of the legend
  p_a1 + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_a1 <- p_a1+theme(legend.position = "none")

# WGD
annodata_2 <- annodata %>% select(Seq,WGD_Status) 
p_a2 <- annodata_2 %>% 
  ggplot(aes("anno2",Seq,fill=WGD_Status))+
  geom_tile()+
  scale_fill_aaas()+
  theme_void()+
  scale_y_continuous(breaks = NULL,expand = c(0,0))
#limits = c(1, nrow(label(hcdata)))


p_a2_legend <- get_legend(
  # create some space to the left of the legend
  p_a2 + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_a2 <- p_a2+theme(legend.position = "none")


# Subtype
annodata_3 <- annodata %>% select(Seq,Subtype) 
p_a3 <- annodata_3 %>% 
  ggplot(aes("anno3",Seq,fill=Subtype))+
  geom_tile()+
  scale_fill_manual(values = subtypecol)+
  theme_void()+
  scale_y_continuous(breaks = NULL,expand = c(0,0))
#limits = c(1, nrow(label(hcdata)))


p_a3_legend <- get_legend(
  # create some space to the left of the legend
  p_a3 + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_a3 <- p_a3+theme(legend.position = "none")


# # Smoker
# annodata_4 <- annodata %>% select(Seq,Smoker) 
# p_a4 <- annodata_4 %>% 
#   ggplot(aes("anno4",Seq,fill=Smoker))+
#   geom_tile()+
#   scale_fill_npg()+
#   theme_void()+
#   scale_y_continuous(breaks = NULL,expand = c(0,0))
# #limits = c(1, nrow(label(hcdata)))
# 
# 
# p_a4_legend <- get_legend(
#   # create some space to the left of the legend
#   p_a4 + theme(legend.box.margin = margin(8, 0, 5, 0))
# )
# 
# p_a4 <- p_a4+theme(legend.position = "none")
# 
# # Population  
# annodata_5 <- annodata %>% select(Seq,Assigned_Population) 
# p_a5 <- annodata_5 %>% 
#   ggplot(aes("anno4",Seq,fill=Assigned_Population))+
#   geom_tile()+
#   scale_fill_npg()+
#   theme_void()+
#   scale_y_continuous(breaks = NULL,expand = c(0,0))
# #limits = c(1, nrow(label(hcdata)))
# 
# 
# p_a5_legend <- get_legend(
#   # create some space to the left of the legend
#   p_a5 + theme(legend.box.margin = margin(8, 0, 5, 0))
# )
# 
# p_a5 <- p_a5+theme(legend.position = "none")
# Test chr19 loss  ---------------------------------------------------------
#annodata %>% left_join(chr19loss) %>% mutate(chr19loss=if_else(is.na(chr19loss),"N",chr19loss) ) %>% select(RTK_Altered_Status,chr19loss) %>% mutate(chr19loss=factor(chr19loss,levels = c('Y','N'))) %>% table() %>% fisher.test(alternative = 'greater')


# Frequency_Plot ----------------------------------------------------------
# use bin window as gene id 

hg38 <- left_join(
  cnvdata %>% group_by(chr) %>% arrange(startpos) %>% slice(1) %>% select(chr,startpos) %>% ungroup(),
  cnvdata %>% group_by(chr) %>% arrange(desc(endpos)) %>% slice(1) %>% select(chr,endpos) %>% ungroup()
) %>% mutate(len=endpos-startpos+1)

bin=1000000
hg38bins <- NULL
for(i in 1:23){
  size=hg38$len[i]
  chr=hg38$chr[i]
  startx=hg38$startpos[i]
  endx=hg38$endpos[i]
  tmp <- tibble(chr=chr,start=seq(from = startx,to=endx,by = bin),end=c(seq(from = startx-1,to=endx,by = bin)[-1],size))
  lastval <- as.integer(tail(tmp,1)[3])
  if( lastval < endx){
    tmp <- bind_rows(tmp,tibble(chr=chr,start=lastval+1,end=endx))
  }
  #tmp <- tmp %>% mutate(start2=start+startx-1,end2=end+startx-1)
  hg38bins <- bind_rows(hg38bins,tmp)
}

hg38info <- 
  hg38bins %>% 
  mutate(geneid=paste(chr,start,end,sep="_"),genename=geneid) %>%
  select(chrom=chr,start,end,geneid,genename) %>% 
  mutate(start=as.integer(start),end=as.integer(end)) %>% 
  as.data.frame()



# Freq_LOH ----------------------------------------------------------------
segdata3 <- cnvdata %>% 
  mutate(LOH=if_else(clone_nMin==0 & clone_nMaj==2,1,0)) %>% 
  mutate(num.mark=1000) %>% 
  select(ID=Tumor_Barcode,chrom=chr,loc.start=startpos,loc.end=endpos,num.mark,seg.mean=LOH) %>%
  mutate(seg.mean=if_else(is.na(seg.mean),0,seg.mean)) %>% 
  as.data.frame()
segdata3 %>% filter_all(any_vars(is.na(.)))
cnseg3 <- CNSeg(segdata3)

rdseg3 <- getRS(cnseg3, by = "gene", imput = FALSE, XY = FALSE, what = "mean",geneMap = hg38info)
reducedseg3 <- rs(rdseg3)
msegdata3 <- as.matrix(reducedseg3[,-(1:5)])
msegdata3 <- apply(msegdata3, 2, as.numeric)


freqdata_loh <- bind_cols(
  reducedseg3[,1:3],
  as_data_frame(msegdata3)
) %>% 
  pivot_longer(cols = -c(chrom,start,end)) %>% 
  as_tibble() %>% 
  mutate(start=as.integer(start),end=as.integer(end))

freqdata_loh <- freqdata_loh %>% 
  filter(value!=0) %>% 
  mutate(calling="LOH")


freqdata_loh <- freqdata_loh  %>% 
  count(chrom,start,end,calling) %>% 
  mutate(freq=-n/totalsamples)

# remove cytoband centro
hg38centro2 <- bed_intersect(hg38centro,freqdata_loh %>% mutate(calling=as.character(calling)),suffix = c('','.y')) %>% select(chrom:Centro)

freqdata_loh <- bed_intersect(freqdata_loh %>% mutate(calling=as.character(calling)),hg38centro,invert = T)
freqdata_loh <- bind_rows(
  freqdata_loh %>% mutate(Centro=0),
  hg38centro2 %>% mutate(calling="LOH",n=0,freq=-1e-36)
) %>% 
  arrange(chrom,start,end) %>% select(-Centro)


freqdata_loh <- bind_rows(
  freqdata_loh %>% rename(pos=start) %>% select(-end),
  freqdata_loh %>% rename(pos=end) %>% select(-start)
) %>%
  filter(!is.na(calling)) %>% 
  unique()

freqdata_loh <- freqdata_loh %>% arrange(chrom,pos) %>% mutate(type="Del",calling=NA_character_) %>% mutate(chrom=factor(chrom,levels = chrlevels))


# Freq_CNV ----------------------------------------------------------------
# segdata <- cnvdata %>% 
#   mutate(num.mark=1000) %>% 
#   select(ID=Tumor_Barcode,chrom=chr,loc.start=startpos,loc.end=endpos,num.mark,seg.mean=relative_copy) %>%
#   mutate(seg.mean=if_else(is.na(seg.mean),0,seg.mean)) %>% 
#   as.data.frame()
# segdata %>% filter_all(any_vars(is.na(.)))
# cnseg <- CNSeg(segdata)
rdseg2 <- getRS(cnseg, by = "gene", imput = FALSE, XY = TRUE, what = "mean",geneMap = hg38info)
reducedseg2 <- rs(rdseg2)
msegdata2 <- as.matrix(reducedseg2[,-(1:5)])
msegdata2 <- apply(msegdata2, 2, as.numeric)


freqdata <- bind_cols(
  reducedseg2[,1:3],
  as_data_frame(msegdata2)
) %>% 
  pivot_longer(cols = -c(chrom,start,end)) %>% 
  as_tibble() %>% 
  mutate(start=as.integer(start),end=as.integer(end))

#version1
# freqdata <- freqdata %>% 
#   mutate(calling=case_when(
#     value < log2(1/2) ~ "-2",
#     value < 0 & value >=log2(1/2) ~ "-1",
#     value ==0 ~ NA_character_,
#     value>0 & value < log2(3/2)  ~ "1",
#     value > log2(3/2) ~ "2"
#   )) %>% 
#   filter(value!=0) %>% 
#   mutate(calling=factor(calling,levels = c("-2","-1","2","1")))%>%
#   mutate(type=if_else(calling %in% c("1","2"),"Amp","Del"))

cplevels <- c("Amplification (Copy gain >=4)","Gain","Loss","Homozygous deletion")

freqdata <- freqdata %>% 
  mutate(calling=case_when(
    value <= -4 ~ "Homozygous deletion",
    value < 0 & value > -4 ~ "Loss",
    value ==0 ~ NA_character_,
    value>0 & value < 4  ~ "Gain",
    value >= 4 ~ "Amplification (Copy gain >=4)"
  )) %>% 
  filter(value!=0) %>% 
  mutate(calling=factor(calling,levels = cplevels)) %>% 
  mutate(type=if_else(calling %in% cplevels[1:2],"Amp","Del"))


freqdata <- freqdata  %>% 
  count(chrom,start,end,type,calling) %>% 
  mutate(freq=if_else(calling %in% cplevels[1:2],n/totalsamples,-n/totalsamples)) 

# add boundary to line up heatmap 
#heatscale <- cnvdata %>% select(chrom=chr,startpos,endpos) %>% group_by(chrom) %>% summarise(start=min(startpos),end=max(endpos)) %>% ungroup() %>% mutate(Centro=0)
# remove cytoband centro
hg38centro2 <- bed_intersect(hg38centro,freqdata %>% mutate(calling=as.character(calling)),suffix = c('','.y')) %>% select(chrom:Centro)

freqdata <- bed_intersect(freqdata %>% mutate(calling=as.character(calling)),hg38centro,invert = T)
freqdata <- bind_rows(
  freqdata %>% mutate(Centro=0),
  hg38centro2 %>% mutate(type="Amp",calling="Amplification (Copy gain >=4)",n=0,freq=1e-36),
  hg38centro2 %>% mutate(type="Amp",calling="Gain",n=0,freq=1e-36),
  hg38centro2 %>% mutate(type="Del",calling="Loss",n=0,freq=-1e-36),
  hg38centro2 %>% mutate(type="Del",calling="Homozygous deletion",n=0,freq=-1e-36)
) %>% mutate(calling=factor(calling,levels = cplevels)) %>% 
  arrange(chrom,start,end,type,calling) %>% select(-Centro)

# %>% bind_rows(
#   heatscale %>% mutate(type="Amp",calling="2",n=0,freq=1e-36),
#   heatscale %>% mutate(type="Amp",calling="1",n=0,freq=1e-36),
#   heatscale %>% mutate(type="Del",calling="-1",n=0,freq=-1e-36),
#   heatscale %>% mutate(type="Del",calling="-2",n=0,freq=-1e-36)
# ) 


freqdata <- bind_rows(
  freqdata %>% rename(pos=start) %>% select(-end),
  freqdata %>% rename(pos=end) %>% select(-start)
) %>%
  filter(!is.na(calling)) %>% unique()
#group_by(chrom,pos,type,calling) %>% 
#arrange(desc(freq)) %>% 
#slice(1) %>%
#ungroup()


freqdata <- freqdata %>% 
  pivot_wider(id_cols = -c(type,n),names_from = "calling",values_from = "freq") %>% 
  pivot_longer(cols = -c(chrom,pos),names_to = "calling",values_to = "freq") %>%
  mutate(sig=if_else(calling %in% cplevels[3:4],-1e-36,1e-36),freq=if_else(is.na(freq),sig,freq)) %>% 
  select(-sig) %>%
  arrange(chrom,pos,calling)


# overlap freqdata and cnvdata
tmp <- cnvdata %>% group_by(chr) %>% summarise(start=min(startpos,na.rm = TRUE),end=max(endpos,na.rm = TRUE)) %>% select(chrom=chr,start,end)
freqdata <- freqdata %>% left_join(tmp) %>% filter(pos>=start,pos<=end) %>% select(-start,-end)
freqdata_loh <- freqdata_loh %>% left_join(tmp) %>% filter(pos>=start,pos<=end) %>% select(-start,-end) %>% mutate(chrom=factor(chrom,levels = chrlevels))


freqdata_extra1 <- 
  freqdata %>% group_by(chrom) %>% summarise(pos=min(pos)) %>% 
  left_join(tmp %>% select(chrom,start)) %>% 
  filter(pos!=start) %>% left_join(freqdata) %>% 
  mutate(pos=start) %>% select(-start)

freqdata_extra2 <-  
  freqdata %>% group_by(chrom) %>% summarise(pos=max(pos)) %>% 
  left_join(tmp %>% select(chrom,end)) %>% 
  filter(pos!=end) %>% left_join(freqdata) %>% 
  mutate(pos=end) %>% select(-end)

freqdata <- bind_rows(freqdata,freqdata_extra1,freqdata_extra2) %>% arrange(chrom,pos,calling)


# freqdata <- 
#   freqdata %>% mutate(pos=(start+end)/2) %>%
#   select(-start,-end) %>%
#   arrange(chrom,pos,type,calling) 

#freqdata <- freqdata %>% complete(nesting(chrom,pos),calling,fill=list(n=0,freq=0.0001))

#freqdata2 <- freqdata %>% filter(chrom=="1") %>% arrange(chrom,pos,freq,calling)
#%>% slice(421:1500)

cnvcolor <- rev(c('#2166ac','#92c5de','#f4a582','#b2182b'))
#names(cnvcolor) <- c('-2','-1','1','2')
names(cnvcolor) <- cplevels

hg38centro_data <- hg38centro %>% select(chrom,pos=end) %>% group_by(chrom) %>% arrange(pos) %>% slice(1) %>% ungroup() %>% mutate(n=0,freq=0,type="Del",chrom=factor(chrom,levels = chrlevels))

hg38centro_data <- hg38centro_data %>% left_join(tmp) %>% filter(pos>=start,pos<=end) %>% select(-start,-end) %>% mutate(chrom=factor(chrom,levels = chrlevels))

#heatscale <- cnvdata %>% select(chrom=chr,startpos,endpos) %>% group_by(chrom) %>% summarise(start=min(startpos),end=max(endpos)) %>% ungroup() 
#scales_ <- list(heatscale$chrom=scale_x_continuous(limits = c(heatscale$start,heatscale$end)))

p_freq <- freqdata %>% 
  mutate(chrom=factor(chrom,levels = chrlevels)) %>% 
  ggplot(aes(x = pos,y=freq,fill=factor(calling,levels = cplevels[c(1,2,4,3)])))+
  geom_area(position = "stack",linetype=2)+
  #geom_bar(data=freqdata2 %>% filter(type=="Del"),stat="identity",position = "stack")+
  scale_fill_manual(values = cnvcolor)+
  labs(fill="SCNA Calling")+
  geom_hline(yintercept = 0,color="gray40")+
  #facet_wrap(~chr,nrow = 1,scales = 'free_x',strip.position = "bottom")+
  facet_grid(~chrom,scales = 'free_x',space = 'free')+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
  scale_y_continuous(expand = expand_scale(add = c(0, 0)), sec.axis = dup_axis(),limits = c(-1,1),breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),labels = function(x) paste0(x*100, "%"))+
  labs(x="",y="% CNV gain/loss, copy neutral LOH")+
  theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0, "lines"),panel.grid.major.y = element_line(colour = 'gray50',linetype = 5),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.placement = "inside",strip.switch.pad.wrap = unit(0, "cm"),strip.text = element_blank(),strip.background = element_blank(),panel.background = element_rect(fill = 'white'),axis.line= element_line(colour = 'black'),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank() )+
  coord_cartesian(clip="off")+panel_border(colour = 'gray80',size = 0.2)+
  geom_line(data=freqdata_loh,aes(pos,freq),colour="black",size=0.4)+
  geom_point(data=hg38centro_data,aes(pos,freq),size=1.5,pch=21,fill="black")

p_freq <- flush_ticks(gg = p_freq)+theme(axis.text.x = element_blank())

p_freq_legend <- get_legend(
  # create some space to the left of the legend
  p_freq + theme(legend.box.margin = margin(8, 0, 5, 0))
)

p_freq <- p_freq+theme(legend.position = "none")


# frequency
tmp1 <- freqdata %>% filter(calling %in% c('Gain')) %>% group_by(chrom) %>% summarise(freq1=median(freq,na.rm=T)) %>% arrange(desc(freq1))
tmp2 <- freqdata %>% filter(calling %in% c('Amplification (Copy gain >=4)')) %>% group_by(chrom) %>% summarise(freq2=median(freq,na.rm=T)) %>% arrange(desc(freq2))
left_join(tmp1,tmp2) %>% 
  mutate(freq=freq1+freq2) %>% 
  arrange(desc(freq))

freqdata %>% filter(calling %in% c('Gain'),pos<33200000,chrom=='12') %>% group_by(chrom) %>% summarise(freq1=median(freq,na.rm=T)) %>% arrange(desc(freq1))
freqdata %>% filter(calling %in% c('Amplification (Copy gain >=4)'),pos<33200000,chrom=='12') %>% group_by(chrom) %>% summarise(freq1=median(freq,na.rm=T)) %>% arrange(desc(freq1))

freqdata %>% filter(calling %in% c('Gain'),pos>37800000,chrom=='12') %>% group_by(chrom) %>% summarise(freq1=median(freq,na.rm=T)) %>% arrange(desc(freq1))
freqdata %>% filter(calling %in% c('Amplification (Copy gain >=4)'),pos>37800000,chrom=='12') %>% group_by(chrom) %>% summarise(freq1=median(freq,na.rm=T)) %>% arrange(desc(freq1))

freqdata %>% filter(calling %in% c('Gain'),pos>13000000,chrom=='21') %>% group_by(chrom) %>% summarise(freq1=median(freq,na.rm=T)) %>% arrange(desc(freq1))
freqdata %>% filter(calling %in% c('Amplification (Copy gain >=4)'),pos>13000000,chrom=='21') %>% group_by(chrom) %>% summarise(freq1=median(freq,na.rm=T)) %>% arrange(desc(freq1))


# Genome_Coverage ---------------------------------------------------------
cnv_freq <- 
  cnvdata %>% 
  select(Tumor_Barcode,chr,startpos,endpos,relative_copy) %>% 
  filter(relative_copy!=0) %>% 
  mutate(size=endpos-startpos) %>% 
  group_by(Tumor_Barcode) %>% 
  summarise(total=sum(size)) %>% 
  mutate(coverage=total/2881033286) %>% 
  right_join(cnvclust) %>% 
  mutate(coverage=if_else(is.na(coverage),0,coverage)) %>% 
  left_join(sample_order) 
cnv_freq %>% 
  ggplot(aes(Seq,coverage,fill=CNV_Clust))+geom_col()+coord_flip()

cnv_coverage <- cnv_freq %>% select(Tumor_Barcode,CNV_Coverage=coverage)
#save(cnv_coverage,file='cnv_coverage.RData')


# Combined figures --------------------------------------------------------
p_dend <- p_dend+theme(plot.margin=margin(r=0,unit="cm"))
p_heatmap <- p_heatmap+theme(plot.margin=margin(l=-1,unit="cm"))+panel_border(colour = 'gray80',size = 0.2)+theme(axis.text.y = element_blank())
p_freq <- p_freq+theme(plot.margin=margin(b=-0.4,t=0.3,l=-1,unit="cm"))

p_com <- align_plots(p_freq,
                     p_heatmap,
                     align = 'v',
                     axis = 'lr'
)

# p_combined_plot <- plot_grid(p_dend,p_com[[2]] ,p_a1,p_a2,p_a3,p_a4,p_a5,align = 'h',axis = "tb",rel_widths = c(0.5,1.5,0.02,0.02,0.02,0.02,0.02),nrow = 1)
# p_combined_plot2 <- plot_grid(NULL,p_com[[1]],NULL,NULL,NULL,NULL,NULL,align = 'h',rel_widths = c(0.5,1.5,0.02,0.02,0.02,0.02,0.02),nrow = 1)
# p_combined_legend <- plot_grid(NULL,NULL,NULL,p_freq_legend,p_heatmap_legend,p_a1_legend,p_a2_legend,p_a3_legend,p_a4_legend,p_a5_legend,NULL,NULL,nrow = 1,align = 'h')

p_combined_plot <- plot_grid(p_dend,p_com[[2]] ,p_a1,p_a2,p_a3,align = 'h',axis = "tb",rel_widths = c(0.5,1.5,0.02,0.02,0.02),nrow = 1)
p_combined_plot2 <- plot_grid(NULL,p_com[[1]],NULL,NULL,NULL,align = 'h',rel_widths = c(0.5,1.5,0.02,0.02,0.02),nrow = 1)
p_combined_legend <- plot_grid(NULL,NULL,NULL,p_freq_legend,p_heatmap_legend,p_a1_legend,p_a2_legend,p_a3_legend,NULL,NULL,nrow = 1,align = 'h')

p_combined <- plot_grid(p_combined_plot2,p_combined_plot,p_combined_legend,nrow = 3,align = "h",rel_heights = c(2,10,1.2))

ggsave(file="heatmap_all_segements_raw_tmp.pdf",plot = p_combined,width = 22,height = 17,device = cairo_pdf)
#ggsave(file="heatmap_all_clone.pdf",plot = p_combined,width = 22,height = 17,device = cairo_pdf)
#ggsave(file="heatmap_all_subclone.pdf",plot = p_combined,width = 20,height = 17,device = cairo_pdf)



 
# #save.image(file = 'bb_heatmap.RData')







# isochromosome 12p -------------------------------------------------------
hg38centro %>% filter(chrom=='12')
seg12 <- BBprofile %>% 
  filter(chr=='12') %>% 
  mutate(seglen=abs(endpos-startpos)) %>% 
  mutate(arm=if_else(startpos>35500000,'chr12q',if_else(endpos<35500000,'chr12p',if_else(abs(35500000 - startpos) > abs(35500000 - endpos), 'chr12p','chr12q')))) %>% 
  select(Tumor_Barcode,arm,chr,startpos,endpos,seglen,BAF,LogR,contains('clone')) %>%
  mutate(clone_total=clone_nMaj + clone_nMin) %>% 
  group_by(Tumor_Barcode,arm) %>% 
  arrange(desc(seglen)) %>% 
  slice(1) %>% 
  ungroup()

sample_order2 <- seg12 %>% 
  select(Tumor_Barcode,arm,clone_total) %>% 
  pivot_wider(names_from=arm,values_from=clone_total,values_fill=0) %>% 
  mutate(isochromsome=if_else((chr12p-chr12q)<1,'No','Yes')) %>% 
  arrange(desc(chr12p-chr12q)) %>% 
  mutate(Seq2=seq_along(Tumor_Barcode)) 
  


## heatmap



# Dendrograms -------------------------------------------------------------
reducedseg_chr12 <- reducedseg[reducedseg$chrom=='12',]
msegdata <- as.matrix(reducedseg_chr12[,-(1:3)])
msegdata <- apply(msegdata, 2, as.numeric)

hc <- hclust(dist(t(msegdata),method = 'euclidean'),method = 'ward.D')
hcdata <- dendro_data_k(hc, 3)

cols <- c("#a9a9a9", "#2ca02c", "#ff7f0e","#1f77b4")

p_dend <- plot_ggdendro(hcdata,
                        direction   = "lr",
                        scale.color = cols,
                        label.size  = 2.5,
                        branch.size = 0.5,
                        expand.y    = 0,nudge.label = 0)

p_dend <- p_dend+
  scale_x_continuous(breaks = NULL,limits = c(1, nrow(label(hcdata))),expand = c(0,0))+
  theme_void()

sample_order <- as_tibble(hcdata$labels) %>% 
  select(Tumor_Barcode=label) %>%
  unique() %>%
  mutate(Seq=seq_along(Tumor_Barcode)) %>% 
  left_join(sample_order2) %>% 
  arrange(isochromsome,Seq) %>% 
  mutate(Seq=seq_along(Tumor_Barcode))
  

#chr19loss <- sample_order %>% tail(n=11) %>% mutate(chr19loss="Y") %>% select(-Seq)
#chr19loss <- cnvdata %>% filter(chr=="19",clone_nMaj==0) %>% mutate(size=endpos-startpos) %>% filter(size>30000000) %>% mutate(chr19loss="Y") %>% select(Tumor_Barcode,chr19loss)
#save(chr19loss,file='chr19loss.RData')

cnvclust <- as_tibble(hcdata$labels) %>% select(Tumor_Barcode=label,CNV_Clust=clust) %>% mutate(CNV_Clust=paste0('C',CNV_Clust))

cnvdata %>% 
  mutate(chr=factor(chr,levels = chrlevels)) %>%
  mutate(startpos=as.integer(startpos),endpos=as.integer(endpos)) %>% 
  left_join(sample_order) %>% 
  arrange(Seq,chr) %>% 
  filter(chr %in% c(12)) %>% 
  ggplot(aes(fill=relative_copy))+
  #geom_segment(aes(x=startpos,xend=endpos,y=Tumor_Barcode,yend=Tumor_Barcode),size=3)+
  geom_rect(aes(xmin=startpos,xmax=endpos,ymin=Seq-0.5,ymax=Seq+0.5))+
  facet_grid(~chr,scales = 'free_x',space = 'free',switch="both")+
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",midpoint = 0)+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
  labs(fill="Relative Copy Number")+
  scale_y_continuous(breaks = sample_order$Seq,labels = sample_order$Tumor_Barcode,expand = expand_scale(add = c(0, 0)), sec.axis = dup_axis())+
  #facet_wrap(~chr,nrow = 1,scales = 'free_x',strip.position = "bottom")+
  labs(x="",y="")+
  #theme_ipsum_rc()+
  theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0, "lines"),axis.text.x = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),strip.placement = "outside",strip.switch.pad.wrap = unit(0, "cm"),strip.text = element_text(face = "bold",size=10),strip.background = element_rect(color = "gray50",fill="white"),panel.background = element_rect(fill = 'white'),axis.line= element_line(colour = 'black'),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank() )+
  coord_cartesian(clip="off")+
  geom_hline(yintercept = 30.5,linetype = 1)

ggsave(filename = 'heatmap_all_segements_isochromosome12p_raw.pdf',width = 5,height = 10,device = cairo_pdf)

isochromosome_12p <- sample_order
save(isochromosome_12p,file='isochromosome_12p.RData')



# Association with KIT mutations ------------------------------------------
load('tgct_maf.RData')

tdata <- wgs_data %>% 
  left_join(
    isochromosome_12p %>% select(Tumor_Barcode, isochromsome)  
  ) %>% 
  left_join(
    tgct_maf %>% filter(Hugo_Symbol == 'KIT') %>% select(Tumor_Barcode,aaChange) %>% mutate(KIT='Mutant')
  ) %>% 
  mutate(KIT = replace_na(KIT,'Wildtype'))


source('~/NIH-Work/R/ZTW_function/Sherlock_functions.R')


barplot_fisher(mdata0 = tdata %>% filter(Subtype=='Seminoma'),var1name = "isochromsome",var1lab = 'Isochromosome 12p', var2name = 'KIT',var2lab = 'KIT',Pcol = ncicolpal[1],filename = 'seminoma_i12p_kit.pdf',width = 4.5,height = 6.1 )

