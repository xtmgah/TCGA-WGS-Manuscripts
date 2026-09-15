# Genome-wide copy-number frequency calculation with an explicit state parameter.
# Uses relative copy-number bins, centromere masks, CNTools means and chr21 filtering.
ptc_cnv_frequency <- function(cnv,ngspurity,cytoband_file,state=c("clone","subclone")) {
state <- match.arg(state)
totalsamples <- cnv %>% count(Tumor_Barcode) %>% pull(Tumor_Barcode) %>% length()
allsamples <- cnv %>% count(Tumor_Barcode) %>% pull(Tumor_Barcode) %>% unique()
cnvdata <- cnv %>%
  filter(chr %in% c(1:22)) %>%
  mutate(clone_total=if (state=='clone') clone_nMin+clone_nMaj else subclone_nMin+subclone_nMaj) %>%
  left_join(
    ngspurity %>% select(Tumor_Barcode,WGD_Status=MCN_WGD,Tumor_Purity=BB_Purity,Tumor_Ploidy=BB_Ploidy)
  ) %>%
  mutate(WGD_Status=if_else(is.na(WGD_Status),'nWGD',WGD_Status)) %>%
  mutate(relative_copy=clone_total-if_else(WGD_Status=="WGD",4,2)) %>%
  mutate(relative_copy=if_else(relative_copy>4,4,relative_copy)) %>%
  mutate(relative_copy=if_else(relative_copy< -4,-4,relative_copy)) %>%
  mutate(relative_copy=if_else(clone_nMaj==0,-4,relative_copy)) %>%
  mutate(startpos=as.integer(startpos),endpos=as.integer(endpos))
cyto <- readr::read_tsv(cytoband_file,col_names=FALSE,show_col_types=FALSE)
colnames(cyto) <- c('Chromosome', 'Start' , 'End', 'Band', 'gieStain')
cyto <- cyto %>% mutate(Chromosome = str_remove(Chromosome,"chr"),Centro=if_else(gieStain=="acen",1,0), Arm=paste0(Chromosome,str_sub(Band,1,1))) %>% filter(Chromosome %in% c(1:22,"X","Y"))
hg38centro <- cyto %>% filter(Centro!=0,Chromosome %in% c(1:22)) %>% select(chrom=Chromosome,start=Start,end=End,Centro)
chrlevels <- seq(1:22)
hg38 <- left_join(
  cnvdata %>% group_by(chr) %>% arrange(startpos) %>% slice(1) %>% select(chr,startpos) %>% ungroup(),
  cnvdata %>% group_by(chr) %>% arrange(desc(endpos)) %>% slice(1) %>% select(chr,endpos) %>% ungroup()
) %>% mutate(len=endpos-startpos+1)
bin=1000000
hg38bins <- NULL
for(i in 1:22){
  size=hg38$len[i]
  chr=hg38$chr[i]
  startx=hg38$startpos[i]
  endx=hg38$endpos[i]
  tmp <- tibble(chr=chr,start=seq(from = startx,to=endx,by = bin),end=c(seq(from = startx-1,to=endx,by = bin)[-1],size))
  lastval <- as.integer(tail(tmp,1)[3])
  if( lastval < endx){
    tmp <- bind_rows(tmp,tibble(chr=chr,start=lastval+1,end=endx))
  }
  hg38bins <- bind_rows(hg38bins,tmp)
}
hg38info <-
  hg38bins %>%
  mutate(geneid=paste(chr,start,end,sep="_"),genename=geneid) %>%
  select(chrom=chr,start,end,geneid,genename) %>%
  mutate(start=as.integer(start),end=as.integer(end)) %>%
  as.data.frame()
segdata <- cnvdata %>%
  mutate(num.mark=1000) %>%
  select(ID=Tumor_Barcode,chrom=chr,loc.start=startpos,loc.end=endpos,num.mark,seg.mean=relative_copy) %>%
  mutate(seg.mean=if_else(is.na(seg.mean),0,seg.mean)) %>%
  as.data.frame()
cnseg <- CNSeg(segdata)
rdseg2 <- getRS(cnseg, by = "gene", imput = FALSE, XY = FALSE, what = "mean",geneMap = hg38info)
reducedseg2 <- rs(rdseg2)
msegdata2 <- as.matrix(reducedseg2[,-(1:5)])
msegdata2 <- apply(msegdata2, 2, as.numeric)
freqdata <- bind_cols(
  reducedseg2[,1:3],
  tibble::as_tibble(msegdata2)
) %>%
  pivot_longer(cols = -c(chrom,start,end)) %>%
  as_tibble() %>%
  mutate(start=as.integer(start),end=as.integer(end))
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
freqdata <- bind_rows(
  freqdata %>% rename(pos=start) %>% select(-end),
  freqdata %>% rename(pos=end) %>% select(-start)
) %>%
  filter(!is.na(calling)) %>% unique()
freqdata <- freqdata %>%
  pivot_wider(id_cols = -c(type,n),names_from = "calling",values_from = "freq") %>%
  pivot_longer(cols = -c(chrom,pos),names_to = "calling",values_to = "freq") %>%
  mutate(sig=if_else(calling %in% cplevels[3:4],-1e-36,1e-36),freq=if_else(is.na(freq),sig,freq)) %>%
  select(-sig) %>%
  arrange(chrom,pos,calling)
tmp <- cnvdata %>% group_by(chr) %>% summarise(start=min(startpos,na.rm = TRUE),end=max(endpos,na.rm = TRUE)) %>% select(chrom=chr,start,end)
freqdata <- freqdata %>% left_join(tmp) %>% filter(pos>=start,pos<=end) %>% select(-start,-end)
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
cnvcolor <- rev(c('#2166ac','#92c5de','#f4a582','#b2182b'))
names(cnvcolor) <- cplevels
hg38centro_data <- hg38centro %>% select(chrom,pos=end) %>% group_by(chrom) %>% arrange(pos) %>% slice(1) %>% ungroup() %>% mutate(n=0,freq=0,type="Del",chrom=factor(chrom,levels = chrlevels))
hg38centro_data <- hg38centro_data %>% left_join(tmp) %>% filter(pos>=start,pos<=end) %>% select(-start,-end) %>% mutate(chrom=factor(chrom,levels = chrlevels))
freqdata_clean <- freqdata %>% filter(!(chrom==21 & (freq>0.05 | freq< -0.05)))
p_freq <- freqdata_clean %>%
  mutate(chrom=factor(chrom,levels = chrlevels)) %>%
  ggplot(aes(x = pos,y=freq,fill=factor(calling,levels = cplevels[c(1,2,4,3)])))+
  geom_area(position = "stack",linetype=2)+
  scale_fill_manual(values = cnvcolor)+
  labs(fill="SCNA Calling")+
  geom_hline(yintercept = seq(-0.10, 0.10, by = 0.05), colour = "gray50", linetype = "dotted", size = 0.3, inherit.aes = FALSE)+
  geom_hline(yintercept = 0,color="gray40",size=0.3)+
  facet_grid(~chrom,scales = 'free_x',space = 'free')+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0),sec.axis = dup_axis())+
  scale_y_continuous(expand = expand_scale(add = c(0, 0)), limits = c(-0.15,0.15),breaks =c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15),labels = percent_format())+
  labs(x="",y="% CNV gain/loss, copy neutral LOH")+
  theme(text = element_text(family = "Roboto Condensed"),panel.spacing = unit(0.001, "lines"),panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.placement = "inside",strip.switch.pad.wrap = unit(0, "cm"),strip.background = element_blank(),panel.background = element_rect(fill = 'white'),axis.line= element_blank(),axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y = element_blank(),axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),legend.position = 'top' )+
  coord_cartesian(clip="off")+ panel_border(colour = 'black',size = 0.15)+
  geom_point(data=hg38centro_data,aes(pos,freq),size=0.5,pch=21,fill="black")
p_freq <- hrbrthemes::flush_ticks(gg = p_freq)+theme(axis.text.x = element_blank())
p_freq + ptc_theme() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing=grid::unit(.02,"lines"),strip.text=element_text(size=11),legend.position="top")
}
