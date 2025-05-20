set_wd()
libztw()
library(maftools)
### for the subject of 1263
load('wgs_data.RData')

load('intogene_drivers.RData')
geneset0 <- intogene$SYMBOL

# Comparision between Smoking and Non-smoking  ----------------------------
load('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/maf_orignal.RData')


# Comparison between Seminoma and Non-Seminoma ----------------------------
samplelist_seminoma <- wgs_data %>% filter(Subtype == 'Seminoma') %>% pull(Tumor_Barcode)
maf_seminoma<- subsetMaf(maf = maf_orignal, tsb = samplelist_seminoma, mafObj = TRUE)
samplelist_nonseminoma <- wgs_data %>% filter(Subtype == 'Non-Seminoma') %>% pull(Tumor_Barcode)
maf_nonseminoma<- subsetMaf(maf = maf_orignal, tsb = samplelist_seminoma, mafObj = TRUE)


pt.vs.rt <- mafCompare(m1 = maf_seminoma, m2 = maf_nonseminoma, m1Name = 'Seminoma', m2Name = 'Non-Seminoma', minMut = 1)
pt.vs.rt$results %>% filter(pval<0.1)
#pt.vs.rt$results <- pt.vs.rt$results %>% filter(Hugo_Symbol %in% geneset0)
pdf('maf_subtype_comparison..pdf',width = 8.5,height = 2.5,useDingbats = FALSE)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05,geneFontSize = 0.8)
dev.off()


## load total MAF
tgct_maf=as_tibble(maf_orignal@data,.name_repair = 'unique') %>% select(Tumor_Sample_Barcode:AAChange.ensGene) %>% rename(Tumor_Barcode=Tumor_Sample_Barcode) %>% filter(Tumor_Barcode %in% wgs_data$Tumor_Barcode)
# frequency
tgct_maf_freq <- tgct_maf %>% 
  select(Tumor_Barcode,Hugo_Symbol) %>% 
  unique() %>%
  count(Hugo_Symbol,name = 'nMutSample') %>% 
  mutate(Mut_Freq=nMutSample/dim(wgs_data)[1]) %>% 
  arrange(desc(Mut_Freq))

tgct_maf_freq_seminoma <- tgct_maf %>% 
  filter(Tumor_Barcode %in% samplelist_seminoma) %>% 
  select(Tumor_Barcode,Hugo_Symbol) %>% 
  unique() %>%
  count(Hugo_Symbol,name = 'nMutSample') %>% 
  mutate(Mut_Freq=nMutSample/length(samplelist_seminoma)) %>% 
  arrange(desc(Mut_Freq))


tgct_maf_freq_nonseminoma <- tgct_maf %>% 
  filter(Tumor_Barcode %in% samplelist_nonseminoma) %>% 
  select(Tumor_Barcode,Hugo_Symbol) %>% 
  unique() %>%
  count(Hugo_Symbol,name = 'nMutSample') %>% 
  mutate(Mut_Freq=nMutSample/length(samplelist_nonseminoma)) %>% 
  arrange(desc(Mut_Freq))

save(tgct_maf_freq,tgct_maf_freq_nonseminoma,tgct_maf_freq_seminoma,tgct_maf,file = 'tgct_maf.RData')

save(maf_orignal,maf_seminoma,maf_nonseminoma,file = 'tgct_maf_orignal.RData')



# Hotspot mutations -------------------------------------------------------





# Proteinpaint ------------------------------------------------------------



qgene <- 'KIT'
#qtranscript <-  'NM_005228'

class_data <- tibble(
  Variant_Classification = c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Translation_Start_Site'),
  class = c('frameshift','frameshift','proteinDel','proteinIns','missense','nonsense','nonsense','splice','Translation_Start_Site')
)


mutated_samples <- tgct_maf %>%
  filter(Hugo_Symbol %in% qgene) %>% 
  pull(Tumor_Barcode) %>% 
  unique()

tdata1 <- tgct_maf %>%
  filter(Hugo_Symbol %in% qgene,Variant_Classification != "Splice_Site") %>% 
  mutate(AAChange.refGene = if_else(is.na(AAChange.refGene),'Splice_Site',AAChange.refGene)) %>% 
  select(Tumor_Barcode,Chromosome,Start_Position,Reference_Allele, Tumor_Seq_Allele2,Variant_Classification,AAChange.refGene) %>% 
  separate_rows(AAChange.refGene,sep = ',') %>% 
  separate(AAChange.refGene,into = c('Gene','TS','Exon','cDNA','AAChange'),sep = ':') %>% 
  mutate(AAPosition = parse_number(str_remove(AAChange,'p.'))) 


tdata2 <- tgct_maf %>%
  filter(Hugo_Symbol %in% qgene, Variant_Classification == "Splice_Site") %>% 
  select(Tumor_Barcode,Chromosome,Start_Position,Reference_Allele, Tumor_Seq_Allele2,Variant_Classification,GeneDetail.refGene) %>% 
  separate_rows(GeneDetail.refGene,sep = ';') %>% 
  separate(GeneDetail.refGene,into = c('TS','Exon','cDNA'),sep = ':') %>% 
  mutate(Gene = qgene) %>% 
  mutate(AAPosition = ceiling(parse_number(str_remove(cDNA,'^c.'))/3)) %>% 
  mutate(AAChange = paste0('X',AAPosition,'_splice'))

tdata0 <- bind_rows(tdata1,tdata2)

tslist <- tdata0 %>% filter(!is.na(TS)) %>% pull(TS) %>% unique()
tslist

#if(!(qtranscript %in% tslist)){print("No transcript found"); break}

tdata <- tdata0 %>% unique() #filter(TS==qtranscript) 

if(length(unique(tdata$Tumor_Barcode)) == length(mutated_samples)){
  tgct_maf %>%
    filter(Hugo_Symbol %in% qgene,!(Tumor_Barcode %in% unique(tdata$Tumor_Barcode))) 
}

tdata %>% 
  left_join(class_data) %>% 
  left_join(wgs_data) %>% 
  mutate(origin = 'somatic') %>% 
  filter(Subtype=='Seminoma') %>% 
  select(disease = Subtype,sample=Tumor_Barcode, origin, gene = Gene, refseq=TS,chromosome=Chromosome, start=Start_Position,aachange=AAChange,class,REF=Reference_Allele, ALT=Tumor_Seq_Allele2) %>% 
  write_delim(paste0('~/Downloads/ProteinPaint_input1_',qgene,'.txt'),delim = '\t',col_names = T)


tdata %>% 
  left_join(class_data) %>% 
  left_join(wgs_data) %>% 
  mutate(origin = 'somatic') %>% 
  filter(Subtype!='Seminoma') %>% 
  select(disease = Subtype,sample=Tumor_Barcode, origin, gene = Gene, refseq=TS,chromosome=Chromosome, start=Start_Position,aachange=AAChange,class,REF=Reference_Allele, ALT=Tumor_Seq_Allele2) %>% 
  write_delim(paste0('~/Downloads/ProteinPaint_input2_',qgene,'.txt'),delim = '\t',col_names = T)
























































## parallele evolution
load('../DP_info_data.RData')
#load('/Volumes/data/NSLC2/NGSpurity2/DP_info_data.RData')
load('/Volumes/data/NSLC/Drive_Gene/DriveGene.RData')


tmpdata <- DP_info_data %>%  filter(str_detect(Variant_Classification,'exonic|splicing|UTR')) 
tmpdata2 <- tmpdata %>% filter(!is.na(AAChange)|str_detect(Variant_Classification,'splicing'),!str_detect(Variant_Classification,'UTR|ncRNA|unknown|@synonymous'))
tmpdata3 <- tmpdata2 %>% select(Tumor_Barcode,Gene_Name,Variant_Classification,Clone,Cluster,CCF,Cluster_CCF) %>% filter(!is.na(Clone))
tmpdata3 <- tmpdata3 %>% count(Tumor_Barcode,Gene_Name) %>% filter(n>1) %>% left_join(tmpdata3)
tmpdata3 <- tmpdata3 %>% select(Tumor_Barcode,Gene_Name,Clone) %>% unique() %>% group_by(Tumor_Barcode,Gene_Name) %>% arrange(Tumor_Barcode,Gene_Name,Clone) %>% summarise(Clone=paste0(Clone,collapse = '_')) %>% ungroup()
tmpdata3 <- tmpdata3 %>% count(Gene_Name,Clone) %>% pivot_wider(names_from = Clone,values_from = n,values_fill = 0) %>% mutate(Total=N+Y+N_Y) %>% arrange(desc(Total))
tmpdata4 <- tmpdata3 %>% filter(Gene_Name %in% drivegene$Gene) %>% filter(Total>1) 


## mutational pattern ## 
tmpdata2 %>% filter(Gene_Name %in% c('EGFR','KRAS','TP53')) %>% select(Tumor_Barcode,context,mutType) %>% mutate(MutationType=str_sub(mutType,3,5)) %>% count(MutationType)

tmpdata2 %>% filter(Gene_Name %in% c('EGFR','KRAS','TP53'))  %>% count(Tumor_Barcode,Gene_Name) %>% filter(n>1) %>% 
  left_join(
    tmpdata2 %>% filter(Gene_Name %in% c('EGFR','KRAS','TP53')) %>% select(Tumor_Barcode,Gene_Name,context,mutType) 
  )%>%
  mutate(MutationType=str_sub(mutType,3,5)) %>% count(MutationType)



## load gene length
genelist <- read_delim('../hg38_refGene.txt',delim = '\t',col_names = F)
genelist <- genelist %>% mutate(Size=abs(X6-X5))%>% select(RefSeq=X2,Gene=X13,Size) %>% 
  group_by(Gene) %>% slice(1) %>% ungroup()


tmpdata4 <- tmpdata4 %>% left_join( genelist %>% rename(Gene_Name=Gene)) %>% mutate(seq=seq_along(Gene_Name))

tmpdata4 %>% ggplot(aes(log2(Size),Total))+geom_point()+ggrepel::geom_text_repel(aes(label=Gene_Name))

# library(scatterpie)
# ggplot()+geom_scatterpie(aes(x=Size,y=Total,group=seq,cols=NA),data=tmpdata4 %>% mutate(Size=Size/10000),cols = colnames(tmpdata4)[2:4])+ coord_equal()


load('tgct_maf.RData')
tmp <- tmpdata3 %>% left_join(tgct_maf_freq %>% rename(Gene_Name=Hugo_Symbol)) 
tmp %>% 
  ggplot(aes(Mut_Freq,Total/1217))+
  geom_point(pch=21,fill='#cccccc',col='black',size=2,stroke=0.5)+
  ggrepel::geom_text_repel(data=tmp %>% filter(Mut_Freq>0.1|Total>25),aes(label=Gene_Name))+
  #geom_smooth(method = 'lm')+
  scale_x_percent()+
  scale_y_percent()+
  labs(x='Overall Mutation Frequency',y='Multiple Mutation Freqeuncy')+
  theme_ipsum_rc(base_size = 14,axis_title_size = 16,axis_title_just = 'm' ,axis = FALSE)+
  panel_border(color = 'black')

ggsave('Multiple_mutation1.pdf',width = 10,height = 8,device = cairo_pdf)

tmp %>% left_join(genelist %>% rename(Gene_Name=Gene)) %>% 
  ggplot(aes(log2(Size/1000),log2(Total)))+
  geom_point(pch=21,fill='#cccccc',col='black',size=2,stroke=0.5)+
  #ggrepel::geom_text_repel(data=tmp %>% filter(Mut_Freq>0.1|Total>25),aes(label=Gene_Name))+
  geom_smooth(method = 'lm')+
  #scale_x_percent()+
  scale_y_percent()+
  labs(x='Overall Mutation Frequency',y='Multiple Mutation Freqeuncy')+
  theme_ipsum_rc(base_size = 14,axis_title_size = 16,axis_title_just = 'm' ,axis = FALSE)+
  panel_border(color = 'black')


source('~/NIH-Work/R/ZTW_function/ztw.R')
## EGFR hotspot ## 
tmp <- tmpdata2 %>% filter(Gene_Name=='EGFR') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:')) %>% count(ID,AAChange) %>% arrange(desc(n)) %>% mutate(Hotspot=if_else(n>10,"Yes",'N')) %>% mutate(Key=if_else(Hotspot=="Yes",AAChange,'Other')) %>% mutate(Key=fct_reorder(Key,n))

PieDonut_ztw(tmp,aes(pies=Key,count=n),mainCol = c('gray50',pal_d3()(5)),showNum = TRUE,showRatioPie = TRUE,showRatioThreshold = 0,pieLabelSize=3.5,titlesize = 5,r0=0,title='Population', showPieName = FALSE,family = 'Roboto Condensed')
ggsave(filename = 'EGFR_hotspot.pdf',width = 8,height = 8)

tmpdata2 %>% filter(Gene_Name=='EGFR') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:'),) %>% count(Clone) %>% arrange(desc(n)) %>% janitor::adorn_percentages(denominator = 'col')

tmp2 <- tmpdata2 %>% filter(Gene_Name=='EGFR') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:')) %>% 
  left_join(tmp) %>% 
  select(Tumor_Barcode,AAChange,Key,Clone,CCF) 

left_join(
  tmp2 %>% filter(Key!="Other")  %>% select(Tumor_Barcode,CCF_hotspot=CCF) %>% group_by(Tumor_Barcode) %>% arrange(desc(CCF_hotspot)) %>% slice(1) %>% ungroup(),
  tmp2 %>% filter(Key=="Other")  %>% select(Tumor_Barcode,CCF_other=CCF)%>% group_by(Tumor_Barcode) %>% arrange(CCF_other) %>% slice(1) %>% ungroup()
) %>% 
  filter(!is.na(CCF_hotspot),!is.na(CCF_other)) %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  mutate(name=str_replace(name,'CCF','EGFR')) %>% 
  ggplot(aes(name,value))+geom_boxplot(width=0.5)+
  labs(x="",y="Cancer cell fraction (CCF)\n")+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = FALSE)+
  scale_y_continuous(breaks = pretty_breaks(n = 10),limits = c(0.5,1.5),expand = c(0,0))+
  panel_border(color = "black")

ggsave(filename = 'EGFR_hotspot2.pdf',width = 3.5,height = 7)


a <- left_join(
  tmp2 %>% filter(Key!="Other")  %>% select(Tumor_Barcode,CCF_hotspot=CCF) %>% group_by(Tumor_Barcode) %>% arrange(desc(CCF_hotspot)) %>% slice(1) %>% ungroup(),
  tmp2 %>% filter(Key=="Other")  %>% select(Tumor_Barcode,CCF_other=CCF)%>% group_by(Tumor_Barcode) %>% arrange(CCF_other) %>% slice(1) %>% ungroup()
) %>% 
  filter(!is.na(CCF_hotspot),!is.na(CCF_other)) 

wilcox.test(a$CCF_hotspot,a$CCF_other)


# TP53
tmp <- tmpdata2 %>% filter(Gene_Name=='TP53') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:')) %>% count(ID,AAChange) %>% arrange(desc(n)) %>% mutate(Hotspot=if_else(n>10,"Yes",'N')) %>% mutate(Key=if_else(Hotspot=="Yes",AAChange,'Other')) %>% mutate(Key=fct_reorder(Key,n))

PieDonut_ztw2(tmp,aes(pies=Key,count=n),mainCol = c('gray50',pal_d3()(5)),showNum = TRUE,showRatioPie = TRUE,showRatioThreshold = 0,pieLabelSize=3.5,titlesize = 5,r0=0,title='Population', showPieName = FALSE,family = 'Roboto Condensed')
ggsave(filename = 'TP53_hotspot.pdf',width = 8,height = 8)

tmpdata2 %>% filter(Gene_Name=='TP53') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:'),) %>% count(Clone) %>% arrange(desc(n)) %>% janitor::adorn_percentages(denominator = 'col')



# KRAS
tmp <- tmpdata2 %>% filter(Gene_Name=='KRAS') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:')) %>% count(ID,AAChange) %>% arrange(desc(n)) %>% mutate(Hotspot=if_else(n>10,"Yes",'N')) %>% mutate(Key=if_else(Hotspot=="Yes",AAChange,'Other')) %>% mutate(Key=fct_reorder(Key,n))

PieDonut_ztw2(tmp,aes(pies=Key,count=n),mainCol = c('gray50',pal_d3()(5)),showNum = TRUE,showRatioPie = TRUE,showRatioThreshold = 0,pieLabelSize=3.5,titlesize = 5,r0=0,title='Population', showPieName = FALSE,family = 'Roboto Condensed')
ggsave(filename = 'KRAS_hotspot.pdf',width = 8,height = 8)

tmpdata2 %>% filter(Gene_Name=='KRAS') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:')) %>% count(Clone) %>% arrange(desc(n)) %>% janitor::adorn_percentages(denominator = 'col')


tmp2 <- tmpdata2 %>% filter(Gene_Name=='KRAS') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:')) %>% 
  left_join(tmp) %>% 
  select(Tumor_Barcode,AAChange,Key,Clone,CCF) 

left_join(
  tmp2 %>% filter(Key!="Other")  %>% select(Tumor_Barcode,CCF_hotspot=CCF) %>% group_by(Tumor_Barcode) %>% arrange(desc(CCF_hotspot)) %>% slice(1) %>% ungroup(),
  tmp2 %>% filter(Key=="Other")  %>% select(Tumor_Barcode,CCF_other=CCF)%>% group_by(Tumor_Barcode) %>% arrange(CCF_other) %>% slice(1) %>% ungroup()
) %>% 
  filter(!is.na(CCF_hotspot),!is.na(CCF_other)) %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  mutate(name=str_replace(name,'CCF','KRAS')) %>% 
  ggplot(aes(name,value))+geom_boxplot(width=0.5)+
  labs(x="",y="Cancer cell fraction (CCF)\n")+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = FALSE)+
  scale_y_continuous(breaks = pretty_breaks(n = 10),limits = c(0.5,1.5),expand = c(0,0))+
  panel_border(color = "black")

ggsave(filename = 'KRAS_hotspot2.pdf',width = 3.5,height = 7)


# MUC19
tmp <- tmpdata2 %>% filter(Gene_Name=='MUC19') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:')) %>% count(ID,AAChange) %>% arrange(desc(n)) %>% mutate(Hotspot=if_else(n>10,"Yes",'N')) %>% mutate(Key=if_else(Hotspot=="Yes",AAChange,'Other')) %>% mutate(Key=fct_reorder(Key,n))

PieDonut_ztw2(tmp,aes(pies=Key,count=n),mainCol = c('gray50',pal_d3()(5)),showNum = TRUE,showRatioPie = TRUE,showRatioThreshold = 0,pieLabelSize=3.5,titlesize = 5,r0=0,title='Population', showPieName = FALSE,family = 'Roboto Condensed')
ggsave(filename = 'KRAS_hotspot.pdf',width = 8,height = 8)

tmpdata2 %>% filter(Gene_Name=='KRAS') %>% mutate(AAChange=str_remove(AAChange,',.*'),AAChange=str_remove(AAChange,'.*:')) %>% count(Clone) %>% arrange(desc(n)) %>% janitor::adorn_percentages(denominator = 'col')
















# EGFR hotspot  -----------------------------------------------------------

tmp <- tgct_maf %>% filter(Hugo_Symbol=='EGFR') %>% select(Tumor_Barcode,AAChange.refGene) %>% mutate(AAChange.refGene = str_replace(AAChange.refGene,'.*NM_001346898:','')) %>% mutate(AAChange.refGene = str_remove(AAChange.refGene, ',.*')) %>% mutate(AAChange.refGene = str_remove(AAChange.refGene, '.*:')) %>% rename(aaChange=AAChange.refGene) %>% unique()










# Multiple Mutations ------------------------------------------------------

tmp <- tgct_maf %>% 
  filter(Hugo_Symbol %in% c('EGFR','TP53','KRAS')) %>% 
  count(Tumor_Barcode,Hugo_Symbol) %>% 
  filter(n>1) %>% 
  left_join(wgs_groups_info %>% select(Tumor_Barcode,Subject)) %>% 
  mutate(Alteration = 'Multiple_Mutations',Type='Multiple_Mutations') %>% 
  select(Subject,Tumor_Barcode,Gene=Hugo_Symbol,Alteration,Type)


tmp %>% 
  select(Tumor_Barcode,Gene) %>%
  unique() %>%
  left_join(tgct_maf %>% select(Tumor_Barcode,Gene=Hugo_Symbol,POS=Start_Position)) %>% 
  group_by(Tumor_Barcode,Gene) %>%
  arrange(Tumor_Barcode,Gene,desc(POS)) %>%
  mutate(Dist=lag(POS)-POS) %>%
  filter(!is.na(Dist)) %>%
  unique() %>% 
  mutate(Gene = factor(Gene,levels = c('EGFR','TP53','KRAS'))) %>% 
  ggplot(aes(Gene,log2(Dist)))+
  geom_boxplot(width=0.7,outlier.shape = NA)+
  geom_point(pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0),fill='#cccccc')+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'yY')+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  labs(x='',y='Distance between two mutations, log2(bp)')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  guides(fill=guide_legend(override.aes = list(size=2,shape=21)),shape=guide_legend(override.aes = list(size=2)))+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'Multiple_Mutation_Distance.pdf',width = 3, height = 6,device = cairo_pdf)




tmp %>% 
  select(Tumor_Barcode,Gene) %>%
  unique() %>%
  left_join(tgct_maf %>% select(Tumor_Barcode,Gene=Hugo_Symbol,POS=Start_Position,aaChange)) %>% 
  group_by(Tumor_Barcode,Gene) %>% 
  summarise(aaChange=paste0(aaChange,collapse = ',')) %>% View()








# Hotspot frequency between smokers and nonsmokers ------------------------
load('tgct_maf.RData')
load('../BBsolution_final3_short.RData')

tmp <- tgct_maf %>% 
  filter(Hugo_Symbol %in% c('EGFR','KRAS','TP53')) %>% 
  select(Tumor_Barcode,Hugo_Symbol, aaChange) %>% 
  unique() %>% 
  filter(!is.na(aaChange)) %>% 
  group_by(Hugo_Symbol,aaChange) %>%
  mutate(n=n()) %>% 
  ungroup() %>% 
  filter(n>10) %>% 
  arrange(desc(n)) 

tdata <- wgs_groups_info %>% select(Tumor_Barcode, Smoking) %>% count(Smoking,name = 'Total') %>% filter(Smoking !='Unknown')


tmp %>% 
  left_join(
    wgs_groups_info %>% select(Tumor_Barcode, Smoking)
    ) %>% 
  filter(Smoking !='Unknown') %>% 
  count(Smoking,Hugo_Symbol,aaChange,sort=T) %>% 
  left_join(tdata) %>% 
  mutate(Freq=n/Total) %>% 
  select(Hugo_Symbol,aaChange,Smoking,Freq) %>% 
  pivot_wider(names_from = 'Smoking',values_from = 'Freq',values_fill = 0) %>% 
  arrange(Hugo_Symbol)
  
  
  
