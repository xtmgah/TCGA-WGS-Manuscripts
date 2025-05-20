set_wd()
libztw()
pdfhr2()

load('wgs_data.RData')
load('~/NIH-Work/R/annotable_ens90.RData',verbose = T)
pc_genelist <- grch38 %>% filter(biotype=='protein_coding') %>% count(symbol,sort=T) %>% filter(n==1, !is.na(symbol),symbol!='.',symbol!='') %>% pull(symbol)

## TCGA data exploring by XENA method
Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)
expdata <- read_delim('../../cBioPortal_Studies/XENA/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz',delim = '\t',col_names = T)
colnames(expdata)[1] <- "ENST"
expdata <- expdata %>% pivot_longer(-ENST,names_to = 'sample',values_to = 'Exp')

geneinfo <- read_delim('../../cBioPortal_Studies/XENA/probeMap_gencode.v23.annotation.gene.probemap',delim = '\t',col_names = T)
colnames(geneinfo)[1] <- "ENST"
geneinfo <- geneinfo %>% filter(gene %in% pc_genelist)

sample2study <- read_delim('../../cBioPortal_Studies/XENA/TCGA_GTEX_category.txt',col_names = T) 
cdata <- read_delim('../../cBioPortal_Studies/TCGA_manifest/clinical.cohort.2024-12-23/clinical.tsv') %>% select(Subject=case_submitter_id,Study=project_id,Sex=gender,primary_diagnosis) %>% unique() %>% mutate(Study=str_remove(Study,"TCGA-")) %>% filter(Study=='TGCT') %>% mutate(Subtype=if_else(is.na(primary_diagnosis)|primary_diagnosis != 'Seminoma, NOS', 'Non-Seminoma','Seminoma')) %>% select(Subject,Subtype,Sex) %>% unique()


# autosome gene 
autoenst <- geneinfo %>% filter(chrom %in% c(paste0('chr',1:22))) %>% pull(ENST)
auto_exp <- expdata %>% filter(ENST %in% autoenst ) %>% group_by(sample) %>% summarise(Auto_mExp=median(Exp,na.rm=T))

chrxenst <- geneinfo %>% filter(chrom %in% 'chrX') %>% pull(ENST)
chrx_exp <- expdata %>% filter(ENST %in% chrxenst ) %>% left_join(geneinfo)


## Compare each gene to median expresison level of autosome genes and test across differnt samples 

tdata <- chrx_exp %>% 
  left_join(auto_exp) %>% 
  mutate(log2FC=log2((2^Exp)/(2^Auto_mExp)))

tdata <- tdata %>% 
  left_join(sample2study) %>% 
  filter(TCGA_GTEX_main_category %in% c('GTEX Testis','TCGA Testicular Germ Cell Tumor')) %>% 
  mutate(Subject=if_else(TCGA_GTEX_main_category == 'TCGA Testicular Germ Cell Tumor',str_sub(sample,1,12),sample)) %>% 
  left_join(cdata) %>% 
  mutate(Group=if_else(TCGA_GTEX_main_category == 'TCGA Testicular Germ Cell Tumor',Subtype,'GTEx Testis'))

pdata <- left_join(
  tdata %>% group_by(Group,gene) %>% do(tidy(wilcox.test(.$Exp,.$uto_mExp,data=.))) %>% ungroup() %>% select(-alternative,-method),
  tdata %>% group_by(Group,gene) %>% summarise(log2FC = median(log2FC)) %>% ungroup()
)


pdata %>% ggplot(aes(log2FC,-log10(p.value)))+geom_point()+facet_wrap(~Group,nrow = 3)


## for each sample, compare the median gene expression between chrx and others and test across different sample

tdata2 <- chrx_exp %>% 
  group_by(sample) %>% 
  summarise(ChrX_mExp=log2(median(2^Exp,na.rm=T))) %>% 
  left_join(auto_exp) %>% 
  mutate(log2FC=log2((2^ChrX_mExp)/(2^Auto_mExp)))

tdata2 <- tdata2 %>% 
  left_join(sample2study) %>% 
  filter(TCGA_GTEX_main_category %in% c('GTEX Testis','TCGA Testicular Germ Cell Tumor')) %>% 
  mutate(Subject=if_else(TCGA_GTEX_main_category == 'TCGA Testicular Germ Cell Tumor',str_sub(sample,1,12),sample)) %>% 
  left_join(cdata) %>% 
  mutate(Group=if_else(TCGA_GTEX_main_category == 'TCGA Testicular Germ Cell Tumor',Subtype,'GTEx Testis'))

my_comparisons <-  list(c('Seminoma','Non-Seminoma'),c('Seminoma','GTEx Testis'),c('GTEx Testis','Non-Seminoma'))

stat.test <- tdata2 %>%
  rstatix::wilcox_test(ChrX_mExp ~ Group, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Group") %>%
  mutate(myformatted.p = sprintf("P = %.2f",p)) %>% 
  mutate(Group = NA)

tdata2 %>% 
  ggplot(aes(Group,ChrX_mExp,fill=Group))+
  ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="white",stroke=0.2)+
  geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.6,size=0.4)+
  scale_fill_manual(values = as.character(c('GTEx Testis' = '#01665e',rev(subtypecol))))+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  labs(y='median chrX mRNA expression\nDESeq2 standardized RSEM (log2)',x='',color='Group')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,-6,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')+
  stat_pvalue_manual(stat.test, label = "myformatted.p",color = ncicolpal[1],size = 3.5)


# DEG between tumor and normal tissues
tdata3 <- chrx_exp
tdata3 <- tdata3 %>% 
  left_join(sample2study) %>% 
  filter(TCGA_GTEX_main_category %in% c('GTEX Testis','TCGA Testicular Germ Cell Tumor')) %>% 
  mutate(Subject=if_else(TCGA_GTEX_main_category == 'TCGA Testicular Germ Cell Tumor',str_sub(sample,1,12),sample)) %>% 
  left_join(cdata) %>% 
  mutate(Group=if_else(TCGA_GTEX_main_category == 'TCGA Testicular Germ Cell Tumor',Subtype,'GTEx Testis'))


pdata1 <- left_join(
  tdata3 %>% filter(Group!='Non-Seminoma') %>% group_by(gene) %>% do(tidy(wilcox.test(Exp~Group,data=.))) %>% select(-method,-alternative) %>% ungroup(),
  tdata3 %>% filter(Group!='Non-Seminoma') %>% group_by(Group,gene) %>% summarise(Exp=median(Exp)) %>% pivot_wider(names_from = Group,values_from = Exp) %>% mutate(log2FC=log2((2^Seminoma)/(2^`GTEx Testis`))) %>% ungroup() %>% rename(TGCT=Seminoma)
) %>% 
  mutate(Group='Seminoma')


pdata2 <- left_join(
  tdata3 %>% filter(Group!='Seminoma') %>% group_by(gene) %>% do(tidy(wilcox.test(Exp~Group,data=.))) %>% select(-method,-alternative) %>% ungroup(),
  tdata3 %>% filter(Group!='Seminoma') %>% group_by(Group,gene) %>% summarise(Exp=median(Exp)) %>% pivot_wider(names_from = Group,values_from = Exp) %>% mutate(log2FC=log2((2^`Non-Seminoma`)/(2^`GTEx Testis`))) %>% ungroup() %>% rename(TGCT=`Non-Seminoma`)
) %>% 
  mutate(Group='Non-Seminoma')

pdata <- bind_rows(pdata1,pdata2)

pdata %>% ggplot(aes(log2FC,-log10(p.value)))+geom_point()+facet_wrap(~Group)

pdata %>% 
  #filter(p.value<0.01) %>% 
  select(gene,log2FC,Group) %>% 
  pivot_wider(values_from = log2FC,names_from = Group) %>% 
  ggplot(aes(Seminoma, `Non-Seminoma`))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)



# Gene expression difference among different SCNA status ------------------
library(valr)
load('scna_chrx.RData',verbose = T)
scna_status <- scna_chrx %>% mutate(Subject=str_sub(Tumor_Barcode,1,12)) %>% select(Subject,chrom,start=startpos,end=endpos,Subtype,nMaj1) %>% mutate(chrom='chrX') 

tdata <- chrx_exp %>% 
  left_join(sample2study) %>% 
  filter(TCGA_GTEX_main_category %in% c('GTEX Testis','TCGA Testicular Germ Cell Tumor')) %>% 
  mutate(Subject=if_else(TCGA_GTEX_main_category == 'TCGA Testicular Germ Cell Tumor',str_sub(sample,1,12),sample)) 

tdata1 <- tdata %>%
  filter(TCGA_GTEX_main_category == 'TCGA Testicular Germ Cell Tumor') %>% 
  rename(start=chromStart,end=chromEnd) %>% 
  bed_intersect(scna_status,suffix = c('','.y')) %>% 
  filter(Subject==Subject.y,abs(start-end)/.overlap >0.8) %>% 
  select(-Subject.y,-start.y,-end.y,-.overlap) %>% 
  rename(Subtype=Subtype.y,nMaj1=nMaj1.y) 

tdata2 <- tdata %>%
  filter(TCGA_GTEX_main_category == 'GTEX Testis') %>% 
  mutate(Subtype='GTEX Testis', nMaj1 = 1)


## by linear regression

pdata1 <- tdata1 %>%
  group_by(gene) %>%
  do(tidy(lm(Exp~nMaj1,data=.))) %>%
  filter(term=='nMaj1') %>%
  ungroup() %>%
  mutate(FDR=p.adjust(p.value,method='bonferroni')) %>%
  arrange(p.value)



pdata1 %>% 
  ggplot(aes(estimate,-log10(FDR)))+
  geom_point(pch=21,size = 2.5,fill = '#71767A',stroke=0.01,col='white')+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[1])+
  ggrepel::geom_text_repel(data = pdata1 %>% filter(FDR<0.05),aes(label=gene),col=ncicolpal[1])+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  labs(x='Linear regression coefficient',x='-log10(FDR)')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'XY')+
  theme(panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')

ggsave(filename = 'chrX_exp_SCNA.pdf',width = 4,height = 4,device = cairo_pdf)


pdata2 <- tdata1 %>% 
  group_by(Subtype,gene) %>% 
  do(tidy(lm(Exp~nMaj1,data=.))) %>% 
  filter(term=='nMaj1') %>% 
  ungroup() %>% 
  group_by(Subtype) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  arrange(p.value)

pdata2 %>% 
  ggplot(aes(estimate,-log10(p.value)))+
  geom_point(pch=21,size = 2.5,fill = '#71767A',stroke=0.01,col='white')+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[1])+
  facet_wrap(~Subtype)+
  ggrepel::geom_text_repel(data = pdata2 %>% filter(p.value<0.05),aes(label=gene),col=ncicolpal[1])+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  labs(x='Linear regression coefficient',x='-log10(p.value)')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'XY')+
  theme(panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')

ggsave(filename = 'chrX_exp_SCNA_Subtype.pdf',width = 8,height = 4,device = cairo_pdf)


tdata1 %>% 
  filter(gene == 'KDM6A') %>% 
  ggplot(aes(as.factor(nMaj1),Exp))+
  geom_boxplot()+
  ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="blue",stroke=0.2)


### by wilcox test

pdata1 <-
  left_join(
    tdata1 %>% 
      mutate(Group=if_else(nMaj1>1,"Yes","No")) %>% 
      group_by(gene) %>% 
      do(tidy(wilcox.test(Exp~Group,data=.))) %>% 
      ungroup() %>% 
      mutate(FDR=p.adjust(p.value,method='bonferroni')) %>% 
      arrange(p.value),
    
    tdata1 %>% 
      mutate(Group=if_else(nMaj1>1,"Yes","No")) %>% 
      group_by(gene,Group) %>% 
      summarise(Exp=median(2^Exp)) %>% 
      ungroup() %>% 
      pivot_wider(names_from = Group,values_from = Exp) %>% 
      mutate(log2FC=log2(Yes/No))
  )

pdata1 %>% 
  ggplot(aes(log2FC,-log10(FDR)))+
  geom_point(pch=21,size = 2.5,fill = '#71767A',stroke=0.1,col='white')+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[1])+
  ggrepel::geom_text_repel(data = pdata1 %>% filter(FDR<0.05),aes(label=gene),col=ncicolpal[1])+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  labs(x='log2(Fold Change)',x='-log10(FDR)')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'XY')+
  theme(panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')

ggsave(filename = 'chrX_exp_SCNA_wilcox.pdf',width = 4,height = 4,device = cairo_pdf)

