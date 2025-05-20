set_wd()
libztw()
pdfhr2()

#https://aacrjournals.org/cancerrescommun/article/2/6/503/705014/Replication-Stress-Defines-Distinct-Molecular


# load weight -------------------------------------------------------------
load('repstress_score_weight.RData')
repstress_score_weight


# RepStree Score calculation -------------------------------------------------

tcga_studies <- str_remove_all(list.files(path = '~/NIH-Work/Earl_Stadtman_Inverstigator/Projects/cBioPortal_Studies/TCGA_PanCancer_Atlas/',pattern = '.tar.gz'),'_.*') %>% unique()

#study <- 'tgct'

cdata_all <- NULL
rdata_all <- NULL

for(study in tcga_studies){
  cdata <- NULL
  rdata <- NULL
  
  cdata <- read_delim(paste0('~/NIH-Work/Earl_Stadtman_Inverstigator/Projects/cBioPortal_Studies/TCGA_PanCancer_Atlas/',study,'_tcga_pan_can_atlas_2018/data_clinical_patient.txt'),skip = 4) %>% select(Subject=PATIENT_ID,Sex=SEX)
  
  rdata <-  read_delim(paste0('~/NIH-Work/Earl_Stadtman_Inverstigator/Projects/cBioPortal_Studies/TCGA_PanCancer_Atlas/',study,'_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt'),delim = '\t',col_names = T) %>% 
    pivot_longer(-c(Hugo_Symbol,Entrez_Gene_Id),names_to = 'TCGA_Barcode',values_to = 'RSEM') %>% 
    group_by(TCGA_Barcode) %>% 
    mutate(Z_score = (RSEM - mean(RSEM, na.rm = TRUE)) / sd(RSEM, na.rm = TRUE)) %>% 
    filter(Hugo_Symbol %in% repstress_score_weight$Hugo_Symbol) %>% 
    left_join(repstress_score_weight) %>% 
    summarise(RepStress_Score=sum(Z_score*Weight)) %>% 
    mutate(Study=study)
  
  cdata_all <- bind_rows(cdata_all,cdata)
  rdata_all <- bind_rows(rdata_all,rdata)
}



rdata_all <- rdata_all %>% 
  mutate(RepStress_ZScore = (RepStress_Score - mean(RepStress_Score, na.rm = TRUE)) / sd(RepStress_Score, na.rm = TRUE)) %>% 
  mutate(Study=toupper(Study))
  
tcga_repstress_score_cbioportal <- rdata_all

mvalue <- tcga_repstress_score_cbioportal %>% group_by(Study) %>% summarise(mvalue=median(RepStress_ZScore))
tmp <- tcga_repstress_score_cbioportal %>% select(Study,TCGA_Barcode) %>% unique() %>% count(Study) %>% mutate(StudyLab=paste0(Study,' (',n,')'))
tcga_repstress_score_cbioportal %>% 
  left_join(tmp) %>% 
  left_join(mvalue) %>% 
  mutate(StudyLab=fct_reorder(StudyLab,RepStress_ZScore)) %>% 
  ggplot(aes(StudyLab,RepStress_ZScore,fill=mvalue))+
  geom_boxplot(width=0.7,size=0.25,outlier.stroke = 0.1,outlier.size = 1,outlier.shape = 21,outlier.fill = NA)+
  scale_y_continuous(breaks = pretty_breaks(n=7),limits = c(-5,5))+
  scale_fill_material("deep-orange")+
  #facet_grid(.~Group,scales = 'free_x',space='free')+
  labs(y='Repstress Z-score',x='',color='Group')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,-6,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')


ggsave(filename = 'tcga_repstress_zscore_cbioportal.pdf',width = 5.5,height = 5,device = cairo_pdf)


# by Subtype
#cdata <- read_delim(paste0('~/NIH-Work/Earl_Stadtman_Inverstigator/Projects/cBioPortal_Studies/TCGA_PanCancer_Atlas/','tgct','_tcga_pan_can_atlas_2018/data_clinical_patient.txt'),skip = 4) %>% select(Subject=PATIENT_ID,Sex=SEX)
#cdata <- cdata %>% select(Subject=PATIENT_ID,SUBTYPE) %>% mutate(Group=if_else(SUBTYPE == 'TGCT_seminoma','Seminoma',if_else(SUBTYPE == 'TGCT_non-seminoma','Non-Seminoma',NA_character_))) %>% select(Subject,Subtype=Group)


tdata <- tcga_repstress_score_cbioportal %>% 
  filter(Study=='TGCT') %>% 
  mutate(Subject=str_sub(TCGA_Barcode,1,12)) %>% 
  left_join(wgs_data) %>% 
  filter(!is.na(Subtype))

my_comparisons <-  list(c('Seminoma','Non-Seminoma'))

stat.test <- tdata %>%
  rstatix::wilcox_test(RepStress_ZScore ~ Subtype, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Subtype") %>%
  mutate(myformatted.p = sprintf("P = %.2f",p)) %>% 
  mutate(Subtype = NA)

stat.test

tdata %>% 
  ggplot(aes(Subtype,RepStress_ZScore,fill=Subtype))+
  ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="white",stroke=0.2)+
  geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.6,size=0.4)+
  scale_fill_manual(values = subtypecol)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  labs(y='Repstress Z-score',x='',color='Group')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,-6,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')+
  stat_pvalue_manual(stat.test, label = "myformatted.p",color = ncicolpal[1],size = 3.5)
  
ggsave(filename = 'tcga_repstress_zscore_cbioportal_subtype.pdf',width = 1.6,height = 5,device = cairo_pdf)



# based on XENA result ----------------------------------------------------
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
xena_exp <- read_delim('../../cBioPortal_Studies/XENA/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz',delim = '\t',col_names = T)
colnames(xena_exp)[1] <- 'ENSG'
#gene-exp-counts.deseq2-normalized.log2
tmp <- read_delim('../../cBioPortal_Studies/XENA/repstress_score_weight_genes.txt',delim = '\t',col_names = F) %>% select(X1,X2)
colnames(tmp) <- c('ENSG','Hugo_Symbol')

sexdata <- read_delim('../../cBioPortal_Studies/TCGA_manifest/clinical.cohort.2024-12-23/clinical.tsv') %>% select(Subject=case_submitter_id,Study=project_id,Sex=gender,primary_diagnosis) %>% unique() %>% mutate(Study=str_remove(Study,"TCGA-"))


xena_score <- xena_exp %>% 
  pivot_longer(-ENSG,names_to = 'Barcode',values_to = 'RSEM') %>% 
  group_by(Barcode) %>% 
  mutate(Z_score = (RSEM - mean(RSEM, na.rm = TRUE)) / sd(RSEM, na.rm = TRUE)) 

xena_score <- xena_score %>% 
  filter(ENSG %in% tmp$ENSG) %>% 
  left_join(tmp) %>% 
  left_join(repstress_score_weight)

xena_score <- xena_score %>% 
  summarise(RepStress_Score=sum(Z_score*Weight))

xena_score <- xena_score %>% 
  mutate(Subject=str_sub(Barcode,1,12)) %>% 
  left_join(sexdata)
  
tcga_repstress_score_xena <- xena_score %>% filter(str_starts(Barcode,'TCGA-')) %>%  mutate(RepStress_ZScore = (RepStress_Score - mean(RepStress_Score, na.rm = TRUE)) / sd(RepStress_Score, na.rm = TRUE)) %>% filter(!is.na(Study))


mvalue <- tcga_repstress_score_xena %>% group_by(Study) %>% summarise(mvalue=median(RepStress_ZScore))
tmp <- tcga_repstress_score_xena %>% select(Study,Barcode) %>% unique() %>% count(Study) %>% mutate(StudyLab=paste0(Study,' (',n,')'))
tcga_repstress_score_xena %>% 
  left_join(tmp) %>% 
  left_join(mvalue) %>% 
  mutate(StudyLab=fct_reorder(StudyLab,RepStress_ZScore)) %>% 
  ggplot(aes(StudyLab,RepStress_ZScore,fill=mvalue))+
  geom_boxplot(width=0.7,size=0.25,outlier.stroke = 0.1,outlier.size = 1,outlier.shape = 21,outlier.fill = NA)+
  scale_y_continuous(breaks = pretty_breaks(n=7),limits = c(-5,5))+
  scale_fill_material("deep-orange")+
  #facet_grid(.~Group,scales = 'free_x',space='free')+
  labs(y='Repstress Z-score',x='',color='Group')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,-6,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')


ggsave(filename = 'tcga_repstress_zscore_xena.pdf',width = 5.5,height = 5,device = cairo_pdf)


save(tcga_repstress_score_cbioportal,tcga_repstress_score_xena,xena_score,file = 'tcga_repstress_zscore.RData')
