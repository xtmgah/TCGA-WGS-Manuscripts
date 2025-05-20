libztw()
set_wd()
pdfhr2()



# load wgs data info ------------------------------------------------------
load('wgs_data.RData',verbose = T)

# load ecDNA result
ecdna_leves <- c("Linear", "ecDNA","BFB","Complex-non-cyclic")

ecdna <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/ecDNA/ecDNA_default/AA_result_table.tsv',delim = '\t')
ecdna <- ecdna %>% filter(!is.na(Classification)) %>% rename(Barcode=`Sample name`) %>% 
  left_join(
    bind_rows(
      wgs_data %>% select(Subject,Barcode=Tumor_Barcode,Subtype),
      wgs_data %>% select(Subject,Barcode=Normal_Barcode,Subtype) %>% mutate(Subtype='Normal')
    ) %>% unique()
  ) %>% 
  filter(!is.na(Subtype)) %>% 
  relocate(Subject) %>% 
  mutate(Classification = factor(Classification,levels = ecdna_leves)) 


ecdna %>% filter(Subtype == 'Normal') %>% filter(Classification == 'ecDNA') %>% pull(`Captured interval length`)

ecdna <- ecdna %>% filter(Subtype != 'Normal')

## AA count distribution
sub_levs <- ecdna %>% 
  count(Subtype,Subject,Classification) %>% 
  pivot_wider(names_from = Classification,values_from = n,values_fill = 0) %>% 
  mutate(n=`Complex-non-cyclic`+ Linear+ ecDNA + BFB) %>% 
  arrange(desc(n)) %>% 
  pull(Subject)

#class_leves <- ecdna %>% count(Classification,sort=T) %>% pull(Classification)

ecdna %>% 
  count(Subtype,Subject,Classification) %>% 
  mutate(Subject=factor(Subject,levels=sub_levs)) %>% 
  mutate(Classification = factor(Classification,levels = ecdna_leves)) %>% 
  ggplot(aes(Subject,n,fill=Classification))+
  geom_bar(position="stack", stat="identity") +
  #facet_grid(Classification~Subtype, scales = 'free_x',space='free_x')+
  facet_wrap(~Subtype,scales = 'free_x')+
  scale_fill_d3()+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid='XY')+
  theme(axis.text.x = element_text(size = 8,angle = 90,hjust = 1,vjust = 0.5),plot.margin = margin(4,4,4,4),panel.spacing = unit('0.1','cm'),strip.text.y = element_text(vjust = 0,hjust = 0.5),legend.position = 'top',axis.ticks.x = element_blank())+
  labs(x='',y = 'Count')+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'ecdna_amplicon_summary.pdf',width = 8,height = 7.2,device = cairo_pdf)

# ecDNA size

ecdna %>% 
  ggplot(aes(`Captured interval length`/1e6,col=Classification))+
  geom_density(size=1)+
  scale_color_d3()+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = "Y",axis = "XY",axis_col = 'black')+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit('0.1','cm'),strip.text.y = element_text(vjust = 0,hjust = 0.5),legend.position = c(0.86,0.8))+
  labs(x='Amplication captured interval length (Mbp)',y = 'Density')
#panel_border(color = 'black',linetype = 1)
ggsave(filename = 'ecdna_amplicon_size.pdf',width = 6,height = 4,device = cairo_pdf)


#overall frequency
ecdna %>% group_by(Classification) %>% count(Subject) %>% count(Classification) %>% ungroup() %>% mutate(freq=n/252) %>% select(Classification,n,freq) %>% write_clip()

ecdna %>% group_by(Classification,Subtype) %>% count(Subject) %>% count(Classification,Subtype) %>% ungroup() %>%  left_join(wgs_data %>% count(Subtype,name = 'total')) %>% mutate(freq=n/total) %>% select(Subtype,Classification,n,freq) %>% arrange(Subtype,Classification) %>% write_clip()

#with oncogene
ecdna %>% filter(Oncogenes != '[]') %>% group_by(Classification) %>% count(Subject) %>% count(Classification) %>% ungroup() %>% mutate(freq=n/252) %>% select(Classification,freq) %>% write_clip()


# oncogene 

pdata <- ecdna %>% 
  select(Subject,Classification,Oncogenes) %>% 
  mutate(Oncogenes = str_remove_all(Oncogenes,"[\\[\\'\\] ]")) %>% 
  separate_rows(Oncogenes,sep = ',') %>% 
  filter(Oncogenes != '') %>% 
  unique() %>% 
  count(Classification,Oncogenes,sort=T) 

gene_levels <-  pdata %>% group_by(Oncogenes) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% pull(Oncogenes)

pdata %>% 
  mutate(Oncogenes = factor(Oncogenes, levels=gene_levels)) %>% 
  ggplot(aes(Oncogenes,n,fill=Classification)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_d3()+
  scale_y_continuous(breaks = pretty_breaks(n=10))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = "Y",axis = "XY",axis_col = 'black')+
  theme(axis.text.x=element_text(size = 6,angle = 90,vjust = 0.5,hjust = 1),plot.margin = margin(4,4,4,4),panel.spacing = unit('0.1','cm'),strip.text.y = element_text(vjust = 0,hjust = 0.5),legend.position = c(0.86,0.8))+
  labs(x=NULL,y = 'Number of subjects')

ggsave(filename = 'ecdna_amplicon_oncogenes.pdf',width = 10,height = 4,device = cairo_pdf)

pdata <- ecdna %>% 
  filter(Classification == 'ecDNA') %>% 
  mutate(Oncogenes = str_remove_all(Oncogenes,"[\\[\\'\\] ]")) %>% 
  separate_rows(Oncogenes,sep = ',') %>% 
  filter(Oncogenes != '') %>% 
  select(Subtype,Oncogenes,Subject) %>% 
  unique() %>% 
  count(Subtype,Oncogenes,sort=T) 

gene_levs <- pdata %>% group_by(Oncogenes) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% pull(Oncogenes)

pdata %>% 
  #filter(Classification == 'ecDNA') %>% 
  mutate(Oncogenes = factor(Oncogenes,levels = gene_levs)) %>% 
  ggplot(aes(Oncogenes,n,fill=Subtype)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = subtypecol)+
  scale_y_continuous(breaks = pretty_breaks(n=2))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = "Y",axis = "XY",axis_col = 'black')+
  theme(axis.text.x=element_text(size = 6,angle = 90,vjust = 0.5,hjust = 1),plot.margin = margin(4,4,4,4),panel.spacing = unit('0.1','cm'),strip.text.y = element_text(vjust = 0,hjust = 0.5),legend.position = c(0.82,0.75))+
  labs(x=NULL,y = 'Number of subjects')

ggsave(filename = 'ecdna_amplicon_oncogenes_ecDNA.pdf',width = 3,height = 2,device = cairo_pdf)


# copy number association
my_comparisons = list(c('ecDNA','Linear'),c('ecDNA','BFB'),c('ecDNA','Complex-non-cyclic'),c('BFB','Linear'),c('BFB','Complex-non-cyclic'),c('Linear','Complex-non-cyclic'))

ecdna %>% 
  #filter(Type=='Tumor') %>% 
  mutate(Group=if_else(Oncogenes=='[]','Amplification with oncogene','Amplification without oncogene')) %>% 
  ggplot(aes(Classification,log2(`Feature median copy number`),fill=Classification))+
  geom_point(pch=21,size=2,position = position_jitter(width = 0.2,height = 0),stroke=0.1,col='white')+
  geom_boxplot(width=0.6,outlier.shape = NA,fill=NA,col=ncicolpal[1],size=0.6)+
  facet_wrap(~Group)+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.2, "lines"),strip.background = element_rect(fill = 'gray95'),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4))+
  labs(x=NULL, y = 'Median copy number')+
  scale_fill_d3()+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'ecdna_amplicon_copynumber.pdf',width = 6,height = 6,device = cairo_pdf)


my_comparisons = list(c('ecDNA','Linear'),c('ecDNA','BFB'),c('ecDNA','Complex-non-cyclic'),c('BFB','Linear'),c('BFB','Complex-non-cyclic'),c('Linear','Complex-non-cyclic'))

ecdna %>% 
  #filter(Type=='Tumor') %>% 
  mutate(Group=if_else(Oncogenes=='[]','Amplification with oncogene','Amplification without oncogene')) %>% 
  ggplot(aes(Classification,log2(`Feature median copy number`),fill=Classification))+
  geom_point(pch=21,size=2,position = position_jitter(width = 0.2,height = 0),stroke=0.1,col='white')+
  geom_boxplot(width=0.6,outlier.shape = NA,fill=NA,col=ncicolpal[1],size=0.6)+
  #facet_wrap(~Group)+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.2, "lines"),strip.background = element_rect(fill = 'gray95'),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4))+
  labs(x=NULL, y = 'Median copy number')+
  scale_fill_d3()+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons,)

ggsave(filename = 'ecdna_amplicon_copynumber2.pdf',width = 3,height = 6,device = cairo_pdf)



ecdna %>% 
  mutate(Group=if_else(Oncogenes!='[]','Amplification with oncogene','Amplification without oncogene')) %>% 
  select(Barcode,Classification,Group) %>% 
  left_join(wgs_data %>% select(Barcode=Tumor_Barcode,Subtype)) %>% 
  unique() %>% 
  count(Subtype,Classification,Group)


# ecdna %>% 
#   left_join(wgs_data %>% select(Barcode=Tumor_Barcode,Subtype)) %>% 
#   filter(Oncogenes!='[]') %>% 
#   write_csv(file = 'ecdna_oncogene.csv')


ecdna2 <- read_csv(file = 'ecdna_oncogene.csv')

pdata <- ecdna2 %>% 
  select(Subject,Subtype,Classification,Locus) %>% 
  unique() %>% 
  count(Subtype,Classification,Locus,sort=T) 

gene_levels <-  pdata %>% group_by(Locus) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% pull(Locus)

pdata %>% 
  mutate(Locus = factor(Locus, levels=gene_levels)) %>% 
  mutate(Subtype=fct_rev(Subtype)) %>% 
  ggplot(aes(Locus,n,fill=Classification)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~Subtype,scales = 'free',ncol = 1)+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks(n=10))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = "Y",axis = "XY",axis_col = 'black')+
  theme(axis.text.x=element_text(size = 12,angle = 90,vjust = 0.5,hjust = 1),plot.margin = margin(4,4,4,4),panel.spacing = unit('0.1','cm'),strip.text.y = element_text(vjust = 0,hjust = 0.5),legend.position = c(0.8,0.9))+
  labs(x=NULL,y = 'Number of tumors')

ggsave(filename = 'ecdna_amplicon_oncogenes2.pdf',width = 4,height = 7,device = cairo_pdf)

