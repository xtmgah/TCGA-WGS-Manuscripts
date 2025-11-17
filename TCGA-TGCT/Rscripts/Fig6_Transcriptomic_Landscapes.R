# ------------------------------------------------------------------------------
# Script: Figure 6 - TCGA-TGCT WGS analysis
# Description: This script reproduces all main analyses and figures for the follwoing TGCT manuscript: 
# A Thiopurine-like Mutagenic Process Defines TGCT Subtypes
# ------------------------------------------------------------------------------

# --- Load Required Libraries and Set Plotting Styles ---------------------------
# (You may need to install some packages if missing)

source('../functions/Sherlock_functions.R')

set_wd()
libztw()
pdfhr()
pdfhr2()



# Fig. 6a -----------------------------------------------------------------
load('../data/ssGSEA.RData')
load('../data/wgs_data.RData')

pdata <- ssgsea_results_h %>%
  left_join(barcode2study) %>%
  left_join(
    wgs_data %>%
      mutate(TCGA_Barcode = str_sub(Tumor_Barcode, 1, 15)) %>%
      select(TCGA_Barcode, Subtype)
  ) %>%
  mutate(
    Group = if_else(
      Study == 'TGCT',
      if_else(is.na(Subtype), 'TGCT', Subtype),
      'Others'
    )
  )


## comparing other TCGA tumors as a whole
pdata2 <- pdata %>%
  mutate(Group = if_else(Study == 'TGCT', 'TGCT', 'Others'))

pdata3 <- pdata2 %>%
  group_by(Pathway) %>%
  do(tidy(wilcox.test(ssGSEA_socre ~ Group, data = .))) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

# #if adjusted by purity
# pdata3 <- pdata2 %>%
#   left_join(pureedata %>% rename(TCGA_Barcode=Barcode)) %>% filter(!is.na(Purity)) %>%
#   group_by(Pathway) %>%
#   do(tidy(lm(ssGSEA_socre ~ Group + Purity, data = .))) %>%
#   ungroup() %>%
#   filter(term=='GroupTGCT') %>%
#   mutate(FDR = p.adjust(p.value, method = 'BH'))
# 

tmp <- pdata2 %>%
  group_by(Pathway, Group) %>%
  summarise(score = median(ssGSEA_socre)) %>%
  ungroup() %>%
  pivot_wider(names_from = Group, values_from = score) %>%
  mutate(FC = TGCT - Others)

pdata3 <- pdata3 %>%
  left_join(tmp) %>%
  mutate(Pathway = str_remove(Pathway, '^HALLMARK_')) %>%
  arrange(FC)

jitter_pos1 <- position_jitter(width = .1, height = 0, seed = 1)
jitter_pos2 <- position_jitternudge(width = .1, height = 0, seed = 1, x = 0.2)

library(ggpp)
library(ggrepel)

pdata3 %>%
  ggplot(aes(x = 1, y = FC)) +
  geom_hline(yintercept = 0, col = 'gray20') +
  geom_point(
    aes(size = -log10(FDR), fill = TGCT),
    pch = 21,
    position = jitter_pos1
  ) +
  scale_x_continuous(limits = c(0.9, 2)) +
  scale_fill_gradient2(
    midpoint = 0,
    low = 'blue',
    high = 'red',
    mid = 'white'
  ) +
  #scale_size_binned_area()+
  labs(
    x = NULL,
    y = 'ssGSEA score difference (TGCT - Others)',
    fill = 'TGCT ssGSEA score'
  ) +
  theme_ipsum_rc(
    base_size = 18,
    axis = F,
    grid = F,
    axis_title_just = 'm',
    axis_text_size = 16,
    axis_title_size = 18
  ) +
  theme(
    plot.margin = margin(4, 4, 4, 4),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(size = 0.5, colour = 'gray20'),
    axis.ticks.y = element_line()
  ) +
  geom_text_repel(
    data = pdata3 %>%
      arrange(desc(abs(FC))) %>%
      group_by(FC > 0) %>%
      slice(1:10),
    aes(label = Pathway),
    color = 'black', #color=if_else(FC>0,'red4','blue4')
    position = jitter_pos2,
    force = 2,
    #nudge_x           = 0.2,
    direction = "y",
    hjust = 0,
    segment.size = 0.2,
    segment.color = 'gray',
    segment.curvature = -0.1
  )

#ggsave('ssGSEA_TCGA_HALLMARK.pdf', width = 8, height = 8, device = cairo_pdf)



# Fig. 6b -----------------------------------------------------------------

pdata <- ssgsea_results_kegg %>%
  left_join(barcode2study) %>%
  left_join(
    wgs_data %>%
      mutate(TCGA_Barcode = str_sub(Tumor_Barcode, 1, 15)) %>%
      select(TCGA_Barcode, Subtype)
  ) %>%
  mutate(
    Group = if_else(
      Study == 'TGCT',
      if_else(is.na(Subtype), 'TGCT', Subtype),
      'Others'
    )
  )


## comparing other TCGA tumors as a whole
pdata2 <- pdata %>%
  mutate(Group = if_else(Study == 'TGCT', 'TGCT', 'Others'))

pdata3 <- pdata2 %>%
  group_by(Pathway) %>%
  do(tidy(wilcox.test(ssGSEA_socre ~ Group, data = .))) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

# #if adjusted by purity
# pdata3 <- pdata2 %>% 
#   left_join(pureedata %>% rename(TCGA_Barcode=Barcode)) %>% filter(!is.na(Purity)) %>% 
#   group_by(Pathway) %>%
#   do(tidy(lm(ssGSEA_socre ~ Group + Purity, data = .))) %>%
#   ungroup() %>%
#   filter(term=='GroupTGCT') %>% 
#   mutate(FDR = p.adjust(p.value, method = 'BH'))
# 

tmp <- pdata2 %>%
  group_by(Pathway, Group) %>%
  summarise(score = median(ssGSEA_socre)) %>%
  ungroup() %>%
  pivot_wider(names_from = Group, values_from = score) %>%
  mutate(FC = TGCT - Others)

pdata3 <- pdata3 %>%
  left_join(tmp) %>%
  mutate(Pathway = str_remove(Pathway, '^KEGG_')) %>%
  arrange(FC)

jitter_pos1 <- position_jitter(width = .1, height = 0, seed = 1)
jitter_pos2 <- position_jitternudge(width = .1, height = 0, seed = 1, x = 0.2)

library(ggpp)
library(ggrepel)

pdata3 %>%
  ggplot(aes(x = 1, y = FC)) +
  geom_hline(yintercept = 0, col = 'gray20') +
  geom_point(
    aes(size = -log10(FDR), fill = TGCT),
    pch = 21,
    position = jitter_pos1
  ) +
  scale_x_continuous(limits = c(0.9, 3)) +
  scale_fill_gradient2(
    midpoint = 0,
    low = 'blue',
    high = 'red',
    mid = 'white'
  ) +
  #scale_size_binned_area()+
  labs(
    x = NULL,
    y = 'ssGSEA score difference (TGCT - Others)',
    fill = 'TGCT ssGSEA score'
  ) +
  theme_ipsum_rc(
    base_size = 18,
    axis = F,
    grid = F,
    axis_title_just = 'm',
    axis_text_size = 16,
    axis_title_size = 18
  ) +
  theme(
    plot.margin = margin(4, 4, 4, 4),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(size = 0.5, colour = 'gray20'),
    axis.ticks.y = element_line()
  ) +
  geom_text_repel(
    data = pdata3 %>%
      arrange(desc(abs(FC))) %>%
      group_by(FC > 0) %>%
      slice(1:10),
    aes(label = Pathway),
    color = 'black', #color=if_else(FC>0,'red4','blue4')
    position = jitter_pos2,
    force = 2,
    #nudge_x           = 0.2,
    direction = "y",
    hjust = 0,
    segment.size = 0.2,
    segment.color = 'gray',
    segment.curvature = -0.1
  )

#ggsave('ssGSEA_TCGA_KEGG.pdf', width = 10, height = 8, device = cairo_pdf)



# Fig. 6c -----------------------------------------------------------------

pdata <- ssgsea_results_h %>%
  left_join(barcode2study) %>%
  left_join(
    wgs_data %>%
      mutate(TCGA_Barcode = str_sub(Tumor_Barcode, 1, 15)) %>%
      select(TCGA_Barcode, Subtype)
  ) %>%
  mutate(
    Group = if_else(
      Study == 'TGCT',
      if_else(is.na(Subtype), 'TGCT', Subtype),
      'Others'
    )
  )

pdata2 <- pdata %>%
  filter(
    Pathway %in%
      c(
        "HALLMARK_MYC_TARGETS_V2",
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_G2M_CHECKPOINT",
        "HALLMARK_SPERMATOGENESIS",
        "HALLMARK_MITOTIC_SPINDLE",
        "HALLMARK_DNA_REPAIR"
      )
  ) %>%
  group_by(Study, Pathway) %>%
  summarise(score = median(ssGSEA_socre)) %>%
  ungroup() %>%
  arrange(score) %>%
  mutate(Study = fct_inorder(Study)) %>%
  mutate(Pathway = str_remove(Pathway, '^HALLMARK_'))

# polar plot

ggplot(
  data = pdata2,
  aes(x = (Study), y = score, group = Pathway, color = Pathway)
) +
  #ylim(0, NA)+
  geom_point(stat = 'identity') +
  geom_line(linewidth = 0.5) +
  coord_polar(start = -pi * 1 / 24 + 0.2) +
  theme_ipsum_rc(base_size = 11, ticks = FALSE, grid = "Xx") +
  scale_color_npg() +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_line(linetype = 2)
  ) +
  geom_hline(
    yintercept = seq(-0.15, 0.45, 0.1),
    color = "gray90",
    size = 0.15
  ) +
  annotate(
    'text',
    x = 1,
    y = seq(-0.15, 0.45, 0.1),
    label = seq(-0.15, 0.45, 0.1),
    size = 3,
    col = '#542788'
  ) +
  #guides(color="none")+
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# ggsave(
#   'ssGSEA_TCGA_HALLMARK_examples.pdf',
#   width = 8,
#   height = 6,
#   device = cairo_pdf
# )



# Fig. 6d -----------------------------------------------------------------
pdata <- ssgsea_results_kegg %>%
  left_join(barcode2study) %>%
  left_join(
    wgs_data %>%
      mutate(TCGA_Barcode = str_sub(Tumor_Barcode, 1, 15)) %>%
      select(TCGA_Barcode, Subtype)
  ) %>%
  mutate(
    Group = if_else(
      Study == 'TGCT',
      if_else(is.na(Subtype), 'TGCT', Subtype),
      'Others'
    )
  )


pdata2 <- pdata %>%
  filter(
    Pathway %in%
      c(
        'KEGG_BASE_EXCISION_REPAIR',
        'KEGG_HOMOLOGOUS_RECOMBINATION',
        'KEGG_MISMATCH_REPAIR',
        'KEGG_DNA_REPLICATION'
      )
  ) %>%
  group_by(Study, Pathway) %>%
  summarise(score = median(ssGSEA_socre)) %>%
  ungroup() %>%
  arrange(score) %>%
  mutate(Study = fct_inorder(Study)) %>%
  mutate(Pathway = str_remove(Pathway, '^KEGG_'))

# polar plot

ggplot(
  data = pdata2,
  aes(x = (Study), y = score, group = Pathway, color = Pathway)
) +
  #ylim(0, NA)+
  geom_point(stat = 'identity') +
  geom_line(linewidth = 0.5) +
  coord_polar(start = -pi * 1 / 24 + 0.42, ) +
  theme_ipsum_rc(base_size = 11, ticks = FALSE, grid = "Xx") +
  scale_color_bmj() +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_line(linetype = 2)
  ) +
  geom_hline(
    yintercept = c(0, 0.1, 0.2, 0.3, 0.35),
    color = "gray90",
    size = 0.15
  ) +
  annotate(
    'text',
    x = 1,
    y = c(0, 0.1, 0.2, 0.3),
    label = c(0, 0.1, 0.2, 0.3),
    size = 3,
    col = '#542788'
  ) +
  #guides(color="none")+
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
# 
# ggsave(
#   'ssGSEA_TCGA_KEGG_examples.pdf',
#   width = 8,
#   height = 6,
#   device = cairo_pdf
# )



# Fig. 6e -----------------------------------------------------------------
load('../data/tcga_expdata_all.RData',verbose = T)

pdata <- expdata %>% filter(Gene=='UNG') 
mvalue <- pdata %>% group_by(Study) %>% summarise(mvalue=median(log2(RSEM+1)))
tmp <- pdata %>% select(Study,TCGA_Barcode) %>% unique() %>% count(Study) %>% mutate(StudyLab=paste0(Study,' (',n,')'))

p1 <- pdata %>% 
  left_join(tmp) %>% 
  left_join(mvalue) %>% 
  mutate(StudyLab=fct_reorder(StudyLab,RSEM)) %>% 
  ggplot(aes(StudyLab,log2(RSEM+1),fill=mvalue))+
  geom_boxplot(outlier.shape = NA,width=0.7,size=0.25)+
  scale_y_continuous(breaks = pretty_breaks(n=7),limits = c(7.5,13.5))+
  scale_fill_viridis_c()+
  labs(y='UNG mRNA expression\nlog2(RNA Seq V2 RSEM + 1)',x='',color='Group')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,-6,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')

# by subtype

pdata2 <- pdata %>% 
  filter(Study=='TGCT') %>% 
  left_join(wgs_data %>% select(Subtype,Subject)) %>% 
  mutate(Study=Subtype) %>% 
  filter(!is.na(Subtype))

my_comparisons <-  list(c('Seminoma','Non-Seminoma'))

stat.test <- pdata2 %>% 
  mutate(RSEM=log2(RSEM+1)) %>% 
  rstatix::wilcox_test(RSEM ~ Subtype, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Subtype") %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(SubtypeLab = NA)

stat.test

tmp <- pdata2 %>% select(Subtype,TCGA_Barcode) %>% unique() %>% count(Subtype) %>% mutate(SubtypeLab=paste0(Subtype,' (',n,')'))

p2 <- pdata2 %>% 
  left_join(tmp) %>% 
  ggplot(aes(SubtypeLab,log2(RSEM+1),fill=SubtypeLab))+
  ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="white",stroke=0.2)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.4)+
  scale_fill_manual(values = rev(as.character(subtypecol)))+
  scale_y_continuous(breaks = pretty_breaks(n = 7),limits = c(7.5,13.5))+
  labs(y=NULL,x='',color='Group')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,-6,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')+
  stat_pvalue_manual(stat.test, label = "myformatted.p",color = ncicolpal[1])

plot_grid(p1,p2,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(5,1))

#ggsave(filename = 'tcga_UNG_exp_tcga_only.pdf',width = 8,height = 6,device = cairo_pdf)



# Fig. 6f -----------------------------------------------------------------
load('../data/tcga_expdata_all.RData',verbose = T)

apobec_genes <- c('AICDA','APOBEC2','APOBEC3A','APOBEC3B','APOBEC3C','APOBEC3D','APOBEC3E','APOBEC3F','APOBEC3G','APOBEC3H','APOBEC4')
tdata <- expdata %>% 
  filter(Gene %in% apobec_genes) %>% 
  left_join(
    expdata %>% filter(Gene == 'UNG') %>% select(TCGA_Barcode,Subject,Study,UNG=RSEM)
  )

pdata <- tdata %>% 
  group_by(Study,Gene) %>% 
  do(tidy(cor.test(log2(.$UNG+1),log2(.$RSEM+1)))) %>% 
  arrange(p.value) %>% 
  group_by(Gene) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  arrange(FDR)

pdata %>% 
  filter(Gene=='APOBEC3B') %>% 
  ggplot(aes(estimate,-log10(FDR),fill=Gene))+
  geom_point(pch=21,stroke=0.2,size=3)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='red',size=0.5)+
  ggrepel::geom_text_repel(aes(label=Study))+
  scale_fill_jco()+
  scale_y_sqrt()+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid="XY",ticks = T)+
  panel_border(color = 'black',size = 0.5)+
  #theme(legend.position = 'none')+
  labs(x='Pearson correlation coefficient', y = '-log10(FDR)')


pdata %>% 
  filter(Gene=='APOBEC3B') %>% 
  mutate(FDR=if_else(FDR<0.05,FDR,NA_real_)) %>% 
  mutate(Study=fct_rev(fct_reorder(Study,estimate))) %>% 
  ggplot(aes(Study,estimate,fill=-log10(FDR),group=1))+
  geom_hline(yintercept = 0,col='gray40',linewidth=0.2)+
  geom_line(colour = '#cccccc',linetype=2)+
  geom_point(pch=21,stroke=0.3,size=4)+
  scale_y_continuous(breaks = pretty_breaks(n = 8))+
  ggrepel::geom_text_repel(aes(label=Study),size=2.5)+
  scale_fill_viridis_c()+
  labs(x=NULL,y='Pearson correlation coefficient',fill="-log10(FDR)\n")+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = F)+
  theme(axis.text.x = element_blank(),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4),axis.ticks.x = element_blank())+
  theme(legend.position = 'top',legend.key.width = unit(1.5,'cm'),legend.key.height = unit(0.3,'cm'))+
  panel_border(color = 'black',size=0.5)+
  theme(legend.position = 'top')

#ggsave(filename = 'tcga_APOBEC3B_UNG_tcga_only.pdf',width = 5,height = 5,device = cairo_pdf)


# Fig. 6g -----------------------------------------------------------------

tdata <- expdata %>% 
  filter(Gene %in% apobec_genes) %>% 
  left_join(
    expdata %>% filter(Gene == 'UNG') %>% select(TCGA_Barcode,Subject,Study,UNG=RSEM)
  )
tdata %>% 
  filter(Gene=='APOBEC3B',Study=='TGCT') %>% 
  left_join(wgs_data) %>% 
  filter(!is.na(Subtype)) %>% 
  ggplot(aes(log2(RSEM+1),log2(UNG+1)))+
  geom_point(pch=21,size=2.5,aes(fill=Subtype))+
  stat_smooth(method = "lm") +
  stat_cor(size=5,col=ncicolpal[1],label.y =10.5) +
  scale_fill_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  labs(x='APOBEC3B mRNA expression\nlog2(RNA Seq V2 RSEM + 1)',y='UNG mRNA expression\nlog2(RNA Seq V2 RSEM + 1)')+
  theme_ipsum_rc(axis_text_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'XY',ticks = T)+
  panel_border(color = 'black')+
  theme(#axis.text.x = element_blank(),
    panel.spacing = unit(0.4,'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4,4,4,4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(hjust = 0.5,face = 'bold',size = 16)
  )+
  guides(fill='none',col='none')

#ggsave(filename = 'tgct_APOBEC3B_UNG_tcga_only.pdf',width = 6,height = 5,device = cairo_pdf)



