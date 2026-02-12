#Fig4a----
# Mutational signature landscape ------------------------------------------

source("./Rcode/ztw_function.R")
source('./Rcode/Sherlock_functions.R')
source("./Rcode/ztw.R")
load('./Rdata/ucec_mutational_signatures.RData',verbose = T)
load('./Rdata/wgsdata.RData',verbose = T)
load('./Rdata/ucec_subtype_data.RData',verbose = T)
library(rstatix)
library(scales)
library(ggpubr)
library(hrbrthemes)
# Mutational signature landscape ID ------------------------------------------
plot_mutational_signature_landscape(
  tumor_activity_all = ucec_activity_all %>% select(Tumor_Barcode,starts_with('ID')),
  tumor_activity_all_ratio = ucec_activity_all_ratio %>% select(Tumor_Barcode,starts_with('ID')),
  tumor_decompsite = ucec_id83_decompsite,
  sigcol = sigcol_id,
  group_data = subtype_data %>% select(Tumor_Barcode,Subtype),
  sortby_signatures = NULL,
  sortby_cutoff = 0.1 ,
  #samlevs = sbs_landscape$samlevs,
  output_file = './Rscripts/Figures/Fig4a.png'
)

#Fig4b----
# By Subtype --------------------------------------------------------------
maxsigdata <- ucec_activity_all %>% select(Tumor_Barcode,starts_with('SBS')) %>% pivot_longer(-Tumor_Barcode) %>% group_by(Tumor_Barcode) %>% arrange(desc(value)) %>% slice(1) %>% ungroup() 
ucec_activity_all %>% filter(ID1==0|ID2==0) ## 6 samples 

tdata <- ucec_activity_all %>%
  left_join(maxsigdata) %>% 
  filter(ID1>0,ID2>0) %>% 
  mutate(Ratio=((ID2)/(ID1)),ID1_ID2=ID1+ID2) %>% 
  mutate(Group=if_else(SBS44<50,'No SBS44',if_else(name=='SBS44','SBS44 dominant','SBS44 non-dominant' ))) %>% 
  select(Tumor_Barcode,SBS44,ID1_ID2,Ratio,Group) %>% 
  mutate(Group=factor(Group, levels = c('No SBS44','SBS44 non-dominant','SBS44 dominant'))) 


tdata <- tdata %>% left_join(subtype_data)

stat.test <- tdata %>%
  pairwise_wilcox_test(Ratio ~ Subtype, p.adjust.method = "BH") %>%
  add_xy_position(x = "Subtype") %>%
  #mutate(p.label = p_format(p.adj, digits = 2, add.p = TRUE))  # e.g., "p = 0.003"
  mutate(p.label = if_else(p.adj<0.001,sprintf("FDR = %.2e",p.adj),sprintf("FDR = %.2f",p.adj))) %>% 
  mutate(col=if_else(p.adj<0.05,ncicolpal[1],'black'))

stat.test$y.position <- c(20,20.5,21,21.5,22,22.5)-2

tdata %>% 
  mutate(Ratio = if_else(Ratio>20,20,Ratio)) %>% 
  ggplot(aes(Subtype, Ratio, fill = Subtype)) +
  # jittered points inside boxes, size by SBS44
  geom_point(aes(size = log2(ID1_ID2+1)),
             position = position_jitter(width = 0.1, height = 0, seed = 1),
             shape = 21, color = "white", alpha = 0.5,stroke=0.2) +
  geom_boxplot(width = 0.5, fill = NA, color = "gray20",
               outlier.shape = NA) +
  geom_hline(yintercept = 1,linetype = 2,col='gray',size=0.5)+
  # p-values
  stat_pvalue_manual(stat.test, label = "p.label", col=ncicolpal[1],
                     tip.length = 0.01, step.increase = 0.05, hide.ns = TRUE) +
  labs(x = NULL, y = 'Signature activity ratio: ID2/ID1', size = "ID1+ID2 mutations") +
  scale_fill_manual(values = subtype_color) +
  scale_size_area() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme_ipsum_rc(axis_text_size = 14, axis_title_just = 'm',
                 axis_title_size = 16, grid = 'Y', ticks = TRUE) +
  panel_border(color = 'black') +
  theme(axis.text.x = element_blank(),
        panel.spacing.x = unit(0.4, 'cm'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(4, 4, 4, 4),
        legend.position = 'right',
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 16)) +
  guides(fill=guide_legend(override.aes = list(size=4)),size=guide_legend(override.aes = list(fill='gray20')))
#ylim(c(0,50))
ggsave(filename = './Rscripts/Figures/Fig4b.png',width = 6,height = 6)# separated by Group

#Fig4c----
# ID1 and ID2 Ratio -------------------------------------------------------
library(rstatix)

ucec_activity_all %>% filter(ID1==0|ID2==0) ## 6 samples 

tdata <- ucec_activity_all %>%
  left_join(maxsigdata) %>% 
  filter(ID1>0,ID2>0) %>% 
  mutate(Ratio=((ID2)/(ID1)),ID1_ID2=ID1+ID2) %>% 
  mutate(Group=if_else(SBS44<50,'No SBS44',if_else(name=='SBS44','SBS44 dominant','SBS44 non-dominant' ))) %>% 
  select(Tumor_Barcode,SBS44,ID1_ID2,Ratio,Group) %>% 
  mutate(Group=factor(Group, levels = c('No SBS44','SBS44 non-dominant','SBS44 dominant'))) 

stat.test <- tdata %>%
  pairwise_wilcox_test(Ratio ~ Group, p.adjust.method = "BH") %>%
  add_xy_position(x = "Group") %>%
  #mutate(p.label = p_format(p.adj, digits = 2, add.p = TRUE))  # e.g., "p = 0.003"
  mutate(p.label = if_else(p.adj<0.001,sprintf("FDR = %.2e",p.adj),sprintf("FDR = %.2f",p.adj))) %>% 
  mutate(col=if_else(p.adj<0.05,ncicolpal[1],'black'))

stat.test$y.position <- c(15.5,17.5,18)


tdata %>% 
  ggplot(aes(Group, Ratio, fill = Group)) +
  # jittered points inside boxes, size by SBS44
  geom_point(aes(size = log2(ID1_ID2+1)),
             position = position_jitter(width = 0.1, height = 0, seed = 1),
             shape = 21, color = "white", alpha = 0.6,stroke=0.2) +
  geom_boxplot(width = 0.5, fill = NA, color = "gray20",
               outlier.shape = NA) +
  geom_hline(yintercept = 1,linetype = 2,col='gray',size=0.5)+
  # p-values
  stat_pvalue_manual(stat.test, label = "p.label", col=ncicolpal[1],
                     tip.length = 0.01, step.increase = 0.05, hide.ns = TRUE) +
  labs(x = NULL, y = 'Signature activity ratio: ID2/ID1', size = "ID1+ID2 mutations") +
  scale_fill_manual(values = ncicolpal[c(3:1)]) +
  scale_size_area(max_size = 6) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  theme_ipsum_rc(axis_text_size = 14, axis_title_just = 'm',
                 axis_title_size = 16, grid = 'Y', ticks = TRUE) +
  panel_border(color = 'black') +
  theme(axis.text.x = element_blank(),
        panel.spacing.x = unit(0.4, 'cm'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(4, 4, 4, 4),
        legend.position = 'right',
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 16)) +
  guides(fill=guide_legend(override.aes = list(size=4)),size=guide_legend(override.aes = list(fill='gray20')))
#ylim(c(0,50))

ggsave(filename = './Rscripts/Figures/Fig4c.png',width = 5.5,height = 6)

#fig4d----
tdata <- tdata %>% left_join(subtype_data)
stat.test <- tdata %>%
  pairwise_wilcox_test(Ratio ~ Subtype, p.adjust.method = "BH") %>%
  add_xy_position(x = "Subtype") %>%
  #mutate(p.label = p_format(p.adj, digits = 2, add.p = TRUE))  # e.g., "p = 0.003"
  mutate(p.label = if_else(p.adj<0.001,sprintf("FDR = %.2e",p.adj),sprintf("FDR = %.2f",p.adj))) %>% 
  mutate(col=if_else(p.adj<0.05,ncicolpal[1],'black'))

stat.test$y.position <- c(20,20.5,21,21.5,22,22.5)-2

tdata=tdata %>% 
  mutate(Ratio = if_else(Ratio>20,20,Ratio)) 

tdata2 <- tdata %>% 
  filter(Subtype == 'MSI') %>% 
  left_join(
    ucec_activity_all_ratio %>% mutate(MSI_Signature=SBS44+SBS54+SBS15+SBS21) %>%  select(Tumor_Barcode,SBS44_Ratio=SBS44,MSI_Signature)
  ) 


plot_scatter(
  data = tdata2,
  x = SBS44_Ratio, y = Ratio,
  #group = name,
  #facet = TRUE, facet_scales = "free", palette = ncicolpal[1],
  point_color = "white",
  point_fill = "#542788",
  x_lab = "SBS44 mutation ratio", y_lab = "Signature activity ratio: ID2/ID1",
  label_x_rel = 0.55,
  label_y_rel=0.05,               # one per facet (recycled if lengths differ)
  alpha_sig = 0.05,
  #label_override = ("β = 0.01, P = 6.70e-04"),
  #label_color_override = ncicolpal[1]
  # p_override = c(GroupA = 0.012, GroupB = 0.18),  # optional override by group
  #add_reg_eq = FALSE                      # set TRUE to show the equation
)

ggsave(filename = './Rscripts/Figures/Fig4d.png',width = 5.1,height = 4)

#Fig4f----
# Neoantigene -------------------------------------------------------------

tdata <- read_delim('./Jian/UCEC_Neoaantigen_Only.txt',delim = '\t',col_names = T)

tdata <- tdata %>% 
  left_join(subtype_data %>% select(Subject,Subtype)) %>% 
  filter(!is.na(Subtype)) %>%
  filter(!is.na(Alteration)) %>% 
  mutate(Alteration = log2(Alteration+1)) %>% 
  mutate(Gene = str_remove(Gene,'Neoantigens')) %>% 
  mutate(Gene=fct_rev(Gene))

res_burden2 <- plot_group_compare(
  tdata       = tdata,
  value       = Alteration,
  group       = "Subtype",
  facet =     "Gene",
  FDR_facet = T,
  facet_ncol = 2,
  ylab        = "Neoantigens",
  hide_ns    = T,
  palette     = subtype_color,
  pcol_sig    = ncicolpal[1],
  #facet_scales = 'free_y',
  output_file = "./Rscripts/Figures/fig4f.png",
  #y_pos_adjust = c(3,60,7,5),
  width       = 7,
  height      = 5,
  p_off = F
)
