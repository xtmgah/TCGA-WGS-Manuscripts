#Fig5----
set_wd()
libztw()
library(rstatix)
pdfhr()
pdfhr2()
myggstyle()
source('./Rcode/Sherlock_functions.R')
source("./Rcode/ztw_function.R")
source("./Rcode/ztw.R")

#Fig5a----
# TMB and replication stress --------------------------------------------
load('./Rdata/tcga_repstress_zscore.RData',verbose = T)

pdata <- ucec_tmb %>% 
  mutate(TCGA_Barcode=str_sub(Tumor_Barcode,1,15)) %>% 
  left_join(tcga_repstress_score_cbioportal) %>% 
  filter(!is.na(RepStress_ZScore)) %>% 
  left_join(subtype_data %>% select(Tumor_Barcode,Subtype)) 

pdata %>% group_by(Subtype) %>% do(tidy(cor.test(.$RepStress_ZScore, .$TMB,data=.))) %>% arrange((p.value))

res_depth <- plot_group_compare(
  tdata       = pdata,
  value       = RepStress_ZScore,
  group       = "Subtype",
  ylab        = "Repstress Z-score",
  hide_ns    = TRUE,
  palette     = subtype_color,
  pcol_sig    = ncicolpal[1],
  output_file = "./Rscripts/Figures/Fig5a.png", 
  width       = 4,
  height      = 5,
  y_pos_adjust = c(5, 6, 5.2, 8, 5.8, 6.2)
)


#Fig5b----
# Proliferation score -----------------------------------------------------
tdata <- read_delim('./Jian/UCEC_Proliferation_Only.txt',delim = '\t')

tdata <- tdata %>% 
  left_join(
    subtype_data %>% select(Subject,Subtype)
  ) %>% 
  dplyr::filter(!is.na(Subtype))

res_burden2 <- plot_group_compare(
  tdata       = tdata,
  value       = Alteration,
  group       = "Subtype",
  #facet =     "Gene",
  #FDR_facet = T,
  #facet_ncol = 2,
  ylab        = "Proliferation Score",
  hide_ns    = T,
  palette     = subtype_color,
  pcol_sig    = ncicolpal[1],
  #facet_scales = 'free_y',
  output_file = "./Rscripts/Figures/Fig5b.png",
  #y_pos_adjust = c(3,60,7,5),
  width       = 5,
  height      = 5,
  p_off = F
)

#Fig5c----
load('./Rdata/ucec_RNAseq.RData')
load('./Rdata/ucec_subtype_data.RData')
# XCI: XIST ---------------------------------------------------------------

markers <- c('XIST')

tdata <- subtype_data %>% 
  left_join(
    rdata %>% filter(Hugo_Symbol %in% markers) %>% mutate(RSEM=log2(RSEM+1)) 
  ) #%>% 
#filter(Subtype == 'MSI')


# box plot show the difference

plot_group_compare(
  tdata       = tdata %>% filter(Hugo_Symbol == 'XIST'),
  value       = RSEM,
  group       = "Subtype",
  facet       = "Hugo_Symbol",facet_ncol = 7,
  ylab        = "XIST mRNA expression\nlog2(RNA Seq V2 RSEM + 1)",
  hide_ns    = TRUE,
  palette     = subtype_color,
  pcol_sig    = ncicolpal[1],
  output_file = "./Rscripts/Figures/Fig5c.png", 
  width       = 4.5,
  height      = 5
  #y_pos_adjust = c(80,150)
)

#Fig5d----
source('./Rcode/Sherlock_functions.R')
source("./Rcode/ztw_function.R")
source("./Rcode/ztw.R")
load('./Rdata/ucec_subtype_data.RData',verbose = T)
# BMI ---------------------------------------------------------------------
bmidata <- read_delim('./Jian/merge_ucec_activity_all_obs_association_subtype_data.txt',delim = '\t') %>%
  select(TCGA_ID:Subtype) %>% 
  mutate(Subtype=factor(Subtype,levels = c('POLE','MSI','CN_HIGH','CN_LOW'))) %>% 
  mutate(BMI=as.numeric(BMI)) %>% 
  mutate(BMI_Class = factor(BMI_Class,levels = c('normal','overweight','obese')))


res_burden2 <- plot_group_compare(
  tdata       = bmidata %>% filter(!is.na(BMI)) %>% mutate(BMI=if_else(BMI>200,90,BMI)),
  value       = BMI,
  group       = "Subtype",
  #facet =     "Gene",
  #FDR_facet = T,
  #facet_ncol = 2,
  ylab        = expression(BMI~(kg/m^2)),
  hide_ns    = T,
  palette     = subtype_color,
  pcol_sig    = ncicolpal[1],
  #facet_scales = 'free_y',
  output_file = "./Rscripts/Figures/Fig5d.png",
  y_pos_adjust = c(88,91,95,99,102,106)-25,
  width       = 4.5,
  height      = 5,
  p_off = F
)
#Fig5e----

load('./Rdata/ucec_burden.RData',verbose = T)
load('./Rdata/wgs_clinical.RData',verbose = T)
load('./Rdata/wgs_covdata.RData',verbose = T)
bmidata <- read_delim('./Jian/merge_ucec_activity_all_obs_association_subtype_data.txt',delim = '\t') %>%
  select(TCGA_ID:Subtype) %>% 
  mutate(Subtype=factor(Subtype,levels = c('POLE','MSI','CN_HIGH','CN_LOW'))) %>% 
  mutate(BMI=as.numeric(BMI)) %>% 
  mutate(BMI_Class = factor(BMI_Class,levels = c('normal','overweight','obese')))

tdata <- ucec_burden %>% 
  left_join(bmidata %>% select(Subject,BMI,BMI_Class)) %>% 
  filter(!is.na(BMI)) %>% 
  left_join(ucec_burden) %>% 
  left_join(wgs_covdata) %>% 
  mutate(Subtype=factor(Subtype,levels = c('CN_LOW','CN_HIGH','MSI','POLE')))

tdata %>% group_by(name) %>% do(tidy(lm(value ~ BMI+age_at_diagnosis+Tumor_Purity+Subtype,data=.))) %>% filter(term=='BMI') %>% arrange(p.value)

pdata <- tdata %>% 
  group_by(name) %>%
  do(tidy(lm(value ~ (BMI_Class)+age_at_diagnosis+Tumor_Purity+Subtype,data=.),conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE)) %>% 
  filter(name=="TMB\nlog2(mutations/Mb)") %>%
  #arrange(p.value) %>% 
  ungroup() %>% 
  mutate(term=str_remove(str_remove(term,'Subtype'),'BMI_Class')) %>% 
  filter(term != "(Intercept)") %>%
  mutate(
    sig = p.value < 0.05,
    term = fct_reorder(term, estimate),      # order rows by estimate (or use p.value if you prefer)
    y = as.numeric(term),
    # put label slightly ABOVE the CI line (nudged upward in y)
    y_lab = y + 0.22,
    lab = sprintf("β = %.3f; P = %s", # [%.3f, %.3f]
                  estimate,# conf.low, conf.high,
                  format.pval(p.value, digits = 2, eps = 1e-300))
  )



ggplot(pdata, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  # CI line (colored by significance)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = sig),
                 height = 0, linewidth = 0.9) +
  # Point (colored by significance)
  geom_point(aes(color = sig), size = 2.6) +
  # Text ABOVE the CI line (same significance color)
  geom_text(aes(y = y_lab, label = lab, color = sig),
            hjust = 0,vjust=0, size = 5, lineheight = 0.95) +
  # leave room on the right for text labels
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red"), guide = "none") +
  labs(
    x = "Regression coefficient (β, 95% CI)",
    y = NULL,
    #title = "Forest plot: TMB log2(mutations/Mb) model terms"
  ) +
  theme_classic(base_size = 12,base_family = 'Roboto Consend') +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size=14)
  )

ggsave(filename = './Rscripts/Figures/Fig5e.png',width = 10,height = 5)

#Fig5f----
load('./Rdata/ucec_subtype_data.RData',verbose = T)
tdata %>% 
  mutate(Subtype=factor(Subtype,levels = c('POLE','MSI','CN_HIGH','CN_LOW'))) %>% 
  filter(name=="TMB\nlog2(mutations/Mb)") %>% 
  ggplot(aes(BMI_Class,value,col=Subtype))+
  geom_boxplot(outlier.color = NA)+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_color_manual(values = subtype_color)+
  geom_smooth(method = "lm", se=TRUE, aes(group=1),col=ncicolpal[2])+
  labs(x='BMI Class', y='TMB\nlog2(mutations/Mb)')+
  facet_wrap(~Subtype,scale='free',nrow=1)+
  theme_ipsum_rc(axis_title_just = "m", axis_title_size = 14, ticks = TRUE, grid = FALSE) +
  theme(
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(4, 4, 0, 4),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.position = "none"
  ) +
  ggpubr::border(color = "black", size = 0.35) +
  coord_cartesian(clip = "off")

ggsave(filename = './Rscripts/Figures/Fig5f.png',width = 8,height = 5)
