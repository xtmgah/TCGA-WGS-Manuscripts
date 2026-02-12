#Fig3----
library(tidyverse)
library(janitor)
library(scales)
library(extrafont)
#font_import()  # Import all fonts from your system
loadfonts()    # Load the fonts
pdfhr()
pdfhr2()
myggstyle()
source('./Rcode/Sherlock_functions.R')
source('./Rcode/Additional_Funcs.R')
source("./Rcode/Sigvisualfunc.R")
#Fig3a----
load('./Rdata/ucec_subtype_data.RData',verbose = T)
load('./Rdata/ucec_mutational_signatures.RData',verbose = T)
# Mutational signature landscape SBS ------------------------------------------
sbs_landscape <- plot_mutational_signature_landscape(
  tumor_activity_all = ucec_activity_all %>% select(Tumor_Barcode,starts_with('SBS')),
  tumor_activity_all_ratio = ucec_activity_all_ratio %>% select(Tumor_Barcode,starts_with('SBS')),
  tumor_decompsite = ucec_sbs96_decompsite,
  sigcol = sigcol_sbs,
  group_data = subtype_data %>% select(Tumor_Barcode,Subtype),
  sortby_signatures = c('SBS2','SBS13'),sortby_cutoff = 0.25,
  legend_sig_nrow = 2,
  output_file = './Rscripts/Figures/Fig3a.png'
)

#Fig3b----
# Mutational signature landscape DBS ------------------------------------------
plot_mutational_signature_landscape(
  tumor_activity_all = ucec_activity_all %>% select(Tumor_Barcode,starts_with('DBS')),
  tumor_activity_all_ratio = ucec_activity_all_ratio %>% select(Tumor_Barcode,starts_with('DBS')),
  tumor_decompsite = ucec_dbs78_decompsite,
  sigcol = sigcol_dbs,
  group_data = subtype_data %>% select(Tumor_Barcode,Subtype),
  sortby_signatures = NULL,
  sortby_cutoff = 0.1 ,
  #samlevs = sbs_landscape$samlevs,
  output_file = './Rscripts/Figures/Fig3b.png'
)

#fig3c----
load("./Rdata/Signature_freqdata.RData")
# by mutations
freqdata %>%
  mutate(Signature = factor(Signature, levels = siglevs)) %>%
  mutate(Type = str_replace(Type, "Prevalence by", "By")) %>%
  mutate(
    lab = if_else(
      Value > 0.05,
      as.character(percent(Value, accuracy = 1)),
      NA_character_
    )
  ) %>%
  mutate(Value = Value * 100) %>%
  mutate(Type = str_remove(Type, 's$')) %>%
  ggplot(aes(Signature, Subtype, fill = Value)) +
  #geom_tile(col='white',linewidth=0.2)+
  geom_point(aes(size = Value), pch = 21, stroke = 0.5) +
  #geom_text(aes(label=lab),col='gray50',size=3)+
  facet_grid(Type ~ Profile, scales = 'free', space = 'free') +
  scale_fill_viridis_c(limit = c(1, 100), breaks = pretty_breaks()) +
  scale_size_area(limit = c(1, 100), breaks = pretty_breaks()) +
  labs(x = NULL, y = NULL, fill = 'Prevalence\n') +
  theme_ipsum_rc(grid = FALSE, plot_margin = margin(5.5, 5.5, 5.5, 5.5)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.spacing = unit(0.2, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 'bold'),
    strip.text.y = element_text(hjust = 0.5, face = 'bold'),
    legend.key.height = unit(0.5, 'cm'),
    legend.position = 'top',
    legend.key.width = unit(2, 'cm')
  ) +
  guides(size = guide_legend(nrow = 1)) +
  panel_border(color = 'black', size = 0.4)

ggsave(
  filename = './Rscripts/Figures/Fig3c.png',
  width = 11.5,
  height = 4
)
