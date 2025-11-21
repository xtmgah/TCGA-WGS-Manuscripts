# ------------------------------------------------------------------------------
# Script: Figure 5 - TCGA-TGCT WGS analysis
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


# Fig. 5a -----------------------------------------------------------------
# Number of Subclones -----------------------------------------------------
load('../data/DP_info_data.RData', verbose = T)
load('../data/ngspurity_data.RData')

library(plotthis)

pdata <- Cluster_info %>%
  filter(Clone == "N") %>%
  group_by(Tumor_Barcode) %>%
  summarise(N = n(), CCF = min(Cluster_CCF)) %>%
  right_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  mutate(CCF = if_else(is.na(N), 0, CCF)) %>%
  mutate(N = if_else(is.na(N), 0, N)) %>%
  arrange(CCF) %>%
  group_by(Subtype) %>%
  mutate(seq = seq_along(CCF)) %>%
  ungroup()

pdata <- pdata %>%
  mutate(
    CCF_group = case_when(
      CCF == 0 ~ "0",
      CCF > 0 & CCF <= 0.05 ~ "(0,0.05]",
      CCF > 0.05 & CCF <= 0.1 ~ "(0.05,0.1]",
      CCF > 0.1 & CCF <= 0.2 ~ "(0.1,0.2]",
      CCF > 0.2 & CCF <= 0.3 ~ "(0.2,0.3]",
      CCF > 0.3 & CCF <= 0.4 ~ "(0.3,0.4]",
      CCF > 0.4 & CCF <= 0.5 ~ "(0.4,0.5]",
      CCF > 0.5 & CCF <= 0.6 ~ "(0.5,0.6]",
      CCF > 0.6 & CCF <= 0.7 ~ "(0.6,0.7]",
      CCF > 0.7 & CCF <= 0.8 ~ "(0.7,0.8]",
      CCF > 0.8 ~ "(0.8,0.9]",
      TRUE ~ "Other" # Catch-all for any unexpected values
    ),
    # Optional: Make it a factor with ordered levels
    CCF_group = factor(
      CCF_group,
      levels = c(
        "0",
        "(0,0.05]",
        "(0.05,0.1]",
        "(0.1,0.2]",
        "(0.2,0.3]",
        "(0.3,0.4]",
        "(0.4,0.5]",
        "(0.5,0.6]",
        "(0.6,0.7]",
        "(0.7,0.8]",
        "(0.8,0.9]",
        "Other"
      )
    )
  )

pdata %>%
  count(CCF_group, Subtype) %>%
  mutate(Subtype = fct_rev(Subtype)) %>%
  RingPlot(
    x = 'Subtype',
    y = 'n',
    group_by = 'CCF_group',
    group_by_sep = "A",
    group_name = 'Minimal CCF of detected subclones'
  ) +
  scale_fill_viridis_d()

# ggsave(
#   'subclone_miniccf_subtype.pdf',
#   width = 6,
#   height = 6,
#   device = cairo_pdf()
# )


# Fig. 5b -----------------------------------------------------------------
load('../data/clsdata.RData',verbose = T)

clsdata2 <- clsdata %>%
  mutate(CLS = if_else(str_starts(CLS, 'clonal'), 'clonal', CLS)) %>%
  group_by(Tumor_Barcode, CLS) %>%
  summarise(n = sum(n), Freq = sum(Freq)) %>%
  ungroup()


clsdata2 %>%
  left_join(wgs_data) %>%
  filter(!is.na(CLS)) %>%
  ggplot(aes(CLS, Freq, fill = Subtype)) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.color = 'white',
    outlier.size = 2,
    outlier.stroke = 0.2
  ) +
  geom_vline(xintercept = c(1.5), linewidth = 0.25, col = '#cccccc') +
  labs(x = "", y = "Proportion of total mutations", fill = 'Subtype') +
  scale_fill_manual(values = subtypecol, breaks = names(subtypecol)) +
  theme_ipsum_rc(
    axis_text_size = 12,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = F
  ) +
  theme(
    axis.ticks.x = element_line(linewidth = .5),
    plot.margin = margin(4, 4, 4, 4)
  ) +
  #theme(legend.position = 'top',legend.direction = 'horizontal')+
  panel_border(color = 'black')

# ggsave(
#   filename = 'CLS_subtype2.pdf',
#   width = 3.5,
#   height = 5,
#   device = cairo_pdf
# )



# Fig. 5c -----------------------------------------------------------------
load('../data/MCN_data.RData')
library(ggtern)
MCN_data %>%
  mutate(MCN_R1 = 1 - MCN_R2 - MCN_R3) %>%
  ggtern(ggtern::aes(MCN_R1, MCN_R2, MCN_R3, fill = Subtype, shape = WGD)) +
  geom_mask() +
  geom_point(size = 3) +
  scale_shape_manual(values = c(24, 21)) +
  scale_fill_manual(values = ncicolpal[c(4, 5)]) +
  ggtern::theme_rgbw(base_size = 14, base_family = 'Roboto') +
  theme_nogrid() +
  theme_hidetitles() +
  theme_clockwise() +
  theme_legend_position('tr') +
  geom_Lline(Lintercept = 0.5) +
  geom_Rline(Rintercept = 0.5) +
  labs(
    xarrow = '% autosomal genome with MCN <2',
    yarrow = '% autosomal genome with MCN >=2 & <3',
    zarrow = '% autosomal genome with MCN >=3'
  )

#ggtern::ggsave(filename = 'Tumor_WGD_tmp.pdf', width = 6.5, height = 6.5)


# Fig. 5d -----------------------------------------------------------------
sample_wgd2 <- MCN_data %>% filter(MCN_R3 > 0.5) %>% pull(Tumor_Barcode)
load('tgct_signature_probability.RData', verbose = T)

DP_info_data <- DP_info_data %>% left_join(tgct_sbs96_probability)

tmp <- DP_info_data %>%
  filter(no.chrs.bearing.mut %in% c(1, 2, 3)) %>%
  rename(CNV_WGD = no.chrs.bearing.mut) %>%
  group_by(Tumor_Barcode, CNV_WGD) %>%
  summarise(CNV = sum(SBS1 + SBS5, na.rm = T))

tmp <- DP_info_data %>%
  count(Tumor_Barcode, name = 'Total') %>%
  left_join(tmp)
tmp <- tmp %>% replace_na(list(CNV = 0))

tmp <- ngspurity %>%
  select(Tumor_Barcode, WGD_Status = MCN_WGD) %>%
  left_join(tmp) %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  mutate(Ratio = CNV / Total)


# tidyHeatmap::heatmap(
#   tmp,.row = Tumor_Barcode,.column = CNV_WGD,.value = Ratio
# )

tmp1 <- left_join(
  tmp %>%
    filter(CNV_WGD == 1) %>%
    select(Tumor_Barcode, Subtype, Total, Ratio1 = Ratio),
  tmp %>%
    filter(WGD_Status == 'WGD', CNV_WGD == 2) %>%
    select(Tumor_Barcode, Ratio2 = Ratio) %>%
    mutate(Group = 'WGD')
) %>%
  filter(!is.na(Group)) %>%
  mutate(R = Ratio2 / Ratio1, Subratio = Ratio1 + Ratio2) %>%
  select(Group, Subtype, Tumor_Barcode, Subratio, Total, R)

tmp2 <- left_join(
  tmp %>%
    filter(CNV_WGD == 2) %>%
    select(Tumor_Barcode, Subtype, Total, Ratio1 = Ratio),
  tmp %>%
    filter(
      WGD_Status == 'WGD',
      CNV_WGD == 3,
      Tumor_Barcode %in% sample_wgd2
    ) %>%
    select(Tumor_Barcode, Ratio2 = Ratio) %>%
    mutate(Group = 'WGD2')
) %>%
  filter(!is.na(Group)) %>%
  mutate(R = Ratio2 / Ratio1, Subratio = Ratio1 + Ratio2) %>%
  select(Group, Subtype, Tumor_Barcode, Subratio, Total, R)


pdata <- bind_rows(tmp1, tmp2)

pdata %>%
  ggplot(aes(Group, R, fill = Subtype)) +
  #geom_violin(trim=T,size=0.2)+
  geom_boxplot(width = 1, color = 'black', outlier.shape = NA, alpha = 0.2) +
  ggbeeswarm::geom_quasirandom(
    pch = 21,
    size = 2,
    width = 0.2,
    color = "white",
    stroke = 0.2,
    dodge.width = 1
  ) +
  geom_hline(yintercept = 1, linetype = 2, col = ncicolpal[2]) +
  geom_vline(xintercept = 1.5, col = '#cccccc', size = 0.2) +
  labs(
    x = NULL,
    y = '(mutations with copy number = a)/(mutations with copy number = b)\n(WGD: a=2, b=1; WGD2: a=3, b=2)'
  ) +
  scale_fill_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 7)) +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    grid = F,
    ticks = T
  ) +
  panel_border(color = 'black') +
  theme(
    #axis.text.x = element_blank(),
    panel.spacing.x = unit(0.4, 'cm'),
    axis.ticks.y = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  scale_color_manual(values = c(ncicolpal[1], 'black')) +
  guides(color = 'none') +
  coord_flip()

# ggsave(
#   filename = 'WGD_timing_subtypes_SBS15.pdf',
#   width = 8,
#   height = 4,
#   device = cairo_pdf
# )


# Fig. 5e -----------------------------------------------------------------

# WGS timming and clonal structure/drivers  -------------------------------
# WGD timming -------------------------------------------------------------
load('../data/wgs_data.RData')
load('../data/ngspurity_data.RData')
load('../data/DP_info_data.RData', verbose = T)
load('../data/clsdata.RData')
#load('../ZTW_functions.RData')
tmp <- DP_info_data %>%
  filter(no.chrs.bearing.mut == 2) %>%
  count(Tumor_Barcode, name = "CNV2")

tmp <- DP_info_data %>%
  count(Tumor_Barcode, name = 'Total') %>%
  left_join(tmp)
tmp <- tmp %>% replace_na(list(CNV2 = 0))

tdata <- ngspurity %>%
  select(Tumor_Barcode, WGD_Status = MCN_WGD) %>%
  left_join(tmp) %>%
  mutate(Ratio = CNV2 / Total)


# driver gene
load('../data/tgct_mtvcf.RData', verbose = T)
load('../data/tgct_maf.RData', verbose = T)
#load('tgct_driver_mutations.RData')
load('../data/intogene_drivers.RData', verbose = T)
genelist <- intogene %>% pull(SYMBOL)

tmp <- tgct_mtvcf %>% select(Tumor_Barcode, MutationID, CLS)

#tmp <- tgct_maf %>%
tmp <- tgct_maf %>%
  mutate(Chromosome = str_remove(Chromosome, "^chr")) %>%
  mutate(
    MutationID = paste(
      Chromosome,
      Start_Position,
      Reference_Allele,
      Tumor_Seq_Allele2,
      sep = ':'
    )
  ) %>%
  select(Tumor_Barcode, MutationID, Hugo_Symbol, Variant_Classification) %>%
  left_join(tmp)

tmpdata <- tmp %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  filter(Hugo_Symbol %in% genelist)


# percentage
left_join(
  tmpdata %>%
    group_by(Subtype, Hugo_Symbol) %>%
    filter(CLS == 'clonal [early]') %>%
    summarise(n1 = n_distinct(Tumor_Barcode)) %>%
    ungroup(),
  tmpdata %>%
    group_by(Subtype, Hugo_Symbol) %>%
    summarise(n2 = n_distinct(Tumor_Barcode)) %>%
    ungroup()
) %>%
  mutate(freq = n1 / n2)


pdata <- tdata %>%
  left_join(clsdata) %>%
  left_join(tmpdata %>% select(-Subtype)) %>%
  #left_join(tgct_maf %>% select(Tumor_Barcode,Hugo_Symbol) %>% unique() %>% filter(Hugo_Symbol %in% genelist)) %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  arrange(!is.na(Hugo_Symbol)) %>%
  filter(WGD_Status == 'WGD', CLS == 'clonal [early]')


#ggrepel::geom_text_repel(aes(label=Hugo_Symbol))
pdata %>%
  ggplot(aes(Ratio, Freq)) +
  geom_point(aes(fill = Hugo_Symbol), pch = 21, size = 4, stroke = 0.1) +
  facet_wrap(~Subtype) +
  scale_fill_d3() +
  stat_smooth(method = "lm") +
  stat_cor(size = 5, col = ncicolpal[1], ) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_d3() +
  #scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  #facet_wrap(~Subtype,scale='free')+
  labs(
    x = 'WGD timming (% mutations with copy number 2) ',
    y = 'Proportion of early clonal mutations',
    fill = 'Early clonal driver mutation'
  ) +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    grid = 'XY',
    ticks = T
  ) +
  panel_border(color = 'black') +
  theme(
    #axis.text.x = element_blank(),
    panel.spacing.x = unit(0.4, 'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# ggsave(
#   filename = 'WGD_timing_early_clonal.pdf',
#   width = 9,
#   height = 6,
#   device = cairo_pdf
# )


# Fig. 5f -----------------------------------------------------------------
load('../data/predup_with_divs.RData')

predup_with_divs %>%
  ggplot(aes(cell_divisions_mid, age_at_diagnosis, , fill = Subtype)) +
  geom_point(pch = 21, size = 2.5) +
  stat_smooth(aes(fill = Subtype, color = Subtype), method = "lm") +
  stat_cor(size = 5, col = ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  facet_wrap(~Subtype, scale = 'free_x') +
  labs(x = 'Number of cell division before WGD', y = 'Age at diagnosis') +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    grid = 'XY',
    ticks = T
  ) +
  panel_border(color = 'black') +
  theme(
    #axis.text.x = element_blank(),
    panel.spacing.x = unit(0.4, 'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  guides(fill = 'none', col = 'none')



# Fig. 5g -----------------------------------------------------------------

predup_with_divs %>%
  ggplot(aes(cell_divisions_mid, age_at_diagnosis, , fill = Subtype)) +
  geom_point(pch = 21, size = 2.5) +
  stat_smooth(aes(fill = Subtype, color = Subtype), method = "lm") +
  stat_cor(size = 5, col = ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  facet_wrap(~Subtype, scale = 'free_x') +
  labs(x = 'Number of cell division before WGD', y = 'Age at diagnosis') +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    grid = 'XY',
    ticks = T
  ) +
  panel_border(color = 'black') +
  theme(
    #axis.text.x = element_blank(),
    panel.spacing.x = unit(0.4, 'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  guides(fill = 'none', col = 'none')

# ggsave(
#   filename = 'Cell_Division_before_WGD_age.pdf',
#   width = 9,
#   height = 4,
#   device = cairo_pdf
# )
# 


