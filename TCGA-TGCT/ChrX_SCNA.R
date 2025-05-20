set_wd()
libztw()
pdfhr2()
library(ggnewscale)

load('wgs_data.RData', verbose = T)
load('ngspurity_data.RData', verbose = T)

scna_chrx <- read_delim(
  '/Volumes/Zhang_Group/ZTW/TCGA-TGCT/NGSpurity/Result/BBprofile_final_chrX.txt',
  delim = '\t',
  col_names = T
)

# switch the nMaj and nMinor
tmp <- scna_chrx %>%
  mutate(size = abs(endpos - startpos)) %>%
  group_by(Tumor_Barcode) %>%
  arrange(desc(size)) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(wgs_data) %>%
  arrange(Subtype, desc(nMaj1)) %>%
  select(Tumor_Barcode, Subtype) %>%
  mutate(seq = seq_along(Tumor_Barcode))

scna_chrx <- left_join(scna_chrx, tmp)

scna_color <- c(pal_gsea()(12)[c(2, 6:12)], 'darkred')
names(scna_color) <- as.character(c(0:8))

p1 <- scna_chrx %>%
  #filter(Tumor_Barcode == 'TCGA-2G-AAEQ-01A-11D-A734-36') %>%
  ggplot() +
  geom_rect(aes(
    xmin = seq - 0.5,
    xmax = seq + 0.5,
    ymin = startpos / 1e6,
    ymax = endpos / 1e6,
    fill = as.factor(nMaj1)
  )) +
  scale_fill_manual(breaks = 0:8, values = scna_color) +
  #geom_hline(yintercept = -1e4)+
  labs(x = NULL, y = 'Chromosome X (Mbp)', fill = 'SCNA') +
  guides(fill = guide_legend(nrow = 1)) +
  new_scale_fill() +
  geom_rect(aes(
    xmin = seq - 0.5,
    xmax = seq + 0.5,
    ymin = -5,
    ymax = 0,
    fill = frac1
  )) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 1)) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 7),
    limits = c(-5, 156040895 / 1e6),
    expand = c(0, 0)
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 7), expand = c(0, 0)) +
  theme_ipsum_rc(
    base_size = 13,
    axis_title_just = 'm',
    axis_title_size = 15,
    plot_margin = margin(5.5, 5.5, 5.5, 5.5),
    ticks = T,
    grid = F
  ) +
  theme(
    panel.spacing.x = unit(0.1, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 'bold'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.direction = 'horizontal',
    legend.key.width = unit(2, 'cm')
  ) +
  panel_border(color = '#cccccc', size = 0.8) +
  facet_wrap(~Subtype, scales = 'free_x') +
  labs(fill = 'Clonal Fraction\n')

p2 <- scna_chrx %>%
  mutate(nMaj2 = if_else(frac2 == 0, NA_real_, nMaj2)) %>%
  #filter(Tumor_Barcode == 'TCGA-2G-AAEQ-01A-11D-A734-36') %>%
  ggplot() +
  geom_rect(aes(
    xmin = seq - 0.5,
    xmax = seq + 0.5,
    ymin = startpos / 1e6,
    ymax = endpos / 1e6,
    fill = as.factor(nMaj2)
  )) +
  scale_fill_manual(breaks = c(0:8), values = scna_color, na.value = 'white') +
  #geom_hline(yintercept = -1e4)+
  labs(x = NULL, y = 'Chromosome X (Mbp)', fill = 'SCNA') +
  guides(fill = guide_legend(nrow = 1)) +
  new_scale_fill() +
  geom_rect(aes(
    xmin = seq - 0.5,
    xmax = seq + 0.5,
    ymin = -5,
    ymax = 0,
    fill = frac2
  )) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 1)) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 7),
    limits = c(-5, 156040895 / 1e6),
    expand = c(0, 0)
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 7), expand = c(0, 0)) +
  theme_ipsum_rc(
    base_size = 13,
    axis_title_just = 'm',
    axis_title_size = 15,
    plot_margin = margin(5.5, 5.5, 5.5, 5.5),
    ticks = T,
    grid = F
  ) +
  theme(
    panel.spacing.x = unit(0.1, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 'bold'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.direction = 'horizontal',
    legend.key.width = unit(2, 'cm')
  ) +
  panel_border(color = '#cccccc', size = 0.8) +
  facet_wrap(~Subtype, scales = 'free_x') +
  labs(fill = 'Subclonal Fraction\n')


# WGD status
p3 <- scna_chrx %>%
  left_join(
    ngspurity %>% select(Tumor_Barcode, MCN_WGD)
  ) %>%
  ggplot() +
  geom_rect(aes(
    xmin = seq - 0.5,
    xmax = seq + 0.5,
    ymin = -1,
    ymax = 1,
    fill = as.factor(MCN_WGD)
  )) +
  scale_fill_manual(values = ncicolpal[c(14, 18)]) +
  #geom_hline(yintercept = -1e4)+
  labs(x = NULL, y = 'Chromosome X (Mbp)', fill = 'SCNA') +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 7),
    limits = c(-20, 20),
    expand = c(0, 0)
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 7), expand = c(0, 0)) +
  theme_ipsum_rc(
    base_size = 13,
    axis_title_just = 'm',
    axis_title_size = 15,
    plot_margin = margin(5.5, 5.5, 5.5, 5.5),
    ticks = T,
    grid = F
  ) +
  theme(
    panel.spacing.x = unit(0.1, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 'bold'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.direction = 'horizontal'
  ) +
  #panel_border(color = '#cccccc',size=0.8)+
  facet_wrap(~Subtype, scales = 'free_x') +
  labs(fill = 'WGD status\n')

plot_grid(p1, p2, p3, align = 'v', axis = 'lr', nrow = 3)

ggsave(
  filename = 'chrx_SCNA_clonality.pdf',
  width = 15,
  height = 8,
  device = cairo_pdf
)


save(scna_chrx, file = 'scna_chrx.RData')


scna_chrx %>%
  mutate(size = abs(endpos - startpos)) %>%
  filter(nMaj1 > 2) %>%
  group_by(Tumor_Barcode) %>%
  summarise(Size = sum(size)) %>%
  filter(Size > 0.5 * 156040895) %>%
  nrow() /
  252


scna_chrx %>%
  filter(Subtype == 'Seminoma') %>%
  mutate(size = abs(endpos - startpos)) %>%
  filter(nMaj1 > 2) %>%
  group_by(Tumor_Barcode) %>%
  summarise(Size = sum(size)) %>%
  filter(Size > 0.5 * 156040895) %>%
  nrow() /
  143

scna_chrx %>%
  filter(Subtype == 'Non-Seminoma') %>%
  mutate(size = abs(endpos - startpos)) %>%
  filter(nMaj1 > 2) %>%
  group_by(Tumor_Barcode) %>%
  summarise(Size = sum(size)) %>%
  filter(Size > 0.5 * 156040895) %>%
  nrow() /
  109


# subclone
scna_chrx %>%
  mutate(size = abs(endpos - startpos)) %>%
  filter(frac2 > 0) %>%
  group_by(Tumor_Barcode) %>%
  summarise(Size = sum(size)) %>%
  filter(Size > 0.5 * 156040895) %>%
  nrow() /
  252

scna_chrx %>%
  filter(Subtype == 'Seminoma') %>%
  mutate(size = abs(endpos - startpos)) %>%
  filter(frac2 > 0) %>%
  group_by(Tumor_Barcode) %>%
  summarise(Size = sum(size)) %>%
  filter(Size > 0.5 * 156040895) %>%
  nrow() /
  143

scna_chrx %>%
  filter(Subtype == 'Non-Seminoma') %>%
  mutate(size = abs(endpos - startpos)) %>%
  filter(frac2 > 0) %>%
  group_by(Tumor_Barcode) %>%
  summarise(Size = sum(size)) %>%
  filter(Size > 0.5 * 156040895) %>%
  nrow() /
  109


# SCNA_X and replication stress -------------------------------------------
load('tcga_repstress_zscore.RData', verbose = T)
load('tcga_xist_exp.RData', verbose = T)

tmp1 <- scna_chrx %>%
  mutate(size = abs(endpos - startpos)) %>%
  filter(nMaj1 > 2) %>%
  group_by(Tumor_Barcode) %>%
  summarise(Size = sum(size)) %>%
  filter(Size > 0.5 * 156040895) %>%
  pull(Tumor_Barcode)

tmp2 <- scna_chrx %>%
  mutate(size = abs(endpos - startpos)) %>%
  filter(nMaj1 > 1) %>%
  group_by(Tumor_Barcode) %>%
  summarise(Size = sum(size)) %>%
  filter(Size > 0.5 * 156040895) %>%
  filter(!(Tumor_Barcode %in% tmp1)) %>%
  pull(Tumor_Barcode)


tmp <- wgs_data %>%
  select(Subject, Tumor_Barcode, Subtype) %>%
  mutate(
    SCNA = case_when(
      Tumor_Barcode %in% tmp1 ~ 'Gain >=2 copies',
      Tumor_Barcode %in% tmp2 ~ 'Gain 1 copy',
      TRUE ~ 'No gain'
    )
  )

# XIST expression
my_comparisons <- list(
  c('Gain >=2 copies', 'Gain 1 copy'),
  c('Gain 1 copy', 'No gain'),
  c('Gain >=2 copies', 'No gain')
)

tdata <- xist_exp_tgct %>%
  left_join(tmp) %>%
  filter(!is.na(SCNA)) %>%
  mutate(RSEM = log2(RSEM + 1)) %>%
  mutate(SCNA = fct_rev(SCNA))

stat.test <- tdata %>%
  group_by(Subtype) %>%
  rstatix::wilcox_test(RSEM ~ SCNA, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "SCNA") %>%
  mutate(myformatted.p = sprintf("P = %.2f", p)) %>%
  mutate(SCNA = NA)

stat.test$y.position <- c(15, 16, 17, 15.2, 16, 17)

tdata %>%
  ggplot(aes(SCNA, RSEM, fill = SCNA)) +
  facet_wrap(~Subtype) +
  ggbeeswarm::geom_quasirandom(
    pch = 21,
    size = 2,
    width = 0.2,
    color = "black",
    stroke = 0.4
  ) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5, size = 0.4) +
  scale_fill_manual(values = as.character(scna_color[c(2, 5, 9)])) +
  scale_y_continuous(breaks = pretty_breaks(n = 7)) +
  labs(
    y = 'XIST mRNA expression\nlog2(RNA Seq V2 RSEM + 1)',
    x = '',
    color = 'Group'
  ) +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    ticks = T,
    grid = 'Y'
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0.5),
    plot.margin = margin(4, 4, -6, 4)
  ) +
  guides(fill = "none") +
  panel_border(color = 'black', linetype = 1) +
  panel_border(color = 'black', size = 0.5) +
  theme(legend.position = 'none') +
  stat_pvalue_manual(stat.test, label = "myformatted.p", color = "black")

ggsave(
  filename = 'chrx_SCNA_XIST_exp.pdf',
  width = 4.5,
  height = 6,
  device = cairo_pdf
)


tdata %>%
  mutate(Group = if_else(SCNA == 'No gain', SCNA, 'Gain')) %>%
  group_by(Subtype) %>%
  do(tidy(wilcox.test(RSEM ~ Group, data = .)))
# P=0.111 and 0.0265

tdata %>%
  mutate(X = as.integer(SCNA)) %>%
  group_by(Subtype) %>%
  do(tidy(lm(RSEM ~ X, data = .)))

tdata %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  mutate(X = as.integer(SCNA)) %>%
  group_by(Subtype) %>%
  do(tidy(lm(RSEM ~ X + BB_Purity, data = .)))


## RepStress Zscore

my_comparisons <- list(
  c('Gain >=2 copies', 'Gain 1 copy'),
  c('Gain 1 copy', 'No gain'),
  c('Gain >=2 copies', 'No gain')
)

tdata <- tcga_repstress_score_cbioportal %>%
  mutate(Subject = str_sub(TCGA_Barcode, 1, 12)) %>%
  left_join(tmp) %>%
  filter(!is.na(SCNA)) %>%
  mutate(SCNA = fct_rev(SCNA))

stat.test <- tdata %>%
  group_by(Subtype) %>%
  rstatix::wilcox_test(
    RepStress_ZScore ~ SCNA,
    comparisons = my_comparisons
  ) %>%
  rstatix::add_xy_position(x = "SCNA") %>%
  mutate(myformatted.p = sprintf("P = %.2f", p)) %>%
  mutate(SCNA = NA)

#stat.test$y.position <- c(4.5,4.6,4.8)

tdata %>%
  ggplot(aes(SCNA, RepStress_ZScore, fill = SCNA)) +
  facet_wrap(~Subtype) +
  ggbeeswarm::geom_quasirandom(
    pch = 21,
    size = 2,
    width = 0.2,
    color = "black",
    stroke = 0.4
  ) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5, size = 0.4) +
  scale_fill_manual(values = as.character(scna_color[c(2, 5, 9)])) +
  scale_y_continuous(breaks = pretty_breaks(n = 7)) +
  labs(y = 'Repstress Z-score', x = '', color = 'Group') +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    ticks = T,
    grid = 'Y'
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0.5),
    plot.margin = margin(4, 4, -6, 4)
  ) +
  guides(fill = "none") +
  panel_border(color = 'black', linetype = 1) +
  panel_border(color = 'black', size = 0.5) +
  theme(legend.position = 'none')
#stat_pvalue_manual(stat.test, label = "myformatted.p", color = "black")

ggsave(
  filename = 'chrx_SCNA_Repstress.pdf',
  width = 4.5,
  height = 6,
  device = cairo_pdf
)


tdata %>%
  mutate(Group = if_else(SCNA == 'No gain', SCNA, 'Gain')) %>%
  group_by(Subtype) %>%
  do(tidy(wilcox.test(RepStress_ZScore ~ Group, data = .)))
# P=0.0485, 0.82

tdata %>%
  mutate(Group = if_else(SCNA == 'No gain', SCNA, 'Gain')) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  mutate(X = as.integer(SCNA)) %>%
  group_by(Subtype) %>%
  do(tidy(lm(RepStress_ZScore ~ X + BB_Purity, data = .)))

# RepStree score vs XIST expression ---------------------------------------

tdata <- left_join(
  tcga_repstress_score_cbioportal %>% filter(Study == 'TGCT'),
  xist_exp_tgct %>%
    select(TCGA_Barcode, Subtype, RSEM) %>%
    mutate(RSEM = log2(RSEM + 1))
) %>%
  filter(!is.na(RSEM))


tdata %>%
  ggplot(aes(RSEM, RepStress_ZScore, fill = Subtype)) +
  geom_point(pch = 21, size = 3, color = "black", stroke = 0.4) +
  geom_smooth(method = 'lm') +
  facet_wrap(~Subtype, scales = 'free') +
  scale_fill_manual(values = subtypecol) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = pretty_breaks(n = 7)) +
  labs(
    x = 'XIST mRNA expression log2(RNA Seq V2 RSEM + 1)',
    x = '',
    color = 'Group'
  ) +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    ticks = T,
    grid = 'Y'
  ) +
  theme(
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0.5),
    plot.margin = margin(4, 4, 0, 4)
  ) +
  guides(fill = "none") +
  panel_border(color = 'black', linetype = 1) +
  panel_border(color = 'black', size = 0.5) +
  theme(legend.position = 'none') +
  stat_cor(label.y = c(2.3, 4.2))
#stat_regline_equation(label.y = c(1.8,3))

ggsave(
  filename = 'XIST_repStressZscore.pdf',
  width = 8,
  height = 4,
  device = cairo_pdf
)
