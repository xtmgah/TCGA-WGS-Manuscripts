libztw()
set_wd()
pdfhr()
pdfhr2()

load('wgs_data.RData')
load('ngspurity_data.RData')
load('wgs_clinical.RData')

## load telseq

telseq <- read_delim(
  "/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Telomere/telseq/telseq_u.out",
  delim = "\t",
  col_names = T
) %>%
  select(-ReadGroup, -Library, -`...44`) %>%
  unique()

telseq <- telseq %>% select(Barcode = Sample, Telseq_TL = LENGTH_ESTIMATE)

telseq <- wgs_data %>%
  select(Subtype, Tumor_Barcode, Normal_Barcode) %>%
  left_join(
    telseq %>% select(Tumor_Barcode = Barcode, Telseq_TL_Tumor = Telseq_TL)
  ) %>%
  left_join(
    telseq %>% select(Normal_Barcode = Barcode, Telseq_TL_Normal = Telseq_TL)
  ) %>%
  mutate(Telseq_TL_Ratio = log2(Telseq_TL_Tumor / Telseq_TL_Normal))


#tmp <- telseq %>% filter_all(any_vars(is.na(.)|is.infinite(.)))

pdata <- telseq %>%
  rename(
    `Normal TL` = Telseq_TL_Normal,
    `Tumor TL` = Telseq_TL_Tumor,
    `Ratio (Tumor_TL/Normal_TL)` = Telseq_TL_Ratio
  ) %>%
  pivot_longer(cols = c(-Subtype, -Tumor_Barcode, -Normal_Barcode)) %>%
  mutate(
    name = factor(
      name,
      levels = c("Normal TL", "Tumor TL", "Ratio (Tumor_TL/Normal_TL)")
    )
  )

my_comparisons <- list(c('Seminoma', 'Non-Seminoma'))

stat.test <- pdata %>%
  group_by(name) %>%
  rstatix::wilcox_test(value ~ Subtype, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Subtype") %>%
  mutate(
    myformatted.p = if_else(
      p < 0.001,
      sprintf("P = %.2e", p),
      sprintf("P = %.2f", p)
    )
  ) %>%
  mutate(Subtype = my_comparisons[[1]][1]) %>%
  mutate(col = if_else(p < 0.05, ncicolpal[1], 'black'))

stat.test


tdata <- telseq %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subject)) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  left_join(wgs_clinical %>% select(Subject, age_at_diagnosis))

a <- tdata %>%
  do(tidy(lm(Telseq_TL_Normal ~ Subtype + age_at_diagnosis, data = .))) %>%
  filter(term == 'SubtypeNon-Seminoma') %>%
  mutate(
    myformatted.p = if_else(
      p.value < 0.001,
      sprintf("P = %.2e", p.value),
      sprintf("P = %.2f", p.value)
    )
  ) %>%
  pull(myformatted.p)

b <- tdata %>%
  do(tidy(lm(
    Telseq_TL_Tumor ~ Subtype + age_at_diagnosis + BB_Purity,
    data = .
  ))) %>%
  filter(term == 'SubtypeNon-Seminoma') %>%
  mutate(
    myformatted.p = if_else(
      p.value < 0.001,
      sprintf("P = %.2e", p.value),
      sprintf("P = %.2f", p.value)
    )
  ) %>%
  pull(myformatted.p)
c <- tdata %>%
  do(tidy(lm(
    Telseq_TL_Ratio ~ Subtype + age_at_diagnosis + BB_Purity,
    data = .
  ))) %>%
  filter(term == 'SubtypeNon-Seminoma') %>%
  mutate(
    myformatted.p = if_else(
      p.value < 0.001,
      sprintf("P = %.2e", p.value),
      sprintf("P = %.2f", p.value)
    )
  ) %>%
  pull(myformatted.p)

stat.test$myformatted.p <- c(a, b, c)


pdata %>%
  ggplot(aes(Subtype, value, fill = Subtype)) +
  #geom_violin(trim=T,size=0.2)+
  geom_boxplot(
    width = 0.8,
    fill = "gray95",
    color = 'black',
    outlier.shape = 21,
    outlier.size = 0.8
  ) +
  ggbeeswarm::geom_quasirandom(
    pch = 21,
    size = 2,
    width = 0.3,
    color = "white",
    stroke = 0.2
  ) +
  facet_wrap(~name, scales = 'free_y') +
  labs(x = NULL, y = 'Telomere value') +
  scale_fill_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 7)) +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    grid = 'XY',
    ticks = T
  ) +
  panel_border(color = 'black') +
  theme(
    axis.text.x = element_blank(),
    panel.spacing.x = unit(0.4, 'cm'),
    axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 5)
  ) +
  stat_pvalue_manual(data = stat.test, col = 'col', label = "myformatted.p") +
  scale_color_manual(values = c(ncicolpal[1], 'black')) +
  guides(color = 'none')

ggsave(
  filename = 'Telseq_subtype_comparisons.pdf',
  width = 4.5,
  height = 6,
  device = cairo_pdf
)


tmp <- telseq %>% count(Subtype, name = 'Total')

telseq %>%
  group_by(Subtype) %>%
  filter(Telseq_TL_Tumor > Telseq_TL_Normal) %>%
  count(Subtype) %>%
  left_join(tmp) %>%
  mutate(Freq = n / Total)

telseq %>%
  group_by(Subtype) %>%
  filter(Telseq_TL_Tumor < Telseq_TL_Normal) %>%
  count(Subtype) %>%
  left_join(tmp) %>%
  mutate(Freq = n / Total)


save(telseq, file = 'tgct_telomere_length.RData')


# Association with Signature ----------------------------------------------
load('ngspurity_data.RData')
load('wgs_clinical.RData')
load('tgct_mutational_signatures.RData', verbose = T)

tmp <- telseq %>%
  rename(
    `Normal TL` = Telseq_TL_Normal,
    `Tumor TL` = Telseq_TL_Tumor,
    `log2(Tumor_TL/Normal_TL)` = Telseq_TL_Ratio
  ) %>%
  pivot_longer(
    cols = c(-Subtype, -Tumor_Barcode, -Normal_Barcode),
    names_to = 'TL',
    values_to = 'TL_value'
  ) %>%
  mutate(
    TL = factor(
      TL,
      levels = c("Normal TL", "Tumor TL", "log2(Tumor_TL/Normal_TL)")
    )
  ) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subject)) %>%
  left_join(wgs_clinical %>% select(Subject, Age = age_at_diagnosis))

tmp %>%
  group_by(TL) %>%
  do(tresult = safely(lm)(TL_value ~ BB_Purity + Subtype + Age, data = .)) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(TL, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup()


tdata <- tgct_activity_all %>%
  pivot_longer(-Tumor_Barcode) %>%
  left_join(tmp)

tresult_all <- tdata %>%
  mutate(value = log2(value + 1)) %>%
  group_by(TL, name) %>%
  do(
    tresult = safely(lm)(TL_value ~ value + BB_Purity + Subtype + Age, data = .)
  ) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(TL, name, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(term == 'value') %>%
  filter(TL == 'log2(Tumor_TL/Normal_TL)') %>%
  arrange(p.value) %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

tresult_all

tresult_all %>%
  ggplot(aes((estimate), -log10(FDR), fill = name)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'red', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, col = '#cccccc', size = 0.5) +
  geom_point(pch = 21, stroke = 0.2, size = 4) +
  scale_fill_manual(values = sigcol) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  ggrepel::geom_text_repel(
    data = tresult_all %>% filter(FDR < 0.05),
    aes(label = name),
    max.overlaps = 30,
    size = 4
  ) +
  labs(
    x = 'Linear regression coefficient',
    y = '-log10(FDR)',
    fill = 'Signature'
  ) +
  guides(fill = 'none') +
  theme_ipsum_rc(
    base_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    grid = F,
    ticks = T
  ) +
  theme(plot.margin = margin(4, 4, 4, 4)) +
  panel_border(color = 'black', size = 0.5) +
  coord_cartesian(clip = 'off')

ggsave(
  filename = 'TL_fdr_signatures.pdf',
  width = 5,
  height = 4,
  device = cairo_pdf
)


tdata %>%
  filter(name == 'SBS87') %>%
  ggplot(aes(log2(value + 1), TL_value, fill = Subtype)) +
  geom_point(pch = 21, size = 2.5) +
  stat_smooth(aes(fill = Subtype, color = Subtype), method = "lm") +
  stat_cor(method = 'spearman', size = 5, col = ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  ggh4x::facet_grid2(Subtype ~ TL, scale = 'free', independent = 'y') +
  labs(x = 'log2(SBS87 + 1)', y = 'Telomere Length (TL) value') +
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
    panel.spacing = unit(0.1, 'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  guides(fill = 'none', col = 'none')

ggsave(
  filename = 'SBS87_subtype_tl.pdf',
  width = 10,
  height = 6,
  device = cairo_pdf
)

tdata %>%
  filter(name == 'SBS1') %>%
  ggplot(aes(log2(value + 1), TL_value, fill = Subtype)) +
  geom_point(pch = 21, size = 2.5) +
  stat_smooth(aes(fill = Subtype, color = Subtype), method = "lm") +
  stat_cor(size = 5, col = ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  ggh4x::facet_grid2(TL ~ Subtype, scale = 'free', independent = 'y') +
  labs(x = 'log2(SBS1 + 1)', y = 'Telomere Length (TL) value') +
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
    panel.spacing = unit(0.4, 'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  guides(fill = 'none', col = 'none')

ggsave(
  filename = 'SBS1_subtype_tl.pdf',
  width = 8,
  height = 9,
  device = cairo_pdf
)


tmp %>%
  ggplot(aes(Age, TL_value, fill = Subtype)) +
  geom_point(pch = 21, size = 2.5) +
  stat_smooth(aes(fill = Subtype, color = Subtype), method = "lm") +
  stat_cor(size = 5, col = ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  ggh4x::facet_grid2(TL ~ Subtype, scale = 'free', independent = 'y') +
  labs(x = 'Age at diagnosis', y = 'Telomere Length (TL) value') +
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
    panel.spacing = unit(0.4, 'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(hjust = 0.5, face = 'bold', size = 13)
  ) +
  guides(fill = 'none', col = 'none')

ggsave(
  filename = 'Age_subtype_tl.pdf',
  width = 8,
  height = 8,
  device = cairo_pdf
)

# Association with genomic features ---------------------------------------
tdata <- telseq %>%
  filter(!is.na(Telseq_TL_Tumor)) %>%
  select(Tumor_Barcode, value = Telseq_TL_Tumor) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subject)) %>%
  left_join(wgs_clinical %>% select(Subject, Age = age_at_diagnosis))

load('data_top.RData')

tgct_data_full <- wgs_data %>%
  select(Subject, Tumor_Barcode) %>%
  left_join(data_top) %>%
  unique() %>%
  mutate(Alteration = 'Yes') %>%
  complete(
    Gene,
    nesting(Subject, Tumor_Barcode, Type),
    fill = list(Alteration = 'No')
  ) %>%
  filter(!is.na(Gene), !is.na(Type)) %>%
  select(Subject, Tumor_Barcode, Gene, Alteration, Type)

## wicox test ##
tdata <- tgct_data_full %>%
  filter(Tumor_Barcode %in% tdata$Tumor_Barcode) %>%
  mutate(Gene = paste0(Gene, "|", Type)) %>%
  select(Tumor_Barcode, Gene, Alteration) %>%
  left_join(tdata)

tdata_all <- tdata %>% left_join(wgs_data %>% select(Tumor_Barcode, Subtype))

#tdata <- tdata %>% filter(Freq > 0.03)

tresult_all <- tdata_all %>%
  group_by(Gene) %>%
  do(
    tresult = safely(lm)(
      value ~ Subtype + Alteration + BB_Purity + Age,
      data = .
    )
  ) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']], exponentiate = TRUE))) %>%
  select(Gene, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(term == 'AlterationYes') %>%
  arrange(p.value) %>%
  mutate(TMP = Gene) %>%
  separate(col = 'TMP', into = c('Gene_short', 'Type'), sep = '\\|') %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))


# visualizaiton
tresult_all %>%
  ggplot(aes(log2(estimate), -log10(FDR), fill = Type)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'red', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, col = '#cccccc', size = 0.5) +
  geom_point(pch = 21, stroke = 0.2, size = 3) +
  scale_fill_nejm() +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  labs(
    x = 'Odds ratio (log2)',
    y = '-log10(FDR)',
    fill = 'Genomic alterations'
  ) +
  #guides(fill = 'none') +
  theme_ipsum_rc(
    base_size = 12,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = F,
    ticks = T
  ) +
  theme(plot.margin = margin(4, 4, 4, 4)) +
  panel_border(color = 'black', size = 0.5) +
  coord_cartesian(clip = 'off')

ggsave(
  filename = 'tl_tumor_fdr_alterations.pdf',
  width = 5,
  height = 4,
  device = cairo_pdf
)
