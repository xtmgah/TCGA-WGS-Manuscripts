set_wd()
libztw()
pdfhr()
pdfhr2()

load('wgs_data.RData')
load('~/NIH-Work/R/annotable_ens90.RData')


grch38 %>% filter(chr == 'Y') %>% filter(start > 21500000, end < 21900000)


gdata <- bind_rows(
  grch38 %>%
    filter(chr == 'Y') %>%
    filter(start > 21500000, end < 21900000) %>%
    filter(biotype == 'protein_coding'),

  grch38 %>%
    filter(chr == 'Y') %>%
    filter(start > 9200000, end < 11800000) %>%
    filter(biotype == 'protein_coding')
) %>%
  mutate(Subtype = 'Non-Seminoma') %>%
  select(gene = symbol, start, end, Subtype) %>%
  mutate(Subtype = factor(Subtype, levels = c('Seminoma', 'Non-Seminoma')))

scna <- read_delim(
  '/Volumes/Zhang_Group/ZTW/TCGA-TGCT/cnvkit/all_call_purity.cns.gz',
  delim = '\t',
  col_names = T
)


scnay <- scna %>%
  filter(chromosome == 'chrY') %>%
  mutate(cn = if_else(cn < 0, 0, cn)) %>%
  mutate(end = if_else(end > 23881271, 23881271, end))

# clustering scnay and get the order
# Assuming your original dataframe is called 'scnay'
scnay_subset <- scnay %>%
  select(Tumor_Barcode, Subtype, start, end, cn)

# Discretize genomic coordinates into bins to create a copy number matrix
bin_size <- 1e5 # Adjust bin size if needed

# Create bins
genomic_bins <- seq(
  min(scnay_subset$start),
  max(scnay_subset$end),
  by = bin_size
)

# Function to expand segments into bins
expand_to_bins <- function(df, bins) {
  expanded <- map_dfr(unique(df$Tumor_Barcode), function(tb) {
    tb_df <- df %>% filter(Tumor_Barcode == tb)
    bins_df <- tibble(bin = bins)
    bins_df <- bins_df %>%
      mutate(
        cn = map_dbl(bin, function(b) {
          segment <- tb_df %>% filter(start <= b, end >= b)
          if (nrow(segment) == 0) return(NA)
          segment$cn[1]
        })
      ) %>%
      mutate(Tumor_Barcode = tb)
    bins_df
  })
  expanded
}

expanded_cnv <- expand_to_bins(scnay_subset, genomic_bins)

# Reshape to wide matrix form
cn_matrix <- expanded_cnv %>%
  pivot_wider(names_from = bin, values_from = cn) %>%
  column_to_rownames(var = "Tumor_Barcode")

# Replace NA with median copy number (or alternative choice)
#cn_matrix[is.na(cn_matrix)] <- median(cn_matrix, na.rm = TRUE)
cn_matrix[is.na(cn_matrix)] <- 1

# Clustering
clustered <- pheatmap::pheatmap(
  cn_matrix,
  scale = "none",
  clustering_method = "ward.D2"
)

# Order of tumor barcodes after clustering
ordered_barcodes <- rownames(cn_matrix)[clustered$tree_row$order]

# Output ordered barcodes
scnay_order <- tibble(Tumor_Barcode = ordered_barcodes) %>%
  left_join(wgs_data) %>%
  mutate(seq1 = seq_along(Tumor_Barcode)) %>%
  arrange(Subtype, rev(seq1)) %>%
  select(Tumor_Barcode, Subtype) %>%
  mutate(seq = seq_along(Tumor_Barcode)) %>%
  select(Tumor_Barcode, Subtype, seq)


# tmp <- scnay %>%
#   mutate(size = abs(end - start)) %>%
#   group_by(Tumor_Barcode) %>%
#   arrange(desc(size)) %>%
#   slice(1) %>%
#   ungroup() %>%
#   left_join(wgs_data) %>%
#   arrange(Subtype, desc(cn)) %>%
#   select(Tumor_Barcode, Subtype) %>%
#   mutate(seq = seq_along(Tumor_Barcode))

scnay <- left_join(scnay, scnay_order)


scna_color <- c(pal_gsea()(12)[c(2, 6:12)], 'darkred')
names(scna_color) <- as.character(c(0:8))

p1 <- scnay %>%
  ggplot() +
  geom_rect(
    aes(
      xmin = seq - 0.5,
      xmax = seq + 0.5,
      ymin = start / 1e6,
      ymax = end / 1e6,
      fill = as.factor(cn)
    ),
    col = 'white',
    linewidth = 0.001
  ) +
  scale_fill_manual(breaks = 0:8, values = scna_color) +
  labs(x = NULL, y = 'Chromosome Y (Mbp)', fill = 'SCNA') +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 7),
    limits = c(NA, 23881271 / 1e6),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 7),
    expand = c(0, 0)
  ) +
  theme_ipsum_rc(
    base_size = 13,
    axis_title_just = 'm',
    axis_title_size = 15,
    plot_margin = margin(5.5, 5.5, 5.5, 5.5),
    ticks = T,
    grid = F
  ) +
  theme(
    panel.spacing.x = unit(0.15, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 'bold'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.direction = 'horizontal',
    legend.key.width = unit(2, 'cm')
  ) +
  panel_border(color = '#cccccc', size = 0.8) +
  facet_wrap(~Subtype, scales = 'free_x') +
  labs(fill = 'Chromosome Y SCNA')
#  geom_hline(yintercept = c(21.5, 21.9), colour = 'black') +
#  geom_hline(yintercept = c(9.2, 11.8), colour = 'black')
# ggrepel::geom_text_repel(data=gdata,aes(label=gene,x=252,y=start/1e6),
#                          force        = 0.5,
#                          nudge_x      = -0.15,
#                          direction    = "y",
#                          hjust        = 0,
#                          segment.size = 0.2)

ggsave(
  filename = 'chrY_SCNA.pdf',
  plot = p1,
  width = 15,
  height = 4,
  device = cairo_pdf
)


# Quantification ----------------------------------------------------------
tmp <- scnay %>%
  mutate(size = abs(end - start)) %>%
  select(Tumor_Barcode, Subtype, size, cn) %>%
  filter(size > 10e6, cn > 1) %>%
  pull(Tumor_Barcode)

wgs_data %>%
  mutate(Group = if_else(Tumor_Barcode %in% tmp, 'A', 'B')) %>%
  count(Subtype, Group) %>%
  group_by(Subtype) %>%
  mutate(freq = n / sum(n))


wgs_data %>%
  mutate(Group = if_else(Tumor_Barcode %in% tmp, 'A', 'B')) %>%
  select(Subtype, Group) %>%
  table() %>%
  fisher.test() %>%
  tidy()


# RBMY genes --------------------------------------------------------------
tmp <- scnay %>%
  filter(start > 21500000, end < 21900000) %>%
  mutate(size = abs(end - start)) %>%
  select(Tumor_Barcode, Subtype, size, cn) %>%
  filter(cn == 0) %>%
  pull(Tumor_Barcode)


wgs_data %>%
  mutate(Group = if_else(Tumor_Barcode %in% tmp, 'A', 'B')) %>%
  count(Subtype, Group) %>%
  group_by(Subtype) %>%
  mutate(freq = n / sum(n))


wgs_data %>%
  mutate(Group = if_else(Tumor_Barcode %in% tmp, 'A', 'B')) %>%
  select(Subtype, Group) %>%
  table() %>%
  fisher.test() %>%
  tidy()


# Asscoiation with other features -----------------------------------------
load('data_top.RData')

p2 <- data_top %>%
  left_join(scnay_order) %>%
  mutate(Key = paste0(Gene, "@", Alteration)) %>%
  filter(!is.na(Subtype)) %>%
  mutate(Subtype = factor(Subtype, levels = c('Seminoma', 'Non-Seminoma'))) %>%
  filter(Type != 'Mutation_Driver') %>%
  ggplot(aes(seq, Key)) +
  geom_tile(fill = 'black') +
  theme(
    panel.spacing.x = unit(0.15, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 'bold'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.direction = 'horizontal',
    legend.key.width = unit(2, 'cm')
  ) +
  panel_border(color = '#cccccc', size = 0.8) +
  facet_wrap(~Subtype, scales = 'free_x')


plot_grid(p1, p2, align = 'v', axis = 'lr', ncol = 1)


# Assocaition test --------------------------------------------------------

tmp1 <- scnay %>%
  filter(start > 21500000, end < 21900000) %>%
  mutate(size = abs(end - start)) %>%
  select(Tumor_Barcode, Subtype, size, cn) %>%
  filter(cn == 0) %>%
  #filter(size>300000) %>%
  pull(Tumor_Barcode)


tmp2 <- data_top %>%
  filter(Gene %in% c('KIT', 'KRAS', 'NRAS'), Type == 'Mutation_Driver') %>%
  select(Tumor_Barcode, Gene, Alteration) %>%
  mutate(Alteration = 'Y') %>%
  pivot_wider(names_from = Gene, values_from = Alteration, values_fill = 'N')


tdata <- wgs_data %>%
  select(Subtype, Tumor_Barcode) %>%
  left_join(tmp2) %>%
  mutate_all(~ replace_na(., "N")) %>%
  mutate(RBMY = if_else(Tumor_Barcode %in% tmp1, 'Y', 'N')) %>%
  arrange(desc(KIT), desc(KRAS), desc(NRAS), desc(RBMY)) %>%
  mutate(Tumor_Barcode = fct_inorder(Tumor_Barcode)) %>%
  filter(Subtype == 'Seminoma')

tdata %>% select(RBMY, KIT) %>% table() %>% fisher.test() %>% tidy()
tdata %>% select(RBMY, KRAS) %>% table() %>% fisher.test() %>% tidy()
tdata %>% select(RBMY, NRAS) %>% table() %>% fisher.test() %>% tidy()


tdata %>%
  pivot_longer(cols = -c(Subtype, Tumor_Barcode)) %>%
  ggplot(aes(name, Tumor_Barcode, fill = value)) +
  geom_tile()

#
# scnay %>%
#   ggplot() +
#   geom_segment(
#     aes(
#       x = start,
#       y = Tumor_Barcode,
#       xend = end,
#       yend = Tumor_Barcode,
#       col = cn
#     ),
#     size = 3
#   ) +
#   scale_color_viridis_c() +
#   coord_flip()
#
#
# barcode <- 'TCGA-S6-A8JX-01A-11D-A748-36'
#
# cnrdata <- read_delim(
#   '/Volumes/Zhang_Group/ZTW/TCGA-TGCT/cnvkit/Result/TCGA-S6-A8JX-01A-11D-A748-36/12a2ad63-fb11-4c52-9115-2819d5de81c4_wgs_gdc_realn.targetcoverage.cnn',
#   delim = '\t',
#   col_names = T
# )
#
# cnrdata <- read_delim(
#   '/Volumes/Zhang_Group/ZTW/TCGA-TGCT/cnvkit/Result/TCGA-S6-A8JX-01A-11D-A748-36/128ee293-58ca-4435-a425-f240097e43e6_wgs_gdc_realn.targetcoverage.cnn',
#   delim = '\t',
#   col_names = T
# )
#
# cnrdata <- read_delim(
#   '/Volumes/Zhang_Group/ZTW/TCGA-TGCT/cnvkit/Result/TCGA-2G-AAH3-01A-11D-A747-36/9e12d4dd-0381-49be-8146-709e7d312985_wgs_gdc_realn.targetcoverage.cnn',
#   delim = '\t',
#   col_names = T
# )
#
# cnrdata <- read_delim(
#   '/Volumes/Zhang_Group/ZTW/TCGA-TGCT/cnvkit/Result/TCGA-2G-AAH3-01A-11D-A747-36/d7905b7b-56ee-4642-b0ac-bbc195e2ab5d_wgs_gdc_realn.targetcoverage.cnn',
#   delim = '\t',
#   col_names = T
# )
#
#
# cnrdata %>%
#   filter(chromosome == 'chrY') %>%
#   mutate(pos = (start + end) / 2) %>%
#   ggplot(aes(pos, log2)) +
#   geom_point()
