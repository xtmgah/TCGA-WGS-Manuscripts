set_wd()
libztw()
pdfhr2()
myggstyle()


# presetting
library(data.table)
conflicts_prefer(cowplot::get_legend)
source(
  '/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/RDS/Oncoplots_functions.R'
)
tmp <- read_csv(
  '/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/RDS/oncoplot_colors.csv'
)
landscape_colors <- tmp$Color
names(landscape_colors) <- tmp$Name
#landscape_colors <- c(landscape_colors,sp_group_color_new)

landscape_colors <- c(
  landscape_colors,
  c('Multiple WGD' = "#14315C", 'EUR_Pop' = "#fc8d62", 'Mix_Pop' = "#542788")
)

# load Genomic alterations data -------------------------------------------------
load('wgs_data.RData', verbose = T)
# driver genes
load('tgct_maf.RData', verbose = T)
load('intogene_drivers.RData', verbose = T)
load('DriveGene.RData')


genelist <- tgct_maf_freq %>%
  filter(
    Hugo_Symbol %in% c(drivegene$Gene, intogene$SYMBOL, 'SRCAP'),
    nMutSample > 2
  ) %>%
  pull(Hugo_Symbol)


data_mutation <- tgct_maf %>%
  filter(Hugo_Symbol %in% genelist) %>%
  left_join(wgs_data %>% select(Subject, Tumor_Barcode)) %>%
  mutate(Alteration = Variant_Classification) %>%
  #mutate(Alteration=paste0(Variant_Classification,":",tx,":",aaChange)) %>%
  select(Subject, Tumor_Barcode, Gene = Hugo_Symbol, Alteration) %>%
  mutate(Type = 'Mutation_Driver') %>%
  unique()


# focal amplification
load('tgct_gistic_gene.RData')
## need to define the gene set ##
genelist2 <- c(
  genelist,
  c(
    'GRID2',
    'JAK2',
    'CHEK1',
    'FANCA',
    'SPRY3',
    'KIT',
    'HIP1',
    'CDKN1B',
    'MDM2',
    'MAPK1'
  )
) %>%
  unique()

data_gistic <- tgct_gistic_gene %>%
  mutate(GeneSymbol = str_remove(GeneSymbol, "\\|.*")) %>%
  filter(GeneSymbol %in% genelist2)

data_gistic <- data_gistic %>%
  filter(CNV %in% c(-2, 2)) %>%
  mutate(CNV = if_else(CNV == 2, "Amp", "Del")) %>%
  left_join(wgs_data %>% select(Subject, Tumor_Barcode)) %>%
  mutate(Gene = GeneSymbol, Alteration = CNV, Type = "SCNA_Focal_Gene") %>%
  select(Subject, Tumor_Barcode, Gene, Alteration, Type)


# Fusion
load('SV_Fusion.RData', verbose = T)
genelist3 <- c(genelist, c('JARID2', 'PARK2', 'GRID2'))
data_fusion <- data_fusion %>% filter(Gene %in% genelist3)

## save R object
data_top <-
  bind_rows(
    data_mutation,
    data_gistic,
    data_fusion
  )

save(data_top, file = 'data_top.RData')

## combined all together
data_top <-
  bind_rows(
    data_mutation %>% mutate(Type = 'Alterations'),
    data_gistic %>% mutate(Type = 'Alterations'),
    data_fusion %>% mutate(Type = 'Alterations'),
  )

# Add data feature
load('MCN_data.RData', verbose = T)
load('VerifyBamID_idata.RData')
load('isochromosome_12p.RData')

data_feature1 <- idata %>%
  select(Subject, WGS_Inferred_Ancestry = POP) %>%
  left_join(wgs_data) %>%
  mutate(Gene = 'Ancestry', Type = 'Feature') %>%
  select(Subject, Tumor_Barcode, Gene, Alteration = WGS_Inferred_Ancestry, Type)

data_feature2 <- MCN_data %>%
  mutate(
    Alteration = if_else(WGD2 == 'WGD2', 'Multiple WGD', WGD),
    Gene = 'WGD',
    Type = 'Feature'
  ) %>%
  select(Subject, Tumor_Barcode, Gene, Alteration, Type)

data_feature3 <- isochromosome_12p %>%
  select(Tumor_Barcode, isochromsome) %>%
  left_join(wgs_data) %>%
  mutate(Gene = 'i(12p)', Type = 'Feature') %>%
  select(Subject, Tumor_Barcode, Gene, Alteration = isochromsome, Type)


data_feature <- data_feature1 %>%
  bind_rows(data_feature2) %>%
  bind_rows(data_feature3)

#
# data_feature <- sp_group_data2 %>%
#   filter(Tumor_Barcode %in% kras_samples) %>%
#   mutate(Type='Genomic_Feature',Gene='Group') %>%
#   left_join(wgs_groups_info) %>%
#   select(Subject,Tumor_Barcode,Gene,Alteration = SP_Group_New,Type)

## oncoplot
sample_level0 <- wgs_data %>%
  filter(Subtype != 'Seminoma') %>%
  pull(Tumor_Barcode)


# select only recurrent events --------------------------------------------
tmp <- data_mutation %>%
  filter(Tumor_Barcode %in% sample_level0) %>%
  select(Subject, Gene) %>%
  unique() %>%
  count(Gene, sort = T) %>%
  filter(n > 1) %>%
  pull(Gene)

data_mutation <- data_mutation %>% filter(Gene %in% tmp)

tmp <- data_gistic %>%
  filter(Tumor_Barcode %in% sample_level0) %>%
  select(Subject, Gene) %>%
  unique() %>%
  count(Gene, sort = T) %>%
  filter(n > 1) %>%
  pull(Gene)

data_gistic <- data_gistic %>% filter(Gene %in% tmp)

tmp <- data_fusion %>%
  filter(Tumor_Barcode %in% sample_level0) %>%
  select(Subject, Gene) %>%
  unique() %>%
  count(Gene, sort = T) %>%
  filter(n > 1) %>%
  pull(Gene)

data_fusion <- data_fusion %>% filter(Gene %in% tmp)


result_mutation <- oncoplot(
  data = data_mutation,
  landscape_colors = landscape_colors,
  sample_level0 = sample_level0,
  GeneSortOnly = TRUE
)


result_gistic <- oncoplot(
  data = data_gistic,
  landscape_colors = landscape_colors,
  sample_level0 = sample_level0,
  GeneSortOnly = TRUE
)


result_fusion <- oncoplot(
  data = data_fusion,
  landscape_colors = landscape_colors,
  sample_level0 = sample_level0,
  GeneSortOnly = TRUE
)

result_feature <- oncoplot(
  data = data_feature,
  landscape_colors = landscape_colors,
  sample_level0 = sample_level0,
  GeneSortOnly = TRUE
)

sample_new_level <- result_mutation$sample_level %>%
  left_join(
    result_gistic$sample_level
  ) %>%
  left_join(
    result_fusion$sample_level
  ) %>%
  left_join(
    result_feature$sample_level
  ) %>%
  arrange(Mutation_Driver, SCNA_Focal_Gene, Fusion, Feature) %>%
  pull(Tumor_Barcode)


result_mutation <- oncoplot(
  data = data_mutation,
  landscape_colors = landscape_colors,
  sample_level0 = sample_level0,
  sample_level = sample_new_level,
  tmar = -0.2,
  bmar = -0.05
  #p2_axis_hidden = TRUE,
  #p2_hidden = TRUE
)

result_gistic <- oncoplot(
  data = data_gistic,
  landscape_colors = landscape_colors,
  sample_level0 = sample_level0,
  sample_level = sample_new_level,
  tmar = -0.2,
  bmar = -0.05
  #p2_axis_hidden = TRUE,
  #p2_hidden = TRUE
)

result_fusion <- oncoplot(
  data = data_fusion,
  landscape_colors = landscape_colors,
  sample_level0 = sample_level0,
  sample_level = sample_new_level,
  GeneSortOnly = TRUE,
  tmar = -0.2,
  bmar = -0.05
  #p2_axis_hidden = TRUE,
  #p2_hidden = TRUE,
)

result_feature <- oncoplot(
  data = data_feature,
  landscape_colors = landscape_colors,
  sample_level0 = sample_level0,
  sample_level = sample_new_level,
  GeneSortOnly = TRUE,
  bmar = 0,
  p2_axis_hidden = TRUE,
  p2_hidden = TRUE,
  removeAlterbyColor = FALSE
)


oncoplot_final <- oncoplot_combined(
  result_mutation,
  result_gistic,
  result_fusion,
  result_feature
)


## for the legned
plegends <- plot_grid(
  result_mutation$oncoplot_legend,
  result_gistic$oncoplot_legend,
  result_fusion$oncoplot_legend,
  result_feature$oncoplot_legend
)

save_plot(
  filename = 'Genome_landscape_legends.pdf',
  plot = plegends,
  base_height = 6,
  base_width = 6,
  device = cairo_pdf
)


# Mutual exclusivity ------------------------------------------------------
load('data_top.RData')

tdata0 <- bind_rows(
  data_top %>%
    filter(Type == 'Mutation_Driver') %>%
    mutate(Alteration = 'Yes', Gene = paste0(Gene, '-Mut')) %>%
    select(Tumor_Barcode, Gene, Alteration) %>%
    unique(),

  data_top %>%
    filter(Type == 'Fusion') %>%
    mutate(Alteration = 'Yes', Gene = paste0(Gene, '-Fus')) %>%
    select(Tumor_Barcode, Gene, Alteration) %>%
    unique(),

  data_top %>%
    filter(Type == 'SCNA_Focal_Gene') %>%
    mutate(Gene = paste0(Gene, '-', Alteration)) %>%
    mutate(Alteration = 'Yes') %>%
    select(Tumor_Barcode, Gene, Alteration) %>%
    unique()
)


samlist <- wgs_data %>% filter(Subtype != 'Seminoma') %>% pull(Tumor_Barcode)

genelist <- tdata0 %>%
  filter(Tumor_Barcode %in% samlist) %>%
  count(Gene) %>%
  mutate(freq = n / length(samlist)) %>%
  filter(n > 5) %>%
  pull(Gene)

tdata <- tdata0 %>% filter(Tumor_Barcode %in% samlist, Gene %in% genelist)


# mutational exclusivity test
library(ggasym)

tdata <- tibble(Tumor_Barcode = samlist) %>%
  left_join(
    tdata %>% pivot_wider(names_from = Gene, values_from = Alteration)
  ) %>%
  mutate(across(-Tumor_Barcode, ~ ifelse(is.na(.) | . != "Yes", 0, 1)))

# Get the list of genes (excluding Tumor_Barcode column)
genes <- colnames(tdata)[-1]

# Create a function to perform Fisher's Exact Test and extract p-value and OR
fisher_test <- function(gene1, gene2, data) {
  tryCatch(
    {
      # Create a contingency table for the two genes
      contingency_table <- table(data[[gene1]], data[[gene2]])

      # Perform Fisher's Exact Test
      fisher_result <- fisher.test(contingency_table)

      # Return the p-value and odds ratio as a tibble
      tibble(
        #gene1 = gene1,
        #gene2 = gene2,
        p_value = fisher_result$p.value,
        odds_ratio = fisher_result$estimate
      )
    },
    error = function(e) {
      # In case of error, return NA for both p-value and odds ratio
      tibble(
        #gene1 = gene1,
        #gene2 = gene2,
        p_value = NA,
        odds_ratio = NA
      )
    }
  )
}
fisher_results <- NULL
# Use purrr to apply the fisher_test function to all unique pairs of genes
fisher_results <- expand.grid(gene1 = genes, gene2 = genes) %>%
  filter(gene1 != gene2) %>% # Exclude comparisons of the same gene
  distinct() %>%
  as_tibble() %>%
  mutate(result = map2(gene1, gene2, ~ fisher_test(.x, .y, tdata))) %>%
  unnest_wider(result) # Use unnest_wider to expand the 'result' column

# View the results with p-values and odds ratios
fisher_results <- fisher_results %>% asymmetrise(gene1, gene2)

fisher_results <- fisher_results %>%
  mutate(odds_ratio = if_else(p_value > 0.05, NA_real_, odds_ratio)) %>%
  mutate(p_value = if_else(p_value < 0.05, p_value, NA_real_)) %>%
  mutate(odds_ratio = if_else(odds_ratio == 0, 1 / 256, odds_ratio)) %>%
  mutate(odds_ratio = if_else(is.infinite(odds_ratio), 256, odds_ratio))

ggplot(fisher_results, aes(x = gene1, y = gene2)) +
  geom_asymmat(
    aes(fill_tl = -log10(p_value), fill_br = log2(odds_ratio)),
    col = 'white'
  ) +
  scale_fill_tl_gradient2(
    low = "white",
    high = "tomato",
    na.value = 'gray90',
    breaks = pretty_breaks(n = 5)
  ) +
  scale_fill_br_gradient2(
    low = 'purple3',
    high = 'orange3',
    midpoint = 0,
    mid = 'white',
    na.value = 'gray90',
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    #title = "Others",
    fill_tl = "Fisher's exact test\n-log10(p-value)",
    fill_br = "Fisher's exact test\nlog2(OR)"
  ) +
  theme_ipsum_rc() +
  theme(
    #panel.background = element_rect(fill = "grey75"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  panel_border(color = 'black', size = 0.5)

ggsave(filename = 'tmp.pdf', width = 6.5, height = 5, device = cairo_pdf)
