# # Install required packages (if not already installed)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("GSVA")
# BiocManager::install("msigdbr")

# Load required libraries
set_wd()
library(GSVA)
library(msigdbr)
libztw()
pdfhr()
pdfhr2()
# Load your RNA-Seq expression data
load('~/NIH-Work/R/annotable_ens90.RData', verbose = T)

protein_coding_list <- bind_rows(
  grch37 %>% filter(biotype == 'protein_coding') %>% select(symbol),
  grch38 %>% filter(biotype == 'protein_coding') %>% select(symbol)
) %>%
  unique() %>%
  pull(symbol)

tcga_studies <- str_remove_all(
  list.files(
    path = '~/NIH-Work/Earl_Stadtman_Inverstigator/Projects/cBioPortal_Studies/TCGA_PanCancer_Atlas/',
    pattern = '.tar.gz'
  ),
  '_.*'
) %>%
  unique()

cdata_all <- NULL
rdata_all <- NULL

for (study in tcga_studies) {
  cdata <- NULL
  rdata <- NULL

  rdata <- read_delim(
    paste0(
      '~/NIH-Work/Earl_Stadtman_Inverstigator/Projects/cBioPortal_Studies/TCGA_PanCancer_Atlas/',
      study,
      '_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt'
    ),
    delim = '\t',
    col_names = T
  )
  tmp <- NULL
  tmp <- rdata %>%
    count(Hugo_Symbol, sort = T) %>%
    filter(n > 1) %>%
    pull(Hugo_Symbol)
  rdata <- rdata %>%
    filter(!(Hugo_Symbol %in% tmp)) %>%
    pivot_longer(
      -c(Hugo_Symbol, Entrez_Gene_Id),
      names_to = 'TCGA_Barcode',
      values_to = 'RSEM'
    ) %>%
    filter(Hugo_Symbol %in% protein_coding_list) %>%
    mutate(Study = toupper(study)) %>%
    relocate(Study)

  rdata_all <- bind_rows(rdata_all, rdata)
  print(paste0(study, '.......Completed'))
}

save(rdata_all, file = 'TCGA-RSEM-protein-coding.RData')

barcode2study <- rdata_all %>% select(Study, TCGA_Barcode) %>% unique()


## pathway analysis by ssgsea
na_genes <- rdata_all %>% filter(is.na(RSEM)) %>% pull(Hugo_Symbol) %>% unique()
gene_expr <- rdata_all %>%
  filter(!(Hugo_Symbol %in% na_genes)) %>%
  select(TCGA_Barcode, Hugo_Symbol, RSEM) %>%
  pivot_wider(names_from = TCGA_Barcode, values_from = RSEM) #%>% filter(Study=='TGCT')

# Convert your data frame into a numeric matrix
gene_expr_mat <- as.matrix(gene_expr[, -1])
rownames(gene_expr_mat) <- gene_expr$Hugo_Symbol

# Obtain Hallmark or KEGG pathways (example with Hallmark)
#msigdbr_collections() %>% View()
m_df <- msigdbr(species = "Homo sapiens", collection = "H") # For Hallmark gene sets
#m_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") # For KEGG gene sets
# Format gene sets
pathway_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)
# Perform ssGSEA analysis
ssgsea_results <- gsva(
  gene_expr_mat,
  pathway_list,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE
)
# View the results
ssgsea_results_h <- as.data.frame(ssgsea_results) %>%
  rownames_to_column(var = 'Pathway') %>%
  as_tibble() %>%
  pivot_longer(-Pathway, names_to = 'TCGA_Barcode', values_to = 'ssGSEA_socre')


#m_df <- msigdbr(species = "Homo sapiens", collection = "H") # For Hallmark gene sets
m_df <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:KEGG_LEGACY"
) # For KEGG gene sets
# Format gene sets
pathway_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)
# Perform ssGSEA analysis
ssgsea_results <- gsva(
  gene_expr_mat,
  pathway_list,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE
)
# View the results
ssgsea_results_kegg <- as.data.frame(ssgsea_results) %>%
  rownames_to_column(var = 'Pathway') %>%
  as_tibble() %>%
  pivot_longer(-Pathway, names_to = 'TCGA_Barcode', values_to = 'ssGSEA_socre')

save(
  barcode2study,
  ssgsea_results_h,
  ssgsea_results_kegg,
  file = 'ssGSEA.RData'
)


# explore the kegg pathways -----------------------------------------------
load('wgs_data.RData')
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

ggsave('ssGSEA_TCGA_KEGG.pdf', width = 10, height = 8, device = cairo_pdf)


## between TGCT study...
pdata2 <- pdata %>% filter(!is.na(Subtype))

pdata3 <- pdata2 %>%
  group_by(Pathway) %>%
  do(tidy(wilcox.test(ssGSEA_socre ~ Subtype, data = .))) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

tmp <- pdata2 %>%
  group_by(Pathway, Subtype) %>%
  summarise(score = median(ssGSEA_socre)) %>%
  ungroup() %>%
  pivot_wider(names_from = Subtype, values_from = score) %>%
  mutate(FC = `Non-Seminoma` - Seminoma)

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
    aes(size = -log10(FDR), fill = `Non-Seminoma`),
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

ggsave('ssGSEA_TGCT_KEGG.pdf', width = 10, height = 8, device = cairo_pdf)

## polar plots across all TCGA study

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

ggsave(
  'ssGSEA_TCGA_KEGG_examples.pdf',
  width = 8,
  height = 6,
  device = cairo_pdf
)


# explore the hallmark pathways -----------------------------------------------

load('wgs_data.RData')
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

ggsave('ssGSEA_TCGA_HALLMARK.pdf', width = 8, height = 8, device = cairo_pdf)


## between TGCT study...
pdata2 <- pdata %>% filter(!is.na(Subtype))

pdata3 <- pdata2 %>%
  group_by(Pathway) %>%
  do(tidy(wilcox.test(ssGSEA_socre ~ Subtype, data = .))) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

tmp <- pdata2 %>%
  group_by(Pathway, Subtype) %>%
  summarise(score = median(ssGSEA_socre)) %>%
  ungroup() %>%
  pivot_wider(names_from = Subtype, values_from = score) %>%
  mutate(FC = `Non-Seminoma` - Seminoma)

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
    aes(size = -log10(FDR), fill = `Non-Seminoma`),
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

ggsave('ssGSEA_TGCT_HALLMARK.pdf', width = 9, height = 8, device = cairo_pdf)

## polar plots across all TCGA study

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

ggsave(
  'ssGSEA_TCGA_HALLMARK_examples.pdf',
  width = 8,
  height = 6,
  device = cairo_pdf
)


# Association with mutaitonal signature -----------------------------------
#hallmark
load('wgs_data.RData')
load('tgct_mutational_signatures.RData', verbose = T)
load('ssGSEA.RData', verbose = T)
load('tgct_clinical.RData')
load('ngspurity_data.RData')


tdata <- wgs_data %>%
  left_join(ngspurity %>% select(Tumor_Barcode, Tumor_Purity = BB_Purity)) %>%
  left_join(tgct_clinical %>% select(Subject, Age = age_at_diagnosis)) %>%
  left_join(
    tgct_activity_all %>%
      pivot_longer(
        -Tumor_Barcode,
        names_to = 'Signature',
        values_to = 'Exposure'
      )
  ) %>%
  mutate(TCGA_Barcode = str_sub(Tumor_Barcode, 1, 15)) %>%
  left_join(
    bind_rows(
      ssgsea_results_h %>% mutate(Group = 'HALLMARK'),
      ssgsea_results_kegg %>% mutate(Group = 'KEGG')
    )
  ) %>%
  filter(!is.na(Pathway))


tresult <- tdata %>%
  group_by(Signature, Group, Pathway) %>%
  do(
    tresult = safely(lm)(
      log2(Exposure + 1) ~ Tumor_Purity + ssGSEA_socre + Age + Subtype,
      data = .
    )
  ) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Signature, Group, Pathway, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(term == 'ssGSEA_socre') %>%
  arrange(p.value)


tresult <- tdata %>%
  group_by(Subtype, Signature, Group, Pathway) %>%
  do(
    tresult = safely(lm)(
      log2(Exposure + 1) ~ Tumor_Purity + ssGSEA_socre + Age,
      data = .
    )
  ) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Subtype, Signature, Group, Pathway, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(term == 'ssGSEA_socre') %>%
  arrange(p.value)


# difference between two subtypes
tdata %>%
  filter(Pathway == 'KEGG_PURINE_METABOLISM', Signature == 'SBS87') %>%
  do(tidy(wilcox.test(ssGSEA_socre ~ Subtype, data = .)))


tdata %>%
  filter(Pathway == 'KEGG_PURINE_METABOLISM', Signature == 'SBS87') %>%
  do(tidy(lm(
    ssGSEA_socre ~ Subtype + Tumor_Purity + ssGSEA_socre + Age,
    data = .
  )))


tdata %>%
  filter(Pathway == 'KEGG_PURINE_METABOLISM', Signature == 'SBS87') %>%
  ggplot(aes(Subtype, ssGSEA_socre)) +
  geom_boxplot() +
  facet_wrap(~Signature)

tdata %>%
  filter(Pathway == 'KEGG_PURINE_METABOLISM', Signature == 'SBS87') %>%
  ggplot(aes(log2(Exposure + 1), ssGSEA_socre)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~Subtype)

h_gene_sets <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:KEGG_LEGACY"
) # For KEGG gene sets
# Format gene sets
head(h_gene_sets)


# Extracted pathway related genes -----------------------------------------

pgenes <- c(
  "TPMT",
  "HPRT",
  "IMPDH1",
  "IMPDH2",
  "GMPS",
  "XDH",
  "XO",
  "AOX1",
  "ITPA",
  "ABCC4",
  "NUDT15",
  "NUDT1",
  "MTHFR",
  "IL6ST",
  "FSTL5",
  "PACSIN2",
  "MOCOS"
)

load('TCGA-RSEM-protein-coding.RData', verbose = T)


tmp <- rdata_all %>%
  filter(Study == 'TGCT', Hugo_Symbol %in% pgenes) %>%
  mutate(RSEM = log2(RSEM + 1))
tdata <- wgs_data %>%
  left_join(ngspurity %>% select(Tumor_Barcode, Tumor_Purity = BB_Purity)) %>%
  left_join(tgct_clinical %>% select(Subject, Age = age_at_diagnosis)) %>%
  mutate(TCGA_Barcode = str_sub(Tumor_Barcode, 1, 15)) %>%
  left_join(
    tgct_activity_all %>%
      pivot_longer(
        -Tumor_Barcode,
        names_to = 'Signature',
        values_to = 'Exposure'
      )
  ) %>%
  left_join(tmp) %>%
  filter(!is.na(RSEM))


tresult <- tdata %>%
  group_by(Signature, Hugo_Symbol) %>%
  do(
    tresult = safely(lm)(
      log2(Exposure + 1) ~ Tumor_Purity + RSEM + Age + Subtype,
      data = .
    )
  ) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Signature, Hugo_Symbol, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(term == 'RSEM') %>%
  arrange(p.value)


tresult <- tdata %>%
  group_by(Subtype, Signature, Hugo_Symbol) %>%
  do(
    tresult = safely(lm)(
      log2(Exposure + 1) ~ Tumor_Purity + RSEM + Age,
      data = .
    )
  ) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Subtype, Signature, Hugo_Symbol, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(term == 'RSEM') %>%
  arrange(p.value)


tdata %>%
  group_by(Signature, Hugo_Symbol) %>%
  mutate(Group = Exposure > 0) %>%
  do(tresult = safely(wilcox.test)(RSEM ~ Group, data = .)) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Signature, Hugo_Symbol, fit) %>%
  unnest(cols = c(fit)) %>%
  filter(Signature == 'SBS87') %>%
  ungroup() %>%
  arrange(p.value)


# # visualization ---------------------------------------------------------

genelev <- tdata %>%
  filter(Signature == 'SBS87') %>%
  group_by(Hugo_Symbol) %>%
  do(tresult = safely(lm)(RSEM ~ Subtype, data = .)) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Hugo_Symbol, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(str_detect(term, 'Subtype')) %>%
  arrange(estimate) %>%
  pull(Hugo_Symbol)


tresult <- tdata %>%
  filter(Signature == 'SBS87') %>%
  group_by(Hugo_Symbol) %>%
  do(tresult = safely(wilcox.test)(RSEM ~ Subtype, data = .)) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Hugo_Symbol, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  arrange(p.value) %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

tresult <- tdata %>%
  filter(Signature == 'SBS87') %>%
  mutate(RSEM = 2^RSEM - 1) %>%
  group_by(Subtype, Hugo_Symbol) %>%
  summarise(RSEM = median(RSEM)) %>%
  ungroup() %>%
  pivot_wider(names_from = Subtype, values_from = RSEM) %>%
  mutate(log2FC = log2(`Non-Seminoma` / Seminoma)) %>%
  left_join(tresult)


tresult %>%
  ggplot(aes(log2FC, -log10(FDR), fill = log2FC > 0)) +
  geom_point(pch = 21, size = 3) +
  ggrepel::geom_text_repel(aes(label = Hugo_Symbol)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = ncicolpal[1]) +
  geom_vline(xintercept = 0, linetype = 2, col = "gray20") +
  scale_fill_manual(values = as.character(subtypecol)) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = pretty_breaks(n = 7)) +
  labs(
    x = 'Non-Seminoma vs. Seminoma, log2(Fold Change)',
    y = 'Significance, -log10(FDR)'
  ) +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = F,
    ticks = T
  ) +
  theme(
    panel.spacing.x = unit(0.2, 'cm'),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'none',
    strip.text.x = element_text(hjust = 0.5, face = 'bold.italic', size = 10),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(linewidth = 0.8, fill = NA)
  )

ggsave(
  filename = 'thiopurine_genes_degs_overview.pdf',
  width = 6,
  height = 5,
  device = cairo_pdf
)


tdata %>%
  filter(Signature == 'SBS87') %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = genelev)) %>%
  ggplot(aes(Subtype, RSEM, fill = Subtype)) +
  #geom_violin(trim = T, size = 0.1) +
  geom_boxplot(
    width = 0.6,
    #fill = NA,
    color = 'black',
    outlier.shape = 21,
    outlier.size = 0.8,
    linewidth = 0.3
  ) +
  #ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="white",stroke=0.2)+
  facet_wrap(~Hugo_Symbol, nrow = 1) +
  labs(x = NULL, y = 'log2(RNA Seq V2 RSEM +1)') +
  scale_fill_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 7)) +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = 'Y',
    ticks = T
  ) +
  theme(
    axis.text.x = element_blank(),
    panel.spacing.x = unit(0.2, 'cm'),
    axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'right',
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    strip.text.x = element_text(hjust = 0.5, face = 'bold.italic', size = 10),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(linewidth = 0.5, fill = NA)
  ) +
  #panel_border(color = 'black',size=0.2) +
  # stat_pvalue_manual(
  #   data = result,
  #   col = 'col',
  #   label = "myformatted.p",
  #   size = 5
  # ) +
  scale_color_manual(values = c(ncicolpal[1], 'black')) +
  guides(color = 'none')

ggsave(
  filename = 'thiopurine_genes_degs.pdf',
  width = 12,
  height = 6,
  device = cairo_pdf
)

### association within Non-Seminoma

tresult <- tdata %>%
  filter(Signature == 'SBS87', Subtype == 'Non-Seminoma') %>%
  group_by(Hugo_Symbol) %>%
  do(
    tresult = safely(lm)(
      log2(Exposure + 1) ~ Tumor_Purity + RSEM + Age,
      data = .
    )
  ) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Hugo_Symbol, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(term == 'RSEM') %>%
  arrange(p.value)


tdata %>%
  filter(Signature == 'SBS87', Subtype == 'Non-Seminoma') %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = genelev)) %>%
  ggplot(aes(RSEM, log2(Exposure + 1))) +
  geom_point(
    pch = 21,
    stroke = 0.1,
    size = 2,
    fill = as.character(subtypecol)[2]
  ) +
  stat_smooth(method = "lm") +
  stat_cor(size = 5, col = ncicolpal[1]) +
  facet_wrap(~Hugo_Symbol, scale = 'free') +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_d3() +
  #scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  labs(x = 'log2(RNA Seq V2 RSEM +1)', y = 'log2(SBS87+1)') +
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
    panel.spacing.x = unit(0.2, 'cm'),
    panel.spacing.y = unit(0.7, 'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 14)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(
  filename = 'thiopurine_genes_SBS87_cor.pdf',
  width = 10,
  height = 10,
  device = cairo_pdf
)
