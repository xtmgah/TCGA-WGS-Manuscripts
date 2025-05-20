set_wd()
libztw()
pdfhr()
pdfhr2()

library(GenomicRanges)
# load PCAWG function
source('PCAWG_funs.R')

#load('../ZTW_functions.RData')

# load dataset ------------------------------------------------------------
load('wgs_data.RData')
load(
  '/Volumes/Zhang_Group/ZTW/TCGA-TGCT/MutationTimeR/tgct_mtvcf.RData',
  verbose = T
)
load('DP_info_data.RData')

barcode2spg <- wgs_data$Subtype
names(barcode2spg) <- wgs_data$Tumor_Barcode

## update tgct mtvcf ## use the dpcluster subclone mutations ##
tmplevel <- levels(tgct_mtvcf$CLS)
#c("clonal [early]", "clonal [late]",  "clonal [NA]","subclonal")
tgct_mtvcf <- tgct_mtvcf %>%
  mutate(CLS_old = CLS) %>%
  left_join(
    DP_info_data %>%
      select(Tumor_Barcode, MutationID = ID, Clone) %>%
      mutate(
        Clone = if_else(
          Clone == "Y",
          "clonal [NA]",
          if_else(Clone == "N", "subclonal", Clone)
        )
      )
  ) %>%
  mutate(
    CLS = if_else(is.na(CLS) | CLS == 'subclonal', Clone, as.character(CLS))
  )

save(tgct_mtvcf, file = 'tgct_mtvcf.RData')
save(tgct_mtvcf, tgct_mtbb, file = 'tgct_mt.RData')

# Summary the Clone and Subclone Mutations --------------------------------
clsdata <- tgct_mtvcf %>% count(Tumor_Barcode, CLS)

clsdata <- clsdata %>%
  group_by(Tumor_Barcode) %>%
  mutate(Freq = n / sum(n)) %>%
  ungroup()


clscolors <- pal_d3()(5)[c(3, 5, 1, 4)]
names(clscolors) <- levels(as.factor(clsdata$CLS))

clsdata %>%
  left_join(wgs_data) %>%
  filter(!is.na(CLS)) %>%
  ggplot(aes(CLS, Freq, fill = Subtype)) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.color = 'white',
    outlier.size = 2,
    outlier.stroke = 0.2
  ) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linewidth = 0.25, col = '#cccccc') +
  labs(x = "", y = "Proportion of total mutations", fill = 'Subtype') +
  scale_fill_manual(values = subtypecol, breaks = names(subtypecol)) +
  theme_ipsum_rc(
    axis_text_size = 12,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = F
  ) +
  theme(
    axis.text.y = element_text(angle = 30, vjust = 1, hjust = 1),
    axis.ticks.x = element_line(linewidth = .5),
    plot.margin = margin(4, 4, 4, -20)
  ) +
  theme(legend.position = 'top', legend.direction = 'horizontal') +
  panel_border(color = 'black') +
  coord_flip()
ggsave(
  filename = 'CLS_subtype.pdf',
  width = 6,
  height = 3.5,
  device = cairo_pdf
)

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

ggsave(
  filename = 'CLS_subtype2.pdf',
  width = 3.5,
  height = 5,
  device = cairo_pdf
)

clsdata %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  group_by(Subtype, Tumor_Barcode, CLS) %>%
  summarise(n = sum(n)) %>%
  mutate(freq = n / sum(n)) %>%
  group_by(Subtype, CLS) %>%
  summarise(median(freq)) %>%
  ungroup()

clsdata %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  filter(CLS == 'clonal [late]') %>%
  do(tidy(wilcox.test(Freq ~ Subtype, data = .)))

clsdata %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  filter(CLS == 'subclonal') %>%
  do(tidy(wilcox.test(Freq ~ Subtype, data = .)))


# adjust purity
load('ngspurity_data.RData')
load('VerifyBamID_idata.RData', verbose = T)
clsdata2 %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype, Subject)) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  left_join(idata %>% select(Subject, PC1, PC2)) %>%
  group_by(CLS) %>%
  #do(tidy(lm(Freq~Subtype+BB_Purity+PC1+PC2,data=.))) %>%
  do(tidy(lm(Freq ~ Subtype + BB_Purity, data = .))) %>%
  filter(term == "SubtypeNon-Seminoma")


#saveRDS(clsdata,file="clsdata.RData")
save(clsdata, clscolors, file = 'clsdata.RData')


# Number of Subclones -----------------------------------------------------
load('DP_info_data.RData', verbose = T)
load('ngspurity_data.RData')
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

ggsave(
  'subclone_miniccf_subtype.pdf',
  width = 6,
  height = 6,
  device = cairo_pdf()
)


# NRPCC>10 ----------------------------------------------------------------
tmp <- ngspurity %>% filter(NRPCC > 10) %>% pull(Tumor_Barcode)
pdata <- Cluster_info %>%
  filter(Clone == "N") %>%
  group_by(Tumor_Barcode) %>%
  summarise(N = n(), CCF = min(Cluster_CCF)) %>%
  right_join(
    wgs_data %>%
      select(Tumor_Barcode, Subtype) %>%
      filter(Tumor_Barcode %in% tmp)
  ) %>%
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

ggsave(
  'subclone_miniccf_subtype_NRPCC10.pdf',
  width = 6,
  height = 6,
  device = cairo_pdf()
)


# Mutational Signature Difference for ID2 paper -----------------------------------------
load('../Signature_ludmil3/Mutation_Signature_Probability_SBS.RData')
load('../Signature_ludmil3/sigcol.RData')

tdata <- tgct_mtvcf %>%
  select(Tumor_Barcode, MutationID, CLS) %>%
  left_join(
    Mutation_Signature_Probability_SBS %>%
      select(Tumor_Barcode, MutationID = ID, starts_with('SBS'))
  )

save(tdata, file = 'CLS2Signature.RData')
load('CLS2Signature.RData')
load('../Subtype_data.RData')

tdata2 <- tdata %>%
  filter(Tumor_Barcode %in% hq_samples2, !is.na(SBS1)) %>%
  left_join(Subtype_data2 %>% select(Tumor_Barcode, Subtype)) %>%
  pivot_longer(cols = -c(Subtype, Tumor_Barcode, MutationID, CLS)) %>%
  group_by(Subtype, CLS, name) %>%
  summarise(value = sum(value, na.rm = T)) %>%
  ungroup()

myggstyle()
tdata2 %>%
  filter(!is.na(CLS)) %>%
  mutate(
    name = factor(name, levels = rev(names(sigcol))),
    CLS = fct_rev(CLS)
  ) %>%
  ggplot(aes(CLS, value, fill = name)) +
  geom_col(position = "fill", width = 0.7, col = 'gray', linewidth = 0.2) +
  facet_wrap(~Subtype, ncol = 1) +
  labs(x = "", y = "Proportion of total mutations", fill = 'Signature') +
  scale_y_continuous(
    breaks = pretty_breaks(),
    expand = expansion(add = c(0, 0))
  ) +
  scale_fill_manual(
    values = sigcol,
    na.value = '#cccccc',
    breaks = names(sigcol)
  ) +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    axis = 'Y',
    axis_col = 'black'
  ) +
  theme(
    panel.spacing = unit(0.15, "lines"),
    strip.text.x = element_text(hjust = 0.5, face = 'bold')
  ) +
  panel_border(color = 'black') +
  coord_flip()

ggsave(
  filename = 'Mutation_clonality_proportion.pdf',
  width = 10,
  height = 6.5,
  device = cairo_pdf
)


# By Chromosome
clsdata <- tgct_mtvcf %>%
  mutate(chr = str_remove(MutationID, ":.*")) %>%
  count(Tumor_Barcode, chr, CLS) %>%
  filter(chr %in% c(1:22, "X", "Y"))

clsdata <- clsdata %>%
  group_by(Tumor_Barcode, chr) %>%
  mutate(Freq = n / sum(n)) %>%
  ungroup()

clscolors <- pal_d3()(5)[c(3, 5, 1, 4)]
names(clscolors) <- levels(clsdata$CLS)

clsdata %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  filter(!is.na(CLS)) %>%
  #mutate(chr = if_else(chr %in% "X",paste0(Gender,"-",chr),chr)) %>%
  filter(chr %in% c(1:22, "X")) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X"))) %>%
  ggplot(aes(CLS, Freq, fill = CLS)) +
  geom_boxplot(outlier.size = 1) +
  facet_grid(Subtype ~ chr) +
  labs(x = "", y = "Proportion of total mutations") +
  scale_fill_manual(values = clscolors) +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16
  ) +
  theme(panel.spacing = unit(0.1, "lines"), axis.text.x = element_blank()) +
  panel_border(color = 'black')

ggsave(
  filename = 'CLS_Subtype_chrs.pdf',
  plot = p,
  width = 16,
  height = 8,
  device = cairo_pdf
)


# Driver Gene Frequency by Clonality --------------------------------------
load('tgct_maf.RData', verbose = T)
#load('tgct_driver_mutations.RData')
load('intogene_drivers.RData', verbose = T)
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

tmp <- tmpdata %>%
  count(Subtype, Hugo_Symbol, CLS)
#filter(Subtype != 'Others')

tmp <- tmp %>%
  left_join(
    wgs_data %>% count(Subtype, name = 'size')
  ) %>%
  mutate(Freq = n / size) %>%
  arrange(desc(Freq))

genelist <- tmp %>%
  group_by(Hugo_Symbol) %>%
  summarise(n = sum(n)) %>%
  arrange(desc(n)) %>%
  filter(n > 1) %>%
  pull(Hugo_Symbol)

tmp %>%
  filter(Hugo_Symbol %in% genelist) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = genelist)) %>%
  ggplot(aes(Hugo_Symbol, Freq, fill = CLS)) +
  geom_bar(position = 'stack', stat = 'identity', size = 0.1) +
  facet_wrap(~Subtype, ncol = 1, scales = 'free') +
  scale_fill_manual(values = clscolors) +
  scale_y_continuous(breaks = pretty_breaks(), labels = percent_format()) +
  labs(x = "", y = "Driver mutation freqeuncy\n", fill = "") +
  theme_ipsum_rc(
    axis_text_size = 12,
    axis_title_just = 'm',
    axis_title_size = 16
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = 'Driver_Gene_Clonality.pdf',
  width = 7,
  height = 6,
  device = cairo_pdf
)

tmp %>%
  filter(Hugo_Symbol %in% genelist) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = genelist)) %>%
  ggplot(aes(Hugo_Symbol, Freq, fill = CLS)) +
  geom_bar(position = 'fill', stat = 'identity', size = 0.1) +
  facet_wrap(~Subtype, ncol = 1, scales = 'free') +
  scale_fill_manual(values = clscolors) +
  scale_y_continuous(breaks = pretty_breaks(), labels = percent_format()) +
  labs(x = "", y = "Proportion\n", fill = "") +
  theme_ipsum_rc(
    axis_text_size = 12,
    axis_title_just = 'm',
    axis_title_size = 16
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = 'Driver_Gene_Clonality2.pdf',
  width = 7,
  height = 6,
  device = cairo_pdf
)

# proportion of clonal mutations vs other mutation by fisher exact test
tmpdata %>%
  mutate(CLS = if_else(CLS == 'clonal [early]', 'Yes', 'No')) %>%
  select(Tumor_Barcode, Hugo_Symbol, CLS, Subtype) %>%
  # unique() %>%
  # count(Tumor_Barcode,Hugo_Symbol,sort=T) %>%
  mutate(CLS = as.factor(CLS), Subtype = as.factor(Subtype)) %>%
  group_by(Hugo_Symbol) %>%
  do(tidy(fisher.test(.$CLS, .$Subtype))) %>%
  arrange(p.value)


# Examples: early to late clonal mutation spectrum changes ----------------
# source('~/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Codes/Sigvisualfunc.R')
# load('~/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Database/Signature/signature_refsets.RData')
# profile0 <- signature_refsets %>% filter(Signature_set_name=='COSMIC_v3.2_Signatures_GRCh37_SBS96') %>% select(Signature_name,MutationType,Contribution) %>%
#   pivot_wider(names_from = Signature_name,values_from = Contribution)
#
# load('../DP_info_data.RData')
# tmp <- tgct_mtvcf %>% select(Tumor_Barcode,MutationID,CLS)
#
# mtype_clonality <- DP_info_data %>%
#   select(Tumor_Barcode,MutationID=ID,MutationType=mutType) %>%
#   left_join(tmp)
#
# # add clonal
# mtype_clonality <- mtype_clonality %>% filter(!is.na(CLS),CLS!='subclonal') %>% mutate(CLS="clonal") %>% bind_rows(mtype_clonality)
#
# mtype_clonality <- mtype_clonality %>% filter(Tumor_Barcode %in% hq_samples2, !is.na(CLS))
# mtype_clonality <- mtype_clonality %>% count(Tumor_Barcode,MutationType,CLS)
#
# profile1 <- mtype_clonality %>% mutate(sample=paste0(Tumor_Barcode,"_",CLS)) %>% select(sample,MutationType,n) %>% pivot_wider(names_from = sample,values_from = n)
# profile1 <- profile1 %>% filter(MutationType %in% profile0$MutationType)
#
# cos_sim_res <- cos_sim_df(profile1,profile0)
#
# cos_sim_res <- cos_sim_res %>% pivot_longer(cols = -rowname) %>% group_by(rowname) %>% arrange(desc(value)) %>% slice(1) %>% ungroup()
#
# mtype_clonality <- mtype_clonality %>% group_by(Tumor_Barcode,CLS) %>% summarise(n=sum(n,na.rm = T)) %>% ungroup()%>% mutate(sample=paste0(Tumor_Barcode,"_",CLS)) %>% left_join(cos_sim_res %>% select(sample=rowname,Signature=name,Cosine=value))
#
# # early clone to later clonal mutation spectrum changes
#
# tmp <- mtype_clonality %>%
#   filter(CLS %in% c('clonal [early]','clonal [late]')) %>%
#   filter(n>100,!is.na(Cosine),Cosine>0.7)
#
# tmplist <- tmp %>% count(Tumor_Barcode) %>% filter(n>1) %>% pull(Tumor_Barcode)
#
# tmp <- tmp %>% filter(Tumor_Barcode %in% tmplist) %>%
#   left_join(wgs_data %>% select(Tumor_Barcode,Subtype)) %>%
#   mutate(Value=paste0(Signature,"_",round(Cosine,3))) %>%
#   select(Tumor_Barcode,Subtype,CLS,Value) %>%
#   pivot_wider(names_from = CLS,values_from = Value)
#
# colnames(tmp)[3:4] <- c('early','late')
#
# tmp %>%
#   separate(col = early,into = c('early_sig','early_cos'),sep = "_",convert = T) %>%
#   separate(col = late,into = c('late_sig','late_cos'),sep = "_",convert = T) %>%
#   filter(early_sig != late_sig) %>%
#   filter(!(early_sig %in% c('SBS5','SBS40') & late_sig %in% c('SBS5','SBS40'))) %>%
#   View()
#
# barcode <- "NSLC-1180-T01"
# barcode <- "NSLC-1132-T01"
# barcode <- "NSLC-1177-T01"
# barcode <- "NSLC-1219-T01"
# barcode <- "NSLC-1107-T01"
# barcode <- "NSLC-1192-T01"
# barcode <- "NSLC-0602-T01"
# samples <- paste0(barcode,"_",c('clonal [early]','clonal [late]'))
#
# profilex <- profile1 %>% select(MutationType,samples[1])
# plot_sbs_96_profile(data = profilex,samplename = samples[1],filename = 'tmp1.pdf')
#
# profilex <- profile1 %>% select(MutationType,samples[2])
# plot_sbs_96_profile(data = profilex,samplename = samples[2],filename = 'tmp2.pdf')
#
#
# # clonal to subclonal mutation spectrum changes
# tmp <- mtype_clonality %>%
#   filter(CLS %in% c('clonal','subclonal', '')) %>%
#   filter(n>100,!is.na(Cosine),Cosine>0.8)
#
# tmplist <- tmp %>% count(Tumor_Barcode) %>% filter(n>1) %>% pull(Tumor_Barcode)
#
# tmp <- tmp %>% filter(Tumor_Barcode %in% tmplist) %>%
#   left_join(wgs_data %>% select(Tumor_Barcode,Subtype)) %>%
#   mutate(Value=paste0(Signature,"_",round(Cosine,3))) %>%
#   select(Tumor_Barcode,Subtype,CLS,Value) %>%
#   pivot_wider(names_from = CLS,values_from = Value)
#
# tmp %>%
#   separate(col = clonal,into = c('clonal_sig','clonal_cos'),sep = "_",convert = T) %>%
#   separate(col = subclonal,into = c('subclonal_sig','subclonal_cos'),sep = "_",convert = T) %>%
#   filter(clonal_sig != subclonal_sig) %>%
#   filter(!(clonal_sig %in% c('SBS5','SBS40') & subclonal_sig %in% c('SBS5','SBS40'))) %>%
#   View()
#
# barcode <- "NSLC-1219-T01"
# barcode <- "NSLC-0331-T01"
#
# barcode <- "NSLC-0405-T01"
#
# barcode <- "NSLC-0509-T01"
#
# barcode <- "NSLC-0680-T01"
#
# barcode <- "16MH0041_P2"
# barcode <- "NSLC-1168-T01"
# samples <- paste0(barcode,"_",c('clonal','subclonal'))
#
# profilex <- profile1 %>% select(MutationType,samples[1])
# plot_sbs_96_profile(data = profilex,samplename = samples[1],filename = 'tmp1.pdf')
#
# profilex <- profile1 %>% select(MutationType,samples[2])
# plot_sbs_96_profile(data = profilex,samplename = samples[2],filename = 'tmp2.pdf')

# Overview of the molecular timing distributions of copy number gain --------
timedGains <- tgct_mtbb %>%
  filter(!is.na(time)) %>%
  filter(chr != "Y") %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  mutate(length = abs(endpos - startpos)) %>%
  filter(length > 10000)

timedGains.all = timedGains
timedGains.all$chr = "all"
timedGains = rbind(
  timedGains,
  timedGains.all,
  timedGains %>% filter(type == "Bi-allelic Gain (WGD)") %>% mutate(chr = "WGD")
)

timedGains <- timedGains %>%
  mutate(prop = NA_real_, start_segment = NA_real_, end_segment = NA_real_)
# Add segment identifier
timedGains$segId = paste0(
  timedGains$Tumor_Barcode,
  "_",
  timedGains$chr,
  "_",
  timedGains$startpos,
  "_",
  timedGains$endpos,
  "_",
  timedGains$length,
  "_",
  timedGains$type
)


for (c in c(1:22, "X", "all", "WGD")) {
  for (h in unique(wgs_data$Subtype)) {
    timedGains_CH <- timedGains %>% filter(chr == c, Subtype == h)
    # Define proportion of events - this is the size of the pie chart
    n_samples = length(unique(timedGains_CH$Tumor_Barcode))
    n_samples_total = wgs_data %>%
      count(Subtype) %>%
      filter(Subtype == h) %>%
      pull(n)
    prop = sqrt((n_samples / n_samples_total) / pi)
    timedGains$prop[timedGains$chr == c & timedGains$Subtype == h] = prop

    # Define start/end boundaries for segments
    timedGains_CH = timedGains_CH[order(-timedGains_CH$time), ]
    timedGains_CH$start_segment = cumsum(rep(
      1 / nrow(timedGains_CH),
      nrow(timedGains_CH)
    )) -
      1 / nrow(timedGains_CH)
    timedGains_CH$end_segment = cumsum(rep(
      1 / nrow(timedGains_CH),
      nrow(timedGains_CH)
    ))

    timedGains$start_segment[
      timedGains$chr == c & timedGains$Subtype == h
    ] = timedGains_CH$start_segment[match(
      timedGains$segId[timedGains$chr == c & timedGains$Subtype == h],
      timedGains_CH$segId
    )]

    timedGains$end_segment[
      timedGains$chr == c & timedGains$Subtype == h
    ] = timedGains_CH$end_segment[match(
      timedGains$segId[timedGains$chr == c & timedGains$Subtype == h],
      timedGains_CH$segId
    )]
  }
}

# Sort cancer types by median pi0
medianTime = aggregate(timedGains$time, list(timedGains$Subtype), median)
timedGains$Subtype = factor(
  timedGains$Subtype,
  medianTime[order(medianTime$x), ]$Group.1
)
# Remove cancer types with fewer than 100 segments, 15 samples, and mean pie size < 0.2
timedGains = subset(
  timedGains,
  Subtype %in% names(which(table(timedGains$Subtype) > 100))
)
# quiet.ttypes = c(names(which(table(summaryTable$Subtype) < 20)))
mean_sizes = aggregate(
  timedGains$prop[timedGains$chr != "all"],
  list(timedGains$Subtype[timedGains$chr != "all"]),
  mean
)
# quiet.ttypes = unique(c(quiet.ttypes,  as.character(mean_sizes[which(mean_sizes$x < 0.2),]$Group.1)))

theme_set(theme_grey(base_family = 'Roboto Condensed'))
p <- ggplot(data = timedGains) + #[!timedGains$Subtype %in% quiet.ttypes,] #%>% sample_frac(siz=0.01)
  facet_grid(Subtype ~ factor(chr, levels = c(1:22, "X", "all", "WGD"))) +
  coord_polar(theta = "x") +
  geom_rect(aes(
    xmin = start_segment,
    xmax = end_segment,
    ymin = 0,
    ymax = prop,
    fill = time
  )) +
  scale_fill_distiller(
    palette = "PRGn",
    name = "Molecular time\n\n",
    breaks = c(0, 0.25, 0.5, 0.75, 1.00),
    limits = c(0, 1)
  ) +
  theme(strip.background = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, size = 20, hjust = 0)) +
  theme(strip.text.x = element_text(size = 15)) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank(), panel.grid = element_blank()) +
  theme(
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "gray", fill = NA, size = 0.25)
  ) +
  theme(panel.background = element_blank()) +
  theme(
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  ) +
  theme(legend.position = "bottom", legend.direction = 'horizontal') +
  guides(fill = guide_colorbar(barwidth = 12, barheight = 0.8)) +
  #ggtitle("Mutational timing estimates - 2 mutation threshold") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  file = 'time_clonal_copy_number_gains.png',
  plot = p,
  width = 25,
  height = 15
)
# #ggsave(file='time_clonal_copy_number_gains.pdf',plot = p,width = 25,height = 4,device = cairo_pdf)
# theme_set(theme_grey(base_family = 'Roboto Condensed'))
# p <- ggplot(data=timedGains) + #[!timedGains$Subtype %in% quiet.ttypes,] #%>% sample_frac(siz=0.01)
#   facet_grid(Subtype~chr) +
#   coord_polar(theta="x") +
#   geom_rect(aes(xmin=start_segment, xmax=end_segment, ymin=0, ymax=prop, fill=time)) +
#   scale_fill_distiller(palette = "PRGn", name="Molecular time\n\n", breaks=c(0, 0.25, 0.5, 0.75, 1.00), limits=c(0,1)) +
#   theme(strip.background = element_blank()) +
#   theme(strip.text.y=element_text(angle=0, size=20, hjust=0)) +
#   theme(strip.text.x=element_text(size=15)) +
#   theme(axis.text=element_blank()) +
#   theme(axis.ticks=element_blank(), panel.grid=element_blank()) +
#   theme(panel.spacing=unit(0, "lines"),
#         panel.border=element_rect(color="gray", fill=NA, size=0.25)) +
#   theme(panel.background=element_blank()) +
#   theme(legend.text=element_text(size=15), legend.title=element_text(size=15)) +
#   theme(legend.position="bottom",legend.direction = 'horizontal') +
#   guides(fill = guide_colorbar(barwidth=12, barheight=0.8)) +
#   #ggtitle("Mutational timing estimates - 2 mutation threshold") +
#   theme(plot.title = element_text(hjust = 0.5))
#
# ggsave(file='time_clonal_copy_number_gains_all.png',plot = p,width = 8,height = 15)

summary(timedGains$prop)

# # chr21/22 vs others
# my_comparisons <- list(c(21,"Others"),c(22,"Others"))
#
# stat.test <- timedGains %>%
#   filter(Subtype!="Others",chr %in% c(1:22,"X")) %>%
#   mutate(Group=if_else(chr %in% c(22,21),chr,'Others')) %>%
#   left_join(Subtype_data2) %>%
#   group_by(Subtype) %>%
#   rstatix::wilcox_test(time ~ Group, comparisons = my_comparisons) %>%
#   rstatix::add_xy_position(x = "Group") %>%
#   mutate(myformatted.p = sprintf("P = %.2e",p)) %>%
#   mutate(Group=NA)
#
# stat.test$y.position <- rep(c(1.05,1.1),3)
#
# timedGains %>%
#   filter(Subtype!="Others",chr %in% c(1:22,"X")) %>%
#   mutate(Group=if_else(chr %in% c(22,21),chr,'Others')) %>%
#   left_join(Subtype_data2) %>%
#   ggplot(aes(Group,time,fill=Group))+
#   #ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="white",stroke=0.2)+
#   geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.6,size=0.4)+
#   facet_wrap(~Subtype)+
#   scale_fill_bmj()+
#   scale_x_discrete(breaks=c(21,22,"Others"),labels=c('Chr21','Chr22','Other Chromosomes'))+
#   scale_y_continuous(breaks = pretty_breaks(n = 7))+
#   labs(y='Molecular time',x=NULL)+
#   theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
#   theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),,panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,0,4),legend.position = 'none',plot.title = element_text(hjust = 0.5,size=14))+
#   guides(fill="none")+
#   panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
#   stat_pvalue_manual(stat.test, label = "myformatted.p",color = ncicolpal[1],size = 3.5)
#
# ggsave(filename = 'time_clonal_copy_number_gains_chr21_22.pdf',width = 7,height = 6,device = cairo_pdf)
#
# # EUS chrx vs others
# my_comparisons <- list(c("X","Others"))
#
# stat.test <- timedGains %>%
#   filter(Subtype=="S_U",chr %in% c(1:22,"X")) %>%
#   mutate(Group=if_else(chr %in% c('X'),chr,'Others')) %>%
#   rstatix::wilcox_test(time ~ Group, comparisons = my_comparisons) %>%
#   rstatix::add_xy_position(x = "Group") %>%
#   mutate(myformatted.p = sprintf("P = %.2e",p)) %>%
#   mutate(Group=NA)
#
# stat.test$y.position <- c(1.05)
#
# timedGains %>%
#   filter(Subtype=="S_U",chr %in% c(1:22,"X")) %>%
#   mutate(Group=if_else(chr %in% c('X'),chr,'Others')) %>%
#   mutate(Group=factor(Group,levels=c('X','Others'))) %>%
#   ggplot(aes(Group,time,fill=Group))+
#   #ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="white",stroke=0.2)+
#   geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.6,size=0.4)+
#   scale_fill_manual(values = pal_bmj()(4)[c(4,3)])+
#   scale_x_discrete(breaks=c("X","Others"),labels=c('ChrX','Other Chromosomes'))+
#   scale_y_continuous(breaks = pretty_breaks(n = 7))+
#   labs(y='Molecular time',x=NULL,title='EU_S')+
#   theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
#   theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),,panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,0,4),legend.position = 'none',plot.title = element_text(face = 1,hjust = 0.5,size=14))+
#   guides(fill="none")+
#   panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
#   stat_pvalue_manual(stat.test, label = "myformatted.p",color = ncicolpal[1],size = 3.5)
#
# ggsave(filename = 'time_clonal_copy_number_gains_chrX.pdf',width = 2,height = 6.5,device = cairo_pdf)

# WGS/ND timming ----------------------------------------------------------
tgct_mtbb0 <- tgct_mtbb
library(GenomicRanges)

# timing_df <- do.call(rbind, lapply(tgct_mtbb0$timing_param, as.data.frame))
# tgct_mtbb0 <- cbind(tgct_mtbb0, timing_df)

finalBB <- GRanges(
  seqnames = tgct_mtbb0$chr,
  ranges = IRanges(tgct_mtbb0$startpos, tgct_mtbb0$endpos),
  major_cn = tgct_mtbb0$clone_nMaj,
  minor_cn = tgct_mtbb0$clone_nMin,
  total_cn = tgct_mtbb0$clone_nMin + tgct_mtbb0$clone_nMaj,
  clonal_frequency = tgct_mtbb0$clonal_frequency
  #timing_param = tgct_mtbb0$timing_param
)

mcols(finalBB) <- cbind(mcols(finalBB), tgct_mtbb0[, c(1, 8, 13:21)])
finalBB <- split(finalBB, mcols(finalBB)$Tumor_Barcode)

# WGD Prelim --------------------------------------------------------------
#Final ploidy, weighted if subclonal CN

finalPloidy <- sapply(finalBB, averagePloidy)
names(finalPloidy) <- names(finalBB)
#Final homozygousity, weighted if subclonal CN # minor_cn == 0
finalHom <- sapply(finalBB, averageHom)
names(finalHom) <- names(finalBB)

isWgd <- .classWgd(finalPloidy, finalHom)
table(isWgd)
plot(
  finalHom,
  finalPloidy,
  col = .classWgd(finalPloidy, finalHom) + 1,
  xlim = c(0, 1)
)

# compared to know WGD method
load('ngspurity_data.RData')
wgdinfo <- tibble(WGD = isWgd, Tumor_Barcode = names(isWgd)) %>%
  left_join(
    ngspurity %>% dplyr::select(Tumor_Barcode, WGD_Status = MCN_WGD)
  ) %>%
  dplyr::rename(PCAWG_WGD = WGD) %>%
  left_join(
    tibble(PCAWG_Hom = finalHom, Tumor_Barcode = names(finalHom))
  ) %>%
  left_join(
    tibble(PCAWG_Ploidy = finalPloidy, Tumor_Barcode = names(finalPloidy))
  )

saveRDS(wgdinfo, file = 'wgdinfo.RData')
#wgdinfo %>% select(PCAWG_WGD,WGD_Status) %>% table()

# WGD_Status
# WGD     nWGD WGD
# FALSE  449   0
# TRUE    21 747

# update with for loop
fracGenomeWgdComp <- t(sapply(finalBB, function(bb) {
  fgw <- try(fractionGenomeWgdCompatible(bb))
  if (class(fgw) != 'try-error') fgw else rep(NA, 10)
}))
rownames(fracGenomeWgdComp) <- names(finalBB)


wgdStar <- factor(
  rep(1, nrow(fracGenomeWgdComp)),
  levels = 0:3,
  labels = c("unlikely", "uninformative", "likely", "very likely")
)
wgdStar[
  fracGenomeWgdComp[, "avg.ci"] <= 0.75 &
    fracGenomeWgdComp[, "nt.total"] / chrOffset["M"] >= 0.33
] <- "likely"
wgdStar[
  fracGenomeWgdComp[, "nt.wgd"] / fracGenomeWgdComp[, "nt.total"] < 0.66
] <- "unlikely"
wgdStar[
  wgdStar == "likely" &
    fracGenomeWgdComp[, "nt.wgd"] / fracGenomeWgdComp[, "nt.total"] > 0.8 &
    fracGenomeWgdComp[, "sd.wgd"] < 0.1 &
    fracGenomeWgdComp[, "nt.total"] / chrOffset["M"] > 0.5
] <- "very likely"
names(wgdStar) <- names(finalBB)
prop.table(table(wgdStar[!isWgd]))

wgdPoss <- !isWgd & 2.5 - 1.5 * finalHom <= finalPloidy

wgdStat <- factor(
  wgdPoss + 2 * isWgd - wgdPoss * isWgd,
  labels = c("absent", "possible", "present")
)
table(wgdStat, wgdStar)


# Example of display very early timing on chromosomes --------------------
# display very early timing on chromosomes 7, 20 and 21. Plot a few.
w <- which(
  fracGenomeWgdComp[, "time.wgd"] < 0.1 &
    fracGenomeWgdComp[, "nt.total"] / chrOffset["M"] > 0.1 &
    !isWgd
)
wgs_data %>%
  filter(Tumor_Barcode %in% names(w)) %>%
  select(Subtype, Tumor_Barcode)

#Explore the early timing of chromosome 7 a bit more deeply. Find all GBMs with +7
w <- w[sapply(
  finalBB[w],
  function(bb)
    sum(
      width(bb)[as.logical(seqnames(bb) == 7) & bb$total_cn >= 3],
      na.rm = TRUE
    ) /
      width(refLengths[7]) >
      0.8
)]
#Tabulate the number of point mutations preceding +7

t <- t(sapply(w, function(ww) {
  finalSnv[[ww]] -> v
  v <- v[seqnames(v) == 7]
  v <- v[which(info(v)$MajCN <= 4 & info(v)$MajCN > 1)]
  n <- sum(info(v)$MutCN == info(v)$MajCN, na.rm = TRUE)
  c(pre_7 = n, total = nrow(v))
}))
rownames(t) <- names(finalSnv)[w]

# Less than 3 early

table(t[, 1] <= 3)


# Temporal distribution of chromosomal gains ------------------------------
aggregatePerChromosome <- function(bb, isWgd = FALSE) {
  .aggregateSegments <- function(m) {
    #m <- mcols(bb)
    t <- weighted.mean(m$time, m$n.snv_mnv, na.rm = TRUE)
    n <- sum(m$n.snv_mnv[!is.na(m$time)], na.rm = TRUE)
    sd <- sd(m$time, na.rm = TRUE)
    ci <- weighted.mean(m$time.up - m$time.lo, m$n.snv_mnv, na.rm = TRUE)
    w <- sum(m$width[!is.na(m$time)], na.rm = TRUE)
    c(time = t, n = n, sd = sd, ci = ci, w = w)
  }
  #   if(!isWgd){
  s <- split(
    as.data.frame(bb)[, c("time", "time.up", "time.lo", "n.snv_mnv", "width")],
    seqnames(bb)
  )
  r <- t(sapply(s, .aggregateSegments))
  r <- r[c(1:22, "X"), ]
  #   }else{
  w <- .aggregateSegments(as.data.frame(bb))
  r <- rbind(r, WGD = w)
  #   }
  return(r)
}


library(parallel)
allChrAgg <- simplify2array(mclapply(
  finalBB,
  aggregatePerChromosome,
  mc.cores = 2
))

t <- allChrAgg[1:23, "time", !isWgd]
t[allChrAgg[1:23, "w", !isWgd] < diff(chrOffset)[1:23] * .33] <- NA

s <- split(as.data.frame(t(t)), barcode2spg[colnames(t)])
n <- 10


at <- function(x, n) {
  if (sum(!is.na(x)) < 3) return(rep(sum(!is.na(x)) / n, n))
  bw = if (sum(!is.na(x)) < 6) 0.5 else "nrd0"
  d <- density(
    x,
    n = n,
    from = 1 / n / 2,
    to = 1 - 1 / n / 2,
    bw = bw,
    na.rm = TRUE
  )
  d$y / sum(d$y) * d$n
}

allChrCancerHist <- sapply(s, apply, 2, at, n = n, simplify = "array")
u <- split(
  data.frame(WGD = allChrAgg["WGD", "time", isWgd]),
  barcode2spg[names(allChrAgg["WGD", "time", isWgd])]
)
wgdCancerHist <- sapply(
  u,
  function(x)
    if (nrow(x) > 0) {
      at(x$WGD, n = n)
    } else {
      rep(0, n)
    },
  simplify = "array"
)
allChrCancerHist <- abind::abind(
  allChrCancerHist,
  All = sapply(sapply(s, as.matrix), at, n = n, simplify = "array") / 23 * 5,
  WGD = wgdCancerHist,
  along = 2
)


prgn <- RColorBrewer::brewer.pal(11, "PRGn")
set1 <- RColorBrewer::brewer.pal(9, "Set1")
col <- colorRampPalette(set1[c(4, 9, 3)])(n)

p <- 0
v <- table(barcode2spg[names(allChrAgg["WGD", "time", ])])
h <- (allChrCancerHist + p) /
  rep(v + p, each = prod(dim(allChrCancerHist)[1:2]))
h <- aperm(h, c(2, 3, 1))

conflicts_prefer(base::setdiff)
a <- colMeans(h[c("All", "WGD"), , ] * c(23 / 5, 1)) %*%
  1:n /
  asum(h * c(23 / 5, 1), c(1, 3))
o <- order(-a)
h <- h[, o, ]
w <- v[o] >= 15 & apply(h, 2, max) > 0.05 * 8 / n
h <- h[, w, ]
h <- h[, -3, ]

m <- 0.02

pdf('time_per_chromosome.pdf', width = 16, height = 8)
layout(
  matrix(1:prod(dim(h)[1:2] + 1), ncol = dim(h)[1] + 1, byrow = TRUE),
  height = c(rev(apply(h, 2, max)) + m, 0.15),
  width = c(5, rep(1, dim(h)[1]))
)
par(mar = c(0.05, 0.1, 0, 0.1), xpd = NA)
for (j in dim(h)[2]:0 + 1)
  for (i in 0:dim(h)[1] + 1) {
    #if(all(h[i,j,]==0))
    if (i == 1 & j != 1) {
      plot(
        NA,
        NA,
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n",
        xlim = c(0, 1),
        ylim = c(0, 1),
        bty = "n"
      )
      text(1, 0, as.character(Subtype_lift[dimnames(h)[[2]]])[j - 1], pos = 2)
      next
    }
    if (j == 1) {
      plot(
        NA,
        NA,
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n",
        xlim = c(0, 1),
        ylim = c(0, 1),
        bty = "n"
      )
      if (i == 1) next
      text(0.5, 1, dimnames(h)[[1]][i - 1], pos = 1)
      next
    }
    r <- c(0, max(h[, j - 1, ] + m))
    par(bty = if (i == 2) "L" else "n")
    barplot(
      h[i - 1, j - 1, ],
      ylim = r,
      width = 1 / n,
      space = 0,
      col = rev(col),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      border = NA,
      xpd = TRUE,
      yaxs = "i",
      xaxs = "i",
      xlim = c(-0.5 / n, 1 + 0.5 / n)
    )
    axis(side = 1, at = c(-0.5 / n, 1 + 0.5 / n), labels = c("", ""), tcl = -.1)
    if (i > 1) abline(v = 0, col = 'lightgrey', lty = 3)
    if (i == 2) {
      abline(h = 0.05 * 8 / n, col = 'lightgrey', lty = 1)
      axis(side = 2, at = c(0, 0.05 * 8 / n), labels = c("", ""), tcl = -.1)
    }
  }
dev.off()


# Synchronous gains -------------------------------------------------------
colTime <- c("#A0C758", "#6B8934", "#BEC6AD", "#CEB299", "#CC6415", "#EF7B00")
names(colTime) <- levels(timingClass)[c(4, 5, 6, 3, 2, 1)]
## timing class by group
timingclass_groups <- NULL
for (igroup in unique(wgs_data$Subtype)) {
  tmplist <- wgs_data %>% filter(Subtype == igroup) %>% pull(Tumor_Barcode)
  d <- fracGenomeWgdComp[rownames(fracGenomeWgdComp) %in% tmplist, ]
  i <- d[, "avg.ci"] <= 0.5 & d[, "chr.all"] > 2 #&  fracGenomeWgdComp[,"nt.total"]/chrOffset["M"] >= 0.1
  timingClass <- paste(
    ifelse(isWgd, "WGD", "ND"),
    ifelse(!i, "uninformative", "")
  )
  timingClass[i] <- paste0(
    timingClass[i],
    ifelse(d[i, "nt.wgd"] / d[i, "nt.total"] > 0.75, "sync", "async")
  )
  #timingClass[i] <- paste0(timingClass[i], cut(fracGenomeWgdComp[i,"nt.wgd"]/fracGenomeWgdComp[i,"nt.total"], c(0,0.5,0.8,1), include.lowest=TRUE))
  tmp <- tibble(Type = timingClass) %>% count(Type) %>% mutate(Subtype = igroup)
  timingclass_groups <- bind_rows(tmp, timingclass_groups)
}


timingclass_groups %>%
  ggplot(aes(Subtype, n, fill = Type)) +
  geom_bar(position = "fill", stat = "identity", width = 0.7) +
  labs(x = "", y = "Proportion of timing class") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_fill_manual(values = colTime, na.value = '#cccccc') +
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    axis = 'Y',
    axis_col = 'black'
  ) +
  theme(
    panel.spacing = unit(0.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text.x = element_text(hjust = 0.5),
    plot.margin = margin(4, 4, 4, 4)
  ) +
  panel_border(color = 'black')

ggsave('TimingClass_group.pdf', width = 4, height = 6, device = cairo_pdf)


d <- fracGenomeWgdComp
i <- d[, "avg.ci"] <= 0.5 & d[, "chr.all"] > 2 #&  fracGenomeWgdComp[,"nt.total"]/chrOffset["M"] >= 0.1
timingClass <- paste(
  ifelse(isWgd, "WGD", "ND"),
  ifelse(!i, "uninformative", "")
)
timingClass[i] <- paste0(
  timingClass[i],
  ifelse(d[i, "nt.wgd"] / d[i, "nt.total"] > 0.75, "sync", "async")
)
#timingClass[i] <- paste0(timingClass[i], cut(fracGenomeWgdComp[i,"nt.wgd"]/fracGenomeWgdComp[i,"nt.total"], c(0,0.5,0.8,1), include.lowest=TRUE))
timingClass <- factor(timingClass)

#timingClass <- factor(timingClass,levels = c(levels(timingClass),"WGD uninformative"))

## repeat this for subtypes

pdf("TimingClass.pdf", 8, 8)
colTime <- c("#A0C758", "#6B8934", "#BEC6AD", "#CEB299", "#CC6415", "#EF7B00")
names(colTime) <- levels(timingClass)[c(4, 5, 6, 3, 2, 1)]
c <- c(RColorBrewer::brewer.pal(9, "Pastel1"), "#DDDDDD")
t <- table(timingClass)[names(colTime)]
t[is.na(t)] <- 0
names(t[t == 0]) <- "WGD uninformative"
pie(t, init.angle = 90, labels = paste0(names(t), ",\nn=", t), col = colTime)
#t <- table(isWgd)
par(new = TRUE)
symbols(x = 0, y = 0, circles = 0.4, inches = FALSE, add = TRUE, bg = "white")
#pie(t, labels=c("",""), col=NA, lwd=5, lty=1, init.angle=90)
dev.off()

#Figure1f <- createSheet(Figure1, "Figure1f")
#addDataFrame(t, Figure1f)

colnames(d) <- c(
  "ntCoamp",
  "ntAmp",
  "timeCoamp",
  "segCoamp",
  "segAmp",
  "chrCoamp",
  "chrAmp",
  "sdTimeCoamp",
  "avgCiSeg",
  "sdAllSeg"
)
timingInfo <- data.frame(
  avgPloidy = finalPloidy,
  avgHom = finalHom,
  isWgd = isWgd,
  d,
  informative = i,
  timingClass = timingClass
)
#tab <- rbind(tab, data.frame(WGD_call=otherStat, WGD_timing=NA, ploidy=otherPloidy, hom=otherHom, nt.wgd=NA, nt.total=NA, time.wgd=NA, sd.wgd=NA,avg.ci=NA, sd.all=NA))
write.table(
  file = paste0(Sys.Date(), "-Timing-info.txt"),
  timingInfo,
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE,
  sep = "\t"
)

# figure 1g
names(timingClass) <- names(finalBB)

cn <- do.call(
  "rbind",
  sapply(
    names(finalBB),
    function(n) {
      bb <- finalBB[[n]]
      data.frame(
        sample = n,
        chr = seqnames(bb),
        width = width(bb),
        M = bb$major_cn,
        m = bb$minor_cn
      )
    },
    simplify = FALSE
  )
)

t <- table(
  pmin(cn$M, 3),
  pmax(3, round(log10(cn$width), 1)),
  timingClass[cn$sample]
)
x <- as.numeric(colnames(t))
plot(
  NA,
  NA,
  type = 'p',
  col = colTime[1],
  pch = 16,
  ylim = c(0, 0.8),
  xlim = range(10^x),
  xlab = "Length of segment",
  ylab = "Proportion >2 allelic copies",
  log = 'x'
)
for (n in dimnames(t)[[3]]) {
  y <- as.numeric(t[4, , n] / colSums(t[3:4, , n]))
  lines(
    10^x[!is.na(y)],
    y[!is.na(y)],
    type = 'p',
    col = paste0(colTime[n], "44"),
    pch = 16,
    cex = 1
  ) #sqrt(colSums(t[3:4,,i]/1000)))
  lines(
    10^x[!is.na(y)],
    predict(loess(y[!is.na(y)] ~ x[!is.na(y)])),
    col = colTime[n],
    lwd = 2
  )
}


pdf('relative_latency.pdf', width = 4, height = 4)
tt <- mg14:::asum(t[, x >= 7, ], 2)
o <- names(colTime)[-3]
p <- tt[4, o] / colSums(tt[3:4, o])
ci <- sapply(
  c(0.025, 0.975),
  qbeta,
  0.025,
  shape1 = tt[4, o] + 1,
  shape2 = tt[3, o] + 1
)
barplot(
  p,
  col = colTime[-3],
  border = NA,
  las = 2,
  ylab = "Proportion >2 allelic copies",
  names = sub("ormative", "", sub("near-diploid", "ND", names(colTime)[-3])),
  ylim = c(0, 0.4)
) -> b
segments(b, ci[, 1], b, ci[, 2])
dev.off()


# Chronological WGD & MRCA ------------------------------------------------
load('wgs_clinical.RData')
tgct_mtvcf0 <- tgct_mtvcf
load('DP_info_data.RData')

tgct_mtvcf0 <- tgct_mtvcf0 %>%
  left_join(
    DP_info_data %>% select(Tumor_Barcode, MutationID = ID, context, mutType)
  ) %>%
  select(-pAllSubclones)

tgct_mtvcf0 %>% count(is.na(mutType) | is.na(CLS))
tgct_mtvcf0 <- tgct_mtvcf0 %>% filter(!is.na(mutType), !is.na(CLS))

finalSnv <- tibble(Tumor_Barcode = names(finalBB)) %>%
  left_join(tgct_mtvcf0)

finalSnv <- split(finalSnv, finalSnv$Tumor_Barcode)


clusters <- Cluster_info %>%
  mutate(Cluster_CCF = if_else(Cluster_CCF > 1.2, 1.2, Cluster_CCF)) %>%
  dplyr::select(
    Tumor_Barcode,
    cluster = Cluster,
    n_ssms = Cluster_nMut,
    proportion = Cluster_CCF
  ) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  mutate(proportion = proportion * BB_Purity)

finalClusters <- tibble(Tumor_Barcode = names(finalBB)) %>%
  left_join(clusters)

finalClusters <- split(finalClusters, finalClusters$Tumor_Barcode)

unique(names(finalBB) == names(finalClusters))


# CpG>TpG branch lengths and mutation rates -------------------------------
tmp <- wgs_clinical %>%
  select(Subject, age_at_diagnosis) %>%
  left_join(wgs_data)
age <- tmp$age_at_diagnosis
names(age) <- tmp$Tumor_Barcode

# Effective (time-averaged) genome size
# tmp = tgct_mtvcf0 %>%
#   filter(CLS!="subclonal",str_detect(mutType,'\\[C>T\\]G')) %>%
#   mutate(cnWeights = MutCN/(MajCN + MinCN)) %>%
#   group_by(Tumor_Barcode) %>%
#   summarise(effGenome = 2/mean(2*cnWeights, na.rm=TRUE))
#
# effGenome <-  tmp$effGenome
# names(effGenome) <-  tmp$Tumor_Barcode
#effGenome_old <- effGenome

cnWeights <- function(vcf) {
  t <- vcf$MajCN + vcf$MinCN
  vcf$MutCN / t
}

avgWeights <- function(vcf, type = c("all", "deamination")) {
  type <- match.arg(type)
  if (type == "deamination") w <- str_detect(vcf$mutType, '\\[C>T\\]G') else
    w <- rep(TRUE, nrow(vcf))
  cnw <- cnWeights(vcf)
  mean(2 * cnw[w], na.rm = TRUE)
}

effGenome <- unlist(mclapply(finalSnv, function(vcf) {
  w <- vcf$CLS != "subclonal"
  w <- w & str_detect(vcf$mutType, '\\[C>T\\]G')
  2 / avgWeights(vcf[na.omit(w), ])
}))
names(effGenome) <- names(finalSnv)


# Power per (sub)clone
# finalPower <- sapply(names(finalBB), function(n) {
#   x <- finalBB[[n]]
#   f <- finalClusters[[n]]$proportion
#   ##Rounding to 1 significant figure
#   #f <- round(f,1)
#   tmp <- tibble(f=f) %>% mutate(fid=seq_along(f)) %>% arrange(f) %>% mutate(oid=seq_along(f))
#
#   for(i in 1:length(x)){
#     t <- x$timing_param[[i]]
#     if(length(unique(t[,"cfi"])) != length(f)) next
#     suppressMessages(
#       tid <- tibble(t=t[,"cfi"]) %>% unique() %>% mutate(tid=seq_along(t)) %>% arrange(t) %>% mutate(oid=seq_along(t)) %>%
#         full_join(tmp) %>%
#         arrange(fid) %>%
#         pull(tid)
#     )
#     #p <- t[match(f, round(t[,"cfi"],1)), "power.s"] # f ordered
#     p <- t[tid, "power.s"] # f ordered
#     if(!is.null(p)) if(all(!is.na(p))) break
#   }
#   if(is.null(p)) return(rep(NA, length(f)))
#   return(p)
# })

finalPower <- sapply(names(finalBB), function(n) {
  x <- finalBB[[n]]
  f <- finalClusters[[n]]$proportion
  ##Rounding to 1 significant figure
  #f <- round(f,1)
  tmp <- tibble(f = f) %>%
    mutate(fid = seq_along(f)) %>%
    arrange(f) %>%
    mutate(oid = seq_along(f)) %>%
    mutate(round = round(f, 1))
  tmp2 <- as_tibble(do.call(
    "rbind",
    tgct_mtbb0[tgct_mtbb0$Tumor_Barcode == n, ]$timing_param
  )) %>%
    select(cfi, power.s) %>%
    unique() %>%
    arrange(cfi) %>%
    mutate(round = round(cfi, 1)) %>%
    filter(!is.na(power.s))
  suppressMessages(
    expr = (p <- tmp %>%
      left_join(tmp2) %>%
      group_by(fid) %>%
      arrange(is.na(power.s), abs(cfi - f)) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      arrange(fid) %>%
      pull(power.s))
  )

  if (length(p[is.na(p)]) == 1) {
    p[is.na(p)] <- tmp2$power.s[which.min(abs(tmp2$cfi - f[is.na(p)]))]
  }

  if (is.null(p)) return(rep(NA, length(f)))
  return(p)
})


# # manual adjusted
# n="NSLC-0993-T01"
# finalPower[[n]] <- c(0.9995582,,0.9597421)
# "NSLC-1030-T01"
# "NSLC-0327-T01"
# #tmp <- finalPower[["NSLC-0778-T01"]]
# #finalPower[["NSLC-0778-T01"]] <- c(tmp[1],NA_real_,tmp[2:3])

# Branch lengths
branchDeam <- t(simplify2array(mclapply(finalSnv, function(vcf) {
  #n <- meta(header(vcf))$META["ID",]
  n <- vcf$Tumor_Barcode[1]
  w <- str_detect(vcf$mutType, '\\[C>T\\]G')
  if (sum(w) == 0) return(c(0, 0))
  p <- vcf$pSub[w]
  n.subclonal <- aggregate(p, list(vcf$CNF[w]), sum, na.rm = TRUE)
  m <- apply(
    abs(outer(n.subclonal$Group.1, finalClusters[[n]]$proportion, `-`)),
    1,
    which.min
  ) # Match to subclones
  p.subclonal <- finalPower[[n]][m] # Power of subclones
  #tm <- (n.subclonal$Group.1 / p.subclonal) / max(n.subclonal$Group.1)
  #tm[is.na(tm)] <- 0
  #b.subclonal<- n.subclonal$x %*% tm
  b.subclonal <- n.subclonal$x %*%
    (n.subclonal$Group.1 / p.subclonal) /
    max(n.subclonal$Group.1) # Subclonal branch, power adjusted & 1/f-scaled

  b.clonal <- sum(1 - p, na.rm = TRUE) / finalPower[[n]][1] # Clonal branch (trunk), power adjusted & 1/f-scaled
  c(b.subclonal, b.clonal)
})))

branchDeam[is.na(branchDeam)]
branchDeam[is.na(branchDeam[, 1]) | branchDeam[, 1] == 0, ]
branchDeam[is.na(branchDeam[, 2]) | branchDeam[, 2] == 0, ]

#d <- droplevels(donor2type[sample2donor[names(finalSnv)]])
#typesSubclones <- setdiff(levels(d), c(typeNa, names(which(table(d)<5))))

nClones <- sapply(finalClusters, nrow)

branchDeamLinear <- t(simplify2array(mclapply(finalSnv, function(vcf) {
  w <- str_detect(vcf$mutType, '\\[C>T\\]G')
  if (sum(w) == 0) return(c(0, 0))
  n <- vcf$Tumor_Barcode[1]
  p <- vcf$pSub[w]
  n.subclonal <- aggregate(p, list(vcf$CNF[w]), sum, na.rm = TRUE)
  m <- apply(
    abs(outer(n.subclonal$Group.1, finalClusters[[n]]$proportion, `-`)),
    1,
    which.min
  ) # Match to subclones
  p.subclonal <- finalPower[[n]][m] # Power of subclones
  b.subclonal <- n.subclonal$x %*% (1 / p.subclonal) # Subclonal branch, power adjusted
  b.clonal <- sum(1 - p, na.rm = TRUE) / finalPower[[n]][1] # Clonal branch (trunk), power adjusted
  c(b.subclonal, b.clonal)
})))


branchDeamLinear[is.na(branchDeamLinear)]
branchDeamLinear[is.na(branchDeamLinear[, 1]) | branchDeamLinear[, 1] == 0, ]
branchDeamLinear[is.na(branchDeamLinear[, 2]) | branchDeamLinear[, 2] == 0, ]


#Compare branching and linear topologies:
f <- (branchDeam[, 1] / finalPloidy) /
  rowSums(branchDeam / cbind(finalPloidy, effGenome))
l <- (branchDeamLinear[, 1] / finalPloidy) /
  rowSums(branchDeamLinear / cbind(finalPloidy, effGenome))
# t <- donor2type[sample2donor[names(finalSnv)]]
#plot(f, l, xlab="Subclonal branch length (branching)", ylab="Subclonal branch length (linear)", pch=21, bg=tissueColors[t], col=tissueBorder[t], cex=tissueCex[t])
plot(
  f,
  l,
  xlab = "Subclonal branch length (branching)",
  ylab = "Subclonal branch length (linear)",
  pch = 21
)
abline(0, 1, lty = 2)

quantile(l / f, na.rm = TRUE)
quantile(l, na.rm = TRUE)
quantile(f, na.rm = TRUE)

unique(names(f) == names(l))

# Compare branching and linear topologies:
tibble(Tumor_Barcode = names(f), f = f, l = l) %>%
  left_join(wgs_data) %>%
  ggplot(aes(f, l, fill = Subtype)) +
  geom_point(pch = 21, size = 3, col = "white") +
  geom_abline(slope = 1, linetype = 2) +
  scale_fill_manual(values = subtypecol, breaks = names(subtypecol)[1:3]) +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_text_size = 14,
    axis_title_size = 16
  ) +
  panel_border(color = 'black') +
  labs(
    x = "Subclonal branch length (branching)",
    y = "Subclonal branch length (linear)"
  )

ggsave(
  filename = 'Compare_branching_linear.pdf',
  width = 6,
  height = 5,
  device = cairo_pdf
)


# Average mutation rates

# Average mutation rates --------------------------------------------------

rateDeam <- cc <- list()
#remove <- "8454fe53-869d-41c8-b0c8-a7929d00eec3" # a cell line, add more samples in the following
remove <- NULL
#par(mfrow=c(6,6), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1, xpd=FALSE)

typesSubclones <- unique(wgs_data$Subtype)
d <- tibble(Tumor_Barcode = names(finalBB)) %>%
  left_join(wgs_data) %>%
  select(Tumor_Barcode, Subtype)

for (n in typesSubclones) {
  i <- d$Tumor_Barcode[d$Subtype == n]
  tt0 <- branchDeam[i, ] / cbind(finalPloidy[i], effGenome[i]) / 3 #cbind(nClones[i]-1, 1)/3 # 3Gb Haploid genome
  tt0[is.infinite(tt0) | is.nan(tt0)] <- 0
  yy <- rowSums(tt0)
  a <- age[i]
  xx <- a
  r <- yy / xx
  # m <- median(r[TiN[names(xx)] <= 0.01 & ! is.na(TiN[names(xx)])],na.rm=TRUE)
  m <- median(r, na.rm = T)
  rateDeam[[n]] <- r
  try({
    #w <- (r-m)^2/m^2 <= 2^2 & TiN[names(xx)] <= 0.01 & ! is.na(TiN[names(xx)])
    w <- (r - m)^2 / m^2 <= 2^2
    remove <- c(remove, names(which(!w)))
    plot(
      xx,
      yy,
      pch = NA,
      log = '',
      xlab = "Age at diagnosis",
      ylab = "SNVs/Gb",
      main = n,
      ylim = c(0, pmin(1000, max(yy, na.rm = TRUE))),
      xlim = c(0, max(age, na.rm = TRUE)),
      cex.main = 1
    )
    #bg=tissueColors[n], col=tissueBorder[n],
    #par(xpd=NA)
    #segments(x0=0,y0=0, xx, yy, col=tissueLines[n], lty=tissueLty[n])
    points(
      xx,
      yy,
      pch = ifelse(w, 21, 4),
      bg = as.character(subtypecol[n]),
      cex = 1.5
    ) # bg=tissueColors[n], col=ifelse(w,tissueBorder[n], tissueColors[n])
    abline(0, m, lty = 3)
    #lines(c(x0,2*x0), c(0,1))
    #print(paste(n,cor(xx[w],yy[w],use='c'), cor(xx[w],tt0[w,] %*% c(0.2,1), use='c'), sep=": "))
    l <- lm(
      yy ~ xx,
      data = data.frame(
        yy = yy[w],
        xx = xx[w],
        row.names = as.character(seq_along(yy[w]))
      )
    )
    f <- (summary(l))
    p <- predict(l, newdata = data.frame(xx = seq(0, 100, 1)), se = TRUE)
    polygon(
      c(seq(0, 100, 1), rev(seq(0, 100, 1))),
      c(p$fit - 2 * p$se.fit, rev(p$fit + 2 * p$se.fit)),
      border = 'black',
      col = alpha(as.character(subtypecol[n]), 0.2)
    ) #border=tissueBorder[n], col=paste0(tissueColors[n],"44")
    v <- which(!is.na(xx[w] * yy[w]))
    cc[[n]] <- cbind(
      f$coefficients,
      f$cov.unscaled * f$sigma^2,
      coef(nnls::nnls(cbind(1, xx[w][v]), yy[w][v]))
    )
  })
}


n <- names(rateDeam)
qRateDeam <- sapply(rateDeam, function(r) {
  m <- median(r, na.rm = TRUE)
  w <- (r - m)^2 / m^2 <= 2^2
  quantile(r[w], na.rm = TRUE)
})
plot(
  sapply(rateDeam, median, na.rm = TRUE),
  pch = NA,
  ylab = "SNVs/Gb/yr",
  main = "CpG>TpG rate",
  ylim = c(0, max(qRateDeam)),
  cex.main = 1,
  xaxt = 'n',
  xlab = "Tumour type"
)
segments(
  seq_along(rateDeam),
  qRateDeam["0%", ],
  seq_along(rateDeam),
  qRateDeam["100%", ],
  col = as.character(subtypecol[n]),
  lty = 1
)
points(
  sapply(rateDeam, median, na.rm = TRUE),
  pch = 21,
  col = '#cccccc',
  bg = as.character(subtypecol[n])
)

length(remove)


#Plot the average rates as barplot across tissues. Relatively little variation.
par(mar = c(6, 3, 1, 1))
o <- order(qRateDeam["50%", ])
barplot(
  qRateDeam["50%", ][o],
  border = "#cccccc",
  col = as.character(subtypecol[colnames(qRateDeam)])[o],
  las = 2,
  names.arg = rep("", ncol(qRateDeam)),
  ylab = "CpG>TpG rate [SNVs/Gb/yr]",
  ylim = c(0, max(qRateDeam))
) -> b
mg14::rotatedLabel(b, labels = colnames(qRateDeam)[o])
segments(
  b,
  qRateDeam["50%", ][o],
  b,
  qRateDeam["100%", ][o],
  col = as.character(subtypecol[colnames(qRateDeam)])[o],
  lwd = 2
)
segments(
  b,
  qRateDeam["0%", ][o],
  b,
  qRateDeam["50%", ][o],
  col = as.character(subtypecol[colnames(qRateDeam)])[o],
  lwd = 2
)


#Plot the copynumber and ITH-adjusted mutation burden versus age and add a loess fit.

tt0 <- branchDeam / cbind(finalPloidy, effGenome) / cbind(nClones - 1, 1) / 3 # 3Gb Haploid genome
tt0[is.infinite(tt0) | is.nan(tt0)] <- 0
m <- sapply(rateDeam, function(r) {
  m <- median(r, na.rm = TRUE)
})
s <- rowSums(tt0) #/m[as.character(donor2type[sample2donor[names(finalSnv)]])]
s[remove] <- NA
t <- barcode2spg[names(finalSnv)]
x <- age[names(finalSnv)]
plot(
  x,
  s,
  bg = as.character(subtypecol[t]),
  pch = 21,
  ylim = c(0, 1000),
  col = "#cccccc",
  cex = 1.5,
  lwd = 0.25,
  xlab = "Age",
  ylab = "SNVs/Gb"
)
p <- predict(loess(s ~ x), newdata = sort(x, na.last = NA), se = TRUE)
r <- function(x) c(x, rev(x))
polygon(
  r(sort(x, na.last = NA)),
  c(p$fit + 2 * p$se, rev(p$fit - 2 * p$se)),
  col = "#00000044",
  border = NA
)
lines(sort(x, na.last = NA), p$fit)

# ExtendedDataFigure8 <- xlsx::createWorkbook()
# ExtendedDataFigure8a <- xlsx::createSheet(ExtendedDataFigure8, "ExtendedDataFigure8a")
mrate <- data.frame(
  Tumor_Barcode = names(x),
  age = x,
  `CpG>TpG/Gb` = s,
  Subtype = t
) %>%
  as_tibble()
# xlsx::addDataFrame(tab, ExtendedDataFigure8a)
a <- simplify2array(cc)
all(a[1, 1, ] + 2 * a[1, 2, ] > 0)

plot(
  a[1, 1, ],
  a[2, 1, ],
  col = as.character(subtypecol[dimnames(a)[[3]]]),
  pch = NA,
  xlab = "Offset",
  ylab = "SNVs/Gb/yr"
)
segments(
  a[1, 1, ],
  a[2, 1, ] - a[2, 2, ],
  a[1, 1, ],
  a[2, 1, ] + a[2, 2, ],
  col = as.character(subtypecol[dimnames(a)[[3]]]),
  pch = 19
)
segments(
  a[1, 1, ] - a[1, 2, ],
  a[2, 1, ],
  a[1, 1, ] + a[1, 2, ],
  a[2, 1, ],
  col = as.character(subtypecol[dimnames(a)[[3]]]),
  pch = 19
)
points(
  a[1, 1, ],
  a[2, 1, ],
  pch = 21,
  bg = as.character(subtypecol[dimnames(a)[[3]]]),
  col = as.character(subtypecol[dimnames(a)[[3]]])
)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)


#The two quantities above can be used to calculate a fraction of mutations due to linear accumulation
# par(mar=c(6,3,1,1))
ma <- sapply(split(age, barcode2spg[names(age)]), median, na.rm = TRUE)
fm <- pmax(a[2, 7, ], 0) *
  ma[dimnames(a)[[3]]] /
  (pmax(0, a[2, 7, ]) * ma[dimnames(a)[[3]]] + pmax(0, a[1, 7, ])) *
  100
o <- order(fm)
fmq <- sapply(names(fm), function(n) {
  aa <- mvtnorm::rmvnorm(10000, mean = a[, 1, n], sigma = a[, 5:6, n])
  aa <- aa[aa[, 1] >= 0 & aa[, 2] >= 0, ]
  quantile(
    pmax(aa[, 2], 0) * ma[n] / (pmax(0, aa[, 2]) * ma[n] + pmax(0, aa[, 1])),
    c(0.025, 0.975)
  )
}) *
  100
barplot(
  fm[o],
  col = as.character(subtypecol[dimnames(a)[[3]]][o]),
  border = as.character(subtypecol[dimnames(a)[[3]]][o]),
  las = 2,
  names.arg = rep("", length(fm)),
  ylab = "Age-attributed mutations [%]"
) -> b
mg14::rotatedLabel(b, labels = names(fm[o]))
segments(
  b,
  fm[o],
  b,
  fmq[2, o],
  col = as.character(subtypecol[dimnames(a)[[3]]][o]),
  lwd = 2
)
segments(
  b,
  fmq[1, o],
  b,
  fm[o],
  col = as.character(subtypecol[dimnames(a)[[3]]][o]),
  lwd = 2
)
abline(h = min(fmq[2, ]))
abline(h = max(fmq[1, ]))


# Hierarchical Bayesian models of deamination rates -----------------------
#
# # prepre data for rstan
# library(rstan)
# y <- Reduce("c",rateDeam)
# y <- y[!names(y) %in% remove] #Remove hypermutated samples
# x <- age[names(y)]
# rownames(classification) <- classification$Tumor_Barcode
# rownames(histology) <- histology$Tumor_Barcode
# tmp <-  tibble(Tumor_Barcode=names(y)) %>% left_join(wgs_data) %>% left_join(clinical_data)
#
# t <- tmp$Subtype
# names(t) <- tmp$Tumor_Barcode
# z <-tmp$Histology
# names(z) <- tmp$Tumor_Barcode
# y <- y*x
# df <- data.frame(x,y,t,z)
# df <- df[rowSums(is.na(df))==0,]
# df <- df[df$y > 0,]
#
# #Visualizing the data. The rate is normalized for ploidy, tumor purity, subclonal architecture. ## Need to split piano in this plot further to carcinoid and LUAD
# ggplot(data = df, aes(x=x, y=y, color=t)) +
#   geom_point(shape=16) +
#   xlab("Age (years)") + ylab("CpG>TpG SNVs/Gb") + xlim(c(0,90)) + labs(color="Subtypes")+ geom_smooth(method = 'lm')+scale_color_d3()
#
#
# df$t <- as.factor(df$t)
# tt <- model.matrix(~droplevels(t)-1, data=df)
#
# data <- list(n = nrow(df),
#              p = ncol(tt),
#              y = df$y,
#              x = tt * df$x,
#              t = tt
# )
#

# MRCA --------------------------------------------------------------------
accel <- c(1, 2.5, 5, 7.5, 10, 20)
names(accel) <- paste0(accel, "x")

set.seed(42)

d <- barcode2spg[names(finalSnv)]

computeSubclonesTimeAbs <- function(l, b) {
  i <- d == l
  tt0 <- b[i, ] / cbind(finalPloidy[i], effGenome[i]) #/ cbind(nClones[i]-1, 1)
  resB <- sapply(
    1:1000,
    function(foo) {
      ## Assess the impact of Poisson fluctuations on numbers
      tt <- matrix(rpois(length(tt0), lambda = tt0), ncol = ncol(tt0))
      res <- sapply(
        accel,
        function(a) tt[, 1] / a / rowSums(tt / rep(c(a, 1), each = nrow(tt)))
      ) *
        age[names(finalSnv)[i]]
      colnames(res) <- paste0(accel, "x")
      #res[res==0] <- NA
      res
    },
    simplify = 'array'
  )
  res <- sapply(
    accel,
    function(a) tt0[, 1] / a / rowSums(tt0 / rep(c(a, 1), each = nrow(tt0)))
  ) *
    age[names(finalSnv)[i]]
  colnames(res) <- paste0(accel, "x")
  resCI <- apply(resB, 1:2, quantile, c(0.1, 0.9), na.rm = TRUE)
  arr <- abind::abind(res, resCI, along = 1)
  rownames(arr)[1] <- "hat"
  arr <- aperm(arr, c(2, 1, 3))
  tt0[is.infinite(tt0) | is.nan(tt0)] <- 0
  r <- which(rowSums(b[i, ]) < 50) ## Exclude samples with less than 50 subs
  arr[r, , ] <- NA
  return(arr)
}

# Timing
subclonesTimeAbs <- mclapply(
  typesSubclones,
  computeSubclonesTimeAbs,
  b = branchDeam
)
subclonesTimeAbsLinear <- mclapply(
  typesSubclones,
  computeSubclonesTimeAbs,
  b = branchDeamLinear
)
names(subclonesTimeAbsLinear) <- names(subclonesTimeAbs) <- typesSubclones

#Buest guess acceleration values (usually 5x, except Ovarian and CNS).
guessAccel <- sapply(subclonesTimeAbs, function(x) "5x")
#guessAccel["Ovary-AdenoCa"] <- guessAccel["Liver-HCC"] <- "7.5x"
#guessAccel[grep('CNS', names(guessAccel))] <- "2.5x

#Plot the inferred latencies per sample with overlaid boxplots.
#u <- setdiff(names(finalSnv)[uniqueSamples], remove)
u <- setdiff(names(finalSnv), remove)
#par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
qSubclone <- sapply(
  subclonesTimeAbs,
  function(x)
    apply(
      x[, "hat", ][rownames(x) %in% u, , drop = FALSE],
      2,
      quantile,
      c(0.05, 0.25, 0.5, 0.75, 0.95),
      na.rm = TRUE
    ),
  simplify = 'array'
)
a <- "5x"
subclonesTimeAbsType <- sapply(names(subclonesTimeAbs), function(n) {
  x <- subclonesTimeAbs[[n]]
  x[,, guessAccel[n]][setdiff(rownames(x), remove), 1:3, drop = FALSE]
})
m <- diag(qSubclone["50%", guessAccel[dimnames(qSubclone)[[3]]], ]) #t[1,3,]
names(m) <- dimnames(qSubclone)[[3]]
nSubclones <- sapply(subclonesTimeAbsType, function(x) sum(!is.na(x[, 1])))
m[nSubclones < 5] <- NA
o <- order(m, na.last = NA)
plot(
  NA,
  NA,
  xlim = c(0.5, length(m[o])),
  ylab = "Years before diagnosis",
  xlab = "",
  xaxt = "n",
  yaxs = "i",
  ylim = c(0, 30)
)
abline(h = seq(10, 20, 10), col = "#DDDDDD", lty = 3)
x <- seq_along(m[o])
mg14::rotatedLabel(x, labels = names(sort(m)))
for (i in seq_along(o))
  try({
    n <- names(m)[o[i]]
    f <- function(x) x / max(abs(x))
    a <- guessAccel[n]
    bwd <- 0.8 / 2
    j <- if (length(na.omit(subclonesTimeAbsType[[o[i]]][, "hat"])) > 1)
      f(mg14::violinJitter(na.omit(subclonesTimeAbsType[[o[i]]][, "hat"]))$y) /
        4 +
        i else i
    tpy <- 2
    segments(
      j,
      na.omit(subclonesTimeAbsType[[o[i]]][, "90%"]),
      j,
      na.omit(subclonesTimeAbsType[[o[i]]][, "10%"]),
      col = mg14::colTrans(as.character(subtypecol[n]), tpy)
    )
    points(
      j,
      na.omit(subclonesTimeAbsType[[o[i]]][, "hat"]),
      pch = 21,
      col = mg14::colTrans(as.character(subtypecol[n]), tpy),
      bg = mg14::colTrans(as.character(subtypecol[n]), tpy),
      cex = 1.5,
      lwd = 1
    )
    rect(
      i - bwd,
      qSubclone["25%", a, n],
      i + bwd,
      qSubclone["75%", a, n],
      border = as.character(subtypecol[n]),
      col = as.character(subtypecol[n])
    )
    segments(
      i - bwd,
      qSubclone["50%", a, n],
      i + bwd,
      qSubclone["50%", a, n],
      col = as.character(subtypecol[n]),
      lwd = 2
    )
    segments(
      i,
      qSubclone["75%", a, n],
      i,
      qSubclone["95%", a, n],
      col = as.character(subtypecol[n]),
      lwd = 1.5
    )
    segments(
      i,
      qSubclone["5%", a, n],
      i,
      qSubclone["25%", a, n],
      col = as.character(subtypecol[n]),
      lwd = 1.5
    )
    # rect(i-bwd,qSubclone["25%",a,n],i+bwd,qSubclone["75%",a,n], border=tissueLines[n],  col=paste(tissueColors[n],"44", sep=""))
    # segments(i-bwd,qSubclone["50%",a,n],i+bwd,qSubclone["50%",a,n],col=tissueLines[n], lwd=2)
    # segments(i,qSubclone["75%",a,n],i,qSubclone["95%",a,n],col=tissueLines[n], lwd=1.5)
    # segments(i,qSubclone["5%",a,n],i,qSubclone["25%",a,n],col=tissueLines[n], lwd=1.5)
    f <- function(x) x / max(abs(x))
  })


#Alternative representation of median acceleration across the range of simulated rate increases, per cancer type:
u <- setdiff(names(finalSnv), remove)
m <- qSubclone["50%", "5x", ] #t[1,3,]
names(m) <- dimnames(qSubclone)[[3]]
m[nSubclones < 5] <- NA
o <- order(m, na.last = NA)
n <- dimnames(qSubclone)[[3]]
par(
  mar = c(7, 3, 1, 1),
  mgp = c(2, .5, 0),
  tcl = 0.25,
  cex = 1,
  bty = "L",
  xpd = FALSE,
  las = 1
)
plot(
  NA,
  NA,
  xlim = c(0.5, length(m[o]) + 0.5),
  ylab = "Latency [yr]",
  xlab = "",
  xaxt = "n",
  yaxs = "i",
  ylim = c(0, 18)
)
abline(h = seq(10, 20, 10), col = "#DDDDDD", lty = 3)
x <- seq_along(m[o])
mg14::rotatedLabel(x, labels = names(sort(m)))
b <- .3
rect(
  seq_along(o) - b,
  qSubclone["50%", "1x", o],
  seq_along(o) + b,
  qSubclone["50%", "10x", o],
  col = alpha(as.character(subtypecol[n[o]]), 0.4),
  border = 1
)
rect(
  seq_along(o) - b,
  qSubclone["50%", "2.5x", o],
  seq_along(o) + b,
  qSubclone["50%", "7.5x", o],
  col = alpha(as.character(subtypecol[n[o]]), 0.9),
  border = 1
)
rect(
  seq_along(o) - b,
  qSubclone["50%", "20x", o],
  seq_along(o) + b,
  qSubclone["50%", "10x", o],
  col = alpha(as.character(subtypecol[n[o]]), 0.2),
  border = 1
)
segments(
  seq_along(o) - b,
  qSubclone["50%", "5x", o],
  seq_along(o) + b,
  qSubclone["50%", "5x", o]
)


### output the table for all the accesslation rate
MRCAdata <- NULL
for (accel_tmp in names(accel)) {
  guessAccel <- sapply(subclonesTimeAbs, function(x) accel_tmp)
  subclonesTimeAbsType <- sapply(names(subclonesTimeAbs), function(n) {
    x <- subclonesTimeAbs[[n]]
    x[,, guessAccel[n]][setdiff(rownames(x), remove), 1:3, drop = FALSE]
  })
  mrca_tmp <- data.frame(do.call("rbind", subclonesTimeAbsType)) %>%
    rownames_to_column(var = 'Tumor_Barcode') %>%
    as_tibble() %>%
    mutate(acceleration = accel_tmp)
  MRCAdata <- bind_rows(MRCAdata, mrca_tmp)
}


MRCAdata <- MRCAdata %>%
  left_join(wgs_data %>% select(Subject, Tumor_Barcode, Subtype)) %>%
  left_join(wgs_clinical %>% select(Subject, Age = age_at_diagnosis)) %>%
  dplyr::rename(Latency = hat, Latency_X10 = X10., Latency_X90 = X90.) %>%
  mutate(MRCA_age = Age - Latency)


# Save object for ploting -------------------------------------------------
save(
  MRCAdata,
  subclonesTimeAbs,
  subclonesTimeAbsLinear,
  mrate,
  rateDeam,
  timingclass_groups,
  file = 'Chronological_timing_short.RData'
) #timingInfo
save(
  mrate,
  MRCAdata,
  subclonesTimeAbs,
  subclonesTimeAbsLinear,
  finalPloidy,
  finalPower,
  effGenome,
  finalBB,
  finalSnv,
  finalClusters,
  remove,
  branchDeam,
  qRateDeam,
  rateDeam,
  timingclass_groups,
  file = 'Chronological_timing.RData'
) #timingInfo
