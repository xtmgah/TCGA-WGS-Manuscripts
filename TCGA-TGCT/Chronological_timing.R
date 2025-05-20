set_wd()

libztw()
#library(ggnewscale)
pdfhr()
pdfhr2()

# Read Mutation time data for chronological timming -------------------------------------------------
load('wgs_data.RData')
load('Chronological_timing_short.RData')

# Mutational Rate  --------------------------------------------------------
#Plot the copynumber and ITH-adjusted mutation burden versus age and add a loess fit.

myggstyle()
# overall
mrate %>%
  drop_na() %>%
  ggplot(aes(age, log2(CpG.TpG.Gb))) +
  geom_point(pch = 21, size = 2, fill = ncicolpal[14]) +
  geom_smooth(method = 'lm') +
  scale_fill_manual(values = subtypecol) +
  theme_ipsum_rc(axis_title_just = 'm', axis_title_size = 14) +
  labs(x = "Age at diagnosis (years)", y = "SNVs/Gb (log2)") +
  guides(fill = "none") +
  scale_x_continuous(breaks = pretty_breaks()) +
  panel_border(color = 'black', linetype = 1)

ggsave(
  filename = 'mutation_rate_age_all.pdf',
  width = 4,
  height = 3.6,
  device = cairo_pdf
)

# seperated by the group
mrate %>%
  drop_na() %>%
  ggplot(aes(age, log2(CpG.TpG.Gb), fill = Subtype)) +
  geom_point(pch = 21, size = 2) +
  geom_smooth(method = 'lm') +
  facet_wrap(~Subtype, nrow = 1, scales = 'free_x') +
  scale_fill_manual(values = subtypecol) +
  theme_ipsum_rc(axis_title_just = 'm', axis_title_size = 14) +
  labs(x = "Age at diagnosis (years)", y = "SNVs/Gb (log2)") +
  guides(fill = "none") +
  scale_x_continuous(breaks = pretty_breaks()) +
  panel_border(color = 'black', linetype = 1)

ggsave(
  filename = 'mutation_rate_age.pdf',
  width = 12,
  height = 3.5,
  device = cairo_pdf
)

mrate %>%
  do(tidy(lm(log2(CpG.TpG.Gb) ~ age, data = .))) %>%
  filter(term != "(Intercept)")

mrate %>%
  group_by(Subtype) %>%
  do(tidy(lm(log2(CpG.TpG.Gb) ~ age, data = .))) %>%
  filter(term != "(Intercept)") %>%
  ungroup() %>%
  arrange(p.value)
#mrate %>% left_join(clinical_data) %>% filter(Histology == 'Adenocarcinoma') %>% group_by(Subtype) %>% do(tidy(lm(log2(CpG.TpG.Gb)~age,data=.))) %>% filter(term!="(Intercept)") %>% ungroup() %>% arrange(p.value)
#sorder <- read_rds('../SmokingVar/sorder.RDS')
#mrate %>% left_join(sorder) %>% filter(!is.na(APOBEC_Signature)) %>% group_by(APOBEC_Signature) %>% do(tidy(lm(CpG.TpG.Gb~age,data=.))) %>% filter(term!="(Intercept)") %>% ungroup() %>% arrange(p.value)

# average rate as barplot
as.data.frame(qRateDeam) %>%
  rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  mutate(rowname = paste0("d", str_remove(rowname, "%"))) %>%
  pivot_wider(names_from = rowname, values_from = value) %>%
  dplyr::rename(Subtype = name) %>%
  ggplot(aes(Subtype, d50, fill = Subtype)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = subtypecol) +
  guides(fill = "none") +
  geom_errorbar(aes(ymin = d0, ymax = d100), width = 0.1, size = 0.7) +
  labs(y = "CpG>TpG rate [SNVs/Gb/yr]", x = "") +
  theme_ipsum_rc(axis_title_just = 'm', axis_title_size = 14) +
  panel_border(color = 'black', linetype = 1)

ggsave(
  filename = 'mutation_rate_age_average.pdf',
  width = 3,
  height = 4,
  device = cairo_pdf
)


#
# # Latency Plot by differnt accelration rate ------------------------------------------------------------
# u <- base::setdiff(names(finalSnv), remove)
# guessAccel <- sapply(subclonesTimeAbs, function(x) "5x")
# qSubclone <- sapply(subclonesTimeAbs, function(x) apply(x[,"hat",][rownames(x)%in%u,,drop=FALSE], 2, quantile, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE), simplify='array')
# subclonesTimeAbsType <- sapply(names(subclonesTimeAbs), function(n) {x <- subclonesTimeAbs[[n]]; x[,,guessAccel[n]][setdiff(rownames(x),remove), 1:3, drop=FALSE]})
# nSubclones <- sapply(subclonesTimeAbsType, function(x) sum(!is.na(x[,1])))
#
# qSubclone <- qSubclone[,,-2]
#
# m <- qSubclone["50%","5x",]#t[1,3,]
# names(m) <- dimnames(qSubclone)[[3]]
# m[nSubclones < 5] <- NA
# o <- rev(order(m, na.last=NA))
# n <- dimnames(qSubclone)[[3]]
#
# cairo_pdf(filename = 'Acceleration_alterative.pdf',width = 3,height = 5)
# par(mar=c(3,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
# plot(NA,NA, xlim=c(0.5,length(m[o])+0.5), ylab="Latency [yr]", xlab="", xaxt="n", yaxs="i", ylim=c(0,18))
# abline(h=seq(10,20,10), col="#DDDDDD", lty=3)
# x <- seq_along(m[o])
# mg14::rotatedLabel(x, labels=as.character(Subtype_lift[names(rev(sort(m)))]))
# b <- .3
# rect(seq_along(o)-b,qSubclone["50%","1x",o],seq_along(o)+b,qSubclone["50%","10x",o], col=alpha(as.character(Subtype_color[n[o]]),0.4), border=1)
# rect(seq_along(o)-b,qSubclone["50%","2.5x",o],seq_along(o)+b,qSubclone["50%","7.5x",o], col=alpha(as.character(Subtype_color[n[o]]),0.9), border=1)
# rect(seq_along(o)-b,qSubclone["50%","20x",o],seq_along(o)+b,qSubclone["50%","10x",o], col=alpha(as.character(Subtype_color[n[o]]),0.2), border=1)
# segments(seq_along(o)-b,qSubclone["50%","5x",o],seq_along(o)+b,qSubclone["50%","5x",o])
# dev.off()

# Subtype Latency difference----------------------------------------------------------------
library(ggpubr)
my_comparisons <- list(c("Seminoma", "Non-Seminoma"))

MRCAdata %>%
  filter(acceleration == "1x") %>%
  group_by(Subtype) %>%
  summarise(
    Latency = median(Latency, na.rm = T),
    Age = median(Age, na.rm = T),
    MRCA_age = median(MRCA_age, na.rm = T)
  )


MRCAdata %>%
  filter(acceleration == "1x") %>%
  filter(!is.na(Latency)) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  do(tidy(lm(Latency ~ Subtype + BB_Purity + Age, data = .)))

my_comparisons <- list(c('Seminoma', 'Non-Seminoma'))

stat.test <- MRCAdata %>%
  filter(acceleration == "1x") %>%
  filter(!is.na(Latency)) %>%
  rstatix::wilcox_test(Latency ~ Subtype, comparisons = my_comparisons) %>%
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


p1 <- MRCAdata %>%
  filter(acceleration == "1x") %>%
  filter(!is.na(Latency)) %>%
  ggplot(aes(Subtype, Latency, fill = Subtype)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray95",
    color = 'black',
    outlier.shape = 21,
    outlier.size = 0.8
  ) +
  ggbeeswarm::geom_quasirandom(
    pch = 21,
    size = 2.5,
    width = 0.2,
    color = "white",
    stroke = 0.2
  ) +
  labs(x = NULL, y = 'Latency, years before diagnosis') +
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


# ggsave(
#   filename = 'latency_spgroup_1x.pdf',
#   width = 3.5,
#   height = 6,
#   device = cairo_pdf
# )

MRCAdata %>%
  filter(acceleration == "1x") %>%
  filter(!is.na(Latency)) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  do(tidy(lm(MRCA_age ~ Subtype + BB_Purity, data = .)))

my_comparisons <- list(c('Seminoma', 'Non-Seminoma'))

stat.test <- MRCAdata %>%
  filter(acceleration == "1x") %>%
  filter(!is.na(Latency)) %>%
  rstatix::wilcox_test(MRCA_age ~ Subtype, comparisons = my_comparisons) %>%
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


p2 <- MRCAdata %>%
  filter(acceleration == "1x") %>%
  filter(!is.na(Latency)) %>%
  ggplot(aes(Subtype, MRCA_age, fill = Subtype)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray95",
    color = 'black',
    outlier.shape = 21,
    outlier.size = 0.8
  ) +
  ggbeeswarm::geom_quasirandom(
    pch = 21,
    size = 2.5,
    width = 0.2,
    color = "white",
    stroke = 0.2
  ) +
  labs(x = NULL, y = 'MRCA age (years)') +
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


#ggsave(filename = 'MRCA_spgroup.pdf', width = 2, height = 6, device = cairo_pdf)

plot_grid(p2, p1, align = 'h', axis = 'tb', nrow = 1)

ggsave(
  filename = 'MRCA_latency_spgroup.pdf',
  width = 4,
  height = 6,
  device = cairo_pdf
)


MRCAdata %>%
  filter(acceleration == "1x") %>%
  #filter((Subtype %in% c('AS_N','EU_N') & acceleration == "1x")|(Subtype %in% c('EU_S') & acceleration == "7.5x")) %>%
  filter(!is.na(Latency)) %>%
  #filter(Subtype!="Others") %>%
  ggplot(aes(Subtype, Age, fill = Subtype)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  geom_point(
    pch = 21,
    size = 1.5,
    position = position_jitter(width = 0.15, height = 0),
    color = "black"
  ) +
  scale_fill_manual(values = subtypecol) +
  theme_ipsum_rc(axis_title_just = 'm', axis_title_size = 14) +
  labs(x = "", y = 'Age at diagnosis (years)') +
  guides(fill = "none") +
  panel_border(color = 'black', linetype = 1) +
  stat_compare_means(comparisons = my_comparisons)

ggsave(
  filename = 'MRCA_spgroup_age.pdf',
  width = 3.5,
  height = 6,
  device = cairo_pdf
)

MRCAdata %>%
  filter(acceleration == "1x") %>%
  #filter((Subtype %in% c('AS_N','EU_N') & acceleration == "1x")|(Subtype %in% c('EU_S') & acceleration == "5x")) %>%
  filter(!is.na(Latency)) %>%
  do(tidy(wilcox.test(MRCA_age ~ Subtype, data = .)))


# WGD latency -------------------------------------------------------------
library(ggpubr)
my_comparisons <- list(c("WGD", "nWGD"))

MRCAdata %>%
  filter(acceleration == "1x") %>%
  filter(!is.na(Latency)) %>%
  filter(Subtype != "Others") %>%
  left_join(ngspurity %>% select(Tumor_Barcode, WGD_Status = MCN_WGD)) %>%
  ggplot(aes(WGD_Status, Latency, fill = WGD_Status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  geom_point(
    pch = 21,
    size = 1.5,
    position = position_jitter(width = 0.15, height = 0.05),
    color = "black"
  ) +
  facet_wrap(~Subtype) +
  scale_fill_jama() +
  theme_ipsum_rc(axis_title_just = 'm', axis_title_size = 14) +
  theme(panel.spacing = unit(0.1, "cm")) +
  labs(x = "", y = 'Latency, years before diagnosis') +
  guides(fill = "none") +
  panel_border(color = 'black', linetype = 1) +
  stat_compare_means(comparisons = my_comparisons)

ggsave(
  filename = 'latency_WGD_spgroup.pdf',
  width = 4,
  height = 6,
  device = cairo_pdf
)


#  Latency difference among all genome alteration or features--------------------------------------------------------------
tdata <- MRCAdata %>%
  filter(!is.na(Latency)) %>%
  select(Tumor_Barcode, Acc = acceleration, Latency) #Histology == "Adenocarcinoma"; acceleration == '1x'
#tdata <- MRCAdata %>% filter(!is.na(Latency),acceleration == '1x',Subtype %in% c('N_A','N_U')) %>% select(Tumor_Barcode,Latency) #Histology == "Adenocarcinoma"

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
  group_by(Subtype, Acc, Gene) %>%
  do(tresult = safely(wilcox.test)(Latency ~ Alteration, data = .)) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Subtype, Acc, Gene, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  arrange(p.value) %>%
  mutate(TMP = Gene) %>%
  separate(col = 'TMP', into = c('Gene_short', 'Type'), sep = '\\|') %>%
  ungroup()


# add diff
tdata_all <- tdata_all %>%
  mutate(
    Alteration = case_when(
      Alteration %in% c("Female", "Y", "Yes", "Non-Smoker", "WGD") ~ "Yes",
      Alteration %in% c("Male", "N", "No", "Smoker", "nWGD") ~ "No",
      TRUE ~ NA_character_
    )
  )
tmp <- tdata_all %>%
  filter(!is.na(Alteration)) %>%
  group_by(Subtype, Acc, Gene, Alteration) %>%
  summarise(value = median(Latency, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'Alteration', values_from = value) %>%
  mutate(diff = Yes - No)

# add freq
tmpsize <- tdata_all %>%
  drop_na() %>%
  select(Subtype, Acc, Tumor_Barcode, Gene) %>%
  unique() %>%
  group_by(Subtype, Acc) %>%
  count(Gene) %>%
  dplyr::rename(size = n)
freq_mrca <- tdata_all %>%
  drop_na() %>%
  count(Subtype, Acc, Gene, Alteration) %>%
  filter(!is.na(Alteration)) %>%
  left_join(tmpsize) %>%
  mutate(Freq = n / size) %>%
  group_by(Subtype, Acc, Gene) %>%
  arrange(Freq) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  select(Subtype, Acc, Gene, Freq) %>%
  mutate(Freq = if_else(Freq == 1, 0, Freq))

# tmpsize <- mdata %>% select(Tumor_Barcode,name) %>% unique()  %>% count(name)
# freq_mrca <- mdata %>% count(name,value) %>% filter(!is.na(value)) %>% left_join(tmpsize)%>% mutate(Freq=n/size) %>% group_by(name) %>% arrange(Freq) %>% dplyr::slice(1) %>% ungroup() %>% select(name,Freq) %>% mutate(Freq=if_else(Freq==1,0,Freq))

tresult_all <- tresult_all %>%
  left_join(tmpsize) %>%
  left_join(tmp) %>%
  left_join(freq_mrca)

#saveRDS(tresult_all,file='latency_tresult.RDS')
tdata_all <- tdata_all %>%
  mutate(TMP = Gene) %>%
  separate(col = 'TMP', into = c('Gene_short', 'Type'), sep = '\\|')

save(tresult_all, tdata_all, file = 'latency_tresult.RData')


# visualizaiton
tresult_all %>%
  filter(Acc == '1x') %>%
  group_by(Type) %>%
  mutate(FDR = p.adjust(p.value, method = 'BH')) %>%
  ungroup() %>%
  ggplot(aes(diff, -log10(FDR), fill = Type)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'red', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, col = '#cccccc', size = 0.5) +
  geom_point(pch = 21, stroke = 0.2, size = 3) +
  scale_fill_nejm() +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  labs(
    x = 'Latency difference (years)',
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
  filename = 'latency_fdr_alterations.pdf',
  width = 5,
  height = 4,
  device = cairo_pdf
)


# Signatures Only ---------------------------------------------------------
load('Chronological_timing_short.RData', verbose = T)
load('tgct_mutational_signatures.RData', verbose = T)
load('ngspurity_data.RData')

tdata0 <- MRCAdata %>%
  filter(!is.na(Latency), acceleration == '1x') %>%
  select(Tumor_Barcode, Acc = acceleration, Latency) #Histology == "Adenocarcinoma"; acceleration == '1x'

tdata <- tgct_activity_all %>%
  pivot_longer(-Tumor_Barcode) %>%
  left_join(tdata0) %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  left_join(ngspurity %>% select(Tumor_Barcode, BB_Purity)) %>%
  filter(!is.na(Latency))

tresult_all <- tdata %>%
  mutate(value = log2(value + 1)) %>%
  group_by(name) %>%
  do(tresult = safely(lm)(Latency ~ value + BB_Purity + Subtype, data = .)) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(name, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  filter(term == 'value') %>%
  arrange(p.value) %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

# tresult_all <- tdata %>%
#   group_by(name,Subtype) %>%
#   do(tresult = safely(cor.test)(.$Latency,log2(.$value+1))) %>%
#   mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
#   filter(!tresult_null) %>%
#   mutate(fit = list(tidy(tresult[['result']]))) %>%
#   select(Subtype,name,fit) %>%
#   unnest(cols = c(fit)) %>%
#   ungroup() %>%
#   group_by(Subtype) %>%
#   mutate(FDR=p.adjust(p.value,method='fdr')) %>%
#   arrange(FDR)

# tresult_all <- tdata %>%
#   group_by(name) %>%
#   do(tresult = safely(cor.test)(.$Latency,.$value)) %>%
#   mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
#   filter(!tresult_null) %>%
#   mutate(fit = list(tidy(tresult[['result']]))) %>%
#   select(name,fit) %>%
#   unnest(cols = c(fit)) %>%
#   ungroup() %>%
#   arrange(p.value) %>%
#   mutate(FDR=p.adjust(p.value,method='BH'))

tresult_all %>%
  ggplot(aes((estimate), -log10(FDR), fill = name)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'red', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, col = '#cccccc', size = 0.5) +
  geom_point(pch = 21, stroke = 0.2, size = 3) +
  scale_fill_manual(values = sigcol) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  ggrepel::geom_text_repel(
    data = tresult_all %>% filter(FDR < 0.05),
    aes(label = name),
    max.overlaps = 30,
    size = 3
  ) +
  labs(
    x = 'Linear regression coefficient',
    y = '-log10(FDR)',
    fill = 'Signature'
  ) +
  guides(fill = 'none') +
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
  filename = 'latency_fdr_signatures.pdf',
  width = 5,
  height = 4,
  device = cairo_pdf
)


tdata %>%
  filter(name == 'SBS87') %>%
  ggplot(aes(log2(value + 1), Latency, fill = Subtype)) +
  geom_point(pch = 21, size = 2.5) +
  stat_smooth(aes(fill = Subtype, color = Subtype), method = "lm") +
  stat_cor(size = 5, method = 'spearman', col = ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  facet_wrap(~Subtype, scale = 'free') +
  labs(x = 'log2(SBS87 + 1)', y = 'Tumor latency years') +
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


ggsave(
  filename = 'SBS87_subtype_latency.pdf',
  width = 9,
  height = 4,
  device = cairo_pdf
)


tdata <- tgct_activity_all_obs %>%
  pivot_longer(-Tumor_Barcode) %>%
  left_join(tdata0) %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  filter(!is.na(Latency))


tresult_all <- tdata %>%
  group_by(Subtype, name) %>%
  do(tresult = safely(wilcox.test)(Latency ~ value, data = .)) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(Subtype, name, fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  arrange(p.value) %>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))


resulttmp <- tresult_all %>%
  filter(Type == 'Signature_Cosmic_final', !str_detect(Gene, 'APOBEC')) %>%
  filter(Acc == '1x') %>%
  filter(Subtype == 'ALL') %>%
  filter(Freq > 0.01, Gene_short != 'SBS5') %>%
  mutate(
    Sig = if_else(
      str_detect(Gene, '^SBS'),
      'SBS',
      if_else(
        str_detect(Gene, '^DBS'),
        'DBS',
        if_else(str_detect(Gene, '^ID'), 'ID', NA_character_)
      )
    )
  ) %>%
  mutate(Sig = factor(Sig, levels = c('SBS', 'ID', 'DBS'))) %>%
  group_by(Sig) %>%
  mutate(FDR = p.adjust(p.value)) %>%
  ungroup()

#filter(p.value<0.1)
#filter(str_detect(Gene_short,'^SBS'))

datatmp <- tdata_all %>%
  filter(Type == 'Signature_Cosmic_final') %>%
  filter(Acc == '1x') %>%
  filter(Subtype == 'ALL') %>%
  filter(Gene_short %in% resulttmp$Gene_short) %>%
  mutate(
    Sig = if_else(
      str_detect(Gene, '^SBS'),
      'SBS',
      if_else(
        str_detect(Gene, '^DBS'),
        'DBS',
        if_else(str_detect(Gene, '^ID'), 'ID', NA_character_)
      )
    )
  ) %>%
  mutate(Sig = factor(Sig, levels = c('SBS', 'ID', 'DBS')))

#tmp <- datatmp %>% count(Gene,Gene_short,Type,Alteration) %>% filter(n>5) %>% count(Gene,Gene_short,Type) %>% filter(n==2) %>% pull(Gene)
#datatmp <- datatmp %>% filter(Gene %in% tmp)
#resulttmp <- resulttmp %>% filter(Gene %in% tmp)

tmplevel <- datatmp %>%
  filter(Alteration == "Yes") %>%
  group_by(Gene_short) %>%
  summarise(value = median(Latency)) %>%
  arrange(value) %>%
  pull(Gene_short)

resulttmp <- resulttmp %>%
  mutate(Gene_short = factor(Gene_short, levels = tmplevel))
datatmp <- datatmp %>%
  mutate(Gene_short = factor(Gene_short, levels = tmplevel))

# barplot
p1 <- datatmp %>%
  #filter(Gene_short=='ID-Novel-C') %>%
  ggplot(aes(Gene_short, Latency, fill = Alteration)) +
  geom_point(
    pch = 21,
    size = 1,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
    color = "black"
  ) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  scale_fill_manual(values = c("Yes" = "#4daf4a", "No" = "#cccccc")) +
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "Y",
    grid_col = "gray90",
    ticks = T
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = 'top',
    panel.spacing.x = unit(0.1, 'cm')
  ) +
  labs(x = "", y = 'Latency, years before diagnosis') +
  scale_y_continuous(breaks = pretty_breaks()) + #,limits = c(-1.5,32)
  #guides(fill="none")+
  facet_grid(. ~ Sig, scales = 'free', space = 'free') +
  panel_border(color = 'black', linetype = 1)

p2 <- resulttmp %>%
  ggplot(aes(Gene_short, "1", fill = -log10(p.value))) +
  geom_tile() +
  facet_grid(. ~ Sig, scales = 'free', space = 'free') +
  scale_fill_viridis_c(option = "D") +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = FALSE,
    grid_col = "gray90",
    ticks = F
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = 'top',
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.3, "cm"),
    panel.spacing.x = unit(0.1, 'cm')
  ) +
  labs(x = "", y = '', fill = "Wilcox test, -log10(P)\n")
#guides(fill="none")+
# panel_border(color = 'black',linetype = 1)

plot_grid(p2, p1, axis = 'v', align = 'lr', ncol = 1, rel_heights = c(1.1, 4))

ggsave(
  filename = 'latency_signatures.pdf',
  width = 16,
  height = 10,
  device = cairo_pdf
)

## vacanol plot

resulttmp %>%
  mutate(Sig = as.factor(Sig)) %>%
  ggplot(aes((diff), -log10(FDR), fill = Sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'red', size = 0.5) +
  geom_vline(xintercept = c(-5, 5), linetype = 2, col = 'blue', size = 0.5) +
  geom_point(aes(size = Freq), pch = 21, stroke = 0.2) +
  scale_size_binned(labels = percent_format()) +
  scale_fill_manual(values = ncicolpal[c(1, 3, 4)]) +
  scale_x_continuous(breaks = pretty_breaks()) +
  ggrepel::geom_text_repel(
    data = resulttmp,
    aes(label = Gene_short),
    max.overlaps = 30,
    size = 3
  ) +
  labs(
    x = 'Latency Change (Years)',
    y = '-log10(FDR)',
    size = 'Frequency',
    fill = 'Signature'
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.5))) +
  theme_ipsum_rc(
    base_size = 12,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = 'XY',
    ticks = T
  ) +
  panel_border(color = 'black', size = 0.5) +
  coord_cartesian(clip = 'off')


ggsave(
  filename = 'latency_fdr_signatures.pdf',
  width = 7,
  height = 5.5,
  device = cairo_pdf
)


# OLD -------------
# Driver Mutations --------------------------------------------------------
load('latency_tresult.RData')
drglist <- readRDS(
  '../../../Collaborators/Nuria/Update2/drivers_intogene.RDS'
) %>%
  pull(symbol)

resulttmp <- tresult_all %>%
  filter(Type == 'Mutation_Driver', Gene_short %in% drglist) %>%
  filter(Acc == '1x') %>%
  filter(Subtype == 'ALL') %>%
  filter(Freq > 0.01) %>%
  mutate(FDR = p.adjust(p.value))

#filter(p.value<0.1)
#filter(str_detect(Gene_short,'^SBS'))

datatmp <- tdata_all %>%
  filter(Type == 'Mutation_Driver', Gene_short %in% drglist) %>%
  filter(Acc == '1x') %>%
  filter(Subtype == 'ALL') %>%
  filter(Gene_short %in% resulttmp$Gene_short)

#tmp <- datatmp %>% count(Gene,Gene_short,Type,Alteration) %>% filter(n>5) %>% count(Gene,Gene_short,Type) %>% filter(n==2) %>% pull(Gene)
#datatmp <- datatmp %>% filter(Gene %in% tmp)
#resulttmp <- resulttmp %>% filter(Gene %in% tmp)

tmplevel <- datatmp %>%
  filter(Alteration == "Yes") %>%
  group_by(Gene_short) %>%
  summarise(value = median(Latency)) %>%
  arrange(value) %>%
  pull(Gene_short)

resulttmp <- resulttmp %>%
  mutate(Gene_short = factor(Gene_short, levels = tmplevel))
datatmp <- datatmp %>%
  mutate(Gene_short = factor(Gene_short, levels = tmplevel))
## vacanol plot

resulttmp %>%
  ggplot(aes((diff), -log10(FDR))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'red', size = 0.5) +
  geom_vline(xintercept = c(-5, 5), linetype = 2, col = 'blue', size = 0.5) +
  geom_point(
    aes(size = Freq),
    pch = 21,
    stroke = 0.1,
    fill = '#3D4551',
    col = 'white'
  ) +
  scale_size_binned(labels = percent_format()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  ggrepel::geom_text_repel(
    data = resulttmp %>% filter(FDR < 0.2),
    aes(label = Gene_short),
    max.overlaps = 30,
    size = 4
  ) +
  labs(
    x = 'Latency Change (Years)',
    y = '-log10(FDR)',
    size = 'Mutation frequency',
    fill = 'Signature'
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.5))) +
  theme_ipsum_rc(
    base_size = 12,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = 'XY',
    ticks = T
  ) +
  theme(legend.position = 'top', legend.key.width = unit(1, 'cm')) +
  panel_border(color = 'black', size = 0.5) +
  coord_cartesian(clip = 'off')


resulttmp %>%
  ggplot(aes((diff), -log10(p.value))) +
  geom_hline(
    yintercept = -log10(0.001),
    linetype = 2,
    col = 'red',
    size = 0.5
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    col = 'orange',
    size = 0.5
  ) +
  geom_vline(xintercept = c(-5, 5), linetype = 2, col = 'blue', size = 0.5) +
  geom_point(
    aes(size = Freq),
    pch = 21,
    stroke = 0.1,
    fill = '#3D4551',
    col = 'white'
  ) +
  scale_size_binned(labels = percent_format()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  ggrepel::geom_text_repel(
    data = resulttmp %>% filter(p.value < 0.05),
    aes(label = Gene_short),
    max.overlaps = 30,
    size = 3.5
  ) +
  labs(
    x = 'Latency Change (Years)',
    y = '-log10(FDR)',
    size = 'Mutation frequency',
    fill = 'Signature'
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.5))) +
  theme_ipsum_rc(
    base_size = 12,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = 'XY',
    ticks = T
  ) +
  panel_border(color = 'black', size = 0.5) +
  coord_cartesian(clip = 'off')

ggsave(
  filename = 'latency_driverGene_signatures.pdf',
  width = 5,
  height = 5.5,
  device = cairo_pdf
)

## for EGFR
tmp <- resulttmp <- tresult_all %>%
  filter(Type == 'Mutation_Driver', Gene_short %in% 'EGFR') %>%
  filter(Acc == '1x') %>%
  filter(Subtype %in% c('AS_N', 'EU_N', 'EU_S', "Others")) %>%
  mutate(label = paste0(Subtype, '\nP=', scientific_format()(p.value))) %>%
  select(Subtype, label)

tdata_all %>%
  filter(Type == 'Mutation_Driver', Gene_short %in% 'EGFR') %>%
  filter(Acc == '1x') %>%
  filter(Subtype %in% c('AS_N', 'EU_N', 'EU_S', "Others")) %>%
  mutate(Aternation = factor(Alteration, levels = c('Yes', 'No'))) %>%
  left_join(tmp) %>%
  ggplot(aes(Alteration, Latency, fill = Alteration)) +
  geom_point(
    pch = 21,
    size = 2,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
    color = "black"
  ) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  scale_fill_manual(values = c("Yes" = "#4daf4a", "No" = "#cccccc")) +
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "Y",
    grid_col = "gray90",
    ticks = T
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = 'bottom',
    panel.spacing.x = unit(0.1, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 2)
  ) +
  labs(x = "", y = 'Latency, years before diagnosis', fill = 'EGFR Mutation') +
  scale_y_continuous(breaks = pretty_breaks()) + #,limits = c(-1.5,32)
  #guides(fill="none")+
  facet_grid(. ~ label, scales = 'free', space = 'free') +
  panel_border(color = 'black', linetype = 1)

ggsave(
  filename = 'latency_driverGene_signatures2.pdf',
  width = 6,
  height = 6,
  device = cairo_pdf
)


## for KRAS
tmp <- resulttmp <- tresult_all %>%
  filter(Type == 'Mutation_Driver', Gene_short %in% 'KRAS') %>%
  filter(Acc == '1x') %>%
  filter(Subtype %in% c('AS_N', 'EU_N', 'EU_S', "Others")) %>%
  mutate(label = paste0(Subtype, '\nP=', scientific_format()(p.value))) %>%
  select(Subtype, label)


tdata_all %>%
  filter(Type == 'Mutation_Driver', Gene_short %in% 'KRAS') %>%
  filter(Acc == '1x') %>%
  filter(Subtype %in% c('AS_N', 'EU_N', 'EU_S', 'Others')) %>%
  mutate(Aternation = factor(Alteration, levels = c('Yes', 'No'))) %>%
  left_join(tmp) %>%
  ggplot(aes(Alteration, Latency, fill = Alteration)) +
  geom_point(
    pch = 21,
    size = 2,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
    color = "black"
  ) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  scale_fill_manual(values = c("Yes" = "#4daf4a", "No" = "#cccccc")) +
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "Y",
    grid_col = "gray90",
    ticks = T
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = 'bottom',
    panel.spacing.x = unit(0.1, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 2)
  ) +
  labs(x = "", y = 'Latency, years before diagnosis', fill = 'KRAS Mutation') +
  scale_y_continuous(breaks = pretty_breaks()) + #,limits = c(-1.5,32)
  #guides(fill="none")+
  facet_grid(. ~ label, scales = 'free', space = 'free') +
  panel_border(color = 'black', linetype = 1)

ggsave(
  filename = 'latency_driverGene_signatures_kras.pdf',
  width = 6,
  height = 6,
  device = cairo_pdf
)

# example of regression
load('../covdata0.RData')
tdata_all %>%
  filter(Type == 'Mutation_Driver', Gene_short %in% 'EGFR') %>%
  filter(Acc == '1x') %>%
  left_join(covdata0) %>%
  filter(Subtype %in% c('AS_N', 'EU_N', 'EU_S')) %>%
  group_by(Subtype) %>%
  do(tidy(lm(
    Latency ~ Alteration + Gender + Histology + Tumor_Purity,
    data = .
  ))) %>%
  filter(term == 'AlterationYes') %>%
  arrange(p.value)
#do(tidy(t.test(Latency~Alteration,data=.)))

# Latency for Gender ------------------------------------------------------

tmp <- tgct_data_full %>%
  filter(Type == 'Mutation_Driver', Gene == 'EGFR') %>%
  select(Tumor_Barcode, EGFR = Alteration)

tdata <- MRCAdata %>%
  filter(Tumor_Barcode %in% hq_samples2) %>%
  filter(acceleration == '1x', Tumor_Barcode %in% hq_samples2) %>%
  left_join(Subtype_data2) %>%
  left_join(covdata0 %>% select(Tumor_Barcode, Gender, Tumor_Purity)) %>%
  filter(Subtype %in% c('AS_N', 'EU_N', 'EU_S', 'Others')) %>%
  left_join(tmp)

tdata %>% do(tidy(wilcox.test(Latency ~ Gender, data = .)))
tmp <- tdata %>%
  group_by(Subtype) %>%
  do(tidy(wilcox.test(Latency ~ Gender, data = .))) %>%
  mutate(label = paste0(Subtype, '\nP=', scientific_format()(p.value))) %>%
  select(Subtype, label)

tdata %>%
  left_join(tmp) %>%
  ggplot(aes(Gender, Latency, fill = Gender)) +
  geom_point(
    pch = 21,
    size = 2,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
    color = "black"
  ) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  scale_fill_manual(values = c("Female" = "#4daf4a", "Male" = "#cccccc")) +
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "Y",
    grid_col = "gray90",
    ticks = T
  ) +
  #theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'top',panel.spacing.x = unit(0.1,'cm'))+
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = 'bottom',
    panel.spacing.x = unit(0.1, 'cm'),
    strip.text.x = element_text(hjust = 0.5, face = 2)
  ) +
  labs(x = "", y = 'Latency, years before diagnosis', fill = 'Sex') +
  scale_y_continuous(breaks = pretty_breaks()) + #,limits = c(-1.5,32)
  facet_grid(. ~ label, scales = 'free', space = 'free') +
  panel_border(color = 'black', linetype = 1)

ggsave(
  filename = 'latency_gender_signatures.pdf',
  width = 6,
  height = 5.8,
  device = cairo_pdf
)


p1 <- tdata %>%
  filter(Subtype == 'EU_N') %>%
  ggplot(aes(Gender, Latency, fill = Gender)) +
  geom_point(
    pch = 21,
    size = 2,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
    color = "black"
  ) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  scale_fill_manual(values = c("Female" = "#984ea3", "Male" = "#cccccc")) +
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "Y",
    grid_col = "gray90",
    ticks = T
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = 'top',
    panel.spacing.x = unit(0.1, 'cm')
  ) +
  labs(x = "", y = 'Latency, years before diagnosis', fill = 'Sex') +
  scale_y_continuous(breaks = pretty_breaks()) + #,limits = c(-1.5,32)
  facet_grid(. ~ EGFR, scales = 'free', space = 'free') +
  panel_border(color = 'black', linetype = 1)


p2 <- tdata %>%
  filter(Subtype == 'EU_N') %>%
  ggplot(aes(EGFR, Latency, fill = EGFR)) +
  geom_point(
    pch = 21,
    size = 2,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
    color = "black"
  ) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  scale_fill_manual(values = c("Yes" = "#4daf4a", "No" = "#cccccc")) +
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "Y",
    grid_col = "gray90",
    ticks = T
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = 'top',
    panel.spacing.x = unit(0.1, 'cm')
  ) +
  labs(x = "", y = 'Latency, years before diagnosis', fill = 'EGFR Mutation') +
  scale_y_continuous(breaks = pretty_breaks()) + #,limits = c(-1.5,32)
  facet_grid(. ~ Gender, scales = 'free', space = 'free') +
  panel_border(color = 'black', linetype = 1)

plot_grid(p1, p2, align = 'v', axis = 'lr')


ggsave(
  filename = 'latency_gender_EGFR_signatures2.pdf',
  width = 7.5,
  height = 5.5,
  device = cairo_pdf
)

tdata %>%
  filter(Subtype == 'EU_N') %>%
  group_by(Gender) %>%
  do(tidy(wilcox.test(Latency ~ EGFR, data = .)))

tdata %>%
  filter(Subtype == 'EU_N') %>%
  group_by(EGFR) %>%
  do(tidy(wilcox.test(Latency ~ Gender, data = .)))


tdata %>%
  group_by(Subtype) %>%
  do(tidy(wilcox.test(Latency ~ Gender, data = .)))

tmpresult <- tdata %>%
  group_by(Subtype) %>%
  do(tidy(
    lm(Latency ~ Gender + EGFR + Tumor_Purity, data = .),
    conf.int = TRUE,
    conf.level = 0.95
  )) %>%
  filter(term != '(Intercept)') %>%
  ungroup() %>%
  arrange(p.value) %>%
  mutate(term = str_remove(term, "EGFR")) %>%
  mutate(term = str_remove(term, "Histology")) %>%
  mutate(term = str_remove(term, "Gender"))


tmpresult <- tmpresult %>%
  mutate(
    label = paste0(
      'β = ',
      round(estimate, 2),
      ', p = ',
      scientific_format(digits = 3)(p.value)
    )
  )
# %>%
#   bind_rows(
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(Subtype='AS_N'),
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(Subtype='EU_N'),
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(Subtype='EU_S')
#   )

# tmplevels <- c('Tumor_Purity','Male','Female','No','Yes')
# tmplables <- c('Tumor purity', 'Male','Female','EGFR wildtype','EGFR mutant')

tmplevels <- c('Tumor_Purity', 'Female', 'Yes')
tmplables <- c('Tumor purity', 'Female\nRef=Male', 'EGFR Mutant\nRef=Wildtype')

tmpresult %>%
  mutate(term = factor(term, levels = tmplevels, labels = tmplables)) %>%
  ggplot(aes(
    estimate,
    term,
    xmin = conf.low,
    xmax = conf.high,
    height = 0,
    color = p.value < 0.05
  )) +
  geom_point(pch = 19, size = 3) +
  geom_vline(xintercept = 0, lty = 4) +
  geom_errorbarh(height = .3) +
  scale_color_manual(values = c("FALSE" = 'black', "TRUE" = 'red')) +
  ggrepel::geom_text_repel(aes(label = label), size = 3.5, nudge_y = 0.2) +
  facet_grid(~Subtype) +
  theme_ipsum_rc(axis_title_just = 'm', axis_title_size = 14) +
  theme(panel.spacing = unit(0.2, "cm")) +
  labs(x = "Regression coefficient for tumor latency", y = NULL) +
  guides(color = "none") +
  scale_x_continuous(breaks = pretty_breaks()) +
  panel_border(color = 'black', linetype = 1)

ggsave(
  filename = 'latency_gender_EGFR_signatures.pdf',
  width = 14,
  height = 3,
  device = cairo_pdf
)


# OLD visualization for all groups and target genes -----------------------

# tresult <- tresult %>% filter(Freq>0.03)
## for trended on only major alterations
tmp <- tresult %>% filter(Yes > 5 | Yes < 2 | (p.value < 0.05 & Freq > 0.2))
tresult %>%
  ggplot(aes((Yes), -log10(p.value), fill = Type)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'red', size = 0.5) +
  geom_point(aes(size = Freq), pch = 21, stroke = 0.2) +
  scale_fill_manual(values = tgct_type_colors[sort(unique(tresult$Type))]) +
  scale_size_binned() +
  #ggrepel::geom_text_repel(data=tmp,aes(label=Gene_short),max.overlaps = 30)+
  labs(
    x = 'Latency, years before diagnosis',
    y = 'Wilcox test, -log10(pvalue)'
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.5))) +
  theme_ipsum_rc(
    base_size = 12,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "XY",
    ticks = T
  ) +
  panel_border(color = 'black', size = 0.5) +
  coord_cartesian(clip = 'off')

tmp <- tresult %>%
  filter(abs(diff) > 5 | (p.value < 0.05 & Freq > 0.2) | p.value < 0.01)
tresult %>%
  ggplot(aes((diff), -log10(p.value), fill = Type)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = 'red', size = 0.5) +
  geom_point(aes(size = Freq), pch = 21, stroke = 0.2) +
  scale_fill_manual(values = tgct_type_colors[sort(unique(tresult$Type))]) +
  scale_size_binned() +
  ggrepel::geom_text_repel(
    data = tmp,
    aes(label = Gene_short),
    max.overlaps = 30,
    size = 3
  ) +
  labs(x = 'Latency Change (Year)', y = 'Wilcox test, -log10(pvalue)') +
  guides(fill = guide_legend(override.aes = list(size = 3.5))) +
  theme_ipsum_rc(
    base_size = 12,
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "XY",
    ticks = T
  ) +
  panel_border(color = 'black', size = 0.5) +
  coord_cartesian(clip = 'off')

ggsave(
  filename = 'latency_tgct_all.pdf',
  width = 12,
  height = 8,
  device = cairo_pdf
)


# barplot #
tmp <- tresult %>%
  filter(qvalue < 0.12) %>%
  select(Gene, p.value, Gene_short, Type) %>%
  left_join(tdata) %>%
  mutate(
    Alteration = if_else(
      Alteration == "Female",
      "Yes",
      if_else(Alteration == "Male", "No", Alteration)
    )
  ) %>%
  mutate(
    Gene_short = if_else(Gene_short == "Gender", "Gender_Female", Gene_short)
  ) %>%
  mutate(Gene_short = paste0(Gene_short, "\n", Type))

tmplevel <- tmp %>%
  filter(Alteration == "Yes") %>%
  group_by(Gene_short) %>%
  summarise(value = median(Latency)) %>%
  arrange(value) %>%
  pull(Gene_short)

tmp <- tmp %>% mutate(Gene_short = factor(Gene_short, levels = tmplevel))
p1 <- tmp %>%
  left_join(wgs_groups_info %>% select(Tumor_Barcode, Subtype)) %>%
  left_join(Subtype_data) %>%
  filter(Alteration == 'No') %>%
  ggplot(aes(Gene_short, Latency, fill = Subtype)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  geom_point(
    pch = 21,
    size = 1,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
    color = "black"
  ) +
  #scale_fill_manual(values = c("Yes" = "#008B45FF","No" = "#cccccc"))+
  scale_fill_manual(values = subtypecol) +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "Y",
    grid_col = "gray90",
    ticks = T
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = 'top'
  ) +
  labs(x = "", y = 'Latency, years before diagnosis') +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(-1.5, 32)) +
  #guides(fill="none")+
  panel_border(color = 'black', linetype = 1)

p2 <- tmp %>%
  filter(Alteration == "Yes") %>%
  ggplot(aes(Gene_short, "1", fill = -log10(p.value))) +
  geom_tile() +
  scale_fill_viridis_c(option = "D") +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = FALSE,
    grid_col = "gray90",
    ticks = F
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = 'top',
    legend.key.width = unit(1, "cm")
  ) +
  labs(x = "", y = '', fill = "Wilcox test, -log10(pvalue)\n")
#guides(fill="none")+
# panel_border(color = 'black',linetype = 1)

plot_grid(p2, p1, axis = 'v', align = 'lr', ncol = 1, rel_heights = c(1.1, 4))

ggsave(
  filename = 'latency_fdr01.pdf',
  width = 16,
  height = 9,
  device = cairo_pdf
)


# targeable alterations
genelist <- c(
  'TP53',
  'EGFR',
  'KRAS',
  'HER2',
  "BRAF",
  "ALK",
  "RET",
  "ROS1",
  "NTRK",
  "MET",
  "ERBB2",
  "NF1",
  "MAPK2K1",
  "KRAS",
  "FRFRG2",
  "3p21.31",
  '6q22.31',
  '9p21.3',
  '12p11.22',
  '13q21.33',
  '16q24.1',
  '17p12',
  '18q22.3',
  '22q12.2',
  '5p15.33',
  '8q24.21',
  '10p11.1',
  '12q15',
  '14q21.1',
  '19q13.2',
  'NRAS',
  'RIT1',
  'FGFR1'
)

tmp <- tresult %>%
  filter(Type %in% c('Mutation_Driver', 'SCNA_Focal_Cytoband', 'Fusion')) %>%
  filter(Gene_short %in% genelist) %>%
  select(Gene, p.value, Gene_short, Type) %>%
  left_join(tdata) %>%
  mutate(
    Alteration = if_else(
      Alteration == "Female",
      "Yes",
      if_else(Alteration == "Male", "No", Alteration)
    )
  ) %>%
  mutate(
    Gene_short = if_else(Gene_short == "Gender", "Gender_Female", Gene_short)
  ) %>%
  mutate(Gene_short = paste0(Gene_short, "\n", Type))

exclude_tmp <- tmp %>%
  filter(Alteration == 'Yes') %>%
  count(Gene) %>%
  filter(n < 5) %>%
  pull(Gene)
tmp <- tmp %>% filter(!(Gene %in% exclude_tmp))

tmplevel <- tmp %>%
  filter(Alteration == "Yes") %>%
  group_by(Gene_short) %>%
  summarise(value = median(Latency, na.rm = T)) %>%
  arrange(value) %>%
  pull(Gene_short)

tmp <- tmp %>% mutate(Gene_short = factor(Gene_short, levels = tmplevel))
p1 <- tmp %>%
  ggplot(aes(Gene_short, Latency, fill = Alteration)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25) +
  geom_point(
    pch = 21,
    size = 1,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
    color = "black"
  ) +
  scale_fill_manual(values = c("Yes" = "#FF7F0EFF", "No" = "#cccccc")) +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = "Y",
    grid_col = "gray90",
    ticks = T
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = 'top'
  ) +
  labs(x = "", y = 'Latency, years before diagnosis') +
  scale_y_continuous(breaks = pretty_breaks()) +
  #guides(fill="none")+
  panel_border(color = 'black', linetype = 1)

p2 <- tmp %>%
  filter(Alteration == "Yes") %>%
  ggplot(aes(Gene_short, "1", fill = -log10(p.value))) +
  geom_tile() +
  scale_fill_viridis_c(option = "D") +
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid = FALSE,
    grid_col = "gray90",
    ticks = F
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = 'top',
    legend.key.width = unit(1, "cm")
  ) +
  labs(x = "", y = '', fill = "Wilcox test, -log10(pvalue)\n")
#guides(fill="none")+
# panel_border(color = 'black',linetype = 1)

plot_grid(p2, p1, axis = 'v', align = 'lr', ncol = 1, rel_heights = c(1.1, 4))

ggsave(
  filename = 'latency_targetable.pdf',
  width = 16,
  height = 9,
  device = cairo_pdf
)


# SBS17b ------------------------------------------------------------------
tgct_data_full %>%
  filter(Gene == 'SBS17b') %>%
  left_join(clinical_data) %>%
  ggplot(aes(Alteration, Age)) +
  geom_boxplot()

tgct_data_full %>%
  filter(Gene == 'SBS17b') %>%
  left_join(clinical_data) %>%
  do(tidy(wilcox.test(Age ~ Alteration, data = .)))

load('../RDS/tgct_variable.RData')

tgct_variable %>%
  filter(name == 'SBS17b') %>%
  left_join(clinical_data) %>%
  filter(value > 0) %>%
  ggplot(aes(log2(value), Age)) +
  geom_point() +
  geom_smooth(method = 'lm')

tgct_variable %>%
  filter(name == 'SBS17b') %>%
  left_join(clinical_data) %>%
  filter(value > 0) %>%
  do(tidy(cor.test(.$Age, .$value)))


# SBS17b ------------------------------------------------------------------

# Germline L1 -------------------------------------------------------------
load('Chronological_timing_short.RData')
load('../RDS/tgct_data_all.RData')
load('../BBsolution_final3_short.RData')

tgct_data_full

tdata <- MRCAdata %>%
  filter(acceleration == '1x') %>%
  left_join(
    tgct_data_full %>%
      filter(str_detect(Gene, 'All_L1')) %>%
      select(Tumor_Barcode, Gene, Alteration)
  ) %>%
  filter(Tumor_Barcode %in% hq_samples2)

tdata %>%
  group_by(Gene) %>%
  do(tidy(wilcox.test(Latency ~ Alteration, data = .)))

tdata %>%
  ggplot(aes(Gene, Latency, fill = Alteration)) +
  geom_point(
    pch = 21,
    size = 2,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
    stroke = 0.05,
    col = 'white'
  ) +
  geom_boxplot(width = 0.5, pch = 21, outlier.shape = NA, alpha = 0.3) +
  scale_fill_manual(values = id2color) +
  theme_ipsum_rc(axis_title_just = 'm', axis_title_size = 14) +
  labs(x = "", y = 'Latency, years before diagnosis') +
  panel_border(color = 'black', linetype = 1)

ggsave(
  file = 'L1_latency_diff.pdf',
  width = 6,
  height = 5,
  device = cairo_pdf()
)


# Multivariate analysis ---------------------------------------------------
load('../ID2_TE/id2data.RData')

tmp1 <- tgct_data_full %>%
  filter(Type == 'Mutation_Driver', Gene == 'EGFR') %>%
  select(Tumor_Barcode, EGFR = Alteration)

tmp2 <- tgct_data_full %>%
  filter(Type == 'Mutation_Driver', Gene == 'KRAS') %>%
  select(Tumor_Barcode, KRAS = Alteration)

tdata <- MRCAdata %>%
  filter(Tumor_Barcode %in% hq_samples2) %>%
  filter(acceleration == '1x', Tumor_Barcode %in% hq_samples2) %>%
  left_join(Subtype_data2) %>%
  left_join(
    covdata0 %>%
      select(Tumor_Barcode, Gender, Smoking, Assigned_Population, Tumor_Purity)
  ) %>%
  left_join(id2data %>% select(Tumor_Barcode, ID2_Present)) %>%
  #filter(Subtype %in% c('AS_N','EU_N','EU_S')) %>%
  left_join(tmp1) %>%
  left_join(tmp2) %>%
  mutate(Smoking = if_else(Smoking == 'Unknown', NA_character_, Smoking)) %>%
  mutate(Smoking = factor(Smoking, levels = c('Smoker', 'Non-Smoker')))


tmpresult <- tdata %>%
  do(tidy(
    lm(
      Latency ~
        Gender +
          Assigned_Population +
          Smoking +
          ID2_Present +
          EGFR +
          KRAS +
          Tumor_Purity,
      data = .
    ),
    conf.int = TRUE,
    conf.level = 0.95
  )) %>%
  filter(term != '(Intercept)') %>%
  ungroup() %>%
  arrange(desc(estimate)) %>%
  #mutate(term=str_remove(term,"EGFR")) %>%
  mutate(term = str_remove(term, "Histology")) %>%
  mutate(term = str_remove(term, "Gender")) %>%
  mutate(
    label = paste0(
      'β = ',
      round(estimate, 2),
      ', p = ',
      scientific_format(digits = 3)(p.value)
    )
  )


# tmpresult <- tmpresult %>% mutate(label=paste0('β = ',round(estimate,2), ', p = ',scientific_format(digits = 3)(p.value))) %>%
#   bind_rows(
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(Subtype='AS_N'),
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(Subtype='EU_N'),
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(Subtype='EU_S')
#   )

# tmplevels <- c('Tumor_Purity','Male','Female','No','Yes')
#
# tmplables <- c('Tumor purity', 'Male','Female','EGFR wildtype','EGFR mutant')

tmplevels <- tmpresult$term
tmplables <- tmpresult$term

#tmplables <- c('EGFR Mutant\nRef=Wildtype','Female\nRef=Male','Smokers\nRef=Nonsmokers','Tumor Purity','Ancestry EAS\nRef=EUR','Ancestry Others\nRef=EUR','KRAS Mutant\nRef=Wildtype','ID2 Signature Present\nRef=Absent')

tmplables <- c(
  'EGFR Mutant\nRef=Wildtype',
  'Female\nRef=Male',
  'Tumor Purity',
  'Ancestry EAS\nRef=EUR',
  'Nonsmokers\nRef=Smokers',
  'Ancestry Others\nRef=EUR',
  'KRAS Mutant\nRef=Wildtype',
  'ID2 Signature Present\nRef=Absent'
)


tmpresult %>%
  mutate(p.value = if_else(p.value > 0.05, 0.5, p.value)) %>%
  mutate(term = factor(term, levels = tmplevels, labels = tmplables)) %>%
  ggplot(aes(
    estimate,
    term,
    xmin = conf.low,
    xmax = conf.high,
    height = 0,
    color = -log10(p.value)
  )) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3) +
  geom_point(pch = 19, size = 4) +
  #scale_color_manual(values = c("FALSE"='black',"TRUE"='red'))+
  scale_color_material(palette = 'red') +
  ggrepel::geom_text_repel(aes(label = label), size = 3.5, nudge_y = 0.4) +
  #facet_grid(~Subtype)+
  theme_ipsum_rc(
    axis_title_just = 'm',
    axis_title_size = 14,
    grid_col = 'gray85'
  ) +
  theme(panel.spacing = unit(0.2, "cm")) +
  labs(x = "Regression coefficient", y = NULL) +
  guides(color = "none") +
  scale_x_continuous(breaks = pretty_breaks()) +
  panel_border(color = 'black', linetype = 1)

ggsave(
  filename = 'latency_multivariables.pdf',
  width = 8,
  height = 5,
  device = cairo_pdf
)
