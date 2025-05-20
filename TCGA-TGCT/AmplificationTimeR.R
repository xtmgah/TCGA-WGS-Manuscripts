set_wd()
libztw()
load('wgs_data.RData')
library(AmplificationTimeR)
library(BSgenome.Hsapiens.UCSC.hg38)

load('ngspurity_data.RData', verbose = T)
load('DP_info_data.RData', verbose = T)
load('tgct_signature_probability.RData', verbose = T)
# Prepare data ------------------------------------------------------------

clock_mutations <- tgct_sbs96_probability %>%
  filter(SBS1 + SBS5 > 0.5) %>%
  select(Tumor_Barcode, ID) %>%
  mutate(clock = TRUE)
DP_info_data <- DP_info_data %>% left_join(clock_mutations) %>% filter(clock)

segment_time_all <- NULL

for (barcode in ngspurity$Tumor_Barcode) {
  tryCatch(
    {
      amp_chrom <- 12
      # amp_start <- 638760
      # amp_stop <- 21216391
      #
      amp_start <- 50000
      amp_stop <- 35000000

      wgd_status <- ngspurity %>%
        filter(Tumor_Barcode == barcode) %>%
        select(MCN_WGD) %>%
        pull(MCN_WGD)

      wgd_status <- if_else(wgd_status == 'WGD', TRUE, FALSE)

      demo_cn <- BBprofile %>%
        filter(
          Tumor_Barcode %in% barcode,
          chr == amp_chrom,
          startpos <= amp_stop,
          endpos >= amp_start
        ) %>%
        mutate(seglength = abs(endpos - startpos)) %>%
        group_by(Tumor_Barcode) %>%
        arrange(desc(seglength)) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        select(chr, startpos, endpos, nMaj1_A, nMin1_A) %>%
        as.data.frame()

      demo_mult <- DP_info_data %>%
        filter(Tumor_Barcode %in% barcode, CHROM %in% amp_chrom) %>%
        select(chr = CHROM, end = POS, no.chrs.bearing.mut) %>%
        as.data.frame()

      demo_muts <- DP_info_data %>%
        filter(Tumor_Barcode %in% barcode, CHROM %in% amp_chrom) %>%
        select(chr = CHROM, start = POS, end = POS, ref = REF, alt = ALT) %>%
        as.data.frame()

      segment_time <- time_amplification(
        cn_data = demo_cn,
        multiplicity_data = demo_mult,
        mutation_data = demo_muts,
        muts_type = "All",
        sample_id = barcode,
        amplification_chrom = amp_chrom,
        amplification_start = amp_start,
        amplification_stop = amp_stop,
        is_WGD = wgd_status,
        genome = "hg38"
      )

      if (!is.null(segment_time)) {
        segment_time <- as_tibble(segment_time) %>%
          mutate(Tumor_Barcode = barcode) %>%
          relocate(Tumor_Barcode)

        segment_time_all <- bind_rows(segment_time_all, segment_time)
      }
    },
    error = function(e) {
      message("Error for barcode ", barcode, ": ", e$message)
    }
  )
}


segment_time_all2 <- segment_time_all %>%
  filter(
    !event_order %in%
      c('Order uncertain and cannot be timed', 'Cannot be timed')
  ) %>%
  left_join(wgs_data %>% select(Tumor_Barcode, Subtype)) %>%
  relocate(Subtype)


# Pivot to long format
## Gain first

pdata <- segment_time_all2 %>%
  # pivot all columns that start with t_
  pivot_longer(
    cols = starts_with("t_"),
    names_to = c("t", "measure"),
    names_pattern = "t_(\\d+)_?(.*)", # split t_1, t_1_mean_bootstrap etc
    values_to = "value"
  ) %>%
  # fill missing measure with "t_value"
  mutate(measure = if_else(measure == "", "t_value", measure)) %>%
  # pivot wider to get t_value, mean_bootstrap, lower_ci, upper_ci in separate columns
  pivot_wider(
    names_from = measure,
    values_from = value
  ) %>%
  mutate(label = paste0('T', t)) %>%
  filter(!is.na(mean_bootstrap)) %>%
  mutate(t = as.integer(t)) %>%
  mutate(t2 = t) %>%
  mutate(mean_bootstrap = if_else(mean_bootstrap > 1, 1, mean_bootstrap)) %>%
  mutate(lower_ci = if_else(lower_ci > 1, 1, lower_ci)) %>%
  mutate(upper_ci = if_else(upper_ci > 1, 1.01, upper_ci)) %>%
  mutate(
    Tumor_Barcode2 = paste0(
      str_sub(Tumor_Barcode, 1, 12),
      ', nMut=',
      num_mutations_used,
      "\nCNV=",
      highest_copy_number,
      ", Event=",
      event_order
    )
  ) %>%
  mutate(WGD_locate = str_locate(event_order, "W")[, 1]) %>%
  filter(
    str_length(event_order) > 1,
    str_detect(event_order, 'W'),
    str_detect(event_order, 'G')
  ) %>%
  filter(!str_starts(event_order, 'W'))


tmplev <- pdata %>%
  group_by(Subtype, Tumor_Barcode, Tumor_Barcode2) %>%
  summarise(n = n(), maxvalue = max(mean_bootstrap)) %>%
  ungroup() %>%
  arrange(Subtype, desc(n), desc(maxvalue)) %>%
  pull(Tumor_Barcode2)

pdata <- pdata %>%
  mutate(Tumor_Barcode2 = factor(Tumor_Barcode2, levels = tmplev))

library(ggforce)
pdata %>%
  ggplot() +
  geom_link(
    data = pdata %>% filter(WGD_locate == t),
    aes(x = t2, xend = t2, y = 0, yend = mean_bootstrap),
    size = 4,
    lineend = "round",
    color = ncicolpal[1]
  ) +
  geom_link(
    data = pdata %>% filter(WGD_locate != t),
    aes(x = t2, xend = t2, y = 0, yend = mean_bootstrap),
    size = 2.5,
    lineend = "round",
    color = "black"
  ) +
  geom_link(
    aes(
      x = t2,
      xend = t2,
      y = 0,
      yend = mean_bootstrap,
      color = mean_bootstrap
    ),
    size = 2,
    lineend = "round"
  ) +
  geom_errorbar(
    aes(x = t2, ymin = lower_ci, ymax = upper_ci),
    width = 0.2,
    size = 0.5,
    lineend = "round"
  ) +
  geom_text(aes(t2, y = 1.05, label = label), hjust = 1.5) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = pretty_breaks()) +
  scale_color_viridis_c(direction = -1) +
  labs(x = NULL, y = NULL, col = 'Time') +
  coord_polar(theta = "y") +
  facet_wrap(~Tumor_Barcode2) +
  #facet_wrap(~Subtype+Tumor_Barcode2)+
  theme_ipsum_rc(grid = "XY", axis = F, axis_title_just = 'm', base_size = 10) +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_text(hjust = 0.5, size = 11),
    panel.spacing = unit(0.1, 'cm'),
    panel.background = element_rect(fill = '#984ea3', linewidth = 0),
    panel.border = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

ggsave(
  filename = 'gain12_WGD_timing_G.pdf',
  width = 8,
  height = 8,
  device = cairo_pdf
)


# Pivot to long format
## W first

pdata <- segment_time_all2 %>%
  # pivot all columns that start with t_
  pivot_longer(
    cols = starts_with("t_"),
    names_to = c("t", "measure"),
    names_pattern = "t_(\\d+)_?(.*)", # split t_1, t_1_mean_bootstrap etc
    values_to = "value"
  ) %>%
  # fill missing measure with "t_value"
  mutate(measure = if_else(measure == "", "t_value", measure)) %>%
  # pivot wider to get t_value, mean_bootstrap, lower_ci, upper_ci in separate columns
  pivot_wider(
    names_from = measure,
    values_from = value
  ) %>%
  mutate(label = paste0('T', t)) %>%
  filter(!is.na(mean_bootstrap)) %>%
  mutate(t = as.integer(t)) %>%
  mutate(t2 = t) %>%
  mutate(mean_bootstrap = if_else(mean_bootstrap > 1, 1, mean_bootstrap)) %>%
  mutate(lower_ci = if_else(lower_ci > 1, 1, lower_ci)) %>%
  mutate(upper_ci = if_else(upper_ci > 1, 1.01, upper_ci)) %>%
  mutate(
    Tumor_Barcode2 = paste0(
      str_sub(Tumor_Barcode, 1, 12),
      ', nMut=',
      num_mutations_used,
      "\nCNV=",
      highest_copy_number,
      ", Event=",
      event_order
    )
  ) %>%
  mutate(WGD_locate = str_locate(event_order, "W")[, 1]) %>%
  filter(
    str_length(event_order) > 1,
    str_detect(event_order, 'W'),
    str_detect(event_order, 'G')
  ) %>%
  filter(str_starts(event_order, 'W'))


tmplev <- pdata %>%
  group_by(Subtype, Tumor_Barcode, Tumor_Barcode2) %>%
  summarise(n = n(), maxvalue = max(mean_bootstrap)) %>%
  ungroup() %>%
  arrange(Subtype, desc(n), desc(maxvalue)) %>%
  pull(Tumor_Barcode2)

pdata <- pdata %>%
  mutate(Tumor_Barcode2 = factor(Tumor_Barcode2, levels = tmplev))

library(ggforce)
pdata %>%
  #filter(Tumor_Barcode %in% c("TCGA-2G-AAFC-01A-11D-A734-36", "TCGA-X3-A8G4-01A-11D-A749-36", "TCGA-XY-A8S3-01B-11D-A749-36", 'TCGA-2G-AAFC-01A-11D-A734-36')) %>%
  #filter(Tumor_Barcode %in% tmp) %>%
  ggplot() +
  geom_link(
    data = pdata %>% filter(WGD_locate == t),
    aes(x = t2, xend = t2, y = 0, yend = mean_bootstrap),
    size = 4,
    lineend = "round",
    color = ncicolpal[1]
  ) +
  geom_link(
    data = pdata %>% filter(WGD_locate != t),
    aes(x = t2, xend = t2, y = 0, yend = mean_bootstrap),
    size = 2.5,
    lineend = "round",
    color = "black"
  ) +
  geom_link(
    aes(
      x = t2,
      xend = t2,
      y = 0,
      yend = mean_bootstrap,
      color = mean_bootstrap
    ),
    size = 2,
    lineend = "round"
  ) +
  geom_errorbar(
    aes(x = t2, ymin = lower_ci, ymax = upper_ci),
    width = 0.2,
    size = 0.5,
    lineend = "round"
  ) +
  geom_text(aes(t2, y = 1.05, label = label), hjust = 1.5) +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = pretty_breaks()) +
  scale_color_viridis_c(direction = -1) +
  labs(x = NULL, y = NULL, col = 'Time') +
  coord_polar(theta = "y") +
  #facet_wrap(~Tumor_Barcode2,nrow = 4)+
  facet_wrap(~ Subtype + Tumor_Barcode2, nrow = 4) +
  theme_ipsum_rc(grid = "XY", axis = F, axis_title_just = 'm', base_size = 12) +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_text(hjust = 0.5, size = 14),
    panel.spacing = unit(0.1, 'cm'),
    panel.background = element_rect(fill = '#984ea3', linewidth = 0),
    panel.border = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

ggsave(
  filename = 'gain12_WGD_timing_W.pdf',
  width = 16,
  height = 16,
  device = cairo_pdf
)
