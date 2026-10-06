# Data-driven oncoplot drawing extracted from scripts/manuscript-figures/figure2/helpers/generate_astro_focused_oncoplot.R
feature_layout <- bind_rows(
  tibble(Type = "Feature", Gene = "Evo_Group"),
  tibble(Type = "Mutation_Driver", Gene = c(target_genes, "TERT")),
  tibble(Type = "Allele_Specific_CN", Gene = allele_cn_genes)
) %>%
  mutate(
    Feature_ID = paste(.data$Type, .data$Gene, sep = "__"),
    row_height = if_else(.data$Type == "Feature", 0.45, 1),
    row_ymax = cumsum(.data$row_height),
    row_ymin = .data$row_ymax - .data$row_height,
    y = .data$row_ymin + .data$row_height / 2,
    tile_height = .data$row_height * 0.96,
    Feature_Label = case_when(
      .data$Type == "Feature" ~ "",
      .data$Type == "Mutation_Driver" & .data$Gene == "TERT" ~ "TERT promoter",
      .data$Type == "Mutation_Driver" ~ paste(.data$Gene, "mutation"),
      .data$Type == "Allele_Specific_CN" ~ paste(.data$Gene, "copy number"),
      TRUE ~ .data$Gene
    )
  )
x_limits <- c(0.5, length(sample_level) + 0.5)
y_limits <- c(max(feature_layout$row_ymax), -0.55)
feature_boundaries <- feature_layout$row_ymax[-nrow(feature_layout)]

sample_index <- tibble(
  Tumor_Barcode = sample_level,
  x = seq_along(sample_level)
)
group_boundary <- wgs_sub %>%
  filter(.data$Tumor_Barcode %in% sample_level) %>%
  mutate(Tumor_Barcode = factor(.data$Tumor_Barcode, levels = sample_level)) %>%
  arrange(.data$Tumor_Barcode) %>%
  count(.data$Evo_Group, name = "n") %>%
  mutate(boundary = cumsum(.data$n)) %>%
  filter(.data$Evo_Group == "ASTRO_Group1") %>%
  pull(.data$boundary)
group_annotation <- sample_index %>%
  left_join(
    wgs_sub %>% select("Tumor_Barcode", "Evo_Group"),
    by = "Tumor_Barcode"
  ) %>%
  group_by(.data$Evo_Group) %>%
  summarise(
    xmin = min(.data$x) - 0.47,
    xmax = max(.data$x) + 0.47,
    x = mean(.data$x),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(group_display_labels[as.character(.data$Evo_Group)], " (n=", .data$n, ")")
  )

oncoplot_data <- bind_rows(data_mutation, data_promoter, data_pl_cn, data_feature) %>%
  filter(.data$Tumor_Barcode %in% sample_level) %>%
  mutate(Feature_ID = paste(.data$Type, .data$Gene, sep = "__")) %>%
  inner_join(feature_layout, by = c("Type", "Gene", "Feature_ID")) %>%
  inner_join(sample_index, by = "Tumor_Barcode") %>%
  group_by(.data$Tumor_Barcode, .data$Feature_ID, .data$x, .data$y, .data$tile_height) %>%
  summarise(
    Alteration = paste(unique(.data$Alteration), collapse = "/"),
    .groups = "drop"
  ) %>%
  separate_rows(Alteration, sep = "/") %>%
  filter(!is.na(.data$Alteration), .data$Alteration != "NA") %>%
  group_by(.data$Tumor_Barcode, .data$Feature_ID, .data$x, .data$y, .data$tile_height) %>%
  mutate(
    n_segments = n(),
    segment = row_number(),
    xmin = .data$x - 0.475,
    xmax = .data$x + 0.475,
    ymin = .data$y - .data$tile_height / 2 + (.data$segment - 1) * .data$tile_height / .data$n_segments,
    ymax = .data$y - .data$tile_height / 2 + .data$segment * .data$tile_height / .data$n_segments
  ) %>%
  ungroup()

oncoplot_background <- crossing(
  Tumor_Barcode = sample_level,
  Feature_ID = feature_layout$Feature_ID
) %>%
  left_join(sample_index, by = "Tumor_Barcode") %>%
  left_join(feature_layout %>% select(Feature_ID, y, tile_height), by = "Feature_ID")

oncoplot_main <- ggplot(oncoplot_background, aes(x = .data$x, y = .data$y)) +
  geom_tile(aes(height = .data$tile_height), width = 0.95, fill = "gray95") +
  geom_rect(
    data = oncoplot_data,
    aes(
      xmin = .data$xmin,
      xmax = .data$xmax,
      ymin = .data$ymin,
      ymax = .data$ymax,
      fill = .data$Alteration
    ),
    inherit.aes = FALSE,
    linewidth = 0
  ) +
  geom_vline(xintercept = seq(1.5, length(sample_level) - 0.5, by = 1), color = "white", linewidth = 0.05) +
  geom_vline(xintercept = group_boundary + 0.5, color = "grey25", linewidth = 0.42) +
  geom_segment(
    data = group_annotation,
    aes(x = .data$xmin, xend = .data$xmax, y = -0.08, yend = -0.08),
    inherit.aes = FALSE,
    color = "grey45",
    linewidth = 0.32
  ) +
  geom_text(
    data = group_annotation,
    aes(x = .data$x, y = -0.34, label = .data$label),
    inherit.aes = FALSE,
    family = base_family,
    fontface = "bold",
    color = "grey15",
    size = 2.45
  ) +
  geom_hline(yintercept = feature_boundaries, color = "white", linewidth = 0.18) +
  scale_fill_manual(values = landscape_colors, breaks = sort(unique(oncoplot_data$Alteration))) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_reverse(
    limits = y_limits,
    breaks = feature_layout$y,
    labels = feature_layout$Feature_Label,
    expand = c(0, 0)
  ) +
  hrbrthemes::theme_ipsum_rc(base_family = base_family, base_size = 8, axis = FALSE, ticks = FALSE, grid = FALSE) +
  theme(
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.3),
    plot.margin = margin(t = 0.5, r = 0, b = 0, l = 0, unit = "pt")
  )

