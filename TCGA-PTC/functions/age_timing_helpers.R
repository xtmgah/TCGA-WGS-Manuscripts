# Shared, minimal support for Fig5, Fig6 and ExtendedDataFig4.
# Historical definitions: TP_Comparision.R and Sherlock_functions.R.
# Input objects remain in a private input directory; no .Rprofile is sourced.
ptc_age_packages <- function(extra = character()) {
  pkgs <- unique(c('dplyr','tidyr','forcats','ggplot2','scales','broom','ggrepel', extra))
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Missing R packages: ', paste(missing, collapse = ', '))
  for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  ptc_font()
}
ptc_age_load <- function(input_dir, file, object) {
  path <- file.path(input_dir, file)
  if (!file.exists(path)) stop('Required input missing: ', file)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object, env, inherits = FALSE)) stop('Object ', object, ' missing in ', file)
  env[[object]]
}
ptc_age_common <- function(input_dir) {
  list(
    cohort = ptc_age_load(input_dir, 'wgsdata_tp.RData', 'wgsdata'),
    covariates = ptc_age_load(input_dir, 'wgs_covdata_tp.RData', 'wgs_covdata'),
    clinical = ptc_age_load(input_dir, 'wgs_clinical.RData', 'wgs_clinical'),
    alterations = ptc_age_load(input_dir, 'Genome_landscape_manual_final.RData', 'data_top0'),
    colors = ptc_age_load(input_dir, 'wgsdata_tp.RData', 'study_color2'))
}
ptc_age_result <- function(id, expression) {
  tryCatch({
    value <- force(expression)
    c(list(panel_id = id, reproducibility_status = 'generated; validation pending'), value)
  }, error = function(e) list(panel_id = id,
    reproducibility_status = 'code identified but reproduction blocked',
    notes = conditionMessage(e)))
}
ptc_age_exposure <- function(data, clinical, replace_clinical_na = FALSE) {
  exposure <- clinical %>% dplyr::select(Subject, age_at_exposure, exposure_dose)
  # Original latency blocks first replace NA in matched clinical rows by zero;
  # the logistic block does not. Unmatched subjects then receive the final defaults.
  if (replace_clinical_na)
    exposure <- exposure %>% mutate(across(-Subject, ~ replace_na(.x, 0)))
  data %>% left_join(exposure, by = 'Subject') %>%
    mutate(exposure_dose = pmin(replace_na(exposure_dose, 0), 1000),
      age_at_exposure = replace_na(age_at_exposure, 100))
}
ptc_age_mrca <- function(input_dir, cohort) {
  ptc_age_load(input_dir, 'thyroid_evolution_analysis/Chronological_timing_short.RData', 'MRCAdata') %>%
    filter(acceleration == '1x') %>%
    right_join(cohort %>% dplyr::select(Study2, Tumor_Barcode), by = 'Tumor_Barcode') %>%
    filter(!is.na(Latency))
}
ptc_age_mutant_group <- function(data, alterations, order_by) {
  braf <- alterations %>% filter(Gene == 'BRAF') %>% pull(Tumor_Barcode) %>% unique()
  tert <- alterations %>% filter(Gene == 'TERT') %>% pull(Tumor_Barcode) %>% unique()
  data %>% mutate(Group = case_when(
    Tumor_Barcode %in% braf & Tumor_Barcode %in% tert ~ 'BRAF + TERT',
    Tumor_Barcode %in% braf ~ 'BRAF', Tumor_Barcode %in% tert ~ 'TERT', TRUE ~ 'Others')) %>%
    mutate(Group = fct_reorder(Group, .data[[order_by]]))
}
ptc_age_box <- function(data, value, ylab, positions) {
  # Extracted statistical logic of historical plot_group_compare(): two-sided
  # pairwise Wilcoxon tests with BH correction, quasirandom dots and seed=1.
  data <- data %>% filter(!is.na(.data[[value]]), !is.na(Group))
  tests <- rstatix::pairwise_wilcox_test(data, as.formula(paste(value, '~ Group')),
    p.adjust.method = 'BH') %>% rstatix::add_xy_position(x = 'Group')
  if (length(positions) == nrow(tests)) tests$y.position <- positions
  else tests$y.position <- tests$y.position + positions
  tests$p.label <- ifelse(tests$p.adj < .001, sprintf('FDR = %.2e', tests$p.adj),
    sprintf('FDR = %.3f', tests$p.adj))
  tests$col <- ifelse(tests$p.adj < .05, '#BB0E3D', 'black')
  set.seed(1)
  plot <- ggplot(data, aes(Group, .data[[value]], fill = Group)) +
    ggbeeswarm::geom_quasirandom(pch = 21, size = 2, width = .4, color = 'white', stroke = .2, alpha = 1) +
    geom_boxplot(width = .5, fill = NA, color = 'gray20', outlier.shape = NA, linewidth = .6) +
    ggpubr::stat_pvalue_manual(tests, label = 'p.label', y.position = 'y.position',
      xmin = 'xmin', xmax = 'xmax', color = 'col', tip.length = .01,
      step.increase = if (nrow(tests) > 1) .05 else 0, hide.ns = TRUE, size = 4.5) +
    scale_fill_brewer(palette = 'Set2') + scale_color_identity() +
    scale_y_continuous(breaks = breaks_pretty(n = 7)) +
    labs(x = NULL, y = ylab, fill = 'Group') + ptc_theme() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    coord_cartesian(clip = 'off')
  list(plot = plot, tests = tests)
}
