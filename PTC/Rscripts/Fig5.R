# Figure 5: age, exposure and driver-event associations.
run_Fig5 <- function(input_dir, output_dir) {
  ptc_age_packages(c('ggsci'))
  x <- ptc_age_common(input_dir)
  cov <- x$covariates %>% left_join(x$cohort %>% dplyr::select(Tumor_Barcode, Study2), by = 'Tumor_Barcode')
  mrca <- ptc_age_mrca(input_dir, x$cohort)
  results <- list()

  results[['Fig. 5a']] <- ptc_age_result('Fig. 5a', {
    # Only tumors with at least one curated alteration enter the model matrix.
    alteration_matrix <- x$alterations %>% filter(Tumor_Barcode %in% x$cohort$Tumor_Barcode) %>%
      transmute(Tumor_Barcode, name = paste(Gene, Type, sep = '|')) %>% distinct() %>%
      mutate(value = TRUE) %>% pivot_wider(names_from = name, values_from = value, values_fill = FALSE) %>%
      pivot_longer(-Tumor_Barcode)
    eligible <- alteration_matrix %>% filter(value) %>% count(name) %>% filter(n >= 10) %>% pull(name)
    dat <- alteration_matrix %>% filter(name %in% eligible) %>% left_join(cov, by = 'Tumor_Barcode') %>%
      ptc_age_exposure(x$clinical)
    fits <- dat %>% group_by(name) %>% group_modify(~ broom::tidy(glm(
      value ~ age_at_diagnosis + Tumor_Purity + Sex + PC1 + PC2 + exposure_dose + age_at_exposure,
      data = .x, family = binomial()))) %>% ungroup() %>%
      filter(term %in% c('exposure_dose','age_at_diagnosis','age_at_exposure')) %>%
      mutate(effect_type = term) %>% group_by(effect_type) %>% mutate(p_adj = p.adjust(p.value, 'BH')) %>%
      ungroup() %>% separate(name, c('Gene','Type'), sep = '\\|') %>%
      mutate(label = if_else(p_adj < .2, Gene, NA_character_))
    p <- ggplot(fits, aes(estimate, -log10(p_adj), fill = effect_type, shape = Type)) +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray50') +
      geom_hline(yintercept = -log10(.05), linetype = 'dotted', color = '#BB0E3D') +
      geom_point(size = 3.5, color = 'black', stroke = .25) +
      ggrepel::geom_text_repel(aes(label = label), size = 4.5, na.rm = TRUE, seed = 1) +
      scale_shape_manual(values = c(22,21,23)) + ggsci::scale_fill_aaas() +
      labs(x = 'Effect size (log-odds estimate)', y = expression(-log[10]('BH-adjusted P')), fill = 'Effect type') +
      ptc_theme() + guides(fill = guide_legend(override.aes = list(shape = 21)))
    list(files = ptc_save(p, output_dir, 'Fig5a', 10, 6), statistics = fits,
      notes = 'Model matrix includes tumors with at least one curated alteration. Dose is capped at 1000; missing age at exposure is 100. BH correction is applied within predictor.')
  })

  frequency_specs <- list('Fig. 5b' = c('TERT','all','TERT promoter mutation frequency','Fig5b'),
    'Fig. 5c' = c('BRAF','Mutation_Driver','BRAF mutation frequency','Fig5c'),
    'Fig. 5d' = c('all','Fusion','Fusion frequency','Fig5d'))
  for (id in names(frequency_specs)) {
    spec <- frequency_specs[[id]]
    results[[id]] <- ptc_age_result(id, {
      selected <- x$alterations
      if (spec[1] != 'all') selected <- selected %>% filter(Gene == spec[1])
      if (spec[2] != 'all') selected <- selected %>% filter(Type == spec[2])
      ids <- unique(selected$Tumor_Barcode)
      counts <- cov %>% count(Study2, age_cat, name = 'total') %>%
        left_join(cov %>% filter(Tumor_Barcode %in% ids) %>% count(Study2, age_cat), by = c('Study2','age_cat')) %>%
        filter(!is.na(age_cat)) %>% mutate(Freq = replace_na(n / total, 0))
      p <- ggplot(counts, aes(age_cat, Freq, color = Study2)) +
        geom_point(aes(size = total, fill = Study2), pch = 21) + geom_line(aes(group = Study2)) +
        scale_color_manual(values = x$colors) + scale_fill_manual(values = x$colors) +
        scale_y_continuous(labels = label_percent()) + scale_size_binned(breaks = breaks_pretty()) +
        labs(x = 'Age at diagnosis (years)', y = spec[3], size = 'Sample size', color = 'Group', fill = 'Group') + ptc_theme()
      list(files = ptc_save(p, output_dir, spec[4], 10, 6), statistics = counts,
        notes = 'Age 70 is included in 60–70. Age groups without altered tumors have frequency zero.')
    })
  }

  results[['Fig. 5e']] <- ptc_age_result('Fig. 5e', {
    # Join on tumor barcode and diagnosis age to align the timing and covariate rows.
    dat <- mrca %>% left_join(x$covariates, by = c('Tumor_Barcode','age_at_diagnosis')) %>%
      ptc_age_exposure(x$clinical, replace_clinical_na = TRUE)
    fit <- lm(Latency ~ Tumor_Purity + Sex + exposure_dose + age_at_exposure + PC1 + PC2, data = dat)
    stats <- broom::tidy(fit, conf.int = TRUE) %>% filter(term != '(Intercept)') %>%
      mutate(term = recode(term, Tumor_Purity = 'Tumor purity', SexFemale = 'Female vs male (ref)',
        exposure_dose = 'Exposure dose', age_at_exposure = 'Age at exposure', PC1 = 'Ancestry PC1', PC2 = 'Ancestry PC2'),
        term = fct_reorder(term, estimate), label = paste0('β = ', round(estimate,2), ', P = ', label_scientific(digits = 3)(p.value)))
    p <- ggplot(stats, aes(estimate, term, xmin = conf.low, xmax = conf.high, color = p.value < .05)) +
      geom_point(size = 3) + geom_errorbar(orientation = 'y', width = .2) +
      ggrepel::geom_text_repel(aes(label = label), size = 4.5, nudge_y = .2, seed = 1) +
      scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = '#BB0E3D')) +
      labs(x = 'Regression coefficient', y = NULL) + guides(color = 'none') + ptc_theme()
    list(files = ptc_save(p, output_dir, 'Fig5e', 10, 6), statistics = stats,
      notes = 'Uses the 1x acceleration MRCA intermediate. Diagnosis age is not included in this latency model.')
  })

  results[['Fig. 5f']] <- ptc_age_result('Fig. 5f', {
    alts <- bind_rows(x$alterations, x$alterations %>% filter(Type == 'Fusion') %>% mutate(Gene = 'All_Gene') %>% distinct())
    eligible <- alts %>% mutate(Key = paste(Gene,Type,sep='@')) %>% count(Key) %>% filter(n > 10) %>% pull(Key)
    pairs <- alts %>% mutate(Key = paste(Gene,Type,sep='@')) %>% dplyr::select(Tumor_Barcode,Key) %>%
      filter(Key %in% eligible) %>% distinct() %>% mutate(value='Yes')
    # Preserve all cohort members and all eligible alteration categories.
    dat <- x$cohort %>% dplyr::select(Tumor_Barcode,Study2) %>% left_join(pairs, by='Tumor_Barcode') %>%
      pivot_wider(names_from=Key, values_from=value, values_fill='No') %>% dplyr::select(-any_of('NA')) %>%
      pivot_longer(-c(Tumor_Barcode,Study2), names_to='name') %>% mutate(value=factor(value,levels=c('No','Yes'))) %>%
      left_join(mrca, by=c('Tumor_Barcode','Study2')) %>%
      left_join(x$covariates, by=c('Tumor_Barcode','age_at_diagnosis')) %>%
      ptc_age_exposure(x$clinical, replace_clinical_na=TRUE)
    stats <- dat %>% group_by(name) %>% group_modify(~ {
      fit <- tryCatch(lm(Latency ~ Tumor_Purity + Sex + exposure_dose + age_at_exposure + PC1 + PC2 + value, data=.x), error=function(e) NULL)
      if(is.null(fit)) tibble() else broom::tidy(fit)
    }) %>% ungroup() %>% filter(term=='valueYes') %>% arrange(p.value) %>%
      mutate(FDR=p.adjust(p.value,method='holm')) %>% separate(name,c('Gene','Type'),sep='@')
    p <- ggplot(stats,aes(estimate,-log10(FDR),fill=Type)) +
      geom_point(pch=21,size=3,stroke=.2) +
      ggrepel::geom_text_repel(data=stats %>% filter(p.value<.05),aes(label=Gene),size=4.5,seed=1) +
      geom_vline(xintercept=0,linewidth=.2) + geom_hline(yintercept=-log10(.05),linetype='dashed',color='#ff7f00') +
      ggsci::scale_fill_d3() + labs(x='Regression coefficient',y=expression(-log[10]('Holm-adjusted P')),fill='Alteration') +
      ptc_theme() + coord_cartesian(clip='off')
    list(files=ptc_save(p,output_dir,'Fig5f',10,6),statistics=stats,
      notes='Holm-adjusted P values. Alteration categories require more than 10 event rows before deduplication.')
  })
  results[['Fig. 5g']] <- list(panel_id='Fig. 5g', status='unavailable',
    notes='The MRCA-age regression and plotting implementation for panel g is not available in this release.')
  results
}
