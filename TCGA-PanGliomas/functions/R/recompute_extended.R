# Executed by recompute.R in the materialized, disposable workspace only.
# Reference tables supply static labels/order, never numerical test results.
coverage <- list(); permutations <- list()
regenerate <- function(path, actual, keys, method, family, inputs) {
  expected <- rd(path); name <- basename(path)
  stopifnot(setequal(names(actual), names(expected)), nrow(actual)==nrow(expected))
  key <- function(x) do.call(paste, c(lapply(x[keys], as.character), sep='\r'))
  ak <- key(actual); ek <- key(expected)
  stopifnot(!anyDuplicated(ak), !anyDuplicated(ek), setequal(ak,ek))
  actual <- actual[match(ek,ak), names(expected), drop=FALSE]
  for (v in names(expected)) {
    a <- actual[[v]]; e <- expected[[v]]
    if (is.numeric(e)) {
      good <- numeric_agreement(a,e,probability_column(v))
      delta <- abs(a-e); delta <- delta[is.finite(delta)]
      checks[[paste(name,v,sep=':')]] <<- list(pass=good, comparisons=length(e),
        maximum_absolute_error=if(length(delta)) max(delta) else 0,
        relative_tolerance=1e-8, absolute_tolerance=if(probability_column(v)) 0 else 1e-10)
    } else good <- identical(as.character(a),as.character(e))
    if (!good) stop(paste('Scientific mismatch:',name,v))
  }
  write_tsv(actual,path)
  write_tsv(actual,file.path(out,paste0('recomputed_',name)))
  coverage[[name]] <<- list(rows=nrow(actual), method=method, correction_family=family,
    measured_inputs=inputs, verified_against='immutable packaged reference table',
    rendering_input=path, status='recomputed_and_verified')
  actual
}
F1 <- 'results/analysis/pooled_814_rerun_2026-09-28/figures/figure1/derived-data'
F2 <- 'results/figures/figure2/derived-data'
F3 <- 'results/figures/figure3/derived-data'
F4 <- file.path(B,'figures/figure4/derived-data')

# Figure 1: three pairs within each metric, retaining the manuscript bracket order.
v <- rd(file.path(F1,'panel_h_chronological_timing_values.tsv'))
a <- pairwise_metrics(v,'Metric','Subtype_Final',list(c('ASTRO','OLIGO'),c('OLIGO','GBM'),c('ASTRO','GBM')),
                      c('MRCA age','Age at diagnosis','Latency'))
a$label <- ifelse(a$padj_bh<.001,'FDR<.001',sprintf('FDR=%.2f',a$padj_bh))
regenerate(file.path(F1,'panel_h_chronological_pairwise_tests.tsv'),a,c('metric','group_1','group_2'),
           'Two-sided Wilcoxon rank sum, exact=FALSE, continuity correction','BH: 3 pairs per metric',
           'panel_h_chronological_timing_values.tsv; 1x chronology')

# Figure 2: latency is part of the three-test family even though only two panels show it.
v <- rd(file.path(F2,'astro_age_mrca_latency_primary_2p5x_sample_table.tsv'))
stopifnot(!anyDuplicated(v$Subject),nrow(v)==254,all(v$acceleration=='2.5x'))
ms <- c(age_at_diagnosis_years='Age at diagnosis',age_mrca_years='MRCA age',latency_years='Latency')
a <- bind_rows(lapply(names(ms),function(m) {
  x <- v[[m]][v$Evo_Group=='ASTRO_Group1' & is.finite(v[[m]])]
  y <- v[[m]][v$Evo_Group=='ASTRO_Group2' & is.finite(v[[m]])]
  tibble(metric=m,label=unname(ms[m]),acceleration='2.5x',n_group1=length(x),n_group2=length(y),
    median_group1=median(x),q1_group1=unname(quantile(x,.25)),q3_group1=unname(quantile(x,.75)),mean_group1=mean(x),
    median_group2=median(y),q1_group2=unname(quantile(y,.25)),q3_group2=unname(quantile(y,.75)),mean_group2=mean(y),
    median_difference_group2_minus_group1=median(y)-median(x),mean_difference_group2_minus_group1=mean(y)-mean(x),
    wilcox_pvalue=wilcox.test(x,y,exact=FALSE)$p.value,welch_pvalue=t.test(x,y)$p.value,
    higher_median_group=if(median(x)==median(y))'tie' else if(median(x)>median(y))'ASTRO_Group1' else 'ASTRO_Group2')
})) |> mutate(q_wilcox_primary=p.adjust(wilcox_pvalue,'BH'),q_welch_primary=p.adjust(welch_pvalue,'BH'))
regenerate(file.path(F2,'astro_age_mrca_latency_primary_tests.tsv'),a,'metric',
           'Two-sided Wilcoxon and Welch; type-7 quartiles','BH separately across all 3 metrics for each test',
           'astro_age_mrca_latency_primary_2p5x_sample_table.tsv; finite observations per metric')
v <- rd(file.path(F2,'diagnostics/tp53_copy_number_analysis_data.tsv'))
stopifnot(nrow(v)==226,!anyDuplicated(v$Subject))
v$Evo_Group <- factor(v$Evo_Group,levels=c('ASTRO_Group1','ASTRO_Group2'))
v$total_cn_centered_at_2 <- v$total_cn_raw-2
model <- lm(mutant_copies_raw ~ total_cn_centered_at_2 * Evo_Group,v)
a <- tidy(model,conf.int=TRUE) |> rename(std_error=std.error,p_value=p.value,conf_low=conf.low,conf_high=conf.high)
regenerate(file.path(F2,'diagnostics/tp53_copy_number_model_terms.tsv'),a,'term',
           'OLS interaction model on uncapped total and mutant copy number; 95% t intervals','None',
           'tp53_copy_number_analysis_data.tsv; 226 unique primary patients')

# Figures 2d and 4e: populate displayed contrasts from the newly fitted Cox models.
for (subtype in c('ASTRO','GBM')) {
  path <- file.path(S,if(subtype=='ASTRO')'Figure2d_display_contrasts.tsv' else 'Figure4e_display_contrasts.tsv')
  a <- rd(path); fit <- tidy(fits[[paste0(subtype,'_primary')]],exponentiate=TRUE,conf.int=TRUE)
  term_map <- if(subtype=='ASTRO')c(Evo_GroupASTRO_Group2='groupCIC/TERT',PC1_scaled='PC1_z',PC2_scaled='PC2_z',sexMale='sexMale',tumor_gradeG3='gradeG3',age_10yr='age10',BB_Purity_10pct='purity10')else setNames(fit$term,fit$term)
  ix <- match(unname(term_map[a$term]),fit$term); stopifnot(!anyNA(ix))
  for (field in intersect(c('estimate','std.error','statistic','p.value','conf.low','conf.high'),names(a))) a[[field]] <- fit[[field]][ix]
  for (field in c('estimate','conf_low','conf_high','statistic')) {
    actual_field <- switch(field,conf_low='conf.low',conf_high='conf.high',field)
    if(paste0('original_',field) %in% names(a)) a[[paste0('original_',field)]] <- a[[actual_field]]
  }
  rev <- a$contrast_reversed
  a$estimate[rev] <- 1/fit$estimate[ix][rev]
  a$conf.low[rev] <- 1/fit$conf.high[ix][rev]; a$conf.high[rev] <- 1/fit$conf.low[ix][rev]
  if('statistic' %in% names(a)) a$statistic[rev] <- -fit$statistic[ix][rev]
  if('sig' %in% names(a)) a$sig <- ifelse(a$p.value<.05,'p < 0.05','NS')
  regenerate(path,a,'term','Efron Cox fit; reciprocal HR, swapped reciprocal CI and negated Z for reversed contrast',
             'None',paste0(subtype,'_cox_data.tsv'))
}

# Figure 3: calculate the full 28-test inventory, its BH family, and the separate promoter test.
metrics <- rd(file.path(F3,'figure3_sample_metrics.tsv'))
stopifnot(nrow(metrics)==263,!anyDuplicated(metrics$Tumor_Barcode))
calc_test <- function(variable,metric,binary) {
  value <- metrics[[variable]]; x <- value[metrics$Evo_Group=='ASTRO_Group1' & !is.na(value)]
  y <- value[metrics$Evo_Group=='ASTRO_Group2' & !is.na(value)]
  if(binary) {
    cells <- c(sum(x),sum(!x),sum(y),sum(!y)); adjusted <- cells + if(any(cells==0)) .5 else 0
    logor <- log((adjusted[1]/adjusted[2])/(adjusted[3]/adjusted[4])); se <- sqrt(sum(1/adjusted))
    tibble(test_family='binary_fisher',variable=variable,metric=metric,n_group1=length(x),n_group2=length(y),
      value_group1=sum(x),value_group2=sum(y),pct_group1=100*mean(x),pct_group2=100*mean(y),
      effect_label='OR 17p-Amp vs CIC/TERT',effect=exp(logor),ci_low=exp(logor-1.96*se),ci_high=exp(logor+1.96*se),
      pvalue=fisher.test(matrix(cells,nrow=2,byrow=TRUE))$p.value)
  } else {
    tibble(test_family='continuous_wilcoxon',variable=variable,metric=metric,n_group1=length(x),n_group2=length(y),
      value_group1=median(x),value_group2=median(y),mean_group1=mean(x),mean_group2=mean(y),
      q25_group1=unname(quantile(x,.25)),q75_group1=unname(quantile(x,.75)),
      q25_group2=unname(quantile(y,.25)),q75_group2=unname(quantile(y,.75)),
      effect_label='Median difference 17p-Amp - CIC/TERT',effect=median(x)-median(y),
      ci_low=NA_real_,ci_high=NA_real_,pvalue=wilcox.test(x,y,exact=FALSE)$p.value)
  }
}
spec <- rd(file.path(F3,'figure3_tests.tsv')) |> select(variable,metric,test_family)
stopifnot(nrow(spec)==28,sum(spec$test_family=='binary_fisher')==2)
a <- bind_rows(lapply(seq_len(nrow(spec)),function(i)calc_test(spec$variable[i],spec$metric[i],spec$test_family[i]=='binary_fisher')))
fmt_p <- function(p) ifelse(p<.001,'p<0.001',paste0('p=',sub('0+$','',sub('\\.$','',ifelse(p<.01,sprintf('%.3f',p),sprintf('%.2f',p))))))
a <- a |> mutate(padj_bh_all_tests=p.adjust(pvalue,'BH'),p_label=fmt_p(pvalue),fdr_label=paste0('FDR=',signif(padj_bh_all_tests,2)))
regenerate(file.path(F3,'figure3_tests.tsv'),a,'variable','Fisher exact or two-sided Wilcoxon; sample OR with 0.5 correction if any zero cell',
           'BH across all 28 tests','figure3_sample_metrics.tsv')
a <- calc_test('TERT_promoter_hotspot','TERT promoter hotspot mutation',TRUE)
regenerate(file.path(F3,'figure3_tert_promoter_test.tsv'),a,'variable','Fisher exact; sample OR and log-Wald interval','None; separate test',
           'figure3_sample_metrics.tsv')
spec <- rd(file.path(F3,'figure3_specificity_fingerprint.tsv')) |> select(variable,metric,metric_class,value_transform)
stopifnot(nrow(spec)==12,all(spec$value_transform %in% c('identity','log10p1')))
a <- bind_rows(lapply(seq_len(nrow(spec)),function(i) {
  v <- metrics[[spec$variable[i]]]; if(spec$value_transform[i]=='log10p1') v <- log10(v+1)
  x <- v[metrics$Evo_Group=='ASTRO_Group1' & is.finite(v)]; y <- v[metrics$Evo_Group=='ASTRO_Group2' & is.finite(v)]
  sd <- sqrt(((length(x)-1)*var(x)+(length(y)-1)*var(y))/(length(x)+length(y)-2))
  effect <- if(is.finite(sd) && sd>0)(mean(x)-mean(y))/sd else 0
  pv <- wilcox.test(x,y,exact=FALSE)$p.value
  bind_cols(spec[i,],tibble(n_group1=length(x),n_group2=length(y),mean_group1=mean(x),mean_group2=mean(y),
    median_group1=median(x),median_group2=median(y),standardized_mean_difference=effect,pvalue=pv,
    significant=pv<.05,significance=if(pv<.05)'p < 0.05' else 'NS',
    enrichment_call=if(pv>=.05 || effect==0)'Not significant' else if(effect>0)'Group1 higher' else 'Group2 higher'))
}))
a <- regenerate(file.path(F3,'figure3_specificity_fingerprint.tsv'),a,'variable','Pooled-SD standardized difference and two-sided Wilcoxon',
                '12-feature BH for panel 3h; nominal labels retained in historical table','figure3_sample_metrics.tsv; recorded transforms')
a$adjusted_pvalue <- p.adjust(a$pvalue,'BH'); write_tsv(a,file.path(out,'figure3_h_BH_adjusted.tsv'))

# Figure 4: all six metric families, including nondisplayed MRCA/age comparisons.
v <- rd(file.path(F4,'publication-panels/figure4_six_metric_values.tsv'))
a <- pairwise_metrics(v,'metric','DN_Group',combn(c('DN1','DN2','DN3'),2,simplify=FALSE),
                      c('MRCA age','Age at diagnosis','Latency','MATH','PGA','Ploidy')) |>
  select(metric,group1=group_1,group2=group_2,n_group1=n_group_1,n_group2=n_group_2,p_value=pvalue,q_value=padj_bh)
regenerate(file.path(F4,'publication-panels/figure4_six_metric_pairwise_tests.tsv'),a,c('metric','group1','group2'),
           'Two-sided Wilcoxon rank sum, exact=FALSE, continuity correction','BH: 3 pairs separately per metric',
           'figure4_six_metric_values.tsv; finite observations per metric')

# Figures 1g, 3d, 4g: reproduce seeded specimen-label permutations from frozen timing estimates.
events <- rd(file.path(F1,'panel_f_pcawg_timing_events.tsv')) |>
  filter(source_type=='mono_allelic_gain',category!='all')
stopifnot(nrow(events)==6665)
sums <- events |> group_by(Tumor_Barcode,Subtype_Final) |>
  summarise(weighted_time_sum=sum(time*weight),segment_weight_sum=sum(weight),n_segments=n(),.groups='drop')
a <- weighted_summary(sums,'Subtype_Final') |> rename(Subtype_Final=group)
regenerate(file.path(F1,'panel_g_gain_timing_summary.tsv'),a,'Subtype_Final','Segment-length-weighted mean, specimen and segment counts',
           'None','panel_f_pcawg_timing_events.tsv; exclude duplicated category=all rows')
p <- weighted_permutation(sums,'Subtype_Final',101,plus_one=TRUE); permutations[['Figure1g']] <- p
regenerate(file.path(F1,'panel_g_gain_timing_global_permutation.tsv'),
           tibble(statistic='range_of_subtype_weighted_mean_times',observed=p$observed,pvalue=p$pvalue,n_permutations=p$n_permutations),
           'statistic','5000 specimen-label permutations; seed 101; range of weighted means; (exceedances+1)/5001','None',
           '729 specimen sums from Figure 1 timing events')
sums3 <- sums |> filter(Subtype_Final=='ASTRO') |> inner_join(metrics |> select(Tumor_Barcode,Evo_Group),by='Tumor_Barcode')
stopifnot(nrow(sums3)==241)
p <- weighted_permutation(sums3,'Evo_Group',11,two_group=TRUE); permutations[['Figure3d']] <- p
a <- weighted_summary(sums3,'Evo_Group') |> rename(Evo_Group=group) |>
  mutate(delta_group1_minus_group2=p$observed,permutation_pvalue=p$pvalue,n_permutations=p$n_permutations)
regenerate(file.path(F3,'figure3_gain_timing_stats.tsv'),a,'Evo_Group',
           '5000 specimen subsets of size 161; seed 11; two-sided difference of weighted means; exceedances/5000','None',
           'Figure 1 ASTRO timing events joined to figure3_sample_metrics.tsv')
sums4 <- rd(file.path(F4,'gain-timing/gbm_dn3_panel_d_gain_timing_sample_sums.tsv'))
stopifnot(nrow(sums4)==392)
p <- weighted_permutation(sums4,'Evo_Group',11); permutations[['Figure4g']] <- p
regenerate(file.path(F4,'gain-timing/gbm_dn3_panel_d_gain_timing_global_permutation.tsv'),
           tibble(statistic='range_of_group_weighted_mean_times',observed=p$observed,pvalue=p$pvalue,n_permutations=p$n_permutations),
           'statistic','5000 specimen-label permutations; seed 11; range of weighted means; exceedances/5000','None',
           'gbm_dn3_panel_d_gain_timing_sample_sums.tsv')
a <- weighted_summary(sums4,'Evo_Group') |> rename(Evo_Group=group) |> mutate(global_range_permutation_pvalue=p$pvalue)
regenerate(file.path(F4,'gain-timing/gbm_dn3_panel_d_gain_timing_stats.tsv'),a,'Evo_Group','Segment-length-weighted group summaries',
           'None','gbm_dn3_panel_d_gain_timing_sample_sums.tsv')

# Route the independently fitted multinomial LRT to the exact row read by panel 5i.
path <- file.path(B,'rna/analysis/GBM_state_global.tsv'); tab <- rd(path)
ix <- which(tab$Purity_Source=='BB_Purity' & tab$Model=='joint' & tab$Test=='Multinomial LRT' & tab$Outcome=='Dominant state')
stopifnot(length(ix)==1)
rna <- rd(file.path(out,'Figure5i_state_recomputed.tsv'))
for(v in c('N','Statistic','DF','P')) {
  check(paste0('Figure5i_rendered_',v),rna[[v]],tab[[v]][ix],probability=v=='P'); tab[[v]][ix] <- rna[[v]]
}
write_tsv(tab,path)
coverage[['Figure5i_displayed_LRT']] <- list(rows=1,method='Center/purity-adjusted multinomial likelihood-ratio test',
  correction_family='None',measured_inputs='GBM_state_analysis_metadata.tsv',rendering_input=path,
  status='displayed row recomputed; unused upstream test rows removed from input')
write_json(list(status='PASS',tables=coverage,permutations=permutations),file.path(stage,'qa/extended_statistics.json'),
           pretty=TRUE,auto_unbox=TRUE,digits=16)
