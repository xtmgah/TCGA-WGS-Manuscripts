# Figure 6: clonal SBS1/SBS5 accumulation with age in TCGA and Chornobyl.
# The two cohort panels are exported together.
run_Fig6 <- function(input_dir, output_dir) {
  ptc_age_packages(c('segmented','cowplot'))
  x <- ptc_age_common(input_dir)
  result <- ptc_age_result('Fig. 6a; Fig. 6b', {
    prob <- ptc_age_load(input_dir,'thyroid_signature_probability_tp_final.RData','thyroid_sbs96_probability')
    dp <- ptc_age_load(input_dir,'DP_info_data.RData','DP_info_data')
    # Sum probability-weighted SBS1+SBS5 assignments for clonal mutations (Clone==Y).
    tdata2 <- prob %>% dplyr::select(Tumor_Barcode,ID,SBS1,SBS5) %>%
      left_join(dp %>% dplyr::select(Tumor_Barcode,ID,Clone),by=c('Tumor_Barcode','ID')) %>%
      filter(Tumor_Barcode %in% x$cohort$Tumor_Barcode,Clone=='Y') %>%
      pivot_longer(c(SBS1,SBS5)) %>% group_by(Tumor_Barcode) %>%
      summarise(value=sum(value,na.rm=TRUE),.groups='drop') %>% left_join(x$covariates,by='Tumor_Barcode')
    data_top0 <- x$alterations
    study_color <- c('TCGA-THCA'='#007BBD','Chornobyl'='#ff7f00')
    ncicolpal <- c('#BB0E3D','#ff7f00')
    select_samples <- data_top0 %>% filter(Gene=='TERT',Type=='Mutation_Driver') %>% pull(Tumor_Barcode)
    select_samples <- data_top0 %>% filter(Gene=='BRAF',Type=='Mutation_Driver',Tumor_Barcode %in% select_samples) %>% pull(Tumor_Barcode)


    tcga <- tdata2 %>% filter(Study=='TCGA-THCA') %>% mutate(Group=if_else(Tumor_Barcode %in% select_samples,'Altered','WT'))

    # ---- 1) TCGA segmented regression ----
    # Adjust the TCGA model for tumor purity and sex.
    set.seed(1) # Stabilize segmented-regression bootstrap restarts.
    fit_tcga <- lm(value ~ age_at_diagnosis + Tumor_Purity + Sex, data = tcga)

    # Initialize the age breakpoint at 50 years.
    seg_tcga <- segmented(fit_tcga, seg.Z = ~ age_at_diagnosis, psi = list(age_at_diagnosis = 50))

    # Breakpoint (estimate and 95% CI)
    bp_est <- as.numeric(summary(seg_tcga)$psi[,"Est."])
    bp_se  <- as.numeric(summary(seg_tcga)$psi[,"St.Err"])
    ci <- confint(seg_tcga)                     # matrix with rownames like "psi1.age_at_diagnosis"
    bp_est <- as.numeric(ci[1, "Est."])
    bp_ci  <- as.numeric(ci[1, c("CI(95%).low","CI(95%).up")])
    bp_ci_low  <- bp_ci[1]
    bp_ci_high <- bp_ci[2]


    # Slopes before/after breakpoint (and their SEs)
    sl <- slope(seg_tcga)$age_at_diagnosis
    slope_pre  <- sl[1,"Est."]
    slope_post <- sl[2,"Est."]
    se_pre     <- sl[1,"St.Err."]
    se_post    <- sl[2,"St.Err."]

    # Create a smooth prediction grid at typical covariate values
    nd_tcga <- data.frame(
      age_at_diagnosis = seq(min(tcga$age_at_diagnosis, na.rm=TRUE),
                             max(tcga$age_at_diagnosis, na.rm=TRUE), length.out = 400),
      Tumor_Purity = mean(tcga$Tumor_Purity, na.rm=TRUE),
      Sex = names(sort(table(tcga$Sex), decreasing = TRUE))[1] # most common sex as reference
    )

    pred_tcga <- cbind(
      nd_tcga,
      fit = as.numeric(predict(seg_tcga, newdata = nd_tcga, se.fit = TRUE)$fit),
      se  = as.numeric(predict(seg_tcga, newdata = nd_tcga, se.fit = TRUE)$se.fit)
    ) %>%
      mutate(lwr = fit - 1.96*se, upr = fit + 1.96*se)

    # ---- 2) TCGA plot with segmented fit ----
    lbl_pre  <- sprintf("slope (pre-bp) = %.1f ± %.1f /yr", slope_pre, se_pre)
    lbl_post <- sprintf("slope (post-bp) = %.1f ± %.1f /yr", slope_post, se_post)
    lbl_bp   <- sprintf("breakpoint = %.1f yrs\n95%% CI: [%.1f, %.1f]", bp_est, bp_ci[1], bp_ci[2])

    p_tcga <- ggplot(tcga, aes(age_at_diagnosis, value)) +
      geom_point(aes(alpha=age_at_diagnosis>52.2),size = 3,pch=21,fill=study_color[1],stroke=0.25) +
      geom_point(data=tcga %>% filter(Group=='Altered'),pch=19,col=ncicolpal[1],size=1)+
      geom_ribbon(data = pred_tcga, aes(x=age_at_diagnosis,ymin = lwr, ymax = upr), inherit.aes = FALSE, alpha = 0.15) +
      geom_line(data = pred_tcga, aes(y = fit), linewidth = 1.2) +
      geom_vline(xintercept = bp_est, linetype = "dashed") +
      scale_alpha_manual(values = c(0.5,1))+
      scale_x_continuous(breaks = pretty_breaks(n=6))+
      scale_y_continuous(breaks = pretty_breaks(n=7))+
      annotate("rect", xmin = bp_ci[1], xmax = bp_ci[2],
               ymin = -Inf, ymax = Inf, alpha = 0.06) +
      annotate("text", x = bp_est, y = max(tcga$value, na.rm=TRUE)*0.95,
               label = lbl_bp, hjust = 0.5, vjust = 1, size = 4.5) +
      annotate("text", x = min(tcga$age_at_diagnosis, na.rm=TRUE)+3,
               y = max(tcga$value, na.rm=TRUE)*0.75, label = lbl_pre, hjust = 0, size = 4.5) +
      annotate("text", x = max(tcga$age_at_diagnosis, na.rm=TRUE)-3,
               y = max(tcga$value, na.rm=TRUE)*1.02, label = lbl_post, hjust = 1, size = 4.5) +
      labs(title = "TCGA-THCA: Segmented regression for SBS1+SBS5 vs. age",
           x = "Age at diagnosis (years)",
           y = "Clonal clock-like mutations (SBS1+SBS5)")+
      guides(alpha='none')+
      theme(axis.text = element_text(size=13),axis.title = element_text(size=14))

    # Chornobyl: adjusted slope label and an unadjusted descriptive line.

    chornobyl <- tdata2 %>% filter(Study=='Chornobyl') %>% mutate(Group=if_else(Tumor_Barcode %in% select_samples,'Altered','WT'))
    fit_ch <- lm(value ~ age_at_diagnosis + Tumor_Purity + Sex, data = chornobyl)
    eq_ch  <- sprintf("slope = %.1f /yr", coef(fit_ch)["age_at_diagnosis"])

    p_ch <- ggplot(chornobyl, aes(age_at_diagnosis, value)) +
      geom_point(alpha = 1, size = 3, pch=21,fill=study_color[2],stroke=0.25) +
      geom_point(data=chornobyl %>% filter(Group=='Altered'),pch=19,col=ncicolpal[1],size=1)+
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, color = "black") +
      annotate("text", x = min(chornobyl$age_at_diagnosis, na.rm=TRUE)+3,
               y = max(chornobyl$value, na.rm=TRUE)*0.92,
               label = eq_ch, hjust = 0, size = 4.5) +
      labs(title = "Chornobyl: Linear regression for SBS1+SBS5 vs. age",
           x = "Age at diagnosis (years)",
           y = "Clonal clock-like mutations (SBS1+SBS5)")+
      scale_x_continuous(breaks = pretty_breaks(n=6))+
      scale_y_continuous(breaks = pretty_breaks(n=7))+
      theme(axis.text = element_text(size=13),axis.title = element_text(size=14))

    # ---- 4) Combine (side-by-side) ----

    combined <- cowplot::plot_grid(p_tcga + ptc_theme(),p_ch + ptc_theme(),align = 'none',nrow = 1)

    statistics <- data.frame(n_TCGA=nrow(tcga),n_Chornobyl=nrow(chornobyl),breakpoint=bp_est,
      breakpoint_low=bp_ci[1],breakpoint_high=bp_ci[2],slope_pre=slope_pre,slope_post=slope_post,
      slope_Chornobyl_adjusted=unname(coef(fit_ch)['age_at_diagnosis']),
      slope_Chornobyl_drawn=unname(coef(lm(value~age_at_diagnosis,data=chornobyl))['age_at_diagnosis']),
      highlight_TCGA=sum(tcga$Group=='Altered'),highlight_Chornobyl=sum(chornobyl$Group=='Altered'))
    list(files=ptc_save(combined,output_dir,'Fig6',14,6),statistics=statistics,
      notes='Highlighted tumors carry both BRAF and TERT mutations. The Chornobyl line is unadjusted; its slope label adjusts for purity and sex. Segmented restart seed is 1.')
  })
  a <- result; a$panel_id <- 'Fig. 6a'
  b <- result; b$panel_id <- 'Fig. 6b'
  list('Fig. 6a'=a,'Fig. 6b'=b)
}
