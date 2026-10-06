# Recompute tractable figure statistics from frozen measured/processed inputs.
suppressPackageStartupMessages({library(dplyr);library(readr);library(tidyr);library(survival);library(broom);library(jsonlite);library(nnet)})
root <- Sys.getenv('PANGLIOMA_WORKSPACE');setwd(root)
stage <- file.path(root,'stage');out <- file.path(stage,'derived')
rd <- function(p)read_tsv(p,show_col_types=FALSE,progress=FALSE)
source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'statistics_helpers.R'))
checks <- list()
check <- function(name,actual,expected,probability=grepl('p.value$',name)) {
 stopifnot(length(actual)==length(expected),all(is.finite(actual)),all(is.finite(expected)))
 delta<-abs(actual-expected);tol<-if(probability)1e-8*abs(expected)else 1e-10+1e-8*abs(expected)
 checks[[name]] <<- list(pass=all(delta<=tol),maximum_absolute_error=max(delta),comparisons=length(delta))
 if(!all(delta<=tol))stop(paste('Scientific mismatch:',name,'max error',max(delta)))
}
B <- 'results/analysis/bcor_reclassification_2026-09-21'
S <- file.path(B,'survival/analysis');terms<-rd(file.path(S,'all_cox_terms.tsv'));kt<-rd(file.path(S,'all_KM_tests.tsv'))
kms<-list();fits<-list()
for(subtype in c('ASTRO','GBM')) {
 id<-paste0(subtype,'_primary');d<-rd(file.path(S,paste0(subtype,'_km_data.tsv')))
 d$grouping<-factor(d$group)
 f<-survfit(Surv(os_time,os_event)~grouping,data=d);lr<-survdiff(Surv(os_time,os_event)~grouping,data=d)
 pv<-pchisq(lr$chisq,length(lr$n)-1,lower.tail=FALSE)
 expected<-filter(kt,analysis==id);check(paste0(id,'_KM'),c(nrow(d),sum(d$os_event),pv),c(expected$n,expected$events,expected$p_value))
 tab<-as.data.frame(summary(f)$table);tab$group<-sub('grouping=','',rownames(tab),fixed=TRUE)
 tab<-tab |> transmute(group,n=records,events=events,median=median,lower=`0.95LCL`,upper=`0.95UCL`)
 kms[[id]]<-list(fit=f,data=d,summary=tab,p=pv)
 cc<-rd(file.path(S,paste0(subtype,'_cox_data.tsv')))
 cc$group<-factor(cc$group,levels=if(subtype=='ASTRO')c('17p-Amp','CIC/TERT')else c('C19','C19/20','TWR'))
 cc$sex<-factor(cc$sex,levels=c('Female','Male'))
 if(subtype=='ASTRO')cc$grade<-factor(cc$grade,levels=c('G2','G3'))
 rhs<-c('group','purity10','PC1_z','PC2_z','sex','age10',if(subtype=='ASTRO')'grade')
 model<-coxph(as.formula(paste('Surv(os_time,os_event)~',paste(rhs,collapse='+'))),data=cc,ties='efron',x=TRUE)
 actual<-tidy(model,exponentiate=TRUE,conf.int=TRUE);expected<-terms |> filter(analysis==id)
 ix<-match(actual$term,expected$term);stopifnot(!anyNA(ix))
 for(v in c('estimate','std.error','statistic','p.value','conf.low','conf.high'))check(paste(id,v),actual[[v]],expected[[v]][ix])
 write_tsv(actual,file.path(out,paste0(id,'_cox_recomputed.tsv')));fits[[id]]<-model
}
saveRDS(list(kms=kms,fits=fits),file.path(S,'survival_results.rds'))

# Figure 3: recalculate the displayed continuous/binary tests and feature effects.
D<-'results/figures/figure3/derived-data';metrics<-rd(file.path(D,'figure3_sample_metrics.tsv'))
tests<-rd(file.path(D,'figure3_tests.tsv'))
for(v in c('high_quality_sv_count','TL_Ratio','TERT_expr','any_hc_ct','multi_hc_ct')) {
 d<-metrics[is.finite(as.numeric(metrics[[v]])),];x<-d[[v]][d$Evo_Group=='ASTRO_Group1'];y<-d[[v]][d$Evo_Group=='ASTRO_Group2']
 pv<-if(v %in% c('any_hc_ct','multi_hc_ct'))fisher.test(table(d$Evo_Group,d[[v]]))$p.value else wilcox.test(x,y,exact=FALSE)$p.value
 ref<-tests$pvalue[tests$variable==v];check(paste('Figure3',v),pv,ref)
}
sp<-rd(file.path(D,'figure3_specificity_fingerprint.tsv'))
for(i in seq_len(nrow(sp))) {
 d<-metrics[is.finite(metrics[[sp$variable[i]]]),];v<-d[[sp$variable[i]]]
 if(sp$value_transform[i]!='identity')v<-log10(v+1)
 x<-v[d$Evo_Group=='ASTRO_Group1'];y<-v[d$Evo_Group=='ASTRO_Group2']
 effect<-(mean(x)-mean(y))/sqrt(((length(x)-1)*var(x)+(length(y)-1)*var(y))/(length(x)+length(y)-2))
 check(paste('Figure3_feature',sp$variable[i]),c(wilcox.test(x,y,exact=FALSE)$p.value,effect),c(sp$pvalue[i],sp$standardized_mean_difference[i]))
}
sp$adjusted_pvalue<-p.adjust(sp$pvalue,'BH');write_tsv(sp,file.path(out,'figure3_h_BH_adjusted.tsv'))

# Figure 5i: center/purity-adjusted multinomial likelihood-ratio test.
d<-rd(file.path(B,'rna/analysis/GBM_state_analysis_metadata.tsv'))
d$RNA_Center<-factor(d$RNA_Center);d$dominant_state<-factor(d$dominant_state);d$DN_Group<-factor(d$DN_Group,levels=c('DN1','DN2','DN3'))
d$Purity<-as.numeric(scale(d$BB_Purity))
full<-multinom(dominant_state~RNA_Center+Purity+DN_Group,d,trace=FALSE,maxit=1000,Hess=TRUE)
reduced<-multinom(dominant_state~RNA_Center+Purity,d,trace=FALSE,maxit=1000)
stat<-2*(as.numeric(logLik(full))-as.numeric(logLik(reduced)));df<-attr(logLik(full),'df')-attr(logLik(reduced),'df');pv<-pchisq(stat,df,lower.tail=FALSE)
ref<-rd(file.path(B,'rna/analysis/GBM_state_global.tsv')) |> filter(Purity_Source=='BB_Purity',Model=='joint',Test=='Multinomial LRT',Outcome=='Dominant state')
check('Figure5i_adjusted_state',c(nrow(d),stat,df,pv),c(ref$N,ref$Statistic,ref$DF,ref$P))
write_tsv(tibble(N=nrow(d),Statistic=stat,DF=df,P=pv),file.path(out,'Figure5i_state_recomputed.tsv'))
source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'recompute_extended.R'))
write_json(list(status='PASS',checks=checks,upstream_frozen=c('ordering and mixture fits','timing estimates','molecular and RNA preprocessing/classification','DESeq2','GSEA')),
 file.path(stage,'qa/statistical_verification.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(stage,'derived/R_session_info.txt'))
