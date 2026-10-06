suppressPackageStartupMessages({library(dplyr);library(readr);library(tidyr);library(purrr)})
root <- normalizePath('.')
out <- file.path(root,'computed_ac2')
groups <- c('DN1','DN2','DN3')
read <- function(p) read_tsv(p,show_col_types=FALSE)
save <- function(x,p) write_tsv(x,p)
fisher <- function(tab) fisher.test(tab,workspace=2e7)$p.value
pairwise <- function(x,y) map_dfr(combn(groups,2,simplify=FALSE),function(g) {
 d<-x[x$DN_Group %in% g,];tibble(group_1=g[1],group_2=g[2],p_value=wilcox.test(d[[y]]~d$DN_Group,exact=FALSE)$p.value)
}) |> mutate(fdr=p.adjust(p_value,'BH'))
for (version in 'current') {
 d<-file.path(out,version,'figure5'); f<-function(n) file.path(d,paste0(n,'.tsv'))
 a<-read(f('panel_a_amplicon_patient_counts')) |> mutate(present=n_amplicons>0)
 ag<-a |> group_by(amplicon_class) |> summarise(p_value=fisher(table(DN_Group,factor(present,levels=c(FALSE,TRUE)))),.groups='drop') |> mutate(fdr=p.adjust(p_value,'BH'))
 save(ag,f('panel_a_amplicon_global_fisher'))
 ap<-map_dfr(unique(a$amplicon_class),function(cl) map_dfr(combn(groups,2,simplify=FALSE),function(g) {
  x<-a[a$amplicon_class==cl&a$DN_Group %in% g,];tibble(amplicon_class=cl,group_1=g[1],group_2=g[2],p_value=fisher(table(x$DN_Group,factor(x$present,levels=c(FALSE,TRUE)))))
 })) |> group_by(amplicon_class) |> mutate(fdr_within_class=p.adjust(p_value,'BH')) |> ungroup() |> mutate(fdr=p.adjust(p_value,'BH'))
 save(ap,f('panel_a_amplicon_pairwise_fisher'))
 # Burden tests include zero counts in all 389 patients; independent classes.
 ab<-a |> group_by(amplicon_class) |> summarise(p_value=if(n_distinct(n_amplicons)>1) kruskal.test(n_amplicons~DN_Group)$p.value else 1,.groups='drop') |> mutate(fdr=p.adjust(p_value,'BH'))
 save(ab,f('panel_a_architecture_burden_global_tests'))
 abp<-map_dfr(unique(a$amplicon_class),function(cl) pairwise(a[a$amplicon_class==cl,],'n_amplicons') |> mutate(amplicon_class=cl)) |> mutate(fdr_within_class=fdr,fdr=p.adjust(p_value,'BH'))
 save(abp,f('panel_a_architecture_burden_pairwise_tests'))
 b<-read(f('panel_b_ecdna_oncogene_burden'))
 save(pairwise(b,'n_distinct_oncogenes') |> mutate(metric='Distinct oncogenes per ecDNA-positive subject',.before=1),f('panel_b_ecdna_oncogene_burden_pairwise_tests'))
 save(tibble(metric='Distinct oncogenes per ecDNA-positive subject',test='Kruskal-Wallis',p_value=kruskal.test(n_distinct_oncogenes~DN_Group,b)$p.value),f('panel_b_ecdna_oncogene_burden_global_test'))
 c<-read(f('panel_c_oncogene_one_vs_rest_fisher')) |> rowwise() |> mutate(p_value=fisher(matrix(c(positive_target,total_target-positive_target,positive_other,total_other-positive_other),nrow=2,byrow=TRUE))) |> ungroup() |> mutate(fdr=p.adjust(p_value,'BH'))
 save(c,f('panel_c_oncogene_one_vs_rest_fisher'))
 cn<-read(f('panel_e_egfr_ecdna_copy_number')); sens<-read(f('panel_e_egfr_copy_number_candidate_sensitivity'))
 expr<-read(f('panel_e_egfr_expression')); fusion<-read(f('panel_e_egfr_fusion_classes'))
 fusion_pairs<-map_dfr(combn(groups,2,simplify=FALSE),function(g) {
  x<-fusion[fusion$DN_Group %in% g,];tibble(group_1=g[1],group_2=g[2],p_value=fisher(xtabs(n_subjects~DN_Group+Architecture,x)))
 }) |> mutate(fdr=p.adjust(p_value,'BH'))
 p<-bind_rows(pairwise(cn,'max_egfr_ecDNA_cn') |> mutate(metric='EGFR ecDNA maximum copy number',.before=1),
              pairwise(expr,'EGFR_expression') |> mutate(metric='EGFR expression',.before=1),
              fusion_pairs |> mutate(metric='EGFR fusion class',.before=1))
 save(p,f('panel_e_pairwise_tests'))
 g<-tibble(metric=c('EGFR ecDNA maximum copy number','EGFR expression','EGFR fusion class'),
           test=c('Kruskal-Wallis','Kruskal-Wallis','Fisher exact'),
           p_value=c(kruskal.test(max_egfr_ecDNA_cn~DN_Group,cn)$p.value,
                     kruskal.test(EGFR_expression~DN_Group,expr)$p.value,
                     fisher(xtabs(n_subjects~DN_Group+Architecture,fusion))))
 save(g,f('panel_e_global_tests'))
 save(pairwise(sens,'max_egfr_ecDNA_cn'),f('panel_e_candidate_sensitivity_pairwise_tests'))
 save(tibble(test='Kruskal-Wallis',p_value=kruskal.test(max_egfr_ecDNA_cn~DN_Group,sens)$p.value),f('panel_e_candidate_sensitivity_global_test'))
 save(b |> group_by(DN_Group) |> summarise(n=n(),median=median(n_distinct_oncogenes),mean=mean(n_distinct_oncogenes),.groups='drop'),f('panel_b_summary'))
 save(cn |> group_by(DN_Group) |> summarise(n=n(),median=median(max_egfr_ecDNA_cn),mean=mean(max_egfr_ecDNA_cn),.groups='drop'),f('panel_e_summary'))
}
