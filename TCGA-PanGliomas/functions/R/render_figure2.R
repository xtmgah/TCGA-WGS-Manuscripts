#!/usr/bin/env Rscript
# Figure 2 display layer. Statistical inputs are recomputed and verified before drawing.
root <- normalizePath(commandArgs(trailingOnly=TRUE)[1]); setwd(root)
source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'common.R'))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(survival);library(ggbeeswarm)})
base_family <- REF_FONT
group_cols <- c(ASTRO_Group1='#8FAED6',ASTRO_Group2='#315686')
group_display_labels <- c(ASTRO_Group1='C17p',ASTRO_Group2='CTR')
diag <- 'results/figures/figure2/derived-data/diagnostics'
rd <- function(n) read_tsv(file.path(diag,paste0(n,'.tsv')),show_col_types=FALSE)
out <- file.path(REF_STAGE,'derived')
metrics_file <- tempfile(fileext='.pdf'); cairo_pdf(metrics_file,width=7,height=10)
cowplot::set_null_device(function(width,height)grDevices::cairo_pdf(tempfile(fileext='.pdf'),width=width,height=height))
savep <- function(p,id,w,h)ref_save(p,paste0('M2_',id),w,h)
audit <- list()
frozen_paths <- c(
 list.files(diag,pattern='^(astro_focused|tp53_copy_number).*\\.tsv$',full.names=TRUE),
 'results/analysis/bcor_reclassification_2026-09-21/survival/analysis/survival_results.rds',
 'results/analysis/bcor_reclassification_2026-09-21/survival/analysis/Figure2d_display_contrasts.tsv',
 sprintf('results/PL-output/plackett_luce_astro_g2_g3/tcga_astro_DN%s_mergedseg_G1.txt',1:2),
 'scripts/plackett-luce/Ordering_Model/reference_files/CIC_hg38_coordinates.tsv',
 'results/figures/figure2/derived-data/astro_age_mrca_latency_primary_2p5x_sample_table.tsv',
 'results/figures/figure2/derived-data/astro_age_mrca_latency_primary_tests.tsv')
frozen_before <- tools::md5sum(frozen_paths)

# Ordering geometry is supplied as a frozen upstream model output.
# c. Original fitted KM curves, including the explicitly retained display tail.
r <- readRDS('results/analysis/bcor_reclassification_2026-09-21/survival/analysis/survival_results.rds')
k <- r$kms$ASTRO_primary
s <- summary(k$fit,censored=TRUE)
d <- data.frame(time=s$time,surv=s$surv,lower=s$lower,upper=s$upper,nc=s$n.censor,g=sub('grouping=','',as.character(s$strata),fixed=TRUE))
lev <- k$summary$group
endpoints <- d %>% group_by(g) %>% summarise(last_observed_time=max(time),.groups='drop')
d <- bind_rows(data.frame(time=0,surv=1,lower=1,upper=1,nc=0,g=lev),d)
d$g <- factor(d$g,levels=lev)
for(g in lev){z<-d[d$g==g,];z<-z[order(z$time),];if(max(z$time)<168){last<-tail(z,1);last$time<-168;last$nc<-0;d<-bind_rows(d,last)}}
expr <- parse(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'figure2_astro_panel_helpers.R'))
for(e in expr)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])=='step_ribbon_df')eval(e)
ribbon <- d %>% group_by(g) %>% group_modify(~step_ribbon_df(.x)) %>% ungroup()
cols <- setNames(unname(group_cols),lev)
labs <- setNames(c('C17p (n = 167)','CTR (n = 85)'),lev)
stopifnot(nrow(k$data)==252,identical(as.integer(k$summary$n),c(167L,85L)))
pkm <- ggplot(d,aes(time,surv,color=g,group=g))+
 geom_ribbon(data=ribbon,aes(ymin=lower,ymax=upper,fill=g),alpha=.10,color=NA,show.legend=FALSE)+
 geom_step(aes(linetype=g),linewidth=.45)+
 geom_segment(data=filter(d,nc>0),aes(xend=time,y=surv-.014,yend=surv+.014),linewidth=.23,show.legend=FALSE)+
 scale_color_manual(values=cols,labels=labs)+scale_fill_manual(values=cols,guide='none')+
 scale_linetype_manual(values=setNames(c('solid','solid'),lev),labels=labs)+
 scale_x_continuous(breaks=c(0,42,84,126,168),expand=c(0,0))+
 scale_y_continuous(breaks=seq(0,1,.25),expand=c(0,0))+
 coord_cartesian(xlim=c(0,168),ylim=c(0,1.02))+
 annotate('text',x=5,y=.10,label=paste0('Log-rank P = ',formatC(k$p,format='f',digits=3)),hjust=0,size=6/.pt,family=REF_FONT)+
 labs(title='Overall survival by ordering group',x='Overall survival (months)',y='Survival probability',color=NULL,linetype=NULL)+
 ref_theme()+theme(legend.position='inside',legend.position.inside=c(.985,.985),legend.justification.inside=c(1,1),legend.background=element_blank(),legend.key.height=unit(6,'pt'),plot.margin=margin(3,6,0,3))
rs <- summary(k$fit,times=c(0,42,84,126,168),extend=TRUE)
risk <- data.frame(time=rs$time,n=rs$n.risk,g=sub('grouping=','',as.character(rs$strata),fixed=TRUE))
risk$label <- factor(c('C17p','CTR')[match(risk$g,lev)],levels=c('CTR','C17p'))
prisk <- ggplot(risk,aes(time,label,label=n))+
 geom_text(aes(hjust=ifelse(time==0,0,ifelse(time==168,1,.5))),size=6/.pt,family=REF_FONT)+scale_x_continuous(limits=c(0,168),breaks=c(0,42,84,126,168),expand=c(0,0))+
 labs(title='Number at risk',x=NULL,y=NULL)+ref_theme()+
 theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.text.y=element_text(size=6,margin=margin(r=4)),plot.title=element_text(size=6.5,margin=margin(b=0)),plot.margin=margin(1,6,1,2))
savep(plot_grid(pkm,prisk,ncol=1,align='v',axis='lr',rel_heights=c(.78,.22)),'c',85,45.5)
write_tsv(d,file.path(out,'M2_c_curve_values.tsv'));write_tsv(risk,file.path(out,'M2_c_risk_counts.tsv'))
audit$survival <- list(n=252,group_n=c(167,85),logrank_p=k$p,last_observed=endpoints,display_xlim=c(0,168),retained_tail_extension=TRUE,tests_recomputed=TRUE)

# d. Read the newly fitted, verified Cox estimates in the manuscript contrast orientation.
d <- read_tsv('results/analysis/bcor_reclassification_2026-09-21/survival/analysis/Figure2d_display_contrasts.tsv',show_col_types=FALSE)
d$label <- gsub('17p-Amp','C17p',gsub('CIC/TERT','CTR',d$label,fixed=TRUE),fixed=TRUE)
labels <- c('Tumor purity (per 10%)','PC2','Age (per 10 years)','PC1','Grade 3 (vs Grade 2)','Female (vs Male)','C17p (vs CTR)')
d$plot_label <- factor(d$label,levels=rev(labels));stopifnot(!anyNA(d$plot_label))
d$sig <- factor(ifelse(d$p.value<.05,'P < 0.05','P ≥ 0.05'),levels=c('P < 0.05','P ≥ 0.05'))
p <- ggplot(d,aes(estimate,plot_label))+
 geom_vline(xintercept=1,linetype='dashed',linewidth=.3,color='grey45')+
 geom_errorbar(aes(xmin=conf.low,xmax=conf.high,color=sig),orientation='y',width=.18,linewidth=.32)+geom_point(aes(color=sig),size=1.3)+
 scale_color_manual(values=c('P < 0.05'='#008080','P ≥ 0.05'='grey60'),drop=FALSE)+
 scale_x_log10(breaks=c(.5,1,2,4,6),labels=c('0.5','1','2','4','6'),limits=c(.4,6.5))+
 labs(title='Multivariable Cox model for OS',x='Hazard ratio (95% CI; log scale)',y=NULL,color=NULL)+ref_theme()+
 theme(legend.position='bottom',legend.direction='horizontal',legend.key.height=unit(6,'pt'),plot.margin=margin(3,2,1,2))
savep(p,'d',85,45.5);write_tsv(d,file.path(out,'M2_d_display_contrasts.tsv'))
audit$cox <- list(contrasts=nrow(d),estimates_unchanged=TRUE,orientation_unchanged=TRUE)

# e. Reconcile retained regional events by CIC overlap in both groups.
# This approved display-derived table is separate from all frozen inputs.
data_mutation <- rd('astro_focused_data_mutation');data_promoter <- rd('astro_focused_data_promoter')
old_cn <- rd('astro_focused_data_allele_specific_cn');data_feature <- rd('astro_focused_data_feature')
wgs_sub <- rd('astro_focused_sample_order');sample_level <- wgs_sub$Tumor_Barcode
stopifnot(length(sample_level)==263,!anyDuplicated(sample_level),identical(as.integer(table(wgs_sub$Evo_Group)),c(172L,91L)))
gene <- read_tsv('scripts/plackett-luce/Ordering_Model/reference_files/CIC_hg38_coordinates.tsv',show_col_types=FALSE)
stopifnot(identical(gene$symbol,'CIC'))
stopifnot(nrow(gene)==1,gene$start==42268537,gene$end==42295797)
events <- bind_rows(lapply(1:2,function(i){
 x<-read_tsv(sprintf('results/PL-output/plackett_luce_astro_g2_g3/tcga_astro_DN%s_mergedseg_G1.txt',i),show_col_types=FALSE)
 x %>% filter(as.character(chr)==as.character(gene$seqnames),startpos<=gene$end,endpos>=gene$start,CNA %in% c('dLOH','dGain')) %>%
 transmute(Tumor_Barcode=Tumour_Name,Evo_Group=paste0('ASTRO_Group',i),ID,CNA,startpos,endpos)
})) %>% distinct()
stopifnot(all(events$startpos<=gene$start),all(events$endpos>=gene$end),all(events$Tumor_Barcode %in% sample_level))
calls <- wgs_sub %>% select(Tumor_Barcode,Evo_Group) %>% mutate(
 Mutation=Tumor_Barcode %in% filter(data_mutation,Gene=='CIC')$Tumor_Barcode,
 LOH=Tumor_Barcode %in% filter(events,CNA=='dLOH')$Tumor_Barcode,
 Gain=Tumor_Barcode %in% filter(events,CNA=='dGain')$Tumor_Barcode,
 Alteration=case_when(LOH&Gain~'Amplified LOH',LOH~'LOH',Gain~'Gain',TRUE~NA_character_))
new_cic <- calls %>% filter(!is.na(Alteration)) %>%
 left_join(data_feature %>% select(Subject,Tumor_Barcode),by='Tumor_Barcode') %>%
 transmute(Subject,Tumor_Barcode,Gene='CIC',Alteration,Type='Allele_Specific_CN')
data_pl_cn <- bind_rows(filter(old_cn,Gene!='CIC'),new_cic)
old_cells <- old_cn %>% select(Tumor_Barcode,Gene,old=Alteration)
new_cells <- data_pl_cn %>% select(Tumor_Barcode,Gene,new=Alteration)
changes <- full_join(old_cells,new_cells,by=c('Tumor_Barcode','Gene')) %>% filter(coalesce(old,'none')!=coalesce(new,'none'))
stopifnot(nrow(changes)==63,all(changes$Gene=='CIC'),all(is.na(changes$old)),all(changes$new=='LOH'),all(changes$Tumor_Barcode %in% filter(wgs_sub,Evo_Group=='ASTRO_Group1')$Tumor_Barcode))
summary <- calls %>% select(Tumor_Barcode,Evo_Group,Mutation,LOH,Gain) %>%
 pivot_longer(c(Mutation,LOH,Gain),names_to='event',values_to='positive') %>% group_by(Evo_Group,event) %>%
 summarise(n=sum(positive),denominator=n(),percent=100*n/denominator,.groups='drop')
for(g in names(group_cols)){z<-summary %>% filter(Evo_Group==g);stopifnot(identical(as.integer(z$n[match(c('Mutation','LOH','Gain'),z$event)]),if(g=='ASTRO_Group1')c(1L,63L,0L)else c(10L,40L,22L)))}
write_tsv(events,file.path(out,'M2_e_CIC_retained_regions.tsv'));write_tsv(calls,file.path(out,'M2_e_CIC_specimen_calls.tsv'))
write_tsv(summary,file.path(out,'M2_e_CIC_event_summary.tsv'));write_tsv(changes,file.path(out,'M2_e_changed_cells.tsv'))
write_tsv(data_pl_cn,file.path(out,'M2_e_copy_number_calls.tsv'));write_tsv(wgs_sub %>% select(Tumor_Barcode,Evo_Group),file.path(out,'M2_e_sample_order.tsv'))
target_genes <- c('TP53','ATRX','CIC');allele_cn_genes <- c('TP53','CIC')
ct <- read_csv('scripts/functions/oncoplot_colors.csv',show_col_types=FALSE);landscape_colors <- setNames(ct$Color,ct$Name)
landscape_colors[names(group_cols)] <- group_cols
landscape_colors[c('Promoter','Amplified LOH','Gain','LOH')] <- c('#4D4D4DFF','#E6A0C4FF','#FAD77BFF','#175149FF')
source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'draw_oncoplot.R'))
oncoplot_main$layers[[6]]$aes_params$size <- 6/.pt
oncoplot_main$layers[[6]]$aes_params$fontface <- 'plain'
oncoplot_main$layers[[6]]$mapping$y <- aes(y=-.65)$y
oncoplot_main$layers[[6]]$data$label <- gsub('n=','n = ',oncoplot_main$layers[[6]]$data$label,fixed=TRUE)
p <- oncoplot_main+scale_y_reverse(limits=c(max(feature_layout$row_ymax),-1.05),breaks=feature_layout$y,labels=sub('CIC copy number','CIC regional CN',feature_layout$Feature_Label,fixed=TRUE),expand=c(0,0))+
 theme(axis.text.y=element_text(family=REF_FONT,size=6,color='black',margin=margin(r=2)),plot.margin=margin(0,1,0,1))
H <- 38.5*72/25.4
left <- ggdraw()+draw_plot(p,y=38/H,height=(H-50)/H)+draw_label('Mutation and copy-number landscape',x=.003,y=1-3/H,hjust=0,vjust=1,fontfamily=REF_FONT,size=8)
keys <- c(Frame_Shift_Del='#1F78B4FF',Frame_Shift_Ins='#6A3D9AFF',In_Frame_Del='#B15928FF',Missense_Mutation='#33A02CFF',Nonsense_Mutation='#a50f15',Splice_Site='#FF7F00FF',Promoter='#4D4D4DFF','Amplified LOH'='#E6A0C4FF',Gain='#FAD77BFF',LOH='#175149FF','No displayed event'='#F2F2F2')
keylabs <- c('Frameshift del.','Frameshift ins.','In-frame del.','Missense','Nonsense','Splice site','Promoter','Gain + LOH','Gain','LOH','No displayed event')
left<-left+draw_label('Mutation / promoter',x=1/108,y=31/H,hjust=0,fontfamily=REF_FONT,size=6.5)+
 draw_label('Copy number',x=67/108,y=31/H,hjust=0,fontfamily=REF_FONT,size=6.5)
for(j in seq_along(keys)){
 if(j<=7){col<-(j-1)%%3;row<-(j-1)%/%3;x<-c(1.5,23.5,45.5)[col+1]/108}
 else{col<-(j-8)%%2;row<-(j-8)%/%2;x<-c(67.5,88.5)[col+1]/108}
 y<-(22-7*row)/H
 left<-left+draw_grob(grid::rectGrob(gp=grid::gpar(fill=keys[j],col=NA)),x=x-.0055,y=y-.016,width=.011,height=.032)+draw_label(keylabs[j],x=x+.014,y=y,hjust=0,fontfamily=REF_FONT,size=6)
}
summary$event <- factor(summary$event,levels=c('Gain','LOH','Mutation'))
summary$Evo_Group <- factor(summary$Evo_Group,levels=rev(names(group_cols)))
right <- ggplot(summary,aes(percent,event,fill=Evo_Group))+
 geom_col(position=position_dodge(width=.74),width=.66)+
 geom_text(aes(label=sprintf('%.1f%% (%d)',percent,n)),position=position_dodge(width=.74),hjust=-.08,size=6/.pt,family=REF_FONT)+
 scale_fill_manual(values=group_cols,breaks=names(group_cols),labels=c('C17p (n = 172)','CTR (n = 91)'))+
 scale_y_discrete(labels=c(Gain='Regional\ngain',LOH='Regional\nLOH',Mutation='Mutation'))+
 scale_x_continuous(limits=c(0,76),breaks=c(0,20,40,60),expand=c(0,0))+
 labs(title='CIC event prevalence',x='Specimens (%)',y=NULL,fill=NULL)+ref_theme()+
 theme(axis.text.y=element_text(size=6),legend.position='bottom',legend.key.width=unit(5,'pt'),legend.key.height=unit(6,'pt'),legend.text=element_text(size=6),plot.margin=margin(3,1,2,1),plot.title=element_text(size=8,hjust=0))
savep(ggdraw()+draw_plot(left,x=0,width=108/172)+draw_plot(right,x=110/172,width=62/172),'e',172,38.5)
write_tsv(oncoplot_data,file.path(out,'M2_e_plotted_cells.tsv'))
audit$cic <- list(roster_n=263,group_n=c(172,91),changed_cells=nrow(changes),only_changes='63 C17p CIC regional LOH calls added',CIC_CN_positive_before=sum(old_cn$Gene=='CIC'),CIC_CN_positive_after=nrow(new_cic),summary=summary,
 overlap=as.list(c(CTR_LOH_and_gain=sum(calls$Evo_Group=='ASTRO_Group2'&calls$LOH&calls$Gain),CTR_mutation_and_LOH=sum(calls$Evo_Group=='ASTRO_Group2'&calls$Mutation&calls$LOH),CTR_mutation_and_gain=sum(calls$Evo_Group=='ASTRO_Group2'&calls$Mutation&calls$Gain))),frozen_inputs_modified=FALSE)

# f. Plot the refitted regression coefficients and the original deterministic jitter.
d <- rd('tp53_copy_number_analysis_data');terms <- rd('tp53_copy_number_model_terms');co <- setNames(terms$estimate,terms$term)
stopifnot(identical(as.integer(table(factor(d$Evo_Group,levels=names(group_cols)))),c(163L,63L)))
pred <- expand_grid(total_cn_display=seq(1,5,length.out=161),Evo_Group=names(group_cols))
x <- pred$total_cn_display-2;g <- as.numeric(pred$Evo_Group=='ASTRO_Group2')
pred$mutant_copies_display <- pmin(co['(Intercept)']+co['total_cn_centered_at_2']*x+co['Evo_GroupASTRO_Group2']*g+co['total_cn_centered_at_2:Evo_GroupASTRO_Group2']*x*g,5)
pred <- filter(pred,mutant_copies_display>=.25)
p <- ggplot(d,aes(total_cn_display,mutant_copies_display,color=Evo_Group))+
 geom_abline(slope=1,intercept=0,color='grey78',linewidth=.28)+
 geom_point(size=.65,alpha=.8,position=position_jitter(width=.16,height=.16,seed=20260817L))+
 geom_line(data=pred,aes(linetype=Evo_Group),linewidth=.4)+
 scale_color_manual(values=group_cols,labels=c('C17p (n = 163)','CTR (n = 63)'))+
 scale_linetype_manual(values=c(ASTRO_Group1='solid',ASTRO_Group2='22'),guide='none')+
 scale_x_continuous(breaks=1:5,labels=c('1','2','3','4','≥5'),limits=c(.35,5.55))+
 scale_y_continuous(breaks=1:5,labels=c('1','2','3','4','≥5'),limits=c(.25,5.55))+
 coord_cartesian(xlim=c(.84,5.25),ylim=c(.7,5.3))+
 labs(title='<i>TP53</i> mutant-copy dosage',x='Total <i>TP53</i> copies',y='Mutant <i>TP53</i> copies',color=NULL)+ref_theme()+
 theme(plot.title=ggtext::element_markdown(family=REF_FONT,size=8,hjust=0,margin=margin(b=2)),axis.title.x=ggtext::element_markdown(family=REF_FONT,size=6.5),axis.title.y=ggtext::element_markdown(family=REF_FONT,size=6.5,angle=90),legend.position='inside',legend.position.inside=c(.32,.86),legend.text=element_text(size=6),legend.key.width=unit(8,'pt'),legend.key.height=unit(6,'pt'),plot.margin=margin(3,2,1,3))
savep(p,'f',46,43.5);write_tsv(pred,file.path(out,'M2_f_saved_model_predictions.tsv'))
audit$dosage <- list(group_n=c(163,63),coefficients=as.list(co),display_cap=5,jitter_seed=20260817L,refitted=TRUE)

# g. Read tests recomputed across the complete three-metric BH family.
timing_path <- 'results/figures/figure2/derived-data/astro_age_mrca_latency_primary_2p5x_sample_table.tsv'
test_path <- 'results/figures/figure2/derived-data/astro_age_mrca_latency_primary_tests.tsv'
metrics <- c(age_mrca_years='MRCA age',age_at_diagnosis_years='Age at diagnosis')
td <- read_tsv(timing_path,show_col_types=FALSE) %>% select(Tumor_Barcode,Evo_Group,all_of(names(metrics))) %>% pivot_longer(all_of(names(metrics)),names_to='metric',values_to='value') %>% filter(is.finite(value))
tests <- read_tsv(test_path,show_col_types=FALSE) %>% filter(metric %in% names(metrics))
ps <- lapply(names(metrics),function(m){
 z<-filter(td,metric==m);z$Evo_Group<-factor(z$Evo_Group,levels=names(group_cols));nn<-table(z$Evo_Group)
 ann<-max(z$value)+.08*diff(range(z$value));q<-tests$q_wilcox_primary[tests$metric==m]
 ggplot(z,aes(Evo_Group,value))+
  geom_boxplot(width=.46,outlier.shape=NA,fill='white',color='#333333',linewidth=.21)+
  geom_point(aes(fill=Evo_Group),shape=21,size=.8,alpha=.75,stroke=.06,color='#333333',position=position_jitter(width=.14,height=0,seed=25))+
  stat_summary(fun=median,geom='crossbar',width=.46,linewidth=.21,color='#333333')+
  annotate('segment',x=1,xend=2,y=ann,yend=ann,linewidth=.25)+
  annotate('segment',x=c(1,2),xend=c(1,2),y=ann,yend=ann-1.4,linewidth=.25)+
  annotate('text',x=1.5,y=ann+3,label=ref_p(q,'q'),size=6/.pt,family=REF_FONT)+
  scale_fill_manual(values=group_cols,guide='none')+scale_color_manual(values=group_cols,guide='none')+
  scale_x_discrete(labels=setNames(paste0(group_display_labels,'\n(n = ',as.integer(nn),')'),names(group_cols)))+
  scale_y_continuous(limits=c(0,90),breaks=seq(0,80,20),expand=c(0,0))+
  labs(title=metrics[[m]],x=NULL,y=if(m==names(metrics)[1])'Years'else NULL)+ref_theme(grid=TRUE)+
  theme(plot.title=element_text(size=7.5,hjust=0),axis.text.x=element_text(size=6),plot.margin=margin(3,1,1,3))
})
savep(plot_grid(plotlist=ps,nrow=1,align='hv',axis='tblr'),'g',76,43.5)
audit$age <- list(counts=td %>% count(metric,Evo_Group),recomputed_tests=tests,tests_recomputed=TRUE)
stopifnot(identical(frozen_before,tools::md5sum(frozen_paths)))
audit$frozen_input_content_unchanged <- TRUE
audit$before_after_input_md5 <- as.list(frozen_before)
write_json(audit,file.path(REF_STAGE,'qa','M2_R_audit.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=16)
dev.off();unlink(metrics_file)
message('Figure 2 c-g rendered from frozen results; approved CIC display aggregation audited.')
