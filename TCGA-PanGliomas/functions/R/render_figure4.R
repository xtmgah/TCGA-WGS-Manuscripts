source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'plot_helpers_34.R'))
suppressPackageStartupMessages(library(survival))
B<-'results/analysis/bcor_reclassification_2026-09-21'
D<-file.path(B,'figures/figure4/derived-data');S<-file.path(B,'survival/analysis')
r<-readRDS(file.path(S,'survival_results.rds'));k<-r$kms$GBM_primary
stopifnot(nrow(k$data)==382,sum(k$data$os_event)==277)
s<-summary(k$fit,censored=TRUE)
surv<-data.frame(time=s$time,surv=s$surv,lower=s$lower,upper=s$upper,nc=s$n.censor,group=sub('grouping=','',as.character(s$strata),fixed=TRUE))
surv<-bind_rows(data.frame(time=0,surv=1,lower=1,upper=1,nc=0,group=c('C19','C19/20','TWR')),surv)
surv$group<-factor(surv$group,levels=c('C19','C19/20','TWR'))
# Preserve the corrected survival panel's existing palette exactly.
kmcols<-c(C19='#F17C7A','C19/20'='#A02723',TWR='#E99CCC')
times<-c(0,50,100,150,200);risk<-summary(k$fit,times=times,extend=TRUE)
riskd<-data.frame(group=sub('grouping=','',as.character(risk$strata),fixed=TRUE),month=risk$time,n=risk$n.risk)
riskd$group<-factor(riskd$group,levels=rev(names(kmcols)))
km<-ggplot(surv,aes(time,surv,color=group,fill=group,group=group))+
 geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.09,color=NA,show.legend=FALSE)+geom_step(linewidth=.55)+
 geom_point(data=filter(surv,nc>0),shape=3,size=.7,stroke=.25,show.legend=FALSE)+
 scale_color_manual(values=kmcols,labels=c('C19 (n = 154)','C19/20 (n = 127)','TWR (n = 101)'),name=NULL)+scale_fill_manual(values=kmcols,guide='none')+
 scale_x_continuous(breaks=times,limits=c(0,max(surv$time)+2),expand=c(0,0))+
 scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=c(0,0))+
 labs(title='Overall survival',subtitle=paste0('382 patients; 277 deaths\nLog-rank P = ',sprintf('%.3f',k$p)),x='Overall survival (months)',y='Survival probability')+
 ref_theme(FALSE)+theme(legend.position='bottom',legend.text=element_text(size=6),legend.key.width=unit(8,'pt'),legend.spacing.x=unit(1,'pt'),plot.subtitle=element_text(size=6,lineheight=.95),axis.text.x=element_text(size=6),axis.title=element_text(size=6.5),plot.margin=margin(1,1,0,1))+
 guides(color=guide_legend(nrow=2,byrow=TRUE))
rt<-ggplot(riskd,aes(month,group,label=n))+geom_text(aes(hjust=ifelse(month==0,0,.5)),size=6/ggplot2::.pt)+
 scale_x_continuous(breaks=times,limits=c(0,max(surv$time)+2),expand=c(0,0))+
 labs(title='At risk',x=NULL,y=NULL)+ref_theme(FALSE)+theme(plot.title=element_text(size=6,hjust=0,margin=margin(b=1)),axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=6),plot.margin=margin(0,1,0,1))
ref_save(km/rt+plot_layout(heights=c(4,1.2)),'M4_d',58,57.5)

cox<-display_read(file.path(S,'Figure4e_display_contrasts.tsv'));cox$label<-factor(cox$label,levels=rev(cox$label));cox$sig<-factor(ifelse(cox$p.value<.05,'P < 0.05','Not significant'),levels=c('P < 0.05','Not significant'))
fp<-ggplot(cox,aes(estimate,label,color=sig))+geom_vline(xintercept=1,linetype='dashed',linewidth=.25,color='#999999')+
 geom_errorbar(aes(xmin=conf.low,xmax=conf.high),orientation='y',width=.14,linewidth=.35)+geom_point(size=1.5)+
 scale_x_log10(breaks=c(.8,1,1.5,2),labels=c('0.8','1.0','1.5','2.0'),limits=c(.78,2.2))+
 scale_color_manual(values=c('P < 0.05'='#008080','Not significant'='grey60'),name=NULL)+
 labs(title='Adjusted Cox model',subtitle='382 patients; 277 deaths',x='Adjusted HR (95% CI; log scale)',y=NULL)+
 ref_theme(FALSE)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6),axis.text.x=element_text(size=6),axis.title.x=element_text(size=6),plot.subtitle=element_text(size=6),legend.position='bottom',legend.text=element_text(size=6),legend.key.width=unit(8,'pt'),legend.spacing.x=unit(1,'pt'))
ref_save(fp,'M4_e',58,46.5)

curve<-display_read(file.path(D,'gain-timing/gbm_dn3_panel_d_gain_timing_curve.tsv'));gs<-display_read(file.path(D,'gain-timing/gbm_dn3_panel_d_gain_timing_stats.tsv'));gt<-display_read(file.path(D,'gain-timing/gbm_dn3_panel_d_gain_timing_global_permutation.tsv'))
curve$group<-factor(unname(dn_names[curve$Evo_Group]),levels=names(dn_cols));gs$group<-factor(unname(dn_names[gs$Evo_Group]),levels=names(dn_cols))
gain<-ggplot(curve,aes(time_display,cum_frac_display,color=group,linetype=group))+geom_hline(yintercept=.5,linetype='dashed',linewidth=.25,color='#999999')+
 geom_vline(data=gs,aes(xintercept=weighted_mean_time,color=group),linetype='dashed',linewidth=.30,show.legend=FALSE)+geom_step(linewidth=.6)+
 scale_color_manual(values=dn_cols,guide='none')+scale_linetype_manual(values=c(C19='solid','C19/20'='longdash',TWR='dotted'),guide='none')+
 scale_x_continuous(breaks=c(0,.5,1),expand=c(0,0))+scale_y_continuous(breaks=c(0,.5,1),labels=label_percent(accuracy=1),expand=c(0,0))+
 coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip='off')+
 annotate('text',x=.03,y=.90,label=paste('Global',ref_p(gt$pvalue[1])),hjust=0,size=6/ggplot2::.pt)+
 annotate('label',x=.98,y=c(.96,.83,.68),label=names(dn_cols),color=dn_cols,fill=alpha('white',.88),linewidth=0,label.padding=unit(.02,'lines'),hjust=1,size=6/ggplot2::.pt)+
 labs(title='Timing of chromosome gains',x='Molecular time',y='Cumulative gain burden')+ref_theme(FALSE)+theme(plot.title=element_text(size=8),axis.text=element_text(size=6),axis.title=element_text(size=6),plot.margin=margin(1,5,1,2))
ref_save(gain,'M4_g',58,28.5)

mv<-display_read(file.path(D,'publication-panels/figure4_six_metric_values.tsv'));mt<-display_read(file.path(D,'publication-panels/figure4_six_metric_pairwise_tests.tsv'))
metricplot<-function(metric,title,ylab) {
 z<-mv|>filter(.data$metric==.env$metric)|>transmute(group=unname(dn_names[DN_Group]),value=value)
 a<-mt|>filter(.data$metric==.env$metric,q_value<.05)|>transmute(g1=unname(dn_names[group1]),g2=unname(dn_names[group2]),label=ref_p(q_value,'q'))
 distribution(z,title,ylab,dn_cols,a,mean_diamond=TRUE,zero=metric%in%c('MRCA age','Age at diagnosis','Latency','MATH','PGA'))
}
hp<-metricplot('MRCA age','MRCA age','Years')|metricplot('Age at diagnosis','Age at diagnosis','Years')|metricplot('Latency','Latency','Years')
ref_save(hp,'M4_h',172,36.5)
ref_save(metricplot('MATH','MATH','Score'),'M4_i',55.333333,33.5)
ref_save(metricplot('PGA','PGA','Genome altered (%)'),'M4_j',55.333333,33.5)
ref_save(metricplot('Ploidy','Ploidy','Copies'),'M4_k',55.333333,33.5)
write_json(scales::gradient_n_pal(c('#117733','#7fbf7b','#f7f7f7','#c9a0cf','#7b2382'),space='Lab')(seq(.005,.995,by=.01)),file.path(REF_STAGE,'qa','M4_time_colors.json'))
ref_export_semantics('M4',list(no_models_or_tests_run=TRUE,KM=list(patients=nrow(k$data),deaths=sum(k$data$os_event),logrank_p=k$p,curve_rows=nrow(surv),risk=riskd),cox=cox,metric_counts=mv|>count(metric,DN_Group),pairwise_tests=mt,gain_curve_rows=nrow(curve),gain_stats=gs,gain_global_test=gt))
