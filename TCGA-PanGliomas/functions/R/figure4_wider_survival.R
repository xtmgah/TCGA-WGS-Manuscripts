source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'plot_helpers_34.R'))
suppressPackageStartupMessages(library(survival))
cairo_pdf(tempfile(fileext='.pdf'),width=7,height=10)
cowplot::set_null_device(function(width,height)grDevices::cairo_pdf(tempfile(fileext='.pdf'),width=width,height=height))
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
times<-seq(0,192,by=24);risk<-summary(k$fit,times=times,extend=TRUE)
riskd<-data.frame(group=sub('grouping=','',as.character(risk$strata),fixed=TRUE),month=risk$time,n=risk$n.risk)
riskd$group<-factor(riskd$group,levels=rev(names(kmcols)))
km<-ggplot(surv,aes(time,surv,color=group,fill=group,group=group))+
 geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.09,color=NA,show.legend=FALSE)+geom_step(linewidth=.55)+
 geom_point(data=filter(surv,nc>0),shape=3,size=.7,stroke=.25,show.legend=FALSE)+
 scale_color_manual(values=kmcols,labels=c('C19 (n = 154)','C19/20 (n = 127)','TWR (n = 101)'),name=NULL)+scale_fill_manual(values=kmcols,guide='none')+
 scale_x_continuous(breaks=times,limits=c(0,max(surv$time)+2),expand=c(0,0))+
 scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=c(0,0))+
 labs(title='Overall survival',subtitle=paste0('382 patients; 277 deaths | Log-rank P = ',sprintf('%.3f',k$p)),x='Overall survival (months)',y='Survival probability')+
 ref_theme(FALSE)+theme(legend.position='inside',legend.position.inside=c(.985,.985),legend.justification.inside=c(1,1),legend.background=element_blank(),legend.text=element_text(size=6),legend.key.width=unit(8,'pt'),legend.spacing.x=unit(1,'pt'),plot.subtitle=element_text(size=6,lineheight=.95),axis.text.x=element_text(size=6),axis.title=element_text(size=6.5),plot.margin=margin(1,1,0,1))+
 guides(color=guide_legend(ncol=1,byrow=TRUE))
rt<-ggplot(riskd,aes(month,group,label=n))+geom_text(aes(hjust=ifelse(month==0,0,.5)),size=6/ggplot2::.pt)+
 scale_x_continuous(breaks=times,limits=c(0,max(surv$time)+2),expand=c(0,0))+
 labs(title='At risk',x=NULL,y=NULL)+ref_theme(FALSE)+theme(plot.title=element_text(size=6,hjust=0,margin=margin(b=1)),axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=6,margin=margin(r=4)),plot.margin=margin(0,1,0,1))
ref_save_original_34(km/rt+plot_layout(heights=c(3,1.3)),'M4_d',105,38.5)


# Saved-fit lookup at the approved display ticks; no model refitting.
write_tsv(surv,file.path(REF_STAGE,'derived','M4_d_survival_values.tsv'))
write_tsv(riskd,file.path(REF_STAGE,'derived','M4_d_24month_risk_counts.tsv'))
# Independently check n.risk directly from the frozen fit's stored event times.
riskcheck<-bind_rows(lapply(names(kmcols),function(g) {
 z<-summary(k$fit,censored=TRUE)
 rows<-sub('grouping=','',as.character(z$strata),fixed=TRUE)==g
 tt<-z$time[rows];nn<-z$n.risk[rows]
 data.frame(group=g,month=times,n=vapply(times,function(t) {
  pos<-which(tt>=t)
  if(length(pos))as.numeric(nn[pos[1]]) else 0
 },numeric(1)))
}))
stopifnot(identical(as.character(riskd$group),riskcheck$group),identical(riskd$month,riskcheck$month),identical(as.numeric(riskd$n),riskcheck$n))
ref_export_semantics('M4_ticks',list(patients=nrow(k$data),deaths=sum(k$data$os_event),logrank_p=k$p,curve_rows=nrow(surv),risk=riskd,ticks=times,x_limits=c(0,max(surv$time)+2),risk_independently_checked=TRUE,no_models_or_tests_run=TRUE))


cox<-display_read(file.path(S,'Figure4e_display_contrasts.tsv'));cox$label<-factor(cox$label,levels=rev(cox$label));cox$sig<-factor(ifelse(cox$p.value<.05,'P < 0.05','Not significant'),levels=c('P < 0.05','Not significant'))
fp<-ggplot(cox,aes(estimate,label,color=sig))+geom_vline(xintercept=1,linetype='dashed',linewidth=.25,color='#999999')+
 geom_errorbar(aes(xmin=conf.low,xmax=conf.high),orientation='y',width=.14,linewidth=.35)+geom_point(size=1.5)+
 scale_x_log10(breaks=c(.8,1,1.5,2),labels=c('0.8','1.0','1.5','2.0'),limits=c(.78,2.2))+
 scale_color_manual(values=c('P < 0.05'='#008080','Not significant'='grey60'),name=NULL)+
 labs(title='Adjusted Cox model',subtitle='382 patients; 277 deaths',x='Adjusted HR (95% CI; log scale)',y=NULL)+
 ref_theme(FALSE)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6),axis.text.x=element_text(size=6),axis.title.x=element_text(size=6),plot.subtitle=element_text(size=6),legend.position='bottom',legend.text=element_text(size=6),legend.key.width=unit(8,'pt'),legend.spacing.x=unit(1,'pt'))
ref_save(fp,'M4_e',65,38.5)
ref_export_semantics('M4_cox',list(cox=cox,no_models_or_tests_run=TRUE))
