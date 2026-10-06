# Display-only adaptation: retain frozen values and regenerate five wider top panels.
source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'plot_helpers_34.R'))

PLOT_LAYERS <- list(); BOX_COMPARISONS <- list()
ref_save <- function(p,id,width_mm,height_mm) {
 map <- c(M3_b='a',M3_c='b',M3_d='c',M3_e='d',M3_f='e',M3_g='f')
 if(!id %in% names(map))return(invisible(NULL))
 letter<-unname(map[id]);width_mm<-c(a=40,b=31,c=31,d=31,e=31,f=64)[[letter]]
 height_mm<-if(letter=='f')46.5 else 48.5
 p<-p+labs(title=NULL)+theme(plot.margin=margin(if(letter=='f')16 else 21,3,if(letter=='f')2 else 10,2),
   axis.text.x=element_text(size=6,lineheight=.95),axis.text.y=element_text(size=6),
   axis.title=element_text(size=6.5),legend.text=element_text(size=6))
 ref_save_original_34(p,paste0('M3_',letter),width_mm,height_mm)
}
D<-'results/figures/figure3/derived-data'
metrics<-display_read(file.path(D,'figure3_sample_metrics.tsv'))
tests<-display_read(file.path(D,'figure3_tests.tsv'))
metrics$group<-unname(astro_names[metrics$Evo_Group])
pv<-function(v)tests$pvalue[match(v,tests$variable)]
bx<-function(v,title,ylab,w,h,test=v,zero=FALSE,breaks=waiver(),labels=waiver(),transform=identity,id) {
 z<-data.frame(group=metrics$group,value=transform(metrics[[v]]))
 a<-data.frame(g1='C17p',g2='CTR',label=ref_p(pv(test)))
 p<-distribution(z,title,ylab,astro_cols,a,zero=zero,breaks=breaks,labels=labels)
 original <- ggplot_build(p)
 p$layers[[1]]$stat_params$width <- .42
 p$layers[[2]]$position <- position_jitter(width=.18,height=0,seed=13)
 updated <- ggplot_build(p)
 fields <- c('ymin','lower','middle','upper','ymax')
 stopifnot(identical(original$data[[1]][,fields],updated$data[[1]][,fields]),identical(original$data[[2]]$y,updated$data[[2]]$y))
 BOX_COMPARISONS[[id]] <<- list(summary_values_identical=TRUE,point_y_values_identical=TRUE,points=nrow(updated$data[[2]]))
 ref_save(p,id,w,h)
 invisible(z)
}
za<-bx('TP53_PL_max_mutant_copies_capped5',expression(italic(TP53)*' amplified LOH'),'Mutant allele copies',26,48.5,test='TP53_PL_max_mutant_copies',zero=TRUE,breaks=0:5,labels=c('0','1','2','3','4','≥5'),id='M3_a')
zc<-bx('high_quality_sv_count','SV burden','log10(SVs + 1)',24,48.5,transform=function(x)log10(x+1),zero=TRUE,breaks=0:3,id='M3_c')
ze<-bx('TL_Ratio','Telomere ratio','log2(T/N)',24,48.5,breaks=c(-1,0,1,2,3),id='M3_d')
zf<-bx('TERT_expr',expression(italic(TERT)*' expression'),'VST expression',24,48.5,breaks=c(2,4,6,8),id='M3_e')

# Prevalence bars retain the original denominators and nominal Fisher tests.
b<-tests|>filter(variable%in%c('any_hc_ct','multi_hc_ct'))
bd<-bind_rows(lapply(seq_len(nrow(b)),function(i)data.frame(metric=c('Any HC CT','Multi-chr CT')[i],group=names(astro_cols),count=c(b$value_group1[i],b$value_group2[i]),n=c(b$n_group1[i],b$n_group2[i]),pct=c(b$pct_group1[i],b$pct_group2[i]))))
bd$metric<-factor(bd$metric,levels=c('Any HC CT','Multi-chr CT'));bd$group<-factor(bd$group,levels=names(astro_cols));bd$x<-as.numeric(bd$metric)+ifelse(bd$group=='C17p',-.20,.20)
ba<-data.frame(x1=c(.8,1.8),x2=c(1.2,2.2),y=c(50,26),tip=1,gap=1,label=ref_p(b$pvalue))
p<-ggplot(bd,aes(x,pct,fill=group))+geom_col(width=.35)+geom_text(aes(label=paste0(formatC(pct,digits=1,format='f'),'%')),vjust=-.35,size=6/ggplot2::.pt)+
 scale_fill_manual(values=astro_cols,labels=c('C17p (n = 172)','CTR (n = 91)'),name=NULL)+
 scale_x_continuous(breaks=1:2,labels=c('Any HC\nCT','Multi-chr\nCT'))+scale_y_continuous(limits=c(0,58),breaks=c(0,20,40),labels=function(x)paste0(x,'%'),expand=c(0,0))+
 labs(title='CT prevalence',x=NULL,y='Specimens (%)')+ref_theme(TRUE)+theme(legend.position='none',legend.direction='vertical',legend.key.height=unit(6,'pt'),legend.text=element_text(size=6),axis.text.x=element_text(size=6.5,lineheight=.9))
ref_save(brackets(p,ba),'M3_b',38,48.5)

# The approved 75%/68% display caps are intentionally retained.
curve<-display_read(file.path(D,'figure3_gain_timing_curve.tsv'));gs<-display_read(file.path(D,'figure3_gain_timing_stats.tsv'))
curve$group<-factor(unname(astro_names[curve$Evo_Group]),levels=names(astro_cols));curve$display<-pmin(curve$cum_frac,ifelse(curve$group=='C17p',.75,.68))
gs$group<-factor(unname(astro_names[gs$Evo_Group]),levels=names(astro_cols))
p<-ggplot(curve,aes(time,display,color=group))+geom_hline(yintercept=.5,linetype='dashed',linewidth=.25,color='#999999')+
 geom_step(linewidth=.6)+geom_vline(data=gs,aes(xintercept=weighted_mean_time,color=group),linetype='dashed',linewidth=.32,show.legend=FALSE)+
 scale_color_manual(values=astro_cols,labels=c('C17p (n = 161)','CTR (n = 80)'),name=NULL)+
 scale_x_continuous(breaks=c(0,.5,1),expand=c(0,0))+scale_y_continuous(breaks=c(0,.25,.5,.75),labels=label_percent(accuracy=1),expand=c(0,0))+
 coord_cartesian(xlim=c(0,1.01),ylim=c(0,.80),clip='off')+
 annotate('text',x=.03,y=.76,label=sprintf('P = %.3f',gs$permutation_pvalue[1]),hjust=0,size=6/ggplot2::.pt)+
 annotate('text',x=.03,y=.035,label='Earlier',hjust=0,size=6/ggplot2::.pt)+annotate('text',x=.98,y=.035,label='Later',hjust=1,size=6/ggplot2::.pt)+
 labs(title='Timing of chromosome gains',x='Molecular time',y='Total gain burden')+ref_theme(FALSE)+theme(legend.position='bottom',legend.key.height=unit(6,'pt'),legend.text=element_text(size=6))
ref_save(p,'M3_g',64,46.5)

tg<-display_read(file.path(D,'figure3_tert_promoter_test.tsv'))
gd<-data.frame(group=factor(names(astro_cols),levels=names(astro_cols)),count=c(tg$value_group1,tg$value_group2),n=c(tg$n_group1,tg$n_group2),pct=c(tg$pct_group1,tg$pct_group2))
p<-ggplot(gd,aes(group,pct,fill=group))+geom_col(width=.42)+geom_text(aes(label=sprintf('%.1f%%\n(%s/%s)',pct,count,n)),vjust=-.4,size=6/ggplot2::.pt)+
 scale_fill_manual(values=astro_cols,guide='none')+scale_x_discrete(labels=c('C17p'='C17p\nn = 172','CTR'='CTR\nn = 91'),expand=expansion(add=.65))+scale_y_continuous(breaks=c(0,10,20),labels=function(x)paste0(x,'%'),limits=c(0,31),expand=c(0,0))+
 labs(title=expression(italic(TERT)*' promoter hotspots'),x=NULL,y='Specimens (%)')+ref_theme(TRUE)
ref_save(brackets(p,data.frame(x1=1,x2=2,y=27,tip=.7,gap=.4,label=ref_p(tg$pvalue))),'M3_f',26,48.5)

sp<-display_read(file.path(D,'figure3_specificity_fingerprint.tsv'))
lev<-c('CIN25','CIN70','FGA','Arm CNA burden','Aneuploidy','Ploidy','Inversion burden','Oscillating CN segments','Maximum clustered SVs','Break load','Translocation burden','HRD-LOH')
sp$metric<-factor(sp$metric,levels=rev(lev));sp$call<-factor(ifelse(sp$significant,'C17p higher (P < 0.05)','Not significant'),levels=c('C17p higher (P < 0.05)','Not significant'))
p<-ggplot(sp,aes(standardized_mean_difference,metric,fill=call))+geom_vline(xintercept=0,linewidth=.25,color='#999999')+geom_point(shape=21,size=1.75,stroke=.3,color='#444444')+
 scale_fill_manual(values=c('C17p higher (P < 0.05)'=unname(astro_cols[1]),'Not significant'='white'),name=NULL)+scale_x_continuous(limits=c(0,.47),breaks=c(0,.2,.4),expand=c(0,0))+
 labs(title='Genomic feature differences',x='Standardized C17p − CTR effect',y=NULL)+ref_theme(FALSE)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6.5),legend.position='bottom',legend.direction='horizontal',legend.text=element_text(size=6),legend.key.height=unit(7,'pt'),plot.margin=margin(1,2,1,2))
ref_save(p,'M3_h',106,46.5)
ref_export_semantics('M3',list(no_models_or_tests_run=TRUE,display_caps=c(C17p=.75,CTR=.68),gain_curve_rows=nrow(curve),gain_statistics=gs,prevalence=bd,TERT_promoter=gd,box_counts=list(a=as.list(table(za$group)),c=as.list(table(zc$group)),e=as.list(table(ze$group)),f=as.list(table(zf$group[is.finite(zf$value)]))),feature_rows=as.character(sp$metric),feature_effects=sp$standardized_mean_difference))

write_json(PLOT_LAYERS,file.path(REF_STAGE,'qa/plot_layers.json'),auto_unbox=TRUE,pretty=FALSE,na='null',digits=16)
write_json(BOX_COMPARISONS,file.path(REF_STAGE,'qa/box_preservation.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
