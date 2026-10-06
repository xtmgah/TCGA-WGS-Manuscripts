# Display-only helpers for Figures 3 and 4. All inferential results are frozen.
source(file.path(Sys.getenv('PANGLIOMA_R_CODE'),'common.R'))
ref_p_scalar <- ref_p
ref_p <- function(x,prefix='P') vapply(x,function(v)ref_p_scalar(v,prefix),character(1))
suppressPackageStartupMessages({library(dplyr);library(readr);library(tidyr);library(patchwork);library(scales)})
# Leave explicit room for Cairo's native glyph bounds at the panel edge.
# This changes the layout margins only; text and marks stay at native size.
ref_save_original_34 <- ref_save
ref_save <- function(p,id,width_mm,height_mm) {
 if(id %in% c(paste0('M3_',letters[1:8]),'M4_e','M4_g','M4_i','M4_j','M4_k'))
  p <- p+theme(plot.margin=margin(3,5,2,3))
 ref_save_original_34(p,id,width_mm,height_mm)
}
astro_cols <- c(C17p='#8FAED6',CTR='#315686')
dn_cols <- c(C19='#E37E78','C19/20'='#8C231D',TWR='#E6A0C4')
dn_names <- c(DN1='C19',DN2='C19/20',DN3='TWR')
astro_names <- c(ASTRO_Group1='C17p',ASTRO_Group2='CTR')
display_read <- function(path) read_tsv(path,show_col_types=FALSE,progress=FALSE)
brackets <- function(p,a) {
 if(!nrow(a)) return(p)
 p+geom_segment(data=a,aes(x=x1,xend=x2,y=y,yend=y),inherit.aes=FALSE,linewidth=.22)+
 geom_segment(data=a,aes(x=x1,xend=x1,y=y,yend=y-tip),inherit.aes=FALSE,linewidth=.22)+
 geom_segment(data=a,aes(x=x2,xend=x2,y=y,yend=y-tip),inherit.aes=FALSE,linewidth=.22)+
 geom_text(data=a,aes(x=(x1+x2)/2,y=y+gap,label=label),inherit.aes=FALSE,size=6/ggplot2::.pt,vjust=0)
}
distribution <- function(dat,title,ylab,cols,annotations=NULL,mean_diamond=FALSE,zero=FALSE,breaks=waiver(),labels=waiver()) {
 dat <- dat[is.finite(dat$value),];lev<-names(cols);dat$group<-factor(dat$group,levels=lev)
 counts<-table(dat$group);xl<-setNames(paste0(lev,'\nn = ',as.integer(counts)),lev)
 ran<-range(dat$value);span<-max(diff(ran),.1)
 p<-ggplot(dat,aes(group,value,fill=group))+
  geom_boxplot(width=.50,outlier.shape=NA,linewidth=.21,color='#333333',fill='white',alpha=1)+
  geom_point(position=position_jitter(width=.13,height=0,seed=13),shape=21,stroke=.06,color='#333333',size=.8,alpha=.75)+
  scale_fill_manual(values=cols,guide='none')+scale_x_discrete(labels=xl,expand=expansion(add=.45))+
  labs(title=title,x=NULL,y=ylab)+ref_theme(TRUE)+theme(axis.text.x=element_text(size=6.5,lineheight=.9),plot.margin=margin(1,2,1,2))
 if(mean_diamond) p<-p+stat_summary(fun=mean,geom='point',shape=23,fill='white',color='#303030',size=1.25,stroke=.25)
 if(!is.null(annotations)&&nrow(annotations)) {
  a<-annotations;a$x1<-match(a$g1,lev);a$x2<-match(a$g2,lev)
  a<-a[order(a$x2-a$x1,a$x1),];a$y<-max(ran)+span*(.15+.22*(seq_len(nrow(a))-1));a$tip<-span*.025;a$gap<-span*.02
  p<-brackets(p,a);upper<-max(a$y)+span*.24
 } else upper<-max(ran)+span*.10
 lower<-if(zero)0 else min(ran)-span*.04
 p+scale_y_continuous(breaks=breaks,labels=labels,expand=expansion(mult=c(0,0)))+coord_cartesian(ylim=c(lower,upper),clip='off')
}
ref_export_semantics <- function(id,obj) write_json(obj,file.path(REF_STAGE,'qa',paste0(id,'_semantics.json')),auto_unbox=TRUE,pretty=TRUE,na='null',digits=16)
