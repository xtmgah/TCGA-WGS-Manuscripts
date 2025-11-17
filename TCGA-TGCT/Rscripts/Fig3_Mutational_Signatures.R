# ------------------------------------------------------------------------------
# Script: Figure 3 - TCGA-TGCT WGS analysis
# Description: This script reproduces all main analyses and figures for the follwoing TGCT manuscript: 
# A Thiopurine-like Mutagenic Process Defines TGCT Subtypes
# ------------------------------------------------------------------------------

# --- Load Required Libraries and Set Plotting Styles ---------------------------
# (You may need to install some packages if missing)

source('../functions/Sherlock_functions.R')

set_wd()
libztw()
pdfhr()
pdfhr2()



# Fig. 2a -----------------------------------------------------------------
load('../data/tgct_mutational_signatures.RData',verbose = T)
# TreeAndLeaf plot --------------------------------------------------------
if(profile_type == 'SBS') { tmp <- tgct_sbs96_profile %>% mutate(Type=str_remove(str_remove(MutationType,'.*\\['),'\\].*')) }
#if(profile_type == 'ID') { tmp <- seqmatrix_refdata %>% mutate(Type=str_remove(MutationType,':[^:]*$')) }
#if(profile_type == 'DBS') { tmp <- seqmatrix_refdata %>% mutate(Type=str_remove(MutationType,'>.*'),Type= paste0(Type,'>NN'))}

dmdata <- tmp %>%
  group_by(Tumor_Barcode,Type) %>% 
  summarise(value=sum(Contribution)) %>%
  mutate(value=value/(sum(value))) %>%
  arrange(desc(value)) %>% 
  slice(1) %>%
  ungroup() %>%
  rename(Dmut=Type,Dmvalue=value)

dsigdata <- tgct_activity_all_ratio %>% pivot_longer(cols = -Tumor_Barcode) %>% group_by(Tumor_Barcode)%>%  arrange(desc(value)) %>% slice(1) %>% ungroup() %>% rename(Dsig=name,Dvalue=value)


# Combined all the covariables as mdata
seqmatrix_refdata_ratio <- 
  bind_rows(
    tgct_sbs96_profile %>% group_by(Tumor_Barcode) %>% mutate(Mutations = Contribution / sum(Contribution)) %>% ungroup(),
    tgct_id83_profile %>% group_by(Tumor_Barcode) %>% mutate(Mutations = Contribution / sum(Contribution)) %>% ungroup()
  )

seqmatrix_refdata_ratio <- tgct_sbs96_profile %>% group_by(Tumor_Barcode) %>% mutate(Mutations = Contribution / sum(Contribution)) %>% ungroup()

mdata <- seqmatrix_refdata_ratio %>% 
  select(MutationType,Mutations,Tumor_Barcode) %>% 
  pivot_wider(names_from = MutationType,values_from = Mutations)

mdatax <- mdata %>% 
  select(Tumor_Barcode) %>% 
  left_join(dsigdata) %>% 
  left_join(dmdata) %>% 
  left_join(
    tgct_activity_all %>% pivot_longer(-Tumor_Barcode) %>% group_by(Tumor_Barcode) %>% summarise(Mutations=sum(value)) %>% ungroup()  ## add mutation as cell size
  ) %>% 
  left_join(
    wgs_data %>% select(Tumor_Barcode,Subtype)
  ) %>% 
  left_join(tgct_activity_all_ratio)


### Perform treeAndleaf visualization
library("TreeAndLeaf")
library("RedeR")
library("igraph")
library("RColorBrewer")

# convert to matrix for clustering
mdata0 <- as.matrix(mdata[,-1])
rownames(mdata0) <- mdata$Tumor_Barcode

hc <- hclust(dist(mdata0), "ward.D") # hc will be used for treeandlef


# TreeAndLeaf -------------------------------------------------------------
tal <- treeAndLeaf(hc)
#--- Map attributes to the tree-and-leaf
#Note: 'refcol = 0' indicates that 'dat' rownames will be used as mapping IDs
tal <- att.mapv(tal, mdatax, refcol = 1)
tal <- att.setv(g = tal, from = "Mutations", to = "nodeSize",nquant = 8)
#--- Set graph attributes using 'att.addv' and 'att.adde' functions
tal <- att.addv(tal, "nodeFontSize", value = 1)
tal <- att.adde(tal, "edgeWidth", value = 10)
# 
# # Update filling Colors, select one of them -----------------------------------------------------------
# # For dominant signature
if(FALSE){
  dsigname <- unique(mdatax$Dsig)
  palcol <- sigcol[dsigname]
  #match color
  palcol <- as.character(palcol[levels(as.factor(mdatax$Dsig))])
  ## if signature name did not included in SBScolor, please random pick one color from pal_d3(palette = 'category20')(20)
  tal <- att.setv(g = tal, from = "Dsig", to = "nodeColor",cols = palcol)
  update_color <- function(){addLegend.color(obj = rdp, tal, title = "Dominant Mutational Signature", position = "topright",vertical=T,dxtitle=50)}
}
# 
# 
# # For dominant mutation types
# if(FALSE){
#   tmp <- sigprofiles %>% select(Type,Color) %>% unique()
#   mutcols <- tmp$Color
#   names(mutcols) <- tmp$Type
#   palcol <- as.character(mutcols[levels(as.factor(mdatax$Dmut))])
#   tal <- att.setv(g = tal, from = "Dmut", to = "nodeColor",cols = palcol)
#   update_color <- function(){addLegend.color(obj = rdp, tal, title = "Dominant MutationType", position = "topright",vertical=T,dxtitle=50)}
#   
# }
# 
# # Cosine Similarity
# if(FALSE){
#   # summary(mdatax$Cosine_similarity)
#   breakstmp <- seq(min(mdatax$Cosine_similarity),1,0.01)
#   palcol <- viridis::viridis(n = length(breakstmp))
#   tal <- att.setv(g = tal, from = "Cosine_similarity", to = "nodeColor",cols = palcol,breaks = breakstmp)
#   update_color <- function(){addLegend.color(obj = rdp, tal, title = "Cosine Similarity", position = "topright",vertical=T,dxtitle=50)}
# }


if(TRUE){
  tal <- att.setv(g = tal, from = "Subtype", to = "nodeColor",cols = rev(as.character(subtypecol)))
  update_color <- function(){addLegend.color(obj = rdp, tal, title = "Subtype", position = "bottomright",vertical=T,dxtitle=50)}
}

if(FALSE){
  # by signature
  breakstmp <- seq(0,max(mdatax$SBS87),0.025)
  pal <- brewer.pal(length(breakstmp)-1, "Greens")
  pal <- c('white',pal)
  
  tal <- att.setv(g = tal, from = "SBS87", to = "nodeColor",cols = pal,breaks = breakstmp)
  update_color <- function(){addLegend.color(obj = rdp, tal, title = "Subtype", position = "bottomright",vertical=T,dxtitle=50)}
}


#--- Call RedeR application
rdp <- RedPort()
calld(rdp)
resetd(rdp)
addGraph(obj = rdp, g = tal, gzoom=8)
#--- Call 'relax' to fine-tune the leaf nodes
relax(rdp, p1=25, p2=200, p3=50, p4=100, p5=10,p6=40,p7=30)
addLegend.size(obj = rdp, tal, title = "Number of Mutations",position = "topleft",vertical=T)
update_color()




# Fig. 3b -----------------------------------------------------------------

# Mutational signature landscape SBS ------------------------------------------
tdata <- bind_rows(
  tgct_activity_all_ratio %>% 
    select(Tumor_Barcode,starts_with('SBS')) %>% 
    pivot_longer(-Tumor_Barcode) %>% 
    mutate(group='Ratio'),
  
  tgct_activity_all %>% 
    select(Tumor_Barcode,starts_with('SBS')) %>% 
    pivot_longer(-Tumor_Barcode) %>% 
    group_by(Tumor_Barcode) %>% 
    summarise(value=sum(value)) %>% 
    mutate(group='Count')
) %>% 
  left_join(
    wgs_data %>% select(Tumor_Barcode,Subtype)
  ) %>% 
  mutate(name = factor(name,levels = sort(names(tgct_activity_all)[-1])))

## ordered by SBS3
samlevs <- tgct_activity_all_ratio %>% left_join(wgs_data) %>% arrange(Subtype,desc(SBS3+SBS18+SBS87)) %>% group_by(Subtype) %>%  mutate(seq=seq_along(Tumor_Barcode)) %>% ungroup() %>% select(Tumor_Barcode,seq)

# ## ordered by clustering
# 
#   ordered_barcodes1 <- get_ordered_barcodes(msa_sbs96_activity %>% filter(Tumor_Barcode %in% (wgs_data %>% filter(Subtype=='Seminoma') %>% pull(Tumor_Barcode))), "Tumor_Barcode")
# ordered_barcodes2 <- get_ordered_barcodes(msa_sbs96_activity %>% filter(Tumor_Barcode %in% (wgs_data %>% filter(Subtype!='Seminoma') %>% pull(Tumor_Barcode))), "Tumor_Barcode")
# 
# samlevs <- bind_rows(
#   tibble(Tumor_Barcode=ordered_barcodes1) %>% mutate(seq=seq_along(Tumor_Barcode)),
#   tibble(Tumor_Barcode=ordered_barcodes2) %>% mutate(seq=seq_along(Tumor_Barcode))
# )


pdata <- tdata %>% left_join(samlevs)
max(pdata$value)
scale_y <- 5000
pdata_ratio <- pdata %>% filter(group == 'Ratio') %>% mutate(value=value * scale_y)
pdata_count <- pdata %>% filter(group == 'Count')

# cosine similarity
pdata_cos <- tgct_sbs96_msa_decompsite %>% 
  select(Tumor_Barcode=Sample_Names,Cosine_similarity) %>% 
  left_join(samlevs) %>% 
  left_join(
    wgs_data %>% select(Tumor_Barcode,Subtype)
  )

# Create the combined plot
p1 <- ggplot() +
  geom_bar(data = pdata_ratio, 
           aes(x = seq, y = value, fill = name), 
           stat = "identity",width = 1) +
  scale_fill_manual(values = sigcol) +
  geom_point(data = pdata_count, 
             aes(x = seq, y = value, color = name), 
             inherit.aes = FALSE,col='#6baed6',size=2) +
  geom_point(data = pdata_count, 
             aes(x = seq, y = value, color = name), 
             inherit.aes = FALSE,col='black',size=0.1) +
  geom_line(data = pdata_count, 
            aes(x = seq, y = value, color = name), 
            inherit.aes = FALSE,col='white',linewidth=0.15) +
  facet_wrap(~Subtype,scales = 'free_x')+
  scale_y_continuous(
    name = "Mutational signature proportion",
    expand = c(0,0),
    breaks = seq(0,5000,length=6),labels = seq(0,1,length=6), #percent_format()(seq(0,1,length=6))
    sec.axis = sec_axis(~ .,breaks = seq(0,5000,length=6),name = "Number of mutations") # Secondary y-axis
  ) +
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0))+
  theme_ipsum(base_size = 12,grid = F,ticks = T,axis_title_just = 'm',axis_title_size = 14)+
  theme(
    axis.title.y.right = element_text(color = "#08519c"), # Customize right y-axis label
    axis.text.y.right = element_text(color = "#08519c"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(hjust = 0.5,face = 2),
    panel.spacing = unit(0.1,"cm"),
    legend.position = 'bottom')+
  labs(x = "Tumor sample", fill = "Signature")


p2 <- ggplot() +
  geom_tile(data=pdata_cos,
            aes(x=seq,y=0.5,fill=Cosine_similarity))+
  facet_wrap(~Subtype,scales = 'free_x')+
  scale_y_continuous(
    name = "Mutational signature proportion",
    expand = c(0,0),
    breaks = seq(0,5000,length=6),labels = seq(0,1,length=6), #percent_format()(seq(0,1,length=6))
    sec.axis = sec_axis(~ .,breaks = seq(0,5000,length=6),name = "Number of mutations") # Secondary y-axis
  ) +
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0))+
  scale_fill_viridis_c(limit=c(0,1))+
  theme_ipsum(base_size = 12,grid = F,ticks = T,axis_title_just = 'm',axis_title_size = 14)+
  theme(
    axis.title.y.right = element_text(color = "#08519c"), # Customize right y-axis label
    axis.text.y.right = element_text(color = "#08519c"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(hjust = 0.5,face = 2),
    panel.spacing = unit(0.1,"cm"),
    legend.position = 'bottom',
    legend.key.width = unit(1.5,'cm'),
    legend.key.height = unit(0.4,'cm')
  )+
  labs(x = "Tumor sample", fill = "Cosine similarity\n")


plot_grid(p1,p2,align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2.24,1))



# Fig. 3c -----------------------------------------------------------------

# Mutational signature landscape ID ------------------------------------------
tdata <- bind_rows(
  tgct_activity_all_ratio %>% 
    select(Tumor_Barcode,starts_with('ID')) %>% 
    pivot_longer(-Tumor_Barcode) %>% 
    mutate(group='Ratio'),
  
  tgct_activity_all %>% 
    select(Tumor_Barcode,starts_with('ID')) %>% 
    pivot_longer(-Tumor_Barcode) %>% 
    group_by(Tumor_Barcode) %>% 
    summarise(value=sum(value)) %>% 
    mutate(group='Count')
) %>% 
  left_join(
    wgs_data %>% select(Tumor_Barcode,Subtype)
  ) %>% 
  mutate(name = factor(name,levels = sort(names(tgct_activity_all)[-1])))

## ordered by clustering

ordered_barcodes1 <- get_ordered_barcodes(msa_id83_activity %>% filter(Tumor_Barcode %in% (wgs_data %>% filter(Subtype=='Seminoma') %>% pull(Tumor_Barcode))), "Tumor_Barcode")
ordered_barcodes2 <- get_ordered_barcodes(msa_id83_activity %>% filter(Tumor_Barcode %in% (wgs_data %>% filter(Subtype!='Seminoma') %>% pull(Tumor_Barcode))), "Tumor_Barcode")

samlevs <- bind_rows(
  tibble(Tumor_Barcode=ordered_barcodes1) %>% mutate(seq=seq_along(Tumor_Barcode)),
  tibble(Tumor_Barcode=ordered_barcodes2) %>% mutate(seq=seq_along(Tumor_Barcode))
)


pdata <- tdata %>% left_join(samlevs)
max(pdata$value)
scale_y <- 250
pdata_ratio <- pdata %>% filter(group == 'Ratio') %>% mutate(value=value * scale_y)
pdata_count <- pdata %>% filter(group == 'Count')

# cosine similarity
pdata_cos <- tgct_id83_msa_decompsite %>% 
  select(Tumor_Barcode=Sample_Names,Cosine_similarity) %>% 
  left_join(samlevs) %>% 
  left_join(
    wgs_data %>% select(Tumor_Barcode,Subtype)
  )

# Create the combined plot
p1 <- ggplot() +
  geom_bar(data = pdata_ratio, 
           aes(x = seq, y = value, fill = name), 
           stat = "identity",width = 1) +
  scale_fill_manual(values = sigcol) +
  geom_point(data = pdata_count, 
             aes(x = seq, y = value, color = name), 
             inherit.aes = FALSE,col='#6baed6',size=2) +
  geom_point(data = pdata_count, 
             aes(x = seq, y = value, color = name), 
             inherit.aes = FALSE,col='black',size=0.1) +
  geom_line(data = pdata_count, 
            aes(x = seq, y = value, color = name), 
            inherit.aes = FALSE,col='white',linewidth=0.15) +
  facet_wrap(~Subtype,scales = 'free_x')+
  scale_y_continuous(
    name = "Mutational signature proportion",
    expand = c(0,0),
    breaks = seq(0,250,length=6),labels = seq(0,1,length=6), #percent_format()(seq(0,1,length=6))
    sec.axis = sec_axis(~ .,breaks = seq(0,250,length=6),name = "Number of mutations") # Secondary y-axis
  ) +
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0))+
  theme_ipsum(base_size = 12,grid = F,ticks = T,axis_title_just = 'm',axis_title_size = 14)+
  theme(
    axis.title.y.right = element_text(color = "#08519c"), # Customize right y-axis label
    axis.text.y.right = element_text(color = "#08519c"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(hjust = 0.5,face = 2),
    panel.spacing = unit(0.1,"cm"),
    legend.position = 'bottom')+
  labs(x = "Tumor sample", fill = "Signature")


p2 <- ggplot() +
  geom_tile(data=pdata_cos,
            aes(x=seq,y=0.5,fill=Cosine_similarity))+
  facet_wrap(~Subtype,scales = 'free_x')+
  scale_y_continuous(
    name = "Mutational signature proportion",
    expand = c(0,0),
    breaks = seq(0,5000,length=6),labels = seq(0,1,length=6), #percent_format()(seq(0,1,length=6))
    sec.axis = sec_axis(~ .,breaks = seq(0,5000,length=6),name = "Number of mutations") # Secondary y-axis
  ) +
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0,0))+
  scale_fill_viridis_c(limit=c(0,1))+
  theme_ipsum(base_size = 12,grid = F,ticks = T,axis_title_just = 'm',axis_title_size = 14)+
  theme(
    axis.title.y.right = element_text(color = "#08519c"), # Customize right y-axis label
    axis.text.y.right = element_text(color = "#08519c"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(hjust = 0.5,face = 2),
    panel.spacing = unit(0.1,"cm"),
    legend.position = 'bottom',
    legend.key.width = unit(1.5,'cm'),
    legend.key.height = unit(0.4,'cm')
  )+
  labs(x = "Tumor sample", fill = "Cosine similarity\n")


plot_grid(p1,p2,align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2.6,1))



# Fig. 3d -----------------------------------------------------------------
# Signature and replication stress --------------------------------------------
load('tcga_repstress_zscore.RData',verbose = T)

pdata <- tgct_activity_all %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  filter(name %in% c('ID1','ID2','ID9')) %>% 
  mutate(TCGA_Barcode=str_sub(Tumor_Barcode,1,15)) %>% 
  left_join(tcga_repstress_score_cbioportal) %>% 
  filter(!is.na(RepStress_ZScore)) %>% 
  left_join(wgs_data %>% select(Tumor_Barcode,Subtype))

pdata %>% group_by(name) %>% do(tidy(lm(RepStress_ZScore~value+Subtype,data=.))) %>% filter(term=='value') %>% arrange(p.value)

pdata2 <- pdata %>% group_by(Subtype,name) %>% do(tidy(cor.test(.$RepStress_ZScore, log2(.$value+1),data=.))) %>% arrange((p.value)) %>% group_by(Subtype) %>% mutate(FDR=p.adjust(p.value)) %>% ungroup() %>% arrange(FDR) %>% 
  mutate(p=if_else(p.value<0.05,p.value,NA_real_))

pdata2 %>% 
  ggplot(aes(name,Subtype,fill=-log10(p),size=estimate)) +
  geom_point(pch=21)


pdata %>% 
  ggplot(aes(log2(value+1),RepStress_ZScore))+
  geom_point(aes(fill=Subtype),pch=21,size=2.5)+
  stat_smooth(method = "lm") +
  stat_cor(size=5,col=ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5))+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  facet_wrap(~Subtype+name,scale='free')+
  #facet_grid(Subtype~name,scales = 'free',space = 'free')+
  labs(x='log2(Mutations+1)',y='Repstress Z-score')+
  theme_ipsum_rc(axis_text_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'XY',ticks = T)+
  panel_border(color = 'black')+
  theme(#axis.text.x = element_blank(),
    panel.spacing.x = unit(0.4,'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4,4,4,4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5,face = 'bold',size = 16)
  )+
  guides(fill='none',col='none')

#ggsave(filename = 'sigIDs_subtype_repstress.pdf',width = 10,height =7,device = cairo_pdf)  


# Fig. 3e -----------------------------------------------------------------
pdata <- tgct_activity_all %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  left_join(wgs_data %>% select(Tumor_Barcode,Subtype)) %>% 
  mutate(name=factor(name,levels=colnames(tgct_activity_all)[-1])) %>% 
  filter(str_starts(name,'SBS'))

pdata %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(log2(value+1)~Subtype,data=.)))

my_comparisons <-  list(c('Seminoma','Non-Seminoma'))

stat.test <- pdata %>% 
  mutate(value=log2(value+1)) %>% 
  group_by(name) %>% 
  rstatix::wilcox_test(value~ Subtype, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Subtype") %>%
  mutate(myformatted.p = if_else(p<0.01,sprintf("P = %.2e",p),sprintf("P = %.2f",p))) %>% 
  mutate(Subtype = my_comparisons[[1]][1]) %>% 
  mutate(col=if_else(p<0.05,ncicolpal[1],'black'))

stat.test

pdata %>% 
  ggplot(aes(Subtype,log2(value+1),fill=Subtype))+
  #geom_boxplot()+
  facet_wrap(~name,scales = 'free_y',nrow = 1)+
  geom_violin(trim=T,size=0.2)+
  geom_boxplot(width=0.15, fill="gray95",color='black',outlier.shape = 21,outlier.size = 0.8)+
  labs(x=NULL,y='log2(Mutation+1)')+
  scale_fill_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_text_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'XY',ticks = T)+
  panel_border(color = 'black')+
  theme(axis.text.x = element_blank(),
        panel.spacing.x = unit(0.4,'cm'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(4,4,4,4),
        legend.position = 'bottom',
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(hjust = 0.5,face = 'bold',size = 16))+
  stat_pvalue_manual(data = stat.test, col='col',label = "myformatted.p")+
  scale_color_manual(values = c(ncicolpal[1],'black'))+
  guides(color='none')

#ggsave(filename = 'SBS_subtype_comparisons.pdf',width = 9,height =6,device = cairo_pdf)  



# Fig. 3f -----------------------------------------------------------------
pdata <- tgct_activity_all %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  left_join(wgs_data %>% select(Tumor_Barcode,Subtype)) %>% 
  mutate(name=factor(name,levels=colnames(tgct_activity_all)[-1])) %>% 
  filter(str_starts(name,'ID'))

pdata %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(log2(value+1)~Subtype,data=.)))

my_comparisons <-  list(c('Seminoma','Non-Seminoma'))

stat.test <- pdata %>% 
  mutate(value=log2(value+1)) %>% 
  group_by(name) %>% 
  rstatix::wilcox_test(value~ Subtype, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Subtype") %>%
  mutate(myformatted.p = if_else(p<0.01,sprintf("P = %.2e",p),sprintf("P = %.2f",p))) %>% 
  mutate(Subtype = my_comparisons[[1]][1]) %>% 
  mutate(col=if_else(p<0.05,ncicolpal[1],'black'))

stat.test

pdata %>% 
  ggplot(aes(Subtype,log2(value+1),fill=Subtype))+
  #geom_boxplot()+
  facet_wrap(~name,scales = 'free_y',nrow = 1)+
  geom_violin(trim=T,size=0.2)+
  geom_boxplot(width=0.15, fill="gray95",color='black',outlier.shape = 21,outlier.size = 0.8)+
  labs(x=NULL,y='log2(Mutation+1)')+
  scale_fill_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_text_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'XY',ticks = T)+
  panel_border(color = 'black')+
  theme(axis.text.x = element_blank(),
        panel.spacing.x = unit(0.4,'cm'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(4,4,4,4),
        legend.position = 'bottom',
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(hjust = 0.5,face = 'bold',size = 16))+
  stat_pvalue_manual(data = stat.test, col='col',label = "myformatted.p")+
  scale_color_manual(values = c(ncicolpal[1],'black'))+
  guides(color='none')

#



