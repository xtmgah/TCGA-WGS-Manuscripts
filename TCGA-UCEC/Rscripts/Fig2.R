#Fig2----

libztw()
pdfhr()
pdfhr2()
myggstyle()
source('./Rcode/Sherlock_functions.R')
source("./Rcode/ztw_function.R")
source("./Rcode/ztw.R")
source("./Rcode/Sigvisualfunc.R")
set_wd()
library(rstatix)
library(tidyverse)
library(hrbrthemes)
library(ggsci)
library(ggpubr)
library(scales)
library(ggrepel)
library(valr)
library(cowplot)
library(grid)

#Fig2a----
load('./Rdata/wgsdata.RData')
load('./Rdata/ucec_subtype_data.RData')
load('./Rdata/ngspurity_data.RData',verbose = T)
trafic <- read_delim("./all_trafic.txt",delim = "\t",col_names = T) %>% rename(Tumor_Barcode=Barcode) %>% filter(Tumor_Barcode %in% subtype_data$Tumor_Barcode) %>% left_join(subtype_data)

#### count ###### 
trafic <- trafic %>% filter(FILTER=="PASS") %>% mutate(TYPE2=case_when(
  CLASS != "L1" ~ "Others",
  TYPE == "TD0" ~ "Solo-L1",
  TYPE == "TD1" ~ "Partnered-3'transduction",
  TYPE == "TD2" ~ "Orphan-3'transduction",
  TYPE == "PSD" ~ "Processed-pseudogene",
  TRUE ~ "Others"
)) %>% mutate(tmp=TYPE,TYPE=TYPE2,TYPE2=tmp) %>% mutate(TYPE=factor(TYPE,levels = (c("Solo-L1","Partnered-3'transduction","Orphan-3'transduction","Processed-pseudogene","Others")))) %>% select(-tmp)
#table(trafic$FILTER,trafic$CLASS)
#table(trafic$FILTER,trafic$SCORE)
#table(trafic$FILTER,trafic$TYPE)

#### length distribution ####

trafic$LEN <- as.integer(trafic$LEN)

mean(trafic$LEN,na.rm = TRUE)

trafic$LEN2 <- cut_width(trafic$LEN/1000,width = 0.1,boundary = 0)

trafic %>% filter(!is.na(LEN2),CLASS=="L1") %>% mutate(TYPE=factor(TYPE,levels = c("Solo-L1","Partnered-3'transduction","Orphan-3'transduction")))%>% ggplot(aes(LEN2,fill=TYPE))+geom_bar(width = 0.7,position = position_stack(reverse = TRUE))+theme_ipsum_rc(base_size = 12,axis = "xy",axis_title_just = "m",axis_title_size = 15)+theme(axis.text.x = element_text(size = 10,angle = 40,hjust = 1,vjust = 1))+scale_fill_manual(values = typecol,drop=FALSE)+labs(x="\nSize of insertions (kb)",y="Number of insertions\n")+geom_vline(xintercept = 12,col="gray30",linetype=2)+geom_text_repel(aes(y=9,label=label),data=data_frame(LEN2="(6,6.1]",label="Full-length\n   Solo-L1    ",TYPE="Solo-L1"),force = 100,nudge_y = 0.01,arrow = arrow(length=unit(0.01,"npc")))

### overview of retrotransposition events identified by TraFiC

trafic %>% select(Tumor_Barcode) %>% unique() %>% 
  left_join(
    trafic %>% select(Tumor_Barcode,CLASS,TYPE) %>% filter(CLASS=="L1") %>% count(Tumor_Barcode) %>% rename(`Total L1`=n)
  ) %>% left_join(
    trafic %>% select(Tumor_Barcode,CLASS,TYPE) %>% filter(TYPE=="Solo-L1") %>% count(Tumor_Barcode) %>% rename(`Solo-L1`=n)
  ) %>% left_join(
    trafic %>% select(Tumor_Barcode,CLASS,TYPE) %>% filter(str_detect(TYPE,"Partnered")) %>% count(Tumor_Barcode) %>% rename(`Partnered`=n)
  ) %>% left_join(
    trafic %>% select(Tumor_Barcode,CLASS,TYPE) %>% filter(str_detect(TYPE,"Orphan")) %>% count(Tumor_Barcode) %>% rename(`Orphan`=n)
  ) %>% left_join(
    trafic %>% select(Tumor_Barcode,CLASS,TYPE) %>% filter(CLASS=="Alu") %>% count(Tumor_Barcode) %>% rename(`Total Alu`=n)
  ) %>% left_join(
    trafic %>% select(Tumor_Barcode,CLASS,TYPE) %>% filter(!(CLASS %in% c("L1","Alu"))) %>% count(Tumor_Barcode) %>% rename(`Total Others`=n)
  ) %>% left_join(
    trafic %>% select(Tumor_Barcode,CLASS,TYPE) %>%  count(Tumor_Barcode) %>% rename(`Total Insertion`=n)
  ) #%>% 
# write_excel_csv('overview_retrotransposition_events.csv',col_names = T,na = "0")


### Somatic 3' transductions originate from a limited repertorie of L1 source element

cytoband <- read_delim('./cytoBand.txt.gz',delim = "\t",col_names = F)
colnames(cytoband) <- c('chrom','start','end','cytoband','location')
cytoband <- cytoband  %>% filter(chrom %in% paste0('chr',c(1:22,"X"))) %>% mutate(chrom=factor(chrom,levels = paste0('chr',c(1:22,"X"))))
cytoclor <- c(gneg="#FFFFFF",
              gpos25="#C8C8C8",
              gpos33="#D2D2D2",
              gpos50="#C8C8C8",
              gpos66="#A0A0A0",
              gpos75="#828282",
              gpos100="#000000",
              gpos="#000000",
              stalk="#647FA4", #repetitive areas
              acen="#D92F27", #centromeres
              gvar="#DCDCDC")

p1 <- cytoband %>% mutate(ystart=0,yend=1) %>% ggplot()+geom_rect(aes(xmin=start,xmax=end,ymin=ystart,ymax=yend,fill=location),col="black",size=0.001)+facet_grid(.~chrom,scales = "free_x",space="free_x",switch = "x")+scale_fill_manual(values = cytoclor)+theme(legend.position = "none",panel.spacing = unit(0.5, "lines"),panel.border = element_rect(colour = "black", fill=NA, size=1),strip.text=element_text(hjust = 0.5,size=12),axis.title = element_blank(),axis.text = element_blank(),strip.background = element_blank(),panel.grid = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))


### sampleset
sampleset <- wgsdata %>% pull(Tumor_Barcode)


library(valr)
trdata <- trafic %>% filter(Tumor_Barcode %in% sampleset) %>% filter(SRC!=".") %>% select(Tumor_Barcode,SRC,SRCTYPE,SRCID) %>% separate(SRC,c("chrom","start","end")) %>% mutate(start=as.integer(start),end=as.integer(end)) %>% select(chrom,start,end,Tumor_Barcode,SRCTYPE,SRCID) 

trdata <- bed_intersect(trdata,cytoband,suffix = c(".x","")) %>% group_by(chrom,start,end,cytoband) %>% summarise(n=n_distinct(Tumor_Barcode.x)) %>% ungroup() %>% mutate(n=n/length(sampleset)) %>% mutate(BP=(start+end)/2) %>% select(chrom,BP,cytoband,n)

trdata <- trdata %>% bind_rows(
  cytoband %>% select(chrom,BP=start,cytoband) %>% mutate(n=0),
  cytoband %>% select(chrom,BP=end,cytoband) %>% mutate(n=0)
)%>% mutate(chrom=factor(chrom,levels = paste0('chr',c(1:22,"X"))))

trdata2 <- trdata %>% filter(n>0.01) 
p2 <- trdata %>% ggplot(aes(BP,n,col=n>0))+geom_bar(stat = "identity",size=1)+facet_grid(.~chrom,scales = "free_x",space="free_x",switch = "x",drop = FALSE)+theme(legend.position = "none",panel.spacing = unit(0.5, "lines"),strip.text=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),strip.background = element_blank(),panel.grid = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),breaks = pretty_breaks())+labs(y="Average transductions per sample")+geom_text_repel(data=trdata2,aes(BP,n,label=cytoband),size=4)+scale_color_manual(values=c("white","green4"))+theme(plot.margin=margin(b=-0.1,t = 1,unit="cm"))

#force = 10,arrow = arrow(length=unit(0.1,"npc")
### combined plot ###
ggarrange(p2,p1,ncol=1,nrow = 2,align = 'v',heights = c(8,1))
ggsave("./Rscripts/Figures/Fig2a_1.png",width = 28,height = 6)


#### list the recurrent #######
library(circlize)

for(i in 1:dim(trdata2)[1]){
  
  chrp <- trdata2$chrom[i]
  cytobandp <- trdata2$cytoband[i]
  cytobandp <- str_remove(cytobandp,"^[0-9XY]*")
  
  trdata <- trafic %>% filter(SRC!=".") %>% select(Tumor_Barcode,SRC,SRCTYPE,SRCID) %>% separate(SRC,c("chrom","start","end")) %>% mutate(start=as.integer(start),end=as.integer(end)) %>% select(chrom,start,end,Tumor_Barcode,SRCTYPE,SRCID) 
  
  srcdata <- bed_intersect(trdata,cytoband,suffix = c(".x","")) %>% filter(chrom==chrp,cytoband==cytobandp) %>% mutate(SRC=paste(chrom,start.x,end.x,sep = "_")) %>% pull(SRC)
  
  bedall <- trafic %>% filter(SRC %in% srcdata) %>% select(CHROM,POS,SRC,SRCTYPE,Tumor_Barcode) %>% separate(SRC,c("chrom","start","end"),convert = T) 
  
  print(bedall)
  
  bed1 <- bedall %>% select(chr=CHROM,start=POS,value1=SRCTYPE) %>% mutate(end=start,start=start-1000000,end=end+1000000) %>% select(chr,start,end,value1) %>% as.data.frame() 
  bed2 <- bedall %>% select(chr=chrom,start,end,value2=SRCTYPE) %>% mutate(end=end,start=start-1000000,end=end+1000000) %>% select(chr,start,end,value2)%>% as.data.frame()
  scols <- c('#e41a1c','#4daf4a')
  names(scols) <- c("SOMATIC","GERMLINE")
  pdf(paste0('./Rscripts/Figures/Hotpot/',chrp,cytobandp,"_Somatic_transductions.pdf"),height = 6,width = 6)
  circos.initializeWithIdeogram()
  circos.genomicLink(bed1, bed2, col =scols[bed1$value1] , border = NA,directional = -1)
  legend("bottomleft", pch = 15, legend = names(scols), col = scols)
  dev.off()
}


#Fig2b----
load('./Rdata/wgsdata.RData')
load('./Rdata/ucec_subtype_data.RData')

chrs=c(1:22,"X")
studyid <- 'POLE'#studyid <- 'MSI'/studyid <- 'CN_LOW'/studyid <- 'CN_HIGH'

samplelist <- subtype_data %>% filter(Subtype==studyid) %>% pull(Tumor_Barcode)

## load final result 
svdata <- read_delim('./manta_meerkat_union_window50_highquality.txt',delim = '\t',col_names = T)
svdata <- svdata %>% filter(SAMPLE %in% subtype_data$Tumor_Barcode) %>% filter(!is.na(CHROM1))

data <- svdata %>% select(barcode=SAMPLE,chrA=CHROM1,posA=POS1,chrB=CHROM2,posB=POS2)
data <- data %>% filter(chrA %in% chrs, chrB %in% chrs) %>% filter(barcode %in% samplelist)

## make reference bin ###
hg38 <- read_delim('./Homo_sapiens_assembly38.fasta.fai',delim = '\t',col_names = F) %>% mutate(X1=str_remove(X1,"chr"))
hg38 <- hg38 %>% 
  select(chr=X1,len=X2) %>%
  slice(1:24) %>% 
  mutate(start=NA_integer_,end=NA_integer_)
hg38$start=cumsum(c(0,hg38$len[1:23]))+1
hg38$end=cumsum(hg38$len[1:24])

hg38bins <- NULL

bin=5000000
for(i in 1:24){
  size=hg38$len[i]
  chr=hg38$chr[i]
  startx=hg38$start[i]
  tmp <- tibble(chr=chr,start=seq(from = 1,to=size,by = bin),end=c(seq(from = 0,to=size,by = bin)[-1],size))
  tmp <- tmp %>% mutate(start2=start+startx-1,end2=end+startx-1)
  hg38bins <- bind_rows(hg38bins,tmp)
}

hg38 <- hg38 %>% rename(start2=start,end2=end) %>% filter(chr %in% chrs)

### 1D data ### 
data1d <- bind_rows(
  data %>% select(chrom=chrA,start=posA),
  data %>% select(chrom=chrB,start=posB)
) %>% 
  mutate(end=start+1)

hg38bins_bed <- hg38bins %>% select(chrom=chr,start,end) 

data1d <- bed_intersect(hg38bins_bed,data1d,suffix=c('','.y')) %>% count(chrom,start,end) %>% rename(chr=chrom)

data1d <- hg38bins %>% left_join(data1d) %>%  mutate(n=if_else(is.na(n),0L,n))

data1d <- data1d %>% mutate(n=n/length(samplelist))

data_chr= hg38 %>% mutate(pos=(start2+end2)/2) %>% select(chr,pos)

ymax <- 2.8

# normalize to sample level

p1 <- data1d %>% 
  filter(chr %in% chrs) %>% 
  mutate(pos=(start2+end2)/2) %>% 
  ggplot(aes(pos,n))+
  geom_area(fill="darkred")+
  geom_line(size=0.25)+
  labs(x="",y="Breakpoints per 5 Mbp per sample",title = "Chromosome")+
  scale_x_continuous(breaks = hg38$end2,expand = c(0,0))+
  scale_y_continuous(expand = c(0,0),breaks = pretty_breaks(),limits = c(0,3))+
  theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,grid="XY",ticks = T,axis = FALSE)+
  theme(axis.ticks.x=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.line = element_line(color = "black",size=0.5),plot.title = element_text(hjust = 0.5,face = "plain"))+
  geom_text_repel(data=data_chr,aes(pos,ymax,label=chr),force = 0,size=3)

ggsave(filename = paste0('./Rscripts/Figures/',studyid,'Fig2b_4.png'),plot = p1,width = 9,height = 2.5)

#Fig2c----
#### count ###### 
load('./Rdata/ecDNA.RData',verbose = T)
## AA count distribution across subtypes
ecdna %>% 
  count(Subject,Classification) %>% 
  left_join(subtype_data) %>%
  count(Subtype,Classification,n) %>% 
  ggplot(aes(nn,Classification,fill=n))+
  geom_col(col='black',linewidth=0.4,width = 0.8)+
  facet_wrap(~Subtype,nrow=1)+
  scale_fill_viridis_c(direction = -1,breaks=pretty_breaks())+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_text_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'XY',ticks = T)+
  panel_border(color = 'black',size=0.35)+
  theme(#axis.text.x = element_text(family = 'Roboto',size = 12,angle = 90,vjust = 0.5,hjust = 1),
    panel.spacing.x = unit(0.2,'cm'),
    axis.ticks.y = element_blank(),
    plot.margin = margin(4,4,0,4),
    legend.position = 'bottom',
    legend.key.width = unit(x = 2,'cm'),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5,face = 'bold',size = 16)
  )+
  labs(x='Number of tumors',fill='Number of amplicons per tumor\n',y=NULL)

ggsave(filename = './Rscripts/Figures/Fig2c.png',width = 12,height = 3)

#Fig2d----
load('./Rdata/ecDNA.RData',verbose = T)
#overall frequency
ecdna %>% group_by(Classification) %>% count(Subject) %>% count(Classification) %>% ungroup() %>% mutate(Frequency=n/dim(wgsdata)[1]) %>% select(Classification,n,Frequency) %>% write_clip()

# ecdna %>% group_by(Classification,Sample_Type) %>% count(Subject)%>% ungroup() %>% count(Classification,Sample_Type)  %>%  left_join(wgsdata %>% count(Sample_Type,name = 'total')) %>% mutate(Frequency=n/total) %>% select(Sample_Type,Classification,n,Frequency) %>% arrange(Sample_Type,Classification) %>% write_clip()

#with oncogene
ecdna %>% filter(Oncogenes != '[]') %>% group_by(Classification) %>% count(Subject) %>% count(Classification) %>% ungroup() %>% mutate(Frequency=n/dim(wgsdata)[1]) %>% select(Classification,n,Frequency) %>% write_clip()
# ecDNA context ####
ecdna_context_colors <- c('Two-foldback'='#4292c6','Simple circular simple background'='#08519c','Simple circular complex background'='#fb6a4a','Heavily rearranged unichromosomal'='#cb181d','Heavily rearranged multichromosomal'='#a50f15','BFB-like'='#67000d','Unknown'='gray50')
ecdna_context_colors2 <- c('Origin unclear'='gray50', 'Simple ecDNA'='#007BBD','Chromothripsis-associated ecDNA'='#BB0E3D')
ecdna_context_map <- tibble(`ecDNA context`=names(ecdna_context_colors),ecDNA_context2=names(ecdna_context_colors2)[c(2,2,3,3,3,3,1)])
tdata <- ecdna %>%
  filter(Classification == 'ecDNA') %>% 
  count(`ecDNA context`) %>% 
  rename(ecDNA_context=`ecDNA context`)

PieDonut_ztw(data = tdata,aes(pies=ecDNA_context,count=n),mainCol = ecdna_context_colors,start=3*pi/2,showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.02,showDonutName = T,pieLabelSize=3.5,titlesize = 5,r0=0.3,title='AS_N', showPieName = TRUE,family = 'Roboto Condensed', labelpositionThreshold =0.02)

ggsave(filename = './Rscripts/Figures/Fig2d_1.png',width = 6,height = 6)


pdata <- ecdna %>% 
  filter(Classification == 'ecDNA') %>% 
  mutate(Oncogenes = str_remove_all(Oncogenes,"[\\[\\'\\] ]")) %>% 
  separate_rows(Oncogenes,sep = ',') %>% 
  filter(Oncogenes != '') %>% 
  select(Sample_Type,Oncogenes,Subject,`ecDNA context`) %>% 
  unique() %>% 
  count(Sample_Type,Oncogenes,`ecDNA context`,sort=T) %>% 
  left_join(ecdna_context_map)

gene_order <- pdata %>% group_by(Oncogenes) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% filter(n>2) %>% pull(Oncogenes)

p1 <- pdata %>% 
  mutate(Oncogenes = factor(Oncogenes,levels = gene_order)) %>% 
  filter(Oncogenes %in% gene_order) %>% 
  ggplot(aes(Oncogenes,n,fill=`ecDNA context`)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = ecdna_context_colors)+
  scale_y_continuous(breaks = pretty_breaks(n=5),expand = expansion(add = c(0,0)))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = "Y",axis = "XY",axis_col = 'black')+
  theme(axis.text.x=element_text(size = 12,angle = 90,vjust = 0.5,hjust = 1),plot.margin = margin(4,4,4,4),panel.spacing = unit('0.1','cm'),strip.text.y = element_text(vjust = 0,hjust = 0.5),legend.position = 'right')+
  labs(x=NULL,y = 'Number of subjects')

## add oncogene co-occurence plot
library(ggforce)
oncogene_features <- ecdna %>% 
  filter(Classification == 'ecDNA') %>% 
  select(`Feature ID`,Oncogenes)

# --- Parse genes per row from your original tibble ---
# oncogene_features: columns `Feature ID` (chr), `Oncogenes` (chr like "['MYC','PVT1']")
parsed <- oncogene_features %>%
  mutate(genes = stringr::str_extract_all(Oncogenes, "'([^']+)'")) %>%
  mutate(genes = purrr::map(genes, ~ stringr::str_replace_all(.x, "'", "")))

# --- Pair co-occurrence table (oncogene1, oncogene2, N_pairs) ---
oncogene_pairs <- parsed %>%
  transmute(`Feature ID`, genes) %>%
  mutate(pairs = purrr::map(genes, ~ {
    g <- unique(.x)
    if (length(g) > 1) combn(sort(g), 2, simplify = FALSE) else NULL
  })) %>%
  tidyr::unnest(pairs, keep_empty = FALSE) %>%
  mutate(
    oncogene1 = purrr::map_chr(pairs, 1),
    oncogene2 = purrr::map_chr(pairs, 2)
  ) %>%
  count(oncogene1, oncogene2, name = "N_pairs")

## only included the genes in order
oncogene_pairs <- oncogene_pairs %>% filter(oncogene1 %in% gene_order,oncogene2 %in% gene_order)

# --- Unique Feature ID counts per gene (for sizing & ordering) ---
gene_counts <- parsed %>%
  select(`Feature ID`, genes) %>%
  tidyr::unnest(genes, names_repair = "minimal") %>%
  distinct(`Feature ID`, genes) %>%
  count(genes, name = "n_features")  # number of unique Feature IDs per gene

gene_counts <- gene_counts %>% filter(genes %in% gene_order)

# --- Linear layout: order genes by n_features (desc), then alphabetically ---
# layout <- gene_counts %>%
#   arrange(desc(n_features), genes) %>%
#   mutate(x = row_number(), y = 0) %>%
#   rename(oncogene = genes)

layout <- gene_counts %>%
  mutate(
    order = match(genes, gene_order)  # get index of each gene in the custom vector
  ) %>%
  #arrange(order) %>%                   # follow that exact order
  arrange(is.na(order), order) %>% 
  mutate(x = row_number(), y = 0) %>%
  rename(oncogene = genes)


# --- Build Bezier edges (downward arcs) ---
edges <- oncogene_pairs %>%
  left_join(layout, by = c("oncogene1" = "oncogene")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(layout, by = c("oncogene2" = "oncogene")) %>%
  rename(x2 = x, y2 = y) %>%
  mutate(
    cx = (x1 + x2) / 2,          # midpoint in x
    cy = -abs(x2 - x1) / 3       # NEGATIVE height => curve below axis
  )

bezier_edges <- edges %>%
  select(oncogene1, oncogene2, N_pairs, x1, y1, cx, cy, x2, y2) %>%
  pmap_dfr(\(oncogene1, oncogene2, N_pairs, x1, y1, cx, cy, x2, y2) {
    tibble(
      oncogene1 = oncogene1,
      oncogene2 = oncogene2,
      N_pairs   = N_pairs,
      x         = c(x1, cx, x2),
      y         = c(y1, cy, y2),
      group     = paste(oncogene1, oncogene2, sep = "_")
    )
  })

# --- Points (nodes) with sizes = n_features ---
nodes <- layout  # has oncogene, x, y
nodes <- nodes %>% left_join(gene_counts, by = c("oncogene" = "genes",'n_features'='n_features'))

# --- Plot ---
p2 <- ggplot() +
  # downward arcs; use linewidth for edge thickness so it won't clash with point size
  geom_bezier2(
    data = bezier_edges,
    aes(x = x, y = y, group = group, linewidth = N_pairs,colour = N_pairs)
    #alpha = 0.6
  ) +
  scale_linewidth_continuous(range = c(0.3, 3), name = "Pair count") +
  scale_color_viridis_c(direction = -1)+
  scale_x_continuous(expand = expansion(add = c(0.6,0.6)))+
  # points sized by unique Feature IDs
  geom_point(
    data = nodes,
    aes(x = x, y = y, size = n_features),
    shape = 21, stroke = 0.4, fill = "white"
  ) +
  scale_size_continuous(range = c(1, 6), name = "Unique Feature IDs") +
  # gene labels just below the axis
  # geom_text(
  #   data = nodes,
  #   aes(x = x, y = y + 0.2, label = oncogene),angle = 90,
  #   vjust = 0.5, hjust=0,size = 3.8, fontface = "bold"
  # ) +
  # clean layout
  coord_cartesian(clip = "off") +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = FALSE,grid = FALSE,axis = FALSE,axis_col = 'black')+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),plot.margin = margin(-4,4,4,4),panel.spacing = unit('0.1','cm'),strip.text.y = element_text(vjust = 0,hjust = 0.5),legend.position = c(1.22,1))+
  labs(x='Oncogene co-occurrence ',y = NULL,color='Pair count',size='Pair count')

p <- plot_grid(p1,p2,align = 'v',axis = 'lr',ncol = 1,rel_heights = c(5,2))

ggsave(filename = './Rscripts/Figures/Fig2d_2.png',plot = p,width = 8,height = 10)

#Fig2e----
# TE vs ecDNA ####
load('./Rdata/trafic.RData',verbose = T)
load('./Rdata/ecDNA.RData',verbose = T)

# TE count vs ecDNA count
tdata <- subtype_data %>% 
  left_join(
    trafic %>% count(Tumor_Barcode,name = 'TE') #filter(SRCID == '22q12.1')
  ) %>% 
  left_join(
    ecdna %>% count(Tumor_Barcode=Barcode,Classification) %>% pivot_wider(names_from = Classification,values_from = n)
  ) %>% 
  mutate(across(-c(Subject:Subtype), ~ replace_na(., 0))) %>% 
  mutate(TE=log2(TE+1)) 


# ecdna vs TE
plot_group_compare(
  tdata       = tdata %>% mutate(ecDNA = if_else(ecDNA>0,'ecDNA+','ecDNA-')), # %>% filter(Subtype=='CN_HIGH'),
  value       = TE,
  group       = "ecDNA",
  ylab        = "log2(TE count +1)",
  pcol_sig    = ncicolpal[1],palette = c("gray20","#D62728FF"),
  output_file = "TE_ecDNA.pdf",   # omit to just get the plot object
  width       = 3,               # optional; if omitted, auto-size by #groups
  height      = 5,hide_ns = F
)

tdata %>% 
  #filter(Subtype=='CN_HIGH') %>% 
  mutate(ecDNA = if_else(ecDNA>2,2,ecDNA)) %>% 
  ggplot(aes(as.factor(ecDNA),TE))+
  ggbeeswarm::geom_quasirandom(aes(fill=as.factor(ecDNA)),pch = 21, size = 3, width = 0.4, color = 'black', stroke = 0.25, alpha = 1) +
  geom_boxplot(width = 0.5, fill = NA, color = "gray20", outlier.shape = NA, size = 0.6) +
  geom_smooth(aes(as.integer(as.factor(ecDNA)),TE),method = 'lm')+
  scale_fill_manual(values = c('gray80',"#FDE725FF","#6DCD59FF"))+
  scale_x_discrete(breaks = c(0,1,2),labels = c("0","1",'2+'))+
  labs(x='Number of ecDNA amplicons',y='log2(TE count +1)')+
  theme_ipsum_rc(
    axis_text_size = 14,
    axis_title_just = 'm',
    axis_title_size = 16,
    grid = 'XY',
    ticks = T
  ) +
  panel_border(color = 'black') +
  theme(
    #axis.text.x = element_blank(),
    panel.spacing.x = unit(0.1, 'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  guides(fill = "none")

ggsave(filename = './Rscripts/Figures/Fig2e.png', width = 3,height = 5)

#Fig2f----
load("./Rdata/ShatterSeek.RData",verbose = T)

tdata0 <- subtype_data %>% 
  left_join(
    chromothripsis_summary %>% 
      filter(HCany) %>% 
      select(Subtype,Tumor_Barcode) %>% 
      unique() %>% 
      mutate(HC="Y")
  ) %>% left_join(
    chromothripsis_summary %>% 
      filter(LC1) %>% 
      select(Subtype,Tumor_Barcode) %>% 
      unique() %>% 
      mutate(LC='Y')
  ) %>% 
  mutate(LC=if_else(!is.na(HC),NA_character_,LC)) %>% 
  pivot_longer(-c(Subject:Subtype)) %>% 
  filter(!is.na(value))

tdata <- tdata0 %>%
  count(Subtype,name) %>% 
  left_join(
    subtype_data %>% count(Subtype,name = 'Total')
  ) %>% 
  mutate(Freq=n/Total) 


tdata %>% 
  mutate(name=factor(name,levels = c('LC','HC'),labels=c('Low','High'))) %>% 
  ggplot(aes(Subtype,Freq,fill=name))+
  geom_col(col='black',linewidth = 0.2)+
  scale_fill_manual(values = c('High'='#007BBD','Low'='#b2df8a'),)+
  labs(x=NULL,y='Percentage of samples with chromothripsis',fill='Confidence')+
  scale_x_discrete(drop=FALSE)+
  scale_y_continuous(breaks = pretty_breaks(n = 7),labels = percent_format(),expand = expansion(mult = c(0,0.05)))+
  theme_ipsum_rc(axis_text_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'Y',ticks = F)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        #panel.spacing.x = unit(0.4,'cm'),
        axis.ticks.x = element_blank(),
        plot.margin = margin(6,4,4,4),
        #legend.position = 'bottom',
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(hjust = 0.5,face = 'bold',size = 16))+
  panel_border(linetype = 1,size=0.5,colour = 'black')

ggsave(filename = './Rscripts/Figures/Fig2f.png',width = 3.5,height = 8)

#Fig2g----
load('./Rdata/ecDNA.RData',verbose = T)

list1 <- tdata0 %>% pull(Tumor_Barcode) %>% unique()
list2 <- ecdna %>% filter(Classification == 'ecDNA') %>% pull(Barcode) %>% unique()

tdata <- subtype_data %>% 
  mutate(chromothripsis=if_else(Tumor_Barcode %in% list1,'Yes','No')) %>% 
  mutate(ecDNA=if_else(Tumor_Barcode %in% list2,'Yes','No'))

tdata %>% filter(Subtype == 'CN_HIGH')%>% select(chromothripsis,ecDNA) %>% table() %>% fisher.test()

barplot_fisher(mdata0 = tdata,var1name = 'ecDNA',var2name = 'chromothripsis',var1lab = "ecDNA",var2lab = "Chromothripsis",var1_num = T,filename = 'chromothripsis_ecDNA_enrichment.pdf')
barplot_fisher(mdata0 = tdata %>% filter(Subtype == 'CN_HIGH'),var1name = 'ecDNA',var2name = 'chromothripsis',var1lab = "ecDNA",var2lab = "Chromothripsis",var1_num = T,filename = 'chromothripsis_ecDNA_enrichment_CN_HIGH.pdf')

# association with trafic
load('./Rdata/trafic.RData',verbose = T)

tdata2 <- subtype_data %>% 
  left_join(
    trafic %>% count(Tumor_Barcode)
  ) %>% 
  mutate(n=replace_na(n,0)) %>% 
  mutate(TE=log2(n+1)) %>% 
  mutate(chromothripsis=if_else(Tumor_Barcode %in% list1,'Chromothripsis+','Chromothripsis-'))


plot_group_compare(
  tdata       = tdata2, # %>% filter(Subtype=='CN_HIGH'),
  value       = TE,
  group       = "chromothripsis",
  ylab        = "log2(TE count +1)",
  pcol_sig    = ncicolpal[1],palette = c("gray20","#D62728FF"),
  output_file = "./Rscripts/Figures/Fig2g.png",   # omit to just get the plot object
  width       = 4,               # optional; if omitted, auto-size by #groups
  height      = 5,hide_ns = F
)
