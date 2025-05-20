set_wd()
libztw()

load('wgs_data.RData')
load('ngspurity_data.RData',verbose = T)

### adjust ccf and clonality
BBsolution4 <- BBprofile %>% select(Tumor_Barcode,Battenberg) %>% unique()

DP_info_data <- NULL
for(i in 1:dim(BBsolution4)[1]){
  barcode <- BBsolution4$Tumor_Barcode[i]
  bbtmp <- BBsolution4$Battenberg[i]
  rdsfile <- paste0('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/NGSpurity/Result/',barcode,'/',bbtmp,'/NGSpurity/',barcode,'_mdata2.RDS')
  DP_info_tmp <- read_rds(rdsfile) %>% mutate(Tumor_Barcode=barcode,Likelihood=as.numeric(Likelihood))
  DP_info_data <- bind_rows(DP_info_data, DP_info_tmp)
}

Cluster_info <- DP_info_data %>% filter(!is.na(Cluster))%>% select(Tumor_Barcode,Cluster,Cluster_CCF,Cluster_nMut) %>% unique()

# adjust the ccf
Cluster_info <- Cluster_info %>%left_join(BBsolution4)

Clone_Cluster_info <- Cluster_info %>% filter(Cluster_CCF>0.95, Cluster_CCF < 1.05) %>% group_by(Tumor_Barcode) %>% arrange(Cluster_CCF) %>% slice(1) %>% ungroup() %>% select(Tumor_Barcode,Minmal_Clone_Cluster_CCF=Cluster_CCF)


tmpdata <- Cluster_info %>%
  filter(!Tumor_Barcode %in% Clone_Cluster_info$Tumor_Barcode) %>%
  group_by(Tumor_Barcode) %>%
  summarise(maxccf=max(Cluster_CCF),minccf=min(Cluster_CCF),n=n()) %>%
  arrange(n,desc(maxccf))

Clone_Cluster_info <- tmpdata %>% filter(n==1,maxccf>0.8) %>% select(Tumor_Barcode,Minmal_Clone_Cluster_CCF=maxccf) %>% bind_rows(Clone_Cluster_info)

Clone_Cluster_info <- tmpdata %>% filter(n==2,maxccf>0.8,minccf<0.6) %>% select(Tumor_Barcode,Minmal_Clone_Cluster_CCF=maxccf) %>% bind_rows(Clone_Cluster_info)

Clone_Cluster_info <- tmpdata %>% filter(n==3,maxccf>0.8) %>% select(Tumor_Barcode,Minmal_Clone_Cluster_CCF=maxccf) %>% bind_rows(Clone_Cluster_info)

Clone_Cluster_info <- Cluster_info %>% filter(Tumor_Barcode %in% (tmpdata %>% filter(n==4) %>% pull(Tumor_Barcode))) %>% filter(Cluster_CCF>0.85) %>% arrange(Cluster_CCF) %>% slice(1) %>% select(Tumor_Barcode,Minmal_Clone_Cluster_CCF=Cluster_CCF) %>% bind_rows(Clone_Cluster_info)
# tmpdata %>% filter(n==2,!(maxccf>0.8 & minccf<0.6)) %>% left_join(BBsolution4 %>% select(Tumor_Barcode,Battenberg)) %>% write_csv('check_clonality.csv',col_names = T)
# tmpdata %>% filter(n==3) %>% left_join(BBsolution4 %>% select(Tumor_Barcode,Battenberg)) %>% write_csv('check_clonality.csv',col_names = T,append = T)
# 
# ## load manually checked clonality
# tmpdata2 <- readxl::read_xlsx('check_clonality.xlsx',sheet = 1,col_names = T)
# 
# tmpdata3 <- Cluster_info %>% filter(Tumor_Barcode %in% (tmpdata %>% filter(n==3) %>% pull(Tumor_Barcode))) %>% group_by(Tumor_Barcode) %>% arrange(Cluster_CCF) %>% slice(2) %>% select(Tumor_Barcode,midccf=Cluster_CCF) %>% ungroup()
# 
# Clone_Cluster_info <- tmpdata2 %>% left_join(tmpdata3) %>%
#   mutate(Minmal_Clone_Cluster_CCF = if_else(Clonality == 'minccf',minccf,if_else(Clonality == 'maxccf', maxccf, midccf))) %>%
#   select(Tumor_Barcode,Minmal_Clone_Cluster_CCF) %>%
#   bind_rows(Clone_Cluster_info)

#Clone_Cluster_info <- Cluster_info %>% filter(Tumor_Barcode %in% (tmpdata %>% filter(n>3) %>% pull(Tumor_Barcode))) %>% filter(Cluster_CCF >0.8) %>%  group_by(Tumor_Barcode) %>% arrange(Cluster_CCF) %>% slice(1) %>% ungroup() %>% select(Tumor_Barcode,Minmal_Clone_Cluster_CCF=Cluster_CCF) %>%  bind_rows(Clone_Cluster_info)


Cluster_info <- Cluster_info %>% left_join(Clone_Cluster_info) %>% mutate(Clone = if_else(Cluster_CCF >= Minmal_Clone_Cluster_CCF, 'Y','N')) %>%
  select(Tumor_Barcode, Battenberg,Cluster,Cluster_CCF,Cluster_nMut,Minmal_Clone_Cluster_CCF,Clone)

DP_info_data <-  DP_info_data %>% select(-Clone, -Cluster_CCF) %>% left_join(Cluster_info) %>% select(Tumor_Barcode,Battenberg,everything()) #mutate(CCF = if_else(is.na(adjust_purity),CCF,CCF/adjust_purity))

save(DP_info_data,Cluster_info, file='DP_info_data.RData')
save(Cluster_info, file='Cluster_info.RData')

