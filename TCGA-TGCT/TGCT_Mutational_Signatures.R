set_wd()
libztw()
pdfhr2()
pdfhr()
load(file = 'wgs_data.RData')

load('~/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Database/Updated/signature_refsets.RData',verbose = T)
source('~/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Codes/Sigvisualfunc.R')


# Mutational Profile ------------------------------------------------------
tgct_sbs96_profile <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/TGCT_spmatrix/output/SBS/TGCT.SBS96.all') %>% pivot_longer(cols = -MutationType,names_to = 'Tumor_Barcode') %>% select(Tumor_Barcode, MutationType, Contribution=value)


tgct_id83_profile <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/TGCT_spmatrix/output/ID/TGCT.ID83.all') %>% pivot_longer(cols = -MutationType,names_to = 'Tumor_Barcode') %>% select(Tumor_Barcode, MutationType, Contribution=value)

tgct_dbs78_profile <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/TGCT_spmatrix/output/DBS/TGCT.DBS78.all') %>% pivot_longer(cols = -MutationType,names_to = 'Tumor_Barcode') %>% select(Tumor_Barcode, MutationType, Contribution=value)

# Mutational Signature Bootstrap Results ------------------------------------
msa_sbs96_activity <- read_csv(file = '/Volumes/data/Ref/MSA/MSA_TGCT/output/TGCT/SBS/output_tables/output_TGCT_SBS_mutations_table.csv',col_names = T) %>%
  rename(Tumor_Barcode=Sample)
msa_sbs96_activity <- msa_sbs96_activity %>% select(Tumor_Barcode,sort(colnames(msa_sbs96_activity)[-1])[c(1,3,4,2,5)])

msa_id83_activity <- read_csv(file = '/Volumes/data/Ref/MSA/MSA_TGCT/output/TGCT/ID/output_tables/output_TGCT_ID_mutations_table.csv',col_names = T) %>%
  rename(Tumor_Barcode=Sample)
msa_id83_activity <- msa_id83_activity %>% select(Tumor_Barcode,sort(colnames(msa_id83_activity)[-1]))

msa_dbs78_activity <- read_csv(file = '/Volumes/data/Ref/MSA/MSA_TGCT/output/TGCT/DBS/output_tables/output_TGCT_DBS_mutations_table.csv',col_names = T) %>%
  rename(Tumor_Barcode=Sample)
msa_dbs78_activity <- msa_dbs78_activity %>% select(Tumor_Barcode,sort(colnames(msa_dbs78_activity)[-1]))

# combined
tgct_activity_all <- left_join(msa_sbs96_activity, msa_id83_activity) %>% left_join(msa_dbs78_activity)
tgct_activity_all_obs <- tgct_activity_all %>% mutate(across(where(is.numeric),~ . > 50)) ## consistence with other

tgct_activity_all_ratio <- msa_sbs96_activity %>% 
  adorn_percentages(denominator = 'row') %>%
  as_tibble() %>% 
  left_join(msa_id83_activity %>% adorn_percentages(denominator = 'row') %>% as_tibble()) %>% 
  left_join(msa_dbs78_activity %>% adorn_percentages(denominator = 'row') %>% as_tibble())


# decompsition based on the msa result
#tmp <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/result_spextractor_cosmic3.4/TGCT.SBS96.all/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Solution_Stats/COSMIC_SBS96_Samples_Stats.txt')

orignal_genome <- tgct_sbs96_profile %>% pivot_wider(names_from = Tumor_Barcode,values_from = Contribution)
refit_signature <- signature_refsets %>% filter(Signature_set_name %in% c('COSMIC_v3.4_Signatures_GRCh38_SBS96')) %>% filter(Signature_name %in% names(msa_sbs96_activity)) %>% select(Signature_name,MutationType,Contribution) %>% pivot_wider(names_from = Signature_name,values_from = Contribution)
refit_genome <- msa_sbs96_activity %>% rename(Sample=Tumor_Barcode)
tgct_sbs96_msa_decompsite <- calculate_similarities(orignal_genomes = orignal_genome,signature = refit_signature,signature_activaties = refit_genome)

orignal_genome <- tgct_id83_profile %>% pivot_wider(names_from = Tumor_Barcode,values_from = Contribution)
refit_signature <- signature_refsets %>% filter(Signature_set_name %in% c('COSMIC_v3.4_Signatures_GRCh37_ID83')) %>% filter(Signature_name %in% names(msa_id83_activity)) %>% select(Signature_name,MutationType,Contribution) %>% pivot_wider(names_from = Signature_name,values_from = Contribution)
refit_genome <- msa_id83_activity %>% rename(Sample=Tumor_Barcode)
tgct_id83_msa_decompsite <- calculate_similarities(orignal_genomes = orignal_genome,signature = refit_signature,signature_activaties = refit_genome)


orignal_genome <- tgct_dbs78_profile %>% pivot_wider(names_from = Tumor_Barcode,values_from = Contribution)
refit_signature <- signature_refsets %>% filter(Signature_set_name %in% c('COSMIC_v3.4_Signatures_GRCh38_DBS78')) %>% filter(Signature_name %in% names(msa_dbs78_activity)) %>% select(Signature_name,MutationType,Contribution) %>% pivot_wider(names_from = Signature_name,values_from = Contribution)
refit_genome <- msa_dbs78_activity %>% rename(Sample=Tumor_Barcode)
tgct_dbs78_msa_decompsite <- calculate_similarities(orignal_genomes = orignal_genome,signature = refit_signature,signature_activaties = refit_genome)




# Saving objects ----------------------------------------------------------

sigcol=c('SBS1'="#4a9855",'SBS3'="#984ea3",'SBS5'="#01665e",'SBS18'="#ff7f00",'SBS87'="#BB0E3D",
         'ID1'="#0073C2FF",'ID2'="#EFC000FF",'ID4'="#868686FF",'ID5'="#CD534CFF",'ID6'="#8F7700FF",ID9="#003C67FF",
         DBS13="#a65628",DBS18="#5aae61"
) 


tgct_tmb <- bind_rows(
  tgct_sbs96_profile,
  tgct_id83_profile,
  tgct_dbs78_profile
) %>% 
  group_by(Tumor_Barcode) %>% 
  summarise(TMB=sum(Contribution)) 

save(tgct_sbs96_profile,tgct_id83_profile,tgct_dbs78_profile,msa_sbs96_activity,msa_id83_activity,msa_dbs78_activity,tgct_activity_all,tgct_activity_all_obs,tgct_activity_all_ratio,tgct_sbs96_msa_decompsite,tgct_id83_msa_decompsite,tgct_dbs78_msa_decompsite,sigcol,tgct_tmb,file='tgct_mutational_signatures.RData')



# Signature Probabilities -------------------------------------------------
library(signature.tools.lib)
load('DP_info_data.RData',verbose = T)


SBSchannels <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", 
                 "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", 
                 "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", 
                 "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", 
                 "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", 
                 "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", 
                 "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", 
                 "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", 
                 "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", 
                 "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", 
                 "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", 
                 "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", 
                 "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", 
                 "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", 
                 "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", 
                 "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", 
                 "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", 
                 "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", 
                 "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", 
                 "T[T>G]G", "T[T>G]T")

DP_info_data <- DP_info_data %>% 
  select(Tumor_Barcode,CHROM,POS,REF,ALT,Gene_Name,AAChange,Variant_Classification,ID,mutType) %>% 
  filter(!is.na(mutType)) %>% 
  rename(context=mutType)

signatures <- read.delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/result_spextractor_cosmic3.4/TGCT.SBS96.all/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures/COSMIC_SBS96_Signatures.txt') 
rownames(signatures) <- signatures$MutationType
signatures <- signatures[,-1]
signatures <- signatures[SBSchannels,]

tgct_sbs96_probability <- NULL

for(barcode in unique(DP_info_data$Tumor_Barcode)){
  
  sampleMutations <- DP_info_data %>% filter(Tumor_Barcode == barcode) %>% as.data.frame()
  sampleSigsExposures <- msa_sbs96_activity %>% filter(Tumor_Barcode == barcode) %>% as.data.frame()
  rownames(sampleSigsExposures) <- sampleSigsExposures$Tumor_Barcode
  sampleSigsExposures <- sampleSigsExposures[,-1]
  
  reassign <- NULL
  
  reassign <- assignSignatureProbabilityToMutations(sampleMutations = sampleMutations,
                                                    sampleSigsExposures = sampleSigsExposures,
                                                    signatures = signatures)
  
  tmpsig <- reassign %>%
    as_tibble() %>% 
    separate(sigsProb, into = c("Signature", "Values"), sep = ";") %>% 
    count(Signature) %>% 
    pull(Signature) %>% 
    str_split(pattern = ":") %>%
    unlist()
  
  reassign <- reassign %>%
    as_tibble() %>% 
    separate(sigsProb, into = c("Signature", "Values"), sep = ";") %>%
    separate_wider_delim(Values, delim = ":", names = tmpsig) %>%
    select(-Signature) %>%
    mutate(across(one_of(tmpsig), as.numeric)) %>% 
    rename(mutType = context)

  tgct_sbs96_probability <-  tgct_sbs96_probability %>% bind_rows(reassign)
  
}

tgct_sbs96_probability <- tgct_sbs96_probability %>%   replace_na(replace = list(SBS1 = 0, SBS3 = 0, SBS5 = 0, SBS18 = 0, SBS87 = 0))


save(tgct_sbs96_probability, file='tgct_signature_probability.RData')


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

ggsave(filename = 'Mutational_signature_landscape_SBS.pdf',width = 16,height = 7,device = cairo_pdf)


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

ggsave(filename = 'Mutational_signature_landscape_ID.pdf',width = 16,height = 8,device = cairo_pdf)



# Mutational signature landscape DBS ------------------------------------------
tdata <- bind_rows(
  tgct_activity_all_ratio %>% 
    select(Tumor_Barcode,starts_with('DBS')) %>% 
    pivot_longer(-Tumor_Barcode) %>% 
    mutate(group='Ratio'),
  
  tgct_activity_all %>% 
    select(Tumor_Barcode,starts_with('DBS')) %>% 
    pivot_longer(-Tumor_Barcode) %>% 
    group_by(Tumor_Barcode) %>% 
    summarise(value=sum(value)) %>% 
    mutate(group='Count')
) %>% 
  left_join(
    wgs_data %>% select(Tumor_Barcode,Subtype)
  ) %>% 
  mutate(name = factor(name,levels = sort(names(tgct_activity_all)[-1]))) %>% 
  mutate(value=if_else(is.na(value),0,value))


pdata <- tdata %>% left_join(samlevs)
max(pdata$value)
scale_y <- 180
pdata_ratio <- pdata %>% filter(group == 'Ratio') %>% mutate(value=value * scale_y)
pdata_count <- pdata %>% filter(group == 'Count')

# cosine similarity
pdata_cos <- tgct_dbs78_msa_decompsite %>% 
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
    breaks = seq(0,180,length=5),labels = seq(0,1,length=5), #percent_format()(seq(0,1,length=6))
    sec.axis = sec_axis(~ .,breaks = seq(0,180,length=5),name = "Number of mutations") # Secondary y-axis
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

ggsave(filename = 'Mutational_signature_landscape_DBS.pdf',width = 16,height = 7,device = cairo_pdf)




# Prevalence --------------------------------------------------------------
conflicts_prefer(dplyr::lag)
nmutation_input <- 50
sigdata <- tgct_activity_all %>% select(Samples=Tumor_Barcode,starts_with('SBS'))
data_sig <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = 'SBS_prevalence.pdf',colset = sigcol,rel_widths = c(1.2,4))
sigdata <- tgct_activity_all %>% select(Samples=Tumor_Barcode,starts_with('ID'))
data_sig <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = 'ID_prevalence.pdf',colset = sigcol,rel_widths = c(2,4))
sigdata <- tgct_activity_all %>% select(Samples=Tumor_Barcode,starts_with('DBS'))
data_sig <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = 'DBS_prevalence.pdf',colset = sigcol,rel_widths = c(1.2,4))






# Comparison --------------------------------------------------------------
pdata <- tgct_activity_all %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  left_join(wgs_data %>% select(Tumor_Barcode,Subtype)) %>% 
  mutate(name=factor(name,levels=colnames(tgct_activity_all)[-1])) %>% 
  filter(str_starts(name,'DBS'))

pdata %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(log2(value+1)~Subtype,data=.)))

my_comparisons <-  list(c('Seminoma','Non-Seminoma'))

stat.test <- pdata %>% 
  mutate(value=log2(value+1)) %>% 
  group_by(name) %>% 
  rstatix::wilcox_test(value~ Subtype, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Subtype") %>%
  mutate(myformatted.p = if_else(p<0.001,sprintf("P = %.2e",p),sprintf("P = %.2f",p))) %>% 
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

ggsave(filename = 'SBS_subtype_comparisons.pdf',width = 9,height =6,device = cairo_pdf)  
ggsave(filename = 'ID_subtype_comparisons.pdf',width = 9,height =6,device = cairo_pdf)  
ggsave(filename = 'DBS_subtype_comparisons.pdf',width = 4,height =6,device = cairo_pdf)  





# PairedWise assocaition ------------------------------------------------
load('tgct_mutational_signatures.RData')
tgct_activity_all
library(ggcorrplot)

tgct_activity_all %>% 
  select(-Tumor_Barcode) %>% 
  correlate(method = 'pearson') %>% 
  rplot()

corr <- round(cor(tgct_activity_all[,-1]), 1)
p.mat <- ggcorrplot::cor_pmat(tgct_activity_all[,-1])

ggcorrplot(corr,p.mat = p.mat,
           hc.order = TRUE,
           type = "upper",sig.level = 0.01,
           outline.color = "white")+
  theme_ipsum_rc()+
  theme(
    #panel.background = element_rect(fill = "grey75"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'signature_correlation.pdf',width = 6.2,height = 5,device = cairo_pdf)



# fisher_results <- fisher_results %>% mutate(odds_ratio = if_else(p_value>0.05,NA_real_,odds_ratio)) %>% mutate(p_value = if_else(p_value < 0.05,p_value,NA_real_))%>% mutate(odds_ratio = if_else(odds_ratio==0,1/64,odds_ratio)) %>%  mutate(odds_ratio = if_else(is.infinite(odds_ratio),64,odds_ratio))
# 
# ggplot(fisher_results, aes(x = gene1, y = gene2)) +
#   geom_asymmat(aes(fill_tl = -log10(p_value), fill_br = log2(odds_ratio)),col='white') +
#   scale_fill_tl_gradient2(low = "white", high = "tomato",na.value = 'gray90',breaks=pretty_breaks(n=5),limits=c(0,16)) +
#   scale_fill_br_gradient2(low = 'purple3',high='orange3',midpoint = 0,mid = 'white',na.value = 'gray90',breaks=pretty_breaks(n=5),limits=c(-6,6)) +
#   labs(
#     title = "Others",
#     fill_tl = "Fisher's exact test\n-log10(p-value)",
#     fill_br = "Fisher's exact test\nlog2(OR)"
#   ) +
#   theme_ipsum_rc()+
#   theme(
#     #panel.background = element_rect(fill = "grey75"),
#     panel.grid = element_blank(),
#     plot.title = element_text(hjust = 0.5),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
#   ) +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete(expand = c(0, 0))+
#   panel_border(color = 'black',size = 0.5)
# 
# ggsave(filename = 'tmp4.pdf',width = 6.2,height = 5,device = cairo_pdf)








# Association between SBS1 and age at diagnosis ---------------------------

pdata <- tgct_activity_all %>% 
  select(Tumor_Barcode,SBS1,SBS5) %>%
  mutate(SBS_age=SBS1+SBS5) %>% 
  left_join(wgs_data %>% select(Subject,Tumor_Barcode)) %>% 
  left_join(
    wgs_clinical %>% select(Subject,Subtype,age_at_diagnosis)
  )

pdata %>% group_by(Subtype) %>% do(tidy(cor.test(.$age_at_diagnosis, .$SBS_age,data=.)))

pdata %>% 
  ggplot(aes(SBS_age,age_at_diagnosis,fill=Subtype))+
  geom_point(pch=21,size=2.5)+
  stat_smooth(aes(fill = Subtype, color = Subtype), method = "lm") +
  stat_cor(size=5,col=ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5))+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  facet_wrap(~Subtype,scale='free')+
  labs(x='Clock-like mutations (SBS1+SBS5)',y='Age at diagnosis')+
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


ggsave(filename = 'SBS1_subtype_age.pdf',width = 9,height =4,device = cairo_pdf)  


# TMB and replication stress --------------------------------------------
load('tcga_repstress_zscore.RData',verbose = T)

pdata <- tgct_tmb %>% 
  mutate(TCGA_Barcode=str_sub(Tumor_Barcode,1,15)) %>% 
  left_join(tcga_repstress_score_cbioportal) %>% 
  filter(!is.na(RepStress_ZScore)) %>% 
  left_join(wgs_data %>% select(Tumor_Barcode,Subtype)) 

pdata %>% group_by(Subtype) %>% do(tidy(cor.test(.$RepStress_ZScore, .$TMB,data=.))) %>% arrange((p.value))

pdata %>% 
  ggplot(aes(RepStress_ZScore,TMB,fill=Subtype))+
  geom_point(pch=21,size=2.5)+
  stat_smooth(aes(fill = Subtype, color = Subtype), method = "lm") +
  stat_cor(size=5,col=ncicolpal[1]) +
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = subtypecol) +
  scale_color_manual(values = subtypecol) +
  scale_y_continuous(breaks = pretty_breaks(n = 5))+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  facet_wrap(~Subtype,scale='free')+
  #facet_grid(Subtype~name,scales = 'free')+
  labs(x='Tumor Mutational Burden (TMB)',y='Repstress Z-score')+
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


ggsave(filename = 'sig_subtyp_tmb_repstress.pdf',width = 9,height =4,device = cairo_pdf)  


# Signature and replication stress --------------------------------------------
load('tcga_repstress_zscore.RData',verbose = T)

pdata <- tgct_activity_all %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  filter(!str_starts(name,'DBS')) %>% 
  #filter(name %in% c('ID1','ID2','ID9')) %>% 
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


ggsave(filename = 'sig_subtype_repstress.pdf',width = 24,height =24,device = cairo_pdf)  

ggsave(filename = 'sigIDs_subtype_repstress.pdf',width = 10,height =7,device = cairo_pdf)  






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
if(TRUE){
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

# for cancer types, only for pancancer data
if(TRUE){
  tal <- att.setv(g = tal, from = "Subtype", to = "nodeColor",cols = rev(as.character(subtypecol)))
  update_color <- function(){addLegend.color(obj = rdp, tal, title = "Subtype", position = "bottomright",vertical=T,dxtitle=50)}
}


if(TRUE){
  tal <- att.setv(g = tal, from = "Subtype", to = "nodeColor",cols = rev(as.character(subtypecol)))
  update_color <- function(){addLegend.color(obj = rdp, tal, title = "Subtype", position = "bottomright",vertical=T,dxtitle=50)}
}

if(TRUE){
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








# Clonality of mutational signature ---------------------------------------
load('tgct_signature_probability.RData',verbose = T)
load('tgct_mtvcf.RData',verbose = T)

load('wgs_data.RData')

tdata <- tgct_mtvcf %>%
  select(Tumor_Barcode,ID=MutationID,CLS) %>% 
  left_join(
    tgct_sbs96_probability %>% select(Tumor_Barcode,ID,SBS1:SBS87)
  ) %>% 
  left_join(
    wgs_data %>% select(Tumor_Barcode,Subtype)
  ) %>% 
  filter(!is.na(SBS1)) %>% 
  pivot_longer(cols = -c(Subtype,Tumor_Barcode,ID,CLS)) %>% 
  group_by(Subtype,name,CLS) %>% 
  summarise(value=sum(value)) %>% 
  ungroup()
  
tdata %>% 
  mutate(name = factor(name,levels=c('SBS1','SBS3','SBS5','SBS18','SBS87'))) %>% 
  ggplot(aes(name,value,fill=CLS))+
  geom_bar(position = "fill", stat = "identity")+
  facet_wrap(~Subtype) +
  labs(x="", y = "Proportion of total mutations",fill='')+
  scale_fill_manual(values = clscolors)+
  scale_y_continuous(breaks = pretty_breaks(),expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  theme_ipsum_rc(axis_text_size = 10,axis_title_just = 'm',axis_title_size = 14)+
  theme_ipsum_rc(strip_text_size = 14,ticks = T)+
  theme(panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5,face = 'bold'),plot.margin = margin(4,4,0,4),axis.title.x = element_text(hjust = 0.5,size = 12),legend.position = 'bottom')+
  panel_border(color = 'black')+
  coord_flip()

ggsave(filename = 'CLS_subtype_signature.pdf',width = 6,height = 2.5,device = cairo_pdf)

  
  

tdata <- tgct_sbs96_probability %>%
  filter(Gene_Name %in% c('KIT','KRAS','NRAS')) %>% 
  select(Tumor_Barcode,Gene_Name,SBS1:SBS87) %>% 
  pivot_longer(-c(Tumor_Barcode,Gene_Name)) %>% 
  group_by(Gene_Name,name) %>% 
  summarise(value=sum(value)) %>% 
  ungroup()


tdata

