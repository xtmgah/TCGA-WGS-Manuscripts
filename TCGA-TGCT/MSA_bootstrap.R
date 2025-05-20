set_wd()
libztw()




exposure_refdata1 <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/result_spextractor_cosmic3.4/TGCT.SBS96.all/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')
#exposure_refdata <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/result_spextractor_cosmic3.4/TGCT.SBS288.all/SBS288/Suggested_Solution/COSMIC_SBS288_Decomposed_Solution/Activities/COSMIC_SBS288_Activities.txt')

exposure_refdata2 <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/result_spextractor_cosmic3.4/TGCT.ID83.all/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Activities/COSMIC_ID83_Activities.txt')

exposure_refdata3 <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/result_spextractor_cosmic3.4/TGCT.DBS78.all/DBS78/Suggested_Solution/COSMIC_DBS78_Decomposed_Solution/Activities/COSMIC_DBS78_Activities.txt')


exposure_refdata <- bind_rows(
  exposure_refdata1 %>% pivot_longer(-Samples,names_to = 'Signature_name',values_to = 'Exposure') %>% mutate(Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_SBS96'),
  exposure_refdata2 %>% pivot_longer(-Samples,names_to = 'Signature_name',values_to = 'Exposure') %>% mutate(Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh37_ID83'),
  exposure_refdata3 %>% pivot_longer(-Samples,names_to = 'Signature_name',values_to = 'Exposure') %>% mutate(Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_DBS78')
)

# Bootstrap mutational signatures -----------------------------------------
load('~/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Database/Updated/signature_refsets.RData')

load('~/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Codes/Sigvisualfunc.RData')
msa_indel <- read_delim('MSA_indel.txt',delim = '\t',col_names = T)

signature_refsets %>% filter(str_detect(Signature_set_name,'COSMIC_v3.4')) %>% count(Signature_set_name)
detected_sigs <- exposure_refdata %>% 
  #filter(Cancer_Type == 'Biliary-AdenoCA') %>% 
  filter(Signature_set_name %in% c('COSMIC_v3.4_Signatures_GRCh38_SBS96','COSMIC_v3.4_Signatures_GRCh37_ID83','COSMIC_v3.4_Signatures_GRCh38_DBS78')) %>% 
  filter(Exposure>20) %>% 
  count(Signature_set_name,Signature_name,sort=T)

dir.create(path = paste0('Bootstrap/signatures'),recursive = T)

tmpdata <- detected_sigs %>% 
  select(Signature_set_name,Signature_name) %>% 
  unique() %>% 
  left_join(signature_refsets)

tmpdata %>% filter(is.na(Contribution))

for(sigset in unique(tmpdata$Signature_set_name)){
  profiletype <- str_extract(sigset, "SBS|DBS|ID")
  
  if(profiletype == 'SBS'){
    tmpdata %>% 
      filter(Signature_set_name == sigset) %>%
      select(Signature_name,MutationType,Contribution) %>% 
      pivot_wider(names_from = Signature_name,values_from = Contribution) %>% 
      left_join(
        sigprofiles %>% select(MutationType,Type,SubType)
      ) %>% 
      select(-MutationType) %>% 
      relocate(Type,SubType) %>% 
      write_csv(file = paste0('Bootstrap/signatures/sigProfiler_',profiletype,'_signatures.csv'),col_names = T)
  }
  
  if(profiletype == 'DBS'){
    tmpdata %>% 
      filter(Signature_set_name == sigset) %>%
      select(Signature_name,MutationType,Contribution) %>% 
      pivot_wider(names_from = Signature_name,values_from = Contribution) %>% 
      rename(`Mutation Type`=MutationType) %>%
      write_csv(file = paste0('Bootstrap/signatures/sigProfiler_',profiletype,'_signatures.csv'),col_names = T)
  }
  
  if(profiletype == 'ID'){
    tmpdata %>% 
      filter(Signature_set_name == sigset) %>%
      select(Signature_name,MutationType,Contribution) %>% 
      pivot_wider(names_from = Signature_name,values_from = Contribution) %>%
      right_join(msa_indel) %>% 
      select(-MutationType) %>% 
      relocate(`Mutation Type`) %>% 
      write_csv(file = paste0('Bootstrap/signatures/sigProfiler_',profiletype,'_signatures.csv'),col_names = T)
  }
  
  
}

# matrix profile
seqmatrix_refdata1 <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/TGCT_spmatrix/output/SBS/TGCT.SBS96.all') %>% pivot_longer(-MutationType,names_to = 'Sample',values_to = 'Mutations') %>% mutate(Profile = 'SBS96')

seqmatrix_refdata2 <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/TGCT_spmatrix/output/ID/TGCT.ID83.all') %>% pivot_longer(-MutationType,names_to = 'Sample',values_to = 'Mutations') %>% mutate(Profile = 'ID83')

seqmatrix_refdata3 <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/SigProfiler/TGCT_spmatrix/output/DBS/TGCT.DBS78.all') %>% pivot_longer(-MutationType,names_to = 'Sample',values_to = 'Mutations') %>% mutate(Profile = 'DBS78')


seqmatrix_refdata <- bind_rows(seqmatrix_refdata1,seqmatrix_refdata2,seqmatrix_refdata3)

dir.create(path = paste0('Bootstrap/TGCT'),recursive = T)
seqmatrix_refdata %>% 
  filter(Profile == 'SBS96') %>% 
  select(Sample,MutationType,Mutations) %>% 
  pivot_wider(names_from = Sample,values_from = Mutations) %>% 
  left_join(
    sigprofiles %>% select(MutationType,Type,SubType)
  ) %>% 
  select(-MutationType) %>% 
  relocate(Type,SubType) %>% 
  write_csv(file = paste0('Bootstrap/TGCT/WGS_TGCT.96.csv'),col_names = T)

seqmatrix_refdata %>% 
  filter(Profile == 'ID83') %>% 
  select(Sample,MutationType,Mutations) %>% 
  pivot_wider(names_from = Sample,values_from = Mutations) %>% 
  right_join(msa_indel) %>% 
  select(-MutationType) %>% 
  relocate(`Mutation Type`) %>% 
  write_csv(file = paste0('Bootstrap/TGCT/WGS_TGCT.indels.csv'),col_names = T)

seqmatrix_refdata %>% 
  filter(Profile == 'DBS78') %>% 
  select(Sample,MutationType,Mutations) %>% 
  pivot_wider(names_from = Sample,values_from = Mutations) %>% 
  rename(`Mutation Type`=MutationType) %>% 
  write_csv(file = paste0('Bootstrap/TGCT/WGS_TGCT.dinucs.csv'),col_names = T)



