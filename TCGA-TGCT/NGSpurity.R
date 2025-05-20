set_wd()
libztw()


# load NGSpurity result ---------------------------------------------------
load('wgs_metrics.RData', verbose = T)
load('ngspurity_data.RData', verbose = T)

ngspurity <- ngspurity %>%
  left_join(
    wgs_metrics %>%
      select(Tumor_Barcode = Sample_Barcode, MEAN_COVERAGE, MEDIAN_COVERAGE)
  ) %>%
  mutate(
    NRPCC = BB_Purity *
      MEAN_COVERAGE /
      (BB_Purity * BB_Ploidy + (1 - BB_Purity) * 2)
  )

save(
  ngspurity0,
  ngspurity,
  BBprofile0,
  BBprofile,
  clone_data,
  subclone_chrx,
  file = 'ngspurity_data.RData'
)
