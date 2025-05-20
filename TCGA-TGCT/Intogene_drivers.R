set_wd()
libztw()
pdfhr2()


intogene <- read_delim('/Volumes/Zhang_Group/ZTW/TCGA-TGCT/Mutations/intogene/output/drivers.tsv',delim = '\t',col_names = T)


## visulaizaiton

intogene %>% 
  mutate(SYMBOL = fct_reorder(SYMBOL,SAMPLES)) %>% 
  mutate(Frequency=SAMPLES/252) %>% 
  ggplot(aes(SYMBOL,Frequency,size=SAMPLES,fill=-log10(QVALUE_COMBINATION)))+
  geom_point(pch=21)+
  scale_y_continuous(breaks = pretty_breaks(),labels = percent_format())+
  scale_fill_viridis_c()+
  scale_size_binned()+
  labs(x=NULL,y='Mutated Frequency',fill="-log10(qvalue)\n",size='nSample')+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'XY',ticks = T)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face =3),
        legend.position = 'bottom',
        legend.key.height = unit(0.3,units = 'cm'),
        legend.key.width = unit(0.7,units = 'cm'),
        plot.margin = margin(4,4,4,4))+
  coord_cartesian(clip = 'off')+
  panel_border(color = 'black',size = 0.4)

ggsave(file='intogene_drivers.pdf',width = 6,height = 5, device = cairo_pdf)

intogene %>% 
  mutate(SYMBOL = fct_reorder(SYMBOL,SAMPLES)) %>% 
  mutate(Frequency=SAMPLES/252) %>% 
  ggplot(aes(SYMBOL,Frequency,size=-log10(QVALUE_COMBINATION),fill=ROLE))+
  geom_point(pch=21)+
  scale_y_continuous(breaks = pretty_breaks(),labels = percent_format())+
  scale_fill_bmj()+
  scale_size_binned()+
  labs(x=NULL,y='Mutated Frequency',size="-log10(qvalue)\n",fill='Role')+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'XY',ticks = T)+
  guides(fill=guide_legend(override.aes = list(size=4)))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face =3),
        legend.position = 'bottom',
        legend.key.height = unit(0.3,units = 'cm'),
        legend.key.width = unit(0.7,units = 'cm'),
        plot.margin = margin(4,4,4,4),
        legend.box.margin = margin(r = 45))+
  coord_cartesian(clip = 'off')+
  panel_border(color = 'black',size = 0.4)

ggsave(file='intogene_drivers2.pdf',width = 6,height = 5, device = cairo_pdf)

save(intogene, file='intogene_drivers.RData')
