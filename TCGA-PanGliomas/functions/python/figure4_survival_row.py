# Adapted from scripts/publishing/figure_standardization/figure4_survival_row.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache
import figures2_4_timing_layout as previous
from figure1_balanced_summary import read,write
OUT=STAGE
QA=STAGE/'qa'
D=previous.D


def timing_panels():
 # The curve renderer already uses physical margins and adapts its data span
 # to native height. All saved coordinates and 50% references are retained.
 previous.OUT=OUT;previous.QA=QA;previous.H=35.5
 previous.curve()
 # Re-render the atlas at native size; never scale a completed panel's text.
 import matplotlib.pyplot as plt
 from matplotlib.patches import Rectangle,Wedge,Circle
 bins=pd.read_csv(D/'pcawg-timing/gbm_dn3_pcawg_timing_pie_bins.tsv',sep='\t');values=pd.read_csv(D/'pcawg-timing/gbm_dn3_pcawg_timing_pie_summary.tsv',sep='\t')
 colors=read(STAGE/'qa/M4_time_colors.json');total={(r.Evo_Group,str(r.category)):r.total_weight_mb*1e6 for r in values.itertuples()};maximum=max(total.values())
 fig=plt.figure(figsize=(121/25.4,35.5/25.4));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,121),ylim=(35.5,0));ax.axis('off')
 ax.text(60.5,.6,'Copy-number gain timing',fontsize=8,ha='center',va='top')
 cats=[str(i) for i in range(1,23)]+['X','all','WGD'];cells=[];sectors=[]
 radius=lambda w:(.13+.32*np.sqrt(w/maximum))*4.4 if w>0 else 0
 for j,cat in enumerate(cats):
  x=10+(j+.5)*4.4;ax.text(x,6.6,cat,ha='center',va='bottom',fontsize=6)
  for g,y in zip(previous.GROUPS,[11,17.5,24]):
   weight=total[g,cat];r=radius(weight);cells.append({'group':g,'category':cat,'weight_bp':weight,'radius_mm':r,'x_mm':x,'y_mm':y})
   ax.add_patch(Rectangle((x-2.2,y-3.25),4.4,6.5,fc='white',ec='#DDDDDD',lw=.35))
   z=bins[(bins.Evo_Group==g)&(bins.category.astype(str)==cat)].sort_values('time_mid');angle=90
   for row in z.itertuples():
    fraction=row.weight/weight;da=360*fraction;ax.add_patch(Wedge((x,y),r,-angle,-angle+da,facecolor=colors[int(row.time_bin)],edgecolor='none'));angle-=da
    sectors.append({'group':g,'category':cat,'time_mid':row.time_mid,'time_bin':int(row.time_bin),'weight':row.weight,'fraction':fraction,'color':colors[int(row.time_bin)]})
   assert abs(angle+270)<.01 or weight==0
 for g,y,n in zip(previous.GROUPS,[11,17.5,24],[160,134,105]):
  ax.text(.8,y-.7,previous.LABEL[g],fontsize=6.5,va='center',color='#000000');ax.text(.8,y+1.7,f'n={n}',fontsize=6,va='center',color='#000000')
 ax.text(10,29.5,'Molecular time (early to late)',fontsize=6.5,va='bottom')
 for i,c in enumerate(colors):ax.add_patch(Rectangle((10+i*54/100,31),54/100+.005,1.4,fc=c,ec='none'))
 for t in [0,.25,.5,.75,1]:ax.text(10+54*t,34.9,f'{t:g}',fontsize=6,ha='center',va='bottom')
 ax.text(76,29.5,'Timed burden (Gb)',fontsize=6.5,va='bottom')
 for x,w in zip([79,94,109],[1e9,1e10,1e11]):
  ax.add_patch(Circle((x,32.3),radius(w),fc='none',ec='#505050',lw=.5));ax.text(x+3,32.3,f'{w/1e9:g}',fontsize=6,va='center')
 fig.savefig(OUT/'panels/M4_f.pdf');plt.close(fig)
 pd.DataFrame(cells).to_csv(OUT/'derived/M4_f_relayout_size_encoding.tsv',sep='\t',index=False)
 pd.DataFrame(sectors).to_csv(OUT/'derived/M4_f_relayout_sector_encoding.tsv',sep='\t',index=False)
