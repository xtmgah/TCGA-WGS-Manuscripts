# Adapted from scripts/publishing/figure_standardization/figure1_tall_lower_panels.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache
setup_matplotlib()
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Wedge,Circle
from matplotlib.colors import LinearSegmentedColormap
OUT=STAGE


DATA=ROOT/'results/analysis/pooled_814_rerun_2026-09-28/figures/figure1/derived-data'


COL={'ASTRO':'#5785C1','OLIGO':'#74A089','GBM':'#BD3027'}


GROUPS=['ASTRO','OLIGO','GBM']


def read(name):return pd.read_csv(DATA/name,sep='\t')


def save(fig,name):fig.savefig(OUT/'panels'/f'M1_{name}.pdf');plt.close(fig)


def ptext(v,prefix='q'):return f'{prefix} < 0.001' if v<.001 else f'{prefix} = {v:.3f}' if v<.01 else f'{prefix} = {v:.2f}'


def chronology():
 d=read('panel_h_chronological_timing_values.tsv');tests=read('panel_h_chronological_pairwise_tests.tsv')
 H=61.5;fig=plt.figure(figsize=(85/25.4,H/25.4));metrics=['MRCA age','Age at diagnosis','Latency'];audit={}
 for k,metric in enumerate(metrics):
  ax=fig.add_axes([.09+k*.307,7.3/H,.24,45.805/H]);z=d[d.Metric==metric];arrays=[z[z.Subtype_Final==g].value.dropna().to_numpy() for g in GROUPS]
  bp=ax.boxplot(arrays,positions=[0,1,2],widths=.55,showfliers=False,patch_artist=True,medianprops={'color':'#202020','linewidth':.75},boxprops={'linewidth':.6},whiskerprops={'linewidth':.5},capprops={'linewidth':.5})
  for p in bp['boxes']:p.set_facecolor('white')
  for j,(group,v) in enumerate(zip(GROUPS,arrays)):
   rng=np.random.default_rng(170+j+k*3);scatter=ax.scatter(j+rng.uniform(-.23,.23,len(v)),v,s=5,facecolor=COL[group],edgecolor='#333333',linewidth=.17,alpha=.75,zorder=3)
   assert np.array_equal(scatter.get_offsets()[:,1],v)
   assert np.allclose(bp['medians'][j].get_ydata(),np.median(v))
   ax.scatter([j],[v.mean()],s=12,marker='D',facecolor='white',edgecolor='#333333',linewidth=.45,zorder=4)
  ymax=max(max(v) for v in arrays)*1.02;zt=tests[tests.metric==metric].sort_values('comparison_order')
  for t,r in enumerate(zt.itertuples()):
   x1=GROUPS.index(r.group_1);x2=GROUPS.index(r.group_2);yy=ymax*(1.06+.19*t)
   ax.plot([x1,x1,x2,x2],[yy-ymax*.03,yy,yy,yy-ymax*.03],color='#444444',lw=.45);ax.text((x1+x2)/2,yy+ymax*.01,ptext(r.padj_bh),fontsize=6,ha='center',va='bottom')
  ax.set_ylim(0,ymax*1.63);ax.set_xlim(-.6,2.6);ax.set_xticks([0,1,2],[f'{g}\nn={len(v)}' for g,v in zip(GROUPS,arrays)],fontsize=6)
  ax.tick_params(axis='both',length=2,pad=1,labelsize=6)
  if k==0:ax.set_ylabel('Years',fontsize=6.5,labelpad=1)
  ax.grid(axis='y',lw=.35,color='#E8E8E8',zorder=0)
  # Current h/i/j titles are individually left aligned in the shared source panel.
  fig.text((k*85/3+1.275)/85,1-.6/H,metric,fontsize=8,va='top')
  audit[metric]={'counts':{g:len(v) for g,v in zip(GROUPS,arrays)},'means':{g:float(v.mean()) for g,v in zip(GROUPS,arrays)},'medians':{g:float(np.median(v)) for g,v in zip(GROUPS,arrays)},'data_area_height_mm':45.805}
 fig.text(.5,.6205/H,'Points: specimens; diamonds: means; boxes: median / IQR',fontsize=6,ha='center')
 save(fig,'hij');return audit


def bic():
 d=read('panel_i_bic_model_selection.tsv');H=61.5;fig=plt.figure(figsize=(85/25.4,H/25.4));fig.text(.015,1-.584/H,'Within-subtype model selection',fontsize=8,va='top');ymax=np.ceil(d.delta_bic.max()/100)*100+100
 for k,group in enumerate(GROUPS):
  z=d[d.Subtype_Final==group];ax=fig.add_axes([.10+k*.304,6.935/H,.24,43.615/H]);ax.vlines(z.G,0,z.delta_bic,color='#AAAAAA',lw=.7)
  points=ax.scatter(z.G,z.delta_bic,c=COL[group],s=13,edgecolor='white',lw=.35,zorder=3);assert np.array_equal(points.get_offsets()[:,1],z.delta_bic.to_numpy())
  selected=z[z.selected.astype(str).str.lower()=='true'];ax.scatter(selected.G,selected.delta_bic,c=COL[group],s=17,marker='D',edgecolor='white',lw=.35,zorder=4)
  for r in selected.itertuples():ax.text(r.G,25,'best',color=COL[group],fontsize=6,ha='center')
  ax.set(xlim=(.5,6.5),ylim=(-10,ymax),xticks=range(1,7),yticks=[0,200,400,600]);ax.set_title(f'{group}\nn={int(z.n_samples.iloc[0])}',fontsize=7.5,pad=2);ax.tick_params(length=2,pad=1,labelsize=6);ax.grid(axis='y',color='#E8E8E8',lw=.35)
  if k==0:ax.set_ylabel('ΔBIC from best model',fontsize=6.5,labelpad=1)
  else:ax.set_yticklabels([])
 fig.text(.55,.9125/H,'Mixture components (G)',fontsize=6.5,ha='center');save(fig,'k')
 return {'rows':len(d),'common_y_limit':float(ymax),'selected':d[d.selected.astype(str).str.lower()=='true'][['Subtype_Final','G','n_samples']].to_dict('records'),'data_area_height_mm':43.615}

