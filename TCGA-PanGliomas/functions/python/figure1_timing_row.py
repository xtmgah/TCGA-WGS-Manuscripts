# Adapted from scripts/publishing/figure_standardization/figure1_timing_row.py; publication side effects removed.
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
import figure1_tall_lower_panels as lower
from figure1_balanced_summary import read,write
OUT=STAGE
QA=STAGE/'qa'


F_WIDTH=121.;G_WIDTH=49.;G_X=127.


F_WIDTH=121.;G_WIDTH=49.;G_X=127.


F_WIDTH=121.;G_WIDTH=49.;G_X=127.


PITCH=4.4


def atlas():
 bins=lower.read('panel_f_pcawg_timing_bins.tsv');events=lower.read('panel_f_pcawg_timing_events.tsv');counts=lower.read('panel_shared_subtype_counts.tsv').set_index('Subtype_Final').n_samples.to_dict()
 fig=plt.figure(figsize=(F_WIDTH/25.4,56.5/25.4));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,F_WIDTH),ylim=(56.5,0));ax.axis('off')
 ax.text(.8,.6,'Copy-number gain timing',fontsize=8,va='top')
 cats=[str(i) for i in range(1,23)]+['X','all','WGD'];cmap=LinearSegmentedColormap.from_list('timing',['#117733','#7FBF7B','#F7F7F7','#C9A0CF','#7B2382']);cells=[];sectors=[]
 for j,cat in enumerate(cats):
  x=10+(j+.5)*PITCH;ax.text(x,8,cat,ha='center',fontsize=6,va='bottom')
  for group,y in zip(lower.GROUPS,[14.5,25.5,36.5]):
   n=events[(events.Subtype_Final==group)&(events.category.astype(str)==cat)].Tumor_Barcode.nunique();prev=n/counts[group];radius=(.09+.38*np.sqrt(prev))*PITCH if n else 0
   cells.append({'subtype':group,'category':cat,'n':n,'denominator':counts[group],'prevalence':prev,'radius_mm':radius,'block':'continuous','x_mm':x,'y_mm':y})
   ax.add_patch(Rectangle((x-PITCH/2,y-5.5),PITCH,11,fc='white',ec='#DDDDDD',lw=.35))
   z=bins[(bins.Subtype_Final==group)&(bins.category.astype(str)==cat)].sort_values('time_mid');angle=90
   for r in z.itertuples():
    da=360*r.weight/z.weight.sum();ax.add_patch(Wedge((x,y),radius,-angle,-angle+da,facecolor=cmap(r.time_mid),edgecolor='none'))
    sectors.append({'subtype':group,'category':cat,'time_mid':float(r.time_mid),'weight':float(r.weight),'fraction':float(r.weight/z.weight.sum())});angle-=da
 for group,y in zip(lower.GROUPS,[14.5,25.5,36.5]):
  ax.text(.8,y-.7,group,fontsize=6.5,va='center');ax.text(.8,y+1.7,f'n={counts[group]}',fontsize=6,va='center')
 ax.text(10,45,'Molecular time (early to late)',fontsize=6.5,va='bottom')
 for i,t in enumerate(np.linspace(0,1,150)):ax.add_patch(Rectangle((10+i*54/150,47),54/150+.005,1.9,fc=cmap(t),ec='none'))
 for t in [0,.25,.5,.75,1]:ax.text(10+54*t,52,f'{t:g}',fontsize=6,ha='center',va='bottom')
 ax.text(76,45,'Prevalence',fontsize=6.5,va='bottom')
 for x,prev in zip([79,94,109],[.25,.5,1]):
  ax.add_patch(Circle((x,49),(.09+.38*np.sqrt(prev))*PITCH,fc='none',ec='#505050',lw=.5));ax.text(x+3,49,f'{prev:.0%}',fontsize=6,va='center')
 fig.savefig(OUT/'panels/M1_f.pdf');plt.close(fig)
 pd.DataFrame(cells).to_csv(OUT/'derived/M1_f_size_encoding.tsv',sep='\t',index=False)
 pd.DataFrame(sectors).to_csv(OUT/'derived/M1_f_sector_encoding.tsv',sep='\t',index=False)
 write(QA/'atlas_audit.json',{'cells':len(cells),'timing_bin_rows':len(bins),'event_rows':len(events),'denominators':counts,'continuous_categories':cats,'radius_scale_mm':PITCH,'radii_unchanged':True,'native_panel_mm':[F_WIDTH,56.5],'models_or_tests_refit':False})


def curve():
 d=lower.read('panel_g_gain_timing_curve.tsv');stats=lower.read('panel_g_gain_timing_summary.tsv');test=lower.read('panel_g_gain_timing_global_permutation.tsv')
 H=56.5;fig=plt.figure(figsize=(G_WIDTH/25.4,H/25.4));ax=fig.add_axes([9.86/G_WIDTH,8.82/H,(G_WIDTH-9.86-2.9)/G_WIDTH,41.38/H]);guides=[]
 for group in lower.GROUPS:
  z=d[d.Subtype_Final==group];x=z.time.to_numpy();y=z.cum_frac_display.to_numpy();assert np.all(np.diff(x)>=0) and np.all(np.diff(y)>=0)
  line=ax.step(x,y,where='post',lw=.9,color=lower.COL[group],label=group)[0];assert np.array_equal(line.get_xdata(),x) and np.array_equal(line.get_ydata(),y)
  k=int(np.flatnonzero(y>=.5)[0]);assert y[k-1]<.5<=y[k];t50=float(x[k])
  guide=ax.vlines(t50,0,.5,ls=(0,(2,2)),lw=.65,color=lower.COL[group]);assert np.allclose(guide.get_segments()[0],[[t50,0],[t50,.5]])
  ax.scatter([t50],[.5],s=8,facecolor='white',edgecolor=lower.COL[group],linewidth=.6,zorder=5)
  guides.append({'subtype':group,'t50':t50,'fraction_before':float(y[k-1]),'fraction_at_crossing':float(y[k]),'curve_row_1based':k+1,'previous_weighted_mean_time':float(stats.loc[stats.Subtype_Final==group,'weighted_mean_time'].iloc[0]),'guide_y_bottom':0.,'guide_y_top':.5})
 horizontal=ax.axhline(.5,lw=.55,ls=(0,(2,2)),color='#777777');assert np.array_equal(horizontal.get_ydata(),[.5,.5])
 ax.set(xlim=(0,1),ylim=(0,1.02),xticks=[0,.5,1],yticks=[0,.5,1],yticklabels=['0%','50%','100%'],xlabel='Molecular time',ylabel='Cumulative burden')
 ax.tick_params(length=2,pad=1);ax.xaxis.labelpad=1;ax.yaxis.labelpad=1
 title=fig.text(.8/G_WIDTH,1-.567/H,'Cumulative gain timing',fontsize=8,va='top')
 legend=fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1-4/H),ncol=3,frameon=False,borderaxespad=0,borderpad=0,handlelength=1,handletextpad=.3,columnspacing=.8,fontsize=6)
 ax.text(.025,.94,lower.ptext(float(test.pvalue.iloc[0]),'Global P'),transform=ax.transAxes,fontsize=6,va='top')
 fig.text(.52,.63/H,'Dashed: 50% cumulative burden',fontsize=6,ha='center')
 fig.canvas.draw();renderer=fig.canvas.get_renderer();lb=legend.get_window_extent(renderer);ab=ax.get_window_extent(renderer);tb=title.get_window_extent(renderer)
 assert not lb.overlaps(ab),('Legend overlaps data',lb.bounds,ab.bounds)
 assert not lb.overlaps(tb),('Legend overlaps title',lb.bounds,tb.bounds)
 fig.savefig(OUT/'panels/M1_g.pdf');plt.close(fig)
 pd.DataFrame(guides).to_csv(OUT/'derived/M1_g_50percent_guides.tsv',sep='\t',index=False)
 write(QA/'curve_audit.json',{'rows':len(d),'pvalue':float(test.pvalue.iloc[0]),'curve_coordinates_preserved':True,'quantile_definition':'First saved time whose displayed cumulative burden is >= 0.5; no interpolation or refitting','guides':guides,'horizontal_guide_y':.5,'legend_outside_data':True,'native_panel_mm':[G_WIDTH,H],'plot_height_mm':41.38,'plot_width_mm':G_WIDTH-9.86-2.9,'models_or_tests_refit':False})

