# Adapted from scripts/publishing/figure_standardization/figures2_4_timing_layout.py; publication side effects removed.
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
from figure1_balanced_summary import read,write
OUT=STAGE
QA=STAGE/'qa'


D=ROOT/'results/analysis/bcor_reclassification_2026-09-21/figures/figure4/derived-data'


GROUPS=['DN1','DN2','DN3'];LABEL=dict(zip(GROUPS,['C19','C19/20','TWR']));COL=dict(zip(GROUPS,['#E37E78','#8C231D','#E6A0C4']))


GROUPS=['DN1','DN2','DN3'];LABEL=dict(zip(GROUPS,['C19','C19/20','TWR']));COL=dict(zip(GROUPS,['#E37E78','#8C231D','#E6A0C4']))


GROUPS=['DN1','DN2','DN3'];LABEL=dict(zip(GROUPS,['C19','C19/20','TWR']));COL=dict(zip(GROUPS,['#E37E78','#8C231D','#E6A0C4']))


H=48.5;FW=121.;GW=49.;PITCH=4.4


H=48.5;FW=121.;GW=49.;PITCH=4.4


H=48.5;FW=121.;GW=49.;PITCH=4.4


H=48.5;FW=121.;GW=49.;PITCH=4.4


def curve():
 d=pd.read_csv(D/'gain-timing/gbm_dn3_panel_d_gain_timing_curve.tsv',sep='\t');stats=pd.read_csv(D/'gain-timing/gbm_dn3_panel_d_gain_timing_stats.tsv',sep='\t');test=pd.read_csv(D/'gain-timing/gbm_dn3_panel_d_gain_timing_global_permutation.tsv',sep='\t')
 fig=plt.figure(figsize=(GW/25.4,H/25.4));ax=fig.add_axes([9.86/GW,8.82/H,(GW-9.86-2.9)/GW,(H-15.12)/H]);guides=[]
 for g,ls in zip(GROUPS,['solid','solid','solid']):
  z=d[d.Evo_Group==g];x=z.time_display.to_numpy();y=z.cum_frac_display.to_numpy();assert np.all(np.diff(x)>=0) and np.all(np.diff(y)>=0)
  line=ax.step(x,y,where='post',lw=.9,color=COL[g],ls=ls,label=LABEL[g])[0];assert np.array_equal(line.get_xdata(),x) and np.array_equal(line.get_ydata(),y)
  k=int(np.flatnonzero(y>=.5)[0]);assert y[k-1]<.5<=y[k];t50=float(x[k])
  v=ax.vlines(t50,0,.5,ls=(0,(2,2)),lw=.65,color=COL[g]);assert np.allclose(v.get_segments()[0],[[t50,0],[t50,.5]])
  ax.scatter([t50],[.5],s=8,facecolor='white',edgecolor=COL[g],linewidth=.6,zorder=5)
  guides.append({'group':g,'label':LABEL[g],'t50':t50,'fraction_before':float(y[k-1]),'fraction_at_crossing':float(y[k]),'curve_row_1based':k+1,'previous_weighted_mean_time':float(stats.loc[stats.Evo_Group==g,'weighted_mean_time'].iloc[0]),'guide_y_bottom':0.,'guide_y_top':.5})
 ax.axhline(.5,lw=.55,ls=(0,(2,2)),color='#777777')
 ax.set(xlim=(0,1),ylim=(0,1.02),xticks=[0,.5,1],yticks=[0,.5,1],yticklabels=['0%','50%','100%'],xlabel='Molecular time',ylabel='Cumulative burden')
 ax.tick_params(length=2,pad=1);ax.xaxis.labelpad=1;ax.yaxis.labelpad=1
 title=fig.text(.5,1-.567/H,'Cumulative gain timing',fontsize=8,va='top',ha='center')
 legend=fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1-4/H),ncol=3,frameon=False,borderaxespad=0,borderpad=0,handlelength=1,handletextpad=.3,columnspacing=.8,fontsize=6)
 from figure1_tall_lower_panels import ptext
 ax.text(.025,.94,ptext(float(test.pvalue.iloc[0]),'Global P'),transform=ax.transAxes,fontsize=6,va='top')
 fig.text(.52,.63/H,'Dashed: 50% cumulative burden',fontsize=6,ha='center')
 fig.canvas.draw();renderer=fig.canvas.get_renderer();lb=legend.get_window_extent(renderer)
 assert not lb.overlaps(ax.get_window_extent(renderer)) and not lb.overlaps(title.get_window_extent(renderer))
 fig.savefig(OUT/'panels/M4_g.pdf');plt.close(fig)
 pd.DataFrame(guides).to_csv(OUT/'derived/M4_g_50percent_guides.tsv',sep='\t',index=False)
 write(QA/'curve_audit.json',{'rows':len(d),'pvalue':float(test.pvalue.iloc[0]),'curve_coordinates_preserved':True,'quantile_definition':'First saved displayed time whose cumulative burden is >= 0.5; no interpolation or refitting','guides':guides,'horizontal_guide_y':.5,'legend_outside_data':True,'native_panel_mm':[GW,H],'plot_height_mm':H-15.12,'plot_width_mm':GW-9.86-2.9})

