# Adapted from scripts/publishing/figure_standardization/figure3_feature_panel.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache


FONTS={k:fitz.Font(fontfile=str(FONT/f'RobotoCondensed-{k}.ttf')) for k in ['Regular','Bold']}


COLS={'C17p':'#8FAED6','CTR':'#315686'}


def rgb(s):return tuple(int(s[j:j+2],16)/255 for j in [1,3,5])


def txt(p,x,y,s,size=6.5,style='Regular',color='#222222',align='left'):
 font=FONTS[style];w=font.text_length(s,fontsize=size)/MM
 if align=='right':x-=w
 elif align=='center':x-=w/2
 p.insert_font(fontname='RC_'+style,fontfile=str(FONT/f'RobotoCondensed-{style}.ttf'))
 p.insert_text((x*MM,y*MM),s,fontname='RC_'+style,fontsize=size,color=rgb(color))
 return w


def save(d,path,title):
 d.set_metadata({'title':title,'creator':'Approved Figure 3h production renderer'})
 d.save(path,garbage=4,deflate=True)
 d.close()


def render(rows,output):
 d=fitz.open();p=d.new_page(width=106*MM,height=46.5*MM)
 txt(p,.6,3.5,'Genomic feature differences',8)
 txt(p,12.7,7.4,'Feature',6,color='#505050')
 txt(p,48,7.4,'Higher in C17p',6,color='#505050')
 p.draw_line((64*MM,6.7*MM),(70*MM,6.7*MM),color=rgb('#505050'),width=.4)
 p.draw_line((69*MM,6.1*MM),(70*MM,6.7*MM),color=rgb('#505050'),width=.4)
 p.draw_line((69*MM,7.3*MM),(70*MM,6.7*MM),color=rgb('#505050'),width=.4)
 txt(p,105,7.4,'BH-adjusted P',6,color='#505050',align='right')
 x0=48;wid=40;yfirst=10.1;pitch=2.37
 family=[(0,1,['CIN','expression']),(2,5,['Copy-number','burden']),(6,6,['HRD scar']),(7,11,['SV / CT','features'])]
 labels={'FGA':'Genome altered','Maximum clustered SVs':'Largest SV cluster','Oscillating CN segments':'Oscillating CN segments'}
 for a,b,ls in family:
  cy=yfirst+(a+b)*pitch/2
  for j,s in enumerate(ls):txt(p,.6,cy+.6+(j-(len(ls)-1)/2)*2.2,s,6,color='#606060')
  if a:
   yy=yfirst+(a-.5)*pitch
   p.draw_line((.6*MM,yy*MM),(105*MM,yy*MM),color=rgb('#E8E8E8'),width=.35)
 for tick in [0,.2,.4]:
  x=x0+wid*tick/.5
  p.draw_line((x*MM,8.8*MM),(x*MM,37.6*MM),color=rgb('#E6E6E6'),width=.35)
  txt(p,x,40.1,f'{tick:g}',6,align='center')
 for i,s in enumerate(rows):
  y=yfirst+i*pitch;v=float(s['standardized_mean_difference'])
  txt(p,12.7,y+.7,labels.get(s['metric'],s['metric']),6.5)
  p.draw_rect(fitz.Rect(x0*MM,(y-.60)*MM,(x0+wid*v/.5)*MM,(y+.60)*MM),fill=rgb(COLS['C17p']),color=None)
  pv=float(s['adjusted_pvalue']);ps='<0.001' if pv<.001 else f'{pv:.3f}'
  txt(p,104.5,y+.7,ps,6.5,align='right')
 p.draw_line((x0*MM,37.6*MM),((x0+wid)*MM,37.6*MM),color=rgb('#707070'),width=.5)
 p.draw_line((x0*MM,8.8*MM),(x0*MM,37.6*MM),color=rgb('#707070'),width=.5)
 txt(p,x0+wid/2,43,'Difference between groups',6,align='center')
 txt(p,.6,45.7,'BH correction across the 12 displayed features.',6,color='#505050')
 save(d,output,'Genomic feature differences: BH-adjusted comparisons')

