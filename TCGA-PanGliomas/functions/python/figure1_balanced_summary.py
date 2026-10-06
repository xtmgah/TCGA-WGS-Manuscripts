# Adapted from scripts/publishing/figure_standardization/figure1_balanced_summary.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache
import group_specific,summary_renderer
from ordering import RFONT,rgb
ORIGINAL=STAGE
OUT=STAGE
QA=STAGE/'qa'


Y=51.833333333333336;EY=99.16666666666669


Y=51.833333333333336;EY=99.16666666666669


WIDTHS=[168*14/38,168*14/38,168*10/38]


XS=[4,4+WIDTHS[0]+2,4+sum(WIDTHS[:2])+4]


COL={'ASTRO':'#5785C1','OLIGO':'#74A089','GBM':'#BD3027','BCOR fusion':'#8D7AAE'}


def read(p):return json.loads(Path(p).read_text())


def write(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')


def rel(p):return str(Path(p).relative_to(ROOT))


def spans(p):return [s for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['text'].strip()]


def text(p,x,y,s,size=6.5,align='left',bold=False,color='#222222'):
 font=fitz.Font(fontfile=str(FONT/f'RobotoCondensed-{"Bold" if bold else "Regular"}.ttf'))
 width=font.text_length(s,fontsize=size)/MM
 if align=='center':x-=width/2
 elif align=='right':x-=width
 name='BALANCED_BOLD' if bold else 'BALANCED_REG';p.insert_font(fontname=name,fontfile=str(FONT/f'RobotoCondensed-{"Bold" if bold else "Regular"}.ttf'))
 p.insert_text((x*MM,y*MM),s,fontname=name,fontsize=size,color=rgb(color))


def selection():
 original,_,_,_=group_specific.selection()
 removed=sorted((r for r in original['M1_c'] if r.event_type!='WGD'),key=lambda r:(r.prevalence_percent,r.display_index))[:6]
 ids={r.display_index for r in removed}
 kept={p:[r for r in rr if p!='M1_c' or r.display_index not in ids] for p,rr in original.items()}
 assert [len(kept[p]) for p in ['M1_b','M1_c','M1_d']]==[14,14,10]
 assert max(r.prevalence_percent for r in removed)<min(r.prevalence_percent for r in kept['M1_c'] if r.event_type!='WGD')
 return kept,removed


def render_pl():
 kept,removed=selection();old=summary_renderer.OUT
 try:
  summary_renderer.OUT=OUT
  for letter,w in zip('bcd',WIDTHS):
   panel=f'M1_{letter}';audit=summary_renderer.render(panel,kept[panel],w,[1.]*(len(kept[panel])-1))
   assert audit['minimum_prevalence_gap_mm']>=1
   records=[]
   for n,(row,a) in enumerate(zip(kept[panel],audit['records']),1):
    records.append({'panel':panel,'display_index':n,'original_display_index':row.display_index,'short_label':a['label'],'full_label':row.full_label,'event_type':row.event_type,'prevalence_percent':row.prevalence_percent,'median_plotpos':row.median_plotpos,'saved_realizations':row.saved_realizations})
   pd.DataFrame(records).to_csv(OUT/f'derived/{panel}_event_labels.tsv',sep='\t',index=False)
 finally:summary_renderer.OUT=old
 pd.DataFrame([r._asdict() for r in removed]).to_csv(OUT/'derived/M1_c_omitted_events.tsv',sep='\t',index=False)


def composition():
 data=pd.read_csv(ORIGINAL/'derived/M1_e_subtype_composition_display.tsv',sep='\t');data.to_csv(OUT/'derived/M1_e_subtype_composition_display.tsv',sep='\t',index=False)
 d=fitz.open();p=d.new_page(width=172*MM,height=8.5*MM);text(p,.8,2.9,'Subtype composition (%)',8)
 for x,label in [(64,'ASTRO'),(78,'OLIGO'),(92,'GBM'),(106,'BCOR fusion')]:
  p.draw_circle((x*MM,1.6*MM),.65*MM,fill=rgb(COL[label]),color=None);text(p,x+1.5,2.4,label,6.5)
 alignment=[];segments=[]
 for i,(x,w) in enumerate(zip(XS,WIDTHS)):
  start=x-4+.8;usable=w-1.6;at=start;z=data[data.PL_Group==f'DN{i+1}']
  for name in ['GBM','ASTRO','OLIGO','BCOR fusion']:
   a=z[z.Subtype_Final==name]
   if len(a):
    row=a.iloc[0];width=usable*float(row.percent_within_group)/100
    p.draw_rect(fitz.Rect(at*MM,3.4*MM,(at+width)*MM,5.05*MM),fill=rgb(COL[name]),color=(1,1,1),width=.2)
    segments.append({'group':f'DN{i+1}','subtype':name,'n':int(row.n),'denominator':int(row.group_total),'percent':float(row.percent_within_group),'x_mm':at,'width_mm':width,'bar_width_mm':usable});at+=width
  assert abs(at-(start+usable))<1e-10
  dominant=z.loc[z.percent_within_group.idxmax()]
  text(p,x-4+w/2,7.75,f'Group {i+1}: {dominant.percent_within_group:.1f}% {dominant.Subtype_Final}',6.5,'center',color=COL[dominant.Subtype_Final])
  alignment.append({'group':f'DN{i+1}','PL_panel':f'M1_{"bcd"[i]}','panel_left_mm':x,'panel_width_mm':w,'PL_plot_left_mm':x+.8,'PL_plot_right_mm':x+w-.8,'composition_bar_left_mm':start+4,'composition_bar_right_mm':at+4})
 d.set_metadata({'title':'Subtype composition aligned to ordering-group summaries'});d.save(OUT/'panels/M1_e.pdf',garbage=4,deflate=True);d.close()
 write(QA/'composition_alignment.json',{'bars':alignment,'segments':segments,'group_denominators':[406,251,157],'each_bar_represents_percent_within_group':True})

