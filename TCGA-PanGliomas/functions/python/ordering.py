# Adapted from scripts/publishing/figure_standardization/refinement/ordering.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache


COLORS={'Mut':'#30C358','Gain':'#E60751','LOH':'#398BE3','HD':'#664499','WGD':'#FFA500','MRCA':'#555555'}


LABELS={'Mut':'Mutation','Gain':'Gain','LOH':'Loss','HD':'Homozygous deletion','WGD':'WGD','MRCA':'MRCA'}


SPEC=json.loads(Path(__file__).with_name('pl_style.json').read_text())


PL_META={
 'M1_a':('Pan-glioma ordering',814),
 'M1_b':('Ordering group 1 (GBM-enriched)',406),
 'M1_c':('Ordering group 2 (ASTRO-enriched)',251),
 'M1_d':('Ordering group 3 (OLIGO-enriched)',157),
 'M2_a':('Astrocytoma ordering: C17p',172),
 'M2_b':('Astrocytoma ordering: CTR',91),
 'M4_a':('GBM ordering: C19',160),
 'M4_b':('GBM ordering: C19/20',134),
 'M4_c':('GBM ordering: TWR',105)}


EXPECTED=dict(zip(PL_META,[51,28,36,23,36,43,33,28,34]))


RFONT=fitz.Font(fontfile=str(FONT/'RobotoCondensed-Regular.ttf'))


def rgb(value):return tuple(int(value[i:i+2],16)/255 for i in [1,3,5])


def concise(label):
 s=re.sub(r'\s*\[[^\]]*\]','',str(label))
 return re.sub(r'\s*\([^)]*\)','',re.sub(r'^.*?:\s*','',s)).strip()


def dimensions(panel):
 p=next(row for row in json.loads(Path(__file__).with_name('layouts.json').read_text())[panel[:2]]['placements'] if row[0]==panel[-1])
 return p[3],p[4]-3.5


def geometry(panel):
 width,height=dimensions(panel)
 source=STAGE/'derived'/f'{panel}_ordering.tsv'
 density=STAGE/'derived'/f'{panel}_violin.tsv'
 assert density.is_file(),f'Frozen violin coordinates required: {density}'
 data=pd.read_csv(source,sep='\t');vd=pd.read_csv(density,sep='\t')
 assert 'display_index' in data, 'Preserve the saved factor order; never resort by medians.'
 ids=data.sort_values('display_index',kind='stable').ID.drop_duplicates().tolist()
 assert len(ids)==EXPECTED[panel]
 lo,hi=float(data.plotpos.min()),float(data.plotpos.max());assert hi>lo
 top=SPEC['event_area_top_mm'];bottom=height-SPEC['event_area_bottom_offset_mm']
 data_top=(top+bottom-SPEC['data_height_mm'])/2
 data_bottom=data_top+SPEC['data_height_mm']
 assert data_top>=top and data_bottom<=bottom
 events=[]
 for i,name in enumerate(ids):
  ev=data[data.ID==name];z=vd[np.isclose(vd.x,i+1)].sort_values('y',kind='stable')
  assert len(z)==515 and len(ev)==950,(panel,name,len(z),len(ev))
  assert abs(z.y.min()-ev.plotpos.min())<1e-12 and abs(z.y.max()-ev.plotpos.max())<1e-12
  assert abs(z.violinwidth.max()-1)<1e-12
  ys=data_bottom-(z.y.to_numpy()-lo)/(hi-lo)*SPEC['data_height_mm']
  vtop,vbottom=float(ys.min()),float(ys.max())
  room_above=vtop-SPEC['density_label_clearance_mm']-(top+SPEC['label_guard_mm'])
  room_below=bottom-SPEC['label_guard_mm']-(vbottom+SPEC['density_label_clearance_mm'])
  events.append(dict(index=i+1,name=name,short=concise(name),ev=ev,z=z,ys=ys,
    vtop=vtop,vbottom=vbottom,room=max(room_above,room_below),side='above' if room_above>=room_below else 'below'))
 return dict(panel=panel,width=width,height=height,source=source,density=density,data=data,
   lo=lo,hi=hi,data_top=data_top,data_bottom=data_bottom,area_top=top,area_bottom=bottom,events=events)


@lru_cache(maxsize=1)
def all_geometry():return {panel:geometry(panel) for panel in PL_META}


@lru_cache(maxsize=1)
def alias_map():
 # An alias needed by any panel is used for that exact locus in every panel.
 frozen=Path(__file__).with_name('event_aliases.json')
 if frozen.is_file():return json.loads(frozen.read_text())
 allg=all_geometry();names=sorted({e['short'] for g in allg.values() for e in g['events']})
 prior={}
 baseline=WORK/'pl_legibility/baseline/staging/derived'
 for p in sorted(baseline.glob('*_event_labels.tsv')):
  for row in pd.read_csv(p,sep='\t').itertuples():
   full=concise(row.full_label)
   if str(row.short_label)!=full:prior[full]=str(row.short_label)
 result={name:name for name in names}
 for g in allg.values():
  for e in g['events']:
   name=e['short'];label=result[name]
   if RFONT.text_length(label,fontsize=SPEC['event_label_pt'])/MM<=e['room']:continue
   chrom=re.match(r'(\d+|X|Y)[pq]',name);prefix=chrom.group(0) if chrom else 'Event'
   peers=[s for s in names if s.startswith(prefix)]
   result[name]=prior.get(name,f'{prefix}[{peers.index(name)+1}]')
 for g in allg.values():
  for e in g['events']:
   assert RFONT.text_length(result[e['short']],fontsize=SPEC['event_label_pt'])/MM<=e['room'],(g['panel'],e['name'],'alias too long')
 assert len(set(result.values()))==len(result),'Ambiguous shared locus alias'
 return result


def text_box(x,y,label,size,rotation=0,align='left'):
 """Exact font-metric bounding box, in native millimeters."""
 length=RFONT.text_length(label,fontsize=size)/MM
 if rotation==90:
  width=(RFONT.ascender-RFONT.descender)*size/MM
  return [x-width/2,y-length,x+width/2,y]
 x0=x-length/2 if align=='center' else x-length if align=='right' else x
 return [x0,y-RFONT.ascender*size/MM,x0+length,y-RFONT.descender*size/MM]


def write_text(page,x,y,label,size,rotation=0,align='left'):
 box=text_box(x,y,label,size,rotation,align)
 px=x+(RFONT.ascender+RFONT.descender)*size/MM/2 if rotation==90 else box[0]
 page.insert_text((px*MM,y*MM),label,fontname='RC',fontsize=size,rotate=rotation,color=rgb('#222222'))
 return box


def render_panel(panel):
 if panel in ['M1_b','M1_c','M1_d']:
  from group_specific import render_one
  return render_one(panel)
 g=all_geometry()[panel];width,height=g['width'],g['height'];events=g['events']
 aliases=alias_map();title,n=PL_META[panel]
 doc=fitz.open();page=doc.new_page(width=width*MM,height=height*MM)
 page.insert_font(fontname='RC',fontfile=str(FONT/'RobotoCondensed-Regular.ttf'))
 title_baseline=.6+RFONT.ascender*SPEC['panel_title_pt']/MM
 write_text(page,.8,title_baseline,title,SPEC['panel_title_pt'])
 write_text(page,width-.8,.8+RFONT.ascender*SPEC['count_pt']/MM,f'n = {n} specimens',SPEC['count_pt'],align='right')
 left=1.8;pitch=(width-2*left)/len(events);barbase=height-SPEC['prevalence_bar_base_offset_mm']
 baseline=height-SPEC['prevalence_baseline_offset_mm'];separator=height-5.15
 for i in range(len(events)):
  if i%2==0:
   page.draw_rect(fitz.Rect((left+i*pitch)*MM,g['area_top']*MM,(left+(i+1)*pitch)*MM,barbase*MM),fill=rgb('#F7F7F7'),color=None)
 wgd=g['data'][g['data'].cna=='WGD'];wgd_median=float(np.median(wgd.plotpos))
 wgd_y=g['data_bottom']-(wgd_median-g['lo'])/(g['hi']-g['lo'])*SPEC['data_height_mm']
 page.draw_line((left*MM,wgd_y*MM),((width-left)*MM,wgd_y*MM),color=rgb('#9B9B9B'),width=.45,dashes='[3 3] 0')
 page.draw_line((left*MM,separator*MM),((width-left)*MM,separator*MM),color=rgb('#777777'),width=.45)
 mappings=[];records=[];number_boxes=[]
 for i,e in enumerate(events):
  x=left+(i+.5)*pitch;ev,z=e['ev'],e['z'];color=rgb(COLORS[ev.cna.iloc[0]])
  half=z.violinwidth.to_numpy()*SPEC['violin_width_mm']/2;ys=e['ys']
  points=np.c_[np.r_[x-half,x+half[::-1]],np.r_[ys,ys[::-1]]]
  shape=page.new_shape();shape.draw_polyline([(a*MM,b*MM) for a,b in points]);shape.finish(fill=color,color=rgb('#333333'),width=SPEC['violin_outline_pt'],closePath=True);shape.commit()
  med=float(np.median(ev.plotpos));medy=g['data_bottom']-(med-g['lo'])/(g['hi']-g['lo'])*SPEC['data_height_mm']
  page.draw_line(((x-SPEC['median_width_mm']/2)*MM,medy*MM),((x+SPEC['median_width_mm']/2)*MM,medy*MM),color=(0,0,0),width=SPEC['median_line_pt'])
  label=aliases[e['short']];length=RFONT.text_length(label,fontsize=SPEC['event_label_pt'])/MM
  anchor=e['vtop']-SPEC['density_label_clearance_mm'] if e['side']=='above' else e['vbottom']+SPEC['density_label_clearance_mm']+length
  box=text_box(x,anchor,label,SPEC['event_label_pt'],90)
  # A background matching its column keeps the dashed reference out of the letters.
  page.draw_rect(fitz.Rect(*[v*MM for v in box])+(-.2,-.2,.2,.2),fill=rgb('#F7F7F7' if i%2==0 else '#FFFFFF'),color=None)
  end=e['vtop'] if e['side']=='above' else e['vbottom']
  near=box[3] if e['side']=='above' else box[1]
  page.draw_line((x*MM,end*MM),(x*MM,near*MM),color=rgb('#B7B7B7'),width=.25)
  write_text(page,x,anchor,label,SPEC['event_label_pt'],90)
  freq=float(ev.freq.iloc[0]);barheight=SPEC['prevalence_bar_height_mm']*freq/100
  page.draw_rect(fitz.Rect((x-SPEC['prevalence_bar_width_mm']/2)*MM,(barbase-barheight)*MM,(x+SPEC['prevalence_bar_width_mm']/2)*MM,barbase*MM),fill=color,color=None)
  display=str(int(np.floor(freq+.5)))
  numbox=write_text(page,x,baseline,display,SPEC['prevalence_label_pt'],align='center');number_boxes.append(numbox)
  reconstructed=g['lo']+(g['data_bottom']-ys)/SPEC['data_height_mm']*(g['hi']-g['lo'])
  error=float(np.max(np.abs(reconstructed-z.y.to_numpy())))
  assert error<1e-12
  assert box[1]>=g['area_top']+SPEC['label_guard_mm']-1e-8 and box[3]<=g['area_bottom']-SPEC['label_guard_mm']+1e-8
  mappings.append(dict(panel=panel,display_index=i+1,short_label=label,full_label=e['name'],event_type=ev.cna.iloc[0],prevalence_percent=freq,median_plotpos=med,saved_realizations=len(ev)))
  records.append(dict(display_index=i+1,label=label,side=e['side'],label_bbox_mm=box,violin_top_mm=e['vtop'],violin_bottom_mm=e['vbottom'],violin_height_mm=e['vbottom']-e['vtop'],violin_max_width_mm=float(2*half.max()),median_plotpos=med,median_y_mm=medy,prevalence_percent=freq,display_prevalence=display,prevalence_baseline_mm=baseline,prevalence_bbox_mm=numbox,prevalence_bar_height_mm=barheight,density_vertices=len(z),inverse_affine_max_error=error))
 gaps=[number_boxes[i+1][0]-number_boxes[i][2] for i in range(len(number_boxes)-1)]
 assert min(gaps)>=1,(panel,'prevalence labels too close',min(gaps))
 output=STAGE/'panels'/f'{panel}.pdf';doc.set_metadata({'title':title,'creator':'Frozen PL geometry; native-size vector display'})
 doc.save(output,garbage=4,deflate=True);doc.close()
 native=fitz.open(output);spans=[s for b in native[0].get_text('dict')['blocks'] if 'lines'in b for line in b['lines'] for s in line['spans'] if s['text'].strip()]
 assert min(s['size'] for s in spans)>=5.99
 assert all(native[0].rect.contains(fitz.Rect(s['bbox'])) for s in spans),(panel,'native text overflow')
 assert all('RobotoCondensed' in s['font'] for s in spans)
 pd.DataFrame(mappings).to_csv(STAGE/'derived'/f'{panel}_event_labels.tsv',sep='\t',index=False)
 audit=dict(panel=panel,status='PASS',page_mm=[width,height],specimens=n,events=len(events),saved_realizations=len(g['data']),source_y_range=[g['lo'],g['hi']],data_top_mm=g['data_top'],data_bottom_mm=g['data_bottom'],data_height_mm=SPEC['data_height_mm'],maximum_violin_width_mm=SPEC['violin_width_mm'],median_violin_height_mm=float(np.median([r['violin_height_mm'] for r in records])),prevalence_baseline_mm=baseline,minimum_prevalence_gap_mm=min(gaps),prevalence_text_rotation=0,wgd_median=wgd_median,wgd_y_mm=wgd_y,minimum_font_pt=min(s['size'] for s in spans),native_text_overflow=False,input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [g['source'],g['density']]},output_sha256=sha(output),records=records)
 (STAGE/'qa'/f'{panel}_pl_geometry.json').write_text(json.dumps(audit,indent=2)+'\n')
 return mappings


def render_pl(input_tsv,id,title,subtitle,width_mm,height_mm):
 """Compatibility call; geometry, headers and labels are now centrally specified."""
 assert Path(input_tsv)==STAGE/'derived'/f'{id}_ordering.tsv'
 assert dimensions(id)==(width_mm,height_mm)
 return render_panel(id)


def render_all():
 return {panel:render_panel(panel) for panel in PL_META}


def draw_figure_key(page):
 page.insert_font(fontname='PLKEY',fontfile=str(FONT/'RobotoCondensed-Regular.ttf'))
 x=SPEC['key_x_mm'];baseline=SPEC['key_baseline_mm'];size=SPEC['legend_pt'];items=[]
 for key in SPEC['key_order']:
  label=LABELS[key]
  page.draw_rect(fitz.Rect(x*MM,(baseline-1.3)*MM,(x+1.25)*MM,(baseline-.05)*MM),fill=rgb(COLORS[key]),color=None)
  page.insert_text(((x+1.8)*MM,baseline*MM),label,fontname='PLKEY',fontsize=size,color=(.13,.13,.13))
  items.append({'event_type':key,'label':label,'x_mm':x,'baseline_mm':baseline,'font_pt':size,'color':COLORS[key]})
  x+=1.8+RFONT.text_length('LOH' if key=='LOH' else label,fontsize=size)/MM+3
 page.insert_text((SPEC['key_x_mm']*MM,SPEC['guide_baseline_mm']*MM),SPEC['guide'],fontname='PLKEY',fontsize=size,color=(.13,.13,.13))
 return {'items':items,'guide':SPEC['guide'],'guide_position_mm':[SPEC['key_x_mm'],SPEC['guide_baseline_mm']],'font_pt':size}

