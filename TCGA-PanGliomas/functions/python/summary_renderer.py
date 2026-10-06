# Adapted from scripts/publishing/figure_standardization/refinement/summary_renderer.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache
from ordering import RFONT,COLORS,PL_META,all_geometry,alias_map,rgb,text_box,concise
from fontTools.ttLib import TTFont
OUT=STAGE


TT=TTFont(FONT/'RobotoCondensed-Regular.ttf');CMAP=TT.getBestCmap();UPM=TT['head'].unitsPerEm


TT=TTFont(FONT/'RobotoCondensed-Regular.ttf');CMAP=TT.getBestCmap();UPM=TT['head'].unitsPerEm


TT=TTFont(FONT/'RobotoCondensed-Regular.ttf');CMAP=TT.getBestCmap();UPM=TT['head'].unitsPerEm


BODY_HEIGHT=42.5


AREA_TOP=7.2


AREA_BOTTOM=BODY_HEIGHT-5.4


DATA_TOP=(AREA_TOP+AREA_BOTTOM-20)/2


DATA_BOTTOM=DATA_TOP+20


def write(path,data):path.write_text(json.dumps(data,indent=2)+'\n')


def key(row):return row.event_type,concise(row.full_label)


def put(page,x,y,label,size=6,align='left',font='PREVIEWREG'):
    length=RFONT.text_length(label,fontsize=size)/MM
    if align=='center':x-=length/2
    elif align=='right':x-=length
    page.insert_text((x*MM,y*MM),label,fontname=font,fontsize=size,color=(.13,.13,.13))


def rotated_label(page,x,y,label,background):
    # TrueType ink extents avoid painting a full line-height mask over neighbors.
    glyphs=[TT['glyf'][CMAP[ord(c)]] for c in label]
    low=min(g.yMin for g in glyphs)/UPM
    high=max(g.yMax for g in glyphs)/UPM
    anchor=x+(low+high)*6/MM/2
    half=(high-low)*6/MM/2
    length=RFONT.text_length(label,fontsize=6)/MM
    box=[x-half,y-length,x+half,y]
    page.draw_rect(fitz.Rect((box[0]-.025)*MM,(box[1]-.025)*MM,(box[2]+.025)*MM,(box[3]+.025)*MM),fill=rgb(background),color=None)
    page.insert_text((anchor*MM,y*MM),label,fontname='PREVIEWREG',fontsize=6,rotate=90,color=(.13,.13,.13))
    return box


def render(panel,rows,width,steps):
    g=all_geometry()[panel];aliases=alias_map()
    doc=fitz.open();page=doc.new_page(width=width*MM,height=BODY_HEIGHT*MM)
    page.insert_font(fontname='PREVIEWREG',fontfile=str(FONT/'RobotoCondensed-Regular.ttf'))
    title,n=PL_META[panel]
    put(page,.8,.6+RFONT.ascender*8/MM,title,8)
    put(page,width-.8,3.9+RFONT.ascender*6.5/MM,f'n = {n} specimens',6.5,align='right')
    centers=np.r_[1.85,1.85+np.cumsum(steps)]
    spare=width-(3.7+sum(steps));centers+=np.linspace(0,spare,len(rows))
    boundaries=np.r_[.8,(centers[:-1]+centers[1:])/2,width-.8]
    base=BODY_HEIGHT-2.9;baseline=BODY_HEIGHT-.6
    for i in range(len(rows)):
        if i%2==0:page.draw_rect(fitz.Rect(boundaries[i]*MM,AREA_TOP*MM,boundaries[i+1]*MM,base*MM),fill=rgb('#F7F7F7'),color=None)
    ymap=lambda y:DATA_BOTTOM-(np.asarray(y)-g['lo'])/(g['hi']-g['lo'])*20
    wg=float(np.median(g['data'].loc[g['data'].cna=='WGD','plotpos']))
    page.draw_line((.8*MM,float(ymap(wg))*MM),((width-.8)*MM,float(ymap(wg))*MM),color=rgb('#9B9B9B'),width=.45,dashes='[3 3] 0')
    page.draw_line((.8*MM,(BODY_HEIGHT-5.15)*MM),((width-.8)*MM,(BODY_HEIGHT-5.15)*MM),color=rgb('#777777'),width=.45)
    records=[]
    for i,(row,x) in enumerate(zip(rows,centers)):
        e=next(e for e in g['events'] if e['index']==row.display_index)
        z=e['z'];ys=ymap(z.y.to_numpy());half=z.violinwidth.to_numpy()*1.05
        points=np.c_[np.r_[x-half,x+half[::-1]],np.r_[ys,ys[::-1]]]
        color=rgb(COLORS[row.event_type]);shape=page.new_shape()
        shape.draw_polyline([(a*MM,b*MM) for a,b in points]);shape.finish(fill=color,color=rgb('#333333'),width=.3,closePath=True);shape.commit()
        med=float(row.median_plotpos)
        page.draw_line(((x-.45)*MM,float(ymap(med))*MM),((x+.45)*MM,float(ymap(med))*MM),width=.4,color=(0,0,0))
        top,bottom=float(ys.min()),float(ys.max())
        above=top-.65-(AREA_TOP+.2);below=AREA_BOTTOM-.2-(bottom+.65)
        side='above' if above>=below else 'below'
        label=aliases[e['short']];length=RFONT.text_length(label,fontsize=6)/MM
        assert length<=max(above,below),(panel,label,length,above,below)
        anchor=top-.65 if side=='above' else bottom+.65+length
        edge=top if side=='above' else bottom;near=anchor if side=='above' else anchor-length
        page.draw_line((x*MM,edge*MM),(x*MM,near*MM),color=rgb('#B7B7B7'),width=.25)
        ink=rotated_label(page,x,anchor,label,'#F7F7F7' if i%2==0 else '#FFFFFF')
        freq=float(row.prevalence_percent)
        page.draw_rect(fitz.Rect((x-1.05)*MM,(base-2*freq/100)*MM,(x+1.05)*MM,base*MM),fill=color,color=None)
        number=str(int(np.floor(freq+.5)));put(page,x,baseline,number,6,align='center')
        numbox=text_box(x,baseline,number,6,align='center')
        records.append({'original_display_index':row.display_index,'preview_display_index':i+1,'label':label,'full_label':row.full_label,'key':key(row),'x_mm':float(x),'side':side,'label_ink_bbox_mm':ink,'prevalence_bbox_mm':numbox,'prevalence_percent':freq,'display_prevalence':number,'median_plotpos':med,'median_y_mm':float(ymap(med)),'violin_height_mm':bottom-top,'maximum_width_mm':2.1,'density_vertices':len(z)})
    for record in records:
        saved=next(e['z'].y.to_numpy() for e in g['events'] if e['index']==record['original_display_index'])
        recovered=g['lo']+(DATA_BOTTOM-ymap(saved))/20*(g['hi']-g['lo'])
        record['inverse_affine_max_error']=float(np.max(np.abs(saved-recovered)))
    gaps=[records[i+1]['prevalence_bbox_mm'][0]-records[i]['prevalence_bbox_mm'][2] for i in range(len(records)-1)]
    ink_gaps=[records[i+1]['label_ink_bbox_mm'][0]-records[i]['label_ink_bbox_mm'][2] for i in range(len(records)-1)]
    assert min(gaps)>=.199 and min(ink_gaps)>0
    # Subtle dividers keep densely packed whole numbers legible as separate cells.
    for left,right in zip(records,records[1:]):
        boundary=(left['prevalence_bbox_mm'][2]+right['prevalence_bbox_mm'][0])/2
        page.draw_line((boundary*MM,(baseline-1.75)*MM),(boundary*MM,(baseline+.1)*MM),color=rgb('#BDBDBD'),width=.2)
    path=OUT/'panels'/f'{panel}.pdf';doc.save(path,garbage=4,deflate=True);doc.close()
    doc=fitz.open(path);spans=[s for b in doc[0].get_text('dict')['blocks'] if 'lines'in b for line in b['lines'] for s in line['spans']]
    assert all(s['size']>=5.99 and 'RobotoCondensed' in s['font'] and doc[0].rect.contains(fitz.Rect(s['bbox'])) for s in spans)
    audit={'status':'PASS','panel':panel,'events':len(rows),'width_mm':width,'body_height_mm':BODY_HEIGHT,'data_span_mm':20,'maximum_violin_width_mm':2.1,'source_y_range':[g['lo'],g['hi']],'wgd_plotpos':wg,'wgd_y_mm':float(ymap(wg)),'prevalence_baseline_mm':baseline,'minimum_prevalence_gap_mm':min(gaps),'minimum_label_ink_x_gap_mm':min(ink_gaps),'records':records,'sha256':sha(path)}
    write(OUT/'qa'/f'{panel}.json',audit)
    doc[0].get_pixmap(dpi=300).save(OUT/'qa'/f'{panel}_300dpi.png')
    return audit

