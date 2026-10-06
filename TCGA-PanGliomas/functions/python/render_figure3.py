# Adapted from scripts/publishing/figure_standardization/refinement/render_figure3.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
for n,f in [('RC','RobotoCondensed-Regular.ttf'),('RCB','RobotoCondensed-Bold.ttf')]:pdfmetrics.registerFont(TTFont(n,str(FONT/f)))


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def yes(value):return str(value).lower() in ('true','1')


def matrix():
    source=ROOT/'results/analysis/bcor_reclassification_2026-09-21/figures/figure3/source-data/panel-i'
    patients=sorted(read(source/'patient_column_key.tsv'),key=lambda r:int(r['column_index']))
    raw=read(source/'patient_gene_matrix.tsv')
    genes=['CDK4','PDGFRA','CDK6','MYCN','MDM2','CDKN2A','RB1','BRCA2']
    labels=['CDK4 / AVIL','PDGFRA','CDK6','MYCN','MDM2','CDKN2A / B','RB1','BRCA2']
    rows=[r for r in raw if r['gene'] in genes]
    assert len(patients)==88 and len(rows)==704
    assert [sum(p['route']==g for p in patients) for g in ['C17p','CTR']]==[69,19]
    assert sum(yes(r['missing_CN']) for r in rows)==2
    lookup={(r['Subject'],r['gene']):r for r in rows};assert len(lookup)==704
    W,H=172,70.5
    c=canvas.Canvas(str(STAGE/'panels/M3_i.pdf'),pagesize=(W*MM,H*MM),initialFontName='RC',initialFontSize=6)
    colors={'C17p':'#8FAED6','CTR':'#315686','no_CT':'#F0F2F4','CT_only':'#BFC8D0','AMP':'#D77753','CN0':'#287FA8','main':'#7755A2','secondary':'#42A79C','missing':'#F4DEA6','G2':'#E5E7EB','G3':'#9AA4AF','G4':'#3B4754'}
    def text(label,x,y,size=6.5,bold=False,color='#222222',align='left'):
        c.setFillColor(color);c.setFont('RCB' if bold else 'RC',size)
        {'left':c.drawString,'center':c.drawCentredString,'right':c.drawRightString}[align](x*MM,(H-y)*MM,str(label))
    def rect(x,y,w,h,color):c.setFillColor(color);c.rect(x*MM,(H-y-h)*MM,w*MM,h*MM,stroke=0,fill=1)
    def line(x,y,xx,yy,color='#D6DDE3',width=.3):c.setStrokeColor(color);c.setLineWidth(width);c.line(x*MM,(H-y)*MM,xx*MM,(H-yy)*MM)
    text('Copy-number alterations at selected genes overlapping chromothripsis',0,3,8)
    x0,gap=28,2.5;dx=(W-x0-gap-.5)/88
    xs=[x0+(i+.5)*dx+(gap if i>=69 else 0) for i in range(88)]
    for start,n,g in [(0,69,'C17p'),(69,19,'CTR')]:
        lo=xs[start]-dx/2;hi=xs[start+n-1]+dx/2
        rect(lo,5.2,hi-lo,3.5,colors[g]);text(f'{g} | {n} CT-positive patients' if g=='C17p' else 'CTR | 19 patients',(lo+hi)/2,7.7,6.5,True,'white' if g=='CTR' else '#202A35','center')
    for top,height,col,maxv,name in [(10.5,5.5,'CT_component_SVs',1500,'CT-component SVs'),(18,5.5,'CT_span_Mb',500,'CT span (Mb)')]:
        text(name,0,top+3,6.5)
        for value in [0,maxv/2,maxv]:
            yy=top+height*(1-value/maxv);line(x0,yy,W-.5,yy)
            text(int(value),x0-1,yy+.6,6,align='right')
        for j,p in enumerate(patients):
            height_value=height*float(p[col])/maxv
            assert height_value<=height+1e-6,(col,height_value)
            rect(xs[j]-dx*.40,top+height-height_value,dx*.80,height_value,colors[p['route']])
    text('Histological grade',0,27,6.5)
    for j,p in enumerate(patients):
        rect(xs[j]-dx*.42,25,dx*.84,2.1,colors[p['hist_grade']])
        if not yes(p['grade_explicit']):text('*',xs[j],27.1,6,True,'white','center')
    for i,(gene,label) in enumerate(zip(genes,labels)):
        top=29+i*2.6;text(label,0,top+1.8,6.5)
        for j,p in enumerate(patients):
            z=lookup[(p['Subject'],gene)]
            key='no_CT' if not yes(z['CT_overlap']) else 'missing' if yes(z['missing_CN']) else 'AMP' if yes(z['AMP']) else 'CN0' if yes(z['CN0']) else 'CT_only'
            left=xs[j]-dx*.43;right=xs[j]+dx*.43
            rect(left,top,dx*.86,2.2,colors[key])
            if yes(z['LOH_main']):rect(left,top,dx*.86,.60,colors['main'])
            if yes(z['LOH_secondary']):rect(left,top+1.6,dx*.86,.60,colors['secondary'])
            if yes(z['missing_CN']):line(left,top,right,top+2.2,'#202A35',.35);line(left,top+2.2,right,top,'#202A35',.35)
    for start,n in [(0,69),(69,19)]:
        for k in sorted(set([1]+list(range(10,n+1,10))+[n])):text(k,xs[start+k-1],51.1,6,align='center')
    text('88 CT-positive primary patients; decreasing CT-component SV count within each group.',0,54.3,6.5)
    text('Cohort CT: C17p 69/168; CTR 19/86.  * G4 derived from documented diagnosis.',0,57.5,6.5)
    legend=[(0,61,'AMP','Focal clonal AMP'),(37,61,'CN0','Whole-gene CN0'),(76,61,'main','Main-state LOH'),(113,61,'secondary','Secondary LOH'),(0,65.4,'missing','CT; CN incomplete'),(44,65.4,'CT_only','CT; none of these CN states'),(104,65.4,'no_CT','No CT overlap at gene')]
    for x,y,key,label in legend:
        rect(x,y-2,2,2,colors['CT_only'] if key in ['main','secondary'] else colors[key])
        if key in ['main','secondary']:rect(x,y-2 if key=='main' else y-.6,2,.6,colors[key])
        if key=='missing':line(x,y-2,x+2,y,'#202A35');line(x,y,x+2,y-2,'#202A35')
        text(label,x+3, y-.2,6.5)
    text('Grade:',0,69,6.5)
    for i,g in enumerate(['G2','G3','G4']):rect(13+i*14,67,2,2,colors[g]);text(g,16+i*14,68.8,6.5)
    c.save()
    return {'patients':88,'gene_rows':8,'cells':704,'group_sizes':{'C17p':69,'CTR':19},'incomplete_CN_cells':2,'patient_order':[p['Subject'] for p in patients],'gene_order':genes,'matrix_source_hashes':{str((source/n).relative_to(ROOT)):sha(source/n) for n in ['patient_column_key.tsv','patient_gene_matrix.tsv']},'minimum_text_pt':6}

