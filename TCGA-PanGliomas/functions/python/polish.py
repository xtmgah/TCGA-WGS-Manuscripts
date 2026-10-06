"""Apply manuscript presentation rules to newly drawn vectors, never to reference artwork."""
from common import *
from legend_geometry import fragment,clear_regions,spans,CONTENT_CACHE
from legend_streams import remove_text
from title_streams import shift
import shutil,re

POSITIONS=json.loads((REPO/'functions/styles/title_positions.json').read_text())
def put(page,x,y,value,size=6.5,style='Regular',color=(.13,.13,.13)):
    name='POLISH_'+style
    page.insert_font(fontname=name,fontfile=str(FONT/f'RobotoCondensed-{style}.ttf'))
    page.insert_text((x,y),value,fontsize=size,fontname=name,color=color)
def moves(n):
    ops=[]
    def move(old,xy):
        box=fitz.Rect(old);new=fitz.Rect(xy[0],xy[1],xy[0]+box.width,xy[1]+box.height)
        ops.append((box*MM,new*MM))
    if n==2:
        for old,xy in [([131.35,149.4,143.1,152.4],[162.3,113]),([144.2,149.4,156,152.4],[162.3,117]),
          ([134.7,193.65,149.75,196.4],[158.7,163.5]),([150.65,193.65,163.7,196.4],[158.7,167.1]),
          ([135.45,240.75,146.2,243.5],[140.5,212.6]),([157.5,240.75,167.4,243.5],[140.5,216.1]),
          ([142.35,242.95,153.65,245.7],[140.5,219.6])]:move(old,xy)
    if n==3:
        for old,xy in [([23.45,114.7,40.1,117.45],[14,81.6]),([41.25,114.7,55.8,117.45],[14,85]),
          ([3.85,183.35,23,186.05],[4,177.65]),([40.85,183.35,59.15,186.05],[43,177.65]),
          ([79.85,183.35,97.2,186.05],[82,177.65]),([116.85,183.35,134.05,186.05],[121,177.65]),
          ([3.85,187.75,24,190.45],[4,181.35]),([47.85,187.75,76.4,190.45],[43,181.35]),
          ([107.85,187.75,131.15,190.45],[82,181.35]),([3.9,191.45,10.05,194.2],[121,181.35]),
          ([16.85,191.35,22.65,194.05],[129,181.35]),([30.85,191.35,36.65,194.05],[138,181.35]),
          ([44.85,191.35,50.65,194.05],[147,181.35]),([3.9,176.95,82.85,179.55],[4,186.35]),
          ([3.9,180.15,75.95,182.75],[4,189.55])]:move(old,xy)
    return ops

def polish(n,path):
    raw=path.with_name(path.stem+'_unpolished.pdf');shutil.copy2(path,raw)
    doc=fitz.open(raw);page=doc[0];ops=moves(n)
    parts=[(dest,fragment(raw,box)) for box,dest in ops]
    if ops:
        clear_regions(doc,raw,[box for box,dest in ops])
        for dest,piece in parts:
            assert spans(piece[0]),(n,'Empty legend',dest)
            page.show_pdf_page(dest,piece,0,keep_proportion=False);piece.close()
    intermediate=path.with_name(path.stem+'_legends.pdf');doc.save(intermediate,garbage=4,deflate=True);doc.close()
    doc=fitz.open(intermediate);page=doc[0];ss=spans(page);targets=[];additions=[];removals=[]
    # Titles are layout text. Their fixed origins/typography are kept separately
    # from all data values, axes and plotted marks.
    for title in POSITIONS[str(n)]['titles']:
        if n==5 and 151*MM<title['origin'][1]<159*MM and title['origin'][0]<101*MM and title['text'] not in ['copy number','expression','classes']:
            additions.append(title);continue
        found=[s for s in ss if s['text']==title['text'] and abs(s['size']-title['size'])<.03]
        if not found:
            resized=[s for s in ss if s['text']==title['text'] and abs(s['size']-title['size'])<1]
            if len(resized)==1:
                removals.extend(resized);additions.append(title);continue
            assert n==3 and (title['origin'][1]<24*MM or title['text']=='Timing of chromosome gains'),(n,title['text'],'Missing title')
            additions.append(title);continue
        s=min(found,key=lambda s:abs(s['origin'][1]-title['origin'][1])*3+abs(s['origin'][0]-title['origin'][0]))
        targets.append({'text':s['text'],'origin_pt':s['origin'],'dx_pt':title['origin'][0]-s['origin'][0],
                        'dy_pt':title['origin'][1]-s['origin'][1]})
    if n==5:
        removals.extend(s for s in ss if abs(s['size']-8)<.03 and 151*MM<s['origin'][1]<159*MM and s['origin'][0]<101*MM and s['text'] not in ['copy number','expression','classes'])
        for s in ss:
            if s['origin'][0]<60*MM and s['origin'][1]<70*MM and re.fullmatch(r'(C19(?:/20)?|TWR) \(n = \d+\)',s['text']):
                removals.append(s);additions.append(dict(s,origin=[s['origin'][0],s['origin'][1]-.25]))
    if n==1:
        for labels,center in [(('ASTRO','n=263'),109.7),(('OLIGO','n=147'),135.54),(('GBM','n=399'),161.38)]:
            for value in labels:
                s=next(s for s in ss if s['text']==value and abs(s['size']-7.5)<.02)
                targets.append({'text':value,'origin_pt':s['origin'],'dx_pt':center*MM-(s['bbox'][0]+s['bbox'][2])/2})
    if n in [1,4]:
        for value in ['Molecular time (early to late)','Prevalence' if n==1 else 'Timed burden (Gb)']:
            s=next(s for s in ss if s['text']==value)
            # The manuscript's local white key backing also clears the adjoining
            # atlas cell border. Keep that presentation geometry at the old site.
            pad=.035*MM
            page.draw_rect(fitz.Rect(s['bbox'])+(-pad,-pad,pad,pad),color=None,fill=(1,1,1),width=0)
            removals.append(s);additions.append(dict(s,origin=[s['origin'][0],s['origin'][1]+(1 if n==1 else .65)*MM]))
    if n in [2,4]:
        value='Number at risk' if n==2 else 'At risk';center=46.5 if n==2 else 56.5
        s=next(s for s in ss if s['text']==value)
        targets.append({'text':value,'origin_pt':s['origin'],'dx_pt':center*MM-(s['bbox'][0]+s['bbox'][2])/2})
    if removals:remove_text(doc,intermediate,[{'text':s['text'],'origin_pt':s['origin']} for s in removals],fitz)
    # Shift in a separate pass because text removal changes stream offsets.
    removed=path.with_name(path.stem+'_removed.pdf');doc.save(removed,garbage=4,deflate=True);doc.close()
    doc=fitz.open(removed);shift(doc,removed,targets,fitz);page=doc[0]
    overlay=fitz.open();o=overlay.new_page(width=page.rect.width,height=page.rect.height)
    for s in additions:
        style='Italic' if 'Italic' in s['font'] else 'Bold' if 'Bold' in s['font'] else 'Regular'
        color=tuple(((s['color']>>k)&255)/255 for k in [16,8,0])
        put(o,*s['origin'],s['text'],s['size'],style,color)
    if n==3:
        import pandas as pd
        values=pd.read_csv(ROOT/'results/figures/figure3/derived-data/figure3_sample_metrics.tsv',sep='\t')
        counts=values.Evo_Group.value_counts()
        value=f'C17p n = {counts["ASTRO_Group1"]}; CTR n = {counts["ASTRO_Group2"]}'
        font=fitz.Font(fontfile=str(FONT/'RobotoCondensed-Regular.ttf'))
        # Centered on the two CT category labels, as in the manuscript.
        labels=[next(s for s in ss if s['text']==v) for v in ['Any HC','Multi-chr']]
        center=sum((s['bbox'][0]+s['bbox'][2])/2 for s in labels)/2
        put(o,center-font.text_length(value,fontsize=6)/2,63.8*MM,value,6)
    if n==4:
        # Both keys are in a combined Cairo text run; preserve their original positions.
        selected=[s for s in spans(page) if s['text'] in ['P < 0.05','Not significant']]
        temp=path.with_name(path.stem+'_titles.pdf');doc.save(temp,garbage=4,deflate=True);doc.close();doc=fitz.open(temp);page=doc[0]
        remove_text(doc,temp,[{'text':s['text'],'origin_pt':s['origin']} for s in selected],fitz)
        for s in selected:put(o,*s['origin'],'P ≥ 0.05' if s['text']=='Not significant' else s['text'],s['size'],color=(34/255,)*3)
    if n==5:
        # Original outside-right legend, translated 15 mm with the enlarged top row.
        from matplotlib import pyplot as plt
        from matplotlib.patches import Rectangle
        fig=plt.figure(figsize=(11.5/25.4,9.75/25.4))
        for j,(value,color) in enumerate([('None','#D5D5D5'),('EGFRvIII','#6FA5B9'),('All others','#00A08A')]):
            y=1.75+j*3.3
            fig.add_artist(Rectangle((.2093/11.5,1-(y+.9)/9.75),1.8/11.5,1.8/9.75,transform=fig.transFigure,facecolor=color,edgecolor='none'))
            fig.text(2.6093/11.5,1-(y+.5732)/9.75,value,fontsize=6.5,va='baseline',color='#222222')
        legend=STAGE/'panels/M5_fusion_key.pdf';fig.savefig(legend,transparent=True);plt.close(fig)
        with fitz.open(legend) as key:o.show_pdf_page(fitz.Rect(100.1*MM,163.305*MM,111.6*MM,173.055*MM),key,0)
    if len(o.get_contents()):page.show_pdf_page(page.rect,overlay,0)
    overlay.close();doc.save(path,garbage=4,deflate=True);doc.close();CONTENT_CACHE.clear()
    (STAGE/f'qa/Figure_{n}_presentation.json').write_text(json.dumps({'legend_translations':len(ops),'title_translations':len(targets),
        'title_additions':len(additions),'reference_artwork_read':False,'scientific_values_changed':False},indent=2)+'\n')
