"""Build all scientific panels anew, then compose the current five page layouts."""
from common import *
import sys, shutil
import numpy as np
import pandas as pd

LAYOUT=json.loads((REPO/'functions/styles/figures.json').read_text())
def compose_current(n,groups):
    spec=LAYOUT[f'M{n}'];doc=fitz.open();page=doc.new_page(width=180*MM,height=spec['height_mm']*MM)
    for source,letters in groups:
        boxes=[spec['panels'][letter] for letter in letters]
        x=min(b[0] for b in boxes);y=boxes[0][1];right=max(b[0]+b[2] for b in boxes);bottom=max(b[1]+b[3] for b in boxes)
        dest=fitz.Rect(x*MM,(y+3.5)*MM,right*MM,bottom*MM)
        with fitz.open(STAGE/f'panels/{source}.pdf') as src:
            from cairo_spaces import normalize
            normalize(src,STAGE/f'panels/{source}.pdf')
            # Pad Cairo's <1-point MediaBox rounding; retain all marks at native size.
            if 0<=dest.width-src[0].rect.width<1.01 and 0<=dest.height-src[0].rect.height<1.01:
                if n==3:
                    src[0].set_mediabox(fitz.Rect(0,0,dest.width,dest.height))
                else:
                    actual=src[0].rect
                    padded=fitz.open();pp=padded.new_page(width=dest.width,height=dest.height)
                    pp.show_pdf_page(actual,src,0)
                    page.show_pdf_page(dest,padded,0);padded.close()
                    continue
            assert abs(src[0].rect.width-dest.width)<.05 and abs(src[0].rect.height-dest.height)<.05,(source,src[0].rect,dest)
            page.show_pdf_page(dest,src,0)
    # Finish grafting panel resources before inserting overlay fonts. Otherwise
    # PyMuPDF's font cache can encode later labels against a grafted subset font.
    text(page,4,8,spec['title'],12,True)
    if n in [1,2,4]:
        from ordering import draw_figure_key
        draw_figure_key(page)
    for letter in spec['letters']:
        x,y,w,h=spec['panels'][letter]
        color=(34/255,)*3 if (n==1 and letter in 'bcde') or (n==5 and letter in 'fghi') else (.13,)*3
        text(page,x,y+2.8,letter,10,True,color)
    if n==3:
        for x,label,color in [(140,'C17p','#8FAED6'),(158,'CTR','#315686')]:
            rgb=tuple(int(color[k:k+2],16)/255 for k in (1,3,5))
            page.draw_rect(fitz.Rect(x*MM,10*MM,(x+1.5)*MM,11.5*MM),fill=rgb,color=None)
            text(page,x+2.3,11.5,label,6.5)
    if n==4:
        text(page,4,246.4,'Points: specimens | Box: median and IQR | Whiskers: 1.5 × IQR | White diamond: mean | Brackets: BH-adjusted q',6,color=(.2,.2,.2))
    target=STAGE/f'main_figures/Figure_{n}.pdf'
    doc.set_metadata({'title':spec['title'],'creator':'TCGA-PanGliomas standalone data-driven reproduction'})
    doc.save(target,garbage=4,deflate=True);doc.close()
    from polish import polish
    polish(n,target)
    (STAGE/f'qa/Figure_{n}_render.json').write_text(json.dumps({'figure':n,'status':'RENDERED',
        'panels':spec['letters'],'page_mm':[180,spec['height_mm']],
        'all_panels_drawn_from_inputs':True,'reference_artwork_used':False,
        'visual_equivalence':'pending comparison'},indent=2)+'\n')

def figure1():
    import ordering,figure1_balanced_summary as b,figure1_timing_row as timing,figure1_tall_lower_panels as lower
    ordering.render_panel('M1_a');b.render_pl();b.composition();timing.atlas();timing.curve();lower.chronology();lower.bic()
    compose_current(1,[(f'M1_{a}',a) for a in 'abcdefg']+[('M1_hij','hij'),('M1_k','k')])
def figure2():
    import ordering,render_figure2 as r
    ordering.render_panel('M2_a');ordering.render_panel('M2_b');r.volcano()
    compose_current(2,[(f'M2_{a}',a) for a in 'abcd']+[('M2_e','ef'),('M2_f','g'),('M2_g','hi'),('M2_h','j')])
def figure3():
    import render_figure3 as r,figure3_feature_panel as feature
    r.matrix();shutil.copy2(STAGE/'panels/M3_i.pdf',STAGE/'panels/M3_h.pdf')
    order=['CIN25','CIN70','FGA','Arm CNA burden','Aneuploidy','Ploidy','HRD-LOH','Inversion burden','Translocation burden','Oscillating CN segments','Maximum clustered SVs','Break load']
    rows=pd.read_csv(STAGE/'derived/figure3_h_BH_adjusted.tsv',sep='\t').to_dict('records')
    feature.render(sorted(rows,key=lambda r:order.index(r['metric'])),STAGE/'panels/M3_g.pdf')
    compose_current(3,[(f'M3_{a}',a) for a in 'abcdefgh'])
def figure4():
    import ordering,figure4_survival_row as r
    for letter in 'abc':ordering.render_panel('M4_'+letter)
    r.timing_panels()
    compose_current(4,[(f'M4_{a}',a) for a in 'abcdefghijk'])
def figure5():
    import render_figure5 as r,figure5_ac2 as ac,figure5_vertical_panels as vertical
    import figure5_space_palette_panels as palette,figure5_row3_panels as row3
    old=r.D;current=ac.D;r.D=current;vertical.D=current;palette.D=current
    ac.panel_a();ac.panel_b();ac.panel_c();vertical.panel_d()
    palette.panel_e()
    # Panel e/f/g is one newly rendered compound image with three current data panels.
    shutil.copy2(STAGE/'panels/M5_e.pdf',STAGE/'panels/M5_efg.pdf')
    row3.D=old;row3.panel_i();shutil.copy2(STAGE/'panels/M5_i.pdf',STAGE/'panels/M5_h.pdf')
    palette.panel_i()
    compose_current(5,[(f'M5_{a}',a) for a in 'abcd']+[('M5_efg','efg'),('M5_h','h'),('M5_i','i')])
if __name__=='__main__':globals()['figure'+sys.argv[1]]()
