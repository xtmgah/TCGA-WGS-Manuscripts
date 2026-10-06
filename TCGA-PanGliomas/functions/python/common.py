"""Portable paths and shared native-vector drawing utilities."""
from pathlib import Path
import hashlib, json, os
import pymupdf as fitz

REPO = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ['PANGLIOMA_WORKSPACE']).resolve()
WORK = ROOT
STAGE = ROOT/'stage'
FONT = REPO/'functions/fonts'
MM = 72/25.4
STYLE = json.loads(Path(__file__).with_name('style.json').read_text())
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def text(page,x,y,label,size=6.5,bold=False,color=(.13,.13,.13)):
    # Reuse the key's regular-font resource name. Registering the same font
    # under both PLKEY and RC after PDF grafting can corrupt glyph encoding.
    name='RCB' if bold else 'PLKEY'
    page.insert_font(fontname=name,fontfile=str(FONT/('RobotoCondensed-Bold.ttf' if bold else 'RobotoCondensed-Regular.ttf')))
    page.insert_text((x*MM,y*MM),label,fontsize=size,fontname=name,color=color)
def setup_matplotlib():
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager
    # Family lookup otherwise depends on system-font discovery/cache order.
    # These are the exact faces used by the reviewed manuscript render.
    font_manager.fontManager.ttflist=[f for f in font_manager.fontManager.ttflist if f.name!='Roboto Condensed']
    for p in [FONT/'system/RobotoCondensed-VariableFont_wght.ttf',FONT/'RobotoCondensed-Bold.ttf',
              FONT/'RobotoCondensed-Italic.ttf',FONT/'system/RobotoCondensed-SemiBoldItalic.ttf']:
        font_manager.fontManager.addfont(str(p))
    matplotlib.rcParams.update({'font.family':'Roboto Condensed','font.size':6.5,'axes.labelsize':6.5,
        'axes.titlesize':8,'axes.titlelocation':'left','xtick.labelsize':6.5,'ytick.labelsize':6.5,
        'legend.fontsize':6.5,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,
        'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
def compose(number,title,height_mm,placements):
    doc=fitz.open();page=doc.new_page(width=180*MM,height=height_mm*MM)
    text(page,4,8,title,12,True)
    for letter,x,y,w,h in placements:
        with fitz.open(STAGE/f'panels/M{number}_{letter}.pdf') as src:
            page.show_pdf_page(fitz.Rect(x*MM,(y+3.5)*MM,(x+w)*MM,(y+h)*MM),src,0)
        text(page,x,y+2.8,letter,10,True)
    out=STAGE/f'main_figures/Figure_{number}.pdf';doc.save(out,garbage=4,deflate=True);doc.close()
    return out
