"""Panel-level visual diagnostics; reference artwork is read only during QA."""
from pathlib import Path
import collections, hashlib, html, json, math
import numpy as np
import pymupdf as fitz
from PIL import Image

MM = 72 / 25.4

def spans(page):
    return [s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]

def text_comparison(page, reference):
    actual, expected = spans(page), spans(reference)
    key = lambda s: (' '.join(s['text'].split()), s['font'], round(s['size'], 2))
    aa, bb = collections.Counter(map(key,actual)), collections.Counter(map(key,expected))
    unused = set(range(len(actual))); offsets=[]; colors=[]
    for s in expected:
        candidates=[i for i in unused if key(actual[i])==key(s)]
        if not candidates:continue
        i=min(candidates,key=lambda j:sum((actual[j]['origin'][k]-s['origin'][k])**2 for k in (0,1)))
        unused.remove(i);t=actual[i]
        offsets.append(max(abs(t['origin'][k]-s['origin'][k]) for k in (0,1)))
        if t['color']!=s['color']:colors.append({'text':s['text'],'actual':t['color'],'reference':s['color']})
    return {'actual_spans':len(actual),'reference_spans':len(expected),
            'extra_text_spans':[{'text':k[0],'font':k[1],'size':k[2],'count':v} for k,v in (aa-bb).items()],
            'missing_text_spans':[{'text':k[0],'font':k[1],'size':k[2],'count':v} for k,v in (bb-aa).items()],
            'max_origin_difference_pt':max(offsets,default=0),'color_differences':colors}

def metric(a,b):
    delta=np.max(np.abs(a.astype(np.int16)-b.astype(np.int16)),axis=2)
    return {'exact_pixels':not bool(delta.any()),'changed_pixel_fraction':float(np.mean(delta>0)),
            'material_pixel_fraction':float(np.mean(delta>16)),
            'maximum_channel_difference':int(delta.max()),'pixels':int(delta.size)}

def compare(page,reference,n,spec,qa):
    result={'text':text_comparison(page,reference),'resolutions':{},'panels':[]}
    images=qa/'panels';images.mkdir(exist_ok=True)
    for dpi in (150,300):
        pix,refpix=[p.get_pixmap(dpi=dpi,alpha=False) for p in (page,reference)]
        a=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,3)
        b=np.frombuffer(refpix.samples,np.uint8).reshape(refpix.height,refpix.width,3)
        result['resolutions'][str(dpi)]=dict(metric(a,b),actual_pixel_sha256=hashlib.sha256(pix.samples).hexdigest(),
                                            reference_pixel_sha256=hashlib.sha256(refpix.samples).hexdigest())
        for index,(letter,(x,y,w,h)) in enumerate(spec['panels'].items()):
            bounds=[math.floor(x*dpi/25.4),math.floor(y*dpi/25.4),math.ceil((x+w)*dpi/25.4),math.ceil((y+h)*dpi/25.4)]
            left,top,right,bottom=bounds;aa=a[top:bottom,left:right];bb=b[top:bottom,left:right]
            if dpi==150:result['panels'].append({'panel':f'{n}{letter}','bounds_mm':[x,y,w,h],'resolutions':{}})
            row=result['panels'][index];row['resolutions'][str(dpi)]=metric(aa,bb)
            if dpi==150:
                for label,pixels in [('reference',bb),('reproduction',aa)]:
                    Image.fromarray(pixels).save(images/f'{n}{letter}_{label}.png')
                diff=np.full_like(aa,255);diff[np.any(aa!=bb,axis=2)]=[220,30,45]
                Image.fromarray(diff).save(images/f'{n}{letter}_difference.png')
        # Whole-page differences include legends and headings outside panel bounds.
        diff=np.full_like(a,255);diff[np.any(a!=b,axis=2)]=[220,30,45]
        Image.fromarray(diff).save(qa/f'Figure_{n}_difference_{dpi}dpi.png')
    return result

def write_report(reports,qa):
    content=['<!doctype html><html lang="en"><meta charset="utf-8"><title>Pan-Glioma panel comparison</title>',
      '<style>body{font:16px system-ui;margin:32px;color:#222}section{border-top:1px solid #ccc;padding:20px 0;break-inside:avoid}.row{display:flex;gap:12px}.row figure{width:33%;margin:0}.row img{width:100%;height:auto;border:1px solid #ddd}figcaption{font-weight:600;margin-bottom:8px}p{max-width:1000px}nav a{margin-right:12px}summary{cursor:pointer}</style>',
      '<h1>Main figures: panel-level comparison</h1><p>Reference manuscript, regenerated standalone figure, and exact pixel differences (red). Images use 150 dpi; metrics include 150 and 300 dpi. Pixel differences measure rendering, not scientific validity. A reviewed status applies only to the recorded pixel fingerprints; changed fingerprints require a new review.</p>',
      '<nav>'+''.join(f'<a href="#figure-{r["figure"]}">Figure {r["figure"]}</a>' for r in reports)+'</nav>']
    for r in reports:
        n=r['figure'];v=r['panel_comparison'];t=v['text']
        content.append(f'<h2 id="figure-{n}">Figure {n}: {html.escape(r["visual_status"])}</h2><p>Maximum matched-text origin difference: {t["max_origin_difference_pt"]:.6f} pt. Missing text spans: {len(t["missing_text_spans"])}; extra spans: {len(t["extra_text_spans"])}; text color differences: {len(t["color_differences"])}.</p>')
        content.append(f'<details><summary>Whole-page difference maps</summary><a href="Figure_{n}_difference_150dpi.png">150 dpi</a> · <a href="Figure_{n}_difference_300dpi.png">300 dpi</a></details>')
        for panel in v['panels']:
            p=panel['panel'];m150=panel['resolutions']['150'];m300=panel['resolutions']['300']
            content.append(f'<section id="panel-{p}"><h3>Panel {p}</h3><p>Changed pixels: {m150["changed_pixel_fraction"]:.5%} at 150 dpi; {m300["changed_pixel_fraction"]:.5%} at 300 dpi. Differences &gt;16/255: {m300["material_pixel_fraction"]:.5%} at 300 dpi.</p><div class="row">')
            for label in ['reference','reproduction','difference']:
                content.append(f'<figure><figcaption>{label.title()}</figcaption><a href="panels/{p}_{label}.png"><img loading="lazy" src="panels/{p}_{label}.png" alt="Panel {p} {label}"></a></figure>')
            content.append('</div></section>')
    content.append('</html>');(qa/'panel_review.html').write_text('\n'.join(content))
