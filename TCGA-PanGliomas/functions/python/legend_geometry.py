from common import *
from legend_streams import remove_text
CONTENT_CACHE={}
def spans(page):return [s for b in page.get_text("dict")["blocks"] for l in b.get("lines",[]) for s in l["spans"]]

def fragment(path, box):
    """Copy selected vector paths and text at 1:1, without hidden page content.

    Reading the display list once avoids expensive redaction of thousands of
    unrelated scatter markers. Key geometry, colors, strokes and font sizes
    come directly from the existing PDF.
    """
    key = str(path)
    if key not in CONTENT_CACHE:
        with fitz.open(path) as src:
            CONTENT_CACHE[key] = (src[0].get_drawings(), spans(src[0]))
    drawings, texts = CONTENT_CACHE[key]
    doc = fitz.open()
    p = doc.new_page(width=box.width, height=box.height)
    delta = box.tl
    for d in drawings:
        if not box.contains(d['rect']):
            continue
        shape = p.new_shape()
        for item in d['items']:
            kind = item[0]
            if kind == 'l':
                shape.draw_line(item[1] - delta, item[2] - delta)
            elif kind == 'c':
                shape.draw_bezier(*(point - delta for point in item[1:]))
            elif kind == 're':
                r = item[1] + (-delta.x, -delta.y, -delta.x, -delta.y)
                shape.draw_rect(r)
            elif kind == 'qu':
                shape.draw_quad(fitz.Quad(*(point - delta for point in item[1])))
            else:
                raise ValueError(('Unsupported legend path', kind))
        shape.finish(width=d.get('width') or 0, color=d.get('color'), fill=d.get('fill'),
            lineCap=int(max(d.get('lineCap') or (0,))), lineJoin=int(d.get('lineJoin') or 0),
            dashes=d.get('dashes'), closePath=d.get('closePath', False),
            even_odd=d.get('even_odd', False),
            fill_opacity=1 if d.get('fill_opacity') is None else d['fill_opacity'],
            stroke_opacity=1 if d.get('stroke_opacity') is None else d['stroke_opacity'])
        shape.commit()
    for s in texts:
        if not box.contains(fitz.Rect(s['bbox'])):
            continue
        bold = 'Bold' in s['font']
        fontname = 'KEY_BOLD' if bold else 'KEY_REG'
        p.insert_font(fontname=fontname, fontfile=str(FONT / f'RobotoCondensed-{"Bold" if bold else "Regular"}.ttf'))
        color = tuple(((s['color'] >> shift) & 255) / 255 for shift in [16, 8, 0])
        p.insert_text(fitz.Point(s['origin']) - delta, s['text'], fontsize=s['size'], fontname=fontname, color=color)
    return doc


def clear_regions(doc, source_path, boxes):
    """Remove selected text-show operators and cover only local vector keys.

    The established text-stream walker validates each original text position
    and subsequent position reset. Removing these show operations preserves
    all unrelated drawing commands byte-for-byte, including dense violins.
    """
    page = doc[0]
    selected = [s for s in spans(page) if any(box.contains(fitz.Rect(s['bbox'])) for box in boxes)]
    targets = [{'text': s['text'], 'origin_pt': list(s['origin']), 'dx_pt': 0} for s in selected]
    remove_text(doc, source_path, targets, fitz)
    for box in boxes:
        page.draw_rect(box, color=None, fill=(1, 1, 1), width=0)

