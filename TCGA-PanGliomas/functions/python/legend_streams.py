"""Remove complete selected text runs, including Cairo's combined legend rows."""
import copy
import unicodedata
import pypdf
from pypdf.generic import TextStringObject, ByteStringObject
from pypdf._cmap import get_encoding
from title_streams import parse


def remove_text(doc, path, targets, fitz):
    reader = pypdf.PdfReader(path)
    to_page = doc[0].transformation_matrix
    font_cache, cuts, found = {}, {}, set()
    norm = lambda s: ''.join(unicodedata.normalize('NFKC', s).split())

    def decode(args, op, font):
        key = font.indirect_reference.idnum
        if key not in font_cache:
            font_cache[key] = get_encoding(font)
        enc, cmap = font_cache[key]
        value = ''
        for item in args[0] if op == b'TJ' else [args[-1]]:
            if not isinstance(item, (TextStringObject, ByteStringObject)):
                continue
            raw = item.original_bytes if isinstance(item, TextStringObject) else bytes(item)
            text = raw.decode(enc, errors='replace') if isinstance(enc, str) else ''.join(enc.get(c, chr(c)) for c in raw)
            value += ''.join(cmap.get(c, c) for c in text)
        return value

    def walk(stream, resources, state, stack=None):
        key = stream.indirect_reference.idnum
        resources = stream.get('/Resources', resources)
        state = copy.copy(state)
        stack = [] if stack is None else stack
        tm, tlm = fitz.Matrix(1, 1), fitz.Matrix(1, 1)
        run = []

        def flush():
            if not run:
                return
            text = norm(''.join(r['text'] for r in run))
            origin = run[0]['origin']
            candidates = sorted([(k, t) for k, t in enumerate(targets)
                                 if abs(t['origin_pt'][1]-origin.y) < .12], key=lambda z: z[1]['origin_pt'][0])
            matched = None
            for start, (k, target) in enumerate(candidates):
                if abs(target['origin_pt'][0]-origin.x) > .12:
                    continue
                value, ids = '', []
                for idx, item in candidates[start:]:
                    value += norm(item['text'])
                    ids.append(idx)
                    if value == text:
                        matched = ids
                        break
                    if len(value) > len(text):
                        break
                if matched is not None:
                    break
            if matched is not None:
                found.update(matched)
                cuts.setdefault(key, set()).update((r['start'], r['end']) for r in run)
            run.clear()

        for args, op, start, end in parse(stream.get_data()):
            if op in [b'q', b'Q', b'cm', b'Do', b'BT', b'ET', b'Tm', b'Td', b'TD', b'T*']:
                flush()
            if op == b'q':
                stack.append(copy.copy(state))
            elif op == b'Q':
                state = stack.pop()
            elif op == b'cm':
                state['ctm'] = fitz.Matrix(*map(float, args)) * state['ctm']
            elif op == b'Do':
                sub = resources['/XObject'][args[0]].get_object()
                if sub.get('/Subtype') == '/Form':
                    child = copy.copy(state)
                    child['ctm'] = fitz.Matrix(*map(float, sub.get('/Matrix', [1, 0, 0, 1, 0, 0]))) * state['ctm']
                    walk(sub, resources, child)
            elif op == b'BT':
                tm, tlm = fitz.Matrix(1, 1), fitz.Matrix(1, 1)
            elif op == b'Tf':
                state['font'] = resources['/Font'][args[0]].get_object()
            elif op == b'Tm':
                tm = fitz.Matrix(*map(float, args))
                tlm = fitz.Matrix(tm)
            elif op in [b'Td', b'TD']:
                tx, ty = map(float, args)
                tlm = fitz.Matrix(1, 0, 0, 1, tx, ty) * tlm
                tm = fitz.Matrix(tlm)
                if op == b'TD':
                    state['leading'] = -ty
            elif op == b'TL':
                state['leading'] = float(args[0])
            elif op == b'T*':
                tlm = fitz.Matrix(1, 0, 0, 1, 0, -state.get('leading', 0)) * tlm
                tm = fitz.Matrix(tlm)
            elif op in [b'Tj', b'TJ']:
                run.append({'text': decode(args, op, state['font']), 'start': start, 'end': end,
                            'origin': fitz.Point(tm.e, tm.f) * state['ctm'] * to_page})
            elif op in [b"'", b'"']:
                raise AssertionError('Unsupported implicit line-positioning text operation')
        flush()
        return state

    state, stack = {'ctm': fitz.Matrix(1, 1), 'font': None}, []
    contents = reader.pages[0]['/Contents']
    for stream in contents if isinstance(contents, list) else [contents]:
        state = walk(stream.get_object(), reader.pages[0]['/Resources'], state, stack)
    assert len(found) == len(targets), ('Unmatched text', [t for k, t in enumerate(targets) if k not in found])
    for key, ranges in cuts.items():
        data = doc.xref_stream(key)
        for start, end in sorted(ranges, reverse=True):
            data = data[:start] + b'\n' + data[end:]
        # Cairo wraps ligatures (for example "fi" in "significant") in
        # ActualText metadata. Remove empty wrappers after deleting that run;
        # leaving them creates a phantom unpositioned extraction warning.
        marked, empty = [], []
        for args, op, start, end in parse(data):
            if op in [b'BDC', b'BMC']:
                marked.append({'range':(start,end),'actual':op==b'BDC' and isinstance(args[-1],dict) and '/ActualText' in args[-1], 'text':False})
            elif op in [b'Tj',b'TJ',b"'",b'"']:
                for item in marked:item['text']=True
            elif op==b'EMC':
                item=marked.pop()
                if item['actual'] and not item['text']:empty.extend([item['range'],(start,end)])
        assert not marked,'Unbalanced marked-content wrappers'
        for start,end in sorted(empty,reverse=True):data=data[:start]+b'\n'+data[end:]
        doc.update_stream(key, data)
    return {'text_spans_removed': len(found), 'text_show_operations_removed': sum(map(len, cuts.values()))}
