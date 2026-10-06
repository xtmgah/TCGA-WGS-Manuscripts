"""Normalize Cairo's blank Type3 space glyphs without changing drawn marks.

Some Cairo builds put spaces in a separate, empty Type3 font. Restoring that
empty glyph to the adjacent embedded TrueType subset keeps text in one font
run for the title/legend layout code. All drawn glyph outlines are unchanged.
"""
from io import BytesIO
import pypdf
from pypdf._cmap import get_encoding
from pypdf.generic import ArrayObject, ByteStringObject, FloatObject, TextStringObject
from fontTools.ttLib import TTFont
from fontTools.pens.ttGlyphPen import TTGlyphPen
from title_streams import parse


def raw_string(value):
    return value.original_bytes if isinstance(value, TextStringObject) else bytes(value)


def blank_spaces(font):
    if font.get('/Subtype') != '/Type3':
        return None
    matrix = list(map(float, font.get('/FontMatrix', [])))
    if len(matrix) != 6 or matrix[0] <= 0 or any(matrix[i] != 0 for i in [1, 2, 4, 5]):
        return None
    # Never replace a Type3 glyph that draws anything, even if labelled a space.
    for proc in font['/CharProcs'].values():
        ops = parse(proc.get_object().get_data())
        if len(ops) != 1 or ops[0][1] not in [b'd0', b'd1']:
            return None
        args, op = ops[0][:2]
        if float(args[1]) != 0 or (op == b'd1' and any(float(v) != 0 for v in args[2:])):
            return None
    enc, cmap = get_encoding(font)
    widths = {}
    for code, width in enumerate(font['/Widths'], int(font['/FirstChar'])):
        char = bytes([code]).decode(enc) if isinstance(enc, str) else enc.get(code, chr(code))
        if cmap.get(char, char) != ' ':
            return None
        widths[code] = float(width) * matrix[0] * 1000
    return widths


def normalize(doc, path):
    if not any(font[2] == 'Type3' for page in doc for font in page.get_fonts(full=True)):
        return 0
    reader = pypdf.PdfReader(path)
    changed = 0
    visited = set()
    patched_fonts = {}

    def add_space(font, width):
        key = font.indirect_reference.idnum
        if key in patched_fonts:
            assert abs(patched_fonts[key] - width) < .001
            return
        assert font['/Subtype'] == '/TrueType' and font['/Encoding'] == '/WinAnsiEncoding'
        first = int(font['/FirstChar'])
        assert first <= 32 <= int(font['/LastChar'])
        program = font['/FontDescriptor']['/FontFile2']
        ttf = TTFont(BytesIO(program.get_data()), recalcTimestamp=False)
        assert ttf.getBestCmap().get(32) is None, 'Expected an omitted space glyph in the Cairo subset'
        name = 'space'
        assert name not in ttf.getGlyphOrder()
        ttf.setGlyphOrder(ttf.getGlyphOrder() + [name])
        ttf['glyf'].glyphs[name] = TTGlyphPen(None).glyph()
        ttf['hmtx'].metrics[name] = (round(width * ttf['head'].unitsPerEm / 1000), 0)
        for table in ttf['cmap'].tables:
            table.cmap[32] = name
        buffer = BytesIO(); ttf.save(buffer); data = buffer.getvalue(); ttf.close()
        doc.update_stream(program.indirect_reference.idnum, data)
        doc.xref_set_key(program.indirect_reference.idnum, 'Length1', str(len(data)))
        widths = ArrayObject(font['/Widths'])
        widths[32 - first] = FloatObject(width)
        buffer = BytesIO(); widths.write_to_stream(buffer)
        doc.xref_set_key(key, 'Widths', buffer.getvalue().decode())
        patched_fonts[key] = width

    def walk(stream, resources):
        nonlocal changed
        key = stream.indirect_reference.idnum
        if key in visited:
            return
        visited.add(key)
        resources = stream.get('/Resources', resources)
        fonts = resources.get('/Font', {})
        spaces = {name: value for name, ref in fonts.items()
                  if (value := blank_spaces(ref.get_object())) is not None}
        for ref in resources.get('/XObject', {}).values():
            child = ref.get_object()
            if child.get('/Subtype') == '/Form':
                walk(child, resources)
        if not spaces:
            return
        source = stream.get_data()
        output, pending = [], []
        active = emitted = None
        converted = 0
        spacing = {b'Tc': 0., b'Tw': 0.}

        def flush():
            if pending:
                buffer = BytesIO()
                ArrayObject(pending).write_to_stream(buffer)
                output.append(buffer.getvalue() + b' TJ\n')
                pending.clear()

        operations = parse(source)
        for index, (args, op, start, end) in enumerate(operations):
            if op == b'Tf':
                active = (args[0], float(args[1]))
                if active[0] in spaces and emitted is not None and emitted[0] not in spaces:
                    target = emitted
                    if fonts[target[0]].get_object().get('/Subtype') != '/TrueType':
                        # Put the space after a Unicode symbol in the following
                        # text font, preserving a real extractable space glyph.
                        for following, operator, begin, finish in operations[index + 1:]:
                            if operator == b'Tf':
                                candidate = (following[0], float(following[1]))
                                if fonts[candidate[0]].get_object().get('/Subtype') == '/TrueType' and candidate[1] == active[1]:
                                    flush()
                                    output.append(source[begin:finish] + b'\n')
                                    target = emitted = candidate
                                break
                            if operator not in [b'Tj', b'TJ']:
                                break
                    if fonts[target[0]].get_object().get('/Subtype') == '/TrueType':
                        continue
                if active != emitted:
                    flush()
                    output.append(source[start:end] + b'\n')
                    emitted = active
            elif op in [b'Tj', b'TJ']:
                values = args[0] if op == b'TJ' else [args[0]]
                if active is not None and active[0] in spaces and emitted != active:
                    assert not any(spacing.values()), 'Unsupported character/word spacing for blank Type3 normalization'
                    assert active[1] == emitted[1], 'Unexpected font-size change around a blank space'
                    target_font = fonts[emitted[0]].get_object()
                    assert target_font.get('/Subtype') == '/TrueType'
                    for value in values:
                        if isinstance(value, (TextStringObject, ByteStringObject)):
                            for code in raw_string(value):
                                width = spaces[active[0]][code]
                                add_space(target_font, width)
                                pending.append(ByteStringObject(b' '))
                                converted += 1
                        else:
                            pending.append(value)
                else:
                    pending.extend(ByteStringObject(raw_string(v)) if isinstance(v, (TextStringObject, ByteStringObject)) else v for v in values)
            else:
                flush()
                # A skipped blank font selection must not leak across a reset.
                if active is not None and active != emitted:
                    output.append(f'{active[0]} {active[1]:.12f} Tf\n'.encode())
                    emitted = active
                output.append(source[start:end] + b'\n')
                if op in spacing:
                    spacing[op] = float(args[0])
                if op in [b'BT', b'ET', b'q', b'Q']:
                    active = emitted = None
        flush()
        if converted:
            doc.update_stream(key, b''.join(output))
            changed += converted

    for page in reader.pages:
        contents = page['/Contents']
        for ref in contents if isinstance(contents, list) else [contents]:
            walk(ref.get_object(), page['/Resources'])
    return changed
