"""Blank-font compatibility must preserve the figure's marks and text."""
from pathlib import Path
import sys
import unittest
import pymupdf as fitz
import pypdf
from pypdf.generic import DecodedStreamObject

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'functions/python'))
from cairo_spaces import blank_spaces, normalize

FIXTURE = ROOT / 'tests/fixtures/cairo_type3_spaces.pdf'


class CairoSpaces(unittest.TestCase):
    def test_blank_spaces_preserve_pixels_and_labels(self):
        with fitz.open(FIXTURE) as doc:
            original = [doc[0].get_pixmap(dpi=dpi, alpha=False).samples for dpi in [150, 300]]
            text = doc[0].get_text()
            self.assertGreater(normalize(doc, FIXTURE), 0)
            # Reopen to test embedded-font serialization as well as live rendering.
            with fitz.open(stream=doc.tobytes(), filetype='pdf') as result:
                self.assertEqual(result[0].get_text(), text)
                spans = [s for b in result[0].get_text('dict')['blocks'] for line in b.get('lines', []) for s in line['spans']]
                self.assertIn('Multivariable Cox model for OS', [s['text'] for s in spans])
                for dpi, pixels in zip([150, 300], original):
                    self.assertEqual(result[0].get_pixmap(dpi=dpi, alpha=False).samples, pixels)
                self.assertEqual(len(result[0].get_drawings()), len(doc[0].get_drawings()))

    def test_painted_type3_glyph_is_not_a_blank_space(self):
        reader = pypdf.PdfReader(FIXTURE)
        fonts = reader.pages[0]['/Resources']['/Font']
        font = next(ref.get_object() for ref in fonts.values() if ref.get_object()['/Subtype'] == '/Type3')
        self.assertIsNotNone(blank_spaces(font))
        proc = DecodedStreamObject()
        proc.set_data(b'0.23 0 0 0 0 0 d1\n0 0 1 1 re f\n')
        key = next(iter(font['/CharProcs']))
        font['/CharProcs'][key] = proc
        self.assertIsNone(blank_spaces(font))

    def test_pdf_without_type3_is_unchanged(self):
        with fitz.open() as doc:
            doc.new_page().insert_text((20, 20), 'No special font')
            before = doc.tobytes(no_new_id=True)
            self.assertEqual(normalize(doc, Path('unused-file.pdf')), 0)
            self.assertEqual(doc.tobytes(no_new_id=True), before)


if __name__ == '__main__':
    unittest.main()
