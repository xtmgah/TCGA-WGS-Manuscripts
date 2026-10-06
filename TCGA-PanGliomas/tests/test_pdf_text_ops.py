"""Regression checks for subpixel PDF translations and removed ligature metadata."""
from pathlib import Path
import sys,tempfile,unittest
import pymupdf as fitz
REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'functions/python'))
from title_streams import shift
from legend_streams import remove_text

class TextOperations(unittest.TestCase):
    def setUp(self):
        root=REPO/'outputs/tmp';root.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=root);self.root=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def test_tiny_translation_is_valid_pdf_and_preserves_other_text(self):
        source=self.root/'source.pdf'
        with fitz.open() as d:
            p=d.new_page(width=200,height=150);p.insert_text((20,30),'Heading',fontsize=8)
            p.insert_text((20,60),'Unchanged',fontsize=8);d.save(source)
        with fitz.open(source) as d:
            shift(d,source,[{'text':'Heading','origin_pt':[20,30],'dx_pt':.00006103515625,'dy_pt':10.00006103515625}],fitz)
            output=self.root/'moved.pdf';d.save(output)
        fitz.TOOLS.mupdf_warnings(reset=True)
        with fitz.open(output) as d:
            spans={s['text']:s for b in d[0].get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']}
            self.assertAlmostEqual(spans['Heading']['origin'][0],20.00006103515625,places=4)
            self.assertAlmostEqual(spans['Heading']['origin'][1],40.00006103515625,places=4)
            self.assertEqual(spans['Unchanged']['origin'],(20,60))
            d[0].get_pixmap()
        self.assertEqual(fitz.TOOLS.mupdf_warnings(reset=True),'')
    def test_removed_ligature_does_not_leave_phantom_actualtext(self):
        source=self.root/'source.pdf'
        with fitz.open() as d:
            p=d.new_page(width=200,height=150);p.insert_text((20,30),'fi',fontsize=8)
            x=p.get_contents()[0];data=d.xref_stream(x)
            d.update_stream(x,b'/Span << /ActualText <feff00660069> >> BDC\n'+data+b'\nEMC\n')
            p.insert_text((20,60),'Retained',fontsize=8);d.save(source)
        with fitz.open(source) as d:
            remove_text(d,source,[{'text':'fi','origin_pt':[20,30]}],fitz)
            output=self.root/'removed.pdf';d.save(output)
        fitz.TOOLS.mupdf_warnings(reset=True)
        with fitz.open(output) as d:self.assertEqual(d[0].get_text().strip(),'Retained')
        self.assertEqual(fitz.TOOLS.mupdf_warnings(reset=True),'')
if __name__=='__main__':unittest.main()
