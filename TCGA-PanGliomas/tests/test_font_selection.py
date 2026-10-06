"""A system-font/cache entry must not override the reviewed Python font files."""
from pathlib import Path
import os
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'functions/python'))
with patch.dict(os.environ, PANGLIOMA_WORKSPACE=str(ROOT / '.local/tmp')):
    from common import setup_matplotlib


class FontSelection(unittest.TestCase):
    def test_competing_font_entry_does_not_change_selected_faces(self):
        from matplotlib import font_manager
        # This other version has the same family and style, and won in one
        # fresh font cache before explicit selection was implemented.
        font_manager.fontManager.addfont(str(ROOT / 'functions/fonts/RobotoCondensed-Regular.ttf'))
        setup_matplotlib()
        expected = {
            ('normal', 'normal'): 'system/RobotoCondensed-VariableFont_wght.ttf',
            ('bold', 'normal'): 'RobotoCondensed-Bold.ttf',
            ('normal', 'italic'): 'RobotoCondensed-Italic.ttf',
            ('bold', 'italic'): 'system/RobotoCondensed-SemiBoldItalic.ttf',
        }
        for (weight, style), relative in expected.items():
            prop = font_manager.FontProperties(family='Roboto Condensed', weight=weight, style=style)
            self.assertEqual(Path(font_manager.findfont(prop)).resolve(), (ROOT / 'functions/fonts' / relative).resolve())


if __name__ == '__main__':
    unittest.main()
