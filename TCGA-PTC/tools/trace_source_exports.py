#!/usr/bin/env python3
"""Read-only source search. Review dynamic/related matches manually; no inferred mappings.

Usage: python3 tools/trace_source_exports.py SOURCE_ROOT OUTPUT_TSV
The source root is the private historical R tree, not this cleaned package.
"""
import csv,re,sys
from pathlib import Path
root=Path(sys.argv[1]);out=Path(sys.argv[2])
ptc=Path(__file__).resolve().parents[1]
with (ptc/'provenance/keynote_panel_map.tsv').open() as h:panels=list(csv.DictReader(h,delimiter='\t'))
files=[p for p in root.rglob('*') if p.suffix.lower() in ('.r','.rmd')]
rows=[]
for p in panels:
 name=p['original_plot_filename']
 if not name:continue
 stem=Path(name).stem
 matched=False
 for f in files:
  for i,line in enumerate(f.read_text(errors='replace').splitlines(),1):
   if name in line or stem in line:
    rows.append(dict(manuscript_id=p['manuscript_figure_panel_id'],original_plot_filename=name,source=f.relative_to(root).as_posix(),line=i,evidence=line.strip(),match='literal' if name in line else 'stem_candidate'))
    matched=True
 if not matched:rows.append(dict(manuscript_id=p['manuscript_figure_panel_id'],original_plot_filename=name,source='',line='',evidence='No literal/stem match. Review dynamic export expressions, upstream object definitions and documented manual mapping.',match='not_found'))
with out.open('w') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(len(files),'scripts searched;',len(rows),'evidence rows')
