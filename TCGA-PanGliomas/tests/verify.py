"""Report scientific checks and visual differences without conflating their status."""
from pathlib import Path
import collections, hashlib, json, re, sys
import pymupdf as fitz
import numpy as np
from PIL import Image
from visual_compare import compare, write_report
REPO=Path(__file__).resolve().parents[1]
def main(out):
    layouts=json.loads((REPO/'functions/styles/figures.json').read_text())
    run=json.loads((out/'run_manifest.json').read_text());figures=run['figures'];reports=[]
    assert run['status']=='RUNNABLE','Wait for the figure build to finish before verification'
    qa=out/'verification';qa.mkdir(exist_ok=True)
    review_path=REPO/'tests/visual_review.json'
    review=json.loads(review_path.read_text()) if review_path.exists() else {'figures':{}}
    fitz.TOOLS.mupdf_warnings(reset=True)
    with fitz.open(out/'figures/main_figures.pdf') as combined:
        assert len(combined)==len(figures)
        for i,n in enumerate(figures):
            with fitz.open(out/f'figures/Figure_{n}.pdf') as actual,fitz.open(REPO/f'reference/figures/Figure_{n}.pdf') as ref:
                assert len(actual)==1
                page=actual[0];base=ref[0]
                fitz.TOOLS.mupdf_warnings(reset=True)
                page.get_text('dict');page.get_pixmap(dpi=150,alpha=False)
                warnings=fitz.TOOLS.mupdf_warnings(reset=True)
                assert not warnings,(n,'Reproduced PDF syntax/rendering warnings',warnings)
                assert '\x00' not in page.get_text(),f'Figure {n} contains invalid font/text encoding'
                assert max(abs(a-b) for a,b in zip(page.rect,base.rect))<.001
                spans=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
                letters=sorted(s['text'] for s in spans if len(s['text'])==1 and abs(s['size']-10)<.01 and 'Bold' in s['font'])
                assert letters==layouts[f'M{n}']['letters'],(n,letters)
                assert not page.get_images(),f'Figure {n} contains embedded raster artwork'
                pix=page.get_pixmap(dpi=150,alpha=False);refpix=base.get_pixmap(dpi=150,alpha=False)
                assert pix.samples==combined[i].get_pixmap(dpi=150,alpha=False).samples
                a=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,3)
                b=np.frombuffer(refpix.samples,np.uint8).reshape(refpix.height,refpix.width,3)
                changed=np.any(a!=b,axis=2);material=np.max(np.abs(a.astype(int)-b.astype(int)),axis=2)>16
                diff=np.full_like(a,255);diff[changed]=[220,30,45]
                Image.fromarray(diff).save(qa/f'Figure_{n}_pixel_difference.png')
                outside=[{'text':s['text'],'bbox':s['bbox']} for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
                numbers=lambda s:collections.Counter(re.findall(r'(?<![\w.])-?\d+(?:\.\d+)?%?',s))
                na,nb=numbers(page.get_text()),numbers(base.get_text())
                assert not outside,(n,'Text outside page',outside)
                assert na==nb,(n,'Numerical text differs',dict(na-nb),dict(nb-na))
                comparison=compare(page,base,n,layouts[f'M{n}'],qa)
                assert hashlib.sha256(combined[i].get_pixmap(dpi=300,alpha=False).samples).hexdigest()==comparison['resolutions']['300']['actual_pixel_sha256'],(n,'Combined PDF differs at 300 dpi')
                approved=review['figures'].get(str(n),{})
                reviewed=review.get('status')=='PASS' and approved.get('status')=='PASS' and approved.get('reference_pdf_sha256')==hashlib.sha256((REPO/f'reference/figures/Figure_{n}.pdf').read_bytes()).hexdigest() and all(
                    approved.get('pixel_sha256',{}).get(dpi)==comparison['resolutions'][dpi]['actual_pixel_sha256'] for dpi in ['150','300'])
                reports.append({'figure':n,'execution':'PASS','one_vector_page':True,'panel_letters':letters,
                    'page_mm':[page.rect.width*25.4/72,page.rect.height*25.4/72],
                    'combined_page_matches':True,'combined_page_matches_dpi':[150,300],'exact_pixels':not changed.any(),
                    'changed_pixel_fraction':float(changed.mean()),'material_pixel_fraction':float(material.mean()),
                    'minimum_text_pt':min(s['size'] for s in spans),'outside_page_text':outside,
                    'extra_numeric_tokens':dict(na-nb),'missing_numeric_tokens':dict(nb-na),
                    'panel_comparison':comparison,'matches_reviewed_visual_fingerprints':reviewed,
                    'visual_status':'REVIEWED_MATCH' if reviewed else 'EXACT' if all(m['exact_pixels'] for m in comparison['resolutions'].values()) else 'DIFFERENCES_REQUIRE_REVIEW'})
    # Manuscript Figure 4 has empty ActualText wrappers left by its original
    # legend edits. They are absent from the regenerated PDF, checked above.
    fitz.TOOLS.mupdf_warnings(reset=True)
    science=json.loads((qa/'statistical_verification.json').read_text());assert science['status']=='PASS'
    extended=json.loads((qa/'extended_statistics.json').read_text())
    assert extended['status']=='PASS' and len(extended['tables'])==15
    assert all(c['pass'] for c in science['checks'].values())
    assert len(extended['permutations'])==3
    ac2=json.loads((qa/'ac2_verification.json').read_text()) if 5 in figures else None
    if ac2 is not None:assert ac2['status']=='PASS' and len(ac2['tables'])>0
    sys.path.insert(0,str(REPO/'functions/python'))
    from release_contract import validate as validate_distribution
    distribution=validate_distribution(REPO)
    result={'execution':'PASS','scientific_checks':science,'figures':reports,
        'ac2_checks':ac2,'extended_statistics':extended,
        'visual_equivalence':'REVIEWED_MATCH' if all(r['matches_reviewed_visual_fingerprints'] for r in reports) else 'EXACT' if all(r['visual_status']=='EXACT' for r in reports) else 'NOT_YET_ESTABLISHED',
        'release_status':'LOCAL_REPRODUCTION_VALIDATED','data_distribution':distribution,'reference_artwork_used_to_render':False,
        'public_data_complete':False,'external_bundle_required':True}
    (qa/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    write_report(reports,qa)
    run['visual_equivalence']=result['visual_equivalence']
    (out/'run_manifest.json').write_text(json.dumps(run,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['figures','scientific_checks','ac2_checks','extended_statistics']},indent=2))
if __name__=='__main__':main(Path(sys.argv[1]).resolve())
