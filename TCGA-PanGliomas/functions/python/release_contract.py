"""Enforce the reviewed Git/request-only split before packaging."""
from pathlib import Path
import csv,gzip,hashlib,io,json,re
PARTICIPANT=re.compile(rb'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}(?:-|\b)')
INDIVIDUAL_COLUMNS={'Subject','Tumor_Barcode','Tumour_Name','Tumor_Barcodes','sample_id','Sample name','Normal_Barcode'}
GIT_CLASSES={'git_figure_support_aggregate','git_project_style'}
REQUEST_CLASSES={'request_measured_data','request_reference_provenance'}
def verify_public_payload(raw,columns):
    if PARTICIPANT.search(raw):raise ValueError('Participant/sample identifier in proposed Git data')
    if set(columns)&INDIVIDUAL_COLUMNS:raise ValueError('Individual-record column in proposed Git data')
def validate(root):
    root=Path(root)
    with (root/'data/manifest.tsv').open() as handle:
        manifest=list(csv.DictReader(handle,delimiter='\t'))
    review=json.loads((root/'data/input_review.json').read_text())
    plans={x['id']:x for x in review['inputs']}
    assert len(plans)==len(review['inputs']),'Duplicate review ID'
    assert {x['id'] for x in manifest}=={i for i,p in plans.items() if p['route']!='retired'},'Unreviewed, missing or retired input in active manifest'
    assert len(manifest)==len({x['path'] for x in manifest})==len({x['id'] for x in manifest})
    public=[];private=[]
    for row in manifest:
        plan=plans[row['id']];route=plan['route'];tier='processed' if route=='git' else 'external'
        assert row['tier']==tier and Path(row['path']).parent==Path('data')/tier/'files',row['id']
        assert row['redistribution']==plan['classification'],row['id']
        assert row['redistribution'] in (GIT_CLASSES if route=='git' else REQUEST_CLASSES),row['id']
        path=root/row['path'];assert not path.is_symlink(),path
        compressed=path.read_bytes();raw=gzip.decompress(compressed)
        assert hashlib.sha256(compressed).hexdigest()==row['compressed_sha256'],row['id']
        assert hashlib.sha256(raw).hexdigest()==row['sha256'],row['id']
        delimiter=',' if row['logical_path'].endswith('.csv') else '\t'
        reader=csv.reader(io.StringIO(raw.decode()),delimiter=delimiter);columns=next(reader)
        assert columns==json.loads(row['columns']),row['id']
        n=sum(1 for _ in reader);assert n==int(row['rows']),row['id']
        if route=='git':verify_public_payload(raw,columns);public.append(path)
        else:private.append(path)
    for tier,allowed in [('processed',public),('external',private)]:
        present={p for p in (root/f'data/{tier}/files').rglob('*') if p.is_file()}
        assert present==set(allowed),f'Unlisted/missing {tier} file(s): {present ^ set(allowed)}'
    return {'status':'PASS','active_inputs':len(manifest),'reviewed_inputs':len(plans),
            'git_inputs':len(public),'request_inputs':len(private),'git_identifier_scan':'PASS',
            'unknown_or_unlisted_inputs':0,'public_storage':'Git only','private_storage':'Author-held Box/Biowulf'}
if __name__=='__main__':
    print(json.dumps(validate(Path(__file__).resolve().parents[2]),indent=2))
