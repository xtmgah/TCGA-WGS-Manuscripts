"""Create separate local code and processed-input archives; never upload or invoke Git."""
from pathlib import Path
import csv,gzip,hashlib,io,json,tarfile
from release_contract import validate
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def archive(path,entries):
    # Stable names, permissions and timestamps make the archive auditable.
    with path.open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as zipped,tarfile.open(fileobj=zipped,mode='w') as tar:
        for source,name in entries:
            if source.is_symlink():raise ValueError(f'Symlinks are not allowed: {source}')
            contents=source.read_bytes();info=tarfile.TarInfo(name);info.size=len(contents);info.mode=0o644;info.mtime=0
            tar.addfile(info,io.BytesIO(contents))
def main():
    distribution=validate(ROOT)
    target=ROOT/'.local/bundles';target.mkdir(parents=True,exist_ok=True)
    paths=[ROOT/x for x in ['README.md','.gitignore','requirements.txt','data/manifest.tsv','data/README.md','data/external/README.md','data/input_review.json','data/redistribution_audit.tsv','data/column_audit.tsv']]
    for folder in ['Rscripts','functions','environment','docs','reference','tests']:
        paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and not any(x in p.parts for x in ['__pycache__','.DS_Store']) and p.suffix not in ['.pyc','.pyo'])
    paths=sorted(set(paths));rows=list(csv.DictReader((ROOT/'data/manifest.tsv').open(),delimiter='\t'))
    paths.extend(ROOT/row['path'] for row in rows if row['tier']=='processed')
    paths=sorted(set(paths))
    external=[]
    for row in rows:
        p=ROOT/row['path']
        if not p.is_file() or digest(p)!=row['compressed_sha256']:raise RuntimeError(f'Missing or changed input: {row["id"]}')
        if row['tier']=='external':external.append((p,str(p.relative_to(ROOT/'data/external'))))
    for n in range(1,6):
        if not (ROOT/f'Rscripts/Fig{n}.html').is_file():raise RuntimeError('Render the five HTML reports before packaging.')
    # Nothing with individual-level data, working files, dependencies or local configuration enters this archive.
    forbidden=['.local/','outputs/','data/external/files/','.venv/','.Rlib/']
    entries=[]
    for p in paths:
        name=str(p.relative_to(ROOT))
        if any(name.startswith(x) for x in forbidden):raise RuntimeError(name)
        if p.stat().st_size>=100_000_000:raise RuntimeError(f'Unexpected large code/repository file: {name}')
        entries.append((p,name))
    manifest=[{'path':n,'bytes':p.stat().st_size,'sha256':digest(p)} for p,n in entries]
    (target/'code_file_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    archive(target/'TCGA-PanGliomas-code.tar.gz',entries)
    archive(target/'panglioma-request-data.tar.gz',external)
    report={'code_files':len(entries),'external_files':len(external),'archives':[
        {'file':p.name,'bytes':p.stat().st_size,'sha256':digest(p)} for p in [target/'TCGA-PanGliomas-code.tar.gz',target/'panglioma-request-data.tar.gz']],
        'distribution_review':distribution,'public_archive_contains_request_data':False,
        'request_bundle_is_private':True,'raw_sequence_data_included':False,
        'data_revision':'main-figures-1-5-data-review-v1','study_specific_access_authorization_certified':False}
    (target/'bundle_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
