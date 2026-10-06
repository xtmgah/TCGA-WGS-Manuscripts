"""Reproduce from the manifest in an isolated workspace; never read the source project."""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, os, shutil, subprocess, sys, time
REPO=Path(__file__).resolve().parents[2]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def manifest():
    with (REPO/'data/manifest.tsv').open() as f:return list(csv.DictReader(f,delimiter='\t'))
def location(row,data):
    return data/Path(row['path']).relative_to('data/external') if row['tier']=='external' else REPO/row['path']
def check_inputs(data,public_only=False):
    errors=[];rows=[r for r in manifest() if not public_only or r['tier']=='processed']
    for row in rows:
        p=location(row,data)
        if not p.is_file():errors.append(f"Missing {row['id']}: {p}");continue
        if sha(p)!=row['compressed_sha256']:errors.append(f"Checksum mismatch: {p}");continue
        if hashlib.sha256(gzip.decompress(p.read_bytes())).hexdigest()!=row['sha256']:errors.append(f"Content mismatch: {p}")
    if errors:raise RuntimeError('\n'.join(errors)+'\nFor request-only files, contact the corresponding authors and follow data/external/README.md. Box/Biowulf are not public download locations.')
    return {'status':'PASS','inputs':len(rows),'external_inputs':sum(r['tier']=='external' for r in rows)}
def run(cmd,cwd,log,env):
    print('Running '+Path(cmd[0]).name+' '+Path(cmd[-1]).name,flush=True)
    with log.open('w') as f:p=subprocess.run(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:
        print(log.read_text()[-6500:],file=sys.stderr)
        raise RuntimeError(f'Command failed ({p.returncode}); see {log}')
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--figures',default='1-5');p.add_argument('--data-dir',type=Path,default=REPO/'data/external')
    p.add_argument('--output-dir',type=Path,default=REPO/'outputs');p.add_argument('--check-inputs',action='store_true')
    p.add_argument('--check-public-inputs',action='store_true',help='Check only the Git-distributed inputs; no private data required')
    p.add_argument('--verify',action='store_true');p.add_argument('--render-reports',action='store_true')
    p.add_argument('--rscript',default=os.environ.get('PANGLIOMA_RSCRIPT','Rscript'))
    a=p.parse_args();a.data_dir=a.data_dir.resolve();a.output_dir=a.output_dir.resolve()
    info=check_inputs(a.data_dir,public_only=a.check_public_inputs)
    if a.check_public_inputs:
        print(json.dumps(dict(info,scope='Git inputs only; full reproduction also requires author-request data'),indent=2));return
    if a.check_inputs:print(json.dumps(info,indent=2));return
    nums=[]
    for item in a.figures.split(','):
        nums.extend(range(int(item.split('-')[0]),int(item.split('-')[1])+1) if '-' in item else [int(item)])
    nums=sorted(set(nums))
    if not nums or not set(nums)<=set(range(1,6)):raise ValueError('Figures must be selected from 1–5')
    out=a.output_dir;work=out/'.work';stage=work/'stage'
    if not out.is_relative_to(REPO):raise ValueError('The output directory must be inside this standalone folder.')
    out.mkdir(parents=True,exist_ok=True)
    (out/'run_manifest.json').write_text(json.dumps({'status':'STARTED','figures':nums},indent=2)+'\n')
    for d in ['panels','derived','qa','main_figures','captions','native_top']:(stage/d).mkdir(parents=True,exist_ok=True)
    for d in ['figures','tables','logs','verification','tmp']:(out/d).mkdir(parents=True,exist_ok=True)
    # Materialize only the manifest inputs. No reference PDF or original code is copied here.
    for row in manifest():
        dest=work/row['logical_path'];dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_bytes(gzip.decompress(location(row,a.data_dir).read_bytes()))
        if '/staging/derived/' in row['logical_path']:
            shutil.copy2(dest,stage/'derived'/dest.name)
    env=os.environ.copy();env.update(PANGLIOMA_WORKSPACE=str(work),PANGLIOMA_R_CODE=str(REPO/'functions/R'),
        PANGLIOMA_PY_CODE=str(REPO/'functions/python'),PANGLIOMA_REPO=str(REPO),
        TCGA_PROJECT_ROOT=str(work),TMPDIR=str(out/'tmp'),MPLCONFIGDIR=str(out/'tmp/matplotlib'),
        PYTHONDONTWRITEBYTECODE='1',XDG_CACHE_HOME=str(out/'tmp/cache'),R_HISTFILE=str(out/'tmp/Rhistory'))
    start=time.time()
    run([a.rscript,'--vanilla',str(REPO/'functions/R/recompute.R')],work,out/'logs/statistics.log',env)
    if 5 in nums:
        run([sys.executable,str(REPO/'functions/python/ac2_aggregate.py')],work,out/'logs/ac2_aggregate.log',env)
        run([a.rscript,'--vanilla',str(REPO/'functions/R/ac2_statistics.R')],work,out/'logs/ac2_statistics.log',env)
        run([sys.executable,str(REPO/'functions/python/ac2_validate.py')],work,out/'logs/ac2_validate.log',env)
    scripts={1:[],2:['render_figure2.R'],3:['render_figure3.R'],
             4:['render_figure4.R','figure4_wider_survival.R','figure4_four_metric_row.R'],5:[]}
    for n in nums:
        for s in scripts[n]:run([a.rscript,'--vanilla',str(REPO/'functions/R'/s),str(work)],work,out/'logs'/f'{n}_{s}.log',env)
        run([sys.executable,str(REPO/'functions/python/render.py'),str(n)],work,out/'logs'/f'Figure_{n}.log',env)
        shutil.copy2(stage/f'main_figures/Figure_{n}.pdf',out/f'figures/Figure_{n}.pdf')
    from pymupdf import open as pdfopen
    combined=pdfopen()
    for n in nums:
        with pdfopen(out/f'figures/Figure_{n}.pdf') as d:
            d[0].get_pixmap(dpi=150,alpha=False).save(out/f'figures/Figure_{n}.png');combined.insert_pdf(d)
    combined.save(out/'figures/main_figures.pdf',garbage=4,deflate=True);combined.close()
    for f in (stage/'derived').glob('*'):
        if f.is_file():shutil.copy2(f,out/'tables'/f.name)
    for f in (stage/'qa').glob('*.json'):shutil.copy2(f,out/'verification'/f.name)
    result={'status':'RUNNABLE','figures':nums,'input_checks':info,'elapsed_seconds':round(time.time()-start,2),
        'visual_equivalence':'not yet established','upstream_models':'frozen ordering/mixture/timing estimates, molecular/RNA classification, DESeq2/GSEA results',
        'source_project_required':False}
    (out/'run_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    if a.verify:
        run([sys.executable,str(REPO/'tests/verify.py'),str(out)],REPO,out/'logs/verification.log',env)
        result=json.loads((out/'run_manifest.json').read_text())
    if a.render_reports:run([a.rscript,'--vanilla',str(REPO/'functions/R/reports.R'),str(out)],REPO,out/'logs/reports.log',env)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
