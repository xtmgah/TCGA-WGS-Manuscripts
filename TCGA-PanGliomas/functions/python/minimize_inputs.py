"""Maintainer migration from a local, checksummed pre-review input snapshot.

Never reads the original analysis project or downloads/uploads data. Retained
cells are selected as strings so floating-point text and row order are preserved.
The normal reproduction entrypoint does not run this migration.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,io,json
ROOT=Path(__file__).resolve().parents[2]
def sha(b):return hashlib.sha256(b).hexdigest()
def table(raw,delimiter='\t'):
    reader=csv.DictReader(io.StringIO(raw.decode()),delimiter=delimiter)
    return reader.fieldnames,list(reader)
def encode(columns,rows,delimiter='\t'):
    stream=io.StringIO(newline='');w=csv.DictWriter(stream,fieldnames=columns,delimiter=delimiter,lineterminator='\n')
    w.writeheader();w.writerows(rows);return stream.getvalue().encode()
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-manifest',type=Path,required=True)
    p.add_argument('--source-root',type=Path,required=True)
    p.add_argument('--cic-interval',type=Path,required=True)
    args=p.parse_args()
    review=json.loads((ROOT/'data/input_review.json').read_text());plans={r['id']:r for r in review['inputs']}
    source=list(csv.DictReader(args.source_manifest.open(),delimiter='\t'))
    assert set(plans)=={r['id'] for r in source}
    output=[];audit=[];column_audit=[];pending=[]
    for old in source:
        relative=Path(old['path'])
        assert relative.parent in [Path('data/processed/files'),Path('data/external/files')],relative
        assert not relative.is_absolute() and '..' not in relative.parts,relative
        plan=plans[old['id']];path=args.source_root/relative;compressed=path.read_bytes()
        assert sha(compressed)==old['compressed_sha256'],old['id']
        raw=gzip.decompress(compressed);assert sha(raw)==old['sha256'],old['id']
        new=dict(old);transform=plan['transform'];new_columns=json.loads(old['columns']);old_count=int(old['rows'] or 0)
        route=plan['route'];count=old_count;removed=[];retained_rows=[]
        if route=='retired':
            new=None;count=0
        elif transform=='CIC_interval_only':
            raw=args.cic_interval.read_bytes();new_columns,data=table(raw);count=len(data)
            assert count==1 and data[0]=={'symbol':'CIC','seqnames':'19','start':'42268537','end':'42295797'}
            new['logical_path']='scripts/plackett-luce/Ordering_Model/reference_files/CIC_hg38_coordinates.tsv'
            new['path']='data/external/files/165_CIC_hg38_coordinates.tsv.gz'
        elif transform!='unchanged':
            delimiter=',' if old['logical_path'].endswith('.csv') else '\t'
            old_columns,data=table(raw,delimiter);assert old_columns==json.loads(old['columns'])
            keep=plan['retain_columns'] or old_columns
            assert len(keep)==len(set(keep)) and set(keep)<=set(old_columns)
            # Retain source column order, including the AC2 expected-table schemas.
            new_columns=[c for c in old_columns if c in keep];removed=[c for c in old_columns if c not in keep]
            seen=set();new_data=[]
            for i,row in enumerate(data):
                if transform=='displayed_RNA_test_only' and not (row['Purity_Source']=='BB_Purity' and row['Model']=='joint' and row['Test']=='Multinomial LRT' and row['Outcome']=='Dominant state'):continue
                if transform=='primary_survival_analyses_only' and row['analysis'] not in ['ASTRO_primary','GBM_primary']:continue
                selected={c:row[c] for c in new_columns}
                if transform=='unique_event_definitions':
                    key=tuple(selected.values())
                    if key in seen:continue
                    seen.add(key)
                assert all(selected[c]==row[c] for c in new_columns)
                new_data.append(selected);retained_rows.append(i+1)
            count=len(new_data);raw=encode(new_columns,new_data,delimiter)
        if new:
            new['tier']='processed' if route=='git' else 'external'
            new['path']=str(Path('data')/new['tier']/'files'/Path(new['path']).name)
            new['redistribution']=plan['classification'];new['rows']=str(count);new['columns']=json.dumps(new_columns)
            new['role']= 'figure_support_aggregate' if route=='git' else 'minimal_author_request_input'
            zipped=compressed if raw==gzip.decompress(compressed) else gzip.compress(raw,mtime=0)
            new.update(sha256=sha(raw),compressed_sha256=sha(zipped),bytes=str(len(raw)),compressed_bytes=str(len(zipped)))
            output.append(new);pending.append((ROOT/new['path'],zipped))
        audit.append({'id':old['id'],'source_logical_path':old['logical_path'],'current_path':new['path'] if new else '',
          'route':route,'classification':plan['classification'],'rationale':plan['rationale'],'consumer':plan['consumer'],
          'transformation':transform,'source_sha256':old['sha256'],'current_sha256':new['sha256'] if new else '',
          'rows_before':old['rows'],'rows_after':count,'columns_before':len(json.loads(old['columns'])),'columns_after':len(new_columns) if new else 0,
          'removed_columns':';'.join(removed),'compressed_bytes_before':old['compressed_bytes'],'compressed_bytes_after':new['compressed_bytes'] if new else 0,
          'retained_values_unchanged':'selected text values and relative row order preserved' if transform not in ['retire','CIC_interval_only'] else 'CIC coordinates checked against source object' if new else 'not used'})
        # Headerless retired files were misdescribed in the initial manifest; do
        # not repeat their first data row (which could contain a TCGA barcode).
        old_columns=json.loads(old['columns'])
        if old['id'] in ['input_166','input_167','input_168']:old_columns=['headerless_field_'+str(i+1) for i in range(len(old_columns))]
        for c in old_columns:
            action='retired_file' if not new else 'retained' if c in new_columns else 'removed'
            column_audit.append({'id':old['id'],'column':c,'action':action,'consumer':plan['consumer'],
              'reason':'required calculation/drawing/verification field or cohort linkage key' if action=='retained' else 'unused by current figures; preserved only in private migration snapshot'})
    assert len(output)==len({x['id'] for x in output})==len({x['path'] for x in output})
    # Nothing is mutated until every original hash and transformation passes.
    for dest,content in pending:
        assert dest.is_relative_to(ROOT/'data')
        dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(content)
    keep={x['path'] for x in output}
    for old in source:
        dest=ROOT/old['path']
        if old['path'] not in keep and dest.is_file():dest.unlink()
    def tsv(path,rows,fields=None):
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    tsv(ROOT/'data/manifest.tsv',output,list(source[0]))
    tsv(ROOT/'data/redistribution_audit.tsv',audit)
    tsv(ROOT/'data/column_audit.tsv',column_audit)
    report={'status':'PASS','source_inputs':len(source),'active_inputs':len(output),'retired_inputs':len(source)-len(output),
       'git_inputs':sum(x['tier']=='processed' for x in output),'request_inputs':sum(x['tier']=='external' for x in output),
       'request_compressed_bytes':sum(int(x['compressed_bytes']) for x in output if x['tier']=='external'),
       'git_compressed_bytes':sum(int(x['compressed_bytes']) for x in output if x['tier']=='processed'),
       'previous_external_compressed_bytes':sum(int(x['compressed_bytes']) for x in source if x['tier']=='external'),
       'columns_removed_from_retained_tables':sum(x['action']=='removed' for x in column_audit),
       'column_entries_removed_with_retired_files':sum(x['action']=='retired_file' for x in column_audit),
       'reduced_tables':sum(x['transformation'] not in ['retire','unchanged'] for x in audit),
       'all_source_inputs_hash_verified':True,'selected_cell_strings_and_relative_row_order_preserved':True}
    (ROOT/'docs/data_minimization_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
