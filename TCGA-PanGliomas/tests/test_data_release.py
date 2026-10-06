"""Failure tests for public export and request-only routing, using synthetic data."""
from pathlib import Path
import csv,gzip,hashlib,json,sys,tempfile,unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'functions/python'))
from release_contract import validate,verify_public_payload
import run

def row(i,tier,columns,values):
    raw=('\t'.join(columns)+'\n'+'\t'.join(values)+'\n').encode();z=gzip.compress(raw,mtime=0)
    return {'id':i,'path':f'data/{tier}/files/{i}.tsv.gz','logical_path':i+'.tsv','tier':tier,
      'sha256':hashlib.sha256(raw).hexdigest(),'compressed_sha256':hashlib.sha256(z).hexdigest(),
      'rows':'1','columns':json.dumps(columns),'redistribution':'git_figure_support_aggregate' if tier=='processed' else 'request_measured_data'},z

class ReleaseContract(unittest.TestCase):
    def setUp(self):
        scratch=ROOT/'outputs/test-data-release';scratch.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=scratch);self.root=Path(self.temp.name)
        self.rows=[];self.plans=[]
        for ident,tier,columns,values in [('public','processed',['group','n'],['A','10']),('private','external',['Subject','value'],['synthetic_subject','5'])]:
            r,z=row(ident,tier,columns,values);self.rows.append(r)
            p=self.root/r['path'];p.parent.mkdir(parents=True);p.write_bytes(z)
            self.plans.append({'id':ident,'route':'git' if tier=='processed' else 'author_request','classification':r['redistribution']})
        self.write()
    def tearDown(self):self.temp.cleanup()
    def write(self):
        with (self.root/'data/manifest.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(self.rows[0]),delimiter='\t');w.writeheader();w.writerows(self.rows)
        (self.root/'data/input_review.json').write_text(json.dumps({'inputs':self.plans}))
    def test_reviewed_split_passes(self):
        result=validate(self.root);self.assertEqual((result['git_inputs'],result['request_inputs']),(1,1))
    def test_identifier_hidden_in_generic_public_field_is_rejected(self):
        # Construct a synthetic barcode without placing an actual study ID in test code.
        with self.assertRaisesRegex(ValueError,'identifier'):
            verify_public_payload(('label\n'+'TCGA'+'-ZZ-'+'0000'+'\n').encode(),['label'])
    def test_individual_column_rejected_even_without_barcodes(self):
        with self.assertRaisesRegex(ValueError,'Individual-record'):
            verify_public_payload(b'Subject\nanonymous\n',['Subject'])
    def test_unreviewed_or_retired_file_in_manifest_rejected(self):
        self.plans[0]['route']='retired';self.write()
        with self.assertRaisesRegex(AssertionError,'Unreviewed'):validate(self.root)
    def test_request_file_cannot_be_promoted_by_manifest_alone(self):
        self.rows[1]['tier']='processed';self.write()
        with self.assertRaises(AssertionError):validate(self.root)
    def test_unlisted_file_under_public_data_rejected(self):
        (self.root/'data/processed/files/stray.tsv.gz').write_bytes(gzip.compress(b'unreviewed'))
        with self.assertRaisesRegex(AssertionError,'Unlisted'):validate(self.root)
    def test_git_only_check_does_not_require_private_files(self):
        (self.root/self.rows[1]['path']).unlink()
        with patch.object(run,'REPO',self.root):
            self.assertEqual(run.check_inputs(self.root/'data/external',public_only=True)['inputs'],1)
            with self.assertRaisesRegex(RuntimeError,'Missing private'):run.check_inputs(self.root/'data/external')
if __name__=='__main__':unittest.main()
