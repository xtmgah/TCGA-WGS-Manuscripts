"""Exercise the real input contract: relocated data, absence and corruption."""
from pathlib import Path
import shutil,sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'functions/python'))
import run

class InputContract(unittest.TestCase):
    def test_relocated_external_bundle_and_corruption(self):
        scratch=ROOT/'.local/tmp';scratch.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as name:
            data=Path(name)/'authorized-inputs'
            shutil.copytree(ROOT/'data/external/files',data/'files')
            self.assertEqual(run.check_inputs(data)['status'],'PASS')
            row=next(r for r in run.manifest() if r['tier']=='external')
            p=run.location(row,data);original=p.read_bytes()
            p.write_bytes(original+b'CORRUPTED')
            with self.assertRaisesRegex(RuntimeError,'Checksum mismatch'):run.check_inputs(data)
            p.write_bytes(original);p.unlink()
            with self.assertRaisesRegex(RuntimeError,'Missing '+row['id']):run.check_inputs(data)

if __name__=='__main__':unittest.main()
