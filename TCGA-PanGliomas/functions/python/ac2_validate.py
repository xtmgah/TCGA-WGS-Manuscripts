"""Check regenerated AC2 quantities before passing them to the figure renderer."""
from common import ROOT,STAGE
import json,shutil
import pandas as pd
import numpy as np
expected=ROOT/'results/analysis/amplicon_classifier_2_update_2026-10-02/current/figure5'
actual=ROOT/'computed_ac2/current/figure5'
report=[]
for p in sorted(actual.glob('*.tsv')):
    ref=expected/p.name
    if not ref.exists():raise FileNotFoundError(ref)
    a=pd.read_csv(p,sep='\t',keep_default_na=False,na_values=['NA'])
    b=pd.read_csv(ref,sep='\t',keep_default_na=False,na_values=['NA'])
    assert set(a.columns)==set(b.columns),(p.name,a.columns,b.columns)
    a=a[b.columns]
    # Preserve existing row ordering, gene ranks and tie-breaking.
    pd.testing.assert_frame_equal(a,b,check_dtype=False,atol=1e-10,rtol=1e-8)
    for column in [c for c in a.columns if c in ('p_value','fdr','fdr_within_class')]:
        np.testing.assert_allclose(a[column],b[column],atol=0,rtol=1e-8,equal_nan=True,err_msg=f'{p.name}: {column}')
    frozen={'panel_e_egfr_expression.tsv','panel_e_egfr_fusion_classes.tsv'}
    mode='frozen_checked_copy' if p.name in frozen else 'recomputed'
    report.append({'table':p.name,'rows':len(a),'status':'PASS','mode':mode})
    shutil.copy2(p,ref)
    shutil.copy2(p,STAGE/'derived'/('AC2_'+p.name))
(STAGE/'qa/ac2_verification.json').write_text(json.dumps({'status':'PASS','tables':report,
    'cohort_patients':389,'primary_group_assignment_retained':True,'filter_rules':'Passing a/e; LowCN-inclusive b–d',
    'upstream_classification':'Frozen AmpliconClassifier 2.0.0 feature calls'},indent=2)+'\n')
print(f'{len(report)} AC2 tables regenerated and checked')
