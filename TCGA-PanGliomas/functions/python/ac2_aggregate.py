"""Current AC2 main-figure aggregations from a compact curated feature table."""
from common import ROOT, STAGE
from pathlib import Path
import ast, itertools
import numpy as np
import pandas as pd
OLD=ROOT/'results/analysis/bcor_reclassification_2026-09-21'
EXPECTED=ROOT/'results/analysis/amplicon_classifier_2_update_2026-10-02/current/figure5'
OUT=ROOT/'computed_ac2'
GROUPS=['DN1','DN2','DN3']
CLASSES=['Linear','ecDNA','BFB','Complex-non-cyclic','FAN']
FIELDS=['n_linear','n_ecDNA','n_BFB','n_complex_non_cyclic','n_FAN']
def read(p):return pd.read_csv(p,sep='\t',keep_default_na=False)
def write(df,p):
 p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,sep='\t',index=False,na_rep='NA')
def genes(s):return [] if s in ['', 'NA', '[]', 'Not provided'] else ast.literal_eval(s)
primary=read(OLD/'analysis/GBM_primary_group_roster.tsv')[['Subject','DN_Group']].drop_duplicates()
assert len(primary)==389 and primary.groupby('DN_Group').size().reindex(GROUPS).tolist()==[159,129,101]
cohort=read(EXPECTED.parent/'curated_feature_calls.tsv')
assert not cohort['Feature ID'].duplicated().any()
version='current'
d=OUT/version/'figure5'; d.mkdir(parents=True,exist_ok=True)
df=cohort.merge(primary,on='Subject',how='inner',validate='many_to_one')
passing=df[df.filter_passing]
counts=passing.groupby(['Subject','Classification'])['Feature ID'].nunique().unstack(fill_value=0).reindex(columns=CLASSES,fill_value=0)
counts.columns=FIELDS
subjects=primary.merge(counts,on='Subject',how='left').fillna(0)
subjects[FIELDS]=subjects[FIELDS].astype(int)
write(subjects,d/'amplicon_subject_table.tsv')
long=subjects.melt(id_vars=['Subject','DN_Group'],value_vars=FIELDS,var_name='field',value_name='n_amplicons')
long['amplicon_class']=long.field.map(dict(zip(FIELDS,CLASSES)))
write(long,d/'panel_a_amplicon_patient_counts.tsv')
rates=long.groupby(['DN_Group','amplicon_class'],sort=False).agg(positive=('n_amplicons',lambda s:sum(s>0)),total=('Subject','size'),n_features=('n_amplicons','sum')).reset_index()
rates['pct']=100*rates.positive/rates.total
write(rates,d/'panel_a_amplicon_rates.tsv')
ec=df[df.Classification.eq('ecDNA')].copy()
write(ec,d/'ecdna_candidate_features.tsv')
patients=ec[['Subject','DN_Group']].drop_duplicates()
den=patients.groupby('DN_Group').size().rename('n_ecdna_subjects').reset_index()
write(den,d/'ecdna_positive_subject_denominators.tsv')
sg=ec.assign(gene=ec.Oncogenes.map(genes)).explode('gene').dropna(subset='gene').rename(columns={'Feature ID':'Feature_ID'})
sg=sg[['Subject','Tumor_Barcode','DN_Group','Feature_ID','context','max_cn','gene']].drop_duplicates()
spg=sg.sort_values(['max_cn','Feature_ID'],ascending=[False,True]).drop_duplicates(['Subject','gene'])
write(spg,d/'ecdna_subject_gene_context.tsv')
burden=patients.merge(spg.groupby('Subject').gene.nunique().rename('n_distinct_oncogenes'),on='Subject',how='left').fillna(0)
burden['n_distinct_oncogenes']=burden.n_distinct_oncogenes.astype(int)
write(burden,d/'panel_b_ecdna_oncogene_burden.tsv')
gc=spg.groupby('gene').Subject.nunique().reset_index(name='n').sort_values(['n','gene'],ascending=[False,True])
top=gc.head(12).gene.tolist(); tests=[]; prevalence=[]
for gene in gc.loc[gc.n.ge(5),'gene']:
    pos=set(spg.loc[spg.gene.eq(gene),'Subject'])
    for g in GROUPS:
        p=set(patients.loc[patients.DN_Group.eq(g),'Subject']); o=set(patients.Subject)-p
        a,c=len(pos&p),len(pos&o); b,dd=len(p)-a,len(o)-c
        aa,bb,cc,d2=np.array([a,b,c,dd])+.5
        lor=np.log(aa*d2/(bb*cc)); se=np.sqrt(1/aa+1/bb+1/cc+1/d2)
        tests.append(dict(gene=gene,DN_Group=g,positive_target=a,total_target=len(p),positive_other=c,total_other=len(o),odds_ratio=np.exp(lor),log2_or=lor/np.log(2),log2_or_low=(lor-1.96*se)/np.log(2),log2_or_high=(lor+1.96*se)/np.log(2)))
for gene,g in itertools.product(top,GROUPS):
    n=int(den.loc[den.DN_Group.eq(g),'n_ecdna_subjects'].iloc[0]); p=sum(spg.gene.eq(gene)&spg.DN_Group.eq(g))
    prevalence.append(dict(DN_Group=g,gene=gene,positive=p,n_ecdna_subjects=n,pct=100*p/n))
write(pd.DataFrame(tests),d/'panel_c_oncogene_one_vs_rest_fisher.tsv')
write(pd.DataFrame(prevalence),d/'panel_c_oncogene_prevalence.tsv')
sden=ec.groupby('DN_Group').Tumor_Barcode.nunique().rename('n_ecdna_tumors').reset_index()
write(sden,d/'panel_d_cooccurrence_tumor_denominators.tsv')
nodes=[]; bars=[]; pairs=[]
for g in GROUPS:
    dat=sg[sg.DN_Group.eq(g)]; n=int(sden.loc[sden.DN_Group.eq(g),'n_ecdna_tumors'].iloc[0])
    nn=dat.groupby('gene').agg(n_tumors=('Tumor_Barcode','nunique'),unique_feature_ids=('Feature_ID','nunique')).reset_index().sort_values(['n_tumors','unique_feature_ids','gene'],ascending=[False,False,True]).head(12)
    nn['rank']=range(1,len(nn)+1); nn['DN_Group']=g; nn['n_ecdna_tumors']=n; nn['pct']=100*nn.n_tumors/n; nn['x']=nn['rank']; nodes.append(nn)
    rank=dict(zip(nn.gene,nn['rank']))
    for (gene,context), cc in dat[dat.gene.isin(rank)].groupby(['gene','context']):
        bars.append(dict(DN_Group=g,gene=gene,context=context,positive=cc.Tumor_Barcode.nunique(),rank=rank[gene],n_ecdna_tumors=n,pct=100*cc.Tumor_Barcode.nunique()/n))
    pair_counts={}
    for _, dd in dat[dat.gene.isin(rank)].groupby('Tumor_Barcode'):
        for pair in itertools.combinations(sorted(set(dd.gene)),2):pair_counts[pair]=pair_counts.get(pair,0)+1
    for (a,b),n_pair in pair_counts.items():
        if n_pair>=3:pairs.append(dict(DN_Group=g,gene_1=a,gene_2=b,pair_count=n_pair,x=min(rank[a],rank[b]),xend=max(rank[a],rank[b])))
nodes=pd.concat(nodes,ignore_index=True)
write(nodes,d/'panel_d_cooccurrence_nodes.tsv');write(nodes,d/'panel_d_cooccurrence_top_genes.tsv')
write(pd.DataFrame(bars),d/'panel_d_cooccurrence_context_prevalence.tsv');write(pd.DataFrame(pairs),d/'panel_d_cooccurrence_pairs.tsv')
egfr=ec[ec['All genes'].map(lambda x:'EGFR' in genes(x))|ec.Oncogenes.map(lambda x:'EGFR' in genes(x))]
# The original EGFR feature table was filter-passing, unlike panels b-d.
# Preserve that analysis and export a candidate-inclusive sensitivity.
cn_all=egfr.groupby(['Subject','DN_Group']).max_cn.max().rename('max_egfr_ecDNA_cn').reset_index()
write(cn_all,d/'panel_e_egfr_copy_number_candidate_sensitivity.tsv')
cn=egfr[egfr.filter_passing].groupby(['Subject','DN_Group']).max_cn.max().rename('max_egfr_ecDNA_cn').reset_index()
write(cn,d/'panel_e_egfr_ecdna_copy_number.tsv')
for stem in ['panel_e_egfr_expression','panel_e_egfr_fusion_classes','panel_e_pairwise_tests','panel_e_global_tests']:
    write(read(EXPECTED/f'{stem}.tsv'),d/f'{stem}.tsv')
audit=primary.copy()
audit['filter_passing_ecDNA']=audit.Subject.isin(set(passing.loc[passing.Classification.eq('ecDNA'),'Subject']))
audit['candidate_ecDNA']=audit.Subject.isin(set(ec.Subject))
audit['candidate_only']=audit.candidate_ecDNA&~audit.filter_passing_ecDNA
write(audit,OUT/version/'ecDNA_patient_filter_audit.tsv')
