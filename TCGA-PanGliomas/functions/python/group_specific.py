# Adapted from scripts/publishing/figure_standardization/refinement/group_specific.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache
HERE=Path(__file__).parent


EXPECTED={'M1_b':14,'M1_c':20,'M1_d':10}


SCHEMA=['panel','display_index','short_label','full_label','event_type','prevalence_percent','median_plotpos','saved_realizations']


def read(path):return json.loads(Path(path).read_text())


def write(path,data):Path(path).write_text(json.dumps(data,indent=2)+'\n')


def full_metadata(panel):
    from ordering import all_geometry,alias_map,concise
    rows=[]
    for event in all_geometry()[panel]['events']:
        ev=event['ev'];name=event['name']
        rows.append({'panel':panel,'display_index':event['index'],'noevent':int(ev.noevent.iloc[0]),'short_label':alias_map()[concise(name)],'full_label':name,'event_type':ev.cna.iloc[0],'prevalence_percent':float(ev.freq.iloc[0]),'median_plotpos':float(np.median(ev.plotpos)),'saved_realizations':len(ev)})
    return pd.DataFrame(rows)


@lru_cache(maxsize=1)
def selection():
    tables={p:full_metadata(p) for p in EXPECTED}
    folder=ROOT/'results/analysis/pooled_814_rerun_2026-09-28/preparation/pooled_G3_display'
    events={};hashes={}
    for i,(panel,table) in enumerate(tables.items(),1):
        path=folder/f'tcga_all_glioma_DN{i}_mergedseg_G1.txt'
        hashes[str(path.relative_to(ROOT))]=sha(path)
        source=pd.read_csv(path,sep='\t',dtype={'chr':str});lookup={}
        for number,rows in source.groupby('noevent',sort=False):
            assert len(rows[['ID','CNA']].drop_duplicates())==1
            r=rows.iloc[0];typ=r.CNA[1:]
            if typ!='Mut':assert len(rows[['chr','startpos','endpos']].drop_duplicates())==1
            lookup[number]={'source_event':int(number),'type':typ,'source_id':r.ID,'chromosome':str(r.chr),'start':int(r.startpos),'end':int(r.endpos)}
        events[panel]={}
        for r in table.itertuples():
            event={'type':'WGD'} if r.event_type=='WGD' else lookup[r.noevent].copy()
            assert event['type']==r.event_type
            event['label']=r.full_label;events[panel][r.display_index]=event
    kept={};rows=[];matches=[]
    for panel,table in tables.items():
        kept[panel]=[]
        for r in table.itertuples():
            event=events[panel][r.display_index];shared=[]
            if event['type']!='WGD':
                for other,other_events in events.items():
                    if other==panel:continue
                    for index,compare in other_events.items():
                        if event['type']!=compare['type']:continue
                        overlap=0
                        if event['type']=='Mut':same=event['source_id']==compare['source_id']
                        else:
                            overlap=max(0,min(event['end'],compare['end'])-max(event['start'],compare['start'])+1) if event['chromosome']==compare['chromosome'] else 0
                            same=overlap>=1
                        if same:
                            shared.append(f'{other}:{index}')
                            matches.append({'panel':panel,'original_display_index':r.display_index,'label':r.full_label,'event_type':r.event_type,'other_panel':other,'other_display_index':index,'other_label':compare['label'],'overlap_bp':overlap})
            selected=event['type']=='WGD' or not shared
            if selected:kept[panel].append(r)
            rows.append({'panel':panel,'original_display_index':r.display_index,'full_label':r.full_label,'event_type':r.event_type,'prevalence_percent':r.prevalence_percent,'selected':selected,'shared_with':';'.join(shared),**{k:v for k,v in event.items() if k in ['source_event','source_id','chromosome','start','end']}})
        assert len(kept[panel])==EXPECTED[panel]
    for panel in EXPECTED:
        for suffix in ['ordering','violin']:
            path=STAGE/f'derived/{panel}_{suffix}.tsv';hashes[str(path.relative_to(ROOT))]=sha(path)
    return kept,rows,matches,hashes

