# Adapted from scripts/publishing/figure_standardization/refinement/render_figure5.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache
setup_matplotlib()
import matplotlib.pyplot as plt
from matplotlib import colormaps,colors,ticker
from matplotlib.patches import Rectangle,PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.lines import Line2D


B = ROOT / 'results/analysis/bcor_reclassification_2026-09-21'


D = B / 'figures/figure5/derived-data'


GROUPS = ['DN1', 'DN2', 'DN3']


NAMES = dict(zip(GROUPS, ['C19', 'C19/20', 'TWR']))


DN = dict(zip(GROUPS, ['#E37E78', '#8C231D', '#E6A0C4']))


STATE_ORDER = ['AC-like', 'NPC-like', 'OPC-like', 'MES-like']


STATE_COLORS = dict(zip(STATE_ORDER, ['#EF8500', '#9143A9', '#3C81BC', '#B05724']))


CONTEXT = {
    'Heavily rearranged multichromosomal': '#8B0000',
    'Heavily rearranged unichromosomal': '#D00000',
    'Simple circular complex background': '#F47C55',
    'Simple circular simple background': '#4169D8',
    'Unknown': '#8A8A8A',
}


FUSION = {'All others': '#00A08A', 'EGFRvIII': '#6FA5B9', 'None': '#D5D5D5'}


AMP = ['#46085C', '#3C4F8A', '#24868E', '#3ABA76']


SOURCES = {}


CHECKS = {}


FONT_CHECKS = {}


def read(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    # “None” is an explicit fusion class, never a missing value.
    return pd.read_csv(path, sep='\t', keep_default_na=False, na_values=['NA'], float_precision='round_trip')


def fig_new(letter, width, height, title):
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))
    fig._mm = (width, height)
    fig._letter = letter
    txt(fig, .2, .60, title, 8)
    return fig


def txt(fig, x, y, value, size=6.5, ha='left', va='top', color='#222222', weight='normal', **kw):
    w, h = fig._mm
    return fig.text(x / w, 1 - y / h, value, fontsize=size, ha=ha, va=va,
                    color=color, fontweight=weight, **kw)


def axes(fig, x, top, width, height):
    w, h = fig._mm
    ax = fig.add_axes([x / w, 1 - (top + height) / h, width / w, height / h])
    ax.tick_params(length=2, width=.6, pad=1.2, labelsize=6.5)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color('#555555')
        spine.set_linewidth(.6)
    return ax


def canvas(fig):
    w, h = fig._mm
    ax = axes(fig, 0, 0, w, h)
    ax.set(xlim=(0, w), ylim=(h, 0))
    ax.axis('off')
    return ax


def swatch(ax, x, y, color, width=1.8, height=1.8):
    ax.add_patch(Rectangle((x, y-height/2), width, height, facecolor=color,
                           edgecolor='none', clip_on=False))


def tidy(ax, grid='y'):
    ax.grid(axis=grid, color='#E6E6E6', linewidth=.45)
    ax.spines[['top', 'right']].set_visible(False)


def fmt_q(value):
    if value < .0001:
        return 'q < 0.0001'
    if value < .001:
        return 'q < 0.001'
    return f'q = {value:.3f}'.rstrip('0').rstrip('.')


def bracket(ax, x1, x2, y, label, dy=.02, text_dy=.035, log=False):
    ylow = y / 1.04 if log else y - dy
    ax.plot([x1, x1, x2, x2], [ylow, y, y, ylow], color='#666666', lw=.45)
    ax.text((x1+x2)/2, y * 1.06 if log else y + text_dy, label,
            fontsize=6, ha='center', va='bottom', color='#555555')


def swarm_offsets(values, log=False, width=.19):
    """Deterministic displacement only; all observations and y values retained."""
    v = np.log10(values) if log else np.asarray(values, dtype=float)
    order = np.argsort(v, kind='stable')
    offsets = np.zeros(len(v))
    extent = max(float(np.ptp(v)), 1)
    bins = np.floor((v-min(v))/extent*30).astype(int)
    density = np.bincount(bins)
    for rank,index in enumerate(order):
        # Deterministic low-discrepancy jitter, tapered by local display density.
        displacement = 2*((rank*.6180339887498949+.5)%1)-1
        offsets[index] = displacement*width*np.sqrt(density[bins[index]]/max(density))
    return offsets


def boxes(ax, data, field, log=False, cap=None, size=5):
    counts = {}
    for index, group in enumerate(GROUPS):
        vals = data.loc[data.DN_Group.eq(group), field].dropna().to_numpy(float)
        counts[group] = len(vals)
        displayed = np.minimum(vals, cap) if cap is not None else vals
        ax.boxplot([displayed], positions=[index], widths=.50, patch_artist=True,
            showfliers=False, manage_ticks=False,
            boxprops={'facecolor':'white', 'edgecolor':'#333333', 'linewidth':.6},
            whiskerprops={'color':'#333333', 'linewidth':.6},
            capprops={'color':'#333333', 'linewidth':.6},
            medianprops={'color':'#222222', 'linewidth':.9})
        ax.scatter(index + swarm_offsets(displayed, log=log), displayed, s=size,
            facecolor=DN[group], edgecolor='#333333', linewidth=.17,
            alpha=.75, zorder=3)
    ax.set_xlim(-.48, 2.48)
    ax.set_xticks(range(3), [f'{NAMES[g]}\nn = {counts[g]}' for g in GROUPS])
    for tick, group in zip(ax.get_xticklabels(), GROUPS):
        tick.set_color('#000000'); tick.set_fontsize(6.5)
    if log:
        ax.set_yscale('log')
        ax.yaxis.set_minor_locator(ticker.NullLocator())
    tidy(ax)
    return counts


def save(fig):
    """Record native font sizes and text bounds for inspection at final size."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    from matplotlib.text import Text
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    minimum = min(t.get_fontsize() for t in texts)
    assert minimum >= 6, (fig._letter, minimum)
    outside = []
    fw, fh = fig.bbox.width, fig.bbox.height
    for t in texts:
        box = t.get_window_extent(renderer)
        if box.x0 < -.75 or box.y0 < -.75 or box.x1 > fw+.75 or box.y1 > fh+.75:
            outside.append({'text':t.get_text(), 'bbox_px':list(box.bounds)})
    FONT_CHECKS[fig._letter] = {'minimum_pt':minimum, 'text_objects':len(texts), 'outside_panel':outside}
    path=STAGE / 'panels' / f'M5_{fig._letter}.pdf'
    fig.savefig(path)
    with fitz.open(path) as doc:
        page=doc[0]
        native_out=[]
        native_sizes=[]
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines',[]):
                for span in line['spans']:
                    native_sizes.append(span['size'])
                    x0,y0,x1,y1=span['bbox']
                    if x0<-.2 or y0<-.2 or x1>page.rect.width+.2 or y1>page.rect.height+.2:
                        native_out.append({'text':span['text'],'bbox_pt':list(span['bbox'])})
        FONT_CHECKS[fig._letter].update(native_pdf_outside_panel=native_out,
            native_pdf_minimum_pt=min(native_sizes),embedded_images=len(page.get_images()))
    plt.close(fig)


def panel_a():
    df = read(B/'ecdna/amplicons/gbm_dn3_amplicon_subject_table.tsv')
    rates = read(D/'panel_a_amplicon_rates.tsv')
    fig = fig_new('a', 56, 45.5, 'Amplicon architecture')
    classes = [('Linear','n_linear'), ('ecDNA','n_ecDNA'), ('Complex non-cyclic','n_complex_non_cyclic')]
    drawn = 0
    all_counts = {}
    for k, group in enumerate(GROUPS):
        sub = df.loc[df.DN_Group.eq(group)].copy()
        n = len(sub); all_counts[group] = n
        top = 6.4 + k*10.3
        txt(fig, 18, top-2.5, f'{NAMES[group]}  n = {n}', 6.5, color=DN[group])
        ax = axes(fig, 18, top, 35, 7.2)
        ax.set(xlim=(0,84), ylim=(2.6,-.6))
        for j, (label, field) in enumerate(classes):
            vals = sub.loc[sub[field].gt(0)].sort_values([field, 'Subject'])
            for i, row in enumerate(vals.itertuples(index=False)):
                value = int(getattr(row, field))
                assert 1 <= value <= 4
                ax.add_patch(Rectangle((i*100/n, j-.35), 100/n, .70,
                             facecolor=AMP[value-1], edgecolor=(1,1,1,.25), lw=.1))
                drawn += 1
            pct = rates.loc[rates.DN_Group.eq(group)&rates.amplicon_class.eq(label),'pct'].iloc[0]
            assert abs(pct-len(vals)*100/n) < 1e-8
            ax.text(pct+1.3, j, f'{pct:.1f}%', fontsize=6, va='center')
        ax.set_yticks(range(3), [p[0] for p in classes], fontsize=6)
        ax.set_xticks([0,25,50,75], ['0%','25%','50%','75%'], fontsize=6)
        ax.tick_params(axis='y', length=0, pad=1.3)
        ax.spines[['top','right','left']].set_visible(False)
        ax.grid(axis='x', color='#E6E6E6', linewidth=.45)
        if k != 2:
            ax.tick_params(axis='x', bottom=False, labelbottom=False)
            ax.spines['bottom'].set_visible(False)
    txt(fig, 35, 37.8, 'Patients (%)', 6.5, ha='center')
    ov = canvas(fig)
    txt(fig, .3, 42.2, 'Amplicons / patient', 6)
    for i, cl in enumerate(AMP):
        x = 28.5+i*6.4
        swatch(ov, x, 43, cl, 1.7,1.7)
        txt(fig, x+2.4, 43, str(i+1), 6, va='center')
    CHECKS['a'] = {'patients_by_group':all_counts, 'individual_tiles':drawn,
        'rates':rates.to_dict('records'), 'amplicon_colors':AMP}
    save(fig)


def panel_b():
    df = read(D/'panel_b_ecdna_oncogene_burden.tsv')
    tests = read(D/'panel_b_ecdna_oncogene_burden_pairwise_tests.tsv')
    fig = fig_new('b',48,45.5,'Oncogene burden')
    ax = axes(fig,8,5.8,39,32.5)
    counts = boxes(ax,df,'n_distinct_oncogenes',cap=10,size=4.5)
    ax.set_ylim(-.2,11.9)
    ax.set_yticks([0,2,4,6,8,10], ['0','2','4','6','8','≥10'])
    ax.set_ylabel('Distinct oncogenes', labelpad=1.5)
    significant = tests.loc[tests.fdr.lt(.05)]
    for i,row in enumerate(significant.itertuples()):
        bracket(ax,GROUPS.index(row.group_1),GROUPS.index(row.group_2),10.45+i*.9,fmt_q(row.fdr),dy=.1,text_dy=.04)
    CHECKS['b'] = {'patients_by_group':counts,'observations':len(df),
        'display_cap':10,'values_at_or_above_cap':int(df.n_distinct_oncogenes.ge(10).sum()),
        'displayed_tests':significant.to_dict('records')}
    save(fig)


def panel_c():
    df = read(D/'panel_c_oncogene_prevalence.tsv')
    tests = read(D/'panel_c_oncogene_one_vs_rest_fisher.tsv')
    genes = df.gene.drop_duplicates().tolist()
    assert len(genes)==12
    table = df.pivot(index='gene',columns='DN_Group',values='pct').loc[genes,GROUPS]
    fig = fig_new('c',60,45.5,'ecDNA oncogene prevalence')
    ax = axes(fig,9.3,5.3,50,26.2)
    cmap = colors.LinearSegmentedColormap.from_list('prevalence',['#F5F7F7','#3B7488'])
    norm = colors.Normalize(0,100)
    ax.pcolormesh(np.arange(4),np.arange(13),table.values,cmap=cmap,norm=norm,
                  edgecolors='white',linewidth=.5,rasterized=False)
    ax.set(ylim=(12,0), xlim=(0,3))
    ax.set_yticks(np.arange(12)+.5, genes, fontsize=6.5)
    ns = df.groupby('DN_Group').n_ecdna_subjects.first()
    ax.set_xticks(np.arange(3)+.5,[f'{NAMES[g]}\nn = {ns[g]}' for g in GROUPS],fontsize=6.5)
    for t,g in zip(ax.get_xticklabels(),GROUPS):t.set_color('#000000')
    ax.tick_params(length=0,pad=1.8)
    for spine in ax.spines.values():spine.set_visible(False)
    marks=[]
    for iy,gene in enumerate(genes):
        for ix,group in enumerate(GROUPS):
            value=table.loc[gene,group]
            ax.text(ix+.40,iy+.5,f'{value:.0f}%',fontsize=6.5,ha='center',va='center')
            t=tests.loc[tests.gene.eq(gene)&tests.DN_Group.eq(group)].iloc[0]
            if t.fdr<.05:
                direction='Higher' if t.log2_or>=0 else 'Lower'
                ax.scatter(ix+.86,iy+.5,marker='^' if direction=='Higher' else 'v',s=11,
                           facecolor='white',edgecolor='#222222',linewidth=.5)
                marks.append({'gene':gene,'group':group,'direction':direction,'q':float(t.fdr)})
    ca=axes(fig,11,40.0,20,1.6)
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=ca,orientation='horizontal',ticks=[0,50,100])
    cb.solids.set_rasterized(False)
    cb.ax.tick_params(labelsize=6,length=1.5,pad=.4)
    cb.outline.set_visible(False)
    txt(fig,10.5,37.5,'Patients (%)',6)
    ov=canvas(fig)
    txt(fig,34.5,38.1,'q < 0.05 vs rest',6)
    ov.scatter([35.5],[42],marker='^',s=11,facecolor='white',edgecolor='#222222',lw=.5)
    ov.scatter([48.6],[42],marker='v',s=11,facecolor='white',edgecolor='#222222',lw=.5)
    txt(fig,37.2,42,'Higher',6,va='center');txt(fig,50.3,42,'Lower',6,va='center')
    CHECKS['c']={'genes':genes,'cells':int(table.size),'patient_denominators':ns.to_dict(),
                 'all_percentages':table.to_dict(),'significance_marks':marks,'colorbar_range':[0,100]}
    save(fig)


def panel_d():
    nodes=read(D/'panel_d_cooccurrence_nodes.tsv')
    pairs=read(D/'panel_d_cooccurrence_pairs.tsv')
    bars=read(D/'panel_d_cooccurrence_context_prevalence.tsv')
    den=read(D/'panel_d_cooccurrence_tumor_denominators.tsv')
    fig=fig_new('d',172,43.5,'ecDNA oncogene co-occurrence')
    pair_max=int(pairs.pair_count.max()); feature_max=int(nodes.unique_feature_ids.max())
    norm=colors.Normalize(3,pair_max); cmap=colormaps['plasma']
    def node_dia(n):return 1.0+2.0*np.sqrt((n-1)/(feature_max-1))
    def edge_lw(n):return .35+(n-3)/(pair_max-3)*1.2
    block=(172-8)/3
    for k,g in enumerate(GROUPS):
        left=k*(block+4)
        nn=nodes.loc[nodes.DN_Group.eq(g)].sort_values('rank')
        bb=bars.loc[bars.DN_Group.eq(g)]
        ee=pairs.loc[pairs.DN_Group.eq(g)]
        n=int(den.loc[den.DN_Group.eq(g),'n_ecdna_tumors'].iloc[0])
        txt(fig,left+5,3.8,f'{NAMES[g]}  (n = {n} specimens)',6.5,color=DN[g],weight='bold')
        ax=axes(fig,left+6.0,7.3,block-6.5,10.2)
        base=np.zeros(len(nn))
        # The original stack puts the final context at the bottom.
        for context in reversed(CONTEXT):
            val=bb.loc[bb.context.eq(context)].set_index('gene').pct.reindex(nn.gene,fill_value=0).to_numpy()
            ax.bar(nn['rank'],val,bottom=base,color=CONTEXT[context],width=.8,
                   edgecolor='#333333',lw=.15)
            base+=val
        assert np.allclose(base,nn.pct.to_numpy())
        ax.set(xlim=(.5,12.5),ylim=(0,100))
        ax.set_xticks(nn['rank'],nn.gene,rotation=90,fontsize=6)
        ax.tick_params(axis='x',length=1.5,pad=1.0)
        ax.set_yticks([0,25,50,75,100],['0%','25%','50%','75%','100%'],fontsize=6)
        if k:ax.tick_params(axis='y',labelleft=False,length=0)
        tidy(ax)
        ar=axes(fig,left+6.0,24.6,block-6.5,7.2)
        ar.set(xlim=(.5,12.5),ylim=(-1.08,.16));ar.axis('off')
        for row in ee.sort_values('pair_count').itertuples():
            dx=row.xend-row.x;depth=.24+dx*.069
            path=MplPath([(row.x,0),(row.x+dx*.22,-depth),(row.xend-dx*.22,-depth),(row.xend,0)],
                         [MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4])
            ar.add_patch(PathPatch(path,fill=False,color=cmap(norm(row.pair_count)),lw=edge_lw(row.pair_count),alpha=.68))
        ar.scatter(nn['rank'],np.zeros(len(nn)),s=(node_dia(nn.unique_feature_ids.to_numpy())*MM)**2,
                   facecolor='white',edgecolor='#222222',lw=.5,zorder=3,clip_on=False)
    ov=canvas(fig)
    entries=list(CONTEXT.items())
    locs=[(0,34),(54,34),(110,34),(0,37.2),(54,37.2)]
    for (label,color),(x,y) in zip(entries,locs):
        swatch(ov,x,y,color,1.8,1.8);txt(fig,x+2.7,y,label,6,va='center')
    txt(fig,0,41.1,'Pair count',6,va='center')
    for x,n in zip([15,28,41],[5,10,15]):
        ov.plot([x,x+5],[41.1,41.1],color=cmap(norm(n)),lw=edge_lw(n))
        txt(fig,x+6.2,41.1,str(n),6,va='center')
    txt(fig,64,41.1,'Unique ecDNA feature IDs',6,va='center')
    for x,n in zip([98,116,134,152],[20,40,60,80]):
        ov.scatter([x],[41.1],s=(node_dia(n)*MM)**2,facecolor='white',edgecolor='#222222',lw=.5)
        txt(fig,x+2.8,41.1,str(n),6,va='center')
    CHECKS['d']={'specimen_denominators':den.to_dict('records'),'nodes':nodes.to_dict('records'),
        'all_edges':pairs.to_dict('records'),'context_values':bars.to_dict('records'),
        'edges_minimum_count':3,'pair_count_range':[3,pair_max],'feature_id_range':[1,feature_max],
        'encodings':'arc color and width = same-specimen pair count; node size = unique feature IDs'}
    save(fig)


def panel_e():
    cn=read(D/'panel_e_egfr_ecdna_copy_number.tsv')
    ex=read(D/'panel_e_egfr_expression.tsv')
    fu=read(D/'panel_e_egfr_fusion_classes.tsv')
    tests=read(D/'panel_e_pairwise_tests.tsv')
    fig=fig_new('e',172,34.5,'EGFR specialization')
    widths=(172-8)/3
    for k,title in enumerate(['EGFR ecDNA copy number','EGFR expression','EGFR fusion classes']):
        txt(fig,k*(widths+4)+widths/2,4.1,title,7.5,ha='center')
    ax=axes(fig,9.5,9.5,44,17.3)
    ncn=boxes(ax,cn,'max_egfr_ecDNA_cn',log=True,size=5.5)
    ax.set_ylim(5,2600);ax.set_yticks([10,100,1000],['10','100','1,000'])
    ax.set_ylabel('Maximum copy number',labelpad=1)
    ctests=tests.loc[tests.metric.eq('EGFR ecDNA maximum copy number')&tests.fdr.lt(.05)]
    for row,y in zip(ctests.itertuples(),[520,1150]):
        bracket(ax,GROUPS.index(row.group_1),GROUPS.index(row.group_2),y,fmt_q(row.fdr),log=True)
    ax2=axes(fig,widths+4+7,9.5,46.5,17.3)
    nex=boxes(ax2,ex,'EGFR_expression',size=5.5)
    ax2.set_ylim(7.8,29);ax2.set_yticks([10,15,20,25]);ax2.set_ylabel('VST',labelpad=1)
    etests=tests.loc[tests.metric.eq('EGFR expression')&tests.fdr.lt(.05)]
    for row,y in zip(etests.itertuples(),[21.2,24.2,27.2]):
        bracket(ax2,GROUPS.index(row.group_1),GROUPS.index(row.group_2),y,fmt_q(row.fdr),dy=.16,text_dy=.05)
    ax3=axes(fig,2*(widths+4)+7,9.5,46.5,17.3)
    base=np.zeros(3)
    for name,color in FUSION.items():
        vals=fu.loc[fu.Architecture.eq(name)].set_index('DN_Group').pct_subjects.reindex(GROUPS).to_numpy()
        ax3.bar(range(3),vals,bottom=base,width=.67,color=color,edgecolor='white',lw=.45)
        for x,(bot,val) in enumerate(zip(base,vals)):
            if val>=10:
                ax3.text(x,bot+val/2,f'{val:.0f}%',fontsize=6,ha='center',va='center')
            else:
                ax3.annotate(f'{val:.0f}%',xy=(x-.335,bot+val/2),
                    xytext=(x-.53,bot+val/2),fontsize=6,ha='center',va='center',
                    arrowprops={'arrowstyle':'-','lw':.45,'color':'#444444',
                                'shrinkA':.5,'shrinkB':0})
        base+=vals
    assert np.allclose(base,100)
    nf=fu.groupby('DN_Group').total_subjects.first().to_dict()
    ax3.set(ylim=(0,100),xlim=(-.5,2.5))
    ax3.set_yticks([0,25,50,75,100],['0%','25%','50%','75%','100%'],fontsize=6)
    ax3.set_xticks(range(3),[f'{NAMES[g]}\nn = {nf[g]}' for g in GROUPS],fontsize=6.5)
    for t,g in zip(ax3.get_xticklabels(),GROUPS):t.set_color('#000000')
    ax3.set_ylabel('Patients',labelpad=1);tidy(ax3)
    ov=canvas(fig)
    for (name,color),x in zip(FUSION.items(),[2*(widths+4)+7,2*(widths+4)+24,2*(widths+4)+41]):
        swatch(ov,x,8.2,color,1.5,1.5);txt(fig,x+2,8.2,name,6,va='center')
    CHECKS['e']={'copy_number_n':ncn,'expression_n':nex,'fusion_n':nf,
        'copy_number_scale':'log10','copy_number_values':cn[['Subject','DN_Group','max_egfr_ecDNA_cn']].to_dict('records'),
        'expression_values':ex[['Subject','DN_Group','EGFR_expression']].to_dict('records'),
        'fusion_counts':fu.to_dict('records'),'displayed_tests':pd.concat([ctests,etests]).to_dict('records')}
    save(fig)


def panels_f_g():
    data=read(B/'rna/analysis/GBM_state_analysis_metadata.tsv')
    composition=read(B/'rna/analysis/GBM_state_composition.tsv')
    tests=read(B/'rna/analysis/GBM_state_global.tsv')
    row_order=['OPC-like','NPC-like','AC-like','MES-like']
    assert len(data)==205
    assert data.groupby('DN_Group').size().reindex(GROUPS).tolist()==[85,66,54]
    data['dominant_score']=data[row_order].max(axis=1)
    data['state_order']=pd.Categorical(data.dominant_state,categories=row_order,ordered=True)
    data=data.sort_values(['DN_Group','state_order','dominant_score'],ascending=[True,True,False])
    z=((data[row_order]-data[row_order].mean())/data[row_order].std()).clip(-2,2)
    selected=tests.loc[tests.Purity_Source.eq('BB_Purity')&tests.Model.eq('joint')&tests.Test.eq('Multinomial LRT')&tests.Outcome.eq('Dominant state')]
    assert len(selected)==1 and int(selected.N.iloc[0])==205
    p=float(selected.P.iloc[0])
    fig=fig_new('f',96,33.5,'Neftel states')
    ax=axes(fig,12.5,7.7,82.5,17.5)
    cmap=colors.LinearSegmentedColormap.from_list('current_state_z',['#4A73B6','#FFFFFF','#D84644'])
    ax.pcolormesh(np.arange(206),np.arange(5),z.to_numpy().T,cmap=cmap,vmin=-2,vmax=2,rasterized=False)
    ax.set(xlim=(0,205),ylim=(4,0))
    ax.set_yticks(np.arange(4)+.5,row_order,fontsize=6.5)
    ax.set_xticks([]);ax.tick_params(length=0,pad=2)
    for sp in ax.spines.values():sp.set_visible(False)
    n=0
    for g in GROUPS:
        count=int(data.DN_Group.eq(g).sum())
        x=12.5+(n+count/2)/205*82.5
        txt(fig,x,4.8,f'{NAMES[g]}  n = {count}',6.5,ha='center',color=DN[g])
        if n:ax.axvline(n,color='white',lw=.7)
        n+=count
    cax=axes(fig,47,28.3,22,1.4)
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(-2,2),cmap=cmap),cax=cax,orientation='horizontal',ticks=[-2,0,2])
    cb.solids.set_rasterized(False)
    cb.ax.tick_params(labelsize=6,length=1.5,pad=.4);cb.outline.set_visible(False)
    txt(fig,25,28,'GSVA z score',6.5)
    CHECKS['f']={'patients':205,'group_n':[85,66,54],'sample_order':data.sample_id.tolist(),
       'row_order':row_order,'z_score_values':z.to_numpy().tolist(),'z_score_display_limits':[-2,2],
       'normalization':'same fixed GSVA scores, cohort means and sample SD as corrected Figure5fg.py'}
    save(fig)

    fig=fig_new('g',72,33.5,'Dominant RNA state')
    txt(fig,.2,4.3,f'Adjusted P = {p:.5f}',6.5)
    ax=axes(fig,16,8.8,52.5,15.0)
    counts=data.groupby('DN_Group').size().reindex(GROUPS)
    left=np.zeros(3)
    recorded=[]
    for st in STATE_ORDER:
        nn=composition.loc[composition.dominant_state.eq(st)].set_index('DN_Group').N.reindex(GROUPS).to_numpy()
        independently=np.array([sum((data.DN_Group==g)&(data.dominant_state==st)) for g in GROUPS])
        assert np.array_equal(nn,independently)
        vals=100*nn/counts.to_numpy()
        ax.barh(range(3),vals,left=left,height=.72,color=STATE_COLORS[st],edgecolor='white',lw=.45)
        for y,(start,value,n) in enumerate(zip(left,vals,nn)):
            ax.text(start+value/2,y,f'{value:.0f}%',fontsize=6.5,ha='center',va='center',
                    color='#222222' if st=='AC-like' else 'white')
            recorded.append({'group':GROUPS[y],'state':st,'n':int(n),'pct':float(value)})
        left+=vals
    assert np.allclose(left,100)
    ax.set(xlim=(0,100),ylim=(2.6,-.6))
    ax.set_yticks(range(3),[f'{NAMES[g]} (n = {counts[g]})' for g in GROUPS],fontsize=6.5)
    ax.set_xticks([0,50,100],['0%','50%','100%'],fontsize=6)
    ax.tick_params(axis='y',length=0,pad=2);ax.spines[['top','right','left']].set_visible(False)
    ax.grid(axis='x',color='#E8E8E8',lw=.45)
    txt(fig,42,27.2,'Patients (%)',6.5,ha='center')
    ov=canvas(fig)
    for st,x in zip(STATE_ORDER,[9,24,40,56]):
        swatch(ov,x,31.1,STATE_COLORS[st],1.7,1.7);txt(fig,x+2.4,31.1,st,6,va='center')
    CHECKS['g']={'patients':205,'counts':recorded,'recomputed_adjusted_test':selected.to_dict('records'),
       'state_order':STATE_ORDER,'state_colors':STATE_COLORS,'orientation':'horizontal 100% stacked'}
    save(fig)


def panel_h():
    data=read(B/'rna/de/pairwise_diverging_de_counts_data.tsv')
    fig=fig_new('h',80,50.5,'Differential expression')
    ov=canvas(fig)
    swatch(ov,16,5.9,colors.to_rgba('#666666',.3),2,1.8)
    txt(fig,19,5.9,'Strong DE',6,va='center')
    swatch(ov,37,5.9,colors.to_rgba('#666666',.95),2,1.8)
    txt(fig,40,5.9,'Replicated subset',6,va='center')
    ax=axes(fig,16,10.1,62.5,31.7)
    contrasts=['DN2 vs DN1','DN3 vs DN1','DN3 vs DN2']
    for iy,contrast in enumerate(contrasts):
        for row in data.loc[data.Contrast.eq(contrast)].itertuples():
            group=row.Higher_Expression_Group
            primary=row.Signed_Strong_DE_Genes;rep=row.Signed_Center_Replicated_Strong_DE_Genes
            ax.barh(iy,primary,height=.62,color=colors.to_rgba(DN[group],.3),edgecolor='#666666',lw=.4)
            ax.barh(iy,rep,height=.25,color=DN[group],edgecolor='none')
            inside=abs(primary)>=800
            x=primary-np.sign(primary)*45 if inside else primary+np.sign(primary)*45
            ha=('right' if primary>0 else 'left') if inside else ('left' if primary>0 else 'right')
            ax.text(x,iy,f'{NAMES[group]}\n{row.Strong_DE_Genes:,}',fontsize=6.5,ha=ha,va='center',linespacing=.95)
    ax.axvline(0,color='#333333',lw=.65)
    ax.set(xlim=(-1450,1450),ylim=(2.65,-.65))
    ax.set_yticks(range(3),['C19/20\nvs C19','TWR\nvs C19','TWR\nvs C19/20'],fontsize=6.5)
    ax.set_xticks([-1000,-500,0,500,1000],['1,000','500','0','500','1,000'],fontsize=6.5)
    ax.tick_params(axis='y',length=0,pad=2)
    ax.spines[['top','right','left']].set_visible(False);ax.grid(axis='x',color='#E6E6E6',lw=.45)
    ax.set_xlabel('DE genes higher in labeled group',labelpad=2)
    CHECKS['h']={'all_counts':data.to_dict('records'),'contrast_order':contrasts,
                 'bar_encodings':'wide pale = strong DE; narrow dark = center-replicated strong DE'}
    save(fig)


def panel_i():
    data=read(D/'row4-source/figure4_hallmark_gsea_top4_nonredundant.tsv')
    assert len(data)==12
    fig=fig_new('i',88,50.5,'Hallmark enrichment')
    contrasts=['DN2 vs DN1','DN3 vs DN1','DN3 vs DN2']
    headers=[('C19','C19/20 vs C19','C19/20'),('C19','TWR vs C19','TWR'),
             ('C19/20','TWR vs C19/20','TWR')]
    shapes={'DN1':'o','DN2':'s','DN3':'^'}
    capped=np.minimum(data.Minus_Log10_FDR_Capped.to_numpy(),35)
    size_lo=float(capped.min());size_hi=float(capped.max())
    def size(v):
        # Same monotone, sqrt-transformed marker-size mapping as scale_size_continuous.
        return (3.0+(7.4-3.0)*np.sqrt(np.clip((np.asarray(v)-size_lo)/(size_hi-size_lo),0,1)))**2
    ov=canvas(fig)
    for k,(contrast,header) in enumerate(zip(contrasts,headers)):
        top=7.1+k*11.5
        ov.add_patch(Rectangle((29,top-2.4),58,2.2,facecolor='#F4F4F4',edgecolor='none'))
        txt(fig,34.2,top-1.3,header[0],6,va='center')
        txt(fig,58.2,top-1.3,header[1],6,ha='center',va='center')
        txt(fig,82.2,top-1.3,header[2],6,ha='right',va='center')
        # Draw arrows as vectors: the embedded Roboto Condensed lacks arrow glyphs.
        ov.annotate('',xy=(29.4,top-1.3),xytext=(33.4,top-1.3),
                    arrowprops={'arrowstyle':'-|>','color':'#444444','lw':.45,'mutation_scale':4})
        ov.annotate('',xy=(87,top-1.3),xytext=(83,top-1.3),
                    arrowprops={'arrowstyle':'-|>','color':'#444444','lw':.45,'mutation_scale':4})
        sub=data.loc[data.Contrast.eq(contrast)].sort_values(['NES','FDR'],ascending=[False,True])
        ax=axes(fig,29,top,58,8.5)
        for iy,row in enumerate(sub.itertuples()):
            ax.scatter(row.NES,iy,s=size(min(row.Minus_Log10_FDR_Capped,35)),marker=shapes[row.Enriched_Group],
                       facecolor=DN[row.Enriched_Group],edgecolor='#333333',lw=.45,alpha=.93)
        ax.set(xlim=(-3.5,3.5),ylim=(3.55,-.55))
        ax.set_yticks(range(4),sub.Pathway_Label,fontsize=6.5)
        ax.set_xticks([-2,0,2],['-2','0','2'],fontsize=6)
        ax.axvline(0,color='#555555',lw=.55)
        ax.grid(axis='x',color='#E6E6E6',lw=.4)
        ax.tick_params(axis='y',length=0,pad=2)
        ax.spines[['top','right','left']].set_visible(False)
        if k<2:
            ax.tick_params(axis='x',labelbottom=False,bottom=False)
            ax.spines['bottom'].set_visible(False)
        else:ax.set_xlabel('NES',labelpad=.5,fontsize=6.5)
    # One compact band retains every graphical encoding key.
    for x,g in zip([1.6,12.6,26.5],GROUPS):
        ov.scatter([x],[48],s=12,marker=shapes[g],facecolor=DN[g],edgecolor='#333333',lw=.45)
        txt(fig,x+1.8,48,NAMES[g],6,va='center')
    txt(fig,39,48,'-log10(FDR)',6,va='center')
    for x,value in zip([59,69,79],[5,20,30]):
        ov.scatter([x],[48],s=size(value),facecolor='white',edgecolor='#444444',lw=.45)
        txt(fig,x+2.2,48,str(value),6,va='center')
    CHECKS['i']={'pathways':data[['Contrast','Pathway','Pathway_Label','NES','FDR','Enriched_Group','Minus_Log10_FDR_Capped']].to_dict('records'),
       'pathway_count':12,'contrast_order':contrasts,'common_nes_limits':[-3.5,3.5],
       'point_size_encoding':'capped -log10(FDR), cap35; diameter uses original sqrt scaling',
       'point_shapes':shapes,'legend_values':[5,20,30]}
    save(fig)


def main():
    panel_a();panel_b();panel_c();panel_d();panel_e();panels_f_g();panel_h();panel_i()
    placements=[('a',4,13,56,49),('b',64,13,48,49),('c',116,13,60,49),
        ('d',4,64,172,47),('e',4,113,172,38),('f',4,153,96,37),('g',104,153,72,37),
        ('h',4,192,80,54),('i',88,192,88,54)]
    title='ecDNA architecture, EGFR specialization and RNA states in GBM'
    out=compose(5,title,250,placements)
    doc=fitz.open(out)
    doc[0].get_pixmap(matrix=fitz.Matrix(3,3),alpha=False).save(STAGE/'qa/M5_216dpi.png')
    audit={'figure':'M5','page_mm':[180,250],'sources':SOURCES,'panels':CHECKS,
        'font_bounds_checks':FONT_CHECKS,'minimum_text_pt':6,
        'scientific_changes':'None; frozen values and saved inference; display-only geometry and deterministic jitter.',
        'panel_letters':list('abcdefghi'),'palette':{'DN':DN,'states':STATE_COLORS,'contexts':CONTEXT,'fusion':FUSION},
        'output_sha256':sha(out)}
    def clean(value):
        if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
        if isinstance(value,list):return [clean(v) for v in value]
        if isinstance(value,float) and not np.isfinite(value):return None
        return value
    (STAGE/'qa/M5_audit.json').write_text(json.dumps(clean(audit),indent=2,allow_nan=False)+'\n')
    caption='''Figure 5. ecDNA architecture, EGFR specialization and RNA states in GBM.
a, Amplicon architecture among patients in each GBM evolutionary group. Each narrow tile is one patient with at least one amplicon of the indicated architecture; fill denotes the number of amplicons per patient. Percentages use the displayed patient denominators.
b, Distinct oncogenes per ecDNA-positive patient. All observations are shown, with values at or above 10 pooled at the labeled ≥10 display boundary. White boxes show median and interquartile range, and whiskers extend to the most extreme observations within 1.5 times the interquartile range. Horizontal displacement separates observations only. Brackets retain the saved significant BH-adjusted pairwise Wilcoxon q values.
c, Prevalence of the same 12 displayed ecDNA oncogenes among 106 C19, 90 C19/20 and 52 TWR ecDNA-positive patients. Cell labels and the continuous color scale denote percentages; upward/downward triangles identify higher/lower prevalence versus the remaining groups at saved BH q < 0.05.
d, Group-specific oncogene co-occurrence among 106, 94 and 53 ecDNA-positive specimens, respectively. Bars give specimen percentages by ecDNA context, using the original group-specific top-gene order. Nodes encode the number of unique ecDNA feature IDs. Arc color and width encode the number of specimens carrying both oncogenes; all saved pairs with at least three co-occurring specimens are shown. Same-specimen co-occurrence does not establish that two genes occur on the same amplicon. Patient and specimen denominators are intentionally different.
e, EGFR ecDNA maximum copy number on the original logarithmic scale, variance-stabilized (VST) EGFR expression, and EGFR fusion classes. Analysis-specific patient counts are shown. Box, point and q-value conventions match b; no copy-number or expression observations are omitted. Fusion percentages and category definitions are unchanged; “None” means no recorded EGFR fusion.
f, Fixed Neftel GSVA signature scores in the corrected cohort of 205 exact primary RNA/WGS matches (85 C19, 66 C19/20 and 54 TWR), standardized across that cohort for display and capped at z = ±2. All patients and the original corrected sample order are retained.
g, Percentage of patients assigned to each dominant bulk-RNA state, shown as horizontal 100% stacked bars. These are dominant bulk signatures, not measured cell proportions. The recomputed joint multinomial likelihood-ratio test adjusts for sequencing center and Battenberg purity (P = 0.000123677618448135). State colors and order are unchanged.
h, Strong differential-expression gene counts from the frozen center-adjusted contrasts. Pale wide bars show the primary counts; narrow saturated bars show the center-replicated subset. Direct labels identify the group with higher expression and the primary count; mirrored axis labels report absolute gene counts.
i, The same four nonredundant FDR-ranked Hallmark gene sets per contrast (12 rows total). Signed normalized enrichment scores (NES) retain their original direction; strip arrows identify the corresponding groups. Fill/shape encode the enriched group; point size encodes −log10(FDR), capped at 35. All statistical results and displayed pathway selections are frozen.
C19, chromosome 19 gain; C19/20, chromosome 19/20 co-gain; TWR, TP53/WGD-rich. AC, astrocyte-like; NPC, neural-progenitor-like; OPC, oligodendrocyte-progenitor-like; MES, mesenchymal-like. OXPHOS, oxidative phosphorylation; EMT, epithelial–mesenchymal transition. The poster informs presentation only: its older PASS-only ecDNA and 207-profile RNA cohorts are not used.
'''
    (STAGE/'captions/M5.txt').write_text(caption)
    print(out)
    print(json.dumps({k:{'render':v['outside_panel'],'native_pdf':v['native_pdf_outside_panel']}
        for k,v in FONT_CHECKS.items() if v['outside_panel'] or v['native_pdf_outside_panel']},indent=2))

