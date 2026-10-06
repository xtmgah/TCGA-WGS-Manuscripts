# Adapted from scripts/publishing/figure_standardization/figure5_space_palette_panels.py; publication side effects removed.
import render_figure5 as source
from render_figure5 import *
from common import *
import math
OUT=STAGE
NATIVE=STAGE


GEOMETRY={}


STATE_COLORS={'AC-like':'#367BB6','MES-like':'#D99A2A',
              'NPC-like':'#7C68A5','OPC-like':'#A8CFBD'}


def panel_e():
    cn=read(D/'panel_e_egfr_ecdna_copy_number.tsv')
    ex=read(D/'panel_e_egfr_expression.tsv')
    fu=read(D/'panel_e_egfr_fusion_classes.tsv')
    tests=read(D/'panel_e_pairwise_tests.tsv')
    fig=fig_new('e',96,66.5,'')
    widths=(96-4)/3
    title_centers=[19.6, widths+2+18.2, 2*(widths+2)+18.4]
    for x,title in zip(title_centers,['EGFR ecDNA'+'\ncopy number','EGFR'+'\nexpression','EGFR fusion'+'\nclasses']):
        txt(fig,x,.6,title,8,ha='center',linespacing=1.0)
    ax=axes(fig,9.5,12,20.2,47.5)
    ncn=boxes(ax,cn,'max_egfr_ecDNA_cn',log=True,size=5.5)
    ax.set_ylim(5,2600);ax.set_yticks([10,100,1000],['10','100','1,000'])
    ax.set_ylabel('Maximum copy number',labelpad=1)
    ctests=tests.loc[tests.metric.eq('EGFR ecDNA maximum copy number')&tests.fdr.lt(.05)]
    for row,y in zip(ctests.itertuples(),[520,1150]):
        bracket(ax,GROUPS.index(row.group_1),GROUPS.index(row.group_2),y,fmt_q(row.fdr),log=True)
    ax2=axes(fig,widths+2+6.8,12,22.8,47.5)
    nex=boxes(ax2,ex,'EGFR_expression',size=5.5)
    ax2.set_ylim(7.8,29);ax2.set_yticks([10,15,20,25]);ax2.set_ylabel('VST',labelpad=1)
    etests=tests.loc[tests.metric.eq('EGFR expression')&tests.fdr.lt(.05)]
    for row,y in zip(etests.itertuples(),[21.2,24.2,27.2]):
        bracket(ax2,GROUPS.index(row.group_1),GROUPS.index(row.group_2),y,fmt_q(row.fdr),dy=.16,text_dy=.05)
    ax3=axes(fig,2*(widths+2)+7,12,22.8,47.5)
    base=np.zeros(3)
    for name,color in FUSION.items():
        vals=fu.loc[fu.Architecture.eq(name)].set_index('DN_Group').pct_subjects.reindex(GROUPS).to_numpy()
        ax3.bar(range(3),vals,bottom=base,width=.67,color=color,edgecolor='white',lw=.45)
        for x,(bot,val) in enumerate(zip(base,vals)):
            ax3.text(x,bot+val/2,f'{val:.0f}%',fontsize=6,ha='center',va='center')
        base+=vals
    assert np.allclose(base,100)
    nf=fu.groupby('DN_Group').total_subjects.first().to_dict()
    ax3.set(ylim=(0,100),xlim=(-.5,2.5))
    ax3.set_yticks([0,25,50,75,100],['0%','25%','50%','75%','100%'],fontsize=6)
    ax3.set_xticks(range(3),[f'{NAMES[g]}\nn = {nf[g]}' for g in GROUPS],fontsize=6.5)
    for t,g in zip(ax3.get_xticklabels(),GROUPS):t.set_color('#000000')
    ax3.set_ylabel('Patients',labelpad=1);tidy(ax3)
    GEOMETRY['egfr'] = {'title_centers_mm': title_centers, 'axis_top_mm': 12,
        'axis_height_mm': 47.5, 'previous_axis_height_mm': 37.5,
        'axis_bottom_mm': 59.5, 'height_increase_percent': 100*(47.5/37.5-1),
        'fusion_legend': 'Existing outside-upper-right vector overlay retained'}
    CHECKS['e']={'copy_number_n':ncn,'expression_n':nex,'fusion_n':nf,
        'copy_number_scale':'log10','copy_number_values':cn[['Subject','DN_Group','max_egfr_ecDNA_cn']].to_dict('records'),
        'expression_values':ex[['Subject','DN_Group','EGFR_expression']].to_dict('records'),
        'fusion_counts':fu.to_dict('records'),'displayed_tests':pd.concat([ctests,etests]).to_dict('records')}
    save(fig)


def panel_i():
    data=read(B/'rna/analysis/GBM_state_analysis_metadata.tsv')
    composition=read(B/'rna/analysis/GBM_state_composition.tsv')
    tests=read(B/'rna/analysis/GBM_state_global.tsv')
    assert len(data)==205
    counts=data.groupby('DN_Group').size().reindex(GROUPS)
    assert counts.tolist()==[85,66,54]
    selected=tests.loc[tests.Purity_Source.eq('BB_Purity')&tests.Model.eq('joint')&tests.Test.eq('Multinomial LRT')&tests.Outcome.eq('Dominant state')]
    assert len(selected)==1 and int(selected.N.iloc[0])==205
    p=float(selected.P.iloc[0])
    fig=fig_new('i',172,25.5,'')
    txt(fig,78,.6,'Dominant RNA state',8,ha='center')
    txt(fig,78,3.85,f'Adjusted P = {p:.5f}',6.5,ha='center')
    # One extra millimetre uses the existing bottom margin, preserving bar size.
    ax=axes(fig,17,7.2,122,12.0)
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
                    color='#222222' if st in ['OPC-like','MES-like'] else 'white')
            recorded.append({'group':GROUPS[y],'state':st,'n':int(n),'pct':float(value)})
        left+=vals
    assert np.allclose(left,100)
    ax.set(xlim=(0,100),ylim=(2.6,-.6))
    ax.set_yticks(range(3),[f'{NAMES[g]} (n = {counts[g]})' for g in GROUPS],fontsize=6.5)
    ax.set_xticks([0,50,100],['0%','50%','100%'],fontsize=6)
    ax.tick_params(axis='y',length=0,pad=2);ax.spines[['top','right','left']].set_visible(False)
    ax.grid(axis='x',color='#E8E8E8',lw=.45)
    txt(fig,78,23,'Patients (%)',6.5,ha='center')
    ov=canvas(fig)
    for st,y in zip(STATE_ORDER,[8.2,12.3,16.4,20.5]):
        swatch(ov,144,y,STATE_COLORS[st],1.7,1.7);txt(fig,146.4,y,st,6,va='center')
    CHECKS['i']={'patients':205,'counts':recorded,'saved_adjusted_test':selected.to_dict('records'),
        'state_order':STATE_ORDER,'state_colors':STATE_COLORS,'orientation':'horizontal 100% stacked'}
    GEOMETRY['rna']={'title_center_x_mm':78,'title_top_mm':.6,'p_center_x_mm':78,'p_top_mm':3.85,
        'axis_top_mm':7.2,'axis_height_mm':12,'native_height_mm':25.5,
        'p_value':p,'p_label':f'Adjusted P = {p:.5f}',
        'percentage_text':'White on blue/purple; dark on orange/mint'}
    save(fig)

