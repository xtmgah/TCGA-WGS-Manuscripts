# Adapted from scripts/publishing/figure_standardization/figure5_vertical_panels.py; publication side effects removed.
import render_figure5 as source
from render_figure5 import *
from common import *
import math
OUT=STAGE
NATIVE=STAGE


def fig_new(letter,width,height,title):
    fig=source.fig_new(letter,width,height,title)
    fig.texts[0].set_x(.5);fig.texts[0].set_ha('center')
    return fig


def panel_d():
    nodes=read(D/'panel_d_cooccurrence_nodes.tsv')
    pairs=read(D/'panel_d_cooccurrence_pairs.tsv')
    bars=read(D/'panel_d_cooccurrence_context_prevalence.tsv')
    den=read(D/'panel_d_cooccurrence_tumor_denominators.tsv')
    fig=fig_new('d',172,63.5,'ecDNA oncogene co-occurrence')
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
        txt(fig,left+5,3.8,f'{NAMES[g]}  (n = {n} specimens)',6.5,color='#000000',weight='bold')
        ax=axes(fig,left+6.0,8.8,block-6.5,24.0)
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
        ar=axes(fig,left+6.0,42.8,block-6.5,10.5)
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
    locs=[(0,54.5),(54,54.5),(110,54.5),(0,57.7),(54,57.7)]
    for (label,color),(x,y) in zip(entries,locs):
        swatch(ov,x,y,color,1.8,1.8);txt(fig,x+2.7,y,label,6,va='center')
    txt(fig,0,61.6,'Pair count',6,va='center')
    for x,n in zip([15,28,41],[5,10,15]):
        ov.plot([x,x+5],[61.6,61.6],color=cmap(norm(n)),lw=edge_lw(n))
        txt(fig,x+6.2,61.6,str(n),6,va='center')
    txt(fig,64,61.6,'Unique ecDNA feature IDs',6,va='center')
    for x,n in zip([98,116,134,152],[20,40,60,80]):
        ov.scatter([x],[61.6],s=(node_dia(n)*MM)**2,facecolor='white',edgecolor='#222222',lw=.5)
        txt(fig,x+2.8,61.6,str(n),6,va='center')
    CHECKS['d']={'specimen_denominators':den.to_dict('records'),'nodes':nodes.to_dict('records'),
        'all_edges':pairs.to_dict('records'),'context_values':bars.to_dict('records'),
        'edges_minimum_count':3,'pair_count_range':[3,pair_max],'feature_id_range':[1,feature_max],
        'encodings':'arc color and width = same-specimen pair count; node size = unique feature IDs'}
    save(fig)

