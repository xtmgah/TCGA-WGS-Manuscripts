# Adapted from scripts/publishing/figure_standardization/figure5_row3_panels.py; publication side effects removed.
import render_figure5 as source
from render_figure5 import *
from common import *
import math
OUT=STAGE
NATIVE=STAGE


GEOMETRY={}


def panel_i():
    data=read(D/'row4-source/figure4_hallmark_gsea_top4_nonredundant.tsv')
    assert len(data)==12
    fig=fig_new('i',72,66.5,'')
    txt(fig,50,.6,'Hallmark\nenrichment',8,ha='center',linespacing=1.0)
    contrasts=['DN2 vs DN1','DN3 vs DN1','DN3 vs DN2']
    headers=[('C19','C19/20 vs C19','C19/20'),('C19','TWR vs C19','TWR'),
             ('C19/20','TWR vs C19/20','TWR')]
    shapes={'DN1':'o','DN2':'o','DN3':'o'}
    capped=np.minimum(data.Minus_Log10_FDR_Capped.to_numpy(),35)
    size_lo=float(capped.min());size_hi=float(capped.max())
    def size(v):
        # Same monotone, sqrt-transformed marker-size mapping as scale_size_continuous.
        return (3.0+(7.4-3.0)*np.sqrt(np.clip((np.asarray(v)-size_lo)/(size_hi-size_lo),0,1)))**2
    ov=canvas(fig)
    for k,(contrast,header) in enumerate(zip(contrasts,headers)):
        header_top=8.4+k*15.3
        top=header_top+3.6
        if k:assert header_top > (8.4+(k-1)*15.3)+3.6+10.5
        ov.add_patch(Rectangle((29,header_top),42,3.2,facecolor='#F4F4F4',edgecolor='none'))
        # Direct negative/positive group labels replace the duplicated contrast title.
        txt(fig,34.2,header_top+1.6,header[0],6.5,va='center')
        txt(fig,66.4,header_top+1.6,header[2],6.5,ha='right',va='center')
        ov.annotate('',xy=(29.4,header_top+1.6),xytext=(33.4,header_top+1.6),
                    arrowprops={'arrowstyle':'-|>','color':'#444444','lw':.45,'mutation_scale':4})
        ov.annotate('',xy=(71,header_top+1.6),xytext=(67,header_top+1.6),
                    arrowprops={'arrowstyle':'-|>','color':'#444444','lw':.45,'mutation_scale':4})
        sub=data.loc[data.Contrast.eq(contrast)].sort_values(['NES','FDR'],ascending=[False,True])
        ax=axes(fig,29,top,42,10.5)
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
    # A dedicated footer contains two aligned, horizontal keys below the axes.
    txt(fig,4,61,'Enriched group',6.5,va='center')
    for x,g in zip([28,44,60],GROUPS):
        ov.scatter([x],[61],s=16,marker=shapes[g],facecolor=DN[g],edgecolor='#333333',lw=.45)
        txt(fig,x+2.2,61,NAMES[g],6.5,va='center')
    txt(fig,4,64.5,'-log10(FDR)',6.5,va='center')
    for x,value in zip([28,44,60],[5,20,30]):
        ov.scatter([x],[64.5],s=size(value),facecolor='white',edgecolor='#444444',lw=.45)
        txt(fig,x+2.2,64.5,str(value),6.5,va='center')
    GEOMETRY['enrichment']={'title_center_mm':50,'title_top_mm':.6,'header_tops_mm':[8.4,23.7,39.0],
        'axis_tops_mm':[12.0,27.3,42.6],'axis_height_mm':10.5,'axis_bottom_mm':53.1,
        'direction_headers':[{'contrast':c,'negative_NES_group':h[0],'positive_NES_group':h[2]} for c,h in zip(contrasts,headers)],
        'legend':{'orientation':'two horizontal rows','row_centers_mm':[61,64.5],'column_centers_mm':[28,44,60]}}
    CHECKS['i']={'pathways':data[['Contrast','Pathway','Pathway_Label','NES','FDR','Enriched_Group','Minus_Log10_FDR_Capped']].to_dict('records'),
       'pathway_count':12,'contrast_order':contrasts,'common_nes_limits':[-3.5,3.5],
       'point_size_encoding':'capped -log10(FDR), cap35; diameter uses original sqrt scaling',
       'point_shapes':shapes,'legend_values':[5,20,30]}
    save(fig)

