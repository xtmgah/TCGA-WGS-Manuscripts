# Adapted from scripts/analysis/amplicon_classifier_2_update_2026-10-02/render.py; publication side effects removed.
import render_figure5 as r
from render_figure5 import *
from common import *
import figure5_vertical_panels as v
OUT=STAGE
D=ROOT/'results/analysis/amplicon_classifier_2_update_2026-10-02/current/figure5'


def new(letter,w,h,title):
    fig=r.fig_new(letter,w,h,title)
    fig.texts[0].set_x(.5);fig.texts[0].set_ha('center')
    return fig


def panel_a():
    data=read(D/'amplicon_subject_table.tsv');rates=read(D/'panel_a_amplicon_rates.tsv')
    fig=new('a',56,60.5,'Amplicon architecture')
    classes=[('Linear','n_linear'),('ecDNA','n_ecDNA'),('BFB','n_BFB'),('Complex non-cyclic','n_complex_non_cyclic'),('FAN','n_FAN')]
    for k,g in enumerate(GROUPS):
        ss=data[data.DN_Group.eq(g)];n=len(ss);top=6.4+k*16.2
        txt(fig,18,top-2.5,f'{NAMES[g]} (n = {n})',6.5,color='#000000',weight='bold')
        ax=axes(fig,18,top,35,13.1);ax.set(xlim=(0,84),ylim=(4.6,-.6))
        for j,(label,field) in enumerate(classes):
            vals=ss[ss[field].gt(0)].sort_values([field,'Subject'])
            for i,row in enumerate(vals.itertuples(index=False)):
                count=int(getattr(row,field));assert 1<=count<=4
                ax.add_patch(Rectangle((i*100/n,j-.35),100/n,.7,facecolor=AMP[count-1],edgecolor=(1,1,1,.25),lw=.1))
            pct=100*len(vals)/n
        ax.set_yticks(range(5),[c[0] for c in classes],fontsize=6)
        ax.set_xticks([0,25,50,75],['0%','25%','50%','75%'],fontsize=6)
        ax.tick_params(axis='y',length=0,pad=1.3);ax.spines[['top','right','left']].set_visible(False)
        ax.grid(axis='x',color='#E6E6E6',linewidth=.45)
        if k<2:ax.tick_params(axis='x',bottom=False,labelbottom=False);ax.spines['bottom'].set_visible(False)
    txt(fig,35,54.7,'Patients (%)',6.5,ha='center')
    ov=canvas(fig);txt(fig,.3,59.1,'Features / patient',6,va='center')
    for i,cl in enumerate(AMP):
        x=19+i*5.8;swatch(ov,x,59.1,cl,1.7,1.7);txt(fig,x+2.4,59.1,str(i+1),6,va='center')
    CHECKS['a']={'rates':rates.to_dict('records'),'classes':[c[0] for c in classes],'independent_classes':True,'height_mm':60.5}
    save(fig)


def panel_b():
    df=read(D/'panel_b_ecdna_oncogene_burden.tsv');tests=read(D/'panel_b_ecdna_oncogene_burden_pairwise_tests.tsv')
    fig=new('b',48,60.5,'Oncogene burden');ax=axes(fig,8,5.8,39,47.5)
    counts=boxes(ax,df,'n_distinct_oncogenes',cap=10,size=4.5)
    ax.set_ylim(-.2,11.9);ax.set_yticks([0,2,4,6,8,10],['0','2','4','6','8','≥10']);ax.set_ylabel('Distinct oncogenes',labelpad=1.5)
    sig=tests[tests.fdr.lt(.05)]
    for i,row in enumerate(sig.itertuples()):bracket(ax,GROUPS.index(row.group_1),GROUPS.index(row.group_2),10.45+i*.9,fmt_q(row.fdr),dy=.1,text_dy=.04)
    CHECKS['b']={'counts':counts,'tests':sig.to_dict('records')};save(fig)


def panel_c():
    # Retain the existing drawing semantics while using the extra top-row height.
    old_new,old_axes,old_txt=r.fig_new,r.axes,r.txt
    def taller(letter,w,h,title):return new(letter,w,h+15,title)
    # Avoid recursive use of r.fig_new in our centered wrapper.
    def taller(letter,w,h,title):
        f=old_new(letter,w,h+15,title);f.texts[0].set_x(.5);f.texts[0].set_ha('center');return f
    def ax_shift(fig,x,top,w,h):return old_axes(fig,x,top+(15 if top>=37 else 0),w,h+(15 if top==5.3 else 0))
    def txt_shift(fig,x,y,value,*a,**k):return old_txt(fig,x,y+(15 if y>=37 else 0),value,*a,**k)
    r.fig_new=taller;r.axes=ax_shift;r.txt=txt_shift
    # Canvas marker positions also move, using a translated data coordinate frame.
    old_canvas=r.canvas
    def lower_canvas(fig):
        ax=old_canvas(fig);ax.set_ylim(fig._mm[1]-15,-15);return ax
    r.canvas=lower_canvas
    try:r.panel_c()
    finally:r.fig_new=old_new;r.axes=old_axes;r.txt=old_txt;r.canvas=old_canvas


def panel_e():
    cn=read(D/'panel_e_egfr_ecdna_copy_number.tsv');tests=read(D/'panel_e_pairwise_tests.tsv')
    fig=new('e',31,66.5,'');txt(fig,19.6,.6,'EGFR ecDNA\ncopy number',8,ha='center',linespacing=1.)
    ax=axes(fig,9.5,12,20.2,47.5);counts=boxes(ax,cn,'max_egfr_ecDNA_cn',log=True,size=5.5)
    ax.set_ylim(5,2600);ax.set_yticks([10,100,1000],['10','100','1,000']);ax.set_ylabel('Maximum copy number',labelpad=1)
    sig=tests[tests.metric.eq('EGFR ecDNA maximum copy number')&tests.fdr.lt(.05)]
    for row,y in zip(sig.itertuples(),[520,1150]):bracket(ax,GROUPS.index(row.group_1),GROUPS.index(row.group_2),y,fmt_q(row.fdr),log=True)
    CHECKS['e']={'counts':counts,'tests':sig.to_dict('records'),'plot_height_mm':47.5};save(fig)

