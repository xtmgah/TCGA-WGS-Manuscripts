# Adapted from scripts/publishing/figure_standardization/refinement/render_figure2.py; publication side effects removed.
from common import *
from pathlib import Path
import json, csv, re, math
import numpy as np
import pandas as pd
from functools import lru_cache
setup_matplotlib()
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Wedge,Circle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


SOURCE_RNA = ROOT/'results/analysis/rna_adjustment_corrections_2026-09-16/analysis/ASTRO_plate_adjusted_DESeq2.tsv'


def volcano():
    r = pd.read_csv(SOURCE_RNA,sep='\t')
    r = r[r.padj.notna()].copy()
    r['y'] = -np.log10(r.padj.clip(lower=1e-300))
    fig = plt.figure(figsize=(46/25.4,43.5/25.4))
    ax = fig.add_axes([.20,.295,.775,.49])
    colors = np.where(r.padj>=.05,'#C7CDD6',np.where(r.log2FoldChange>0,'#315686','#8FAED6'))
    ax.scatter(r.log2FoldChange,r.y,c=colors,s=1.35,alpha=.65,edgecolors='none')
    threshold = -np.log10(.05)
    ax.axhline(threshold,color='#73869E',lw=.45,ls='--')
    ax.axvline(0,color='.5',lw=.4)
    ax.set(xlabel='Log2FC (CTR vs C17p)',ylabel='−log10(FDR)')
    ax.xaxis.labelpad=1
    ax.set_ylim(-.2,r.y.max()*1.18)
    ax.set_yticks([0,2,4,6,8,10,12]);ax.set_xticks([-1,0,1])
    ax.tick_params(length=2,width=.6,pad=1.5)
    top = r[(r.gene_type=='protein_coding')&(r.padj<.05)].nsmallest(3,'padj')
    assert list(top.Symbol)==['ETV4','EDA2R','CREB3L1']
    # Fixed collision-free positions preserve every original labeled gene.
    label_positions={'ETV4':(1.66,12.15),'EDA2R':(.97,9.15),'CREB3L1':(.73,7.60)}
    for _,row in top.iterrows():
        ax.annotate(row.Symbol,(row.log2FoldChange,row.y),xytext=label_positions[row.Symbol],
                    ha='right',va='center',fontsize=6,
                    arrowprops=dict(arrowstyle='-',color='#727272',lw=.35,shrinkA=1,shrinkB=1))
    ax.text(.02,threshold+.35,'FDR = 0.05',transform=ax.get_yaxis_transform(),fontsize=6,ha='left',va='bottom',color='#555555',bbox=dict(facecolor='white',edgecolor='none',pad=.2,alpha=.9))
    fig.text(.015,.975,'Differentially expressed genes',ha='left',va='top',fontsize=8)
    fig.text(.015,.90,'Plate, age, purity and grade adjusted',ha='left',va='top',fontsize=6)
    # Two compact rows retain every color meaning at a 46-mm panel width.
    for x,y,color,label in [(.14,.085,'#8FAED6','C17p-high'),(.62,.085,'#315686','CTR-high'),(.29,.035,'#C7CDD6','FDR ≥ 0.05')]:
        fig.add_artist(Line2D([x-.006],[y],marker='o',markersize=2.5,color=color,linestyle='none',transform=fig.transFigure))
        fig.text(x+.025,y,label,fontsize=6,ha='left',va='center')
    target = STAGE/'panels/M2_h.pdf'
    fig.savefig(target)
    stats = dict(source=str(SOURCE_RNA.relative_to(ROOT)),source_sha256=sha(SOURCE_RNA),plotted_genes=len(r),
                 labels=list(top.Symbol),x_range=list(ax.get_xlim()),y_range=list(ax.get_ylim()),
                 significant_CTR_high=int(((r.padj<.05)&(r.log2FoldChange>0)).sum()),
                 significant_C17p_high=int(((r.padj<.05)&(r.log2FoldChange<0)).sum()),
                 values_unchanged=True,refitted=False)
    plt.close(fig)
    (STAGE/'qa/M2_h_audit.json').write_text(json.dumps(stats,indent=2)+'\n')
    return stats

