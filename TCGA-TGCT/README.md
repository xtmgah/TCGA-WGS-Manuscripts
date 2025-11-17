
## A Thiopurine-like Mutagenic Process Defines TGCT Subtypes

This repository contains R code used to generate primary and supplementary visualizations for the TCGA-TGCT project. It supports the manuscript "A Thiopurine-like Mutagenic Process Defines TGCT Subtypes".

The repository includes:
- R scripts for main figures
- Reference data files for demo
- R Markdown files and rendered HTML outputs for reproducibility

**Note: To properly view the HTML files with full interactivity and formatting, download and open them locally on your computer.**

Software Dependencies:
- R version ≥ 4.0
- All required R packages are listed within each R script, R Markdown (.Rmd), and rendered HTML output files.

Hardware Requirements:
The R scripts are designed to run on standard desktop or laptop machines.

## Figures included ## 

**Fig. 1: Distinct Genomic Landscapes of Seminomas and Non-Seminomas in Testicular Germ Cell Tumors (TGCTs).** 
(a), Comparison of key genomic features between seminomas (blue) and non-seminomas (purple), including tumor mutational burden (TMB; mutations/Mb), percentage genome altered (PGA; WGD-adjusted), number of structural variants (SVs), and number of transposable element insertions (TEs). P-values were obtained from linear regression models adjusted for age and tumor purity. For TMB, tumor ploidy was also included as a covariate. Multiple testing was corrected using the Benjamini–Hochberg FDR method.
(b), Significantly mutated genes identified by the IntOGen pipeline across all TCGT samples. Dot size reflects combined -log10 q-value based on different driver gene algorithms. Colors denote functional annotation (Act: activating; LoF: loss-of-function).
(c), Copy number heatmap of chromosome 12 showing frequent 12p amplification and isochromosome 12p [i(12p)] in TGCTs.
(d), Enrichment of i(12p) in seminomas with KIT mutations. P-value and odds ratio from two-sided Fisher’s exact test shown above the barplot.
(e), GISTIC2.0 analysis of recurrent focal copy number alterations. Significant amplifications (red) and deletions (blue) are annotated with potential target genes. Y-axis represents the significance by q value (left) and normalized amplification signals (G-Score; right). The green line represents the significance cutoff at Q value=0.25.
(f), Oncoprints depicting recurrent driver mutations, focal amplifications/deletions, and genes disrupted by SVs in seminomas and non-seminomas. Mutation frequencies for each gene are shown on the right. For mutations, only recurrent nonsynonymous mutations in known pan-cancer driver genes were included. For focal somatic copy number alterations (SCNAs), only high-level amplifications (copy number >= 2) and deep deletions (copy number = –2) were considered. For genes disrupted by SVs, we included only recurrent events involving known pan-cancer driver genes.
All boxplots display median, interquartile range (IQR), and whiskers extending to 1.5× IQR.


**Fig. 2: Widespread Chromosome X Amplification and X Chromosome Inactivation (XCI) in TGCTs.**
(a) Tumor purity-adjusted sequencing depth ratio of chromosome X versus autosomes in matched normal blood, seminomas, and non-seminomas. P-values from two-sided Wilcoxon rank-sum tests are shown above boxplots.
(b) Copy number profiles of chromosome X in TGCTs, separated by seminoma (left) and non-seminoma (right), ordered by total copy number. The y-axis indicates chromosome X positions. Top panel: clonal somatic copy number alterations (SCNAs); bottom panel: subclonal SCNAs. Subclonal fractions per sample are shown in the barplot below each panel.
(c) XIST mRNA expression across TCGA tumor types (male subjects), with TGCTs showing the highest levels.
(d) XIST mRNA expression in normal testis (GTEx), non-seminomas, and seminomas (TCGA). P-values from two-sided Wilcoxon rank-sum tests are shown above boxplots.
(e) Negative correlation between XIST expression and chromosome X median DNA methylation in seminomas (left) and non-seminomas (right). Pearson correlation coefficient (R) and P-value are shown below.
(f) Replication stress scores across TCGA tumor types based on replication stress gene signatures, with TGCTs showing the highest levels.
(g) Correlation between XIST expression and replication stress scores in seminomas (left) and non-seminomas (right), with significant positive correlation in seminomas only. Pearson R and P-value are shown above. 
All boxplots display median, IQR, and whiskers extending to 1.5× IQR.


**Fig. 3: Distinct Mutational Processes in TGCT Subtypes.**
(a) Unsupervised hierarchical clustering of TGCTs based on collective single-base substitution (SBS) and insertion-deletion (ID) mutational profiles, separating seminomas (blue) and non-seminomas (purple).
(b) SBS mutational signature decomposition, with stacked barplots showing relative contributions per signature (left y-axis). A line with dots indicates the number of mutations (right y-axis). Cosine similarity between original and decomposed profiles is shown below.
(c) ID mutational signature decomposition, with stacked barplots showing relative contributions per signature (left y-axis). A line with dots indicates the number of mutations (right y-axis). Cosine similarity between original and decomposed profiles is shown below.
(d) Correlation between mutational burden of ID signatures (ID1, ID2, ID9) and replication stress scores, with significant positive correlations in seminomas only. Pearson R and P-value are shown above scatterplots.
(e) Comparison of SBS mutational signature burden between seminomas and non-seminomas. Statistical significance (FDR) was assessed by linear regression, adjusted for tumor purity, with Benjamini–Hochberg correction.
(f) Comparison of ID mutational signature burden between subtypes, with non-seminomas showing higher ID1 and ID6 activity. Statistical significance (FDR) was assessed by linear regression, adjusted for tumor purity, with Benjamini–Hochberg correction.
All boxplots display median, IQR, and whiskers extending to 1.5× IQR.


**Fig. 4: Impact of Thiopurine-Associated Mutational Signature SBS87 on Tumor Latency and Telomere Length in TGCTs.**
(a) Comparison of estimated age of the most recent common ancestor (MRCA, left) and tumor latency (time from MRCA to diagnosis, right) between seminomas and non-seminomas. P-values from two-sided Wilcoxon rank-sum tests are shown above boxplots.
(b) Correlation between SBS87 activity (log2-transformed) and tumor latency in seminomas (left) and non-seminomas (right). Pearson R and P-value are shown above scatterplots.
(c) Linear regression analysis of mutational signatures associated with tumor latency in all TGCT samples, adjusted for tumor purity. FDR from Benjamini–Hochberg correction is shown, with FDR<0.05 thresholds marked by a red dashed line.
(d) Comparison of telomere length (TL) in normal blood, tumor, and tumor/normal TL ratio (log2) between seminomas and non-seminomas. P-values above boxplots are from linear regression, adjusted for tumor purity and age, comparing seminomas and non-seminomas; P-values below are from two-sided Wilcoxon rank-sum tests comparing tumor versus normal TL. 
(e) Linear regression analysis of mutational signatures associated with tumor TL in non-seminomas, adjusted for tumor purity and age. FDR from Benjamini–Hochberg correction is shown, with FDR<0.05 thresholds marked by a red dashed line.
(f) Correlation between SBS87 activity and TL in normal tissues (left), tumors (middle), and tumor/normal TL ratio (right), stratified by subtype. Significant positive correlations with tumor TL are observed in non-seminomas only. Pearson R and P-value are shown above scatterplots.
All boxplots display median, IQR, and whiskers extending to 1.5× IQR.


**Fig. 5: Early and Widespread Whole-Genome Doubling (WGD) and Clonal Architecture in TGCTs.**
(a) Distribution of tumors with minimal cancer cell fraction (CCF) of subclones in seminomas and non-seminomas, with non-seminomas showing higher intratumoral heterogeneity.
(b) Comparison of clonal versus subclonal mutation proportions between seminomas and non-seminomas.
(c) Ternary plot illustrating genome-wide major copy number (MCN) composition to infer WGD status. Each point represents a tumor, colored by subtype and shaped by WGD category. Tumors are classified as non-WGD, single WGD, or multiple WGD based on the relative proportions of the autosomal genome with MCN = 1, 2, or ≥3.
(d) Timing of WGD events relative to the accumulation of clock-like mutations (SBS1 and SBS5). Boxplots show the ratio of mutations in genomic segments with major copy number (MCN) = 2 versus MCN = 1 (indicative of single WGD), or MCN = 3 versus MCN = 2 (indicative of multiple WGDs). High mutation ratios in duplicated segments suggest that WGD occurred early during tumor evolution in both seminomas and non-seminomas.
(e) Correlation between WGD timing and early clonal mutation proportion in seminomas (left) and non-seminomas (right). Tumors harboring early clonal driver mutations prior to WGD are annotated by color, with counts and frequencies indicated above the plots. Pearson R and P-value are shown above scatterplots. 
(f) Estimated cell divisions before WGD in seminomas versus non-seminomas, inferred from pre-WGD mutation burden. P-values from two-sided Wilcoxon rank-sum tests are shown above boxplots.
(g) Correlation between pre-WGD cell divisions and patient age at diagnosis in seminomas (left) and non-seminomas (right). Pearson R and P-value are shown above scatterplots.
All boxplots display median, IQR, and whiskers extending to 1.5× IQR.


**Fig. 6: Germ Cell-Like Transcriptomic Program and Unique APOBEC Mutagenesis Suppression in TGCTs.**
(a) Differential enrichment of Hallmark gene sets in TGCTs versus other TCGA tumor types using ssGSEA scores. Enriched pathways include cell cycle, DNA repair, and spermatogenesis; depleted pathways include metabolic and immune-related processes. Circle size indicates −log10 FDR from two-sided Wilcoxon rank-sum tests; color denotes median ssGSEA score enrichment in TGCTs. Benjamini–Hochberg correction was applied.
(b) Radar plot of TGCT-specific enrichment for select Hallmark pathways (SPERMATOGENESIS, DNA_REPAIR, E2F_TARGETS, G2M_CHECKPOINT, MITOTIC_SPINDLE, MYC_TARGETS_V2) across TCGA tumor types.
(c) Differential enrichment of KEGG gene sets in TGCTs, with top pathways including HOMOLOGOUS_RECOMBINATION, MISMATCH_REPAIR, and BASE_EXCISION_REPAIR. Circle size indicates −log10 FDR from two-sided Wilcoxon rank-sum tests; color denotes median ssGSEA score enrichment in TGCTs. Benjamini–Hochberg correction was applied.
(d) Radar plot of ssGSEA enrichment for core DNA repair KEGG pathways across TCGA tumor types, with TGCTs ranking highest.
(e) UNG (uracil DNA glycosylase) mRNA expression across TCGA tumors, with TGCTs showing the highest levels and seminomas higher than non-seminomas. P-values from two-sided Wilcoxon rank-sum tests are shown above boxplots. Boxplots display median, IQR, and whiskers extending to 1.5× IQR.
(f) Pearson correlation between APOBEC3B and UNG expression across TCGA tumors, with TGCTs showing a unique negative correlation. FDR from Benjamini–Hochberg correction is shown.
(g) Scatterplot of negative association between APOBEC3B and UNG mRNA expression in TGCTs. Pearson R and P-value are shown below.

---

## Instructions for Reproducibility ##

All R scripts are written in a generic way to allow users to reproduce the figures presented in this tumor evolution paper. To run the code successfully, please follow the folder structure and guidelines below:

Input Data: Place all required input data files in the data/ directory.

Functions: The functions/ folder contains additional helper functions required for visualization. These will be sourced automatically by the scripts.

**Place all data, Rscripts, and functions in the same directory as the R script or R Markdown file.** 

Note:
The data/ folder is expected to contain most of the whole-genome sequencing (WGS) analysis results from the TCGA-TGCT study. These are large files and are not included in this repository due to data access restrictions. To obtain the raw WGS data of TCGA-TGCT study, users must apply through dbGaP (application details are provided in our manuscript). All sequencing raw data have been deposited in dbGaP in accordance with NIH guidelines and data use restrictions. For any questions, please contact the Principal Investigator (Dr. Tongwu Zhang) of the TCGA-TGCT WGS analysis directly.


## Data format
You can find all the data format in the html file "data_format.html"

