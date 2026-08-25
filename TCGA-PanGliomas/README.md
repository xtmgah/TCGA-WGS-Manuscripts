
## Endogenous mutational mechanisms and metabolic context shape endometrial cancer

This repository contains R code used to generate primary and supplementary visualizations for the TCGA-UCEC project. It supports the manuscript "Deep whole-genome sequencing reveals distinct mutational mechanisms shaping endometrial cancer".

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

**Fig. 1: Whole-genome sequencing defines distinct genomic architectures across UCEC subtypes.** 
a, Distribution of UCEC tumors classified into four molecular subtypes based on WGS-derived genomic features: POLE (n = 41), MSI (n = 136), CN-High (n = 173), and CN-Low (n = 90).
b, Genome-wide alteration metrics across UCEC molecular subtypes, including tumor mutational burden (TMB; log₂ mutations per megabase), proportion of genome altered (PGA; WGD-adjusted), structural variant (SV) burden (log₂ number of events + 1), and transposable element (TE) insertion burden (log₂ number of events + 1). Each point represents one tumor. Boxes show median and interquartile range; whiskers denote 1.5× interquartile range. Two-sided Wilcoxon rank-sum tests with FDR correction were used.
c, Genome-wide tumor mutational burden stratified by mutation class, showing separate contributions from single-nucleotide variants (SNVs) and small insertions and deletions (indels) across UCEC molecular subtypes (shown similarly to b). 
d, Fraction of WGD positive and WGD negative tumors within each molecular subtype.
e, Relative mitochondrial DNA (mtDNA) copy number across UCEC subtypes, expressed as log₂(tumor/normal mtDNA ratio). Statistical significance was assessed using two-sided Wilcoxon rank-sum tests with FDR correction.
f, Significantly mutated driver genes identified using an integrative driver discovery framework. The upper panel shows the number of tumors harboring mutations in each gene, with subtype-specific frequencies indicated. Middle annotations summarize inferred functional roles, novelty, and combined statistical significance (−log₁₀ q value). The lower panel shows mutation class composition and support across multiple driver detection methods.
g, Subtype-specific properties of driver gene alterations, including the proportion of genes affected by single versus multiple-hit mutations, distribution of mutation classifications, and fraction of clonal versus subclonal events.


**Fig. 2: LINE-1 activity associated with structural variation, chromothripsis, and ecDNA formation in UCEC.**
a, Genome-wide distribution of somatic LINE-1-mediated transductions across autosomes and sex chromosomes, shown as the average number of transduction events per sample in fixed genomic bins. Circos plots depict transduction events originating from recurrent germline LINE-1 source elements, including chr22q12.1 and chrXp22.2, with arcs indicating source-to-insertion relationships.
b, Genome-wide profiles of structural variant (SV) breakpoints across UCEC molecular subtypes. SV breakpoint density is shown as the number of breakpoints per 5 Mb window, averaged across tumors within each subtype (POLE, MSI, CN-High, and CN-Low), and plotted along chromosomal coordinates.
c, Summary of focal amplification architectures detected across molecular subtypes. Bars indicate the number of tumors harboring focal amplification events classified as linear, ecDNA, BFB, or complex non-cyclic structures, stratified by subtype.
d, Characterization of oncogene-containing ecDNA amplicons. Stacked bar plots show the number of subjects harboring ecDNA-associated oncogenes, stratified by genomic context of ecDNA formation. Curved links illustrate pairwise co-occurrence of oncogenes detected on the same ecDNA amplicon, with line width and color corresponding to the number of shared events and node size indicating the number of subjects.
e, Relationship between ecDNA burden and transposable element insertion counts. Tumors are grouped by the number of detected ecDNA amplicons (0, 1, or ≥2), and TE insertion burden is shown as log₂(TE count + 1). Boxes indicate median and interquartile range; whiskers represent 1.5× interquartile range. Statistical significance was assessed using a linear regression.  
f, Frequency of chromothripsis events across molecular subtypes. Stacked bars indicate the percentage of tumors with chromothripsis, stratified by confidence level (high or low) based on established detection criteria.
g, Comparison of transposable element (TE) insertion burden between tumors with and without chromothripsis. TE burden is shown as log₂(TE count + 1) per tumor. Each dot represents an individual tumor. Boxes indicate median and interquartile range; whiskers represent 1.5× interquartile range. Statistical significance was assessed using a two-sided Wilcoxon rank-sum test.


**Fig. 3: Subtype-specific mutational signature landscapes in UCEC.**
a, SBS mutational signature decomposition across UCEC tumors, ordered by molecular subtype. Stacked bars show the relative contribution of each SBS signature per tumor (left axis), with total SBS mutation counts shown as a line plot (right axis). Cosine similarities between the original and reconstructed mutational profiles are shown below.
b, DBS mutational signature decomposition across UCEC tumors, displayed similar to a, with stacked bars indicating relative signature contributions and total DBS mutation counts shown as a line plot. Cosine similarity values are shown below.
c, Prevalence of SBS, ID, and DBS mutational signatures across UCEC molecular subtypes. Bubble size and color indicate the proportion of mutations (top) or the proportion of tumors (bottom) attributed to each signature within each subtype. A mutational signature was considered present in a tumor if at least 50 mutations were assigned to the signature and it contributed more than 5% of total mutations.
d, De novo mutational signature profile of DBS78C identified in mismatch repair-deficient tumors, shown as relative contributions across doublet substitution contexts.


**Fig. 4: Indel mutational signatures distinguish UCEC molecular subtypes.**
a, ID mutational signature decomposition across UCEC tumors, ordered by molecular subtype. Stacked bars show the relative contribution of each ID signature per tumor (left axis), with total ID mutation counts shown as a line plot (right axis). Cosine similarities between the original and reconstructed ID profiles are shown below.
b, Comparison of ID2-to-ID1 signature activity ratios across UCEC molecular subtypes. Each point represents an individual tumor. Boxplots indicate the median and interquartile range. Statistical significance was assessed using two-sided Wilcoxon rank-sum tests with FDR correction.
c, ID2-to-ID1 ratios stratified by dominance of the MSI-associated single-base substitution signature SBS44. Tumors were grouped as SBS44-dominant (>50% contribution), SBS44 non-dominant, or SBS44-absent. Statistical significance was assessed using two-sided Wilcoxon rank-sum tests with FDR correction.
d, Correlation between ID2-to-ID1 ratio and SBS44 mutation burden across UCEC tumors. Pearson correlation coefficient (R) and two-sided P value are shown; the shaded region denotes the 95% confidence interval of the fitted linear regression.
e, Replication timing distribution of ID2 mutations across UCEC molecular subtypes. Bars indicate normalized densities of indel mutations across replication timing bins from early to late replicating regions; dashed lines represent smoothed trends.
f, Comparison of indel-derived neoantigen burden across UCEC molecular subtypes. Each point represents an individual tumor. Boxplots summarize distributions, with statistical significance assessed using two-sided Wilcoxon rank-sum tests and Benjamini-Hochberg correction; adjusted P values are shown.
g, Schematic illustration of replication slippage mechanisms underlying ID1 and ID2 mutations. ID1 mutations arise from slippage of the nascent strand, whereas ID2 mutations result from slippage of the template strand. The relative contributions of ID1 and ID2 in MSI and non-MSI tumors are summarized.



**Fig. 5: Transcriptomic features of CN-Low tumors and their association with body mass index**
a, Replication stress Z-scores across UCEC molecular subtypes. Each point represents an individual tumor. Boxplots indicate the median and interquartile range. Pairwise statistical comparisons were performed using two-sided Wilcoxon rank-sum tests with FDR correction.
b, Proliferation scores across UCEC molecular subtypes, derived from proliferation-related gene expression signatures. Each point represents one tumor. Boxplots summarize the distribution within each subtype. FDR-adjusted P values are shown.
c, Expression of XIST across UCEC molecular subtypes. Each point represents one tumor. Boxplots indicate median and interquartile range. FDR-adjusted P values for pairwise comparisons are shown.
d, Distribution of body mass index (BMI) across UCEC molecular subtypes. Each point represents an individual patient. Boxplots summarize the distribution within each subtype. FDR-adjusted P values for pairwise comparisons are indicated.
e-f, Multivariable regression analysis assessing associations with tumor mutational burden in TCGA-UCEC (e; N=409) and CPTAC-UCEC (f; N=91). Variables include UCEC molecular subtype (reference = CN-Low), age at diagnosis, tumor purity, and BMI class (reference = normal). Points represent regression coefficients (β), and horizontal lines indicate 95% confidence intervals. Two-sided P values are shown.

---

## Instructions for Reproducibility ##

All R scripts are written in a generic way to allow users to reproduce the figures presented in this tumor evolution paper. To run the code successfully, please follow the folder structure and guidelines below:

Input Data: Place all required input data files in the data/ directory.

Functions: The functions/ folder contains additional helper functions required for visualization. These will be sourced automatically by the scripts.

**Place all data, Rscripts, and functions in the same directory as the R script or R Markdown file.** 

Note:
The data/ folder is expected to contain most of the whole-genome sequencing (WGS) analysis results from the TCGA-UCEC study. These are large files and are not included in this repository due to data access restrictions. To obtain the raw WGS data of TCGA-UCEC study, users must apply through dbGaP (application details are provided in our manuscript). All sequencing raw data have been deposited in dbGaP in accordance with NIH guidelines and data use restrictions. For any questions, please contact the Principal Investigator (Dr. Tongwu Zhang) of the TCGA-UCEC WGS analysis directly.


## Data format
You can find all the data format in the html file "data_format.html"

## WGS analysis pipeline

The Whole Genome Sequencing (WGS) pipeline used in this study is available at the following repository: [WGS-Pipeline](https://github.com/xtmgah/Sherlock-Lung/tree/main/Sherlock-Lung%20Pipeline)

