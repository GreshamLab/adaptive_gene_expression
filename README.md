# Supplemental Code for "Changes in gene expression at multiple levels in CNV containing adapted populations" 

---

## Introduction
This is an orienting document intended to describe the codes used in the analysis of data for the publication of "Changes in gene expression at multiple levels in CNV containing adapted populations". All files are available in this repository, download unzip, and move all files into a 'data' directory to use the commands provided here.

### Expression Alignment to Ensembl reference genome and reads per gene

We used cutadapt for trimming, followed by two-rounds of STAR, first aligning to ncRNA to filter out, then aligning to Ensembl for alignment to reference genome. 
```{}
Process_STAR_Carolino_v7_RNA_vSGD_CDS.sh
Process_STAR_Carolino_v7_RPF_vSGD_CDS.sh
```

This was followed by the reads per gene calculated by bedtools coverage -counts, using strand specificity.
```{}
run_bedtools_coverage_RNA_RPF_Ensembl_counts.sh
```

### DESeq2 performed on RNA, RPF 
To identify genes with different relative abundances DESeq2 was ran on untransformed, un-normalized reads as per documentation.

```{}
DESeq_difObs_RNA_v1.1.R
DESeq_difObs_RPF_v1.1.R

# these rely on: 
analyses/chemostat_expression/RNA_unnormed_read_counts_by_gene.tsv
analyses/chemostat_expression/RPF_unnormed_read_counts_by_gene.tsv

# these produce:
analyses/deseq/DESeq_Obs_RNA_DGY1657_DGY1726.txt
analyses/deseq/DESeq_Obs_RNA_DGY1657_DGY1735.txt
analyses/deseq/DESeq_Obs_RNA_DGY1657_DGY1741.txt
analyses/deseq/DESeq_Obs_RNA_DGY1657_DGY1743.txt
#
analyses/deseq/DESeq_Obs_RPF_DGY1657_DGY1726.txt
analyses/deseq/DESeq_Obs_RPF_DGY1657_DGY1735.txt
analyses/deseq/DESeq_Obs_RPF_DGY1657_DGY1741.txt
analyses/deseq/DESeq_Obs_RPF_DGY1657_DGY1743.txt

```
### Build "Unit" object
Certain modes of analysis require or are improved by the expression data being similarly scaled. To acheive this we first transform the data using Box-Cox Power Transform followed by a MinMax transform into the Unit range [0,1], (see Figure S5). Unit transformed data are used for Figures 2A and Supplemental Figure S3.

```{}
scale_data_v0.13_Unit.py

# input
analyses/chemostat_expression/expression_dict.tsv
metadata/aa_protein_lengths.tsv
metadata/chemostat_gene_relative_copy_number.tsv
analyses/ms/DGY1726_v_DGY1657_wCON_wBatch_QN_p0.01.txt
analyses/ms/DGY1735_v_DGY1657_wCON_wBatch_QN_p0.01.txt
analyses/ms/DGY1741_v_DGY1657_wCON_wBatch_QN_p0.01.txt
analyses/ms/DGY1743_v_DGY1657_wCON_wBatch_QN_p0.01.txt

# output
analyses/efficiency/ms_counts_expression.tsv
analyses/efficiency/ratio_df_Unit_v13.tab
```

### Perform Unit Transform Test
Unit transformed data are evaluated using the Unit transform test, described in Supplemental Figure 6, 7, 8, 9. The results of which inform the Protein Expression Efficieny analysis (Figure 4A).
```{}
unit_differential_v16.py

# input
analyses/efficiency/ratio_df_Unit_v13.tab

# output
analyses/efficiency/unit_object_level_1v13_Unit_10pct.tab
analyses/efficiency/unit_object_level_2v13v16_Unit_10pct.tab

```

### Make Heatmap and cluster (Figure 2)
Takes Unit transformed expression data, builds 

```{}
build_expression_object_l1_v1.0.py
heatmap_v2_l1.r

# input
data/analyses/efficiency/unit_object_level_2v13v16_Unit_10pct.tab
data/analyses/ms/Perseus_DA_Welchs_t-test_wFDR_edit.txt

# output
data/output/unit_object_level_Unit_l1_tests.tab
figures/unit_object_select_24k_clust.pdf
figures/_pheatmap_all_24_score.pdf
```

### Calculate Exp_RNA from Obs_RNA
To calculate change in transcription efficiency we first calculate the expected RNA abundance given the observed abundance in the ancestor, multiplied by the copy_number in the evolved strain.

```{}
make_exp_rna_from_obs_rna.py

# these rely on: 
metadata/chemostat_gene_relative_copy_number.tsv
metadata/Transposable_elements_rDNA.txt
data/analyses/chemostat_expression/RNA_unnormed_read_counts_by_gene.tsv

#this produces: 
analyses/chemostat_expression/combined_coverage_table_expected.tsv
```

### DESeq2 performed on Obs_RNA versus Exp_RNA 
To identify genes with different relative abundances DESeq2 was ran on copy-number corrected, un-normalized reads.

```{}
DESeq_difExp_RNA_v1.1.R

# input
analyses/chemostat_expression/combined_coverage_table_expected.tsv

# output
analyses/deseq/DESeq_Exp_RNA_DGY1657_DGY1726.txt
analyses/deseq/DESeq_Exp_RNA_DGY1657_DGY1735.txt
analyses/deseq/DESeq_Exp_RNA_DGY1657_DGY1741.txt
analyses/deseq/DESeq_Exp_RNA_DGY1657_DGY1743.txt
```

### Perform Fisher Exact Test on translation efficiency ratios
Calculate mediansof replicates for TPM normalized RPF, RNA in both Evolved and Ancestor. Use subsequent ratios to calculate FET and p-values.

```{}
calc_teff_v2.py

# input
analyses/efficiency/ratio_df_Unit_v13.tab

#output
analyses/efficiency/Unit_v13_teff_results.tab
```

### Make Figure 2B, 2C
Generates plots for the 'copy number correction' figures (Figure 2B, 2c). Also makes supplementary files that summarize DESeq2 results for Observed RNA, Observed RPF, and copy-number corrected Expected RNA. 
```{}
plot_fig2B_Obs_Exp.py

# input 
analyses/efficiency/ratio_df_Unit_v13.tab
analyses/efficiency/Unit_v13_teff_results.tab
metadata/chemostat_gene_relative_copy_number.tsv
analyses/chemostat_deseq/DESeq_Obs_RNA_DGY1657_DGY1726.txt
analyses/chemostat_deseq/DESeq_Obs_RNA_DGY1657_DGY1735.txt
analyses/chemostat_deseq/DESeq_Obs_RNA_DGY1657_DGY1741.txt
analyses/chemostat_deseq/DESeq_Obs_RNA_DGY1657_DGY1743.txt
analyses/chemostat_deseq/DESeq_Exp_RNA_DGY1657_DGY1726.txt
analyses/chemostat_deseq/DESeq_Exp_RNA_DGY1657_DGY1735.txt
analyses/chemostat_deseq/DESeq_Exp_RNA_DGY1657_DGY1741.txt
analyses/chemostat_deseq/DESeq_Exp_RNA_DGY1657_DGY1743.txt

# output
figures/_obs_rna_deseq_pval.pdf
figures/_exp_rna_deseq_pval.pdf
analyses/chemostat_deseq/DESeq_Obs_RNA_results.txt
analyses/chemostat_deseq/DESeq_Obs_RPF_results.txt
analyses/chemostat_deseq/DESeq_Exp_RNA_results.txt
```
### Calculate linear regression and standardized residuals

```{}
plot_unit_level_2_sres_peff.r

input:
analyses/efficiency/unit_object_level_2v13v16_Unit_10pct.tab

output:
figures/_ms_rpf_l2_peff.pdf
analyses/efficiency/unit_object_level2_ms_rpf_16_Unit_Sres_fdr.csv

### Plot translation efficiency and protein expression efficiency 

```{}
plot_teff_peff_v5.py

# input:
metadata/chemostat_gene_relative_copy_number.tsv
analyses/efficiency/Unit_v13_teff_results.tab

# output:
figures/_sig_teff_efficiency_pval.pdf
figures/_sig_peff_efficiency_pval.pdf
```

### SSD1 binding site motif analysis

```{}
SSD1_motif_analysis_v2.py

# input:
analyses/ssd1/TL_from_SGD.tsv
metadata/McManus_2018_saccharomyces_cerevisiae.gff * (not included)
SGD/S288C_reference_sequence_R64-2-1_20150113.fsa * (not included)

# output
analyses/ssd1/TL_from_SN.bed
analyses/ssd1/TL_from_SGD_SN.fa
analyses/ssd1/_SSD1_hits_in_TL_from_SGD_SN.txt
analyses/ssd1/_SSD1_in_TL_from_SGD_SN.bed
analyses/ssd1/_SSD1_random_random_control_in_TL_from_SGD_SN.bed
```
### Plot Translation efficiency and SSD1 binding motif enrichment

```{}
plot_teff_ssd1.py

# input
analyses/efficiency/unit_object_level_2v13v16_Unit_10pct.tab
_predictions_0.5.bed * [from uorfish]
analyses/efficiency/Unit_v13_teff_results.tab
analyses/ssd1/SSD1_hits_in_TL_from_SGD_SN.txt

# output
figures/teff_ssd1_efficiency_pval.pdf
figures/all_strain_teff_cnv_ssd1_efficiency_pval.pdf
```

### Plot Protein expression efficiency and protein complex enrichment

```{}
plot_peff_complexes.py

# input:
analyses/efficiency/unit_object_level2_ms_rpf_16_Unit_Sres_fdr.csv
metadata/chemostat_gene_relative_copy_number.tsv
analyses/complexes/table_of_complexes.csv

# output:
figures/_background_peff_efficiency_pval.pdf
```

### Generic stats tests
Code contains small standalone tests like FET, HGM
```{}
HGM_calc.py

```

