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
data/RNA_unnormed_read_counts_by_gene.tsv
data/RPF_unnormed_read_counts_by_gene.tsv

# these produce:
data/DESeq_Obs_RNA_DGY1657_DGY1726.txt
data/DESeq_Obs_RNA_DGY1657_DGY1735.txt
data/DESeq_Obs_RNA_DGY1657_DGY1741.txt
data/DESeq_Obs_RNA_DGY1657_DGY1743.txt
#
data/DESeq_Obs_RPF_DGY1657_DGY1726.txt
data/DESeq_Obs_RPF_DGY1657_DGY1735.txt
data/DESeq_Obs_RPF_DGY1657_DGY1741.txt
data/DESeq_Obs_RPF_DGY1657_DGY1743.txt

```
### Build "Unit" object
Certain modes of analysis require or are improved by the expression data being similarly scaled. To acheive this we first transform the data using Box-Cox Power Transform followed by a MinMax transform into the Unit range [0,1], (see Figure S5). Unit transformed data are used for Figures 2A and Supplemental Figure S3.

```{}
scale_data_v0.13_Unit_public.py

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
unit_differential_v16_public.py

# input
analyses/efficiency/ratio_df_Unit_v13.tab

# output
analyses/efficiency/unit_object_level_1v13_Unit_10pct.tab
analyses/efficiency/unit_object_level_2v13v16_Unit_10pct.tab

```



```

### Make Heatmap and cluster (Figure 2)
Takes 

```{}
build_expression_object_l1_v1.0.py

# prod

### Calculate Exp_RNA from Obs_RNA
To calculate change in transcription efficiency we first calculate the expected RNA abundance given the observed abundance in the ancestor, multiplied by the copy_number in the evolved strain.

```{}
make_exp_rna_from_obs_rna.py

# these rely on: 
data/chemostat_gene_relative_copy_number.tsv
data/Transposable_elements_rDNA.txt
data/RNA_unnormed_read_counts_by_gene.tsv

#this produces: 
data/combined_coverage_table_expected.tsv
```

### DESeq2 performed on Obs_RNA versus Exp_RNA 
To identify genes with different relative abundances DESeq2 was ran on copy-number corrected, un-normalized reads.

```{}
DESeq_difExp_RNA_v1.1.R

# this relies on: 
data/combined_coverage_table_expected.tsv

#this produces: 
data/DESeq_Exp_RNA_DGY1657_DGY1726.txt
data/DESeq_Exp_RNA_DGY1657_DGY1735.txt
data/DESeq_Exp_RNA_DGY1657_DGY1741.txt
data/DESeq_Exp_RNA_DGY1657_DGY1743.txt
```

### Make robust scale heatmap
```{}
make_multilevel_heatmaps.py

# these rely on: 
data/chemostat_gene_relative_copy_number.tsv
data/DGY1726_global_expression.tsv
data/DGY1735_global_expression.tsv
data/DGY1741_global_expression.tsv
data/DGY1743_global_expression.tsv

# these produce:
heatmap_CNV_map.pdf
Log2FC_map_GR_RNA_chrXI.pdf
Log2FC_map_GR_RPF_chrXI.pdf
Log2FC_map_GR_MS_chrXI.pdf
```

### Efficiency analysis 
One of the core metrics in the analysis is the efficiency ratios. A conceptual extension of translational efficiency (ie. RPF/RNA) to include post-translational ratios like (Protein/RPF).

```{}
rpf_tpm_and_ms_FET_v3.py

# this relies on:
data/aa_protein_lengths.tsv
data/Transposable_elements_rDNA.txt
data/DGY1726_v_DGY1765_wCON_wBatch_QN_p0.01.txt
data/DGY1735_v_DGY1765_wCON_wBatch_QN_p0.01.txt
data/DGY1741_v_DGY1765_wCON_wBatch_QN_p0.01.txt
data/DGY1743_v_DGY1765_wCON_wBatch_QN_p0.01.txt
#
data/DESeq_Obs_RNA_DGY1657_DGY1726.txt
data/DESeq_Obs_RNA_DGY1657_DGY1735.txt
data/DESeq_Obs_RNA_DGY1657_DGY1741.txt
data/DESeq_Obs_RNA_DGY1657_DGY1743.txt
#
data/DESeq_Obs_RPF_DGY1657_DGY1726.txt
data/DESeq_Obs_RPF_DGY1657_DGY1735.txt
data/DESeq_Obs_RPF_DGY1657_DGY1741.txt
data/DESeq_Obs_RPF_DGY1657_DGY1743.txt
#
data/DESeq_Exp_RNA_DGY1657_DGY1726.txt
data/DESeq_Exp_RNA_DGY1657_DGY1735.txt
data/DESeq_Exp_RNA_DGY1657_DGY1741.txt
data/DESeq_Exp_RNA_DGY1657_DGY1743.txt
#
Perseus_DA_Welchs_t-test_wFDR.txt

#this produces:
data/ms_counts_expression.tsv
data/ms_to_rna_ratio_1e2.tsv
data/ms_to_rpf_ratio_1e2.tsv
data/rpf_to_rna_ratio_1e2.tsv
#
data/dosage_compensation_FET_rna_v_rna.csv
data/dosage_compensation_FET_rpf_v_rna.csv
data/dosage_compensation_FET_ms_v_rna.csv
data/dosage_compensation_FET_ms_v_rpf.csv
#
temp_transcript_eff_cnn.csv

#Figures:
data/ms_to_rna_efficiency_pval_ratio.pdf
data/ms_to_rpf_efficiency_pval_ratio.pdf
data/rpf_to_rna_efficiency_pval_ratio.pdf
data/rna_rna_efficiency_pval_ratio.pdf
data/DGY1726_rna_v_DGY1657_rna.pdf
```

### Generic stats tests
Code contains small standalone tests like FET, HGM
```{}
generic_stats_tests.py
```

### Convert RPF to P-site fractions
Ribosome reading frames can be estimated by calculating the P-site (Ingolia et al. 2009) However, not all RPFs are the ideal 28nt fragment size, be it do to altered confirmation, over- or incomplete digestion, or enzymatic addition of nucleotides during library construction. As such for non-ideal fragment sizes a fraction of the ribosome is assigned to possible positions (Spealman and Naik et al. 2018). 
```{}
convert_to_psite.py

# this is invoked as:
strain=DGY1657
	rep=R1
		python convert_to_psite.py -i RPF_${strain}_${rep}.sorted.sam -o RPF_${strain}_${rep}

# this relies on:
RPF_${strain}_${rep}.sorted.sam STAR aligned RPFs, samtools sorted

#this produces:
RPF_${strain}_${rep}.bedgrahp
```

### SSD1 motif analysis
```{}
motif_analysis.py

# this relies on:
data/Spealman_Naik_2018_saccharomyces_cerevisiae.gff
data/TL_from_SGD.tsv
data/TL_from_SN.fa

#this produces:
data/TL_from_SN.bed
data/TL_from_SGD_SN.fa
data/SSD1_hits_in_TL_from_SGD_SN.txt
```

### Violin plot for Figure 5C 
Make violin plot of SSD1 motif Observed-Expected

```{}
violin_ssd1_uorf.py

# this relies on: 
data/hist_xeff_1e2_uORF.txt
data/hist_xeff_1e2_not_sig_uORF.txt
data/hist_xeff_1e2_no_uORF.txt
data/hist_xeff_1e2_sig.txt
data/hist_xeff_1e2.txt
data/hist_Bayne_SSD1_uORFs_targets.txt
data/hist_Bayne_SSD1_targets.txt
data/hist_all_SSD1_targets.txt
data/hist_YKL_all.txt

#Figures: 
hist_SSD1_uORF.pdf
```

