# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:52:05 2023

@author: pspea
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 12:05:39 2022

@author: pspea
"""
import pandas as pd
import numpy as np

output_figure_name = ('figures/DESeq_CNV_chemostat_map.pdf')
copy_number_filename = ('metadata/chemostat_gene_relative_copy_number.tsv')
df = pd.read_table(copy_number_filename, index_col=0)
cn_dict = df.to_dict('index')

strain_list = list(cn_dict['YKR039W'].keys())
strain_list.sort()

strain_list = ['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']
#We need to populate the gene list with those genes that are detected in every strain - 
# otherwise there will be misalignement on between genes between strains on the global heatmap

genes_to_remove_filename = ('metadata/Transposable_elements_rDNA.txt')
df = pd.read_table(genes_to_remove_filename, index_col=0)
genes_to_remove_dict = df.to_dict('index')

counts_filename = ('data/analyses/chemostat_expression/RNA_unnormed_read_counts_by_gene.tsv')
df = pd.read_table(counts_filename, index_col=0)
counts_dict = df.to_dict('index')

exp_counts_dict = {}

for gene in counts_dict:
    if (gene in cn_dict) and (gene not in genes_to_remove_dict):
        
        if gene not in exp_counts_dict:
            exp_counts_dict[gene] = {}
        
        for strain in strain_list:
            if strain in cn_dict[gene]:
                for i in range(1,3):
                    obs_col_name = ("{}_{}").format(strain, i)
                    
                    if obs_col_name in counts_dict[gene]:
                        
                        if obs_col_name not in exp_counts_dict[gene]:
                            exp_counts_dict[gene][obs_col_name] = 0
                            
                        exp_counts_dict[gene][obs_col_name] += counts_dict[gene][obs_col_name]
                        
                        for j in range(1,3):
                            anc_col_name = ("DGY1657_{}").format(j)
                            exp_col_name = ("{}_exp_{}").format(strain, j)
                            exp_count = round(counts_dict[gene][anc_col_name]*cn_dict[gene][strain])
                    
                            exp_counts_dict[gene][exp_col_name] = exp_count

exp_filename = ('data/analyses/chemostat_expression/combined_coverage_table_expected.tsv')
exp_counts_df = pd.DataFrame.from_dict(exp_counts_dict, orient='index')

exp_counts_df = exp_counts_df.reindex(sorted(exp_counts_df.columns), axis=1)

exp_counts_df.to_csv(path_or_buf=exp_filename, na_rep = np.nan)