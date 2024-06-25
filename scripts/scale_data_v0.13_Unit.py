# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 10:28:50 2024



@author: pspea
"""

import pandas as pd
import numpy as np

# from sklearn.preprocessing import QuantileTransformer
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PowerTransformer

DGY1657_sample_to_rep_dict = {'2.A':'rep1', '8.A':'rep2', '14.A':'rep3', '5.B':'rep4', '11.B': 'rep5'}
DGY1726_sample_to_rep_dict = {'3.A':'rep1', '9.A':'rep2', '15.A':'rep3', '6.B':'rep4', '12.B': 'rep5'}
DGY1735_sample_to_rep_dict = {'5.A':'rep1', '11.A':'rep2', '2.B':'rep3', '8.B':'rep4', '14.B': 'rep5'}
DGY1741_sample_to_rep_dict = {'6.A':'rep1', '12.A':'rep2', '3.B':'rep3', '9.B':'rep4', '15.B': 'rep5'}
DGY1743_sample_to_rep_dict = {'7.A':'rep1', '13.A':'rep2', '4.B':'rep3', '10.B':'rep4', '16.B': 'rep5'}

sample_to_rep_dict = {'DGY1657': DGY1657_sample_to_rep_dict,
                      'DGY1726': DGY1726_sample_to_rep_dict,
                      'DGY1735': DGY1735_sample_to_rep_dict,
                      'DGY1741': DGY1741_sample_to_rep_dict,
                      'DGY1743': DGY1743_sample_to_rep_dict}


def record_sample_ms(gene, strain, replicate, value, observed_ms_dict):

    obs_column_name = ('ms_{strain}_{replicate}').format(
        strain = strain, replicate = replicate) 
            
    if gene not in observed_ms_dict:
        observed_ms_dict[gene] = {}
    
    if strain not in observed_ms_dict[gene]:
        observed_ms_dict[gene][obs_column_name] = value
    else:
        print('error')
        1/0

    return(observed_ms_dict)

def calc_median_ms(gene, strain, obs_list, observed_ms_dict):

    obs_column_name = ('ms_{}_median').format(strain) 
            
    if gene not in observed_ms_dict:
        observed_ms_dict[gene] = {}
    
    if strain not in observed_ms_dict[gene]:
        observed_ms_dict[gene][obs_column_name] = np.median(obs_list)
    else:
        print('error')
        1/0

    return(observed_ms_dict)

def seperate_rep_and_median(strain, index, gene, ms_dict, observed_ms_dict):
    obs_list = []
    
    strain_sample_to_rep = sample_to_rep_dict[strain]
    
    
    for sep_sample in strain_sample_to_rep:
        replicate = strain_sample_to_rep[sep_sample]
        
        col_name = ('Reporter.intensity.corrected.{}').format(sep_sample)
        value = ms_dict[index][col_name]
        obs_list.append(value)
        
        observed_ms_dict = record_sample_ms(gene, strain, replicate, value, observed_ms_dict)
        
    observed_ms_dict = calc_median_ms(gene, strain, obs_list, observed_ms_dict)        
    
    return(observed_ms_dict)

def make_observed_ms():
    # zscr_scalar = zscoreScaler()
    observed_ms_dict = {}
    
    genes_to_remove_filename = ('metadata/Transposable_elements_rDNA.txt')
    genes_to_remove_df = pd.read_table(genes_to_remove_filename, index_col=0)
    genes_to_remove_dict = genes_to_remove_df.to_dict('index')
        
    for strain in ['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']:

        infile = ('analyses/ms/{}_v_DGY1657_wCON_wBatch_QN_p0.01.txt').format(strain)
        df = pd.read_table(infile)
        ms_dict = df.to_dict('index')
        
        for index in ms_dict:
            gene = ms_dict[index]['Majority.protein.IDs']
            
            if gene[0] == 'Y':
                if ';' in gene:
                    gene_list = gene.split(';')
            
                    for gene in gene_list:
                        if gene not in genes_to_remove_dict:
                            observed_ms_dict = seperate_rep_and_median(strain, index, gene, ms_dict, observed_ms_dict)
                else:
                    if gene not in genes_to_remove_dict:
                        observed_ms_dict = seperate_rep_and_median(strain, index, gene, ms_dict, observed_ms_dict)
                        
    for strain in set(['DGY1657']):
        infile = ('analyses/ms/DGY1735_v_{}_wCON_wBatch_QN_p0.01.txt').format(strain)
        df = pd.read_table(infile)
        ms_dict = df.to_dict('index')
        
        for index in ms_dict:
            gene = ms_dict[index]['Majority.protein.IDs']
            
            if gene[0] == 'Y':
                if ';' in gene:
                    gene_list = gene.split(';')
            
                    for gene in gene_list:
                        if gene not in genes_to_remove_dict:
                            observed_ms_dict = seperate_rep_and_median(strain, index, gene, ms_dict, observed_ms_dict)
                else:
                    if gene not in genes_to_remove_dict:
                        observed_ms_dict = seperate_rep_and_median(strain, index, gene, ms_dict, observed_ms_dict)

    outfile_name = ('analyses/efficiency/ms_counts_expression.tsv')
    df = pd.DataFrame.from_dict(observed_ms_dict, orient='index')
    df = df.fillna(0)

    df.to_csv(outfile_name, sep = '\t')
    
    return(df)

def make_tpm(df, sample_name, strain):
    """
    convert read counts to TPM (transcripts per million)
    :param df: a dataFrame contains the result coming from featureCounts
    :param sample_name: a list, all sample names, same as the result of featureCounts
    :return: TPM
    """
    #result = df
    
    genes_to_remove_filename = ('metadata/Transposable_elements_rDNA.txt')
    genes_to_remove_df = pd.read_table(genes_to_remove_filename, index_col=0)
    genes_to_remove_dict = genes_to_remove_df.to_dict('index')
    
    result_dict = df.to_dict('index')
    
    for rem_gene in genes_to_remove_dict:
        if rem_gene in result_dict:
            del result_dict[rem_gene]
    
    total_rate = 0
        
    for gene in result_dict:
        if sample_name in result_dict[gene]:
            nt_length = result_dict[gene]['nt_length']
            counts = result_dict[gene][sample_name]
            total_rate += counts / nt_length
                
    tpm_name = ('{}_tpm').format(sample_name)
    
    for gene in result_dict:
        if sample_name in result_dict[gene]:
            nt_length = result_dict[gene]['nt_length']
            counts = result_dict[gene][sample_name]
            rate = counts / nt_length
            result_dict[gene][tpm_name] = 1e6*(rate / total_rate)
            
    df = pd.DataFrame.from_dict(result_dict, orient='index')
    
    return(df)

def make_PT_transform(istype, strain):
    global ratio_df
    
    pt = PowerTransformer()
    mms = MinMaxScaler()
    
    rep_set = set(['R1', 'R2'])
        
    for rep in rep_set:
        rep_raw = ('CDS_{istype}_{strain}_{rep}').format(istype = istype, strain = strain, rep = rep)
        rep_tpm = ('{rep_raw}_tpm').format(rep_raw = rep_raw)
        
        ratio_df = make_tpm(ratio_df, rep_raw, strain)
              
        rep_upt_scale = ('{rep_tpm}_upt').format(rep_tpm = rep_tpm)
        ratio_df[rep_upt_scale] = (pt.fit_transform(ratio_df[[rep_tpm]]))
        
        rep_pt_scale = ('{rep_tpm}_pt').format(rep_tpm = rep_tpm)
        ratio_df[rep_pt_scale] = mms.fit_transform((pt.fit_transform(ratio_df[[rep_tpm]])))
            
        # correction = abs(min(ratio_df[rep_upt_scale]))
        # rep_pt_scale = ('{rep_tpm}_pt').format(rep_tpm = rep_tpm)
        # ratio_df[rep_pt_scale] = ratio_df[rep_upt_scale] + correction
        
        
    median_colname = ('CDS_{istype}_{strain}_median_tpm').format(istype = istype, strain = strain)
    rep_1_tpm_colname = ('CDS_{istype}_{strain}_R1_tpm').format(istype = istype, strain = strain)
    rep_2_tpm_colname = ('CDS_{istype}_{strain}_R2_tpm').format(istype = istype, strain = strain)
    ratio_df[median_colname] = ratio_df[[rep_1_tpm_colname, rep_2_tpm_colname]].median()    

def make_ms_PT_transform(strain):
    global ratio_df
    
    pt = PowerTransformer()
    mms = MinMaxScaler()
    
    rep_set = set(['rep1', 'rep2', 'rep3', 'rep4', 'rep5'])
    
    # name_set = set()
    
    for rep in rep_set:
        rep_raw = ('ms_{strain}_{rep}').format(strain = strain, rep = rep)
        #rep_tpm = ('{rep_raw}_tpm').format(rep_raw = rep_raw)
        
        #ratio_df = make_tpm(ratio_df, rep_raw, strain)
        
        rep_upt_scale = ('{rep_raw}_upt').format(rep_raw = rep_raw)
        ratio_df[rep_upt_scale] = (pt.fit_transform(ratio_df[[rep_raw]])) 

        rep_pt_scale = ('{rep_raw}_pt').format(rep_raw = rep_raw)
        ratio_df[rep_pt_scale] = mms.fit_transform((pt.fit_transform(ratio_df[[rep_raw]])))
        
        # correction = abs(min(ratio_df[rep_upt_scale]))
        
        # rep_pt_scale = ('{rep_raw}_pt').format(rep_raw = rep_raw)
        # ratio_df[rep_pt_scale] = ratio_df[rep_upt_scale] + correction
        
        #rep_tpm_scale = ('{rep_tpm}_robust').format(rep_tpm = rep_tpm)
        #ratio_df[rep_tpm_scale] = (pt.fit_transform(ratio_df[[rep_tpm]]))
        
        # name_set.add(rep_raw)
        # #name_set.add(rep_tpm)
        # name_set.add(rep_raw_scale)
        # #name_set.add(rep_tpm_scale)
      



#import copy_number:
copy_number_filename = ('metadata/chemostat_gene_relative_copy_number.tsv')
cn_df = pd.read_table(copy_number_filename, index_col=0)


infile_name = ('metadata/aa_protein_lengths.tsv')

length_df = pd.read_table(infile_name, index_col=1)
length_df['nt_length'] = length_df["aa_length"]*3
#sizes_dict = df.to_dict('index')

#import read abundances:
infile_name = ('analyses/chemostat_expression/expression_dict.tsv') 

counts_df = pd.read_table(infile_name, index_col=0)

to_transfer_df = pd.merge(left = counts_df,
                           right = length_df,
                           left_index=True,
                           right_index=True)

to_transform_df = pd.merge(left = to_transfer_df,
                           right = cn_df,
                           left_index=True,
                           right_index=True)
            
to_transform_df = to_transform_df[to_transform_df["status"] == "Verified"]

expected_ms_df = make_observed_ms()

ratio_df = pd.merge(left = to_transform_df,
                    right = expected_ms_df,
                    left_index = True,
                    right_index = True)
#
for strain in set(['DGY1657', 'DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
    for istype in set(['RNA', 'RPF']):
        make_PT_transform(istype, strain)
    
    make_ms_PT_transform(strain)


outfile_name = ('analyses/efficiency/ratio_df_Unit_v13.tab')    
#ratio_df = pd.DataFrame.from_dict(ratio_dict, orient='index')
ratio_df.to_csv(outfile_name, sep = '\t') 