# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 10:41:37 2024

@author: pspea
"""

import random
import plotly.io as pio
pio.renderers.default = "browser"

#from scipy.stats import fisher_exact
import pandas as pd
import numpy as np

#
#from scipy.stats import mannwhitneyu
import plotly.graph_objects as go

from scipy.stats import fisher_exact

if False:
    infile_name = ('analyses/efficiency/ratio_df_Unit_v13.tab')    
    ratio_df = pd.read_table(infile_name, index_col=0)
    gene_set = set(ratio_df.index)
        
    infile = ('analyses/efficiency/Unit_v13_teff_results.tab')
    df = pd.read_table(infile, index_col=0)
    
    copy_number_filename = ('metadata/chemostat_gene_relative_copy_number.tsv')
    cn_df = pd.read_table(copy_number_filename, index_col=0)
    cn_df = cn_df.add_suffix('_copy_number')
    
    df = pd.merge(left = df,
                  right = cn_df,
                  left_index=True,
                  right_index=True)
    
    #rna_dict = df.to_dict('index')
        
    
       
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        
        
        suffix_string = ('_deseq_{}').format(strain)
        
        
        deseq_file_name = ("analyses/chemostat_deseq/DESeq_Obs_RNA_DGY1657_{}.txt").format(strain)
        deseq_df = pd.read_table(deseq_file_name, index_col=0)
        deseq_df = deseq_df.add_suffix(suffix_string)
        
        df = pd.merge(left = df,
                      right = deseq_df,
                      left_index=True,
                      right_index=True)   
        
    output_figure_name = ('C:/Gresham/Project_Carolino/figures/'
                          '/_obs_rna_deseq_pval.pdf')    
    fig = go.Figure()
        
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        anc_col_name = ('CDS_RNA_DGY1657_median_tpm')
        evo_col_name = ('CDS_RNA_{}_median_tpm').format(strain)
        
        cn_col_name = ('{}_copy_number').format(strain)
        padj_col_name = ('padj_deseq_{}').format(strain) 

        background_df = df[ (df[padj_col_name] <= 0.05) & (df[cn_col_name] == 1)]
        index_vals = background_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(background_df[anc_col_name]),
                y=np.log2(background_df[evo_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 100,
                    color = 'Grey', #'Blue',
                    opacity = 0.3,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        sig_2_df = df[ (df[padj_col_name] <= 0.05) & (df[cn_col_name] == 2)]
        index_vals = sig_2_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(sig_2_df[anc_col_name]),
                y=np.log2(sig_2_df[evo_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 0,
                    color = 'Blue', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        insig_2_df = df[ (df[padj_col_name] > 0.05) & (df[cn_col_name] == 2)]
        index_vals = insig_2_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(insig_2_df[anc_col_name]),
                y=np.log2(insig_2_df[evo_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 1,
                    color = 'Purple', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        sig_3_df = df[ (df[padj_col_name] <= 0.05) & (df[cn_col_name] == 3)]
        index_vals = sig_3_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(sig_3_df[anc_col_name]),
                y=np.log2(sig_3_df[evo_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 0,
                    color = 'Orange', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        insig_3_df = df[ (df[padj_col_name] > 0.05) & (df[cn_col_name] == 3)]
        index_vals = insig_3_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(insig_2_df[anc_col_name]),
                y=np.log2(insig_2_df[evo_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 1,
                    color = 'Yellow', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        sig_4_df = df[ (df[padj_col_name] <= 0.05) & (df[cn_col_name] == 4)]
        index_vals = sig_4_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(sig_4_df[anc_col_name]),
                y=np.log2(sig_4_df[evo_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 0,
                    color = 'Red', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        insig_4_df = df[ (df[padj_col_name] > 0.05) & (df[cn_col_name] == 4)]
        index_vals = insig_4_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(insig_4_df[anc_col_name]),
                y=np.log2(insig_4_df[evo_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 1,
                    color = 'Green', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )

    plot_title = ('Evolved versus Ancestor RNA')
    xaxis_title = ("Ancestor strain log2(RNA)")
    yaxis_title = ("Evolved strain log2(RNA)")
    
    fig.update_layout(
        title=plot_title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        font=dict(
            size=18,
            color="Black"
        ),
        autosize=False,
        width=800,
        height=800,
    
    )
    
    fig.update_traces(marker=dict(size=12,
                              line=dict(width=2,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))

    # fig.update_xaxes(range=[2.5, 12.5])
    # fig.update_yaxes(range=[2.5, 12.5])
    
    # fig.update_xaxes(range=[min_axis, max_axis])
    # fig.update_yaxes(range=[min_axis, max_axis])
    
    fig.show()
    fig.write_image(output_figure_name)
                
#
if True:
    infile_name = ('analyses/efficiency/ratio_df_Unit_v13.tab')    
    ratio_df = pd.read_table(infile_name, index_col=0)
    gene_set = set(ratio_df.index)
        
    infile = ('analyses/efficiency/Unit_v13_teff_results.tab')
    df = pd.read_table(infile, index_col=0)
    
    copy_number_filename = ('metadata/chemostat_gene_relative_copy_number.tsv')
    cn_df = pd.read_table(copy_number_filename, index_col=0)
    cn_df = cn_df.add_suffix('_copy_number')
    
    df = pd.merge(left = df,
                  right = cn_df,
                  left_index=True,
                  right_index=True)
    
    #rna_dict = df.to_dict('index')
               
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        
        
        suffix_string = ('_deseq_{}').format(strain)
        
        
        deseq_file_name = ("analyses/chemostat_deseq/DESeq_Exp_RNA_DGY1657_{}.txt").format(strain)
        deseq_df = pd.read_table(deseq_file_name, index_col=0)
        deseq_df = deseq_df.add_suffix(suffix_string)
        
        df = pd.merge(left = df,
                      right = deseq_df,
                      left_index=True,
                      right_index=True)   
        

    output_figure_name = ('C:/Gresham/Project_Carolino/figures/'
                          '/_exp_rna_deseq_pval.pdf')    
    fig = go.Figure()
    
    
        
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        anc_col_name = ('CDS_RNA_DGY1657_median_tpm')
        evo_col_name = ('CDS_RNA_{}_median_tpm').format(strain)
        
        exp_col_name = ('Exp_RNA_{}_median_tpm').format(strain)
        cn_col_name = ('{}_copy_number').format(strain)
        
        df[exp_col_name] = df[evo_col_name]/df[cn_col_name]
            
        padj_col_name = ('padj_deseq_{}').format(strain) 
        
        background_df = df[ (df[padj_col_name] <= 0.05) & (df[cn_col_name] == 1)]
        index_vals = background_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(background_df[anc_col_name]),
                y=np.log2(background_df[exp_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 100,
                    color = 'Grey', #'Blue',
                    opacity = 0.3,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        sig_2_df = df[ (df[padj_col_name] <= 0.05) & (df[cn_col_name] == 2)]
        index_vals = sig_2_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(sig_2_df[anc_col_name]),
                y=np.log2(sig_2_df[exp_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 0,
                    color = 'Blue', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        insig_2_df = df[ (df[padj_col_name] > 0.05) & (df[cn_col_name] == 2)]
        index_vals = insig_2_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(insig_2_df[anc_col_name]),
                y=np.log2(insig_2_df[exp_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 1,
                    color = 'Purple', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        sig_3_df = df[ (df[padj_col_name] <= 0.05) & (df[cn_col_name] == 3)]
        index_vals = sig_3_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(sig_3_df[anc_col_name]),
                y=np.log2(sig_3_df[exp_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 0,
                    color = 'Orange', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        insig_3_df = df[ (df[padj_col_name] > 0.05) & (df[cn_col_name] == 3)]
        index_vals = insig_3_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(insig_3_df[anc_col_name]),
                y=np.log2(insig_3_df[exp_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 1,
                    color = 'Yellow', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        sig_4_df = df[ (df[padj_col_name] <= 0.05) & (df[cn_col_name] == 4)]
        index_vals = sig_4_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(sig_4_df[anc_col_name]),
                y=np.log2(sig_4_df[exp_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 0,
                    color = 'Red', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )
        
        insig_4_df = df[ (df[padj_col_name] > 0.05) & (df[cn_col_name] == 4)]
        index_vals = insig_4_df[cn_col_name].astype('category').cat.codes
        
        # many teff
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=np.log2(insig_4_df[anc_col_name]),
                y=np.log2(insig_4_df[exp_col_name]),
                
                text=index_vals.index + ' ' + strain,
                marker=dict(
                    symbol = 1,
                    color = 'Green', #'Blue',
                    opacity = 0.5,
                    line=dict(
                        color='Black',
                        width=1
                    )
                ),
                showlegend=False
            )
        )

    plot_title = ('Evolved versus Ancestor RNA')
    xaxis_title = ("Ancestor strain log2(RNA)")
    yaxis_title = ("Evolved strain log2(RNA)")
    
    fig.update_layout(
        title=plot_title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        font=dict(
            size=18,
            color="Black"
        ),
        autosize=False,
        width=800,
        height=800,
    
    )
    
    fig.update_traces(marker=dict(size=12,
                              line=dict(width=2,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))

    # fig.update_xaxes(range=[2.5, 12.5])
    # fig.update_yaxes(range=[2.5, 12.5])
    
    # fig.update_xaxes(range=[min_axis, max_axis])
    # fig.update_yaxes(range=[min_axis, max_axis])
    
    #fig.show()
    #fig.write_image(output_figure_name)
   
    
cnn_set = set()             
cnv_set = set()
is_diff_cnn_set = set()
is_diff_cnv_set = set()
    
gene_counter = {}
for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        
    anc_col_name = ('CDS_RNA_DGY1657_median_tpm')
    evo_col_name = ('CDS_RNA_{}_median_tpm').format(strain)
    
    exp_col_name = ('Exp_RNA_{}_median_tpm').format(strain)
    cn_col_name = ('{}_copy_number').format(strain)
    
    df[exp_col_name] = df[evo_col_name]/df[cn_col_name]
        
    padj_col_name = ('padj_deseq_{}').format(strain) 
    
    difference_col_name = ('{}_RNA_median_difference').format(strain)
    df[difference_col_name] = (df[exp_col_name] - df[anc_col_name])
    
    subset_df = df[(df[padj_col_name] > 0.05) & 
                   (df[cn_col_name] == 1)]
    cnn_set = cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] > 0.05) &
                   (df[cn_col_name] != 1)]
    cnv_set = cnv_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] == 1)]
    is_diff_cnn_set = is_diff_cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] != 1)]
    is_diff_cnv_set = is_diff_cnv_set.union(set(subset_df.index))
            
print(len(cnn_set), len(is_diff_cnn_set))
print(len(cnv_set), len(is_diff_cnv_set))
print('cnn pct', round(100*len(is_diff_cnn_set)/len(cnn_set),3))
print('cnv pct', round(100*len(is_diff_cnv_set)/len(cnv_set),3))
    
odds, pval = fisher_exact([[len(is_diff_cnv_set), 
                            len(cnv_set)],
                           [len(is_diff_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('diff', odds, pval)
    

    
     
    
odds, pval = fisher_exact([[(11+9), 
                            (69+84)],
                            [(440+257), 
                            (4222+4207)]], alternative='two-sided')
print('g150', odds, pval)

odds, pval = fisher_exact([[(11+9+4+19), 
                            (69+84+17+81)],
                            [(440+257+1112+916), 
                            (4222+4207+4274+4210)]], alternative='two-sided')
print('all', odds, pval)

gene_counter = {}
for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
    
    cnn_set = set()             
    cnv_set = set()
    is_high_cnn_set = set()
    is_low_cnn_set = set()
    is_high_cnv_set = set()
    is_low_cnv_set = set()
    
    anc_col_name = ('CDS_RNA_DGY1657_median_tpm')
    evo_col_name = ('CDS_RNA_{}_median_tpm').format(strain)
    
    exp_col_name = ('Exp_RNA_{}_median_tpm').format(strain)
    cn_col_name = ('{}_copy_number').format(strain)
    
    df[exp_col_name] = df[evo_col_name]/df[cn_col_name]
        
    padj_col_name = ('padj_deseq_{}').format(strain) 
    
    difference_col_name = ('{}_RNA_median_difference').format(strain)
    df[difference_col_name] = (df[exp_col_name] - df[anc_col_name])
    
    subset_df = df[(df[padj_col_name] > 0.05) & 
                   (df[cn_col_name] == 1)]
    cnn_set = cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] > 0.05) &
                   (df[cn_col_name] != 1)]
    cnv_set = cnv_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] == 1) &
                   (df[difference_col_name] > 0) ]
    is_high_cnn_set = is_high_cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] == 1) &
                   (df[difference_col_name] < 0) ]
    is_low_cnn_set = is_low_cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] != 1) &
                   (df[difference_col_name] > 0) ]
    is_high_cnv_set = is_high_cnv_set.union(set(subset_df.index))
    
    for gene in set(subset_df.index):
        if gene not in gene_counter:
            gene_counter[gene] = {'high':set(), 'low':set()}
            
        gene_counter[gene]['high'].add(strain)
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] != 1) &
                   (df[difference_col_name] < 0) ]
    is_low_cnv_set = is_low_cnv_set.union(set(subset_df.index))
    
    for gene in set(subset_df.index):
        if gene not in gene_counter:
            gene_counter[gene] = {'high':set(), 'low':set()}
            
        gene_counter[gene]['low'].add(strain)
    
    print(strain)
    # cnn_total = len(cnn_set) + len(is_high_cnn_set) + len(is_low_cnn_set)
    # print(len(cnn_set), len(is_high_cnn_set), len(is_low_cnn_set))
    # print('cnn high pct', round(100*len(is_high_cnn_set)/(cnn_total),3))
    # print('cnn low pct', round(100*len(is_low_cnn_set)/(cnn_total),3))
    # print('cnn both pct', round(100*len(is_high_cnn_set)/(cnn_total),3)+round(100*len(is_low_cnn_set)/(cnn_total),3))
    
    # cnv_total = len(cnv_set) + len(is_high_cnv_set) + len(is_low_cnv_set)
    # print(len(cnv_set), len(is_high_cnv_set), len(is_low_cnv_set))
    # print('cnv high pct', round(100*len(is_high_cnv_set)/(cnv_total),3))
    # print('cnv low pct', round(100*len(is_low_cnv_set)/(cnv_total),3))
    # print('cnv both pct', round(100*len(is_high_cnv_set)/(cnv_total),3)+round(100*len(is_low_cnv_set)/(cnv_total),3))
    
    
    odds, pval = fisher_exact([[len(is_high_cnv_set), 
                                len(cnv_set)],
                               [len(is_high_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('high', odds, pval)
    
    odds, pval = fisher_exact([[len(is_low_cnv_set), 
                                len(cnv_set)],
                               [len(is_low_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('low', odds, pval)
    
    odds, pval = fisher_exact([[len(is_high_cnv_set)+len(is_low_cnv_set), 
                                len(cnv_set)],
                               [len(is_high_cnn_set)+len(is_low_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('both', odds, pval)
    

    
    # median pct teff for CNNs
    np.median([47.286, 14.377, 35.2729, 8.32])
    
    # median pct teff for CNVs
    np.median([52.941, 21.739, 41.976, 13.095])
    
    # median pct teff for CNNs
    np.median([26.018, 10.422, 21.758, 6.109])
    
    # median pct teff for CNVs
    np.median([23.529, 15.942, 23.457, 10.714])
    
cnn_set = set()             
cnv_set = set()
is_high_cnn_set = set()
is_low_cnn_set = set()
is_high_cnv_set = set()
is_low_cnv_set = set()

gene_counter = {}
for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']): 
#for strain in set(['DGY1726', 'DGY1735']): 
#for strain in set(['DGY1741', 'DGY1743']): 
    anc_col_name = ('CDS_RNA_DGY1657_median_tpm')
    evo_col_name = ('CDS_RNA_{}_median_tpm').format(strain)
    
    exp_col_name = ('Exp_RNA_{}_median_tpm').format(strain)
    cn_col_name = ('{}_copy_number').format(strain)
    
    df[exp_col_name] = df[evo_col_name]/df[cn_col_name]
        
    padj_col_name = ('padj_deseq_{}').format(strain) 
    
    difference_col_name = ('{}_RNA_median_difference').format(strain)
    df[difference_col_name] = (df[exp_col_name] - df[anc_col_name])
    
    subset_df = df[(df[padj_col_name] > 0.05) & 
                   (df[cn_col_name] == 1)]
    cnn_set = cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] > 0.05) &
                   (df[cn_col_name] != 1)]
    cnv_set = cnv_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] == 1) &
                   (df[difference_col_name] > 0) ]
    is_high_cnn_set = is_high_cnn_set.union(set(subset_df.index))
    
    # for gene in set(subset_df.index):
    #     if gene not in gene_counter:
    #         gene_counter[gene] = {'high':set(), 'low':set()}
            
    #     gene_counter[gene]['high'].add(strain)
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] == 1) &
                   (df[difference_col_name] < 0) ]
    is_low_cnn_set = is_low_cnn_set.union(set(subset_df.index))
    
    # for gene in set(subset_df.index):
    #     if gene not in gene_counter:
    #         gene_counter[gene] = {'high':set(), 'low':set()}
            
    #     gene_counter[gene]['low'].add(strain)
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] != 1) &
                   (df[difference_col_name] > 0) ]
    is_high_cnv_set = is_high_cnv_set.union(set(subset_df.index))
    
    # for gene in set(subset_df.index):
    #     if gene not in gene_counter:
    #         gene_counter[gene] = {'high':set(), 'low':set()}
            
    #     gene_counter[gene]['high'].add(strain)
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] != 1) &
                   (df[difference_col_name] < 0) ]
    is_low_cnv_set = is_low_cnv_set.union(set(subset_df.index))
    
    # for gene in set(subset_df.index):
    #     if gene not in gene_counter:
    #         gene_counter[gene] = {'high':set(), 'low':set()}
            
    #     gene_counter[gene]['low'].add(strain)

ct = 0  
print_gene_set = set() 
for gene in gene_counter:
    
    
    hsize = len(gene_counter[gene]['high'])
    
    lsize = len(gene_counter[gene]['low'])

    # if (hsize > 1):
    #     print(gene, gene_counter[gene])
        
    # if (lsize > 1):
    #     ct+=1
    #     print(gene, gene_counter[gene])

    if (hsize > 1) and (lsize == 0):
        ct+=1
        print(gene, gene_counter[gene])
        print_gene_set.add(gene)
            
    # if (hsize == 0) and (lsize > 0):
    #     ct+=1
    #     print(gene, gene_counter[gene])
        
    
print(ct)
print(print_gene_set)
    
    
    # print(strain)
    # cnn_total = len(cnn_set) + len(is_high_cnn_set) + len(is_low_cnn_set)
    # print(len(cnn_set), len(is_high_cnn_set), len(is_low_cnn_set))
    # print('cnn high pct', round(100*len(is_high_cnn_set)/(cnn_total),3))
    # print('cnn low pct', round(100*len(is_low_cnn_set)/(cnn_total),3))
    # print('cnn both pct', round(100*len(is_high_cnn_set)/(cnn_total),3)+round(100*len(is_low_cnn_set)/(cnn_total),3))
    
    # cnv_total = len(cnv_set) + len(is_high_cnv_set) + len(is_low_cnv_set)
    # print(len(cnv_set), len(is_high_cnv_set), len(is_low_cnv_set))
    # print('cnv high pct', round(100*len(is_high_cnv_set)/(cnv_total),3))
    # print('cnv low pct', round(100*len(is_low_cnv_set)/(cnv_total),3))
    # print('cnv both pct', round(100*len(is_high_cnv_set)/(cnv_total),3)+round(100*len(is_low_cnv_set)/(cnv_total),3))
    
odds, pval = fisher_exact([[len(is_high_cnv_set), 
                            len(cnv_set)],
                           [len(is_high_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('high', odds, pval)
#high 1.000167881005943 1.0

odds, pval = fisher_exact([[len(is_low_cnv_set), 
                            len(cnv_set)],
                           [len(is_low_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('low', odds, pval)
#low 0.9162102907873725 0.8206495839743074

#restricted to DGY1726, DGY1735
odds, pval = fisher_exact([[len(is_high_cnv_set), 
                            len(cnv_set)],
                           [len(is_low_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('g150', odds, pval)
#high 1.3088332234673699 0.5938002275502852
        
#restricted to DGY1726, DGY1735
odds, pval = fisher_exact([[len(is_low_cnv_set), 
                            len(cnv_set)],
                           [len(is_high_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('g150', odds, pval)
#g150 1.6968272620446534 0.06456671086193456

#restricted to DGY1741, DGY1743
odds, pval = fisher_exact([[len(is_high_cnv_set), 
                            len(cnv_set)],
                           [len(is_low_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('g250', odds, pval)
#g250 0.8062650558928878 0.4215946198504915
        
#restricted to DGY1726, DGY1735
odds, pval = fisher_exact([[len(is_low_cnv_set), 
                            len(cnv_set)],
                           [len(is_high_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('g250', odds, pval)
#g250 1.1365544102340261 0.5508981358531247
    
#
cnn_set = set()             
cnv_set = set()
is_high_cnn_set = set()
is_low_cnn_set = set()
is_high_cnv_set = set()
is_low_cnv_set = set()

for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
    
    anc_col_name = ('CDS_RNA_DGY1657_median_tpm')
    evo_col_name = ('CDS_RNA_{}_median_tpm').format(strain)
    
    exp_col_name = ('Exp_RNA_{}_median_tpm').format(strain)
    cn_col_name = ('{}_copy_number').format(strain)
    
    df[exp_col_name] = df[evo_col_name]/df[cn_col_name]
        
    padj_col_name = ('padj_deseq_{}').format(strain) 
    
    difference_col_name = ('{}_RNA_median_difference').format(strain)
    df[difference_col_name] = (df[exp_col_name] - df[anc_col_name])
    
    subset_df = df[(df[cn_col_name] == 1)]
    cnn_set = cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[cn_col_name] != 1)]
    cnv_set = cnv_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] == 1) &
                   (df[difference_col_name] > 0) ]
    is_high_cnn_set = is_high_cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] == 1) &
                   (df[difference_col_name] < 0) ]
    is_low_cnn_set = is_low_cnn_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] != 1) &
                   (df[difference_col_name] > 0) ]
    is_high_cnv_set = is_high_cnv_set.union(set(subset_df.index))
    
    subset_df = df[(df[padj_col_name] <= 0.05) & 
                   (df[cn_col_name] != 1) &
                   (df[difference_col_name] < 0) ]
    is_low_cnv_set = is_low_cnv_set.union(set(subset_df.index))
    
    
print(len(cnn_set), len(is_high_cnn_set), len(is_low_cnn_set))
print(len(cnv_set), len(is_high_cnv_set), len(is_low_cnv_set))

odds, pval = fisher_exact([[len(is_high_cnv_set), 
                            len(cnv_set)],
                           [len(is_high_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('high', odds, pval)

odds, pval = fisher_exact([[len(is_low_cnv_set), 
                            len(cnv_set)],
                           [len(is_low_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('low', odds, pval)

###
# Output supplemental
###
deseq_file_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Obs_RNA_DGY1657_DGY1726.txt")
deseq_df = pd.read_table(deseq_file_name, index_col=0)
suffix_string = ('_deseq_obs_RNA_DGY1726')
df = deseq_df.add_suffix(suffix_string)


for strain in set(['DGY1735', 'DGY1741', 'DGY1743']):
        
        
        suffix_string = ('_deseq_obs_RNA_{}').format(strain)
        
        
        deseq_file_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Obs_RNA_DGY1657_{}.txt").format(strain)
        deseq_df = pd.read_table(deseq_file_name, index_col=0)
        deseq_df = deseq_df.add_suffix(suffix_string)
        
        df = pd.merge(left = df,
                      right = deseq_df,
                      left_index=True,
                      right_index=True) 
        
outfile_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Obs_RNA_results.txt")
df.to_csv(outfile_name, sep = '\t') 
#
deseq_file_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Obs_RPF_DGY1657_DGY1726.txt")
deseq_df = pd.read_table(deseq_file_name, index_col=0)
suffix_string = ('_deseq_obs_RPF_DGY1726')
df = deseq_df.add_suffix(suffix_string)


for strain in set(['DGY1735', 'DGY1741', 'DGY1743']):
        
        
        suffix_string = ('_deseq_obs_RPF_{}').format(strain)
        
        
        deseq_file_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Obs_RPF_DGY1657_{}.txt").format(strain)
        deseq_df = pd.read_table(deseq_file_name, index_col=0)
        deseq_df = deseq_df.add_suffix(suffix_string)
        
        df = pd.merge(left = df,
                      right = deseq_df,
                      left_index=True,
                      right_index=True) 
        
outfile_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Obs_RPF_results.txt")
df.to_csv(outfile_name, sep = '\t')
#
deseq_file_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Exp_RNA_DGY1657_DGY1726.txt")
deseq_df = pd.read_table(deseq_file_name, index_col=0)
suffix_string = ('_deseq_exp_RNA_DGY1726')
df = deseq_df.add_suffix(suffix_string)


for strain in set(['DGY1735', 'DGY1741', 'DGY1743']):
        
        
        suffix_string = ('_deseq_exp_RNA_{}').format(strain)
        
        
        deseq_file_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Exp_RNA_DGY1657_{}.txt").format(strain)
        deseq_df = pd.read_table(deseq_file_name, index_col=0)
        deseq_df = deseq_df.add_suffix(suffix_string)
        
        df = pd.merge(left = df,
                      right = deseq_df,
                      left_index=True,
                      right_index=True) 
        
outfile_name = ("C:/Gresham/Project_Carolino/analyses/chemostat_deseq/DESeq_Exp_RNA_results.txt")
df.to_csv(outfile_name, sep = '\t')