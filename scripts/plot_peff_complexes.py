# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 12:59:23 2024

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

def two_way(total, left, right, obs, sim=1000):
    random_bag = range(total)
    sig_fig = len(str(int(sim)))-1

    by_chance_list = []
    exp_list = []

    for i in range(int(sim)):
        
        left_set = set(random.sample(random_bag, left))
        
        right_set = set(random.sample(random_bag, right))
        
        intersect = left_set.intersection(right_set)
        
        exp = len(intersect)
        
        if exp >= obs:
            by_chance_list.append(1)
        else:
            by_chance_list.append(0)
            
        exp_list.append(exp)
    
    exp_median = np.median(exp_list)    
    
    if exp_median != 0:
        ratio = round(obs/exp_median, 3) 
    else:
        ratio = obs
        
    outline = ('Median exp: {exp}, Num obs: {obs}, ratio: {ratio}, pval:{pval}\n').format(
        exp = exp_median,
        obs = obs,
        ratio = ratio,
        pval = round(sum(by_chance_list)/len(by_chance_list), sig_fig)
        )
    
    print(outline)
    
def two_way_deenrichment(total, left, right, obs, sim=10000):
    random_bag = range(total)
    sig_fig = len(str(int(sim)))-1

    by_chance_list = []
    exp_list = []

    for i in range(int(sim)):
        
        left_set = set(random.sample(random_bag, left))
        
        right_set = set(random.sample(random_bag, right))
        
        intersect = left_set.intersection(right_set)
        
        exp = len(intersect)
        
        if exp < obs:
            by_chance_list.append(1)
        else:
            by_chance_list.append(0)
            
        exp_list.append(exp)
    
    exp_median = np.median(exp_list)    
    
    if obs != 0:
        ratio = round(exp_median/(obs), 3) 
    else:
        ratio = exp_median
        
    outline = ('Median exp: {exp}, Num obs: {obs}, denrichment ratio: {ratio}, pval:{pval}\n').format(
        exp = exp_median,
        obs = obs,
        ratio = ratio,
        pval = round(sum(by_chance_list)/len(by_chance_list), sig_fig)
        )
    
    print(outline)


if False:
    infile = ('analyses/efficiency/unit_object_level2_ms_rpf_16_Unit_Sres_fdr.csv')
    df = pd.read_table(infile, index_col=0, sep = ',')
    
    copy_number_filename = ('metadata/chemostat_gene_relative_copy_number.tsv')
    cn_df = pd.read_table(copy_number_filename, index_col=0)
    cn_df = cn_df.add_suffix('_copy_number')
    
    df = pd.merge(left = df,
                  right = cn_df,
                  left_index=True,
                  right_index=True)
        
    rpf_col_name = ('DGY1657_ms_pt_median')
    rna_col_name = ('DGY1657_rpf_pt_median')
    df['anc_ratio'] = np.log2(df[rpf_col_name]/df[rna_col_name])
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        col_name = ('{}_fet_ratio_median').format(strain)
        rpf_col_name = ('{}_ms_pt_median').format(strain)
        rna_col_name = ('{}_rpf_pt_median').format(strain)
        evo_ratio = ('{}_evo_ratio').format(strain)
        df[evo_ratio] = np.log2(df[rpf_col_name]/df[rna_col_name])
        
    max_axis = -1*np.inf
    min_axis = np.inf
    
    for strain in set(['DGY1657', 'DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain == 'DGY1657':
            col_name = ('anc_ratio')
        else:
            col_name = ('{}_evo_ratio').format(strain)
            
        if max(df[col_name]) > max_axis:
            max_axis = max(df[col_name])
        if min(df[col_name]) < min_axis:
            min_axis = min(df[col_name])
            
    axis_wiggle = abs(max_axis-min_axis)*0.1
            
    max_axis+=axis_wiggle
    min_axis-=axis_wiggle
    
    each_cnv_dict = {}
    unit_object = df.to_dict('index')
    
    ''' '''
    complex_filename = ('analyses/complexes/table_of_complexes.csv')
    complex_file = open(complex_filename)

    complex_set = set()

    for line in complex_file:
        if line[0] != '#':
            #Complex	Subunit ID	Subunit Name	Subunit Stoichiometry
            gene = line.split(',')[1]
            
            if gene in unit_object:
                complex_set.add(gene)
            
            # if '_' in gene:
            #     gene_list = gene.split('_')
            #     for gene in gene_list:
                    
            #         if gene in unit_object:
            #             complex_set.add(gene)

    complex_list = list(complex_set)
    '''    '''

    sig_dict = {}    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            evo_ratio = ('{}_evo_ratio').format(strain)
                    
            pval_name = ("{}_ms_rpf_fdr").format(strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
            
            sres_col_name = ('sres_{}').format(strain)
            
            sn_col_name = ('{}_ms_rpf_signal_to_noise').format(strain)
            
            for gene in unit_object:
                # ratio = (unit_object[gene][evo_ratio] / 
                #          unit_object[gene]['anc_ratio'])
                
                #if unit_object[gene][evo_ratio] < unit_object[gene]['anc_ratio']:
                if True:
                    if ((unit_object[gene][pval_name] <= 0.05) or
                        (abs(unit_object[gene][sres_col_name]) >=2 )):
                        if unit_object[gene][sn_col_name]:
                            
                            if gene not in sig_dict:
                                sig_dict[gene] = set()
                                
                            sig_dict[gene].add(strain)


    many_sig_set = set()
    for gene in sig_dict:
        if len(sig_dict[gene]) > 1:
            many_sig_set.add(gene)
            
    multi_sig_list = list(many_sig_set)
    print(len(multi_sig_list))
    
    many_complexes = list(many_sig_set.intersection(complex_set))
    print(len(many_complexes))
    
    
    '''
    #neg select
    select_list = list(set(['YJL096W', 'YNR022C', 'YLR413W', 'YKR054C', 'YLR466W', 'YMR322C',
                            'YOR391C', 'YPL280W', 'YIL160C', 'YMR323W', 'YOR393W', 'YPL281C']))
    '''
    '''
    #neg mito select
    select_list = list(set(['YBR268W', 'YCR046C', 'YDR041W', 'YDR405W', 'YJL096W', 'YKR006C',
                            'YKR085C', 'YNL081C', 'YNR022C', 'YNR037C', ]))
    
    '''
    #YRF_SELECT
    #select_list = list(set(['YER190W', 'YGR296W', 'YLR466W', 'YNL339C', 'YPL283C']))    
    
    
    '''
    #pos_ YBR268W
    select_list = list(set(['YOR028C', 'YAL028W', 'YMR121C',
                            'YCR031C', 'YNL072W']))
    '''
    '''
    #ERR
    select_list = list(set(['YMR323W', 'YPL281C', 'YOR393W']))
    '''
    
    #select_list = multi_sig_list


    # many_sig_set = set()
    # for gene in sig_dict:
    #     if len( sig_dict[gene]) == 1:
    #         many_sig_set.add(gene)
            
    # one_sig_list = list(many_sig_set)
    # print(len(one_sig_list))
    
    many_sig_set = set()
    for gene in sig_dict:
        if len( sig_dict[gene]) == 2:
            many_sig_set.add(gene)
            
    two_sig_list = list(many_sig_set)
    print(len(two_sig_list))
    
    many_sig_set = set()
    for gene in sig_dict:
        if len(sig_dict[gene]) == 3:
            many_sig_set.add(gene)
            
    three_sig_list = list(many_sig_set)
    print(len(three_sig_list))
    
    many_sig_set = set()
    for gene in sig_dict:
        if len(sig_dict[gene]) == 4:
            many_sig_set.add(gene)
            
    four_sig_list = list(many_sig_set)
    print(len(four_sig_list))
    
    output_figure_name = ('figures/_background_peff_efficiency_pval.pdf')

    fig = go.Figure()
            
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            
            evo_ratio = ('{}_evo_ratio').format(strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
            cnv_df = df[df[cn_col_name] > 1]
                        
            sres_col_name = ('sres_{}').format(strain)
            pval_name = ("{}_ms_rpf_fdr").format(strain)
            
            sn_col_name = ('{}_ms_rpf_signal_to_noise').format(strain)
            
            #zero_sig_df = df.loc[complex_list]
            #zero_sig_df = df.loc[['YNL081C']]
            zero_sig_df = df.loc[many_complexes]
            zero_sig_df = zero_sig_df[(zero_sig_df[sn_col_name] == True)]
            zero_sig_df = zero_sig_df[(zero_sig_df[pval_name] <= 0.05) | (abs(zero_sig_df[sres_col_name]) >=2 )]
            zero_sig_index_vals = zero_sig_df[cn_col_name].astype('category').cat.codes
            
            ##Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=zero_sig_df['anc_ratio'],
                    y=zero_sig_df[evo_ratio],
                    text=zero_sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 100,
                        color = 'Green',
                        opacity = 0.5,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            '''
            sig_df = df.loc[two_sig_list]
            sig_df = sig_df[(sig_df[sn_col_name] == True)]
            sig_df = sig_df[(sig_df[pval_name] <= 0.05) | (abs(sig_df[sres_col_name]) >=2 )]

            sig_index_vals = sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=sig_df['anc_ratio'],
                    y=sig_df[evo_ratio],
                    text=sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 100,
                        color = '#6666FF',
                        opacity = 0.4,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            
            sig_df = df.loc[three_sig_list]
            sig_df = sig_df[(sig_df[sn_col_name] == True)]
            sig_df = sig_df[(sig_df[pval_name] <= 0.05) | (abs(sig_df[sres_col_name]) >=2 )]

            sig_index_vals = sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=sig_df['anc_ratio'],
                    y=sig_df[evo_ratio],
                    text=sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 100,
                        color = '#9966FF',
                        opacity = 0.5,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            
            sig_df = df.loc[three_sig_list]
            sig_df = sig_df[(sig_df[sn_col_name] == True)]
            sig_df = sig_df[(sig_df[pval_name] <= 0.05) | (abs(sig_df[sres_col_name]) >=2 )]

            sig_index_vals = sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=sig_df['anc_ratio'],
                    y=sig_df[evo_ratio],
                    text=sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 100,
                        color = '#FF66CC',
                        opacity = 0.6,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            '''
    plot_title = ('Protein efficiencies')
    xaxis_title = ("Ancestor strain log2(MS / RPF)")
    yaxis_title = ("Evolved strain log2(MS / RPF)")
    
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
    
    #fig.update_xaxes(range=[-4.5, 4.5])
    #fig.update_yaxes(range=[-4.5, 4.5])
    
    fig.update_xaxes(range=[min_axis, max_axis])
    fig.update_yaxes(range=[min_axis, max_axis])
    
    fig.show()
    #fig.write_image(output_figure_name)
    
    
sig_peff_set = set()
for gene in sig_dict:
    if len(sig_dict[gene]) > 0:
        print(sig_dict[gene])
        
        sig_peff_set.add(gene)
        
#ssd1_set

#high_qc_uorfs_set = set(list(high_qc_uorfs.keys()))  

''' are sig_peff enriched in complex_set? '''     
total_genes = 4289
sig_peff = len(sig_peff_set) 
complexes = len(complex_set)
overlap = len(sig_peff_set.intersection(complex_set))
two_way(total_genes, sig_peff, complexes, overlap, 10000) 
#Yes Median exp: 66.0, Num obs: 89, ratio: 1.348, pval:0.0015

''' are low sig_peff enriched in complex_set? '''     
total_genes = 4289
sig_peff = len(sig_peff_set) 
complexes = len(complex_set)
overlap = len(sig_peff_set.intersection(complex_set))
two_way(total_genes, sig_peff, complexes, overlap, 100000) 
#Yes! Median exp: 26.0, Num obs: 66, ratio: 2.538, pval:0.0



''' are sig_peff enriched in uorfs?'''
total_genes = 4289
sig_peff = len(sig_peff_set) 
uorf = len(high_qc_uorfs_set)
overlap = len(sig_peff_set.intersection(high_qc_uorfs_set))
two_way(total_genes, sig_peff, uorf, overlap, 10000) 
#No Median exp: 76.0, Num obs: 67, ratio: 0.882, pval:0.9062

two_way_deenrichment(total_genes, sig_peff, uorf, overlap, 10000)
#close to yes Median exp: 76.0, Num obs: 67, denrichment ratio: 1.134, pval:0.0932 

''' are sig_teff enriched in ssd1? '''     
total_genes = 4289
sig_peff = len(sig_peff_set)
ssd1 = len(ssd1_set)
overlap = len(sig_peff_set.intersection(ssd1_set))
two_way(total_genes, sig_peff, ssd1, overlap, 10000) 
#Almost Median exp: 22.0, Num obs: 29, ratio: 1.318, pval:0.0717

''' are complex_set enriched in uorfs? '''     
total_genes = 4289
uorfs = len(high_qc_uorfs_set) 
complexes = len(complex_set)
overlap = len(high_qc_uorfs_set.intersection(complex_set))
two_way(total_genes, uorfs, complexes, overlap, 10000) 
#no Median exp: 64.0, Num obs: 50, ratio: 0.781, pval:0.9842

''' are complex_set enriched in ssd1? '''     
total_genes = 4289
ssd1 = len(ssd1_set) 
complexes = len(complex_set)
overlap = len(ssd1_set.intersection(complex_set))
two_way(total_genes, ssd1, complexes, overlap, 10000) 
#no Median exp: 18.0, Num obs: 16, ratio: 0.889, pval:0.7727

select_list = list(sig_peff_set.intersection(complex_set))
print(select_list)

sig_dict = {}    
for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
    if strain != 'DGY1657':
        evo_ratio = ('{}_evo_ratio').format(strain)
                
        pval_name = ("{}_ms_rpf_fdr").format(strain)
        
        cn_col_name = ('{}_copy_number').format(strain)
        
        sres_col_name = ('sres_{}').format(strain)
        
        sn_col_name = ('{}_ms_rpf_signal_to_noise').format(strain)
        
        for gene in unit_object:
            # ratio = (unit_object[gene][evo_ratio] / 
            #          unit_object[gene]['anc_ratio'])
            
            if unit_object[gene][evo_ratio] < unit_object[gene]['anc_ratio']:
            #if True:
                if ((unit_object[gene][pval_name] <= 0.05) or
                    (abs(unit_object[gene][sres_col_name]) >=2 )):
                    if unit_object[gene][sn_col_name]:
                        
                        if gene not in sig_dict:
                            sig_dict[gene] = set()
                            
                        sig_dict[gene].add(strain)


many_sig_set = set()
for gene in sig_dict:
    if len(sig_dict[gene]) > 2:
        many_sig_set.add(gene)
        
multi_sig_list = list(many_sig_set)
print(len(multi_sig_list))

many_complexes = list(many_sig_set.intersection(complex_set))
print(len(many_complexes))

overlap = (sig_peff_set.intersection(complex_set))
protein_complex_degraded_set = set()
for gene in overlap:
    if len(sig_dict[gene]) > 2 or (gene[0:3] == 'YKR' and len(sig_dict[gene]) > 1):
        print(gene, sig_dict[gene])
        protein_complex_degraded_set.add(gene)
protein_complex_degraded_list = list(protein_complex_degraded_set)

if False:
    output_figure_name = ('figures/_degraded_complex_peff_efficiency_pval.pdf')

    fig = go.Figure()
            
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            
            evo_ratio = ('{}_evo_ratio').format(strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
            cnv_df = df[df[cn_col_name] > 1]
                        
            pval_name = ("{}_ms_rpf_fdr").format(strain)
            
            zero_sig_df = df.loc[protein_complex_degraded_list]
            #zero_sig_df = df.loc[many_complexes]
            zero_sig_df = zero_sig_df[(zero_sig_df[sn_col_name] == True)]
            zero_sig_df = zero_sig_df[(zero_sig_df[pval_name] <= 0.05) | (abs(zero_sig_df[sres_col_name]) >=2 )]
            zero_sig_index_vals = zero_sig_df[cn_col_name].astype('category').cat.codes
            
            ##Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=zero_sig_df['anc_ratio'],
                    y=zero_sig_df[evo_ratio],
                    text=zero_sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 0,
                        color = 'Red',
                        opacity = .8,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )

    plot_title = ('Protein efficiencies')
    xaxis_title = ("Ancestor strain log2(MS / RPF)")
    yaxis_title = ("Evolved strain log2(MS / RPF)")
    
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

    #fig.update_xaxes(range=[-4.5, 4.5])
    #fig.update_yaxes(range=[-4.5, 4.5])
    
    fig.update_xaxes(range=[min_axis, max_axis])
    fig.update_yaxes(range=[min_axis, max_axis])
    
    fig.show()
    fig.write_image(output_figure_name)