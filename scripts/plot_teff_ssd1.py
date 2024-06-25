# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:28:20 2024

Combine uORFs

@author: pspea
"""
import random
import plotly.io as pio
pio.renderers.default = "browser"

from scipy.stats import fisher_exact
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

unit_file_name = ('analyses/efficiency/unit_object_level_2v13v16_Unit_10pct.tab')
unit_object_df = pd.read_table(unit_file_name, index_col=0)
unit_object_dict = unit_object_df.to_dict('index')

genes_in_analysis = set(list(unit_object_dict.keys()))


global_uorf_dict = {}

for strain in set(['DGY1657', 'DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
    for rep in set(['rep1', 'rep2']):
        uorf_file_name = ('C:/Gresham/tiny_projects/quorf/'
                          '{strain}_{rep}/{strain}_{rep}_predictions_0.5.bed').format(
                              strain = strain, rep = rep)
                              
        sample_name = ('{strain}_{rep}').format(strain = strain, rep = rep)
        
        uorf_file = open(uorf_file_name)
                              
        for line in uorf_file:
            print(line)
            #II 504806   504852    YBR135W.504806.504852.xx  0.996429979801178  +
            uid = line.split('\t')[3]
            pred_score = float(line.split('\t')[4])
            
            if pred_score >= 0.9:
                gene = uid.split('.')[0]
                
                if gene in genes_in_analysis:
                    if gene not in global_uorf_dict:
                        global_uorf_dict[gene] = set()
                        
                    global_uorf_dict[gene].add(sample_name)
            
        uorf_file.close()
        
        
present_in_control = set()
evo_strains = {}
evo_genes = {}

for gene in global_uorf_dict:
    samples = global_uorf_dict[gene]
    
    if ('DGY1657_rep1' in samples) and ('DGY1657_rep2' in samples):
        present_in_control.add(gene)
        
    else:
        for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
            sample_name_r1 = ('{strain}_rep1').format(strain = strain)
            sample_name_r2 = ('{strain}_rep2').format(strain = strain)
            
            if (sample_name_r1 in samples) and (sample_name_r2 in samples):
            
                if strain not in evo_strains:
                    evo_strains[strain] = set()
                    
                evo_strains[strain].add(gene)
                
                if gene not in evo_genes:
                    evo_genes[gene] = set()
                    
                evo_genes[gene].add(strain)
                
print(list(evo_genes.keys()))

print(list(present_in_control))

len(present_in_control)
    
                    
high_qc_uorfs = {}

for gene in global_uorf_dict:
    samples = global_uorf_dict[gene]
    
    for strain in set(['DGY1657', 'DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        sample_name_r1 = ('{strain}_rep1').format(strain = strain)
        sample_name_r2 = ('{strain}_rep2').format(strain = strain)
        
        if (sample_name_r1 in samples) and (sample_name_r2 in samples):
                    
            if gene not in high_qc_uorfs:
                high_qc_uorfs[gene] = set()
                
            high_qc_uorfs[gene].add(strain)
          
print(list(high_qc_uorfs))

len(high_qc_uorfs)


'''
'''
if True:    
    infile = ('analyses/efficiency/Unit_v13_teff_results.tab')
    df = pd.read_table(infile, index_col=0)
    
    copy_number_filename = ('metadata/chemostat_gene_relative_copy_number.tsv')
    cn_df = pd.read_table(copy_number_filename, index_col=0)
    cn_df = cn_df.add_suffix('_copy_number')
    
    df = pd.merge(left = df,
                  right = cn_df,
                  left_index=True,
                  right_index=True)
    
    rpf_col_name = ('CDS_RPF_DGY1657_median_tpm')
    rna_col_name = ('CDS_RNA_DGY1657_median_tpm')
    df['anc_ratio'] = np.log2(df[rpf_col_name]/df[rna_col_name])
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        col_name = ('{}_fet_ratio_median').format(strain)
        rpf_col_name = ('CDS_RPF_{}_median_tpm').format(strain)
        rna_col_name = ('CDS_RNA_{}_median_tpm').format(strain)
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
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            pval_name = ("{strain}_fet_pval_median").format(
                strain = strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
                        
            for gene in unit_object:
                if unit_object[gene][cn_col_name] != 1:
                    if gene not in each_cnv_dict:
                        each_cnv_dict[gene] = {}
                    
                    pval = unit_object[gene][pval_name]
                    each_cnv_dict[gene][strain] = pval

    for gene in each_cnv_dict:
        if len(each_cnv_dict[gene]) > 1:
            show = True
            for strain in each_cnv_dict[gene]:
                if each_cnv_dict[gene][strain] > 0.05:
                    show = False
            if show:
                print(gene, each_cnv_dict[gene])
    
    many_sig_dict = {}
    unit_object = df.to_dict('index')
    
    infile = ('analyses/ssd1/SSD1_hits_in_TL_from_SGD_SN.txt')
    ssd1_df = pd.read_table(infile, index_col=0)
    ssd1_dict = ssd1_df.to_dict('index')
    
    ct = 0
    ssd1_set = set()
    for gene in ssd1_dict:
        if (ssd1_dict[gene]['pval'] <= 0.05) and (ssd1_dict[gene]['hits'] > 0):
            ct+=1
            
            if gene in unit_object:
                ssd1_set.add(gene)
    ssd1_list = list(ssd1_set)
    
    ssd1_sig_teff = []
    
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            pval_name = ("{strain}_fet_pval_median").format(
                strain = strain)
                        
            for gene in unit_object:
                if gene not in many_sig_dict:
                    many_sig_dict[gene] = set()
                
                if unit_object[gene][pval_name] <= 0.05:
                    many_sig_dict[gene].add(strain)
                  
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            
            # output_figure_name = ('C:/Gresham/Project_Carolino/figures/efficiency'
            #                       '/_{strain}_teff_uorf_efficiency_pval.pdf').format(strain = strain)
            
            output_figure_name = ('figures/_{strain}_teff_ssd1_efficiency_pval.pdf').format(strain = strain)

            fig = go.Figure()
            
            # select_list = []
            # for gene in high_qc_uorfs:
            #     if (strain in high_qc_uorfs[gene]) or ('DGY1657' in high_qc_uorfs[gene]):
            #         select_list.append(gene)
                    
            select_list = ssd1_list
            
            evo_ratio = ('{}_evo_ratio').format(strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
            cnv_df = df[df[cn_col_name] > 1]
                        
            pval_name = ("{strain}_fet_pval_median").format(
                strain = strain)
                        
            cn_col_name = ('{}_copy_number').format(strain)
            
            zero_sig_df = df.loc[select_list]
            zero_sig_df = zero_sig_df[zero_sig_df[pval_name] <= 0.05]
            zero_sig_index_vals = zero_sig_df[cn_col_name].astype('category').cat.codes
            
            # many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=zero_sig_df['anc_ratio'],
                    y=zero_sig_df[evo_ratio],
                    text=zero_sig_index_vals.index + ' ' + strain,
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

            plot_title = ('Translation efficiencies_ssd1_{strain}').format(strain = strain)
            xaxis_title = ("Ancestor strain log2(RPF / RNA)")
            yaxis_title = ("Evolved strain log2(RPF / RNA)")
            
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
            
sig_teff_set = set()
for gene in many_sig_dict:
    if len(many_sig_dict[gene]) > 0:
        print(many_sig_dict[gene])
        
        sig_teff_set.add(gene)
        
#ssd1_set


high_qc_uorfs_set = set(list(high_qc_uorfs.keys()))  
'''

'''
uorf_cnv_set = set()
cnv_set = set() 
uorf_cnn_set = set() 
cnn_set = set()
for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
    select_list = list(high_qc_uorfs_set)
    
    evo_ratio = ('{}_evo_ratio').format(strain)
    
    cn_col_name = ('{}_copy_number').format(strain)
                
    pval_name = ("{strain}_fet_pval_median").format(
        strain = strain)
                
    cn_col_name = ('{}_copy_number').format(strain)
    
    uorf_df = df.loc[select_list]
    temp_df = uorf_df[(uorf_df[pval_name] <= 0.05) & (uorf_df[cn_col_name] != 1)]
    uorf_cnv_set = uorf_cnv_set.union(set(temp_df.index))
    
    temp_df = df[(df[pval_name] <= 0.05) & (df[cn_col_name] != 1)]
    cnv_set = cnv_set.union(set(temp_df.index))
    
    temp_df = uorf_df[(uorf_df[pval_name] <= 0.05) & (uorf_df[cn_col_name] == 1)]
    uorf_cnn_set = uorf_cnn_set.union(set(temp_df.index))
    
    temp_df = df[(df[pval_name] <= 0.05) & (df[cn_col_name] == 1)]
    cnn_set = cnn_set.union(set(temp_df.index))
    
odds, pval = fisher_exact([[len(uorf_cnv_set), 
                            len(cnv_set)],
                           [len(uorf_cnn_set), 
                            len(cnn_set)]], alternative='two-sided')
print('cnv_uorfs', odds, pval)

cnv_uorf_select_list = list(uorf_cnv_set)

ssd1_cnv_set = set()
for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
    select_list = list(ssd1_set)
    
    evo_ratio = ('{}_evo_ratio').format(strain)
    
    cn_col_name = ('{}_copy_number').format(strain)
                
    pval_name = ("{strain}_fet_pval_median").format(
        strain = strain)
                
    cn_col_name = ('{}_copy_number').format(strain)
    
    ssd1_df = df.loc[select_list]
    temp_df = ssd1_df[(ssd1_df[pval_name] <= 0.05) & (ssd1_df[cn_col_name] != 1)]
    ssd1_cnv_set = uorf_cnv_set.union(set(temp_df.index))
    
cnv_ssd1_select_list = list(ssd1_cnv_set)

###
name_list = ['YDR293C','YHR052W', 'YHR143W', 'YGR279C', 'YLR110C', 'YLR286C', 'YMR305C', 'YNL066W', 'YOR247W', 'YPL256C', 'YBR162C']
for name in name_list:
    print(name, name in high_qc_uorfs_set, name in sig_teff_set)


                             
''' are sig_teff enriched in uorfs? '''     
total_genes = 4289
sig_teff = len(sig_teff_set) 
uorf = len(high_qc_uorfs_set)
overlap = len(sig_teff_set.intersection(high_qc_uorfs_set))
two_way(total_genes, sig_teff, uorf, overlap, 1000) 
#Yes Median exp: 135.0, Num obs: 216, ratio: 1.6, pval:0.0

''' are sig_teff enriched in ssd1? '''     
total_genes = 4289
sig_teff = len(sig_teff_set) 
ssd1 = len(ssd1_set)
overlap = len(sig_teff_set.intersection(ssd1_set))
two_way(total_genes, sig_teff, ssd1, overlap, 1000) 
#Yes Median exp: 39.0, Num obs: 56, ratio: 1.436, pval:0.0021

''' are uorfs enriched in ssd1? '''     
total_genes = 4289
uorfs = len(high_qc_uorfs_set) 
ssd1 = len(ssd1_set)
overlap = len(high_qc_uorfs_set.intersection(ssd1_set))
two_way(total_genes, uorfs, ssd1, overlap, 10000) 
#Yes Median exp: 21.0, Num obs: 50, ratio: 2.381, pval:0.0

''' are sig teff associate uorfs enriched in ssd1? '''     
sig_teff = len(sig_teff_set)
sig_teff_uorfs = (sig_teff_set.intersection(high_qc_uorfs_set))
sig_teff_ssd1 = (sig_teff_set.intersection(ssd1_set))
overlap = len(sig_teff_uorfs.intersection(sig_teff_ssd1))
two_way(sig_teff, len(sig_teff_uorfs), len(sig_teff_ssd1), overlap, 100000) 
#Yes Median exp: 12.0, Num obs: 28, ratio: 2.333, pval:0.0

ssd1_select_list = list(sig_teff_uorfs.intersection(sig_teff_ssd1))


if False:    
    infile = ('analyses/efficiency/Unit_v13_teff_results.tab')
    df = pd.read_table(infile, index_col=0)
    
    copy_number_filename = ('metadata/chemostat_gene_relative_copy_number.tsv')
    cn_df = pd.read_table(copy_number_filename, index_col=0)
    cn_df = cn_df.add_suffix('_copy_number')
    
    df = pd.merge(left = df,
                  right = cn_df,
                  left_index=True,
                  right_index=True)
    
    rpf_col_name = ('CDS_RPF_DGY1657_median_tpm')
    rna_col_name = ('CDS_RNA_DGY1657_median_tpm')
    df['anc_ratio'] = np.log2(df[rpf_col_name]/df[rna_col_name])
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        col_name = ('{}_fet_ratio_median').format(strain)
        rpf_col_name = ('CDS_RPF_{}_median_tpm').format(strain)
        rna_col_name = ('CDS_RNA_{}_median_tpm').format(strain)
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
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            pval_name = ("{strain}_fet_pval_median").format(
                strain = strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
                        
            for gene in unit_object:
                if unit_object[gene][cn_col_name] != 1:
                    if gene not in each_cnv_dict:
                        each_cnv_dict[gene] = {}
                    
                    pval = unit_object[gene][pval_name]
                    each_cnv_dict[gene][strain] = pval

    for gene in each_cnv_dict:
        if len(each_cnv_dict[gene]) > 1:
            show = True
            for strain in each_cnv_dict[gene]:
                if each_cnv_dict[gene][strain] > 0.05:
                    show = False
            if show:
                print(gene, each_cnv_dict[gene])
    
    many_sig_dict = {}
    unit_object = df.to_dict('index')
        
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            pval_name = ("{strain}_fet_pval_median").format(
                strain = strain)
                        
            for gene in unit_object:
                if gene not in many_sig_dict:
                    many_sig_dict[gene] = set()
                
                if unit_object[gene][pval_name] <= 0.05:
                    many_sig_dict[gene].add(strain)
                    
    output_figure_name = ('figures/_all_strain_teff_cnv_ssd1_efficiency_pval.pdf').format()

    fig = go.Figure()
                  
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            
            # output_figure_name = ('C:/Gresham/Project_Carolino/figures/efficiency'
            #                       '/_{strain}_teff_uorf_efficiency_pval.pdf').format(strain = strain)
            

                        
            evo_ratio = ('{}_evo_ratio').format(strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
            cnv_df = df[df[cn_col_name] > 1]
                        
            pval_name = ("{strain}_fet_pval_median").format(
                strain = strain)
                        
            cn_col_name = ('{}_copy_number').format(strain)
            
            # zero_sig_df = df.loc[ssd1_select_list]
            # zero_sig_df = df.loc[cnv_uorf_select_list]
            zero_sig_df = df.loc[cnv_ssd1_select_list]
            zero_sig_df = zero_sig_df[zero_sig_df[pval_name] <= 0.05]
            zero_sig_index_vals = zero_sig_df[cn_col_name].astype('category').cat.codes
            
            # many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=zero_sig_df['anc_ratio'],
                    y=zero_sig_df[evo_ratio],
                    text=zero_sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 1,
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

    plot_title = ('Translation efficiencies_ssd1_uorfs_{strain}').format(strain = strain)
    xaxis_title = ("Ancestor strain log2(RPF / RNA)")
    yaxis_title = ("Evolved strain log2(RPF / RNA)")
    
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


        
