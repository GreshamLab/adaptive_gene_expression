# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 08:48:17 2023

@author: pspea
"""
import plotly.io as pio
pio.renderers.default = "browser"

from scipy.stats import fisher_exact
import pandas as pd
import numpy as np

#
#from scipy.stats import mannwhitneyu
import plotly.graph_objects as go
#import plotly.express as px

#from scipy import stats


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
        #col_name = ('{}_fet_ratio_median').format(strain)
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

    # many_sig_set = set()
    # for gene in many_sig_dict:
    #     if len( many_sig_dict[gene]) == 0:
    #         many_sig_set.add(gene)
            
    # zero_sig_list = list(many_sig_set)
    # print(len(zero_sig_list))
    
    select_list = list(set(['YKR039W', 'YAL005C', 'YBL075C', 'YEL009C', 'YER103W', 'YKR042W', 'YPR016C', 'YOR375C', 'YFL014W', 
                   'YDR171W', 'YJL159W', 'YJR070C', 'YMR037C', 'YKL062W', 'YLL024C', 'YDL229W', 'YNL209W', 'YJR045C',
                   'YBR082C', 'YEL012W', 'YMR022W', 'YBR173C', 'YKR046C', 'YKR057W', 'YKR065C', 'YKR080W', 'YKR094C']))
    
    select_list = list(set(['YAL003W', 'YOR375C', 'YJR109C', ]))    

    many_sig_set = set()
    for gene in many_sig_dict:
        if len( many_sig_dict[gene]) == 1:
            many_sig_set.add(gene)
            
    one_sig_list = list(many_sig_set)
    print(len(one_sig_list))
    
    many_sig_set = set()
    for gene in many_sig_dict:
        if len( many_sig_dict[gene]) == 2:
            many_sig_set.add(gene)
            
    two_sig_list = list(many_sig_set)
    print(len(two_sig_list))
    
    many_sig_set = set()
    for gene in many_sig_dict:
        if len(many_sig_dict[gene]) == 3:
            many_sig_set.add(gene)
            
    three_sig_list = list(many_sig_set)
    print(len(three_sig_list))
    
    many_sig_set = set()
    for gene in many_sig_dict:
        if len(many_sig_dict[gene]) == 4:
            many_sig_set.add(gene)
            
    four_sig_list = list(many_sig_set)
    print(len(four_sig_list))
    
    output_figure_name = ('figures/_count_more_sig_teff_efficiency_pval.pdf')

    fig = go.Figure()
            
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            
            evo_ratio = ('{}_evo_ratio').format(strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
            cnv_df = df[df[cn_col_name] > 1]
                        
            pval_name = ("{strain}_fet_pval_median").format(
                strain = strain)
                        
            cn_col_name = ('{}_copy_number').format(strain)
            
            zero_sig_df = df.loc[select_list]
            #zero_sig_df = zero_sig_df[zero_sig_df[pval_name] <= 0.05]
            zero_sig_index_vals = zero_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            # fig.add_trace(
            #     go.Scatter(
            #         mode='markers',
            #         x=zero_sig_df['anc_ratio'],
            #         y=zero_sig_df[evo_ratio],
            #         text=zero_sig_index_vals.index + ' ' + strain,
            #         marker=dict(
            #             symbol = 100,
            #             color = 'Green',
            #             opacity = 0.5,
            #             line=dict(
            #                 color='Black',
            #                 width=1
            #             )
            #         ),
            #         showlegend=False
            #     )
            # )
            '''
            one_sig_df = df.loc[one_sig_list]
            one_sig_df = one_sig_df[one_sig_df[pval_name] <= 0.05]
            one_sig_index_vals = one_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=one_sig_df['anc_ratio'],
                    y=one_sig_df[evo_ratio],
                    text=one_sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 0,
                        color = '#66CCFF',
                        opacity = 0.2,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            '''
            two_sig_df = df.loc[two_sig_list]
            two_sig_df = two_sig_df[two_sig_df[pval_name] <= 0.05]
            two_sig_index_vals = two_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=two_sig_df['anc_ratio'],
                    y=two_sig_df[evo_ratio],
                    text=two_sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 100,
                        color = '#6666FF',
                        opacity = 0.3,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            
            three_sig_df = df.loc[three_sig_list]
            three_sig_df = three_sig_df[three_sig_df[pval_name] <= 0.05]
            three_sig_index_vals = three_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=three_sig_df['anc_ratio'],
                    y=three_sig_df[evo_ratio],
                    text=three_sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 100,
                        color = '#9966FF',
                        opacity = 0.4,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            
            four_sig_df = df.loc[four_sig_list]
            four_sig_df = four_sig_df[four_sig_df[pval_name] <= 0.05]
            four_sig_index_vals = four_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=four_sig_df['anc_ratio'],
                    y=four_sig_df[evo_ratio],
                    text=four_sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 100,
                        color = '#FF66CC',
                        opacity = 0.5,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            
    plot_title = ('Translation efficiencies_{strain}').format(strain = strain)
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
    
    #fig.show()
    #fig.write_image(output_figure_name)
    ''' ''' ''' '''    
    gene_counter = {}
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        cnn_set = set()             
        cnv_set = set()
        is_high_cnn_set = set()
        is_low_cnn_set = set()
        is_high_cnv_set = set()
        is_low_cnv_set = set()
        
        pval_name = ("{strain}_fet_pval_median").format(
            strain = strain)
        
        cn_col_name = ('{}_copy_number').format(strain)
        difference_col_name = ('{}_difference').format(strain)
        
        evo_ratio = ('{}_evo_ratio').format(strain)
        
        df[difference_col_name] = df[evo_ratio] > df['anc_ratio']
        '''
        '''
        subset_df = df[(df[pval_name] > 0.05) & 
                       (df[cn_col_name] == 1)]
        cnn_set = cnn_set.union(set(subset_df.index))
        
        subset_df = df[(df[pval_name] > 0.05) & 
                       (df[cn_col_name] != 1)]
        cnv_set = cnv_set.union(set(subset_df.index))
        '''
        '''
        subset_df = df[(df[pval_name] <= 0.05) & 
                       (df[cn_col_name] == 1) &
                       (df[difference_col_name] == True) ]
        is_high_cnn_set = is_high_cnn_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['high'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                       (df[cn_col_name] == 1) &
                       (df[difference_col_name] == False) ]
        is_low_cnn_set = is_low_cnn_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['low'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                        (df[cn_col_name] != 1) &
                        (df[difference_col_name] == True) ]
        is_high_cnv_set = is_high_cnv_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['high'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                        (df[cn_col_name] != 1) &
                        (df[difference_col_name] == False) ]
        is_low_cnv_set = is_low_cnv_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['low'].add(strain)
        
        '''
        '''
        print(strain)
        # print(len(cnn_set) + len(is_high_cnn_set) + len(is_low_cnn_set), len(is_high_cnn_set), len(is_low_cnn_set))
        # print(len(cnv_set) + len(is_high_cnv_set) + len(is_low_cnv_set), len(is_high_cnv_set), len(is_low_cnv_set))

        print(len(cnn_set) + len(is_high_cnn_set) + len(is_low_cnn_set), len(is_high_cnn_set) + len(is_low_cnn_set))
        print(len(cnv_set) + len(is_high_cnv_set) + len(is_low_cnv_set), len(is_high_cnv_set) + len(is_low_cnv_set))
    
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

    # median pct teff for both CNNs
    100*np.median([610/4274, 332/4222, 656/4210, 336/4207])
    
    # median pct teff for both CNVs
    100*np.median([5/17, 11/69, 14/81, 9/84])

    # median pct teff for down CNNs
    100*np.median([131/4274, 92/4222, 218/4210, 112/4207])
    
    # median pct teff for down CNVs
    100*np.median([2/17, 4/69, 7/81, 5/84])
    
    cnn_set = set()             
    cnv_set = set()
    is_high_cnn_set = set()
    is_low_cnn_set = set()
    is_high_cnv_set = set()
    is_low_cnv_set = set()
    
    gene_counter = {}
    
    #for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
    #for strain in set(['DGY1726', 'DGY1735']):
    for strain in set(['DGY1741', 'DGY1743']):
        pval_name = ("{strain}_fet_pval_median").format(
            strain = strain)
        
        cn_col_name = ('{}_copy_number').format(strain)
        difference_col_name = ('{}_difference').format(strain)
        
        evo_ratio = ('{}_evo_ratio').format(strain)
        
        df[difference_col_name] = df[evo_ratio] > df['anc_ratio']
        '''
        '''
        subset_df = df[(df[cn_col_name] == 1)]
        cnn_set = cnn_set.union(set(subset_df.index))
        
        subset_df = df[(df[cn_col_name] != 1)]
        cnv_set = cnv_set.union(set(subset_df.index))
        '''
        '''
        subset_df = df[(df[pval_name] <= 0.05) & 
                       (df[cn_col_name] == 1) &
                       (df[difference_col_name] == True) ]
        is_high_cnn_set = is_high_cnn_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['high'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                       (df[cn_col_name] == 1) &
                       (df[difference_col_name] == False) ]
        is_low_cnn_set = is_low_cnn_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['low'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                        (df[cn_col_name] != 1) &
                        (df[difference_col_name] == True) ]
        is_high_cnv_set = is_high_cnv_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['high'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                        (df[cn_col_name] != 1) &
                        (df[difference_col_name] == False) ]
        is_low_cnv_set = is_low_cnv_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['low'].add(strain)
        
        '''
        '''
    print(len(cnn_set), len(is_high_cnn_set), len(is_low_cnn_set))
    print(len(cnv_set), len(is_high_cnv_set), len(is_low_cnv_set))

    odds, pval = fisher_exact([[len(is_high_cnv_set), 
                                len(cnv_set)],
                               [len(is_high_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('high', odds, pval)
    #high 0.7615823235923022 0.48165203243264676

    odds, pval = fisher_exact([[len(is_low_cnv_set), 
                                len(cnv_set)],
                               [len(is_low_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('low', odds, pval)
    #low 1.6964460185637518 0.09376665950622794
    
    odds, pval = fisher_exact([[len(is_high_cnv_set), 
                                len(cnv_set)],
                               [len(is_high_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('g150 high', odds, pval)
    #g150 high 1.2359484777517564 0.5598926217439615

    odds, pval = fisher_exact([[len(is_low_cnv_set), 
                                len(cnv_set)],
                               [len(is_low_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('g150 low', odds, pval)
    #g150 low 2.226793248945148 0.04965051461096953
    
    odds, pval = fisher_exact([[len(is_high_cnv_set), 
                                len(cnv_set)],
                               [len(is_high_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('g250 high', odds, pval)
    #g250 high 0.6544549717676332 0.33885006456284494

    odds, pval = fisher_exact([[len(is_low_cnv_set), 
                                len(cnv_set)],
                               [len(is_low_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('g250 low', odds, pval)
    #g250 low 1.445628276678505 0.2797721237938388
    
    ct = 0   
    to_print_set = set()
    for gene in gene_counter:
        
        
        hsize = len(gene_counter[gene]['high'])
        
        lsize = len(gene_counter[gene]['low'])

        # if (hsize > 1):
        #     print(gene, gene_counter[gene])
            
        # if (lsize > 1):
        #     ct+=1
        #     print(gene, gene_counter[gene])

        # if (hsize > 0) and (lsize == 0):
        #     ct+=1
        #     print(gene, gene_counter[gene])
        
        # if (hsize > 1) and (lsize == 0):
        #     ct+=1
        #     print(gene, gene_counter[gene])
        #     to_print_set.add(gene)
                
        # if (hsize == 0) and (lsize > 0):
        #     ct+=1
        #     print(gene, gene_counter[gene])
        
        # if (hsize == 0) and (lsize > 1):
        #     ct+=1
        #     print(gene, gene_counter[gene])
        #     to_print_set.add(gene)
            
        if (hsize > 0) and (lsize > 0):
            ct+=1
            print(gene, gene_counter[gene])
            to_print_set.add(gene)
        
    print(ct)
   
if False:
        
        
    infile = ('analyses/efficiency/unit_object_level2_ms_rpf_14_Unit_Sres_fdr.csv')
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
        if len(sig_dict[gene]) > 2:
            many_sig_set.add(gene)
            
    multi_sig_list = list(many_sig_set)
    print(len(multi_sig_list))
    
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
    select_list = list(set(['YER190W', 'YGR296W', 'YLR466W', 'YNL339C', 'YPL283C']))    
    
    
    '''
    #pos_ YBR268W
    select_list = list(set(['YOR028C', 'YAL028W', 'YMR121C',
                            'YCR031C', 'YNL072W']))
    '''
    '''
    #ERR
    select_list = list(set(['YMR323W', 'YPL281C', 'YOR393W']))
    '''



    many_sig_set = set()
    for gene in sig_dict:
        if len( sig_dict[gene]) == 1:
            many_sig_set.add(gene)
            
    one_sig_list = list(many_sig_set)
    print(len(one_sig_list))
    
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
    
    output_figure_name = ('figures/_select_cn_peff_efficiency_pval.pdf')

    fig = go.Figure()
            
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        if strain != 'DGY1657':
            
            evo_ratio = ('{}_evo_ratio').format(strain)
            pval_name = ("{}_ms_rpf_fdr").format(strain)
            
            cn_col_name = ('{}_copy_number').format(strain)
            cnv_df = df[df[cn_col_name] > 1]
                        
            cnv_sig_df = cnv_df[cnv_df[pval_name] <= 0.05]
            cnv_sig_index_vals = cnv_sig_df[cn_col_name].astype('category').cat.codes
            
            ##Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=cnv_sig_df['anc_ratio'],
                    y=cnv_sig_df[evo_ratio],
                    text=cnv_sig_index_vals.index + ' ' + strain,
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
            
            
            
            # zero_sig_df = df.loc[select_list]
            # zero_sig_df = zero_sig_df[zero_sig_df[pval_name] <= 0.05]
            # zero_sig_index_vals = zero_sig_df[cn_col_name].astype('category').cat.codes
            
            # ##Add many teff
            # fig.add_trace(
            #     go.Scatter(
            #         mode='markers',
            #         x=zero_sig_df['anc_ratio'],
            #         y=zero_sig_df[evo_ratio],
            #         text=zero_sig_index_vals.index + ' ' + strain,
            #         marker=dict(
            #             symbol = 100,
            #             color = 'Green',
            #             opacity = 0.5,
            #             line=dict(
            #                 color='Black',
            #                 width=1
            #             )
            #         ),
            #         showlegend=False
            #     )
            # )
            '''
            one_sig_df = df.loc[one_sig_list]
            one_sig_df = one_sig_df[one_sig_df[pval_name] <= 0.05]
            one_sig_index_vals = one_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=one_sig_df['anc_ratio'],
                    y=one_sig_df[evo_ratio],
                    text=one_sig_index_vals.index + ' ' + strain,
                    marker=dict(
                        symbol = 100,
                        color = '#66CCFF',
                        opacity = 0.3,
                        line=dict(
                            color='Black',
                            width=1
                        )
                    ),
                    showlegend=False
                )
            )
            
            two_sig_df = df.loc[two_sig_list]
            two_sig_df = two_sig_df[two_sig_df[pval_name] <= 0.05]
            two_sig_index_vals = two_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=two_sig_df['anc_ratio'],
                    y=two_sig_df[evo_ratio],
                    text=two_sig_index_vals.index + ' ' + strain,
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
            
            three_sig_df = df.loc[three_sig_list]
            three_sig_df = three_sig_df[three_sig_df[pval_name] <= 0.05]
            three_sig_index_vals = three_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=three_sig_df['anc_ratio'],
                    y=three_sig_df[evo_ratio],
                    text=three_sig_index_vals.index + ' ' + strain,
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
            
            four_sig_df = df.loc[four_sig_list]
            four_sig_df = four_sig_df[four_sig_df[pval_name] <= 0.05]
            four_sig_index_vals = four_sig_df[cn_col_name].astype('category').cat.codes
            
            #Add many teff
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=four_sig_df['anc_ratio'],
                    y=four_sig_df[evo_ratio],
                    text=four_sig_index_vals.index + ' ' + strain,
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
    
    gene_counter = {}
    
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        cnn_set = set()             
        cnv_set = set()
        is_high_cnn_set = set()
        is_low_cnn_set = set()
        is_high_cnv_set = set()
        is_low_cnv_set = set()
        
        pval_name = ("{strain}_ms_rpf_fdr").format(
            strain = strain)
        
        cn_col_name = ('{}_copy_number').format(strain)
        difference_col_name = ('{}_difference').format(strain)
        
        evo_ratio = ('{}_evo_ratio').format(strain)
        
        df[difference_col_name] = df[evo_ratio] > df['anc_ratio']
        '''
        '''
        subset_df = df[(df[pval_name] > 0.05) & 
                       (df[cn_col_name] == 1)]
        cnn_set = cnn_set.union(set(subset_df.index))
        
        subset_df = df[(df[pval_name] > 0.05) & 
                       (df[cn_col_name] != 1)]
        cnv_set = cnv_set.union(set(subset_df.index))
        '''
        '''
        subset_df = df[(df[pval_name] <= 0.05) & 
                       (df[cn_col_name] == 1) &
                       (df[difference_col_name] == True) ]
        is_high_cnn_set = is_high_cnn_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['high'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                       (df[cn_col_name] == 1) &
                       (df[difference_col_name] == False) ]
        is_low_cnn_set = is_low_cnn_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['low'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                        (df[cn_col_name] != 1) &
                        (df[difference_col_name] == True) ]
        is_high_cnv_set = is_high_cnv_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['high'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                        (df[cn_col_name] != 1) &
                        (df[difference_col_name] == False) ]
        is_low_cnv_set = is_low_cnv_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['low'].add(strain)
        
        '''
        '''
        print(len(cnn_set) + len(is_high_cnn_set) + len(is_low_cnn_set), len(is_high_cnn_set), len(is_low_cnn_set))
        print(len(cnv_set) + len(is_high_cnv_set) + len(is_low_cnv_set), len(is_high_cnv_set), len(is_low_cnv_set))

        #print(len(cnn_set) + len(is_high_cnn_set) + len(is_low_cnn_set), len(is_high_cnn_set) + len(is_low_cnn_set))
        #print(len(cnv_set) + len(is_high_cnv_set) + len(is_low_cnv_set), len(is_high_cnv_set) + len(is_low_cnv_set))
    
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

    # median pct peff for both CNNs
    100*np.median([174/4272, 147/4220, 184/4208, 161/4205])
    
    # median pct peff for both CNVs
    100*np.median([3/17, 10/69, 6/81, 7/84])

    # median pct peff for down CNNs
    100*np.median([32/4272, 37/4220, 19/4208, 38/4205])
    
    # median pct peff for down CNVs
    100*np.median([3/17, 8/69, 4/81, 6/84])
    
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
        pval_name = ("{strain}_ms_rpf_fdr").format(
            strain = strain)
        
        cn_col_name = ('{}_copy_number').format(strain)
        difference_col_name = ('{}_difference').format(strain)
        
        evo_ratio = ('{}_evo_ratio').format(strain)
        
        df[difference_col_name] = df[evo_ratio] > df['anc_ratio']
        '''
        '''
        subset_df = df[(df[cn_col_name] == 1)]
        cnn_set = cnn_set.union(set(subset_df.index))
        
        subset_df = df[(df[cn_col_name] != 1)]
        cnv_set = cnv_set.union(set(subset_df.index))
        '''
        '''
        # subset_df = df[(df[pval_name] <= 0.05) & 
        #                (df[cn_col_name] == 1) &
        #                (df[difference_col_name] == True) ]
        # is_high_cnn_set = is_high_cnn_set.union(set(subset_df.index))
        
        # for gene in set(subset_df.index):
        #     if gene not in gene_counter:
        #         gene_counter[gene] = {'high':set(), 'low':set()}
                
        #     gene_counter[gene]['high'].add(strain)
        
        # subset_df = df[(df[pval_name] <= 0.05) & 
        #                 (df[cn_col_name] == 1) &
        #                 (df[difference_col_name] == False) ]
        # is_low_cnn_set = is_low_cnn_set.union(set(subset_df.index))
        
        # for gene in set(subset_df.index):
        #     if gene not in gene_counter:
        #         gene_counter[gene] = {'high':set(), 'low':set()}
                
        #     gene_counter[gene]['low'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                        (df[cn_col_name] != 1) &
                        (df[difference_col_name] == True) ]
        is_high_cnv_set = is_high_cnv_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['high'].add(strain)
        
        subset_df = df[(df[pval_name] <= 0.05) & 
                        (df[cn_col_name] != 1) &
                        (df[difference_col_name] == False) ]
        is_low_cnv_set = is_low_cnv_set.union(set(subset_df.index))
        
        for gene in set(subset_df.index):
            if gene not in gene_counter:
                gene_counter[gene] = {'high':set(), 'low':set()}
                
            gene_counter[gene]['low'].add(strain)
        
        '''
        '''
    print(len(cnn_set), len(is_high_cnn_set), len(is_low_cnn_set))
    print(len(cnv_set), len(is_high_cnv_set), len(is_low_cnv_set))

    odds, pval = fisher_exact([[len(is_high_cnv_set), 
                                len(cnv_set)],
                               [len(is_high_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('high', odds, pval)
    #high 0.4820219647961486 0.28308356071766677

    odds, pval = fisher_exact([[len(is_low_cnv_set), 
                                len(cnv_set)],
                               [len(is_low_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('low', odds, pval)
    # low 7.429565217391304 5.763689233669496e-07
    
    
    odds, pval = fisher_exact([[len(is_high_cnv_set), 
                                len(cnv_set)],
                               [len(is_high_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('g150 high', odds, pval)
    #g150 high 0.6089466089466089 0.7716125394542112

    odds, pval = fisher_exact([[len(is_low_cnv_set), 
                                len(cnv_set)],
                               [len(is_low_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('g150 low', odds, pval)
    #g150 low 10.627289377289378 6.966467925466662e-08
    
    odds, pval = fisher_exact([[len(is_high_cnv_set), 
                                len(cnv_set)],
                               [len(is_high_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('g250 high', odds, pval)
    #g250 high 0.4606178230632379 0.4416190976760535

    odds, pval = fisher_exact([[len(is_low_cnv_set), 
                                len(cnv_set)],
                               [len(is_low_cnn_set), 
                                len(cnn_set)]], alternative='two-sided')
    print('g250 low', odds, pval)
    #g250 low 7.718157181571816 0.0002815512745766273
    
    
    ct = 0   
    to_print_set = set()
    for gene in gene_counter:
        
        
        hsize = len(gene_counter[gene]['high'])
        
        lsize = len(gene_counter[gene]['low'])

        # if (hsize > 1):
        #     print(gene, gene_counter[gene])
            
        # if (lsize > 1):
        #     ct+=1
        #     print(gene, gene_counter[gene])

        # if (hsize > 0) and (lsize == 0):
        #     ct+=1
        #     print(gene, gene_counter[gene])
        
        # if (hsize > 1) and (lsize == 0):
        #     ct+=1
        #     print(gene, gene_counter[gene])
        #     to_print_set.add(gene)
                
        # if (hsize == 0) and (lsize > 0):
        #     ct+=1
        #     print(gene, gene_counter[gene])
        
        if (hsize == 0) and (lsize > 1):
            ct+=1
            print(gene, gene_counter[gene])
            to_print_set.add(gene)
            
        # if (hsize > 0) and (lsize > 0):
        #     ct+=1
        #     print(gene, gene_counter[gene])
        #     to_print_set.add(gene)
        
    print(ct)
   