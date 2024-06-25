# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 11:01:22 2024

v0.12 - added expectation calculation

@author: pspea
"""
import pandas as pd
import numpy as np
import random

def populate_quantile_lists(selective_category, positive_scored_category, unit_object, names_dict, level):
    if level == 1:
        quantile_dict = {}
    
        for strain in names_dict:        
            if strain != 'DGY1657':
                for istype in names_dict[strain]:                    
                    sel_val_dict_pos = {}
                    sel_val_list_pos = []
                    
                    sel_val_dict_neg = {}
                    sel_val_list_neg = []
                    
                    for gene in unit_object:
                        
                        col_name = ('{strain}_{istype}_{selective_category}').format(
                            strain = strain, istype = istype, 
                            selective_category = selective_category)
        
                        sel_val =  unit_object[gene][col_name]
            
                        col_name = ('{strain}_{istype}_signal_to_noise').format(
                            strain = strain, istype = istype)
        
                        process =  unit_object[gene][col_name]
                            
                        if process:
                            
                            if sel_val not in sel_val_dict_pos:
                                sel_val_dict_pos[sel_val] = set()
                            sel_val_dict_pos[sel_val].add(gene)
                            
                            sel_val_list_pos.append(sel_val)
        

                            
                    if strain not in quantile_dict:
                        quantile_dict[strain] = {}
                        
                    if istype not in quantile_dict[strain]:
                        quantile_dict[strain][istype] = {'sel_val_dict_pos': sel_val_dict_pos,
                                                       'sel_val_list_pos': sel_val_list_pos,
                                                       'sel_val_dict_neg': sel_val_dict_neg,
                                                       'sel_val_list_neg': sel_val_list_neg}
        return(quantile_dict)
                        
                        
    if level == 2:
        quantile_dict = {}

        numer_denom_pair = {'ms':'rpf',
                            'rpf':'rna'}
        
        #for strain in set(['DGY1726']):
        for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):  
            for numer in numer_denom_pair:
                denom = numer_denom_pair[numer]
                
                pair = ('{numer}_{denom}').format(numer = numer, 
                                                  denom = denom)
                
                sel_val_dict_pos = {}
                sel_val_list_pos = []
                
                sel_val_dict_neg = {}
                sel_val_list_neg = []
                
                for gene in unit_object:

                    sel_col_name = ('{strain}_{numer}_{denom}_{selective_category}').format(
                        strain = strain,  numer = numer,  denom = denom, 
                        selective_category = selective_category)
                    
                    sel_val = unit_object[gene][sel_col_name]
                    

                    col_name = ('{strain}_{numer}_{denom}_signal_to_noise').format(
                        strain = strain, numer = numer,  denom = denom)
    
                    process =  unit_object[gene][col_name]
                        
                    if process:
                        if sel_val not in sel_val_dict_pos:
                            sel_val_dict_pos[sel_val] = set()
                        sel_val_dict_pos[sel_val].add(gene)
                        
                        sel_val_list_pos.append(sel_val)
    
                    
                        
                if strain not in quantile_dict:
                    quantile_dict[strain] = {}
                    
                if pair not in quantile_dict[strain]:
                    quantile_dict[strain][pair] = {'sel_val_dict_pos': sel_val_dict_pos,
                                                   'sel_val_list_pos': sel_val_list_pos,
                                                   'sel_val_dict_neg': sel_val_dict_neg,
                                                   'sel_val_list_neg': sel_val_list_neg}
        return(quantile_dict)


def half_quant(index_val, needed_genes, temp_set, sum_list, sum_dict, runmode):

    if runmode == 'up':
        sum_list.sort()
    
        while (len(temp_set) <= needed_genes):
            
            index = len([x for x in sum_list if x <= index_val])
            
            if index < len(sum_list):
                index_val = sum_list[index]
                genes_at_val = sum_dict[index_val]
                temp_set = temp_set.union(genes_at_val) 
                #print(index, index_val, len(temp_set))
            else:
                return(temp_set)
            
    if runmode == 'down':
        sum_list.sort(reverse=True)
    
        while (len(temp_set) <= needed_genes):
            
            index = len([x for x in sum_list if x >= index_val])
            
            if index < len(sum_list):
                index_val = sum_list[index]
                genes_at_val = sum_dict[index_val]
                temp_set = temp_set.union(genes_at_val) 
                #print(index, index_val, len(temp_set))
            else:
                return(temp_set)
            
    return(temp_set)

def directed_half_quant(index_val, needed_genes, temp_set, sum_list, sum_dict, runmode):

    if runmode == 'up':
        sum_list.sort()
    
        while (len(temp_set) <= needed_genes):
            
            index = len([x for x in sum_list if x <= index_val])
            
            if index < len(sum_list):
                index_val = sum_list[index]
                genes_at_val = sum_dict[index_val]
                temp_set = temp_set.union(genes_at_val) 
                #print(index, index_val, len(temp_set))
            else:
                return(temp_set)
            
    if runmode == 'down':
        sum_list.sort(reverse=True)
    
        while (len(temp_set) <= needed_genes):
            
            index = len([x for x in sum_list if x >= index_val])
            
            if index < len(sum_list):
                index_val = sum_list[index]
                genes_at_val = sum_dict[index_val]
                temp_set = temp_set.union(genes_at_val) 
                #print(index, index_val, len(temp_set))
            else:
                return(temp_set)
            
    return(temp_set)

def score_features(positive_scored_category, negative_scored_category, selective_category, quantile_dict, unit_object, level, pct):
    
    if level == 1:
        
        one_side = round((len(unit_object)/pct)/2) 
        
        for gene in unit_object:
            for strain in quantile_dict:
                if strain != 'DGY1657':
                    for istype in quantile_dict[strain]:       
                                                
                        sel_col_name = ('{strain}_{istype}_{selective_category}').format(
                            strain = strain, istype = istype, selective_category = selective_category)
                        
                        index_val =  unit_object[gene][sel_col_name]

                        scored_col_name = ('{strain}_{istype}_{positive_scored_category}').format(
                            strain = strain, istype = istype, positive_scored_category = positive_scored_category)
                        
                        score_val = unit_object[gene][scored_col_name] #1

                        sel_val_dict = quantile_dict[strain][istype]['sel_val_dict_pos']
                        sel_val_list = quantile_dict[strain][istype]['sel_val_list_pos']
                        
                        if score_val < 1:                            
                            scored_col_name = ('{strain}_{istype}_{negative_scored_category}').format(
                                strain = strain, istype = istype, negative_scored_category = negative_scored_category)
                
                        #print(gene)
                        quantile_set = set()
                        
                        quantile_set = half_quant(index_val, one_side, quantile_set, sel_val_list, sel_val_dict, 'up')
                        
                        quantile_set = half_quant(index_val, 2*one_side, quantile_set, sel_val_list, sel_val_dict, 'down')
                
                        if len(quantile_set) < (2*one_side):
                            quantile_set = half_quant(index_val, 2*one_side, quantile_set, sel_val_list, sel_val_dict, 'up')
                            
                        if len(quantile_set) < (2*one_side):
                            quantile_set = half_quant(index_val, 2*one_side, quantile_set, sel_val_list, sel_val_dict, 'down')
                        
                        col_name = ('{strain}_{istype}_quantile_set_number').format(
                            strain = strain, istype = istype)
                        unit_object[gene][col_name] = len(quantile_set)
                        
                        score_list = []
                        
                        for sim in range(int(one_side)): 
                            for quantile_gene in quantile_set:
                                if random.random() <= 0.01:
                                    #print(quantile_gene, strain, istype)
                                    score_list.append(unit_object[quantile_gene][scored_col_name])
                                                        
                        score = unit_object[gene][scored_col_name]
                        
                        if len(score_list) > 1:
                            fdr = len([x for x in score_list if x >= score])/len(score_list)
                        else:
                            print(gene, istype, scored_col_name)
                            print(len(score_list))
                            

                        
                        fdr_col_name = ('{strain}_{istype}_fdr').format(
                            strain = strain, istype = istype)
                        
                        unit_object[gene][fdr_col_name] = fdr
                        
                        if fdr <= 0.05:
                            print(gene, strain, istype, score, fdr)

                            
    if level == 2:
        #one_side = round((len(unit_object)/10)/2) #215 (TP: 55, FP: 280, FN: 170: TN: 3786), (Sensitivity [195/(195+254)]: 0.24, Specificity [3813]: 0.93 )
        one_side = round((len(unit_object)/pct)/2)
        numer_denom_pair = set(['ms_rpf', 'rpf_rna'])

        for gene in unit_object:
            #print(gene)
            for pair in numer_denom_pair:
                numer, denom = pair.split('_')
                #print(pair)
                for strain in quantile_dict:
                    if strain != 'DGY1657':
                        #print(strain)
                        if pair in quantile_dict[strain]:
                            
                            sel_col_name = ('{strain}_{numer}_{denom}_{selective_category}').format(
                                strain = strain,  numer = numer, denom = denom,
                                selective_category = selective_category)
                        
                            index_val =  unit_object[gene][sel_col_name]
                            #print(index_val, sel_col_name)
                            
                            scored_col_name = ('{strain}_{numer}_{denom}_{positive_scored_category}').format(
                                strain = strain, numer = numer, denom = denom,  
                                positive_scored_category = positive_scored_category)
                            
                            score_val = unit_object[gene][scored_col_name]
                            
                            sel_val_dict = quantile_dict[strain][pair]['sel_val_dict_pos']
                            sel_val_list = quantile_dict[strain][pair]['sel_val_list_pos']
                            

                            
                            quantile_set = set()
                            
                            quantile_set = half_quant(index_val, one_side, quantile_set, sel_val_list, sel_val_dict, 'up')
                            
                            quantile_set = half_quant(index_val, 2*one_side, quantile_set, sel_val_list, sel_val_dict, 'down')
                    
                            if len(quantile_set) < (2*one_side):
                                quantile_set = half_quant(index_val, 2*one_side, quantile_set, sel_val_list, sel_val_dict, 'up')
                            
                            score_list = []
                            
                            for sim in range(int(one_side)): 
                                for quantile_gene in quantile_set:
                                    if quantile_gene != gene:
                                        if random.random() <= 0.1:
                                            score_list.append(unit_object[quantile_gene][scored_col_name])
                                
                            score = unit_object[gene][scored_col_name]
                            #print(score, len(score_list), index_val)
                            
                            fdr = len([x for x in score_list if x >= score])/len(score_list)
                            
                                                       
                            fdr_col_name = ('{strain}_{numer}_{denom}_fdr').format(
                                strain = strain, numer = numer, denom = denom)
                            
                            unit_object[gene][fdr_col_name] = fdr
                            
                            if fdr <= 0.05:
                                print(gene, strain, pair, score, fdr)

    return(unit_object)

def build_names_dict(ratio_df):
    names_dict = {}
    
    '''
    for each molecule (istype), for each strain, 
    '''
    for istype in ['RNA', 'RPF', 'ms']:
        for strain in set(['DGY1657', 'DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        #for strain in set(['DGY1657', 'DGY1726']):
        #for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']): 
            
            if strain not in names_dict:
                names_dict[strain] = {'rna':{}, 'rpf':{}, 'ms':{}}
                
            if istype != 'ms':
                for rep in range(1,3):
                    colname = ("CDS_{istype}_{strain}_R{rep}_tpm_pt").format(
                        istype = istype, strain = strain, rep  = rep)
                    
                    lcase_istype = istype.lower()
                    
                    names_dict[strain][lcase_istype][rep] = {'colname': colname}
            else:
                for rep in range(1,6):
                    
                    colname = ("ms_{strain}_rep{rep}_pt").format(
                        strain = strain, rep  = rep)
                                        
                    names_dict[strain][istype][rep] = {'colname': colname}
                    
    for strain in names_dict:
        for istype in names_dict[strain]:
            for rep in names_dict[strain][istype]:
                colname = names_dict[strain][istype][rep]['colname']
                
                exp_list = ratio_df[colname].tolist()
                #tpm_pt_df['CDS_RNA_DGY1657_R1_tpm_pt'].tolist()
                
                names_dict[strain][istype][rep]['exp_list'] = exp_list
                
    return(names_dict)

def build_unit_object(ratio_df):
        
    unit_object = ratio_df.to_dict('index')
    
    ### base unit for strain and istype
    for strain in names_dict:
        for istype in names_dict[strain]:
            for gene in unit_object:
                
                #for each gene we find expression abundance (and rank
                gene_exp_list = []
                gene_rank_list = []
                
                for rep in names_dict[strain][istype]:
                    col_name = ('{strain}_{istype}_{rep}_rank').format(strain = strain, istype = istype, rep = rep)
                    colname = names_dict[strain][istype][rep]['colname']
                    
                    gene_val = unit_object[gene][colname]
                    gene_exp_list.append(gene_val)
                    
                    
                    
                    exp_list = names_dict[strain][istype][rep]['exp_list']
                    rank = len([x for x in exp_list if x <= gene_val])
                    unit_object[gene][col_name] = (rank)     
                
                for rep in names_dict[strain][istype]:
                    colname = names_dict[strain][istype][rep]['colname']
                    
                    gene_val = unit_object[gene][colname]
                    gene_exp_list.append(gene_val)
                    
                    exp_list = names_dict[strain][istype][rep]['exp_list']
                    rank = len([x for x in exp_list if x <= gene_val])
                    gene_rank_list.append(rank)
                    
                #note that these are all post-QT transform
                col_name = ('{strain}_{istype}_pt_median').format(strain = strain, istype = istype)
                unit_object[gene][col_name] = np.median(gene_exp_list)
                
                col_name = ('{strain}_{istype}_pt_std').format(strain = strain, istype = istype)
                unit_object[gene][col_name] = np.std(gene_exp_list)
                

    #first degree ration evo / anc using same istype
    for strain in names_dict:
        for istype in names_dict[strain]:
            for gene in unit_object:
                
                col_name = ('{strain}_{istype}_pt_median').format(strain = strain, istype = istype)
                evo_median = unit_object[gene][col_name]
                
                col_name = ('{strain}_{istype}_pt_std').format(strain = strain, istype = istype)
                evo_std = unit_object[gene][col_name]

                
                if strain != 'DGY1657':
                    col_name = ('DGY1657_{istype}_pt_median').format(istype = istype)
                    anc_median = unit_object[gene][col_name]
                    
                    col_name = ('DGY1657_{istype}_pt_std').format(istype = istype)
                    anc_std = unit_object[gene][col_name]
                    

                    
                    col_name = ('{strain}_{istype}_sum_of_pt_medians').format(strain = strain, istype = istype)
                    sum_of_pt_medians = evo_median + anc_median
                    unit_object[gene][col_name] = sum_of_pt_medians
                    
                    new_col_name = ('{strain}_{istype}_rel_pt_evo_over_anc').format(
                        strain = strain, istype = istype)
                    ratio_level_score = 1
                    if (evo_median != 0) and (anc_median != 0):
                        ratio_level_score = (evo_median/anc_median)
                    unit_object[gene][new_col_name] = ratio_level_score
                    
                    new_col_name = ('{strain}_{istype}_rel_pt_anc_over_evo').format(
                        strain = strain, istype = istype)
                    ratio_level_score = 1
                    if (evo_median != 0) and (anc_median != 0):
                        ratio_level_score = (anc_median/evo_median)
                    unit_object[gene][new_col_name] = ratio_level_score
                    
                    col_name = ('{strain}_{istype}_sum_of_pt_std').format(strain = strain, istype = istype)                    
                    unit_object[gene][col_name] = abs(evo_std + anc_std)
                    
                    col_name = ('{strain}_{istype}_signal_to_noise').format(strain = strain, istype = istype)   
                    process = False
                    if abs(evo_median - anc_median) > abs(evo_std + anc_std):
                        process = True
                        
                    unit_object[gene][col_name] = process
                    


    return(unit_object)

def score_object(unit_object, pct):
    
    level = 1                           
    selective_category = 'sum_of_pt_medians'
    positive_scored_category = 'rel_pt_evo_over_anc'
    negative_scored_category = 'rel_pt_anc_over_evo'
    
    quantile_dict = populate_quantile_lists(selective_category, positive_scored_category, unit_object, names_dict, level)
    
    unit_object = score_features(positive_scored_category, negative_scored_category, selective_category, quantile_dict, unit_object, level, pct)
    
    return(unit_object)


def build_unit_object_level_2(unit_object, names_dict, pct):   
    
    #2 level analysis 
    numer_denom_pair = set(['ms_rpf', 'rpf_rna'])
    
    for gene in unit_object:
        
        for pair in numer_denom_pair:
            numer, denom = pair.split('_')
            # pull values
            col_name = ('DGY1657_{istype}_pt_median').format(istype = denom)
            anc_denom_median = unit_object[gene][col_name]
            
            col_name = ('DGY1657_{istype}_pt_std').format(istype = denom)
            anc_denom_std = unit_object[gene][col_name]
            
            col_name = ('DGY1657_{istype}_pt_median').format(istype = numer)
            anc_numer_median = unit_object[gene][col_name]
            
            col_name = ('DGY1657_{istype}_pt_std').format(istype = numer)
            anc_numer_std = unit_object[gene][col_name]
            
            
            #for strain in set(['DGY1726']):
            for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
                #pull values
                col_name = ('{strain}_{istype}_pt_median').format(strain = strain, istype = denom)
                evo_denom_median = unit_object[gene][col_name]
                
                col_name = ('{strain}_{istype}_pt_std').format(strain = strain, istype = denom)
                evo_denom_std = unit_object[gene][col_name]
                                
                col_name = ('{strain}_{istype}_pt_median').format(strain = strain, istype = numer)
                evo_numer_median = unit_object[gene][col_name]
                
                col_name = ('{strain}_{istype}_pt_std').format(strain = strain, istype = numer)
                evo_numer_std = unit_object[gene][col_name]
                                                                
                new_col_name = ('{strain}_{numer}_{denom}_rel_pt_sum_of_medians').format(
                    strain = strain, numer = numer, denom = denom)
                rel_pt_sum_of_medians = evo_numer_median + anc_numer_median + evo_denom_median + anc_denom_median
                unit_object[gene][new_col_name] = rel_pt_sum_of_medians
                
                new_col_name = ('{strain}_{numer}_{denom}_rel_pt_diff_of_ratios').format(
                    strain = strain, numer = numer, denom = denom)
                
                rel_pt_diff_of_ratios = 0
                if (evo_denom_median != 0) and (evo_numer_median != 0):
                    evo_ratio = (evo_numer_median/evo_denom_median)
                    if (anc_denom_median != 0) and (anc_numer_median != 0):
                        anc_ratio = (anc_numer_median/anc_denom_median)
                        if anc_ratio != 0:
                            rel_pt_diff_of_ratios = evo_ratio - anc_ratio
                unit_object[gene][new_col_name] = rel_pt_diff_of_ratios
                                

                #signal_to_noise
                new_col_name = ('{strain}_{numer}_{denom}_signal_to_noise').format(
                    strain = strain, numer = numer, denom = denom)
                process = False
                
                if abs(evo_numer_median - anc_numer_median) > abs(evo_numer_std + anc_numer_std):
                    process = True
                if abs(evo_denom_median - anc_denom_median) > abs(evo_denom_std + anc_denom_std):
                    process = True
                
                unit_object[gene][new_col_name] = process
                
                new_col_name = ('{strain}_{numer}_{denom}_rel_pt_evo_over_anc').format(
                    strain = strain, numer = numer, denom = denom)
                
                two_level_score = 1
                if (evo_denom_median != 0) and (evo_numer_median != 0):
                    evo_ratio = (evo_numer_median/evo_denom_median)
                    if (anc_denom_median != 0) and (anc_numer_median != 0):
                        anc_ratio = (anc_numer_median/anc_denom_median)
                        if anc_ratio != 0:
                            two_level_score = (evo_ratio / anc_ratio)
                
                unit_object[gene][new_col_name] = two_level_score
                
                new_col_name = ('{strain}_{numer}_{denom}_rel_pt_anc_over_evo').format(
                    strain = strain, numer = numer, denom = denom)
                
                two_level_score = 1
                if (evo_denom_median != 0) and (evo_numer_median != 0):
                    evo_ratio = (evo_numer_median/evo_denom_median)
                    if (anc_denom_median != 0) and (anc_numer_median != 0):
                        anc_ratio = (anc_numer_median/anc_denom_median)
                        if anc_ratio != 0:
                            two_level_score = (anc_ratio / evo_ratio)
                            
                unit_object[gene][new_col_name] = two_level_score
                

    level = 2                           
    #selective_category = 'rel_pt_numer_diff'
    selective_category = 'rel_pt_diff_of_ratios'
    positive_scored_category = 'rel_pt_evo_over_anc'
    negative_scored_category = 'rel_pt_anc_over_evo'
    
    #selective_category = 'pt_diff_of_medians'
    quantile_dict = populate_quantile_lists(selective_category, positive_scored_category, unit_object, names_dict, level)
    
    unit_object = score_features(positive_scored_category, negative_scored_category, selective_category, quantile_dict, unit_object, level, pct)
    #next try pt_sub as scored_category
    return(unit_object)


infile_name = ('analyses/efficiency/ratio_df_Unit_v13.tab')    
ratio_df = pd.read_table(infile_name, index_col=0)
names_dict = build_names_dict(ratio_df)
# ratio_df = pd.DataFrame.from_dict(ratio_dict, orient='index')
# ratio_df.to_csv(outfile_name, sep = '\t') 

for pct in set([10]):
#for pct in set([5,10,15,20,25,30,35,40,45,50]):
    unit_object = build_unit_object(ratio_df)
    unit_object = score_object(unit_object, pct)

    outfile_name = ('analyses/efficiency/unit_object_level_1v13_Unit_{pct}pct.tab').format(pct=pct)
    unit_object_df = pd.DataFrame.from_dict(unit_object, orient='index')
    unit_object_df.to_csv(outfile_name, sep = '\t')

    unit_object_name = ('analyses/efficiency/unit_object_level_1v13_Unit_{pct}pct.tab').format(pct=pct)
    unit_object_df = pd.read_csv(unit_object_name,
                                  sep='\t', index_col = 0)
    unit_object = unit_object_df.to_dict('index')
    
    unit_object = build_unit_object_level_2(unit_object, names_dict, pct)
    
    outfile_name = ('analyses/efficiency/unit_object_level_2v13v16_Unit_{pct}pct.tab').format(pct=pct)    
    unit_object_df = pd.DataFrame.from_dict(unit_object, orient='index')
    unit_object_df.to_csv(outfile_name, sep = '\t') 
    
#
