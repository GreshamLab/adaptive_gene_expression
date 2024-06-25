# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 13:26:45 2024

@author: pspea
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:33:49 2024

#CDS_RPF_DGY1726_R1_tpm_qt
#CDS_RPF_DGY1726_R2_tpm_qt

@author: pspea
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
#teff_check 



infile_name = ('analyses/efficiency/ratio_df_Unit_v13.tab')    
ratio_df = pd.read_table(infile_name, index_col=0)
ratio_dict = ratio_df.to_dict('index') 

anc_rpf_rep1_str = ('CDS_RPF_DGY1657_R1_tpm')
anc_rpf_rep2_str = ('CDS_RPF_DGY1657_R2_tpm')
anc_rna_rep1_str = ('CDS_RNA_DGY1657_R1_tpm')
anc_rna_rep2_str = ('CDS_RNA_DGY1657_R2_tpm')

anc_rpf_median_str = ('CDS_RPF_DGY1657_median_tpm')
anc_rna_median_str = ('CDS_RNA_DGY1657_median_tpm')

sig_dict = {}

for gene in ratio_dict:
    
    anc_rpf_rep1 = ratio_dict[gene][anc_rpf_rep1_str]
    anc_rpf_rep2 = ratio_dict[gene][anc_rpf_rep2_str]
    anc_rna_rep1 = ratio_dict[gene][anc_rna_rep1_str]
    anc_rna_rep2 = ratio_dict[gene][anc_rna_rep2_str]
    
    if gene not in sig_dict:
        sig_dict[gene] = {anc_rpf_rep1_str: anc_rpf_rep1,
                          anc_rpf_rep2_str: anc_rpf_rep2,
                          anc_rna_rep1_str: anc_rna_rep1,
                          anc_rna_rep2_str: anc_rna_rep2}
    
    #for evo_strain in set(['DGY1726']):
    for evo_strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        evo_rpf_rep1_str = ('CDS_RPF_{}_R1_tpm').format(evo_strain)
        evo_rpf_rep2_str = ('CDS_RPF_{}_R2_tpm').format(evo_strain)
        evo_rna_rep1_str = ('CDS_RNA_{}_R1_tpm').format(evo_strain)
        evo_rna_rep2_str = ('CDS_RNA_{}_R2_tpm').format(evo_strain)
        #
        fetr_str_r1 = ('{}_fet_ratio_r1').format(evo_strain)      
        pval_str_r1 = ('{}_fet_pval_r1').format(evo_strain)
        fetr_str_r2 = ('{}_fet_ratio_r2').format(evo_strain)      
        pval_str_r2 = ('{}_fet_pval_r2').format(evo_strain)
        #
        evo_rpf_rep1 = ratio_dict[gene][evo_rpf_rep1_str]
        evo_rpf_rep2 = ratio_dict[gene][evo_rpf_rep2_str]
        evo_rna_rep1 = ratio_dict[gene][evo_rna_rep1_str]
        evo_rna_rep2 = ratio_dict[gene][evo_rna_rep2_str]
        #

        _odds, p1 = fisher_exact([[evo_rpf_rep1, 
                                    evo_rna_rep1],
                                   [anc_rpf_rep1, 
                                    anc_rna_rep1]], alternative='two-sided')
    
        fet_pval = round(p1,5)
        if evo_rna_rep1 != 0:
            evo_ratio = round(evo_rpf_rep1/evo_rna_rep1,3)
        else:
            evo_ratio = round(evo_rpf_rep1,3)
            
        if anc_rna_rep1 != 0:
            anc_ratio = round(anc_rpf_rep1/anc_rna_rep1,3)
        else:
            anc_ratio = round(anc_rpf_rep1,3)
            
        if anc_ratio != 0:
            fet_ratio_1 = (evo_ratio/anc_ratio)
        else:
            fet_ratio_1 = (evo_ratio)
        
        
        _odds, p2 = fisher_exact([[evo_rpf_rep2, 
                                    evo_rna_rep2],
                                   [anc_rpf_rep2, 
                                    anc_rna_rep2]], alternative='two-sided')
            
        fet_pval = round(p2,5)
        if evo_rna_rep2 != 0:
            evo_ratio = round(evo_rpf_rep2/evo_rna_rep2,3)
        else:
            evo_ratio = round(evo_rpf_rep2,3)
            
        if anc_rna_rep2 != 0:
            anc_ratio = round(anc_rpf_rep2/anc_rna_rep2,3)
        else:
            anc_ratio = round(anc_rpf_rep2,3)
            
        if anc_ratio != 0:
            fet_ratio_2 = (evo_ratio/anc_ratio)
        else:
            fet_ratio_2 = (evo_ratio)
            
        sig_dict[gene][evo_rpf_rep1_str] = evo_rpf_rep1
        sig_dict[gene][evo_rpf_rep2_str] = evo_rpf_rep2
        sig_dict[gene][evo_rna_rep1_str] = evo_rna_rep1
        sig_dict[gene][evo_rna_rep2_str] = evo_rna_rep2
        
        sig_dict[gene][fetr_str_r1] = fet_ratio_1
        sig_dict[gene][fetr_str_r2] = fet_ratio_2
        sig_dict[gene][pval_str_r1] = p1
        sig_dict[gene][pval_str_r2] = p2
        
        median_evo_rpf = np.median([evo_rpf_rep1, evo_rpf_rep2])
        median_anc_rpf = np.median([anc_rpf_rep1, anc_rpf_rep2])
        
        median_evo_rna = np.median([evo_rna_rep1, evo_rna_rep2])
        median_anc_rna = np.median([anc_rna_rep1, anc_rna_rep2])
        
        _odds, pm = fisher_exact([[median_evo_rpf, 
                                    median_evo_rna],
                                   [median_anc_rpf, 
                                    median_anc_rna]], alternative='two-sided')
            
        fet_pval = round(pm,5)
        if evo_rna_rep2 != 0:
            evo_ratio = round(median_evo_rpf/median_evo_rna,3)
        else:
            evo_ratio = round(median_evo_rpf,3)
            
        if anc_rna_rep2 != 0:
            anc_ratio = round(median_anc_rpf/median_anc_rna,3)
        else:
            anc_ratio = round(median_anc_rpf,3)
            
        if anc_ratio != 0:
            fet_ratio_median = (evo_ratio/anc_ratio)
        else:
            fet_ratio_median = (evo_ratio)
            
        evo_rpf_median_str = ('CDS_RPF_{}_median_tpm').format(evo_strain)
        evo_rna_median_str = ('CDS_RNA_{}_median_tpm').format(evo_strain)
        fetr_str_median = ('{}_fet_ratio_median').format(evo_strain)      
        pval_str_median = ('{}_fet_pval_median').format(evo_strain)
            
        sig_dict[gene][evo_rpf_median_str] = median_evo_rpf
        sig_dict[gene][anc_rpf_median_str] = median_anc_rpf
        sig_dict[gene][evo_rna_median_str] = median_evo_rna
        sig_dict[gene][anc_rna_median_str] = median_anc_rna        
        sig_dict[gene][fetr_str_median] = fet_ratio_median
        sig_dict[gene][pval_str_median] = pm

        
        
outfile_name = ('analyses/efficiency/Unit_v13_teff_results.tab')   
sig_df = pd.DataFrame.from_dict(sig_dict, orient='index')
sig_df.to_csv(outfile_name, sep = '\t') 