# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 12:47:18 2023

@author: pspea
"""

import random
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu

#DESEQ2 mrna for RSC4
cnv_list = [-0.633482966159613]
cnn_list = [-0.438083922, -0.585934693, -0.326769896]

U1, p = mannwhitneyu(cnv_list, cnn_list, method="exact")
ratio = np.median(cnv_list)/np.median(cnn_list) #1.446031078400574 0.5

print(ratio, p)

def fet_life(a,b,c,d):
        
    odds, p = fisher_exact([[a, b],
                               [c, d]], alternative='two-sided')
    
    ratio = round((a/b)/(c/d), 5)
    pval = round(p, 5)
    
    print(ratio, pval, p)
    



#change in number of significantly different RNA genes between generations?
sig_g250 = np.median([1635, 1126])
sig_g150 = np.median([255, 467])
insig_g250 = np.median([4236-1635, 4236-1126])
insig_g150 = np.median([4236-255, 4236-467])

fet_life(sig_g250, insig_g250, sig_g150, insig_g150)

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


#tEFF
total_genes = 4291
tEFF = 1027
uorfs = 526 #526 - May, SN - #615, uorfish - #361
overlap = (29+33+7+21) #May
two_way(total_genes, tEFF, uorfs, overlap, 1000) #Median exp: 126.0, Num obs: 90, ratio: 0.714, pval:1.0


#tEFF
total_genes = 4291
tEFF = 1027
uorfs = 615 #526 - May, SN - #615, uorfish - #361
#overlap = (29+33+7+21) #May
overlap = (83+33+55+21) #SN
two_way(total_genes, tEFF, uorfs, overlap, 10000) #Median exp: 147.0, Num obs: 192, ratio: 1.306, pval:0.0

#tEFF
total_genes = 4291
tEFF = 1027
uorfs = 361 #526 - May, SN - #615, uorfish - #361
#overlap = (29+33+7+21) #May
#overlap = (83+33+55+21) #SN
overlap = (60+7+55+21)
two_way(total_genes, tEFF, uorfs, overlap, 10000) #Median exp: 86.0, Num obs: 143, ratio: 1.663, pval:0.0

'''
# Evaluated against the 'only r1 teff results'
'''
#tEFF
total_genes = 4291
tEFF = 410
uorfs = 526 #526 - May, SN - #615, uorfish - #361
overlap = (9+14+6+4) #May
two_way(total_genes, tEFF, uorfs, overlap, 1000) #Median exp: 50.0, Num obs: 33, ratio: 0.66, pval:0.9986

#tEFF
total_genes = 4291
tEFF = 410
uorfs = 615 #526 - May, SN - #615, uorfish - #361
#overlap = (29+33+7+21) #May
overlap = (14+29+17+6) #SN
two_way(total_genes, tEFF, uorfs, overlap, 1000) #Median exp: 59.0, Num obs: 66, ratio: 1.119, pval:0.157

#tEFF
#NB - the uorfs here are score > 0.95, there are 446 in raw - but only 361 in the Unit set
total_genes = 4291
tEFF = 410
uorfs = 361 #526 - May, SN - #615, uorfish - #361
#overlap = (29+33+7+21) #May
#overlap = (83+33+55+21) #SN
overlap = (30+6+4+17)
two_way(total_genes, tEFF, uorfs, overlap, 1000) #Median exp: 34.0, Num obs: 57, ratio: 1.676, pval:0.0


'''
uorf and ssd1
'''
total_genes = 4291
ssd1_targets = 83 #ssd1_targets (bayne, both stress and not)
uorfs = 361 
overlap = (34)
two_way(total_genes, ssd1_targets, uorfs, overlap, 1000) #Median exp: 7.0, Num obs: 34, ratio: 4.857, pval:0.0

total_genes = 4291
ssd1_targets = 162 #ssd1_sig (having significantly FET more than expected)
uorfs = 361 
overlap = (34)
two_way(total_genes, ssd1_targets, uorfs, overlap, 10000) #Median exp: 14.0, Num obs: 34, ratio: 2.429, pval:0.0

total_genes = 4291
ssd1_targets = 1165 #ssd1_moe (having more than expected)
uorfs = 361 
overlap = (361)
two_way(total_genes, ssd1_targets, uorfs, overlap, 1000) #Median exp: 35.0, Num obs: 361, ratio: 10.314, pval:0.0

total_genes = 4291
ssd1_targets = 303 #ssd1_moe (>= 2)
uorfs = 361 
overlap = (52)
two_way(total_genes, ssd1_targets, uorfs, overlap, 1000) #Median exp: 35.0, Num obs: 361, ratio: 10.314, pval:0.0

total_genes = 4291
ssd1_targets = 1521 #ssd1 at least one
uorfs = 361 
overlap = (200)
two_way(total_genes, ssd1_targets, uorfs, overlap, 10000) #Median exp: 35.0, Num obs: 361, ratio: 10.314, pval:0.0


'''
uORF and 
'''
total_genes = 6000
ssd1_targets = 125 #ssd1_targets (bayne, both stress and not, not limited to UNIT)
uorfs = 820 # from May (not limited to UNIT)
overlap = (13)
two_way(total_genes, ssd1_targets, uorfs, overlap, 1000) #Median exp: 17.0, Num obs: 13, ratio: 0.765, pval:0.894

'''
uORF and 
'''
total_genes = 5275
ssd1_targets = 76 #ssd1_targets (bayne, both stress and not, not limited to UNIT)
uorfs = 426 # Unique to DGY
overlap = (6)
two_way(total_genes, ssd1_targets, uorfs, overlap, 1000) #Median exp: 17.0, Num obs: 13, ratio: 0.765, pval:0.894

total_genes = 5275
ssd1_targets = 217 #ssd1_motifs (sig positive)
uorfs = 426 # Unique to DGY
overlap = (35)
two_way(total_genes, ssd1_targets, uorfs, overlap, 10000) #Median exp: 17.0, Num obs: 35, ratio: 2.059, pval:0.0001

'''
#peff and protein complex
'''
total_genes = 4289
protein_complex = 487 #
sig_peff = 496 # Unique to DGY
overlap = (87)
two_way(total_genes, protein_complex, sig_peff, overlap, 10000) #Median exp: 17.0, Num obs: 35, ratio: 2.059, pval:0.0001




