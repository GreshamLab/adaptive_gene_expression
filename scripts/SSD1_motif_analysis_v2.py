# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 20:59:02 2023

@author: pspea

_v2 - now generates a bedfile of SSD1 motif hits.
_v3 - 
    _x_ updated for both strands
    _x_ output 'at random' control for bedtools reldist
    
    
### limited analysis to the SN defined TL
"""

import numpy as np
from scipy.stats import fisher_exact

import random

chromo_dict = {}

call_sign = {'W': '+',
             'C': '-'}

n_list = ['A', 'T', 'C', 'G']
y_list = ['C', 'T']

def rev_comp(seq):
    
    compliment_dict = {'A':'T',
                       'T':'A',
                       'C':'G',
                       'G':'C'}
    
    rev_seq = seq[::-1]
    
    rev_comp_seq = ''
    
    for nt in rev_seq:
        if nt in compliment_dict:
            rev_comp_seq += compliment_dict[nt]
        else:
            rev_comp_seq += nt

    return(rev_comp_seq)

def calc_TL_fasta():
    tl_dict = {}
    tl_coord_dict = {}
    
    '''
    #tab seperated
    chrVII	2778	2789	1	YGL263W	ATCGCCATAGGC
    chrVII	2781	2789	1	YGL263W	GCCATAGGC
    '''
    tl_file = open('analyses/ssd1/TL_from_SGD.tsv')
    
    for line in tl_file:
        line = line.strip()
        chromo, start, stop, strand, gene, seq = line.split('\t')
        start = int(start)
        stop = int(stop)
        
        if gene not in tl_dict:
            tl_dict[gene] = seq
            
            sign = call_sign[gene.split('-')[0][-1]]
            
            if gene not in tl_coord_dict:
                tl_coord_dict[gene] = {'chromo': chromo,
                    'start': start,
                    'stop': stop,
                    'sign': sign, 
                    'ssd1_hits': {}}
            
        else:
            if len(seq) > len(tl_dict[gene]):
                tl_dict[gene] = seq
            
                tl_coord_dict[gene] = {'chromo': chromo,
                    'start': start,
                    'stop': stop,
                    'sign': sign, 
                    'ssd1_hits': {}}
                    
                    
    tl_file.close()  

    '''
    # {"source_attrib": "AWN", "MASK_OUT_ORFS": ["YOR302W"], "genome_files": ["/usr1/home/JUMBO/ANNOTATIONS/SGD_saccharomyces_cerevisiae.gff"], "positive_strand_end_peaks": ["/usr1/home/JUMBO/pA_bedGraph/peak_calls/total-Sc_rep1_pA_pos.bedGraph-putative_peaks-outlier_adjusted.bedGraph", "/usr1/home/JUMBO/pA_bedGraph/peak_calls/total-Sc_rep2_pA_pos.bedGraph-putative_peaks-outlier_adjusted.bedGraph"], "EFFECTIVE_LENGTH_CUTOFF": 1000, "prefix": "saccharomyces_cerevisiae", "output_dir": "/usr1/home/JUMBO/flanking_caller/20160712_generated/", "positive_strand_start_peaks": ["/usr1/home/JUMBO/TSS_bedGraph/peak_calls/total-Scer-hap-rep1_pos.bedGraph-putative_peaks-outlier_adjusted.bedGraph", "/usr1/home/JUMBO/TSS_bedGraph/peak_calls/total-Scer-hap-rep2_pos.bedGraph-putative_peaks-outlier_adjusted.bedGraph"], "negative_strand_end_peaks": ["/usr1/home/JUMBO/pA_bedGraph/peak_calls/total-Sc_rep1_pA_neg.bedGraph-putative_peaks-outlier_adjusted.bedGraph", "/usr1/home/JUMBO/pA_bedGraph/peak_calls/total-Sc_rep2_pA_neg.bedGraph-putative_peaks-outlier_adjusted.bedGraph"], "negative_strand_start_peaks": ["/usr1/home/JUMBO/TSS_bedGraph/peak_calls/total-Scer-hap-rep1_neg.bedGraph-putative_peaks-outlier_adjusted.bedGraph", "/usr1/home/JUMBO/TSS_bedGraph/peak_calls/total-Scer-hap-rep2_neg.bedGraph-putative_peaks-outlier_adjusted.bedGraph"], "exclusion_genome_files": [], "WINSORIZATION": 0.95}
    chrI	AWN	mRNA	136914	137510	.	+	.	ID=YAL008W_mRNA;PARENT=YAL008W
    chrI	AWN	gene	136914	137510	.	+	.	ID=YAL008W
    chrI	AWN	five_prime_UTR	136874	136913	.	+	.	PARENT=YAL008W_mRNA
    '''
        
    tl_file = open('metadata/McManus_2018_saccharomyces_cerevisiae.gff')
    tl_bed_file_name = ('analyses/ssd1/TL_from_SN.bed')
    tl_bed_file = open(tl_bed_file_name, 'w')
    
    for line in tl_file:
        line = line.strip()
        
        if 'five_prime_UTR' in line:
            chromo = line.split('\t')[0]
            start = int(line.split('\t')[3])
            stop = int(line.split('\t')[4])
            sign = line.split('\t')[6]
            
            gene = line.split('\t')[8].split('=')[1].split('_')[0]
            
            if gene[0] == 'Y':
            
                outline = ('{chromo}\t{start}\t{stop}\t{gene}\t1\t{sign}\n').format(
                    chromo = chromo, start = start, stop = stop, gene = gene, sign = sign)
                
                if gene not in tl_dict:
                    tl_bed_file.write(outline)
                    
                if gene not in tl_coord_dict:
                    tl_coord_dict[gene] = {'chromo': chromo,
                        'start': start,
                        'stop': stop,
                        'sign': sign, 
                        'ssd1_hits': {}}

                if gene in tl_dict:
                    if abs(stop-start) > len(tl_dict[gene]):
                        tl_bed_file.write(outline)                    
                        tl_coord_dict[gene] = {'chromo': chromo,
                                                'start': start,
                                                'stop': stop,
                                                'sign': sign, 
                                                'ssd1_hits': {}}
                        
    tl_file.close()
    tl_bed_file.close()
                    
    '''
    #fasta format
    >YGR131W::chrVII:754206-754725(+)
    AGCAATTGAAAGAACTGTGATTTATTTCCGCTTGTTCGAAATTATTGATGTTTAGCACTTTGCAGTAGCGACAATACAATATATGTGCTTTTAGTGCTGGGATAGTTCGTAGCTCCATTTCGGGGCGCTTGTTACATTTATTGTATATGCGCGGATGTGGCACATGCTGTTGAGATCTCACTCCTTTGGTATCTCTTTCCTGCGCCGCATTGTGCCGGCAGAATGTCGCGCTTGTATTCTCATGAACTTTTCCTCTTTACGAACCCTTTGGCGGCATGCCGTTTAAAATCTGTTGAAGATTTCCTTTACGAACAATGAGCAATGTTTTGCACAGGCAGGTGGGAAGTAGGGCCTATCGCGCCTTGGATGCAGATATAAGTATAAATATAAATTATAATAATTGGCTGTATCAGTAAATCCTTCTTGCGATGGGAGGAAGCACGATAGAGTATGTTAAGCTTTTGAGAGGCTTCATATTCATTGGAATTTTAAATAACAATAAAGCAACAACAATAATAA
    '''
    tl_file = open('analyses/ssd1/TL_from_SN.fa')
    
    for line in tl_file:
        line = line.strip()
        
        if line[0] == '>':
            gene = line.split('>')[1].split(':')[0]
        else:
            tl_dict[gene] = line
            
    tl_file.close()  
    
    tl_fasta_file_name = ('analyses/ssd1/TL_from_SGD_SN.fa')
    tl_fasta_file = open(tl_fasta_file_name, 'w')
                
    for gene in tl_dict:
        outline = ('>{gene}\n{seq}\n').format(
            gene = gene,
            seq = tl_dict[gene])
        tl_fasta_file.write(outline)
        
    tl_fasta_file.close()
    
    total_hit = 0
    total_length = 0
    hits_per_gene = {}
    
    missing_gene_set = set()

    for gene in tl_dict:
        if gene not in tl_coord_dict:
            missing_gene_set.add(gene)

        if gene in tl_coord_dict:
            #print('missing gene', gene)

            seq = tl_dict[gene]
            
            total_length += len(seq)
            
            if gene not in hits_per_gene:
                hits_per_gene[gene] = 0
                            
            for nt2 in n_list:
                for nt3 in y_list:
                    for nt6 in n_list:
                        for nt7 in y_list:
                            pos_search = ('C{nt2}{nt3}TC{nt6}{nt7}T').format(
                                nt2=nt2, nt3=nt3, nt6=nt6, nt7=nt7)
                            pos_count = seq.count(pos_search)
                            
                            if pos_count > 0:
                                start = 0
                                results = seq.split(pos_search)
                            
                                for each in results[:-1]:
                                    start+=len(each)

                                    tl_coord_dict[gene]['ssd1_hits'][start] = pos_search
                                    start+=len(pos_search)
                            
                            rev_search = rev_comp(pos_search)
                            rev_count = seq.count(rev_search)
                            
                            if rev_count > 0:
                                start = 0
                                results = seq.split(rev_search)
                            
                                for each in results[:-1]:
                                    start+=len(each)

                                    tl_coord_dict[gene]['ssd1_hits'][start] = rev_search
                                    start+=len(rev_search)
                            
                            count = (pos_count + rev_count)
                            hits_per_gene[gene] += count
                            total_hit += count



    hits_list = []                    
    for gene in hits_per_gene:
        hits_list.append(hits_per_gene[gene])
        
    print(np.median(hits_list))
    print(np.mean(hits_list))
        
    background_density = total_hit/total_length
    #genome wide per nucleotide rate: 0.00238590192991225674
    #TL per nucleotide rate: 0.002911805177351201
    
    sig_set = set()                  
      
    SSD1_hits_file = open('analyses/ssd1/_SSD1_hits_in_TL_from_SGD_SN.txt', 'w')
    
    outline = ('#gene\thits\tlength\tobserved_density\texp_hits\tpval\n')
    
    SSD1_hits_file.write(outline)
    
    gene_list = []

    ssd1_bed_file_name = ('analyses/ssd1/_SSD1_in_TL_from_SGD_SN.bed')
    ssd1_bed_file = open(ssd1_bed_file_name, 'w')
    
    print('missing_gene_set', len(missing_gene_set))

    for gene in tl_coord_dict:
        ssd1_hits = tl_coord_dict[gene]['ssd1_hits']

        if len(ssd1_hits) > 0:
            ct = 0
            chromo = tl_coord_dict[gene]['chromo']
            sign = tl_coord_dict[gene]['sign']
            start = tl_coord_dict[gene]['start']
            stop = tl_coord_dict[gene]['stop']

            for ssd1_nt in ssd1_hits:
                seq = ssd1_hits[ssd1_nt]
                name = ('{gene}_{seq}_{ct}').format(
                    gene = gene, seq = seq, ct = ssd1_nt)

                if sign == '+': 
                    ssd1_start = start + ssd1_nt
                    ssd1_stop = ssd1_start + len(seq)
                else:
                    ssd1_stop = stop - ssd1_nt 
                    ssd1_start = ssd1_stop - len(seq) 

                outline = ('{chromo}\t{ssd1_start}\t{ssd1_stop}\t{name}\t0\t{sign}\n').format(
                    chromo = chromo, ssd1_start = ssd1_start, 
                    ssd1_stop = ssd1_stop, name = name, sign = sign)

                ssd1_bed_file.write(outline)

                # outline = ('{chromo}\t{ssd1_start}\t{ssd1_stop}\t{name}_support\t0\t{sign}\n').format(
                #     chromo = chromo, ssd1_start = start, 
                #     ssd1_stop = stop, name = name, sign = sign)

                # ssd1_bed_file.write(outline)
                #print(outline)

                ct+=1
    ssd1_bed_file.close()         
    
    rK_bed_file_name = ('analyses/ssd1/_SSD1_random_random_control_in_TL_from_SGD_SN.bed')
    rK_bed_file = open(rK_bed_file_name, 'w')
        
    for gene in tl_coord_dict:
        ssd1_hits = tl_coord_dict[gene]['ssd1_hits']

        if len(ssd1_hits) > 0:
            
            gene_length = abs(tl_coord_dict[gene]['stop']-tl_coord_dict[gene]['start'])
            
            random_set = set()
            for potential_gene in tl_coord_dict:
                if potential_gene != gene:
                    potential_length = abs(tl_coord_dict[potential_gene]['stop']-tl_coord_dict[potential_gene]['start'])
                    
                    if (potential_length >= (gene_length*0.5)) and (potential_length >= (gene_length*1.5)):
                        random_set.add(potential_gene)
                                
            ssd1_hits = tl_coord_dict[gene]['ssd1_hits']


            if len(random_set) != 0:
                random_gene = random.sample(list(random_set),k=1)[0]
            else:
                print('no near match', gene, gene_length)
                
            ct = 0
            chromo = tl_coord_dict[random_gene]['chromo']
            sign = tl_coord_dict[random_gene]['sign']
            start = tl_coord_dict[random_gene]['start']
            stop = tl_coord_dict[random_gene]['stop']

            for ssd1_nt in ssd1_hits:
                seq = ssd1_hits[ssd1_nt]
                name = ('{gene}_{seq}_{ct}').format(
                    gene = random_gene, seq = seq, ct = ssd1_nt)

                if sign == '+': 
                    ssd1_start = start + ssd1_nt
                    ssd1_stop = ssd1_start + len(seq)
                else:
                    ssd1_stop = stop - ssd1_nt 
                    ssd1_start = ssd1_stop - len(seq) 

                outline = ('{chromo}\t{ssd1_start}\t{ssd1_stop}\t{name}\t0\t{sign}\n').format(
                    chromo = chromo, ssd1_start = ssd1_start, 
                    ssd1_stop = ssd1_stop, name = name, sign = sign)

                #rK_bed_file.write(outline)

                # outline = ('{chromo}\t{ssd1_start}\t{ssd1_stop}\t{name}_support\t0\t{sign}\n').format(
                #     chromo = chromo, ssd1_start = start, 
                #     ssd1_stop = stop, name = name, sign = sign)

                # ssd1_bed_file.write(outline)
                #print(outline)

                ct+=1
    rK_bed_file.close()
    
    for gene in hits_per_gene:
        length = len(tl_dict[gene])
        
        hits = hits_per_gene[gene]
        
        odds, pval = fisher_exact([[(hits), 
                                    (length)],
                                   [total_hit, total_length]], alternative='two-sided')
        
        gene_list.append(gene)
        
        if pval < 0.05:
            print(gene, hits, length)
            sig_set.add(gene)
            
        outline = ('{gene}\t{hits}\t{length}\t{density}\t{exp_hits}\t{pval}\n').format(
            gene = gene, hits = hits, length = length, 
            density = round(hits/length, 3),
            exp_hits = round(background_density*length),
            pval = pval
            )
        
        SSD1_hits_file.write(outline)
        
    SSD1_hits_file.close()
        
    
    
        
        
        
calc_TL_fasta()
        
        
    
    

def calc_genome_fasta():
        

    infast_file = open("C:/Gresham/genomes/SGD/S288C_reference_sequence_R64-2-1_20150113.fsa")
    
    
    
    for line in infast_file:
        '''
        >ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]
        CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACA
        CATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTT
        ACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC
        CACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCTCCAATTACCCATATC
        CAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCACCATAC
        TGTTCTTCTACCCACCATATTGAAACGCTAACAAATGATCGTAAATAACACACACGTGCT
        TACCCTACCACTTTATACCACCACCACATGCCATACTCACCCTCACTTGTATACTGATTT
        '''
        line = line.strip()
        
        if line[0] == '>':
            print(line)
            if 'chromosome=' in line:
                chromo = line.split('chromosome=')[1].split(']')[0]
            else:
                chromo = 'nope'
                
            chromo_dict[chromo] = ''
            
        else:
            chromo_dict[chromo] += line
            
    infast_file.close()
    
    length = 0
    
    for chromo in chromo_dict:
        if chromo != 'nope':
            length += len(chromo_dict[chromo])
            print(length)
    
    hit_counter = {}
            
    for nt2 in n_list:
        for nt3 in y_list:
            for nt6 in n_list:
                for nt7 in y_list:
                    pos_search = ('C{nt2}{nt3}TC{nt6}{nt7}T').format(
                        nt2=nt2, nt3=nt3, nt6=nt6, nt7=nt7)
                    neg_search = rev_comp(pos_search)
                    
                    if pos_search not in hit_counter:
                        hit_counter[pos_search] = 0
                        
                    for chromo in chromo_dict:
                        if chromo != 'nope':
                            count = 0
                            count += chromo_dict[chromo].count(pos_search)
                            count += chromo_dict[chromo].count(neg_search)
                        
                            hit_counter[pos_search] += count
    
    total_hit = 0
    
    for seq in hit_counter:
        total_hit += hit_counter[seq]
        
    freq_per_nt = total_hit/(length)
    
    print(freq_per_nt)
    
    


                
        
        