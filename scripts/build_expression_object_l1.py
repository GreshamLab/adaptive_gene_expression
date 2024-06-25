# -*- coding: utf-8 -*-
"""
Created on Wed May 29 11:29:30 2024

@author: pspea
"""

import pandas as pd

import plotly.express as px
#from umap import UMAP
import plotly.io as pio
pio.renderers.default = "browser"

filename = ('data/analyses/efficiency/unit_object_level_2v13v16_Unit_10pct.tab')
ratio_df = pd.read_table(filename, index_col=0)

for istype in ['RNA', 'RPF']:
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        deseq_filename = ('data/analyses/deseq/'
                          'DESeq_Obs_{istype}_DGY1657_{strain}.txt').format(strain = strain, istype = istype)
        deseq_df = pd.read_table(deseq_filename, index_col=0)
        prefix = ('deseq_obs_{istype}_{strain}').format(istype = istype, strain = strain)
        deseq_df = deseq_df.add_prefix(prefix)
        
        ratio_df = pd.merge(left = ratio_df,
                                   right = deseq_df,
                                   left_index=True,
                                   right_index=True)
        
#
wttest_filename = ('data/analyses/ms/Perseus_DA_Welchs_t-test_wFDR_edit.txt')
wttest_df = pd.read_table(wttest_filename, index_col='T: Majority.protein.IDs')
#wttest_df = 
       
ratio_df = pd.merge(left = ratio_df,
                           right = wttest_df,
                           left_index=True,
                           right_index=True)

outfile_name = ('data/output/unit_object_level_Unit_l1_tests.tab')
ratio_df.to_csv(outfile_name, sep = '\t') 

#unit_object_v3$DGY1726_rna_diff<-(unit_object_v3[c("DGY1726_rna_pt_median")])-(unit_object_v3[c("DGY1657_rna_pt_median")])
diff_colname_list = []
for istype in ['rna', 'rpf', 'ms']:
    for strain in set(['DGY1726', 'DGY1735', 'DGY1741', 'DGY1743']):
        new_colname = ('{strain}_{istype}_diff').format(strain = strain, istype = istype)
        evo_colname = ('{strain}_{istype}_pt_median').format(strain = strain, istype = istype)
        anc_colname = ('DGY1657_{istype}_pt_median').format(istype = istype)
        
        ratio_df[new_colname] = ratio_df[evo_colname]-ratio_df[anc_colname]
        diff_colname_list.append(evo_colname)
        diff_colname_list.append(anc_colname)
        
ratio_df["cn_count"] = ratio_df['DGY1726'] + ratio_df['DGY1735'] + ratio_df['DGY1741'] + ratio_df['DGY1743']
    
diff_colname_list = ['DGY1726_rna_diff', 'DGY1726_rpf_diff', 'DGY1726_ms_diff', 
                      'DGY1735_rna_diff', 'DGY1735_rpf_diff', 'DGY1735_ms_diff',
                        'DGY1741_rna_diff', 'DGY1741_rpf_diff', 'DGY1741_ms_diff',
                        'DGY1743_rna_diff', 'DGY1743_rpf_diff', 'DGY1743_ms_diff']

data = ratio_df[diff_colname_list]

# umap_2d = UMAP(random_state=0)
# umap_2d.fit(data)

# projections = umap_2d.transform(data)

# fig = px.scatter(
#     projections, x=0, y=1,
#     color=ratio_df.cn_count.astype(str), labels={'color': 'CN count'}
# )
# fig.show()

# import plotly.express as px


# def KMeans_cluster(data, k):    # 2. initialize the model    
#     my_kmeans = KMeans(n_clusters= k)    # 3. fit the model to the data    
#     my_kmeans.fit(data) # pass your scaled data here    # 4. obtain the cluster output    
#     clusters = my_kmeans.predict(data) # pass your scaled data here    
#     centroids = my_kmeans.cluster_centers_    
#     return clusters,  pd.DataFrame(centroids)

# clusters, centroids = KMeans_cluster(transformed_data, 9)

clust = "unit_object_all_24k_which.clust.tab"
clust_df = pd.read_table(clust, index_col=0)

ratio_df = pd.merge(left = ratio_df,
                           right = clust_df,
                           left_index=True,
                           right_index=True)

limited_df = ratio_df.loc[
    (
    (
     (ratio_df['deseq_obs_RNA_DGY1726padj'] < 0.01) |
     (ratio_df['deseq_obs_RNA_DGY1735padj'] < 0.01) |
     (ratio_df['deseq_obs_RNA_DGY1741padj'] < 0.01) |
     (ratio_df['deseq_obs_RNA_DGY1743padj'] < 0.01) 
     ) & (
         (ratio_df['deseq_obs_RPF_DGY1726padj'] < 0.01) |
         (ratio_df['deseq_obs_RPF_DGY1735padj'] < 0.01) |
         (ratio_df['deseq_obs_RPF_DGY1741padj'] < 0.01) |
         (ratio_df['deseq_obs_RPF_DGY1743padj'] < 0.01) 
         ) & (
             (ratio_df["N: Welch's T-test q-value DGY1726_DGY1657"] < 0.01) |
             (ratio_df["N: Welch's T-test q-value DGY1735_DGY1657"] < 0.01) |
             (ratio_df["N: Welch's T-test q-value DGY1741_DGY1657"] < 0.01) |
             (ratio_df["N: Welch's T-test q-value DGY1743_DGY1657"] < 0.01) 
             ) | (
                 (ratio_df['DGY1726'] != 1) |
                 (ratio_df['DGY1735'] != 1) |
                 (ratio_df['DGY1741'] != 1) |
                 (ratio_df['DGY1743'] != 1) 
                 )
             )
    ]
 

limited_dict = limited_df.to_dict('index')

cluster_dict = dict()

for gene in limited_dict:
    cluster = limited_dict[gene]['cluster']
    if cluster not in cluster_dict:
        cluster_dict[cluster] = 0
    cluster_dict[cluster] += 1
    
select_cluster_set = set()
total_gene = 0
for cluster in cluster_dict:
    if cluster_dict[cluster] > 10:
        print(cluster)
        total_gene+=cluster_dict[cluster]
        select_cluster_set.add(cluster)
        
limited_data = [[-0.31, 0.31, -0.31, 0.31,-0.31, 0.31, -0.31, 0.31, -0.31, 0.31, -0.31, 0.31]]
gene_list = ['_ismin,ismax']

blank_line = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


for select_cluster in select_cluster_set:
    for i in range(3):
        limited_data.append(blank_line)
        gene_blank = ('_blank_{}').format(len(gene_list))
        gene_list.append(gene_blank)
    
    for gene in limited_dict:
        cluster = limited_dict[gene]['cluster']
        
        if cluster == select_cluster:
            gene_difference = []
            
            gene_cluster = ('{}_{}').format(gene, cluster)
            gene_list.append(gene_cluster)
            
            for diff_colname in diff_colname_list:
                gene_difference.append(limited_dict[gene][diff_colname])
                
            limited_data.append(gene_difference)
    
# ismin = 10
# ismax = 0
# for gene_exp in limited_data:
#     for exp in gene_exp:
#         if exp < ismin:
#             ismin = exp
#         if exp > ismax:
#             ismax = exp

#data=[[1, 25, 30, 50, 1], [20, 1, 60, 80, 30], [30, 60, 1, 5, 20]]
fig = px.imshow(limited_data,
                labels=dict(x="Strain and Level", y="Gene", color="Expression"),
                x=diff_colname_list,
                y=gene_list,
                color_continuous_scale='RdBu_r',
                width=2000, height=2000, aspect="auto"
               )
#fig.update_xaxes(side="top")
fig.show()
output_figure_name = "figures/unit_object_select_24k_clust.pdf"
fig.write_image(output_figure_name)