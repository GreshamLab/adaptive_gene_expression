library(pheatmap) 
library(RColorBrewer)

unit_object_v3 <- read.delim("data/output/unit_object_level_Unit_l1_tests")

unit_object_v3 <- subset(unit_object_v3, (( unit_object_v3$deseq_obs_RNA_DGY1726padj < 0.01 |
                                            unit_object_v3$deseq_obs_RNA_DGY1735padj < 0.01 |
                                            unit_object_v3$deseq_obs_RNA_DGY1741padj < 0.01 |
                                            unit_object_v3$deseq_obs_RNA_DGY1743padj < 0.01) & (
                                            unit_object_v3$deseq_obs_RPF_DGY1726padj < 0.01 |
                                            unit_object_v3$deseq_obs_RPF_DGY1735padj < 0.01 |
                                            unit_object_v3$deseq_obs_RPF_DGY1741padj < 0.01 |
                                            unit_object_v3$deseq_obs_RPF_DGY1743padj < 0.01) & (
                                            unit_object_v3$N..Welch.s.T.test.q.value.DGY1726_DGY1657 < 0.01 |
                                            unit_object_v3$N..Welch.s.T.test.q.value.DGY1735_DGY1657 < 0.01 |
                                            unit_object_v3$N..Welch.s.T.test.q.value.DGY1741_DGY1657 < 0.01 |
                                            unit_object_v3$N..Welch.s.T.test.q.value.DGY1743_DGY1657 < 0.01) | (
                                            unit_object_v3$DGY1726 != 1 |
                                            unit_object_v3$DGY1735 != 1 |
                                            unit_object_v3$DGY1741 != 1 |
                                            unit_object_v3$DGY1743 != 1 )
                                          ))


unit_object_v3$DGY1726_rna_diff<-(unit_object_v3[c("DGY1726_rna_pt_median")])-(unit_object_v3[c("DGY1657_rna_pt_median")])
unit_object_v3$DGY1735_rna_diff<-(unit_object_v3[c("DGY1735_rna_pt_median")])-(unit_object_v3[c("DGY1657_rna_pt_median")])
unit_object_v3$DGY1741_rna_diff<-(unit_object_v3[c("DGY1741_rna_pt_median")])-(unit_object_v3[c("DGY1657_rna_pt_median")])
unit_object_v3$DGY1743_rna_diff<-(unit_object_v3[c("DGY1743_rna_pt_median")])-(unit_object_v3[c("DGY1657_rna_pt_median")])

unit_object_v3$DGY1726_rpf_diff<-(unit_object_v3[c("DGY1726_rpf_pt_median")])-(unit_object_v3[c("DGY1657_rpf_pt_median")])
unit_object_v3$DGY1735_rpf_diff<-(unit_object_v3[c("DGY1735_rpf_pt_median")])-(unit_object_v3[c("DGY1657_rpf_pt_median")])
unit_object_v3$DGY1741_rpf_diff<-(unit_object_v3[c("DGY1741_rpf_pt_median")])-(unit_object_v3[c("DGY1657_rpf_pt_median")])
unit_object_v3$DGY1743_rpf_diff<-(unit_object_v3[c("DGY1743_rpf_pt_median")])-(unit_object_v3[c("DGY1657_rpf_pt_median")])

unit_object_v3$DGY1726_ms_diff<-(unit_object_v3[c("DGY1726_ms_pt_median")])-(unit_object_v3[c("DGY1657_ms_pt_median")])
unit_object_v3$DGY1735_ms_diff<-(unit_object_v3[c("DGY1735_ms_pt_median")])-(unit_object_v3[c("DGY1657_ms_pt_median")])
unit_object_v3$DGY1741_ms_diff<-(unit_object_v3[c("DGY1741_ms_pt_median")])-(unit_object_v3[c("DGY1657_ms_pt_median")])
unit_object_v3$DGY1743_ms_diff<-(unit_object_v3[c("DGY1743_ms_pt_median")])-(unit_object_v3[c("DGY1657_ms_pt_median")])

r_data <- unit_object_v3[c("X",
                           "DGY1726_rna_diff", "DGY1726_rpf_diff", "DGY1726_ms_diff",
                           "DGY1735_rna_diff", "DGY1735_rpf_diff", "DGY1735_ms_diff",
                           "DGY1741_rna_diff", "DGY1741_rpf_diff", "DGY1741_ms_diff",
                           "DGY1743_rna_diff", "DGY1743_rpf_diff", "DGY1743_ms_diff")]


#set specific column as row names
rownames(r_data) <- r_data$X

#remove original column from data frame
r_data$X <- NULL

paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(r_data), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(r_data)/paletteLength, max(r_data), length.out=floor(paletteLength/2)))

myBreaks = myBreaks + seq_along(myBreaks) * .Machine$double.eps

clustered_data <- pheatmap(r_data,
                           main = "unscaled heatmap",
                           breaks = myBreaks,
                           color = myColor
)
pdf('C:/Gresham/Project_Carolino/figures/l1_analysis/_pheatmap_all_24_score.pdf')

clustered_data <- pheatmap(r_data,
                           scale = 'none',
                           cluster_cols = FALSE,
                           main = "all", 
                           breaks = myBreaks,
                           color = myColor,
                           cutree_rows = 24)
dev.off()


clust_24 <- cbind(r_data, 
               cluster = cutree(clustered_data$tree_row, 24))

which_clust_24 <- clust[c("cluster")]

write.table(which_clust,"data/output/unit_object_all_24k_which.clust.tab", sep = '\t', row.names = TRUE)

