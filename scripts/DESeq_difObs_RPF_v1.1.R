#DESeq package: analysis of differentially expressed genes
library(DESeq2)
library(RColorBrewer)
library(yaml)
library(gridExtra)
library(ggplot2)
#install.packages("Rcpp", repos="http://cran.rstudio.com/", dependencies=TRUE)

#1726
#Make dataframe of the raw gene reads 
setwd('data/')
file <- 'RPF_unnormed_read_counts_by_gene.tsv'
countData = read.table(file = file, header = TRUE, row.names = 1, sep = '\t')
countData[is.na(countData)] <- 0 # Make NA values == 0
head(countData) #check the file organization
dim(countData) #check dimensions (rows vs. columns)
typeof(countData)
comp <- countData[c(1, 2, 3, 4)]
head(comp)

#Make the metadata dataframe:
sample_names <- c('DGY1657_1', 'DGY1657_2', 'DGY1726_1', 'DGY1726_2')
condition <- c('Control', 'Control', 'Evolved', 'Evolved')
colData <- data.frame(row.names=sample_names, condition=factor(condition, levels=c('Control','Evolved')))
colData
#Create DESeqDataSet object:
dataset <- DESeqDataSetFromMatrix(countData = comp,
                                  colData = colData,
                                  design = ~condition)
dataset

# Now run the DESeq2 algorithm and extract results for the two-class comparison
dds <- DESeq(dataset)
result <- results(dds, contrast=c('condition','Evolved','Control'))
result <- result[complete.cases(result),]  #remove any rows with NA
head(result)
write.table(result, "DESeq_Obs_RPF_DGY1657_DGY1726.txt", sep = '\t')

summary(result) # summary of DESeq2 results
# Top 50 up-regulated and down-regulated genes by p-value:
n = 50
resOrdered <- result[order(result$padj),]
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                     resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]
#Plot log fold change vs. mean expression. Genes with p < 0.1 colored red:
plotMA(result, main='DGY1657 vs DGY1315', ylim=c(-2,2))
# Plot of dispersion:
plotDispEsts(dds, ylim = c(1e-6, 1e1))
# Histogram of p-values:
hist(result$pvalue, breaks=10, col="grey")
# Plot Control vs. Evolved PCA: 
#Data is logtransformed
rld <- rlogTransformation(dds, blind=TRUE)
plotPCA(rld) + geom_text(aes(label=name),vjust=0.5) #geom adds the labels
# Plot counts for a single gene (with the lowest p-value):
plotCounts(dds, gene=which.min(result$padj), intgroup='condition', pch = 19)
# Top genes normalized counts in heatmap:
hmcol <- brewer.pal(11,'RdBu') #Determine color palette (max: 11)
nCounts <- counts(dds, normalized=TRUE)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(8,2))


###1735
#Make dataframe of the raw gene reads 
file <- 'RPF_unnormed_read_counts_by_gene.tsv'
countData = read.table(file = file, header = TRUE, row.names = 1, sep = '\t')
countData[is.na(countData)] <- 0 # Make NA values == 0
head(countData) #check the file organization
dim(countData) #check dimensions (rows vs. columns)
typeof(countData)
comp <- countData[c(1, 2,7,8)]
head(comp)

#Make the metadata dataframe:
sample_names <- c('DGY1657_1', 'DGY1657_2', 'DGY1735_1', 'DGY1735_2')
condition <- c('Control', 'Control', 'Evolved', 'Evolved')
colData <- data.frame(row.names=sample_names, condition=factor(condition, levels=c('Control','Evolved')))
colData
#Create DESeqDataSet object:
dataset <- DESeqDataSetFromMatrix(countData = comp,
                                  colData = colData,
                                  design = ~condition)
dataset

# Now run the DESeq2 algorithm and extract results for the two-class comparison
dds <- DESeq(dataset)
result <- results(dds, contrast=c('condition','Evolved','Control'))
result <- result[complete.cases(result),]  #remove any rows with NA
head(result)
write.table(result, "DESeq_Obs_RPF_DGY1657_DGY1735.txt", sep = '\t')

summary(result) # summary of DESeq2 results
# Top 50 up-regulated and down-regulated genes by p-value:
n = 50
resOrdered <- result[order(result$padj),]
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                     resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]
#Plot log fold change vs. mean expression. Genes with p < 0.1 colored red:
plotMA(result, main='DGY1657 vs DGY1315', ylim=c(-2,2))
# Plot of dispersion:
plotDispEsts(dds, ylim = c(1e-6, 1e1))
# Histogram of p-values:
hist(result$pvalue, breaks=10, col="grey")
# Plot Control vs. Evolved PCA: 
#Data is logtransformed
rld <- rlogTransformation(dds, blind=TRUE)
plotPCA(rld) + geom_text(aes(label=name),vjust=0.5) #geom adds the labels
# Plot counts for a single gene (with the lowest p-value):
plotCounts(dds, gene=which.min(result$padj), intgroup='condition', pch = 19)
# Top genes normalized counts in heatmap:
hmcol <- brewer.pal(11,'RdBu') #Determine color palette (max: 11)
nCounts <- counts(dds, normalized=TRUE)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(8,2))

###1741
#Make dataframe of the raw gene reads, filtering replicate two
file <- 'RPF_unnormed_read_counts_by_gene.tsv'
countData = read.table(file = file, header = TRUE, row.names = 1, sep = '\t')
countData[is.na(countData)] <- 0 # Make NA values == 0
head(countData) #check the file organization
dim(countData) #check dimensions (rows vs. columns)
typeof(countData)
comp <- countData[c(1, 2, 9, 10)]
head(comp)

#Make the metadata dataframe:
# filtering replicate 2 (ie "quad_2" from analysis du to outlier expression)
sample_names <- c('DGY1657_1', 'DGY1657_2', 'DGY1741_1', 'DGY1741_2')
condition <- c('Control', 'Control', 'Evolved', 'Evolved')
colData <- data.frame(row.names=sample_names, condition=factor(condition, levels=c('Control','Evolved')))
colData
#Create DESeqDataSet object:
dataset <- DESeqDataSetFromMatrix(countData = comp,
                                  colData = colData,
                                  design = ~condition)
dataset

# Now run the DESeq2 algorithm and extract results for the two-class comparison
dds <- DESeq(dataset)
result <- results(dds, contrast=c('condition','Evolved','Control'))
result <- result[complete.cases(result),]  #remove any rows with NA
head(result)
write.table(result, "DESeq_Obs_RPF_DGY1657_DGY1741.txt", sep = '\t')

summary(result) # summary of DESeq2 results
# Top 50 up-regulated and down-regulated genes by p-value:
n = 50
resOrdered <- result[order(result$padj),]
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                     resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]
#Plot log fold change vs. mean expression. Genes with p < 0.1 colored red:
plotMA(result, main='DGY1657 vs DGY1315', ylim=c(-2,2))
# Plot of dispersion:
plotDispEsts(dds, ylim = c(1e-6, 1e1))
# Histogram of p-values:
hist(result$pvalue, breaks=10, col="grey")
# Plot Control vs. Evolved PCA: 
#Data is logtransformed
rld <- rlogTransformation(dds, blind=TRUE)
plotPCA(rld) + geom_text(aes(label=name),vjust=0.5) #geom adds the labels
# Plot counts for a single gene (with the lowest p-value):
plotCounts(dds, gene=which.min(result$padj), intgroup='condition', pch = 19)
# Top genes normalized counts in heatmap:
hmcol <- brewer.pal(11,'RdBu') #Determine color palette (max: 11)
nCounts <- counts(dds, normalized=TRUE)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(8,2))

###1743
#Make dataframe of the raw gene reads 
file <- 'RPF_unnormed_read_counts_by_gene.tsv'
countData = read.table(file = file, header = TRUE, row.names = 1, sep = '\t')
countData[is.na(countData)] <- 0 # Make NA values == 0
head(countData) #check the file organization
dim(countData) #check dimensions (rows vs. columns)
typeof(countData)
comp <- countData[c(1, 2, 11, 12)]
head(comp)

#Make the metadata dataframe:
sample_names <- c('DGY1657_1', 'DGY1657_2', 'DGY1743_1', 'DGY1743_2')
condition <- c('Control', 'Control', 'Evolved', 'Evolved')
colData <- data.frame(row.names=sample_names, condition=factor(condition, levels=c('Control','Evolved')))
colData
#Create DESeqDataSet object:
dataset <- DESeqDataSetFromMatrix(countData = comp,
                                  colData = colData,
                                  design = ~condition)
dataset

# Now run the DESeq2 algorithm and extract results for the two-class comparison
dds <- DESeq(dataset)
result <- results(dds, contrast=c('condition','Evolved','Control'))
result <- result[complete.cases(result),]  #remove any rows with NA
head(result)
write.table(result, "DESeq_Obs_RPF_DGY1657_DGY1743.txt", sep = '\t')

summary(result) # summary of DESeq2 results
# Top 50 up-regulated and down-regulated genes by p-value:
n = 50
resOrdered <- result[order(result$padj),]
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                     resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]
#Plot log fold change vs. mean expression. Genes with p < 0.1 colored red:
plotMA(result, main='DGY1657 vs DGY1315', ylim=c(-2,2))
# Plot of dispersion:
plotDispEsts(dds, ylim = c(1e-6, 1e1))
# Histogram of p-values:
hist(result$pvalue, breaks=10, col="grey")
# Plot Control vs. Evolved PCA: 
#Data is logtransformed
rld <- rlogTransformation(dds, blind=TRUE)
plotPCA(rld) + geom_text(aes(label=name),vjust=0.5) #geom adds the labels
# Plot counts for a single gene (with the lowest p-value):
plotCounts(dds, gene=which.min(result$padj), intgroup='condition', pch = 19)
# Top genes normalized counts in heatmap:
hmcol <- brewer.pal(11,'RdBu') #Determine color palette (max: 11)
nCounts <- counts(dds, normalized=TRUE)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(8,2))

