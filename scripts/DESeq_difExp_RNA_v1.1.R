#DESeq package: analysis of differentially expressed genes
library(DESeq2)

#install.packages("Rcpp", repos="http://cran.rstudio.com/", dependencies=TRUE)

#1726
#Make dataframe of the raw gene reads 
setwd('data/')
file <- 'combined_coverage_table_expected.tsv'
countData = read.table(file = file, header = TRUE, row.names = 1, sep = ',')
countData[is.na(countData)] <- 0 # Make NA values == 0
head(countData) #check the file organization
dim(countData) #check dimensions (rows vs. columns)
typeof(countData)
comp <- countData[c(1, 2, 3, 4)]
head(comp)

#Make the metadata dataframe:
sample_names <- c('DGY1726_1', 'DGY1726_2','DGY1726_exp_1', 'DGY1726_exp_2')
condition <- c('Evolved', 'Evolved', 'Control', 'Control')
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
write.table(result, "DESeq_Exp_RNA_DGY1657_DGY1726_v2.txt", sep = '\t')

#1735
#Make dataframe of the raw gene reads 
file <- 'combined_coverage_table_expected.tsv'
countData = read.table(file = file, header = TRUE, row.names = 1, sep = ',')
countData[is.na(countData)] <- 0 # Make NA values == 0
head(countData) #check the file organization
dim(countData) #check dimensions (rows vs. columns)
typeof(countData)
comp <- countData[c(5, 6, 7, 8)]
head(comp)

#Make the metadata dataframe:
sample_names <- c('DGY1735_1', 'DGY1735_2','DGY1735_exp_1', 'DGY1735_exp_2')
condition <- c('Evolved', 'Evolved', 'Control', 'Control')
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
write.table(result, "DESeq_Exp_RNA_DGY1657_DGY1735_v2.txt", sep = '\t')


#1741
#Make dataframe of the raw gene reads 
file <- 'combined_coverage_table_expected.tsv'
countData = read.table(file = file, header = TRUE, row.names = 1, sep = ',')
countData[is.na(countData)] <- 0 # Make NA values == 0
head(countData) #check the file organization
dim(countData) #check dimensions (rows vs. columns)
typeof(countData)
comp <- countData[c(9, 10, 11, 12)]
head(comp)

#Make the metadata dataframe:
sample_names <- c('DGY1741_1', 'DGY1741_2','DGY1741_exp_1', 'DGY1741_exp_2')
condition <- c('Evolved', 'Evolved', 'Control', 'Control')
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
write.table(result, "DESeq_Exp_RNA_DGY1657_DGY1741_v2.txt", sep = '\t')



#1743
#Make dataframe of the raw gene reads 
file <- 'combined_coverage_table_expected.tsv'
countData = read.table(file = file, header = TRUE, row.names = 1, sep = ',')
countData[is.na(countData)] <- 0 # Make NA values == 0
head(countData) #check the file organization
dim(countData) #check dimensions (rows vs. columns)
typeof(countData)
comp <- countData[c(13, 14, 15, 16)]
head(comp)

#Make the metadata dataframe:
sample_names <- c('DGY1743_1', 'DGY1743_2','DGY1743_exp_1', 'DGY1743_exp_2')
condition <- c('Evolved', 'Evolved', 'Control', 'Control')
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
write.table(result, "DESeq_Exp_RNA_DGY1657_DGY1743_v2.txt", sep = '\t')
