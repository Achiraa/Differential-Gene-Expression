setwd("C:/geo_data/Stromal")

#Installing libraries
library(DESeq2)
library(tidyverse)
library(BiocGenerics)
library(S4Vectors)
library(stats4)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cummeRbund)
library(RColorBrewer)

#reading the sample info of the experiment containing samplenames, cell info
sdata=read.csv("Stromaltomatosource.csv")
dim(sdata)
sdata
head(sdata)
View(sdata)

#reading the csv file with gene counts and ensembleID
rw=read.csv("Stromaltomato.csv")
View(rw)
rw1=rw[,-1]
rw1
row.names(rw1)=make.names(rw[,1],unique = T)
head(rw1)
dim(rw1)
View(rw1)

#checking whether the rownames and column names in both the files are similar
all(colnames(rw1) %in% rownames(sdata))
all(colnames(rw1) == rownames(sdata))

#converting the count data into different forms like class, assays etc
#used to store the input values, intermediate calculations and results of an analysis of differential expression
dds_gene=DESeqDataSetFromMatrix(countData = rw1,colData = sdata, design = ~ Treatment)
dds_gene

#only keeping the genes which have at least 10 counts
keepgene= rowSums(counts(dds_gene)) >= 10
keepgene
dds_gene1=dds_gene[keepgene,]
dds_gene1

#Generate the normalized counts, it uses he median of ratios method of normalization
dds_genenew=estimateSizeFactors(dds_gene1)
dds_genenew

#checking the normalized counts
sizeFactors(dds_genenew)

#To  improve  the  distances/clustering  for  the  PCA  and  hierarchical  clustering  visualization  methods,  
#we  need  to  moderate the variance across the mean by applying the rlog transformation to the normalized counts.
rld=rlog(dds_genenew, blind=TRUE)
rld

#Differential expression analysis with DESeq2
newdds_gene=DESeq(dds_genenew)
newdds_gene

newdds_gene$Treatment=relevel(dds_genenew$Treatment,ref = "Untreated")
newdds_gene$Treatment

#gives us the different column i.e. log2foldchange, baseMean, Pval, padj
resultss=results(newdds_gene)
resultss

#gives us the info regarding how many genes are up and down regulated
summary(resultss)

#taking FDR>0.01 or psdj>0.01
result.05=results(newdds_gene, alpha=0.05)
result.05
summary(result.05)

#getting the most differentiated expressed genes and saving in CSV
resultsort=result.05[order(result.05$pvalue),]
write.csv(resultsort,file = "ALLstromaltomatoGenes.csv")

# Filter out rows with NA p-values, and then subset by p-value
resultsort = result.05[order(result.05$pvalue),]
significant_genes = resultsort[complete.cases(resultsort) & resultsort$pvalue < 0.05, ]
significant_genes
dim(significant_genes)

write.csv(significant_genes,file="Significant tomato genes.csv")
