setwd("C:/geo_data/Luminal")

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
source=read.csv("sourcedata.csv")
dim(source)
source
head(source)
View(source)

#reading the csv file with gene counts and ensembleID
raw=read.csv("luminalrawdata.csv")
View(raw)
raw1=raw[,-1]
raw1
row.names(raw1)=make.names(raw[,1],unique = T)
head(raw1)
dim(raw1)
View(raw1)

#checking whether the rownames and column names in both the files are similar
all(colnames(raw1) %in% rownames(source))
all(colnames(raw1) == rownames(source))

#converting the count data into different forms like class, assays etc
#used to store the input values, intermediate calculations and results of an analysis of differential expression
dds=DESeqDataSetFromMatrix(countData = raw1,colData = source, design = ~ Treatment)
dds

#only keeping the genes which have at least 10 counts
keep= rowSums(counts(dds)) >= 10
keep
dds1=dds[keep,]
dds1

#Generate the normalized counts, it uses he median of ratios method of normalization
ddsnew=estimateSizeFactors(dds1)
ddsnew

#checking the normalized counts
sizeFactors(ddsnew)

#saving the normalized counts
#retrieve  the  normalized  counts  matrix  from  dds,  we  use  the  counts()  function  and  
#add  the  argument  normalized=TRUE.
nm1=counts(ddsnew,normalized=TRUE)
nm1
wt1=write.table(nm1, file=" normalized_countsluminal.txt", sep="\t", quote=F,col.names=NA)
wt1

#To  improve  the  distances/clustering  for  the  PCA  and  hierarchical  clustering  visualization  methods,  
#we  need  to  moderate the variance across the mean by applying the rlog transformation to the normalized counts.
rld=rlog(ddsnew, blind=TRUE)
rld

#Differential expression analysis with DESeq2
newdds=DESeq(ddsnew)
newdds

newdds$Treatment=relevel(ddsnew$Treatment,ref = "Untreated")
newdds$Treatment

#gives us the different column i.e. log2foldchange, baseMean, Pval, padj
resultss=results(newdds)
resultss

#gives us the info regarding how many genes are up and down regulated
summary(resultss)
xy=write.csv(resultss,file = "resuts1.csv",row.names = TRUE,append = TRUE)
xy

#taking FDR>0.01 or psdj>0.01
result.05=results(newdds, alpha=0.05)
result.05
write.csv(result.05,file = "resuts0.5.csv",row.names = TRUE,append = TRUE)
summary(result.05)

#getting the most differentiated expressed genes and saving in CSV
resultsort=result.05[order(result.05$pvalue),]
write.csv(resultsort,file = "All Genes.csv")

# Filter out rows with NA p-values, and then subset by p-value
resultsort = result.05[order(result.05$pvalue),]
significant_genes = resultsort[complete.cases(resultsort) & resultsort$pvalue < 0.05, ]
significant_genes
dim(significant_genes)

# Save the significant genes to a CSV file
write.csv(significant_genes, file = "significant_genes.csv")

#PCA plot
plotPCA(rld, intgroup="CellType")
plotPCA(rld, intgroup="Treatment")

#MA Plot
plotMA(resultss)
plotMA(result.05)
#idx=identify(resultss$baseMean, resultss$log2FoldChange)
#rownames(resultss)[idx]

#volcano
v=EnhancedVolcano(result.05,lab=rownames(newdds),x="log2FoldChange",y="padj")
v

#heatmap
topgenes=head(rownames(topgene),3714)
topgenes
mat1=assay(newdds)[topgenes,]
mat1
mat2=mat1 -rowMeans(mat1)
View(mat2)

#with annotation
annodf=as.data.frame(source)
annodf

annot=HeatmapAnnotation(Cell.Type=annodf$CellType, col = list(CellType =c("eGFP+"="darkgreen","tdtomato+"="darkred")))
annot
annot2=HeatmapAnnotation(Treatment=annodf$Treatment,col=list(Treatment=c("Untreated"="yellow","Treated"="purple")))

heatmap(mat2, name="Gene Expression for Luminal", bottom_annotation = annot, top_annotation = annot2, show_column_names = F, cluster_rows = T)
Heatmap(mat2, name = "Gene expression for Luminal", bottom_annotation = annot,top_annotation = annot2, show_row_names = F, show_column_names = F, cluster_rows = T)

#pheatmap
newq=assay(rld)
select= order(rowMeans(counts(newdds, normalized=T)),decreasing = T)
pp=newq[select,]
pp

a1=data.frame(colData(dds)[,c("CellType","Treatment")])
a1
heat_colors <- brewer.pal(11, "RdYlBu")
heat_colors

pheatmap(mat2,
         main = "Luminal Heatmap", 
         color = heat_colors,
         cluster_rows = T, 
         show_rownames = F,
         show_colnames = T,
         cluster_cols = F,
         annotation_col = a1, 
         border_color = NA, 
         fontsize = 12, 
         fontsize_col = 10,
         scale = "row", 
         fontsize_row = 4)





