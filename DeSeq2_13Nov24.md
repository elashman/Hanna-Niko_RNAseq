# Overview of samples

![Screenshot 2024-11-13 at 11 42 14](https://github.com/user-attachments/assets/e8e1b831-93c6-47cf-84b2-dd5384735c2a)

# load libraries in RStudio

```r
library(dplyr) 
library(tximport) 
library(DESeq2) 
library(rhdf5) 
library(ggplot2) 
library(plotly) 
library(tibble) 
library(biomaRt)

gene_map <- read.csv ("map_genes_2.csv")
```

# Copy files

scp -r elashman@bea81.hpc.ista.ac.at:/nfs/scistore13/bonogrp/hschoen/RNAseq/kallistoKatya /Users/elashman/Documents/Niko_Hanna_RNASeq_13Nov24

```r
sample_table = read.csv("Info_samples_DEseq2.csv", sep=";") 
sample_files = paste0("./kallistoKatya/", pull(sample_table, "Sample") , "/abundance.h5") 
names(sample_files) = pull(sample_table, "Sample")
names(sample_files) 

count_data <- tximport(files = sample_files,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE) 
count_data[['counts']]
```

#  specify what samples are control
```r
conditions = sample_table[,3] 
conditions = factor(conditions, levels = c('WT_1h', 'ST7:AID_1h', 'MACO-1:AID_1h', 'WT_4h', 'ST7:AID_4h', 'MACO-1:AID_4h' )) 
sample_table$conditions = conditions

deseq_dataset = DESeqDataSetFromTximport(txi = count_data, colData =sample_table, design = ~conditions)
```

# 3 STEPS in DESeq2 analysis
# STEP 1: estimate size factors (normalizatiion)
```r
deseq_dataset = estimateSizeFactors(deseq_dataset) 
normalizationFactors(deseq_dataset) 
counts(deseq_dataset, normalized=TRUE)

vst = varianceStabilizingTransformation(deseq_dataset) 
plotPCA(vst, intgroup='conditions')
```
# PCA plot with samples names
```r
plotPCA(vst, intgroup='conditions') + geom_label(aes(label=name))
```
# FOR WT_1h versus ST7_1h

# Exclude all the samples that are not needed
```r
sample_files_ST7_1h = sample_files[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
count_data_ST7_1h <- tximport(files = sample_files_ST7_1h,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE) 
count_data_ST7_1h[['counts']]
```
# specify what samples are control
```r
sample_table_ST7_1h = read.csv("Info_samples_DEseq2_ST7_1h.csv", sep=";")
conditions = sample_table_ST7_1h[,3] 
conditions = factor(conditions, levels = c('WT_1h', 'ST7:AID_1h'))
sample_table_ST7_1h$conditions = conditions

deseq_dataset_ST7_1h = DESeqDataSetFromTximport(txi = count_data_ST7_1h, colData = sample_table_ST7_1h, design = ~conditions)
```
# estimate size factors (normalizatiion)
```r
deseq_dataset_ST7_1h = estimateSizeFactors(deseq_dataset_ST7_1h) 
normalizationFactors(deseq_dataset_ST7_1h) 
counts(deseq_dataset_ST7_1h, normalized=TRUE)
```
# estimate dispersions
```r
deseq_dataset_ST7_1h = estimateDispersions(deseq_dataset_ST7_1h) 
plotDispEsts(deseq_dataset_ST7_1h)
```
# Run statistics
```r
deseq_dataset_ST7_1h = nbinomWaldTest(deseq_dataset_ST7_1h) 
result_table_ST7_1h = results(deseq_dataset_ST7_1h) 
summary(result_table_ST7_1h)
```
Result:
out of 18854 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 608, 3.2%
LFC < 0 (down)     : 428, 2.3%
outliers [1]       : 7, 0.037%
low counts [2]     : 4704, 25%
(mean count < 9)

# MA plot
```r
plotMA(result_table_ST7_1h)
```
# Make results to data.frame
```r
result_df_ST7_1h = as.data.frame(result_table_ST7_1h)

sum(complete.cases(result_df_ST7_1h))

filter_df_ST7_1h = result_df_ST7_1h[complete.cases(result_df_ST7_1h),] 
View(filter_df_ST7_1h)
```
# in filter_df_ST7_1h 14143 genes

# padj < 0.05
```r
filter_df_ST7_1h$padj < 0.05
filter_df2_ST7_1h = filter_df_ST7_1h[filter_df_ST7_1h$padj < 0.05, ]
```
# in filter_df2_ST7_1h 579 genes

# log2FoldChange >1 <-1
```r
abs(filter_df2_ST7_1h$log2FoldChange) > 1 
filter_df3_ST7_1h = filter_df2_ST7_1h[abs(filter_df2_ST7_1h$log2FoldChange) > 1, ] 
View(filter_df3_ST7_1h)
```
# in filter_df3_ST7_1h 120 genes


# Assign gene names
```r
library("tidyverse")

filter_df_ST7_1h = rownames_to_column(filter_df_ST7_1h, var = "ensgene")

listMarts() 
ensembl_110 = useEnsembl(biomart = "ensembl", version = 110) 
ensembl_110 = useDataset("celegans_gene_ensembl", mart = ensembl_110) 

annotation_df_ST7_1h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df_ST7_1h$ensgene, mart = ensembl_110) 
View(annotation_df_ST7_1h) 
annotated_df_ST7_1h = left_join(filter_df_ST7_1h, annotation_df_ST7_1h, by = c("ensgene" = "ensembl_gene_id"))

filter_df3_ST7_1h = rownames_to_column(filter_df3_ST7_1h, var = "ensgene") 
annotation_filter_df3_ST7_1h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df3_ST7_1h$ensgene, mart = ensembl_110) 
annotated_df3_ST7_1h = left_join(filter_df3_ST7_1h, annotation_filter_df3_ST7_1h, by = c("ensgene" = "ensembl_gene_id"))

write_tsv(annotated_df3_ST7_1h, "annotated_df3_ST7_1h")
```





