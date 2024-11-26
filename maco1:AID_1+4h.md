# FOR WT_1+4h versus ST7_1+4h

setwd("/Users/elashman/Documents/Niko_Hanna_RNASeq_13Nov24")

# Exclude all the samples that are not needed

```R
library(dplyr)
library(tximport)
library(DESeq2)
library(rhdf5)
library(ggplot2)
library(plotly)
library(tibble)
library(biomaRt)


gene_map <- read.csv ("map_genes_2.csv")

sample_table = read.csv("Info_samples_DEseq2_1+4h.csv", sep=";")
sample_table_maco1 = sample_table[c(1,2,3,4,5,6,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 31, 32, 33, 34, 35, 36), ] 
sample_files = paste0("./kallistoKatya/", pull(sample_table, "Sample") , ".1/abundance.h5")
names(sample_files) = pull(sample_table, "Sample")
names(sample_files)

sample_files_maco1 = sample_files[c(1,2,3,4,5,6,13,14,15,16,17,18,19, 20, 21, 22, 23, 24, 31, 32, 33, 34, 35, 36)]
count_data_maco1 <- tximport(files = sample_files_maco1,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE)
count_data_maco1[['counts']]
```

# specify what samples are control
```R
conditions = sample_table_maco1[,3]
conditions = factor(conditions, levels = c('WT', 'MACO-1:AID'))
sample_table_maco1$conditions = conditions

deseq_dataset_maco1 = DESeqDataSetFromTximport(txi = count_data_maco1, colData = sample_table_maco1, design = ~conditions)
```
# estimate size factors (normalizatiion)

```R
deseq_dataset_maco1 = estimateSizeFactors(deseq_dataset_maco1)
normalizationFactors(deseq_dataset_maco1)
counts(deseq_dataset_maco1, normalized=TRUE)

vst = varianceStabilizingTransformation(deseq_dataset_maco1)
plotPCA(vst, intgroup='conditions') + geom_label(aes(label=name))

d = assay(vst) 
d = t(d) 
d = dist(d) 
h = hclust(d) 
plot(h)
```
![PCAplot_maco1_1+4h](https://github.com/user-attachments/assets/c56403e4-0675-4156-b6d1-702a40a387af)

![Cluster_Dendrogram_maco1_1+4h](https://github.com/user-attachments/assets/6a021f64-7207-4dc2-bb47-04ab2c82dd51)


# estimate dispersions

```R
deseq_dataset_maco1 = estimateDispersions(deseq_dataset_maco1)
plotDispEsts(deseq_dataset_maco1)
```
![Disperson_maco1_1+4h](https://github.com/user-attachments/assets/37f6aacc-d0a9-4069-bef6-c711e4e88deb)

# Run statistics
```R
deseq_dataset_maco1 = nbinomWaldTest(deseq_dataset_maco1) 
result_table_maco1 = results(deseq_dataset_maco1) 
summary(result_table_maco1)
```
out of 19255 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1272, 6.6%
LFC < 0 (down)     : 1567, 8.1%
outliers [1]       : 134, 0.7%
low counts [2]     : 2574, 13%
(mean count < 1)


# MA plot
```R
plotMA(result_table_maco1)
```
![MAplot_maco1_1+4h](https://github.com/user-attachments/assets/56a25b1f-832f-44ed-8288-de351b00df3a)

# Make results to data.frame

```R
result_df_maco1 = as.data.frame(result_table_maco1)

sum(complete.cases(result_df_maco1))

filter_df_maco1 = result_df_maco1[complete.cases(result_df_maco1),]
View(filter_df_maco1)
```
# number of genes 16547

# padj < 0.05

```R
filter_df_maco1$padj < 0.05
filter_df2_maco1 = filter_df_maco1[filter_df_maco1$padj < 0.05, ]
```
# number of genes 2242

# log2FoldChange >1 <-1
```R
abs(filter_df2_maco1$log2FoldChange) > 1
filter_df3_maco1 = filter_df2_maco1[abs(filter_df2_maco1$log2FoldChange) > 1, ]
View(filter_df3_maco1)
```

# number of genes 776

# Assign gene names

```R
library("tidyverse")

filter_df_maco1 = rownames_to_column(filter_df_maco1, var = "ensgene")

listMarts()
ensembl_110 = useEnsembl(biomart = "ensembl", version = 110)
ensembl_110 = useDataset("celegans_gene_ensembl", mart = ensembl_110)

annotation_df_maco1 = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df_maco1$ensgene, mart = ensembl_110)
View(annotation_df_maco1)
annotated_df_maco1 = left_join(filter_df_maco1, annotation_df_maco1, by = c("ensgene" = "ensembl_gene_id"))

filter_df3_maco1 = rownames_to_column(filter_df3_maco1, var = "ensgene")
annotation_filter_df3_maco1 = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df3_maco1$ensgene, mart = ensembl_110)
annotated_df3_maco1 = left_join(filter_df3_maco1, annotation_filter_df3_maco1, by = c("ensgene" = "ensembl_gene_id"))

write_tsv(annotated_df3_maco1, "annotated_df3_maco1_1+4h")
```



