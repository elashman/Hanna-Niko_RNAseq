# FOR WT_1h versus maco1_1h
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

sample_table = read.csv("Info_samples_DEseq2.csv", sep=";")
sample_table_maco1_1h = sample_table[c(1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18), ] 
sample_files = paste0("./kallistoKatya/", pull(sample_table, "Sample") , ".1/abundance.h5")
names(sample_files) = pull(sample_table, "Sample")
names(sample_files)

sample_files_maco1_1h = sample_files[c(1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18)]
count_data_maco1_1h <- tximport(files = sample_files_maco1_1h,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE)
count_data_maco1_1h[['counts']]
```

# specify what samples are control
```R
conditions = sample_table_maco1_1h[,3]
conditions = factor(conditions, levels = c('WT_1h', 'MACO-1:AID_1h'))
sample_table_maco1_1h$conditions = conditions

deseq_dataset_maco1_1h = DESeqDataSetFromTximport(txi = count_data_maco1_1h, colData = sample_table_maco1_1h, design = ~conditions)
```
# estimate size factors (normalizatiion)

```R
deseq_dataset_maco1_1h = estimateSizeFactors(deseq_dataset_maco1_1h)
normalizationFactors(deseq_dataset_maco1_1h)
counts(deseq_dataset_maco1_1h, normalized=TRUE)
```

# estimate dispersions

```R
deseq_dataset_maco1_1h = estimateDispersions(deseq_dataset_maco1_1h)
plotDispEsts(deseq_dataset_maco1_1h)
```
# Run statistics
```R
deseq_dataset_maco1_1h = nbinomWaldTest(deseq_dataset_maco1_1h) 
result_table_maco1_1h = results(deseq_dataset_maco1_1h) 
summary(result_table_maco1_1h)
```
out of 18929 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1823, 9.6%
LFC < 0 (down)     : 2182, 12%
outliers [1]       : 50, 0.26%
low counts [2]     : 2546, 13%
(mean count < 1)


# MA plot
```R
plotMA(result_table_maco1_1h)
```

# Make results to data.frame

```R
result_df_maco1_1h = as.data.frame(result_table_maco1_1h)

sum(complete.cases(result_df_maco1_1h))

filter_df_maco1_1h = result_df_maco1_1h[complete.cases(result_df_maco1_1h),]
View(filter_df_maco1_1h)
```
# in filter_df_maco1_1h 16333 genes

# padj < 0.05

```R
filter_df_maco1_1h$padj < 0.05
filter_df2_maco1_1h = filter_df_maco1_1h[filter_df_maco1_1h$padj < 0.05, ]
```
# in filter_df2_maco1_1h 3274 genes

# log2FoldChange >1 <-1
```R
abs(filter_df2_maco1_1h$log2FoldChange) > 1
filter_df3_maco1_1h = filter_df2_maco1_1h[abs(filter_df2_maco1_1h$log2FoldChange) > 1, ]
View(filter_df3_maco1_1h)
```
# in filter_df3_maco1_1h 1153 genes


# Assign gene names

```R
library("tidyverse")

filter_df_maco1_1h = rownames_to_column(filter_df_maco1_1h, var = "ensgene")

listMarts()
ensembl_110 = useEnsembl(biomart = "ensembl", version = 110)
ensembl_110 = useDataset("celegans_gene_ensembl", mart = ensembl_110)

annotation_df_maco1_1h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df_maco1_1h$ensgene, mart = ensembl_110)
View(annotation_df_maco1_1h)
annotated_df_maco1_1h = left_join(filter_df_maco1_1h, annotation_df_maco1_1h, by = c("ensgene" = "ensembl_gene_id"))

filter_df3_maco1_1h = rownames_to_column(filter_df3_maco1_1h, var = "ensgene")
annotation_filter_df3_maco1_1h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df3_maco1_1h$ensgene, mart = ensembl_110)
annotated_df3_maco1_1h = left_join(filter_df3_maco1_1h, annotation_filter_df3_maco1_1h, by = c("ensgene" = "ensembl_gene_id"))

write_tsv(annotated_df3_maco1_1h, "annotated_df3_maco1_1h")
