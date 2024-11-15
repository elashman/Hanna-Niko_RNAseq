# FOR WT_1h versus maco1_4h
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
sample_table_maco1_4h = sample_table[c(19, 20, 21, 22, 23, 24, 31, 32, 33, 34, 35, 36), ] 
sample_files = paste0("./kallistoKatya/", pull(sample_table, "Sample") , ".1/abundance.h5")
names(sample_files) = pull(sample_table, "Sample")
names(sample_files)

sample_files_maco1_4h = sample_files[c(19, 20, 21, 22, 23, 24, 31, 32, 33, 34, 35, 36)]
count_data_maco1_4h <- tximport(files = sample_files_maco1_4h,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE)
count_data_maco1_4h[['counts']]
```

# specify what samples are control
```R
conditions = sample_table_maco1_4h[,3]
conditions = factor(conditions, levels = c('WT_4h', 'MACO-1:AID_4h'))
sample_table_maco1_4h$conditions = conditions

deseq_dataset_maco1_4h = DESeqDataSetFromTximport(txi = count_data_maco1_4h, colData = sample_table_maco1_4h, design = ~conditions)
```
# estimate size factors (normalizatiion)

```R
deseq_dataset_maco1_4h = estimateSizeFactors(deseq_dataset_maco1_4h)
normalizationFactors(deseq_dataset_maco1_4h)
counts(deseq_dataset_maco1_4h, normalized=TRUE)

vst = varianceStabilizingTransformation(deseq_dataset_maco1_4h)
plotPCA(vst, intgroup='conditions') + geom_label(aes(label=name))
```

# estimate dispersions

```R
deseq_dataset_maco1_4h = estimateDispersions(deseq_dataset_maco1_4h)
plotDispEsts(deseq_dataset_maco1_4h)
```
# Run statistics
```R
deseq_dataset_maco1_4h = nbinomWaldTest(deseq_dataset_maco1_4h) 
result_table_maco1_4h = results(deseq_dataset_maco1_4h) 
summary(result_table_maco1_4h)
```
out of 18618 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 955, 5.1%
LFC < 0 (down)     : 1176, 6.3%
outliers [1]       : 19, 0.1%
low counts [2]     : 3564, 19%
(mean count < 3)


# MA plot
```R
plotMA(result_table_maco1_4h)
```

# Make results to data.frame

```R
result_df_maco1_4h = as.data.frame(result_table_maco1_4h)

sum(complete.cases(result_df_maco1_4h))

filter_df_maco1_4h = result_df_maco1_4h[complete.cases(result_df_maco1_4h),]
View(filter_df_maco1_4h)
```
# in filter_df_maco1_4h 15035 genes

# padj < 0.05

```R
filter_df_maco1_4h$padj < 0.05
filter_df2_maco1_4h = filter_df_maco1_4h[filter_df_maco1_4h$padj < 0.05, ]
```
# in filter_df2_maco1_4h 1633 genes

# log2FoldChange >1 <-1
```R
abs(filter_df2_maco1_4h$log2FoldChange) > 1
filter_df3_maco1_4h = filter_df2_maco1_4h[abs(filter_df2_maco1_4h$log2FoldChange) > 1, ]
View(filter_df3_maco1_4h)
```
# in filter_df3_maco1_4h 505 genes


# Assign gene names

```R
library("tidyverse")

filter_df_maco1_4h = rownames_to_column(filter_df_maco1_4h, var = "ensgene")

listMarts()
ensembl_110 = useEnsembl(biomart = "ensembl", version = 110)
ensembl_110 = useDataset("celegans_gene_ensembl", mart = ensembl_110)

annotation_df_maco1_4h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df_maco1_4h$ensgene, mart = ensembl_110)
View(annotation_df_maco1_4h)
annotated_df_maco1_4h = left_join(filter_df_maco1_4h, annotation_df_maco1_4h, by = c("ensgene" = "ensembl_gene_id"))

filter_df3_maco1_4h = rownames_to_column(filter_df3_maco1_4h, var = "ensgene")
annotation_filter_df3_maco1_4h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df3_maco1_4h$ensgene, mart = ensembl_110)
annotated_df3_maco1_4h = left_join(filter_df3_maco1_4h, annotation_filter_df3_maco1_4h, by = c("ensgene" = "ensembl_gene_id"))

write_tsv(annotated_df3_maco1_4h, "annotated_df3_maco1_4h")
```

# Common genes for maco-1 between 1h and 4h
```R
annotated_df3_maco1_1h <- read.csv("/Users/elashman/Documents/Niko_Hanna_RNASeq_13Nov24/annotated_df3_maco1_1h.csv", sep="\t")
genes_maco1_1h = annotated_df3_maco1_1h$ensgene
genes_maco1_4h = annotated_df3_maco1_4h$ensgene
common_genes_maco1_1h_4h = intersect(genes_maco1_1h, genes_maco1_4h)

list_common_genes_maco1_1h = annotated_df3_maco1_1h[annotated_df3_maco1_1h$ensgene %in% common_genes_maco1_1h_4h, c('ensgene','log2FoldChange', 'padj')]
list_common_genes_maco1_4h = annotated_df3_maco1_4h[annotated_df3_maco1_4h$ensgene %in% common_genes_maco1_1h_4h, c('ensgene','log2FoldChange', 'padj', 'external_gene_name', 'description')]

colnames(list_common_genes_maco1_1h)[2] = "log2FoldChange_1h"
colnames(list_common_genes_maco1_4h)[2] = "log2FoldChange_4h"
colnames(list_common_genes_maco1_1h)[3] = "padj_1h"
colnames(list_common_genes_maco1_4h)[3] = "padj_4h"
List_common_genes_maco1_padj_FoldChange = left_join(list_common_genes_maco1_1h, list_common_genes_maco1_4h, by = c("ensgene" = "ensgene")) 
write_tsv(List_common_genes_maco1_padj_FoldChange, "List_common_genes_maco1_padj_FoldChange")
```
