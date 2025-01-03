# FOR WT_1h versus ST7_4h

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

sample_table = read.csv("Info_samples_DEseq2.csv", sep=";")
sample_table_ST7_4h = sample_table[c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30), ] 
sample_files = paste0("./kallistoKatya/", pull(sample_table, "Sample") , ".1/abundance.h5")
names(sample_files) = pull(sample_table, "Sample")
names(sample_files)

sample_files_ST7_4h = sample_files[c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)]
count_data_ST7_4h <- tximport(files = sample_files_ST7_4h,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE)
count_data_ST7_4h[['counts']]
```

# specify what samples are control
```R
conditions = sample_table_ST7_4h[,3]
conditions = factor(conditions, levels = c('WT_4h', 'ST7:AID_4h'))
sample_table_ST7_4h$conditions = conditions

deseq_dataset_ST7_4h = DESeqDataSetFromTximport(txi = count_data_ST7_4h, colData = sample_table_ST7_4h, design = ~conditions)
```
# estimate size factors (normalizatiion)

```R
deseq_dataset_ST7_4h = estimateSizeFactors(deseq_dataset_ST7_4h)
normalizationFactors(deseq_dataset_ST7_4h)
counts(deseq_dataset_ST7_4h, normalized=TRUE)

vst = varianceStabilizingTransformation(deseq_dataset_ST7_4h)
plotPCA(vst, intgroup='conditions') + geom_label(aes(label=name))

d = assay(vst) 
d = t(d) 
d = dist(d) 
h = hclust(d) 
plot(h)
```

# estimate dispersions

```R
deseq_dataset_ST7_4h = estimateDispersions(deseq_dataset_ST7_4h)
plotDispEsts(deseq_dataset_ST7_4h)
```
# Run statistics
```R
deseq_dataset_ST7_4h = nbinomWaldTest(deseq_dataset_ST7_4h) 
result_table_ST7_4h = results(deseq_dataset_ST7_4h) 
summary(result_table_ST7_4h)
```
out of 18830 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 868, 4.6%
LFC < 0 (down)     : 274, 1.5%
outliers [1]       : 107, 0.57%
low counts [2]     : 3939, 21%
(mean count < 4)


# MA plot
```R
plotMA(result_table_ST7_4h)
```

# Make results to data.frame

```R
result_df_ST7_4h = as.data.frame(result_table_ST7_4h)

sum(complete.cases(result_df_ST7_4h))

filter_df_ST7_4h = result_df_ST7_4h[complete.cases(result_df_ST7_4h),]
View(filter_df_ST7_4h)
```
# number of genes 14784

# padj < 0.05

```R
filter_df_ST7_4h$padj < 0.05
filter_df2_ST7_4h = filter_df_ST7_4h[filter_df_ST7_4h$padj < 0.05, ]
```
# number of genes 648

# log2FoldChange >1 <-1
```R
abs(filter_df2_ST7_4h$log2FoldChange) > 1
filter_df3_ST7_4h = filter_df2_ST7_4h[abs(filter_df2_ST7_4h$log2FoldChange) > 1, ]
View(filter_df3_ST7_4h)
```

# number of genes 258

# Assign gene names

```R
library("tidyverse")

filter_df_ST7_4h = rownames_to_column(filter_df_ST7_4h, var = "ensgene")

listMarts()
ensembl_110 = useEnsembl(biomart = "ensembl", version = 110)
ensembl_110 = useDataset("celegans_gene_ensembl", mart = ensembl_110)

annotation_df_ST7_4h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df_ST7_4h$ensgene, mart = ensembl_110)
View(annotation_df_ST7_4h)
annotated_df_ST7_4h = left_join(filter_df_ST7_4h, annotation_df_ST7_4h, by = c("ensgene" = "ensembl_gene_id"))

filter_df3_ST7_4h = rownames_to_column(filter_df3_ST7_4h, var = "ensgene")
annotation_filter_df3_ST7_4h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df3_ST7_4h$ensgene, mart = ensembl_110)
annotated_df3_ST7_4h = left_join(filter_df3_ST7_4h, annotation_filter_df3_ST7_4h, by = c("ensgene" = "ensembl_gene_id"))

write_tsv(annotated_df3_ST7_4h, "annotated_df3_ST7_4h")
```

# Common genes for ST7 between 1h and 4h
```R
annotated_df3_ST7_1h <- read.csv("/Users/elashman/Documents/Niko_Hanna_RNASeq_13Nov24/annotated_df3_ST7_1h.csv", sep="\t")
genes_ST7_1h = annotated_df3_ST7_1h$ensgene
genes_ST7_4h = annotated_df3_ST7_4h$ensgene
common_genes_ST7_1h_4h = intersect(genes_ST7_1h, genes_ST7_4h)

list_common_genes_ST7_1h = annotated_df3_ST7_1h[annotated_df3_ST7_1h$ensgene %in% common_genes_ST7_1h_4h, c('ensgene','log2FoldChange', 'padj')]
list_common_genes_ST7_4h = annotated_df3_ST7_4h[annotated_df3_ST7_4h$ensgene %in% common_genes_ST7_1h_4h, c('ensgene','log2FoldChange', 'padj', 'external_gene_name', 'description')]

colnames(list_common_genes_ST7_1h)[2] = "log2FoldChange_1h"
colnames(list_common_genes_ST7_4h)[2] = "log2FoldChange_4h"
colnames(list_common_genes_ST7_1h)[3] = "padj_1h"
colnames(list_common_genes_ST7_4h)[3] = "padj_4h"
List_common_genes_ST7_padj_FoldChange = left_join(list_common_genes_ST7_1h, list_common_genes_ST7_4h, by = c("ensgene" = "ensgene")) 
write_tsv(List_common_genes_ST7_padj_FoldChange, "List_common_genes_ST7_padj_FoldChange")
```
