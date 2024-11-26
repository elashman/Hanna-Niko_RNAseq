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
sample_table_ST7 = sample_table[c(1,2,3,4,5,6,7,8,9,10,11,12,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30), ] 
sample_files = paste0("./kallistoKatya/", pull(sample_table, "Sample") , ".1/abundance.h5")
names(sample_files) = pull(sample_table, "Sample")
names(sample_files)

sample_files_ST7 = sample_files[c(1,2,3,4,5,6,7,8,9,10,11,12,19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)]
count_data_ST7 <- tximport(files = sample_files_ST7,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE)
count_data_ST7[['counts']]
```

# specify what samples are control
```R
conditions = sample_table_ST7[,3]
conditions = factor(conditions, levels = c('WT', 'ST7:AID'))
sample_table_ST7$conditions = conditions

deseq_dataset_ST7 = DESeqDataSetFromTximport(txi = count_data_ST7, colData = sample_table_ST7, design = ~conditions)
```
# estimate size factors (normalizatiion)

```R
deseq_dataset_ST7 = estimateSizeFactors(deseq_dataset_ST7)
normalizationFactors(deseq_dataset_ST7)
counts(deseq_dataset_ST7, normalized=TRUE)

vst = varianceStabilizingTransformation(deseq_dataset_ST7)
plotPCA(vst, intgroup='conditions') + geom_label(aes(label=name))

d = assay(vst) 
d = t(d) 
d = dist(d) 
h = hclust(d) 
plot(h)
```
![PCAplot_ST7_1+4h](https://github.com/user-attachments/assets/dd4cad29-0747-42f7-a7c6-6c6667472fbe)

![Cluster_Dendrogram_ST7_1+4h](https://github.com/user-attachments/assets/0b7df487-4077-45ae-b963-0787ec50edb0)



# estimate dispersions

```R
deseq_dataset_ST7 = estimateDispersions(deseq_dataset_ST7)
plotDispEsts(deseq_dataset_ST7)
```
![Dispersion_ST7_1+4h](https://github.com/user-attachments/assets/3fa85be3-b58c-4b38-a50f-32aa50c5ce03)

# Run statistics
```R
deseq_dataset_ST7 = nbinomWaldTest(deseq_dataset_ST7) 
result_table_ST7 = results(deseq_dataset_ST7) 
summary(result_table_ST7)
```
out of 19340 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 362, 1.9%
LFC < 0 (down)     : 155, 0.8%
outliers [1]       : 207, 1.1%
low counts [2]     : 3286, 17%
(mean count < 2)


# MA plot
```R
plotMA(result_table_ST7)
```
![MAplot_ST7_1+4h](https://github.com/user-attachments/assets/0fc525f5-fa8b-4366-b20a-b89fbb1565b5)

# Make results to data.frame

```R
result_df_ST7 = as.data.frame(result_table_ST7)

sum(complete.cases(result_df_ST7))

filter_df_ST7 = result_df_ST7[complete.cases(result_df_ST7),]
View(filter_df_ST7)
```
# number of genes 15847

# padj < 0.05

```R
filter_df_ST7$padj < 0.05
filter_df2_ST7 = filter_df_ST7[filter_df_ST7$padj < 0.05, ]
```
# number of genes 326

# log2FoldChange >1 <-1
```R
abs(filter_df2_ST7$log2FoldChange) > 1
filter_df3_ST7 = filter_df2_ST7[abs(filter_df2_ST7$log2FoldChange) > 1, ]
View(filter_df3_ST7)
```

# number of genes 89

# Assign gene names

```R
library("tidyverse")

filter_df_ST7 = rownames_to_column(filter_df_ST7, var = "ensgene")

listMarts()
ensembl_110 = useEnsembl(biomart = "ensembl", version = 110)
ensembl_110 = useDataset("celegans_gene_ensembl", mart = ensembl_110)

annotation_df_ST7 = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df_ST7$ensgene, mart = ensembl_110)
View(annotation_df_ST7)
annotated_df_ST7 = left_join(filter_df_ST7, annotation_df_ST7, by = c("ensgene" = "ensembl_gene_id"))

filter_df3_ST7 = rownames_to_column(filter_df3_ST7, var = "ensgene")
annotation_filter_df3_ST7 = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df3_ST7$ensgene, mart = ensembl_110)
annotated_df3_ST7 = left_join(filter_df3_ST7, annotation_filter_df3_ST7, by = c("ensgene" = "ensembl_gene_id"))

write_tsv(annotated_df3_ST7, "annotated_df3_ST7_1+4h")
```


