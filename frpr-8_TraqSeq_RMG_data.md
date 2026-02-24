# The goal is to check the mRNA levels of the gene frpr-8 (WBGene00018728) in TrapSeq data for RMG from Hanna and Niko

# TRAP elution

# WT
273769 (1 File)
273768 (1 File)
273767 (1 File)

# ST7

273772 (1 File)
273771 (1 File)
273770 (1 File)

# Maco

273775 (1 File)
273774 (1 File)
273773 (1 File)

# load libraries in RStudio

```bash
setwd("/Users/elashman/Documents/NikoHanna_TrapSeq")

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
# load data 
```bash
 sample_table = read.csv("Info_samples_DEseq2.csv", sep=";")
 sample_files = paste0("./kallisto/", pull(sample_table, "Sample") , "/abundance.h5")
 names(sample_files) = pull(sample_table, "Sample") 

 count_data <- tximport(files = sample_files,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE)
 count_data[['counts']]
```
# specify what samples are control

```bash
 conditions = sample_table[,3]
 conditions = factor(conditions, levels = c('WT_elution', 'ST7_elution', 'Macoilin_elution', 'WT_input', 'ST7_input', 'Macoilin_input' )) 
 sample_table$conditions = conditions

 deseq_dataset = DESeqDataSetFromTximport(txi = count_data, colData =sample_table, design = ~conditions)
```
  ## 3 STEPS in DESeq2 analysis

# STEP 1: estimate size factors (normalizatiion)

```bash
deseq_dataset = estimateSizeFactors(deseq_dataset)
normalizationFactors(deseq_dataset)
counts(deseq_dataset, normalized=TRUE)

vst = varianceStabilizingTransformation(deseq_dataset)
plotPCA(vst, intgroup='conditions')
```

# FOR WT versus maco-1
# Exclude all the samples that are not needed

```bash
sample_files_maco1_elution = sample_files[c(1, 2, 3, 7, 8, 9)]
count_data_maco1_elution <- tximport(files = sample_files_maco1_elution,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE)
```

# specify what samples are control
```bash
 sample_table_maco1_elution = read.csv("Info_samples_DEseq2_maco1_elution.csv", sep=";")
 conditions =  sample_table_maco1_elution[,3]
 conditions = factor(conditions, levels = c('WT_elution', 'Macoilin_elution')) 
 sample_table_maco1_elution$conditions = conditions

  deseq_dataset_maco1_elution = DESeqDataSetFromTximport(txi = count_data_maco1_elution, colData = sample_table_maco1_elution, design = ~conditions)
```

# estimate size factors (normalizatiion)
```bash
deseq_dataset_maco1_elution = estimateSizeFactors(deseq_dataset_maco1_elution)
normalizationFactors(deseq_dataset_maco1_elution)
counts(deseq_dataset_maco1_elution, normalized=TRUE)

vst_maco1_elution = varianceStabilizingTransformation(deseq_dataset_maco1_elution)
plotPCA(vst_maco1_elution, intgroup='conditions')
```
# estimate dispersions
```bash
deseq_dataset_maco1_elution = estimateDispersions(deseq_dataset_maco1_elution)
plotDispEsts(deseq_dataset_maco1_elution)
```

# Run statistics 
```bash
deseq_dataset_maco1_elution = nbinomWaldTest(deseq_dataset_maco1_elution)
result_table_maco1_elution = results(deseq_dataset_maco1_elution)
summary(result_table_maco1_elution)

out of 20010 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 8, 0.04%
LFC < 0 (down)     : 4, 0.02%
outliers [1]       : 114, 0.57%
low counts [2]     : 3823, 19%
(mean count < 16)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

# result_table_analysis2 is DataFrame not data.frame! To make it data.frame
```bash
result_df_maco1_elution = as.data.frame(result_table_maco1_elution)
write.csv(result_df_maco1_elution, "FC_all_genes_maco1_elution.csv", row.names = TRUE)
```

# plot counts for  WBGene00003807 (npr-1)
```bash
dds = deseq_dataset_maco1_elution
plotCounts(dds, 'WBGene00003807', intgroup ='conditions', normalized = TRUE)
plotCounts(dds, 'WBGene00003807', intgroup ='conditions', normalized = FALSE)
```

# plot counts for frpr-8 (WBGene00018728)
```bash
dds = deseq_dataset_maco1_elution
plotCounts(dds, 'WBGene00018728', intgroup ='conditions', normalized = TRUE)
plotCounts(dds, 'WBGene00018728', intgroup ='conditions', normalized = FALSE)
```

# FOR WT versus ST7

# Exclude all the samples that are not needed
```bash
sample_files_ST7_elution = sample_files[c(4, 5, 6, 7, 8, 9)]
count_data_ST7_elution <- tximport(files = sample_files_ST7_elution,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE)
```
# specify what samples are control

```bash
 sample_table_ST7_elution = read.csv("Info_samples_DEseq2_ST7input.csv", sep=";")
 conditions = sample_table_ST7_elution[,3]
 conditions = factor(conditions, levels = c('WT_elution', 'ST7_elution')) 
 sample_table_ST7_elution$conditions = conditions

  deseq_dataset_ST7_elution = DESeqDataSetFromTximport(txi = count_data_ST7_elution, colData = sample_table_ST7_elution, design = ~conditions)
```

# estimate size factors (normalizatiion)
```bash
deseq_dataset_ST7_elution = estimateSizeFactors(deseq_dataset_ST7_elution)
normalizationFactors(deseq_dataset_ST7_elution)
counts(deseq_dataset_ST7_elution, normalized=TRUE)

vst_ST7_elution = varianceStabilizingTransformation(deseq_dataset_ST7_elution)
plotPCA(vst_ST7_elution, intgroup='conditions')
```
# estimate dispersions
```bash
deseq_dataset_ST7_elution = estimateDispersions(deseq_dataset_ST7_elution)
plotDispEsts(deseq_dataset_ST7_elution)
```
# Run statistics 
```bash
deseq_dataset_ST7_elution = nbinomWaldTest(deseq_dataset_ST7_elution)
result_table_ST7_elution = results(deseq_dataset_ST7_elution)
summary(result_table_ST7_elution)

out of 19719 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 0, 0%
LFC < 0 (down)     : 2, 0.01%
outliers [1]       : 202, 1%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

# result_table_analysis2 is DataFrame not data.frame! To make it data.frame
```bash
result_df_ST7_elution = as.data.frame(result_table_ST7_elution)
write.csv(result_df_ST7_elution, "FC_all_genes_ST7_elution.csv", row.names = TRUE)
```
# plot counts for  WBGene00003807
```bash
dds = deseq_dataset_ST7_elution
plotCounts(dds, 'WBGene00003807', intgroup ='conditions', normalized = TRUE)
plotCounts(dds, 'WBGene00003807', intgroup ='conditions', normalized = FALSE)
```

# plot counts for frpr-8 (WBGene00018728)
```bash
dds = deseq_dataset_ST7_elution
plotCounts(dds, 'WBGene00018728', intgroup ='conditions', normalized = TRUE)
plotCounts(dds, 'WBGene00018728', intgroup ='conditions', normalized = FALSE)
```
