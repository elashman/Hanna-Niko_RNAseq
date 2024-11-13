FOR WT_1h versus ST7_1h
Exclude all the samples that are not needed

sample_files_ST7_1h = sample_files[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)] count_data_ST7_1h <- tximport(files = sample_files_ST7_1h,type = "kallisto",tx2gene = gene_map, ignoreAfterBar = TRUE) count_data_ST7_1h[['counts']]
specify what samples are control

sample_table_ST7_1h = read.csv("Info_samples_DEseq2_ST7_1h.csv", sep=";") conditions = sample_table_ST7_1h[,3] conditions = factor(conditions, levels = c('WT_1h', 'ST7:AID_1h')) sample_table_ST7_1h$conditions = conditions

deseq_dataset_ST7_1h = DESeqDataSetFromTximport(txi = count_data_ST7_1h, colData = sample_table_ST7_1h, design = ~conditions)
estimate size factors (normalizatiion)

deseq_dataset_ST7_1h = estimateSizeFactors(deseq_dataset_ST7_1h) normalizationFactors(deseq_dataset_ST7_1h) counts(deseq_dataset_ST7_1h, normalized=TRUE)
estimate dispersions

deseq_dataset_ST7_1h = estimateDispersions(deseq_dataset_ST7_1h) plotDispEsts(deseq_dataset_ST7_1h)
Run statistics

deseq_dataset_ST7_1h = nbinomWaldTest(deseq_dataset_ST7_1h) result_table_ST7_1h = results(deseq_dataset_ST7_1h) summary(result_table_ST7_1h)

Result: out of 18854 with nonzero total read count adjusted p-value < 0.1 LFC > 0 (up) : 608, 3.2% LFC < 0 (down) : 428, 2.3% outliers [1] : 7, 0.037% low counts [2] : 4704, 25% (mean count < 9)
MA plot

plotMA(result_table_ST7_1h)
Make results to data.frame

result_df_ST7_1h = as.data.frame(result_table_ST7_1h)

sum(complete.cases(result_df_ST7_1h))

filter_df_ST7_1h = result_df_ST7_1h[complete.cases(result_df_ST7_1h),] View(filter_df_ST7_1h)
in filter_df_ST7_1h 14143 genes
padj < 0.05

filter_df_ST7_1h$padj < 0.05 filter_df2_ST7_1h = filter_df_ST7_1h[filter_df_ST7_1h$padj < 0.05, ]
in filter_df2_ST7_1h 579 genes
log2FoldChange >1 <-1

abs(filter_df2_ST7_1h$log2FoldChange) > 1 filter_df3_ST7_1h = filter_df2_ST7_1h[abs(filter_df2_ST7_1h$log2FoldChange) > 1, ] View(filter_df3_ST7_1h)
in filter_df3_ST7_1h 120 genes
Assign gene names

library("tidyverse")

filter_df_ST7_1h = rownames_to_column(filter_df_ST7_1h, var = "ensgene")

listMarts() ensembl_110 = useEnsembl(biomart = "ensembl", version = 110) ensembl_110 = useDataset("celegans_gene_ensembl", mart = ensembl_110)

annotation_df_ST7_1h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df_ST7_1h$ensgene, mart = ensembl_110) View(annotation_df_ST7_1h) annotated_df_ST7_1h = left_join(filter_df_ST7_1h, annotation_df_ST7_1h, by = c("ensgene" = "ensembl_gene_id"))

filter_df3_ST7_1h = rownames_to_column(filter_df3_ST7_1h, var = "ensgene") annotation_filter_df3_ST7_1h = getBM (attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "external_gene_name", "description"), filters = c("ensembl_gene_id"), values = filter_df3_ST7_1h$ensgene, mart = ensembl_110) annotated_df3_ST7_1h = left_join(filter_df3_ST7_1h, annotation_filter_df3_ST7_1h, by = c("ensgene" = "ensembl_gene_id"))

write_tsv(annotated_df3_ST7_1h, "annotated_df3_ST7_1h")
