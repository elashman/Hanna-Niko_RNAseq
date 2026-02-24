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



