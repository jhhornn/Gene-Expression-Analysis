# Vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow
# Data: Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications (GSE52778)

# Study Design
# RNA-Sequencing performed on 4 primary human airway smooth muscle
# cell lines treated with 1 micromolar dexamethasone for 18 hours

# For each of the 4 cell lines, study has a treated and untreated sample

# Goal: To understand the transcriptional changes occurring
# due to treatment with dexamethasone((used by asthma pateient to reduce
# inflammation in the airways).


# Requirements
## R
## RStudio

# R packages
## DESeq2
## tidyverse
## airway

# Script to perform differential gene expression analysis using DESeq2 package
# setwd("~/Documents/Bioinformatics Projects/Bioinformatics_101/Gene Expression")

# load libraries
library(DESeq2)
library(tidyverse)
library(airway) # the package we load the data we want to use from

# Script to load dataset to be used
data(airway)
airway
colData(airway)

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

## Step 1: preparing count data ...

# read in counts data
counts_data <- read.csv("counts_data.csv")
head(counts_data)

# read in sample info
colData <- read.csv("sample_info.csv")

# making sure the row names in colData matches to column names in counts.data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order
all(colnames(counts_data) == rownames(colData))


## Step 2: construct a DESeqDataset object ...
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# set the factor level
# set untreated as the reference level so we can compare it 
# with the treated
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")
# reference level is written first. If we don't explicitly state ref,
# ref is chosen in alphabetical order

# NOTE: if our dataset has technical replicates, we need to
# collapse them before we perform Differential Expression Analysis. Using a
# function called collapse replicates from DESeq2 package
# ONLY COLLAPSE TECHNICAL REPLICATES AND NOT BIOLOGICAL REPLICATES


## Step 3: Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

# res breakdown
# baseMean: the mean of each counts of a gene across all samples.
# average of the normalized counts taken over all the sample

# log2FoldChange: the fold change of the gene in the treated condition when compared with the untreated
# The poitive values stand for the unregulated gene in the treated condition while the
# negative values stand for the down regulated genes in the treated condition

# lfcSE: provides the standard error estimates for the log2FoldChange

# stat: the wald test values for the genes

# pvalue: the pvalue of the test statistic for the gene

# padj: p adjusted value is the corrected p value for multiple testing. To avoid
# the detection of false positive genes, the DESeq2 does adjustment on the p values


## Explore Results ...
summary(res)

res0.01 <- results(dds, alpha = 0.01)
# alpha = pvalue
summary(res0.01)


# contrasts
resultsNames(dds)

# e.g, if we have treated_4hrs, treated_8hrs as multiple levels of comparison
results(dds, constrast = c("dexamethasone", "treated_4hrs", "untreated"))


# MA plot
plotMA(res)


