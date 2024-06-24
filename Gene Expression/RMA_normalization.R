# Script to perform RMA normalization
# setwd("/Users/oluwaseun/Documents/Bioinformatics Projects/Bioinformatics_101/Gene Expression")

# How to read and normalize(RMA) microarray data in R

## Requirements
# R
# RStudio

## R packages:
# tidyverse
# GEOquery (fetch data from GEO DB)
# affy (function that performs the RMA normalization)


# load libraries
library(affy)
library(GEOquery)
library(tidyverse)

# get supplementary files
#getGEOSuppFiles("GSE148537")

# untar files
untar("GSE148537_RAW.tar", exdir = "data/")

# reading in .cel files
raw.data <- ReadAffy(celfile.path = "data/")
raw.data

# performing RMA normalization
normalized.data <- rma(raw.data)
normalized.data

# get expression estimate
normalized.expr <- as.data.frame(exprs(normalized.data))

# map probe IDs to gene symbols
gse <- getGEO("GSE148537", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE148537_series_matrix.txt.gz@featureData@data
# subset
feature.data <- feature.data[,c(1,11)]

normalized.expr <- normalized.expr %>% 
  rownames_to_column(var = "ID") %>% 
  inner_join(., feature.data, by = "ID")

