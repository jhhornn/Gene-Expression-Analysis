# This script manipulates gene expression data
# setwd("~/Documents/Bioinformatics Projects/Bioinformatics_101")

# Tasks
# Read in expression matrix (FPKM)
# Fetching metadata/clinical data using GEOquery R package
# Demonstrate useful functions to tidy, reshape, join and explore data

#Packages Requires
# dplyr(v1.0.7)
# tidyverse(v1.3.1)
# GEOquery(v2.60.0)

# load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

# read in the data
gse_data <- read.csv(file = "./GSE183947_fpkm.csv")
dim(gse_data)

# get metadate
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# metadata_subset <- select(metadata, c(1,10,11,17))


metadata.modified <- metadata %>% 
  select(1,10,11,17) %>% 
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>% 
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(gse_data)


# reshaping data
data.long <- gse_data %>% 
  rename(gene = X) %>% 
  gather(key = "samples", value = "FPKM", -gene)

# join dataframes = data.long + metadata.modified
data.long <- data.long %>% 
  left_join(., metadata.modified, by = c("samples" = "description"))


# explore data
data.long %>% 
  filter(gene == "BRCA1" | gene == "BRCA2") %>% 
  group_by(gene, tissue) %>% 
  summarize(mean_FPKM = mean(FPKM))

