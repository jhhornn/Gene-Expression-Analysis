# scripts to visualise gene expression data (GSE183947)
# setwd("~/Documents/Bioinformatics Projects/Bioinformatics_101")

# load libraries
library(tidyverse)
library(ggplot2)

# basic format for ggplot
#gplot(data, aes(x = variable,  y = variable1)) + 
#  geom_col()

# Tasks
# Visualise gene expression (RNA-Seq) data as -
# 1. Bar plots
# 2. Density plots
# 3. Boxplots
# 4. Scatterplots
# 5. Heatmap


# barplot
data.long %>% 
  filter(gene == "BRCA1") %>% 
  ggplot(., aes(samples, FPKM, fill = tissue)) +
  geom_col()

# densityplot
data.long %>% 
  filter(gene == "BRCA1") %>% 
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.3)

# boxplot
data.long %>% 
  filter(gene == "BRCA1") %>% 
  ggplot(., aes(x = metastasis, y = FPKM)) +
  #geom_boxplot()
  geom_violin()

# scatterplot
# check for correlation between the expression of 2 genes
data.long %>% 
  filter(gene == "BRCA1" | gene == "BRCA2") %>% 
  spread()
  
  