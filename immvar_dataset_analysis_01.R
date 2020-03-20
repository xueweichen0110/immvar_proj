#installation of required package variancePartition
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("variancePartition")
#
rm(list = ls())

library(BiocManager)
library(variancePartition)
library(dplyr)
#read datasets
exprobject <- readRDS('exprObj.RDS')
info_exprobject <- readRDS('info.RDS')
####info_exprobject overview
head(info_exprobject)
#984 genes, 6 batches, 

###exprobject overview
head(exprobject)
featuredata_table <- exprobject@featureData@data
