#cytoTRACE analysis
#using mnn corrected data
library(Seurat)
library(tidyverse)
library(CytoTRACE)
source("function.r")
setwd("/data/khaoxian/project/Stagecrc")

epi.mnn <- anndata::read_h5ad("/data/khaoxian/project/Stagecrc/save/MNNcorrect/epi.mnncorrect.h5ad")
epi.mnn.data <- epi.mnn$X
epi.mnn.data <- t(epi.mnn.data)

epi_cytotracez <- CytoTRACEz(epi.mnn.data)
epi.mnn.data.pheno <- as.character(epi.mnn$obs$Cluster)
names(epi.mnn.data.pheno) <- rownames(epi.mnn$obs)

dir.create("/data/khaoxian/project/Stagecrc/output/cytotrace/epi/1/",recursive = T)
plotCytoTRACE(epi_cytotracez, phenotype = epi.mnn.data.pheno,
              outputDir ="/data/khaoxian/project/Stagecrc/output/cytotrace/epi/1/")
plotCytoGenes(epi_cytotracez, numOfGenes = 20,
              outputDir = "/data/khaoxian/project/Stagecrc/output/cytotrace/epi/1/")
saveRDS(epi_cytotracez,"/data/khaoxian/project/Stagecrc/output/cytotrace/epi/1/epi_cytotracez.rds")