library(reticulate)
library(Matrix)
library(Seurat)
library(tidyverse)
library(harmony)
source("functions.r")
#1. merge data-----
GSE188711_seu_data <- readRDS("/data/zdw/GSE188711/CreateSeuratObject/GSE188711_seu_data.RDS")
GSE161277_seu_data <- readRDS("/data/zdw/GSE161277/CreateSeuratObject/GSE161277_seu_data.RDS")
GSE132465_seu_data_T <- readRDS("/data/zdw/GSE132465/CreateSeuratObject/GSE132465_seu_data_T.RDS")
GSE144735_T <- readRDS("/data/zdw/GSE144735/CreateSeuratObject/GSE144735_T.RDS")
GSE146771_seu_data <- readRDS("/data/zdw/GSE146771/CreateSeuratObject/GSE146771_seu_data.RDS")
GSE164522_CRLM_PT_seu_data <- readRDS("/data/zdw/GSE164522/CreateSeuratObject/GSE164522_CRLM_PT_seu_data.RDS")
GSE178318_seu_data_T <- readRDS("/data/zdw/GSE178318/CreateSeuratObject/GSE178318_seu_data_T.RDS")
GSE200997_seu_data_T <- readRDS("/data/zdw/GSE178318/CreateSeuratObject/GSE200997_seu_data_T.RDS")

Integrate_seu_obj_1 <- merge(GSE188711_seu_data, y = c(GSE161277_seu_data,GSE132465_seu_data_T,
                                                       GSE144735_T,GSE146771_seu_data,GSE164522_CRLM_PT_seu_data,
                                                       GSE178318_seu_data_T,GSE200997_seu_data_T))

#some data have been normalized and transformed by lo2p1 and need re-transformation with logep1
norm.orig <- c("P0104", "P0123", "P0202", "P0305", "P0323", "P0408", "P0410", "P0613",
               "P1025", "P1026","patient08", "patient09", "patient10", "patient11", "patient12",
               "patient13", "patient14", "patient15", "patient16", "patient17")
Integrate_seu_obj_1 <- NormalizeData(Integrate_seu_obj_1)
Integrate_seu_obj_1 <- Correctnormdata(Integrate_seu_obj_1,normalized.samples=norm.orig)

Integrate_seu_obj_2 <- SplitObject(Integrate_seu_obj_1, split.by = "orig.ident")

#2. remove doublets-----
#py_install("scrublet", pip = T)
#py_install("numpy", pip = T)

Integrate_seu_obj_2_list = list()
for (i in 1:length(Integrate_seu_obj_2)) {
  matrixs <- GetAssayData(Integrate_seu_obj_2[[i]], assay="RNA", slot="counts")
  matrixs <- t(matrixs)
  barcodes <- rownames(matrixs)
  Integrate_seu_obj_2_list[[i]]=list(matrixs,barcodes)
  names(Integrate_seu_obj_2_list[[i]])[1]="matrix"
  names(Integrate_seu_obj_2_list[[i]])[2]="barcodes"
  names(Integrate_seu_obj_2_list)[i]=names(Integrate_seu_obj_2)[i]
}


colnames = c("sample", "number","rate")
Multiplet_rates <- matrix(nrow = 81,ncol = 3,dimnames = list(NULL, colnames))

for (i in 1:length(Integrate_seu_obj_2)) {
  Multiplet_rates[i,1] <- names(Integrate_seu_obj_2_list)[i]
  Multiplet_rates[i,2] <- as.numeric(nrow(Integrate_seu_obj_2_list[[i]][["matrix"]]) )
}
a=as.numeric(Multiplet_rates[,2])
for (i in 1:length(a)) {
  Multiplet_rates[i,3]=ifelse(a[i] < 1000,0.008,ifelse(a[i] <2000,0.015,ifelse(a[i] <3000,0.023,ifelse(a[i] <4000,0.030,
                                                                                                       ifelse(a[i] <5000,0.038,ifelse(a[i] <6000,0.046,ifelse(a[i] <7000,0.053,ifelse(a[i] <8000,0.062,
                                                                                                                                                                                      ifelse(a[i] <9000,0.068,ifelse(a[i] <10000,0.080,0.08))))))))))
}  

#do not run these codes in cycles, and set i = 1~81 step by step
for (i in 1:length(Integrate_seu_obj_2_list)) {
  matrixs=Integrate_seu_obj_2_list[[i]][["matrix"]]
  barcodes=Integrate_seu_obj_2_list[[i]][["barcodes"]]
  names=names(Integrate_seu_obj_2_list)[i]
  a=Multiplet_rates[i,3]
  rate=as.numeric(a)
  
  #run in python console
  repl_python()
  # import matplotlib.pyplot as plt
  # import scrublet as scr
  # import scipy.io
  # import numpy as np
  # import os
  # import pandas as pd
  # scrub = scr.Scrublet(r.matrixs, expected_doublet_rate=r.rate) 
  # doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
  # scrub.plot_histogram()
  # plt.show()
  # 
  # #run in console, and threshold could be adjusted according to the results of scrub.plot_histogram
  # predicted_doublets_new =scrub.call_doublets(threshold=0.38)
  # scrub.plot_histogram()
  # plt.show()
  # 
  # sample_id=r.names
  # input_dir='/home/zdw/Scrublet_Doublet'
  # out_df=pd.DataFrame(r.barcodes)
  # out_df['doublet_scores'] = doublet_scores
  # out_df['predicted_doublets'] = predicted_doublets_new
  # os.mkdir(input_dir + '/'+ sample_id)
  # out_df.to_csv(input_dir + '/' + sample_id + '/' + 'doublet.txt',index=False,header=True)
  # exit
  
  #Run in R: subset singlets
  list <- py$out_df
  names(list) <- c("cell_id","doublet_scores","predicted_doublets")
  metadata <- FetchData(Integrate_seu_obj_2_Doubletfiltered[[names]], "orig.ident")
  metadata$cell_id <- rownames(metadata)
  metadata <- left_join(x = metadata, y = list, by = "cell_id")
  rownames(metadata) <- metadata$cell_id
  metadata$cell_id=NULL
  Integrate_seu_obj_2_Doubletfiltered[[names]] <- AddMetaData(Integrate_seu_obj_2_Doubletfiltered[[names]], metadata = metadata)
  A <- paste(names,"_Doublet_list", sep = "")
  saveRDS(list,file = paste("/data/zdw/1_Integrate/Scrublet_Doublet/",A,".RDS",sep="")) 
}

Integrate_seu_obj_2_Doubletfiltered <- Integrate_seu_obj_2
saveRDS(Integrate_seu_obj_2_Doubletfiltered, file = "/data/zdw/1_Integrate/CreateSeuratObject/Integrate_seu_obj_2_Doubletfiltered.RDS")

#3. quality control-----
nFeature_lower <- 250
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 20

A <- Integrate_seu_obj_2_with_Doubletfiltered
for (i in 1:length(A)) {
  A[[i]] <- PercentageFeatureSet(A[[i]], pattern = "^MT-", col.name = "pMT")
  A[[i]] <- subset(A[[i]], subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper 
                   & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & pMT < pMT_upper)
}

Integrate_seu_obj_2_with_Doubletfiltered_QC <- A
saveRDS(Integrate_seu_obj_2_with_Doubletfiltered_QC, file = "/data/zdw/1_Integrate/CreateSeuratObject/Integrate_seu_obj_2_with_Doubletfiltered_QC.RDS")

#4. clustering cells to 3 clusters by scores-----
epithelial_cell_marker <- list(c("EPCAM", "KRT8", "KRT18"))
stromal_cell_marker <- list(c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B","HAND2"))
immune_cell_marker <- list(c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14",  "CD68", "CD83", "CSF1R", "FCER1G","LYZ"))

Integrate_seu_obj_2_with_Doubletfiltered_QC <- map(Integrate_seu_obj_2_with_Doubletfiltered_QC,function(x){
  x <- AddModuleScore(x,features = epithelial_cell_marker,name = "epithelial_cell",nbin = 5)
  x <- AddModuleScore(x,features = stromal_cell_marker,name = "stromal_cell",nbin = 5)
  x <- AddModuleScore(x,features = immune_cell_marker,name = "immune_cell",nbin = 5)
  return(x)
})

Integrate_seu_obj_2_with_Doubletfiltered_QC <- map(Integrate_seu_obj_2_with_Doubletfiltered_QC,function(x){
  metas <- x@meta.data %>% as.data.frame() %>% select(epithelial_cell1, stromal_cell1, immune_cell1)
  metasx <- apply(metas,1,function(x){
    colnames(metas)[match(max(x),x)]
  })
  x@meta.data$cell_type <- metasx
  return(x)
})

Integrate_seu_obj_2_with_Doubletfiltered_QC_Score <- Integrate_seu_obj_2_with_Doubletfiltered_QC
saveRDS(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score, file = "/data/zdw/1_Integrate/CreateSeuratObject/Integrate_seu_obj_2_with_Doubletfiltered_QC_Score.RDS")

Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- merge(x=Integrate_seu_obj_2_with_Doubletfiltered_QC_Score[[1]],
                                                                 y=Integrate_seu_obj_2_with_Doubletfiltered_QC_Score[2:81])

Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$AAA <- Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@active.ident
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge$AAA <- gsub(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$AAA, pattern = "-", replacement = "_", fixed = TRUE)
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$orig.ident <- Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$AAA
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$orig.ident <- as.factor(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$orig.ident)

metas <- readxl::read_excel("/data/zdw/1_Integrate/Metadata/metadata.xlsx")
metadata <- FetchData(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$Patient <- metadata$orig.ident
metadata <- left_join(x = metadata, y = metas, by = "Patient")
rownames(metadata) <- metadata$cell_id
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- AddMetaData(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, metadata = metadata)
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- subset(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, subset = Treatment.history %in% "Treatment-naïve")

#scale data and reduce dimension
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- FindVariableFeatures(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
score_cc <- function(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge) {
  Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- CellCycleScoring(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, s.genes, g2m.genes)
  Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$CC.Diff <- Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$S.Score - Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$G2M.Score
  return(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge)
}
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- score_cc(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge)

#remove mitochondiral, ribosomal genes, TCR and BCR
tcr.gene <- read.csv("/data/zdw/1_Integrate/TCR_gene.csv")
tcr.gene <- tcr.gene$Approved.symbol
hvgs <- VariableFeatures(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge)
hvgs <- hvgs[-grep("^IG[JKLH]|^RNA|^MT-|^RPS|^RPL",hvgs)]
hvgs <- hvgs[!hvgs %in% tcr.gene]
VariableFeatures(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge) <- hvgs

Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- ScaleData(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, vars.to.regress = c("pMT", "CC.Diff"))
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- RunPCA(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, npcs =50)
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$orig.ident <- as.character(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$orig.ident)
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- RunHarmony(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, ndims = 50)
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- FindNeighbors(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, reduction = "harmony", dims = 1:30)
Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- RunUMAP(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, reduction = "harmony", dims = 1:30)

saveRDS(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, file = "/data/zdw/1_Integrate/CreateSeuratObject/Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")

#3. subset major cluster
Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- subset(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, subset = cell_type %in% "stromal_cell1")
Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- subset(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, subset = cell_type %in% "epithelial_cell1")
Imm_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- subset(Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, subset = cell_type %in% "immune_cell1")

saveRDS(Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, file = "/data/zdw/1_Integrate/CreateSeuratObject/Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")
saveRDS(Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, file = "/data/zdw/1_Integrate/CreateSeuratObject/Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")
saveRDS(Imm_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")








