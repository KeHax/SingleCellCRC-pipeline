library(Seurat)
library(tidyverse)
library(monocle3)
source("function.r")

cd4T <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/CD4_T_harmony.RDS") 
cd4T$Cluster <- factor(cd4T$Cluster, levels = sort(unique(cd4T$Cluster)))
cd4t.cds <- seurat2monocle3(cd4T)
cd4t.cds <- preprocess_cds(cd4t.cds, num_dim = 50)
cd4t.cds <- reduce_dimension(cd4t.cds, preprocess_method = "PCA")

cds.embed <- cd4t.cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(cd4T, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cd4t.cds@int_colData$reducedDims$UMAP <- int.embed

cds.embed <- cd4t.cds@int_colData$reducedDims$PCA
int.embed <- Embeddings(cd4T, reduction = "pca")
int.embed <- int.embed[rownames(cds.embed),]
cd4t.cds@int_colData$reducedDims$PCA <- int.embed

cd4t.cds <- cluster_cells(cd4t.cds)
levels(clusters(cd4t.cds,reduction_method = "UMAP"))

cd4t.cds@clusters@listData[["UMAP"]][["clusters"]] <- cd4T$Cluster
levels(clusters(cd4t.cds,reduction_method = "UMAP"))

cd4t.cds <- learn_graph(cd4t.cds,close_loop=F)
cd4t.cds <- order_cells(cd4t.cds)   #choose naive CD4+ T cell

cd4T$pseudotime <- pseudotime(cd4t.cds)

