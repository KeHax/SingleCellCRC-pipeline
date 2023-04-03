library(Seurat)
library(tidyverse)
source("functions.r")

Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")
Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")
Imm_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Imm_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")

#epithelial cells-----
Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- FindVariableFeatures(Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge)
Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- Filtervariablefeatures(Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge,
                                                                                      pattern = "NBT")

Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- ScaleData(Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, vars.to.regress = c("pMT", "CC.Diff"))
Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- RunPCA(Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, npcs =50)
Epi_Integrate_harmony <- RunHarmony(Epi_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Epi_Integrate_harmony, ndims = 50)

#several runs of clustering to purify epithelial cells.
#clustering 1
Epi_Integrate_harmony <- FindNeighbors(Epi_Integrate_harmony, dims = 1:30)
Epi_Integrate_harmony <- RunUMAP(Epi_Integrate_harmony, dims = 1:30)
for (i in seq(0.1,1.2,length.out=12)) {
  Epi_Integrate_harmony <- FindClusters(Epi_Integrate_harmony, resolution = i)
}

#annotation
new.cluster.ids <- c("Epi","Epi","imm","Epi","Epi",
                     "Epi","Epi","Epi","Epi","Epi",
                     "imm","Epi","Epi","Epi","Epi",
                     "Epi","imm","Epi","imm","Epi",
                     "Epi","Epi","unkown","Epi","imm",
                     "Epi","Epi","Endo","Epi","Epi","imm","Epi")
Idents(Epi_Integrate_harmony) <- Epi_Integrate_harmony$RNA_snn_res.0.4
names(new.cluster.ids) <- levels(Epi_Integrate_harmony)
Epi_Integrate_harmony <- RenameIdents(Epi_Integrate_harmony,new.cluster.ids)
Epi_Integrate_harmony@meta.data$epi_cell_type_RNA_snn_res.0.4 <- Idents(Epi_Integrate_harmony)

Epi <- subset(Epi_Integrate_harmony, subset = epi_cell_type_RNA_snn_res.0.4 %in% "Epi")

#clustering 2
Epi_1 <- Epi
Epi_1 <- FindVariableFeatures(Epi_1)
Epi_1 <- Filtervariablefeatures(Epi_1, pattern = "NBT")

Epi_1 <- ScaleData(Epi_1, vars.to.regress = c("pMT", "CC.Diff"))
Epi_1 <- RunPCA(Epi_1, npcs =50)
Epi_1_harmony <- RunHarmony(Epi_1, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Epi_1_harmony, ndims = 50)

Epi_1_harmony <- FindNeighbors(Epi_1_harmony, reduction = "harmony", dims = 1:30)
Epi_1_harmony <- RunUMAP(Epi_1_harmony, reduction = "harmony", dims = 1:30)
Epi_1_harmony <- FindClusters(Epi_1_harmony, resolution = 1.0)


new.cluster.ids <- c("Epi_1","Epi_2","Epi_3","Epi_4","Epi_5","Epi_6",
                     "Epi_7","Epi_8","Epi_9","Epi_10","Epi_11",
                     "Epi_12","Epi_13","Epi_14","Epi_15","Epi_16",
                     "unkown","Epi_17","Epi_18","Epi_19","Epi_20",
                     "Epi_21","Epi_22","Epi_23","Epi_24","Epi_25")
Idents(Epi_1_harmony) <- Epi_1_harmony$RNA_snn_res.1
names(new.cluster.ids) <- levels(Epi_1_harmony)
Epi_1_harmony <- RenameIdents(Epi_1_harmony,new.cluster.ids)
Epi_1_harmony@meta.data$epi_cell_type_RNA_snn_res.1 <- Idents(Epi_1_harmony)

#clustering 3
Epi_2 <- Epi_1_harmony
Epi_2 <- subset(Epi_2, subset = `epi_cell_type_RNA_snn_res.1` != "unkown" )

Epi_2 <- FindVariableFeatures(Epi_2)
Epi_2 <- Filtervariablefeatures(Epi_2, pattern = "NBT")

Epi_2 <- ScaleData(Epi_2, vars.to.regress = c("pMT", "CC.Diff"))
Epi_2 <- RunPCA(Epi_2, npcs =50)
Epi_2_harmony <- RunHarmony(Epi_2, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Epi_2_harmony, ndims = 50)

Epi_2_harmony <- FindNeighbors(Epi_2_harmony, reduction = "harmony", dims = 1:30)
Epi_2_harmony <- RunUMAP(Epi_2_harmony, reduction = "harmony", dims = 1:30)
Epi_2_harmony <- FindClusters(Epi_2_harmony, resolution = 1.0)

new.cluster.ids <- c("Epi_1","Epi_2","Epi_3","Epi_4","Epi_5","Epi_6",
                     "Epi_7","Epi_8","Epi_9","Epi_10","Epi_11",
                     "Epi_12","Epi_13","Epi_14","Epi_15","Epi_16",
                     "Epi_17","Epi_18","Epi_19","Epi_20","Epi_21",
                     "Epi_22","Epi_23","Epi_24","Epi_25","Epi_26","Epi_27")
Idents(Epi_2_harmony) <- Epi_2_harmony$RNA_snn_res.1
names(new.cluster.ids) <- levels(Epi_2_harmony)
Epi_2_harmony <- RenameIdents(Epi_2_harmony,new.cluster.ids)
Epi_2_harmony@meta.data$epi_cell_type_RNA_snn_res.1 <- Idents(Epi_2_harmony)

#annotation
Epi_3 <- Epi_2_harmony
Epi_3 <- FindNeighbors(Epi_3, reduction = "harmony", dims = 1:30)
Epi_3 <- RunUMAP(Epi_3, reduction = "harmony", dims = 1:30)
Epi_3 <- FindClusters(Epi_3, resolution = 0.15)

Idents(Epi_3) <- Epi_3$RNA_snn_res.0.15
markers_0.15 <- FindAllMarkers(Epi_3,only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)
top50_markers_0.15 <- markers_0.15 %>% group_by(cluster) %>% top_n(50, wt = avg_log2FC)

Epi_3 <- subset(Epi_3,subset = RNA_snn_res.0.15 != "6")
Epi_3@meta.data$cell_type <- Epi_3@meta.data$RNA_snn_res.0.15

new.cluster.ids <- c("Malignant_epithelial_cells_1_TM4SF1",
                     "Malignant_epithelial_cells_2_SOX4",
                     "Malignant_epithelial_cells_3_proliferative",
                     "Malignant_epithelial_cells_4_ASCL2",
                     "Malignant_epithelial_cells_5_FCGBP",
                     "Malignant_epithelial_cells_6_OLFM4",
                     "Malignant_epithelial_cells_5_FCGBP")
Idents(Epi_3) <- Epi_3$cell_type
names(new.cluster.ids) <- levels(Epi_3)
Epi_3 <- RenameIdents(Epi_3,new.cluster.ids)
Epi_3@meta.data$cell_type <- Idents(Epi_3)

saveRDS(Epi_3,file = "/data/zdw/1_Integrate/CreateSeuratObject/Epi_6.RDS")

#stromal cells------
Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <-  readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")
Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$orig.ident <- as.character(Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge@meta.data$orig.ident)

Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- FindVariableFeatures(Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge)
Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- Filtervariablefeatures(Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge,
                                                                                      pattern="NBT")

Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- ScaleData(Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, vars.to.regress = c("pMT", "CC.Diff"))
Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <- RunPCA(Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, npcs =50)
Str_Integrate_harmony <- RunHarmony(Str_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Str_Integrate_harmony, ndims = 50)


Str_Integrate_harmony <- FindNeighbors(Str_Integrate_harmony, reduction = "harmony", dims = 1:30)
Str_Integrate_harmony <- RunUMAP(Str_Integrate_harmony, reduction = "harmony", dims = 1:30)

#several runs of clustering to purify stromal cells.
#clustering 1
for (i in seq(0.1,1.2,length.out=12)) {
  Str_Integrate_harmony <- FindClusters(Str_Integrate_harmony, resolution = i)
}

new.cluster.ids <- c("imm","Str","Str","Str","unknow","Str",
                     "imm","imm","Str","imm","imm",
                     "imm","unknow","imm","imm","imm",
                     "Str","imm","unknow","Str","unknow")
Idents(Str_Integrate_harmony) <- Str_Integrate_harmony$RNA_snn_res.0.4
names(new.cluster.ids) <- levels(Str_Integrate_harmony)
Str_Integrate_harmony <- RenameIdents(Str_Integrate_harmony,new.cluster.ids)
Str_Integrate_harmony@meta.data$Str_cell_type_RNA_snn_res.0.4 <- Idents(Str_Integrate_harmony)

Str <- subset(Str_Integrate_harmony, subset = Str_cell_type_RNA_snn_res.0.4 %in% "Str")

#clustering 2
Str_1 <- Str
Str_1 <- FindVariableFeatures(Str_1)
Str_1 <- Filtervariablefeatures(Str_1, pattern = "NBT")

Str_1 <- ScaleData(Str_1, vars.to.regress = c("pMT", "CC.Diff"))
Str_1 <- RunPCA(Str_1, npcs =50)
Str_1_harmony <- RunHarmony(Str_1, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Str_1_harmony, ndims = 50)

Str_1_harmony <- FindNeighbors(Str_1_harmony, reduction = "harmony", dims = 1:30)
Str_1_harmony <- RunUMAP(Str_1_harmony, reduction = "harmony", dims = 1:30)
Str_1_harmony <- FindClusters(Str_1_harmony, resolution = 1.0)

new.cluster.ids <- c("Str_1","Str_2","Str_3","Str_4","Str_5","Str_6",
                     "Str_7","Str_8","Str_9","Str_10","Str_11",
                     "Str_12","Str_13","unkown","Str_14","Str_15",
                     "Str_16","Str_17","Str_18","Str_19","unkown_1")
Idents(Str_1_harmony) <- Str_1_harmony$RNA_snn_res.1
names(new.cluster.ids) <- levels(Str_1_harmony)
Str_1_harmony <- RenameIdents(Str_1_harmony,new.cluster.ids)
Str_1_harmony@meta.data$Str_cell_type_RNA_snn_res.1 <- Idents(Str_1_harmony)

#clustering 3
Str_2 <- Str_1_harmony
Str_2 <- subset(Str_2, subset = `Str_cell_type_RNA_snn_res.1` != "unkown" )
Str_2 <- subset(Str_2, subset = `Str_cell_type_RNA_snn_res.1` != "unkown_1" )
Str_2@meta.data[["Str_cell_type_RNA_snn_res.1"]]=NULL

Str_2 <- FindVariableFeatures(Str_2)
Str_2 <- Filtervariablefeatures(Str_2, pattern = "NBT")

Str_2 <- ScaleData(Str_2, vars.to.regress = c("pMT", "CC.Diff"))
Str_2 <- RunPCA(Str_2, npcs =50)
Str_2_harmony <- RunHarmony(Str_2, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Str_2_harmony, ndims = 50)

Str_2_harmony <- FindNeighbors(Str_2_harmony, reduction = "harmony", dims = 1:25)
Str_2_harmony <- RunUMAP(Str_2_harmony, reduction = "harmony", dims = 1:25)
Str_2_harmony <- FindClusters(Str_2_harmony, resolution = 1.0)

new.cluster.ids <- c("Str_1","Str_2","Str_3","Str_4","Str_5","Str_6",
                     "Str_7","Str_8","Str_9","Str_10","Str_11",
                     "Str_12","Str_13","Str_14","Str_15","Str_16",
                     "Str_17","Str_18","Str_19")
Idents(Str_2_harmony) <- Str_2_harmony$RNA_snn_res.1
names(new.cluster.ids) <- levels(Str_2_harmony)
Str_2_harmony <- RenameIdents(Str_2_harmony,new.cluster.ids)
Str_2_harmony@meta.data$Str_cell_type_RNA_snn_res.1 <- Idents(Str_2_harmony)

#clustering 4
Str_4 <- Str_2_harmony

Str_4 <- FindVariableFeatures(Str_4)
Str_4 <- Filtervariablefeatures(Str_4, pattern = "NBT")

Str_4 <- ScaleData(Str_4, vars.to.regress = c("pMT", "CC.Diff"))
Str_4 <- RunPCA(Str_4, npcs =50)
Str_4_harmony <- RunHarmony(Str_4, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Str_4_harmony, ndims = 50)

Str_4_harmony <- FindNeighbors(Str_4_harmony, reduction = "harmony", dims = 1:30)
Str_4_harmony <- RunUMAP(Str_4_harmony, reduction = "harmony", dims = 1:30)
Str_4_harmony <- FindClusters(Str_4_harmony, resolution = 1.0)
Str_4_harmony <- FindClusters(Str_4_harmony, resolution = 0.1)

new.cluster.ids <- c("MSCs","Endothelial_cell","MSCs","MSCs","MSCs","MSCs","Endothelial_cell",
                     "MSCs","MSCs","MSCs","MSCs","MSCs","Endothelial_cell","MSCs")
Idents(Str_4_harmony) <- Str_4_harmony$RNA_snn_res.0.4
names(new.cluster.ids) <- levels(Str_4_harmony)
Str_4_harmony <- RenameIdents(Str_4_harmony,new.cluster.ids)
Str_4_harmony@meta.data$Str_cell_type_RNA_snn_res.0.4 <- Idents(Str_4_harmony)
Str_4_harmony@meta.data$cell_type <- Str_4_harmony@meta.data$Str_cell_type_RNA_snn_res.0.4

Str_Endo <- subset(Str_4_harmony, subset = Str_cell_type_RNA_snn_res.0.4 %in% c("Endo_1","Endo_2","Endo_3"))
Str_MSC <- subset(Str_4_harmony, subset = Str_cell_type_RNA_snn_res.0.4 %in% c("MSC_1","MSC_2","MSC_3","MSC_4","MSC_5","MSC_6","MSC_7","MSC_8","MSC_9","MSC_10","MSC_11"))

saveRDS(Str_Endo, file = "/data/zdw/1_Integrate/CreateSeuratObject/Str_Endo.RDS")
saveRDS(Str_MSC, file = "/data/zdw/1_Integrate/CreateSeuratObject/Str_MSC.RDS")

#clustering 5: endothelial cells
Str_Endo <- FindVariableFeatures(Str_Endo)
Str_Endo <- Filtervariablefeatures(Str_Endo, pattern = "NBT")

Str_Endo <- ScaleData(Str_Endo, vars.to.regress = c("pMT", "CC.Diff"))
Str_Endo <- RunPCA(Str_Endo, npcs =50)
Str_Endo_harmony <- RunHarmony(Str_Endo, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Str_Endo_harmony, ndims = 50)

Str_Endo_harmony <- FindNeighbors(Str_Endo_harmony, reduction = "harmony", dims = 1:25)
Str_Endo_harmony <- RunUMAP(Str_Endo_harmony, reduction = "harmony", dims = 1:25)
Str_Endo_harmony <- FindClusters(Str_Endo_harmony, resolution = 1)
#Endothelial cells only consisted of one cluster.
saveRDS(Str_Endo_harmony, file = "/data/zdw/1_Integrate/CreateSeuratObject/Str_Endo_harmony.RDS")

#clustering 6: MSC
Str_MSC <- FindVariableFeatures(Str_MSC)
Str_MSC <- Filtervariablefeatures(Str_MSC, pattern = "NBT")

Str_MSC <- ScaleData(Str_MSC, vars.to.regress = c("pMT", "CC.Diff"))
Str_MSC <- RunPCA(Str_MSC, npcs =50)
Str_MSC_harmony <- RunHarmony(Str_MSC, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Str_MSC_harmony, ndims = 50)

Str_MSC_harmony <- FindNeighbors(Str_MSC_harmony, reduction = "harmony", dims = 1:25)
Str_MSC_harmony <- RunUMAP(Str_MSC_harmony, reduction = "harmony", dims = 1:25)
Str_MSC_harmony <- FindClusters(Str_MSC_harmony, resolution = 0.02)
Str_MSC_harmony <- FindClusters(Str_MSC_harmony, resolution = 0.04)
Str_MSC_harmony <- FindClusters(Str_MSC_harmony, resolution = 0.06)
Str_MSC_harmony <- FindClusters(Str_MSC_harmony, resolution = 0.1)
Str_MSC_harmony <- FindClusters(Str_MSC_harmony, resolution = 0.4)

Str_MSC_harmony_1 <- subset(Str_MSC_harmony,subset = `RNA_snn_res.0.04` !="6")
new.cluster.ids <- c("h01_CAF_FAP",
                     "h02_PVL_MCAM",
                     "h03_CAF_CXCL14",
                     "h04_iCAF_CXCL12",
                     "h05_telocytes_ICAM1",
                     "h06_myoCAF_ACTA2")
Idents(Str_MSC_harmony_1) <- Str_MSC_harmony_1$RNA_snn_res.0.04
names(new.cluster.ids) <- levels(Str_MSC_harmony_1)
Str_MSC_harmony_1 <- RenameIdents(Str_MSC_harmony_1,new.cluster.ids)
Str_MSC_harmony_1@meta.data$Str_cell_type_RNA_snn_res.0.4 <- Idents(Str_MSC_harmony_1)

saveRDS(Str_MSC_harmony_1, file = "/data/zdw/1_Integrate/CreateSeuratObject/Str_MSC_harmony_1.RDS")

#immune cells------
Imm_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge <-  readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Imm_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge.RDS")

#Immune cells mixed in stromal cells and epithelial cells were also integrated
Str_Integrate_harmony <-readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Str_Integrate_harmony.RDS")
Epi_Integrate_harmony <-readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Epi_Integrate_harmony.RDS")

Imm_from_Str <- subset(Str_Integrate_harmony, subset = Str_cell_type_RNA_snn_res.0.4 %in% c("imm_1" ,"imm_2","imm_3","imm_4","imm_5","imm_6","imm_7","imm_8","imm_9","imm_10") )
Imm_from_Epi <- subset(Epi_Integrate_harmony, subset = epi_cell_type_RNA_snn_res.0.4 %in% c("imm_1" ,"imm_2","imm_3","imm_4","imm_5","imm_6") )

Imm_seu_obj <- merge(Imm_Integrate_seu_obj_2_with_Doubletfiltered_QC_Score_merge,y=c(Imm_from_Str,Imm_from_Epi))

Imm_seu_obj <- FindVariableFeatures(Imm_seu_obj)
tcr.gene <- read.csv("/data/zdw/1_Integrate/TCR_gene.csv")
tcr.gene <- tcr.gene$Approved.symbol
hvgs <- VariableFeatures(Imm_seu_obj)
hvgs <- hvgs[-grep("^IG[JKL]|^IGH[VDJ]|^RNA|^MT-|^RPS|^RPL",hvgs)]
hvgs <- hvgs[!hvgs %in% tcr.gene]
VariableFeatures(Imm_seu_obj) <- hvgs

Imm_seu_obj <- ScaleData(Imm_seu_obj, vars.to.regress = c("pMT", "CC.Diff"))
Imm_seu_obj <- RunPCA(Imm_seu_obj, npcs =50)
Imm_Integrate_harmony <- RunHarmony(Imm_seu_obj, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Imm_Integrate_harmony, ndims = 50)

Imm_Integrate_harmony <- FindNeighbors(Imm_Integrate_harmony, reduction = "harmony", dims = 1:25)
Imm_Integrate_harmony <- RunUMAP(Imm_Integrate_harmony, reduction = "harmony", dims = 1:25)

#clustering 1
for (i in seq(0.1,1.2,length.out=12)) {
  Str_Integrate_harmony <- FindClusters(Str_Integrate_harmony, resolution = i)
}

new.cluster.ids <- c("T/NKT_1","T/NKT_2","Mast","T/NKT_3","Epi","unknow_1","unknow_2",
                     "B_1","B_2","Plasma_1","T/NKT_4","Mye_1",
                     "Mye_2","T/NKT_5","Plasma_2","Mye_3")
Idents(Imm_Integrate_harmony) <- Imm_Integrate_harmony$RNA_snn_res.0.4
names(new.cluster.ids) <- levels(Imm_Integrate_harmony)
Imm_Integrate_harmony <- RenameIdents(Imm_Integrate_harmony,new.cluster.ids)
Imm_Integrate_harmony@meta.data$Imm_cell_type_RNA_snn_res.0.4 <- Idents(Imm_Integrate_harmony)

#clustering 2
Imm <- subset(Imm_Integrate_harmony, subset = `Imm_cell_type_RNA_snn_res.0.4` != "unknow_1")
Imm <- subset(Imm, subset = `Imm_cell_type_RNA_snn_res.0.4` != "unknow_2")
Imm <- subset(Imm, subset = `Imm_cell_type_RNA_snn_res.0.4` != "Epi")
Imm_1 <- Imm

Imm_1 <- FindVariableFeatures(Imm_1)
tcr.gene <- read.csv("/data/zdw/1_Integrate/TCR_gene.csv")
tcr.gene <- tcr.gene$Approved.symbol
hvgs <- VariableFeatures(Imm_1)
hvgs <- hvgs[-grep("^IG[JKL]|^IGH[VDJ]|^RNA|^MT-|^RPS|^RPL",hvgs)]
hvgs <- hvgs[!hvgs %in% tcr.gene]
VariableFeatures(Imm_1) <- hvgs

Imm_1 <- ScaleData(Imm_1, vars.to.regress = c("pMT", "CC.Diff"))
Imm_1 <- RunPCA(Imm_1, npcs =50)
Imm_1_harmony <- RunHarmony(Imm_1, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Imm_1_harmony, ndims = 50)

Imm_1_harmony <- FindNeighbors(Imm_1_harmony, reduction = "harmony", dims = 1:25)
Imm_1_harmony <- RunUMAP(Imm_1_harmony, reduction = "harmony", dims = 1:25)
Imm_1_harmony <- FindClusters(Imm_1_harmony, resolution = 1.0)

new.cluster.ids <- c("T/NKT_1","T/NKT_2","T/NKT_3","T/NKT_4","Mye_1","T/NKT_5",
                     "T/NKT_6","Plasma_1","Mast_1","Mye_2","B_1","Plasma_2",
                     "B_2","Mye_3","T/NKT_7","Mye_4","T/NKT_8","T/NKT_9",
                     "B_3","Plasma_3","T/NKT_10","T/NKT_11","Mye_5","Mye_6",
                     "T/NKT_10","T/NKT_11")
Idents(Imm_1_harmony) <- Imm_1_harmony$RNA_snn_res.1
names(new.cluster.ids) <- levels(Imm_1_harmony)
Imm_1_harmony <- RenameIdents(Imm_1_harmony,new.cluster.ids)
Imm_1_harmony@meta.data$Imm_cell_type_RNA_snn_res.1 <- Idents(Imm_1_harmony)

#clustering 3
Imm_3 <- Imm_1_harmony
Imm_3 <- FindVariableFeatures(Imm_3)
tcr.gene <- read.csv("/data/zdw/1_Integrate/TCR_gene.csv")
tcr.gene <- tcr.gene$Approved.symbol
hvgs <- VariableFeatures(Imm_3)
hvgs <- hvgs[-grep("^IG[JKL]|^IGH[VDJ]|^RNA|^MT-|^RPS|^RPL",hvgs)]
hvgs <- hvgs[!hvgs %in% tcr.gene]
VariableFeatures(Imm_3) <- hvgs

Imm_3 <- ScaleData(Imm_3, vars.to.regress = c("pMT", "CC.Diff"))
Imm_3 <- RunPCA(Imm_3, npcs =50)
Imm_3_harmony <- RunHarmony(Imm_3, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Imm_3_harmony, ndims = 50)

Imm_3_harmony <- FindNeighbors(Imm_3_harmony, reduction = "harmony", dims = 1:25)
Imm_3_harmony <- RunUMAP(Imm_3_harmony, reduction = "harmony", dims = 1:25)
Imm_3_harmony <- FindClusters(Imm_3_harmony, resolution = 0.02)

new.cluster.ids <- c("T_NK_cells","Myeloid_cells","B_cells","Plasma_cell","Mast_cell")
Idents(Imm_3_harmony) <- Imm_3_harmony$RNA_snn_res.0.02
names(new.cluster.ids) <- levels(Imm_3_harmony)
Imm_3_harmony <- RenameIdents(Imm_3_harmony,new.cluster.ids)
Imm_3_harmony@meta.data$Imm_cell_type_RNA_snn_res.0.02 <- Idents(Imm_3_harmony)

#subset immune cells
Imm_T_NK <- subset(Imm_3_harmony, subset = Imm_cell_type_RNA_snn_res.0.02 %in% "T_NK_cells")
Imm_Myeloid <- subset(Imm_3_harmony, subset = Imm_cell_type_RNA_snn_res.0.02 %in% "Myeloid_cells")
Imm_B_cells <- subset(Imm_3_harmony, subset = Imm_cell_type_RNA_snn_res.0.02 %in% "B_cells")
Imm_Plasma_cell <- subset(Imm_3_harmony, subset = Imm_cell_type_RNA_snn_res.0.02 %in% "Plasma_cell")
Imm_Mast_cell <- subset(Imm_3_harmony, subset = Imm_cell_type_RNA_snn_res.0.02 %in% "Mast_cell")

saveRDS(Imm_T_NK, file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_T_NK.RDS")
saveRDS(Imm_Myeloid, file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_Myeloid.RDS")
saveRDS(Imm_B_cells, file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_B_cells.RDS")
saveRDS(Imm_Plasma_cell, file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_Plasma_cell.RDS")
saveRDS(Imm_Mast_cell, file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_Mast_cell.RDS")


#T cell and NK cells--------------
Imm_T_NK <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Imm_T_NK.RDS")
Imm_T_NK <- FindVariableFeatures(Imm_T_NK)
Imm_T_NK <- Filtervariablefeatures(Imm_T_NK, pattern = "T_cell")

Imm_T_NK <- ScaleData(Imm_T_NK, vars.to.regress = c("pMT", "CC.Diff"))
Imm_T_NK <- RunPCA(Imm_T_NK, npcs =50)
Imm_T_NK_harmony <- RunHarmony(Imm_T_NK, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Imm_T_NK_harmony, ndims = 50)

Imm_T_NK_harmony <- FindNeighbors(Imm_T_NK_harmony, reduction = "harmony", dims = 1:20)
Imm_T_NK_harmony <- RunUMAP(Imm_T_NK_harmony, reduction = "harmony", dims = 1:20)
Imm_T_NK_harmony <- FindClusters(Imm_T_NK_harmony, resolution = 2.0)

new.cluster.ids <- c("CD4","CD8","CD8","CD8","CD8","CD8","CD4","NK_KLRD1","Pro_T","CD8",
                     "CD4","CD8","CD4","CD8","CD4","CD4","CD4","NK_KLRB1","CD4","CD4",
                     "CD4","CD4","CD4_CD8","CD4","CD4","CD4","CD8","CD4","CD4","NK_FCGR3A","CD4")
Idents(Imm_T_NK_harmony) <- Imm_T_NK_harmony$RNA_snn_res.2
names(new.cluster.ids) <- levels(Imm_T_NK_harmony)
Imm_T_NK_harmony <- RenameIdents(Imm_T_NK_harmony,new.cluster.ids)
Imm_T_NK_harmony@meta.data$Imm_cell_type_RNA_snn_res.2 <- Idents(Imm_T_NK_harmony)

Imm_T_NK_harmony@meta.data$cell_type <- Imm_T_NK_harmony@meta.data$Imm_cell_type_RNA_snn_res.2
new.cluster.ids <- c("CD4_cells",
                     "CD8_cells",
                     "NK_cells_KLRD1",
                     "T_cells_proliferative",
                     "NK_cells_KLRB1",
                     "CD4_CD8_cells",
                     "NK_cells_FCGR3A")
Idents(Imm_T_NK_harmony) <- Imm_T_NK_harmony$cell_type
names(new.cluster.ids) <- levels(Imm_T_NK_harmony)
Imm_T_NK_harmony <- RenameIdents(Imm_T_NK_harmony,new.cluster.ids)
Imm_T_NK_harmony@meta.data$cell_type_1 <- Idents(Imm_T_NK_harmony)

CD4_T <- subset(Imm_T_NK_harmony, cell_type_1 == "CD4_cells")
CD8_T <- subset(Imm_T_NK_harmony, cell_type_1 == "CD8_cells")

saveRDS(CD4_T,file="/data/zdw/1_Integrate/CreateSeuratObject/CD4_T.RDS")
saveRDS(CD8_T,file="/data/zdw/1_Integrate/CreateSeuratObject/CD8_T.RDS")


#CD4 T cells---------------
CD4_T <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/CD4_T.RDS")

CD4_T <- FindVariableFeatures(CD4_T)
CD4_T <- Filtervariablefeatures(CD4_T, pattern = "T_cell")

CD4_T <- ScaleData(CD4_T, vars.to.regress = c("pMT", "CC.Diff"))
CD4_T <- RunPCA(CD4_T, npcs =50)
CD4_T_harmony <- RunHarmony(CD4_T, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(CD4_T_harmony, ndims = 50)

CD4_T_harmony <- FindNeighbors(CD4_T_harmony, reduction = "harmony", dims = 1:25)
CD4_T_harmony <- RunUMAP(CD4_T_harmony, reduction = "harmony", dims = 1:25)
CD4_T_harmony <- FindClusters(CD4_T_harmony, resolution = 2.0)
CD4_T_harmony <- FindClusters(CD4_T_harmony, resolution = 1.0)
CD4_T_harmony <- FindClusters(CD4_T_harmony, resolution = 0.4)
CD4_T_harmony <- FindClusters(CD4_T_harmony, resolution = 0.1)

Idents(CD4_T_harmony) <- CD4_T_harmony$RNA_snn_res.2
new.cluster.ids <- c("hT01_CD4_Treg_FOXP3","hT01_CD4_Treg_FOXP3","hT02_CD4_Th1_CXCL13","hT03_CD4_Tn_CCR7","hT04_CD4_Tcm_ANXA1",
                     "hT05_CD4_Th17_IL17A","hT06_CD4_Tem_GZMK","hT03_CD4_Tn_CCR7","hT03_CD4_Tn_CCR7","hT03_CD4_Tn_CCR7",
                     "hT06_CD4_Tem_GZMK","hT03_CD4_Tn_CCR7","hT03_CD4_Tn_CCR7","hT03_CD4_Tn_CCR7","hT01_CD4_Treg_FOXP3",
                     "hT03_CD4_Tn_CCR7","hT03_CD4_Tn_CCR7","hT03_CD4_Tn_CCR7","hT01_CD4_Treg_FOXP3","hT03_CD4_Tn_CCR7",
                     "hT03_CD4_Tn_CCR7","hT03_CD4_Tn_CCR7","hT03_CD4_Tn_CCR7","hT07_CD4_Trm_CXCR6","hT03_CD4_Tn_CCR7",
                     "hT01_CD4_Treg_FOXP3","hT03_CD4_Tn_CCR7","hT01_CD4_Treg_FOXP3","hT07_CD4_Trm_CXCR6")

Idents(CD4_T_harmony) <- CD4_T_harmony$RNA_snn_res.2
names(new.cluster.ids) <- levels(CD4_T_harmony)
CD4_T_harmony <- RenameIdents(CD4_T_harmony, new.cluster.ids)
CD4_T_harmony@meta.data$cell_type_RNA_snn_res.2 <- Idents(CD4_T_harmony)

saveRDS(CD4_T_harmony,file="/data/zdw/1_Integrate/CreateSeuratObject/CD4_T_harmony.RDS")

#CD8 T cells---------------

CD8_3 <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/CD8_T.RDS")

CD8_3 <- FindVariableFeatures(CD8_3)
CD8_3 <- Filtervariablefeatures(CD8_3, pattern = "T_cell")

CD8_3 <- ScaleData(CD8_3, vars.to.regress = c("pMT", "CC.Diff"))
CD8_3 <- RunPCA(CD8_3, npcs =50)
CD8_3_harmony <- RunHarmony(CD8_3, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(CD8_3_harmony, ndims = 50)

CD8_3_harmony <- FindNeighbors(CD8_3_harmony, reduction = "harmony", dims = 1:20)
CD8_3_harmony <- RunUMAP(CD8_3_harmony, reduction = "harmony", dims = 1:20)
CD8_3_harmony <- FindClusters(CD8_3_harmony, resolution = 1.0)
CD8_3_harmony <- FindClusters(CD8_3_harmony, resolution = 0.4)

new.cluster.ids <- c("h01_CD8_TEMRA_TEFF",
                     "h02_CD8_HSPA1A",
                     "h03_CD8_Tem_CD161",
                     "h04_CD8_Tem_GZMK",
                     "h05_CD8_Tem_FGFBP2",
                     "h06_IEL_CD160",
                     "h07_CD8_Tn_CCR7",
                     "h08_CD8_Tex_PDCD1",
                     "h09_CD8_proliferative",
                     "h10_CD8_Trm_CXCL1",
                     "h11_CD8_SELENOK")
Idents(CD8_3_harmony) <- CD8_3_harmony$RNA_snn_res.0.4
names(new.cluster.ids) <- levels(CD8_3_harmony)
CD8_3_harmony <- RenameIdents(CD8_3_harmony, new.cluster.ids)
CD8_3_harmony@meta.data$cell_type_1 <- Idents(CD8_3_harmony)
CD8_3_harmony@meta.data$cell_type <- CD8_3_harmony@meta.data$cell_type_1

saveRDS(CD8_3_harmony,file="/data/zdw/1_Integrate/CreateSeuratObject/CD8_3_harmony.RDS")

#Myeloid cells-------------------
Imm_Myeloid <- readRDS(file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_Myeloid.RDS")
Imm_Myeloid <- FindVariableFeatures(Imm_Myeloid)
Imm_Myeloid <- Filtervariablefeatures(Imm_Myeloid, pattern = "NBT")

Imm_Myeloid <- ScaleData(Imm_Myeloid, vars.to.regress = c("pMT", "CC.Diff"))
Imm_Myeloid <- RunPCA(Imm_Myeloid, npcs =50)
Imm_Myeloid_harmony <- RunHarmony(Imm_Myeloid, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(Imm_Myeloid_harmony, ndims = 50)

Imm_Myeloid_harmony <- FindNeighbors(Imm_Myeloid_harmony, reduction = "harmony", dims = 1:25)
Imm_Myeloid_harmony <- RunUMAP(Imm_Myeloid_harmony, reduction = "harmony", dims = 1:25)
Imm_Myeloid_harmony <- FindClusters(Imm_Myeloid_harmony, resolution = 0.2)

Imm_Myeloid_harmony <- subset(Imm_Myeloid_harmony,subset = RNA_snn_res.0.2 !="8")
new.cluster.ids <- c("hM01_Mph_C1QC",
                     "hM02_cDC2_FCN1",
                     "hM03_Mono_CX3CR1",
                     "hM04_Mono_NLRP3",
                     "hM05_Mono_IL1B",
                     "hM06_cDC2_CD1C",
                     "hM07_Mph_CD68",
                     "hM08_Mph_MKI67",
                     "hM9_Mph_SPP1",
                     "hM10_cDC_LAMPP3")
Idents(Imm_Myeloid_harmony) <- Imm_Myeloid_harmony$RNA_snn_res.0.2
names(new.cluster.ids) <- levels(Imm_Myeloid_harmony)
Imm_Myeloid_harmony <- RenameIdents(Imm_Myeloid_harmony,new.cluster.ids)
Imm_Myeloid_harmony@meta.data$cell_type_0.2 <- Idents(Imm_Myeloid_harmony)

Imm_Myeloid_harmony@meta.data$cell_type <- Imm_Myeloid_harmony@meta.data$cell_type_0.2
saveRDS(Imm_Myeloid_harmony, file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_Myeloid_harmony.RDS")

#DC is heterogeneous, recluster DC
mcell <- Imm_Myeloid_harmony
DC <- subset(mcell, Cluster %in% c("hM06_cDC2_CD1C", "hM10_cDC_LAMPP3"))
Seuratflow2 <- function(object,
                        pattern,
                        group.by.vars="orig.ident",
                        max.iter.harmony=30,
                        dims=15,
                        resolution=0.6){
  cat("harmony workflow.\ndata must bave been normalized")
  DefaultAssay(object) <- "RNA"
  object <- object %>% FindVariableFeatures() %>% ScaleData(features=rownames(object))
  object <- Filtervariablefeatures(object,pattern=pattern)
  object <- RunPCA(object)
  object <- RunHarmony(object, group.by.vars = group.by.vars, max.iter.harmony = max.iter.harmony)
  object <- object %>% FindNeighbors(dims=1:dims, reduction="harmony") %>% 
    FindClusters(resolution=resolution) %>% RunUMAP(dims=1:dims, reduction="harmony")
  return(object)
}
DC <- Seuratflow2(DC, pattern = "NBT")

mcell$Cluster <- gsub("M02_cDC2_FCN1","M02_Mph_FCN1",mcell$Cluster)
mcell$Cluster <- gsub("M03_Mono_CX3CR1","M01_Mono_CX3CR1",mcell$Cluster)
mcell$Cluster <- gsub("M04_Mono_NLRP3","M02_Mono_NLRP3",mcell$Cluster)
mcell$Cluster <- gsub("M05_Mono_IL1B","M03_Mono_IL1B",mcell$Cluster)
mcell$Cluster <- gsub("M01_Mph_C1QC","M04_Mph_C1QC",mcell$Cluster)
mcell$Cluster <- gsub("M02_Mph_FCN1","M05_Mph_FCN1",mcell$Cluster)
mcell$Cluster <- gsub("M07_Mph_CD68","M06_Mph_CD68",mcell$Cluster)
mcell$Cluster <- gsub("M08_Mph_MKI67","M07_Mph_MKI67",mcell$Cluster)
mcell$Cluster <- gsub("M09_Mph_SPP1","M08_Mph_SPP1",mcell$Cluster)
#M09_cDC1_CLEC9A
mcell$Cluster <- gsub("M06_cDC2_CD1C","M10_cDC2_CD1C",mcell$Cluster)
mcell$Cluster <- gsub("M10_cDC_LAMPP3","M11_cDC_LAMPP3",mcell$Cluster)
#M12_pDC_IRF7
mcell$Cluster <- gsub("M11_Mast_cell","M13_Mast_cell",mcell$Cluster)
mcell$Cluster <- gsub("M11_cDC_LAMPP3","M11_cDC_LAMP3",mcell$Cluster)

mcell@meta.data$Cluster[match(colnames(DC)[DC$seurat_clusters==7], colnames(mcell))] <- rep("M09_cDC1_CLEC9A", length(colnames(DC)[DC$seurat_clusters==7]))
mcell@meta.data$Cluster[match(colnames(DC)[DC$seurat_clusters==13], colnames(mcell))] <- rep("M12_pDC_IRF7", length(colnames(DC)[DC$seurat_clusters==13]))

saveRDS(mcell, "/data/khaoxian/project/Stagecrc/data/final_data/myeloid.corrected.rds")

#B cells---------------
Imm_B_cells <- readRDS(file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_B_cells.RDS")
Imm_Plasma_cell <- readRDS(file = "/data/zdw/1_Integrate/CreateSeuratObject/Imm_Plasma_cell.RDS")

B_Plasma <- merge(Imm_B_cells,Imm_Plasma_cell)
B_Plasma <- FindVariableFeatures(B_Plasma)
B_Plasma <- Filtervariablefeatures(B_Plasma, pattern = "B_cell")

B_Plasma <- ScaleData(B_Plasma, vars.to.regress = c("pMT", "CC.Diff"))
B_Plasma <- RunPCA(B_Plasma, npcs =50)
B_Plasma_harmony <- RunHarmony(B_Plasma, group.by.vars = "orig.ident", max.iter.harmony = 20)
ElbowPlot(B_Plasma_harmony, ndims = 50)

B_Plasma_harmony <- FindNeighbors(B_Plasma_harmony, reduction = "harmony", dims = 1:30)
B_Plasma_harmony <- RunUMAP(B_Plasma_harmony, reduction = "harmony", dims = 1:30)
B_Plasma_harmony <- FindClusters(B_Plasma_harmony, resolution = 0.4)

B_Plasma_harmony_1 <- subset(B_Plasma_harmony,subset = RNA_snn_res.0.4 !="5")
B_Plasma_harmony_1 <- subset(B_Plasma_harmony_1,subset = RNA_snn_res.0.4 !="9")

Idents(B_Plasma_harmony_1) <- B_Plasma_harmony_1$RNA_snn_res.0.4
new.cluster.ids <- c("hB03_FollicularB_MS4A1","hB01_PlasmaB_IgA","hB02_PlasmaB_IgG","hB01_PlasmaB_IgA","hB04_naiveB_IGHD","hB02_PlasmaB_IgG","hB03_FollicularB_MS4A1","hB01_GCBCell_AICDA","hB02_GCBCell_LRMP")
names(new.cluster.ids) <- levels(B_Plasma_harmony_1)
B_Plasma_harmony_1 <- RenameIdents(B_Plasma_harmony_1,new.cluster.ids)
B_Plasma_harmony_1@meta.data$cell_type <- Idents(B_Plasma_harmony_1)

saveRDS(B_Plasma_harmony_1, file = "/data/zdw/1_Integrate/CreateSeuratObject/B_Plasma_harmony_1.RDS")









