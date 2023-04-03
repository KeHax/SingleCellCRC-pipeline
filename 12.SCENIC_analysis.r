library(SeuratDisk)
library(Seurat)
library(tidyverse)
setwd("/data/khaoxian/project/Stagecrc")

#myeloid cells-----------------------
mcell <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Imm_Myeloid_harmony.RDS")

#remove datasets which did not provide counts data
norm.orig <- c("P0104", "P0123", "P0202", "P0305", "P0323", "P0408", "P0410", "P0613",
               "P1025", "P1026","patient08", "patient09", "patient10", "patient11", "patient12",
               "patient13", "patient14", "patient15", "patient16", "patient17")
mcellx <- subset(mcell, orig.ident %in% setdiff(unique(mcell$orig.ident), norm.orig))

mcells <- Samplingseurat(mcellx, "Cluster", 500)
mcells <- RunTSNE(mcells, reduction = "harmony")
mcells@assays$RNA@scale.data <- matrix()
mcells@assays$RNA@data <- mcells@assays$RNA@counts

SaveH5Seurat(mcells, filename = "./save/senic_data/mcells.h5Seurat", overwrite = T)
Convert("./save/senic_data/mcells.h5Seurat", dest = "h5ad", overwrite = T)
system("/home/khaoxian/anaconda3/envs/pyscenic/bin/python /home/khaoxian/pycharm/scenic/scenic.py /data/khaoxian/project/Stagecrc/output/newsenic/myeloid /data/khaoxian/project/Stagecrc/save/senic_data/mcells.h5ad Cluster")
#done

#str-----------------------
str <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/Str_MSC_harmony_1.RDS")

#remove datasets which did not provide counts data
norm.orig <- c("P0104", "P0123", "P0202", "P0305", "P0323", "P0408", "P0410", "P0613",
               "P1025", "P1026","patient08", "patient09", "patient10", "patient11", "patient12",
               "patient13", "patient14", "patient15", "patient16", "patient17")
strx <- subset(str, orig.ident %in% setdiff(unique(str$orig.ident), norm.orig))

strs <- Samplingseurat(strx, "Cluster", 500)
strs <- RunTSNE(strs, reduction = "harmony")
strs@assays$RNA@scale.data <- matrix()
strs@assays$RNA@data <- strs@assays$RNA@counts

SaveH5Seurat(strs, filename = "./save/senic_data/strs.h5Seurat", overwrite = T)
Convert("./save/senic_data/strs.h5Seurat", dest = "h5ad", overwrite = T)
system("/home/khaoxian/anaconda3/envs/pyscenic/bin/python /home/khaoxian/pycharm/scenic/scenic.py /data/khaoxian/project/Stagecrc/output/newsenic/str /data/khaoxian/project/Stagecrc/save/senic_data/strs.h5ad Cluster")



