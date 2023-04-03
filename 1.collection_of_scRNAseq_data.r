library(data.table)
library(Seurat)
library(dplyr)
library(tidyverse)
library(magrittr)
library(reshape2)
library(readxl)
library(readr)
#Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer. Nat Genet. 2020----------------
#GSE132465-------------------
#1. read data
GSE132465_data <- fread("/data/zdw/GSE132465/Data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
data = data.frame(GSE132465_data,check.names = F, row.names = 1)
data1 <- as(as.matrix(data),"dgCMatrix")
GSE132465_seu_data  <- CreateSeuratObject(counts = data1, min.cells = 3)

#2. select tumor samples
GSE132465_seu_data_T <- subset(x = GSE132465_seu_data, idents = c("SMC01-T","SMC02-T","SMC03-T","SMC04-T","SMC05-T","SMC06-T","SMC07-T","SMC08-T","SMC09-T","SMC10-T",
                                                                  "SMC11-T","SMC14-T","SMC15-T","SMC16-T","SMC17-T","SMC18-T","SMC19-T","SMC20-T","SMC21-T",
                                                                  "SMC22-T","SMC23-T","SMC24-T","SMC25-T"))
GSE132465_seu_data_T@meta.data$orig.ident <- GSE132465_seu_data_T@active.ident

#3. save data
saveRDS(GSE132465_seu_data_T, file = "/data/zdw/GSE132465/CreateSeuratObject/GSE132465_seu_data_T.RDS")

#GSE144735-------------------
#1. read data
# 读取txt格式的矩阵
GSE144735_data <- fread("/data/zdw/GSE144735/Data/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz")
data = data.frame(GSE144735_data,check.names = F, row.names = 1)
data1 <- as(as.matrix(data),"dgCMatrix")
GSE144735_seu_data  <- CreateSeuratObject(counts = data1, min.cells = 3)

#2. select tumor samples
GSE144735_T <- subset(x = GSE144735_T, idents = c("KUL01-T","KUL19-T","KUL21-T","KUL28-T","KUL30-T","KUL31-T"))
GSE144735_T@meta.data$orig.ident <- GSE144735_T@active.ident

#3. save data
saveRDS(GSE144735_T, file = "/data/zdw/GSE144735/CreateSeuratObject/GSE144735_T.RDS")

#Immune phenotypic linkage between colorectal cancer and liver metastasis. Cancer Cell. 2022--------------------
#GSE146771-------------------
#1. read data
GSE146771_data <- fread("/data/zdw/GSE146771/Data/CRC.Leukocyte.10x.TPM.txt")
GSE146771_Metadata <- fread("/data/zdw/GSE146771/Data/CRC.Leukocyte.10x.Metadata.txt")
data = data.frame(GSE146771_data,check.names = F, row.names = 1)
data1 <- as(as.matrix(data),"dgCMatrix")
GSE146771_seu_data  <- CreateSeuratObject(counts = data1, min.cells = 3)

#2. select tumor samples
GSE146771_seu_data@meta.data$patient <- sapply(rownames(GSE146771_seu_data@meta.data),function(x){
  str_split(x,"_")[[1]][3]
})
GSE146771_seu_data@meta.data$orig.ident <- GSE146771_seu_data@meta.data$patient
GSE146771_seu_data@meta.data$patient=NULL

#3. save data
saveRDS(GSE146771_seu_data, file = "/data/zdw/GSE146771/CreateSeuratObject/GSE146771_seu_data.RDS")

#Single-cell transcriptomic proﬁling unravels the adenoma-initiation role of protein tyrosine kinases during colorectal tumorigenesis-----------------
#GSE161277-------------------
#1. read data
samples <- read_excel("/data/zdw/GSE161277/Data/metadata-GSE161277.xlsx")%>% .$Sample

for (i in seq_along(samples)){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0("/data/zdw/GSE161277/Data/", samples[i], "/filtered_feature_bc_matrix")))
}

#create seurat objects from cellranger files
for (i in seq_along(samples)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i], min.cells = 3))
}
rm(c)

GSE161277_seu_data <- merge(seu_obj1, y = c(seu_obj2, seu_obj3, seu_obj4), add.cell.ids = samples)

#2. save data
saveRDS(GSE161277_seu_data,"/data/zdw/GSE161277/CreateSeuratObject/GSE161277_seu_data.RDS")

#Immune phenotypic linkage between colorectal cancer and liver metastasis-----------
#GSE164522-------------------
#1. read data and select tumor samples
GSE164522_CRLM_PT_expression_data <- fread("/data/zdw/GSE164522/Data/CRLM_PT_expression.csv")
GSE164522_CRLM_metadata <- fread("/data/zdw/GSE164522/Data/GSE164522_CRLM_metadata.csv")
GSE164522_CRLM_metadata_PT <- subset(x= GSE164522_CRLM_metadata,subset = tissue == "primary tumor")
data = data.frame(GSE164522_CRLM_PT_expression_data,check.names = F, row.names = 1)
data1 <- as(as.matrix(data),"dgCMatrix")
GSE164522_CRLM_PT_expression_seu_data  <- CreateSeuratObject(counts = data1, min.cells = 3)

#2. save data
metadata <- FetchData(GSE164522_CRLM_PT_expression_seu_data, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$V1 <- metadata$cell_id
metadata <- left_join(x = metadata, y = GSE164522_CRLM_metadata_PT, by = "V1")
rownames(metadata) <- metadata$V1
GSE164522_CRLM_PT_expression_seu_data <- AddMetaData(GSE164522_CRLM_PT_expression_seu_data, metadata = metadata)
GSE164522_CRLM_PT_expression_seu_data@meta.data$orig.ident <- GSE164522_CRLM_PT_expression_seu_data@meta.data$patient
GSE164522_CRLM_PT_expression_seu_data@meta.data$celltype_sub= NULL
saveRDS(GSE164522_CRLM_PT_expression_seu_data, file = "/data/zdw/GSE164522/CreateSeuratObject/GSE164522_CRLM_PT_seu_data.RDS")

# A single-cell atlas of liver metastases of colorectal cancer reveals reprogramming of the tumor microenvironment in response to preoperative chemotherapy------
#GSE178318-------------------
#1. read data
A <- read_excel("/data/zdw/GSE178318/Metadata/metadata-178318.xlsx",range = cell_cols("B:B")) %>% .$Sample
patients_metadata <- read_excel("/data/zdw/GSE178318/Metadata/metadata-178318.xlsx")

data.data <- Read10X(data.dir ="/data/zdw/GSE178318/data")
GSE178318_seu_data <- CreateSeuratObject(counts = data.data, min.cells = 3)

#2. select tumor samples
GSE178318_seu_data@meta.data$patient <- sapply(rownames(GSE178318_seu_data@meta.data),function(x){
  str_split(x,"_")[[1]][2]
})

GSE178318_seu_data@meta.data$sample <- sapply(rownames(GSE178318_seu_data@meta.data),function(x){
  str_split(x,"_")[[1]][3]
})

GSE178318_seu_data$sample_id <- paste(GSE178318_seu_data@meta.data$patient, GSE178318_seu_data@meta.data$sample, sep = "_")
GSE178318_seu_data$orig.ident <- GSE178318_seu_data$sample_id
GSE178318_seu_data@meta.data$patient=NULL
GSE178318_seu_data$sample_id=NULL
GSE178318_seu_data@meta.data$samples <- GSE178318_seu_data@meta.data$orig.ident

Idents(GSE178318_seu_data) <- GSE178318_seu_data$samples
new_ids_main <- c("COL07_CRC","COL07_LM" , "COL12_CRC",   "COL12_LM", "COL12_PBMC",  "COL15_CRC",   "COL15_LM",  "COL16_CRC",   
                  "COL16_LM",  "COL17_CRC",   "COL17_LM", "COL17_PBMC",  "COL18_CRC", "COL18_LM", "COL18_PBMC")
names(new_ids_main) <- levels(GSE178318_seu_data)
GSE178318_seu_data <- RenameIdents(GSE178318_seu_data, new_ids_main)
GSE178318_seu_data@meta.data$main_cell_type <- Idents(GSE178318_seu_data)
GSE178318_seu_data@meta.data$samples = NULL

GSE178318_seu_data_T <- subset(x = GSE178318_seu_data, idents = c("COL07_CRC","COL12_CRC","COL15_CRC","COL16_CRC","COL17_CRC","COL18_CRC"))
GSE178318_seu_data_T@meta.data$orig.ident <- GSE178318_seu_data_T@active.ident

#3. save data
saveRDS(GSE178318_seu_data_T, file = "/data/zdw/GSE178318/CreateSeuratObject/GSE178318_seu_data_T.RDS")

#Resolving the difference between leftsided and right-sided colorectal cancer by single-cell sequencing------
#GSE188711-------------------
#1. read data
Patient_metadata <- read_excel("/data/zdw/GSE188711/Metadata/metadata-GSE188711.xlsx")
samples <- Patient_metadata %>% .$Sample


for (i in seq_along(samples)){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0("/data/zdw/GSE188711/Data/", samples[i], "/filtered_feature_bc_matrix")))
}

### create seurat objects from cellranger files
for (i in seq_along(samples)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i], min.cells = 3))
}


seu_obj <- merge(seu_obj1, y = c(seu_obj2, seu_obj3, seu_obj4, seu_obj5, seu_obj6), add.cell.ids = samples)

#2. select tumor samples
GSE188711_seu_data@meta.data$samples <- GSE188711_seu_data@meta.data$orig.ident

Idents(GSE188711_seu_data) <- GSE188711_seu_data@meta.data$samples
new_ids_main <- c("Left_sided_CRC_Patient1_T", "Left_sided_CRC_Patient2_T", "Left_sided_CRC_Patient3_T", "Right_sided_CRC_Patient1_T", "Right_sided_CRC_Patient2_T", "Right_sided_CRC_Patient3_T")
names(new_ids_main) <- levels(GSE188711_seu_data)
GSE188711_seu_data <- RenameIdents(GSE188711_seu_data, new_ids_main)
GSE188711_seu_data@meta.data$orig.ident <- Idents(GSE188711_seu_data)
GSE188711_seu_data@meta.data$samples =NULL

#3. save data
saveRDS(GSE188711_seu_data,"/data/zdw/GSE188711/CreateSeuratObject/GSE188711_seu_data.RDS")

#Refining colorectal cancer classification and clinical stratification through a single-cell atlas------------------
#GSE200997-------------------
#1. read data
GSE200997_GEO_processed_CRC_10X_raw_UMI_count_data <- fread("/data/zdw/GSE200997/Data/GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv")
GSE200997_GEO_processed_CRC_10X_cell_annotation <- fread("/data/zdw/GSE200997/Data/GSE200997_GEO_processed_CRC_10X_cell_annotation.csv")
data = data.frame(GSE200997_GEO_processed_CRC_10X_raw_UMI_count_data,check.names = F, row.names = 1)
data1 <- as(as.matrix(data),"dgCMatrix")

GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data  <- CreateSeuratObject(counts = data1, min.cells = 3)

#2. select tumor samples
GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data$patient <- sapply(rownames(GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data),function(x){
  str_split(x,"_")[[1]][2]
})
GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data$patient <- toupper(GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data$patient)
GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data$sample <- sapply(rownames(GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data),function(x){
  str_split(x,"_")[[1]][1]
})

GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data$sample_id <- paste(GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data$sample, GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data$patient, sep = "_")
GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data$orig.ident <- GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data$sample_id

GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data$patient=NULL
GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data$sample_id=NULL
GSE200997_GEO_processed_CRC_10X_raw_UMI_count_seu_data@meta.data$sample=NULL

GSE200997_seu_data@meta.data$orig.ident <- as.factor(GSE200997_seu_data@meta.data$orig.ident)

GSE200997_seu_data@meta.data$samples <- GSE200997_seu_data@meta.data$orig.ident

Idents(GSE200997_seu_data) <- GSE200997_seu_data$samples
new_ids_main <- c("B_CAC10", "B_CAC11", "B_CAC14", "B_CAC15",  "B_CAC4",  "B_CAC6",  "B_CAC7",  "T_CAC1", "T_CAC10", 
                  "T_CAC11", "T_CAC12", "T_CAC13", "T_CAC14", "T_CAC15", "T_CAC16",  "T_CAC2",  "T_CAC3",  "T_CAC4",  
                  "T_CAC5", "T_CAC6",  "T_CAC7",  "T_CAC8",  "T_CAC9")
names(new_ids_main) <- levels(GSE200997_seu_data)
GSE200997_seu_data <- RenameIdents(GSE200997_seu_data, new_ids_main)
GSE200997_seu_data@meta.data$samples = NULL


GSE200997_seu_data_T <- subset(x = GSE200997_seu_data, idents = c("T_CAC1", "T_CAC10", "T_CAC11", "T_CAC12", "T_CAC13", "T_CAC14", "T_CAC15", "T_CAC16",  "T_CAC2",  "T_CAC3",  "T_CAC4",  
                                                                  "T_CAC5", "T_CAC6",  "T_CAC7",  "T_CAC8",  "T_CAC9"))
GSE200997_seu_data_T@meta.data$orig.ident <- GSE200997_seu_data_T@active.ident

#3. save data
saveRDS(GSE200997_seu_data_T, file = "/data/zdw/GSE200997/CreateSeuratObject/GSE200997_seu_data_T.RDS")


