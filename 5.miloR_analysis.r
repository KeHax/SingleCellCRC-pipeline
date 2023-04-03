library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

#representative cluster: Myeloid cells
Myeloid_harmony <-readRDS("/data/khaoxian/project/Stagecrc/data/final_data/myeloid.corrected.rds")
Myeloid_Stage1_2 <- subset(Myeloid_harmony,subset = Stage %in% c("1","2"))
Myeloid_Stage2_3 <- subset(Myeloid_harmony,subset = Stage %in% c("2","3"))
Myeloid_Stage3_4 <- subset(Myeloid_harmony,subset = Stage %in% c("3","4"))
Myeloid_Stage1_4 <- subset(Myeloid_harmony,subset = Stage %in% c("1","4"))

#find k and d value
A <- c(15,20,25,30,35,40,45,50)
B <- c(15,20,25,30,35,40,45,50)
for (j in seq_along(A)) {
  a=A[j]
  for (i in seq_along(B)) {
    b=B[i]
    {
      Myeloid_milo <- as.SingleCellExperiment(Myeloid_harmony)
      Myeloid_milo <- Milo(Myeloid_milo)
      Myeloid_milo <- buildGraph(Myeloid_milo, k = a, d = b, reduced.dim = "HARMONY")
      Myeloid_milo <- makeNhoods(Myeloid_milo, prop = 0.1, k = a, d = b, refined = TRUE, reduced_dims = "HARMONY")
      p <- plotNhoodSizeHist(Myeloid_milo)
      if(mean(p[["data"]][["nh_size"]]) > 77 * 5){
        print(paste("Myeloid_milo",a,B[i],sep = " "))
      }
    }
  }
}


A <- c(15,20,25,30,35,40,45,50)
B <- c(15,20,25,30,35,40,45,50)
for (j in seq_along(A)) {
  a=A[j]
  for (i in seq_along(B)) {
    b=B[i]
    {
      Myeloid_Stage1_2_milo <- as.SingleCellExperiment(Myeloid_Stage1_2)
      Myeloid_Stage1_2_milo <- Milo(Myeloid_Stage1_2_milo)
      Myeloid_Stage1_2_milo <- buildGraph(Myeloid_Stage1_2_milo, k = a, d = b, reduced.dim = "HARMONY")
      Myeloid_Stage1_2_milo <- makeNhoods(Myeloid_Stage1_2_milo, prop = 0.1, k = a, d = b, refined = TRUE, reduced_dims = "HARMONY")
      p <- plotNhoodSizeHist(Myeloid_Stage1_2_milo)
      if(mean(p[["data"]][["nh_size"]]) > 35 * 5){
        print(paste("Myeloid_Stage1_2_milo",a,B[i],sep = " "))
      }
    }
  }
}


A <- c(15,20,25,30,35,40,45,50)
B <- c(15,20,25,30,35,40,45,50)
for (j in seq_along(A)) {
  a=A[j]
  for (i in seq_along(B)) {
    b=B[i]
    {
      Myeloid_Stage2_3_milo <- as.SingleCellExperiment(Myeloid_Stage2_3)
      Myeloid_Stage2_3_milo <- Milo(Myeloid_Stage2_3_milo)
      Myeloid_Stage2_3_milo <- buildGraph(Myeloid_Stage2_3_milo, k = a, d = b, reduced.dim = "HARMONY")
      Myeloid_Stage2_3_milo <- makeNhoods(Myeloid_Stage2_3_milo, prop = 0.1, k = a, d = b, refined = TRUE, reduced_dims = "HARMONY")
      p <- plotNhoodSizeHist(Myeloid_Stage2_3_milo)
      if(mean(p[["data"]][["nh_size"]]) > 50 * 5){
        print(paste("Myeloid_Stage2_3_milo",a,B[i],sep = " "))
      }
    }
  }
}


A <- c(15,20,25,30,35,40,45,50)
B <- c(15,20,25,30,35,40,45,50)
for (j in seq_along(A)) {
  a=A[j]
  for (i in seq_along(B)) {
    b=B[i]
    {
      Myeloid_Stage3_4_milo <- as.SingleCellExperiment(Myeloid_Stage3_4)
      Myeloid_Stage3_4_milo <- Milo(Myeloid_Stage3_4_milo)
      Myeloid_Stage3_4_milo <- buildGraph(Myeloid_Stage3_4_milo, k = a, d = b, reduced.dim = "HARMONY")
      Myeloid_Stage3_4_milo <- makeNhoods(Myeloid_Stage3_4_milo, prop = 0.1, k = a, d = b, refined = TRUE, reduced_dims = "HARMONY")
      p <- plotNhoodSizeHist(Myeloid_Stage3_4_milo)
      if(mean(p[["data"]][["nh_size"]]) > 42 * 5){
        print(paste("Myeloid_Stage3_4_milo",a,B[i],sep = " "))
      }
    }
  }
}


A <- c(15,20,25,30,35,40,45,50)
B <- c(15,20,25,30,35,40,45,50)
for (j in seq_along(A)) {
  a=A[j]
  for (i in seq_along(B)) {
    b=B[i]
    {
      Myeloid_Stage1_4_milo <- as.SingleCellExperiment(Myeloid_Stage1_4)
      Myeloid_Stage1_4_milo <- Milo(Myeloid_Stage1_4_milo)
      Myeloid_Stage1_4_milo <- buildGraph(Myeloid_Stage1_4_milo, k = a, d = b, reduced.dim = "HARMONY")
      Myeloid_Stage1_4_milo <- makeNhoods(Myeloid_Stage1_4_milo, prop = 0.1, k = a, d = b, refined = TRUE, reduced_dims = "HARMONY")
      p <- plotNhoodSizeHist(Myeloid_Stage1_4_milo)
      if(mean(p[["data"]][["nh_size"]]) > 27 * 5){
        print(paste("Myeloid_Stage1_4_milo",a,B[i],sep = " "))
      }
    }
  }
}

Myeloid_milo <- as.SingleCellExperiment(Myeloid_harmony)
Myeloid_milo <- Milo(Myeloid_milo)

#1. ConMyeloiduct KNN graph
Myeloid_milo <- buildGraph(Myeloid_milo, k = 30, d = 30, reduced.dim = "HARMONY")

#2. Defining representative neighbourhoods on the KNN graph
Myeloid_milo <- makeNhoods(Myeloid_milo, prop = 0.05, k = 30, d= 30, refined = TRUE, reduced_dims = "HARMONY")

#3. average number of neighborhood should be more than 5 x number of sample
p <- plotNhoodSizeHist(Myeloid_milo)

#4. Counting cells in neighbourhoods
Myeloid_milo <- countCells(Myeloid_milo, meta.data = as.data.frame(colData(Myeloid_milo)), sample="orig.ident")
#head(nhoodCounts(Myeloid_milo))

#5. Defining experimental design
Myeloid_milo_design <- data.frame(colData(Myeloid_milo))[,c("orig.ident", "Stage")]
Myeloid_milo_design$Stage <- factor(Myeloid_milo_design$Stage, levels = c("4", "3","2","1"))
Myeloid_milo_design <- distinct(Myeloid_milo_design)
rownames(Myeloid_milo_design) <- Myeloid_milo_design$orig.ident

#6. Computing neighbourhood connectivity
Myeloid_milo <- calcNhoodDistance(Myeloid_milo, d=30, reduced.dim = "HARMONY")
saveRDS(Myeloid_milo,file="/data/zdw/1_Integrate/CreateSeuratObject/Myeloid_milo.RDS")

#7. Testing
da_results <- testNhoods(Myeloid_milo, design = ~ Stage, design.df = Myeloid_milo_design)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

Myeloid_milo <- buildNhoodGraph(Myeloid_milo)

da_results <- annotateNhoods(Myeloid_milo, da_results, coldata_col = "cell_type")
da_results$cell_type <- ifelse(da_results$cell_type_fraction < 0.7, "Mixed", da_results$cell_type)

plotDAbeeswarms(da_results, group.by = "cell_type")



#stage 1 vs stage 2
Myeloid_Stage1_2_milo <- as.SingleCellExperiment(Myeloid_Stage1_2)
Myeloid_Stage1_2_milo <- Milo(Myeloid_Stage1_2_milo)

#1. Construct KNN graph
Myeloid_Stage1_2_milo <- buildGraph(Myeloid_Stage1_2_milo, k = 45, d = 50, reduced.dim = "HARMONY")

#2. Defining representative neighbourhoods on the KNN graph
Myeloid_Stage1_2_milo <- makeNhoods(Myeloid_Stage1_2_milo, prop = 0.1, k = 45, d = 50, refined = TRUE, reduced_dims = "HARMONY")

#3. average number of neighborhood should be more than 5 x number of sample
p <- plotNhoodSizeHist(Myeloid_Stage1_2_milo)

#4. Counting cells in neighbourhoods
Myeloid_Stage1_2_milo <- countCells(Myeloid_Stage1_2_milo, meta.data = as.data.frame(colData(Myeloid_Stage1_2_milo)), sample="orig.ident")

#5. Defining experimental design
Myeloid_Stage1_2_milo_design <- data.frame(colData(Myeloid_Stage1_2_milo))[,c("orig.ident", "Stage")]
Myeloid_Stage1_2_milo_design$Stage <- factor(Myeloid_Stage1_2_milo_design$Stage, levels = c("2", "1"))
Myeloid_Stage1_2_milo_design <- distinct(Myeloid_Stage1_2_milo_design)
rownames(Myeloid_Stage1_2_milo_design) <- Myeloid_Stage1_2_milo_design$orig.ident

#6. Computing neighbourhood connectivity
Myeloid_Stage1_2_milo <- calcNhoodDistance(Myeloid_Stage1_2_milo, d= 50, reduced.dim = "HARMONY")

#7. Testing
da_results <- testNhoods(Myeloid_Stage1_2_milo, design = ~ Stage, design.df = Myeloid_Stage1_2_milo_design)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

Myeloid_Stage1_2_milo <- buildNhoodGraph(Myeloid_Stage1_2_milo)

da_results <- annotateNhoods(Myeloid_Stage1_2_milo, da_results, coldata_col = "cell_type")

da_results$cell_type <- ifelse(da_results$cell_type_fraction < 0.7, "Mixed", da_results$cell_type)
p2 <- plotDAbeeswarms(da_results, group.by = "cell_type")
p1 <- p2 + labs(title="stage1_2", subtitle ="k = 45_d = 50")
print(p1)



#stage 2 vs stage 3
Myeloid_Stage2_3_milo <- as.SingleCellExperiment(Myeloid_Stage2_3)
Myeloid_Stage2_3_milo <- Milo(Myeloid_Stage2_3_milo)

#1. Construct KNN graph(找不到最适k，d值)
Myeloid_Stage2_3_milo <- buildGraph(Myeloid_Stage2_3_milo, k = 45, d = 50, reduced.dim = "HARMONY")

#2. Defining representative neighbourhoods on the KNN graph
Myeloid_Stage2_3_milo <- makeNhoods(Myeloid_Stage2_3_milo, prop = 0.1, k = 45, d= 50, refined = TRUE, reduced_dims = "HARMONY")

#3. average number of neighborhood should be more than 5 x number of sample
p <- plotNhoodSizeHist(Myeloid_Stage2_3_milo)
#4. Counting cells in neighbourhoods
Myeloid_Stage2_3_milo <- countCells(Myeloid_Stage2_3_milo, meta.data = as.data.frame(colData(Myeloid_Stage2_3_milo)), sample="orig.ident")

#5. Defining experimental design
Myeloid_Stage2_3_milo_design <- data.frame(colData(Myeloid_Stage2_3_milo))[,c("orig.ident", "Stage")]
Myeloid_Stage2_3_milo_design$Stage <- factor(Myeloid_Stage2_3_milo_design$Stage, levels = c("3", "2"))
Myeloid_Stage2_3_milo_design <- distinct(Myeloid_Stage2_3_milo_design)
rownames(Myeloid_Stage2_3_milo_design) <- Myeloid_Stage2_3_milo_design$orig.ident

#6. Computing neighbourhood connectivity
Myeloid_Stage2_3_milo <- calcNhoodDistance(Myeloid_Stage2_3_milo, d= 50, reduced.dim = "HARMONY")
saveRDS(Myeloid_Stage2_3_milo,file="/data/zdw/1_Integrate/CreateSeuratObject/Myeloid_Stage2_3_milo.RDS")

#7. Testing
da_results <- testNhoods(Myeloid_Stage2_3_milo, design = ~ Stage, design.df = Myeloid_Stage2_3_milo_design)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

Myeloid_Stage2_3_milo <- buildNhoodGraph(Myeloid_Stage2_3_milo)

da_results <- annotateNhoods(Myeloid_Stage2_3_milo, da_results, coldata_col = "cell_type")

da_results$cell_type <- ifelse(da_results$cell_type_fraction < 0.7, "Mixed", da_results$cell_type)
p2 <- plotDAbeeswarm(da_results, group.by = "cell_type")
p1 <- p2 + labs(title="stage2_3", subtitle ="-------k = 45_d = 50")
print(p1)



#stage3 vs stage 4
Myeloid_Stage3_4_milo <- as.SingleCellExperiment(Myeloid_Stage3_4)
Myeloid_Stage3_4_milo <- Milo(Myeloid_Stage3_4_milo)

#1. ConMyeloid_Stage3_4uct KNN graph(找不到最适k，d值)
Myeloid_Stage3_4_milo <- buildGraph(Myeloid_Stage3_4_milo, k = 45, d = 50, reduced.dim = "HARMONY")

#2. Defining representative neighbourhoods on the KNN graph
Myeloid_Stage3_4_milo <- makeNhoods(Myeloid_Stage3_4_milo, prop = 0.1, k = 45, d = 50, refined = TRUE, reduced_dims = "HARMONY")

#3. average number of neighborhood should be more than 5 x number of sample
p <- plotNhoodSizeHist(Myeloid_Stage3_4_milo)

#4. Counting cells in neighbourhoods
Myeloid_Stage3_4_milo <- countCells(Myeloid_Stage3_4_milo, meta.data = as.data.frame(colData(Myeloid_Stage3_4_milo)), sample="orig.ident")

#5. Defining experimental design
Myeloid_Stage3_4_milo_design <- data.frame(colData(Myeloid_Stage3_4_milo))[,c("orig.ident", "Stage")]
Myeloid_Stage3_4_milo_design$Stage <- factor(Myeloid_Stage3_4_milo_design$Stage, levels = c("4", "3"))
Myeloid_Stage3_4_milo_design <- distinct(Myeloid_Stage3_4_milo_design)
rownames(Myeloid_Stage3_4_milo_design) <- Myeloid_Stage3_4_milo_design$orig.ident

#6. Computing neighbourhood connectivity
Myeloid_Stage3_4_milo <- calcNhoodDistance(Myeloid_Stage3_4_milo, d= 50, reduced.dim = "HARMONY")
saveRDS(Myeloid_Stage3_4_milo,file="/data/zdw/1_Integrate/CreateSeuratObject/Myeloid_Stage3_4_milo.RDS")

#7. Testing
da_results <- testNhoods(Myeloid_Stage3_4_milo, design = ~ Stage, design.df = Myeloid_Stage3_4_milo_design)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

Myeloid_Stage3_4_milo <- buildNhoodGraph(Myeloid_Stage3_4_milo)

da_results <- annotateNhoods(Myeloid_Stage3_4_milo, da_results, coldata_col = "cell_type")

da_results$cell_type <- ifelse(da_results$cell_type_fraction < 0.7, "Mixed", da_results$cell_type)
pdata <- plotDAbeeswarms(da_results, group.by = "cell_type", return.data=T)

p2 <- plotDAbeeswarms(da_results, group.by = "cell_type")
p1 <- p2 + labs(title="stage3_4", subtitle ="-----k = 45_d = 50")
print(p1)
