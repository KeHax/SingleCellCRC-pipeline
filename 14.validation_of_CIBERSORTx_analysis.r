#consensus between results of CIBERSORTx and consensusTMEAnalysis
library(ConsensusTME)
library(tidyverse)
library(ggplot2)
library(patchwork)

#1. read cohort data------------
cohort.DEList <- readRDS("/data/khaoxian/project/Stagecrc/data/rds/CRCcohorts_RNAseq.rds")

#transformation: CPM from counts
a=cohort.DEList$COAD$counts
a=t(t(a)*1000000/c(colSums(a)))
a=log(a+1,base=2)
cohort.DEList$COAD$counts <- as.data.frame(a)

a=cohort.DEList$READ$counts
a=t(t(a)*1000000/c(colSums(a)))
a=log(a+1,base=2)
cohort.DEList$READ$counts <- as.data.frame(a)

a=cohort.DEList$GSE33113$counts
a=log(a+1,base=2)
cohort.DEList$GSE33113$counts <- as.data.frame(a)

#using COAD for validation
#input CIBERSORTx analysis result
coad.ciber <- read.table("/data/khaoxian/project/Stagecrc/data/cibersortx/webout/coad.cibersortx.result.txt", sep = "\t",
                         header = T,row.names = 1)
coad.ciber <- as.data.frame(t(coad.ciber))
cluster <- rownames(coad.ciber)
cohort.DEList$COAD$counts <- rbind(cohort.DEList$COAD$counts, coad.ciber)

consensus.DEList <- cohort.DEList

#2. do ConsensusTME analysis------------
cohort.DEList <- purrr::map(cohort.DEList, function(x){
  x$counts <- ConsensusTME::consensusTMEAnalysis(as.matrix(x$counts), cancer = "COAD", statMethod = "gsva")
  return(x)
})
for(i in 1:length(consensus.DEList)){
  consensus.DEList[[i]]$counts <- rbind(consensus.DEList[[i]]$counts, cohort.DEList[[i]]$counts)
}

#3. analyze relationship between results of CIBERSORTx and ConsensusTME------------
cluster <- tail(rownames(consensus.DEList$COAD$counts),100)[26:78]
pcluster <- tail(rownames(consensus.DEList$COAD$counts),100)[82:100]

#integrate B DC CTL Macro MAST NK CD4 CD8 TREG ENdo Fibro Mono Plasma subtypes
bcell.list <- list()
bcell.list[[1]] <- cluster[14:17]
bcell.list[[2]] <- "B_cells"

DC.list <- list()
DC.list[[1]] <- cluster[c(25,29)]
DC.list[[2]] <- "Dendritic_cells"

Macro.list <- list()
Macro.list[[1]] <- cluster[c(20,21,26,27,28)]
Macro.list[[2]] <- "Macrophages"

mono.list <- list()
mono.list[[1]] <- cluster[c(22,24)]
mono.list[[2]] <- "Monocytes"

mast.list <- list()
mast.list[[1]] <- cluster[c(30)]
mast.list[[2]] <- "Mast_cells"

mast.list <- list()
mast.list[[1]] <- cluster[c(30)]
mast.list[[2]] <- "Mast_cells"

nk.list <- list()
nk.list[[1]] <- cluster[c(51,53)] 
nk.list[[2]] <- "NK_cells"

cd4.list <- list()
cd4.list[[1]] <- cluster[c(31,33,34,35,37)] 
cd4.list[[2]] <- "T_cells_CD4"

cd8.list <- list()
cd8.list[[1]] <- cluster[c(45,46,40,47,43)] 
cd8.list[[2]] <- "T_cells_CD8" 

ctl.list <- list()
ctl.list[[1]] <- cluster[c(38)]
ctl.list[[2]] <- "Cytotoxic_cells"

treg.list <- list()
treg.list[[1]] <- cluster[c(32)]
treg.list[[2]] <- "T_regulatory_cells"

endo.list <- list()
endo.list[[1]] <- cluster[c(7)]
endo.list[[2]] <- "Endothelial"

fibro.list <- list()
fibro.list[[1]] <- cluster[c(8:13)]
fibro.list[[2]] <- "Fibroblasts"

plasma.list <- list()
plasma.list[[1]] <- cluster[c(18,19)]
plasma.list[[2]] <- "Plasma_cells"

imm.list <- list()
imm.list[[1]] <- cluster[c(14:53)] 
imm.list[[2]] <- "Immune_Score"

cluster
pcluster

Formconsensus <- function(consensus.DEList,
                          cell.pair){
  genelist <- purrr::map(consensus.DEList, list("counts"))
  genelist1 <- lapply(genelist, function(x){
    x <- x[cell.pair[[1]],]
    x <- colSums(x)
    return(x)
  })
  genelist2 <- lapply(genelist, function(x){
    x <- x[cell.pair[[2]],]
    x <- colSums(x)
    return(x)
  })
  genelist <- purrr::map(1:length(genelist2), function(x){
    x <- data.frame(Cibersortx=genelist1[[x]],
                    ConsensusTME=genelist2[[x]])
    return(x)
  })
  names(genelist) <- names(genelist2)
  return(genelist)
}

bcell.consen <- Formconsensus(consensus.DEList,cell.pair = bcell.list)
dc.consen <- Formconsensus(consensus.DEList,cell.pair = DC.list)
mac.consen <- Formconsensus(consensus.DEList,cell.pair = Macro.list)
mono.consen <- Formconsensus(consensus.DEList,cell.pair = mono.list)
mast.consen <- Formconsensus(consensus.DEList,cell.pair = mast.list)
nk.consen <- Formconsensus(consensus.DEList,cell.pair = nk.list)
cd4.consen <- Formconsensus(consensus.DEList,cell.pair = cd4.list)
cd8.consen <- Formconsensus(consensus.DEList,cell.pair = cd8.list) 
ctl.consen <- Formconsensus(consensus.DEList,cell.pair = ctl.list)
treg.consen <- Formconsensus(consensus.DEList,cell.pair = treg.list)##
endo.consen <- Formconsensus(consensus.DEList,cell.pair = endo.list)
fibro.consen <- Formconsensus(consensus.DEList,cell.pair = fibro.list)
plasma.consen <- Formconsensus(consensus.DEList,cell.pair = plasma.list)
imm.consen <- Formconsensus(consensus.DEList,cell.pair = imm.list)

consen.list <- list()
consen.list[["B cell"]] <- bcell.consen
consen.list[["DC"]] <- dc.consen
consen.list[["Macrophage"]] <- mac.consen
consen.list[["Monocyte"]] <- mono.consen
consen.list[["Mast cell"]] <- mast.consen
consen.list[["NK cell"]] <- nk.consen
consen.list[["CD4 cell"]] <- cd4.consen
consen.list[["CD8 cell"]] <- cd8.consen
consen.list[["Cytotoxic cell"]] <- ctl.consen
consen.list[["Treg"]] <- treg.consen
consen.list[["Endothelium"]] <- endo.consen
consen.list[["Fibroblast"]] <- fibro.consen
consen.list[["Plasma cell"]] <- plasma.consen
consen.list[["Immune cell"]] <- imm.consen

#4. plotting------------
plot.list <- list()
for (i in 1:length(consen.list)) {
  cluster_name <- names(consen.list)[i]
  plot.list[[i]] <- ggplot(consen.list[[i]][[1]], 
                           aes_string(x = "Cibersortx", 
                                      y = "ConsensusTME")) + 
    geom_point(size=0.1, color="orange") 
    geom_smooth(method = 'lm', formula = y ~ x, color="red", lwd=0.3) + 
    ggpmisc::stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label..,sep = '~~~~')), 
                          formula = y ~ x, parse = T, size=1.9)+
    labs(x=paste("CIBERSORTx:",cluster_name,sep = " "),
         y=paste("ConsensusTME:",cluster_name,sep = " "))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title = element_text(size=6),
          axis.text = element_text(size=5),
          axis.line = element_blank(),
          axis.ticks = element_line(size=0.25),
          panel.border=element_rect(fill = "transparent",size = 0.2),
          plot.title = element_text(hjust = 0.5,size = 7))
}
p <- cowplot::plot_grid(plotlist = plot.list, ncol = 4)
ggsave("/data/khaoxian/project/Stagecrc/output/deconvolution/consenesus.pearson.pdf", p ,
       height = 21,width = 21, units = "cm")

plot.list <- list()
for (i in 1:length(consen.list)) {
  cluster_name <- names(consen.list)[i]
  plot.list[[i]] <- ggplot(consen.list[[i]][[1]], 
                           aes_string(x = "Cibersortx", 
                                      y = "ConsensusTME")) + 
    geom_point(size=0.1, color="orange") + 
    geom_smooth(method = 'lm', formula = y ~ x, color="red", lwd=0.3) + 
    stat_cor(method = "spearman", label.x = 0,label.y = max(consen.list[[i]][[1]]$ConsensusTME)*1.5)+
    labs(x=paste("CIBERSORTx:",cluster_name,sep = " "),
         y=paste("ConsensusTME:",cluster_name,sep = " "))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title = element_text(size=6),
          axis.text = element_text(size=5),
          axis.line = element_blank(),
          axis.ticks = element_line(size=0.25),
          panel.border=element_rect(fill = "transparent",size = 0.2),
          plot.title = element_text(hjust = 0.5,size = 7))
}

p <- cowplot::plot_grid(plotlist = plot.list, ncol = 4)
ggsave("/data/khaoxian/project/Stagecrc/output/deconvolution/consenesus.pearson.R.extracted.pdf", p ,
       height = 21,width = 21, units = "cm")









