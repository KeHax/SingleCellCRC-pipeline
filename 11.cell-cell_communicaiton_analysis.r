library(CellChat)
library(Seurat)
library(ggplot2)                  
library(patchwork)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(tidyverse)
#1. make Cellchat object-----
total_object <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/total_object.RDS")
for (i in 1:4) {
  
  total_stage1 <- subset(total_object,subset = Stage %in% as.character(i))
  
  #Part I: Data input & processing and initialization of CellChat object
  #Load data
  data.input = total_stage1@assays[["RNA"]]@data
  meta = total_stage1@meta.data
  unique(meta$cell_type)
  
  #Create a CellChat object
  cellchat_stage1 <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
  
  
  #Add cell information into meta slot of the object (Optional)
  cellchat_stage1 <- addMeta(cellchat_stage1, meta = meta)
  cellchat_stage1 <- setIdent(cellchat_stage1, ident.use = "cell_type")
  groupSize <- as.numeric(table(cellchat_stage1@idents))
  
  #Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat_stage1@DB <- CellChatDB.use
  
  #Preprocessing the expression data for cell-cell communication analysis
  cellchat_stage1 <- subsetData(cellchat_stage1)
  future::plan("multiprocess", workers = 4)
  cellchat_stage1 <- identifyOverExpressedGenes(cellchat_stage1)
  cellchat_stage1 <- identifyOverExpressedInteractions(cellchat_stage1)
  cellchat_stage1 <- projectData(cellchat_stage1, PPI.human)
  
  #Part II: Inference of cell-cell communication network
  #Compute the communication probability and infer cellular communication network
  cellchat_stage1 <- computeCommunProb(cellchat_stage1,raw.use = FALSE,population.size = TRUE)
  cellchat_stage1 <- filterCommunication(cellchat_stage1, min.cells = 10)
  
  #Extract the inferred cellular communication network as a data frame
  #subsetCommunication：access the inferred cell-cell communications of interest
  df.net <- subsetCommunication(cellchat_stage1)
  df.net.1 <- subsetCommunication(cellchat_stage1, sources.use = c(1,2), targets.use = c(4,5))
  df.net.2 <- subsetCommunication(cellchat_stage1, signaling = c("WNT", "TGFb"))
  
  #Infer the cell-cell communication at a signaling pathway level
  cellchat_stage1 <- computeCommunProbPathway(cellchat_stage1)
  
  #Calculate the aggregated cell-cell communication network
  cellchat_stage1 <- aggregateNet(cellchat_stage1)
  
  groupSize <- as.numeric(table(cellchat_stage1@idents))
  
  #Part III: Visualization of cell-cell communication network
  #Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
  
  #Part IV: Systems analysis of cell-cell communication network
  #Compute and visualize the network centrality scores
  # Compute the network centrality scores
  cellchat_stage1 <- netAnalysis_computeCentrality(cellchat_stage1, slot.name = "netP") 
  selectK(cellchat_stage1, pattern = "outgoing")
  
  nPatterns = 4
  cellchat_stage1 <- identifyCommunicationPatterns(cellchat_stage1, pattern = "outgoing", k = nPatterns)
  
  
  #Identify and visualize incoming communication pattern of target cells
  selectK(cellchat_stage1, pattern = "incoming")
  
  nPatterns = 4
  cellchat_stage1 <- identifyCommunicationPatterns(cellchat_stage1, pattern = "incoming", k = nPatterns)
  
  
  #Manifold and classification learning analysis of signaling networks
  #Identify signaling groups based on their functional similarity
  cellchat_stage1 <- computeNetSimilarity(cellchat_stage1, type = "functional")
  cellchat_stage1 <- netEmbedding(cellchat_stage1, type = "functional")
  cellchat_stage1 <- netClustering(cellchat_stage1, type = "functional")
  netVisual_embedding(cellchat_stage1, type = "functional", label.size = 3.5)
  
  #Identify signaling groups based on structure similarity
  cellchat_stage1 <- computeNetSimilarity(cellchat_stage1, type = "structural")
  cellchat_stage1 <- netEmbedding(cellchat_stage1, type = "structural")
  cellchat_stage1 <- netClustering(cellchat_stage1, type = "structural")
  netVisual_embedding(cellchat_stage1, type = "structural", label.size = 3.5)
  
  #Save the cellchat_stage1 object
  saveRDS(cellchat_stage1, file = paste("/data/zdw/1_Integrate/CreateSeuratObject/cellchat_stage",i,".RDS", sep = ""))
}

#2. DE analysis-----
cellchat_stage1 <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/cellchat_stage1.RDS")
cellchat_stage2 <- readRDS("/data/zdw/1_Integrate/CreateSeuratObject/cellchat_stage2.RDS")

group.new = levels(cellchat_stage1@idents)
cellchat_stage2 <- liftCellChat(cellchat_stage2, group.new)

object.list <- list(stage1 = cellchat_stage1, stage2 = cellchat_stage2)
cellchat <- mergeCellChat(object.list, 
                          add.names = names(object.list), 
                          cell.prefix = FALSE,
                          merge.data=TRUE)

#Part I: Predict general principles of cell-cell communication
#Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#Part II: Identify the conserved and context-specific signaling pathways
#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)

#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")

#Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
# dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
gg1 + gg2

#Identify dysfunctional signaling by using differential expression analysis
pos.dataset = "stage1"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "stage1",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "stage2",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


# visualize the enriched ligands in the first condition
net.down_stage1_2 <- net.down
computeEnrichmentScore(net.down_stage1_2, species = 'human')

# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human')
net.up <- subset(net.up,net.up$receptor.pct.2 != "NA")
net.up_stage1_2 <- net.up
computeEnrichmentScore(net.up_stage1_2, species = 'human')

#Save the merged CellChat object
cellchat_stage1_2 <- cellchat
saveRDS(cellchat_stage1_2, file = "/data/zdw/1_Integrate/CreateSeuratObject/cellchat_stage1_2.RDS")

#3. DE analysis visualization----
cellchat_stage1 <- readRDS("/data/khaoxian/project/Stagecrc/data/lps_cellchat/cellchat_stage1.RDS")
cellchat_stage2 <- readRDS("/data/khaoxian/project/Stagecrc/data/lps_cellchat/cellchat_stage2.RDS")
cellchat_stage3 <- readRDS("/data/khaoxian/project/Stagecrc/data/lps_cellchat/cellchat_stage3.RDS")
cellchat_stage4 <- readRDS("/data/khaoxian/project/Stagecrc/data/lps_cellchat/cellchat_stage4.RDS")

cc <- levels(cellchat_stage1@idents)
ccs <- c(8:33)
ccs <- c(8:14, 28:33)

group.cellType <- cc[ccs]

object.list <- list(
  stage1=cellchat_stage1,
  stage2=cellchat_stage2,
  stage3=cellchat_stage3,
  stage4=cellchat_stage4
)
object.list2 <- list(
  stage1=cellchat_stage1,
  stage2=cellchat_stage2,
  stage3=cellchat_stage3,
  stage4=cellchat_stage4
)
object.list <- lapply(object.list, function(x) {
  subsetCellChat(
    x,
    cells.use = NULL,
    idents.use = cc[ccs],
    group.by = NULL,
    invert = FALSE,
    thresh = 0.05
  )})

pathways.show <- c("LIFR")  
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf("/data/khaoxian/project/Stagecrc/output/LR/lifr.stage1-4.15cm_updated.pdf",height = 15, width=15)
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      idents.use = cc[ccs])
}
dev.off()

#4. pathway analysis of enriched ligands and receptors----
load("/data/khaoxian/project/Stagecrc/data/lps_cellchat/net.down_stage1_2.Rdata")
computeEnrichmentScore(net.down_stage1_2, species = 'human')

load("/data/khaoxian/project/Stagecrc/data/lps_cellchat/net.up_stage1_2.Rdata")
computeEnrichmentScore(net.up_stage1_2, species = 'human')

load("/data/khaoxian/project/Stagecrc/data/lps_cellchat/net.down_stage2_3.Rdata")
computeEnrichmentScore(net.down_stage2_3, species = 'human')

load("/data/khaoxian/project/Stagecrc/data/lps_cellchat/net.up_stage2_3.Rdata")
computeEnrichmentScore(net.up_stage2_3, species = 'human')

load("/data/khaoxian/project/Stagecrc/data/lps_cellchat/net.down_stage3_4.Rdata")
computeEnrichmentScore(net.down_stage3_4, species = 'human')

load("/data/khaoxian/project/Stagecrc/data/lps_cellchat/net.up_stage3_4.Rdata")
computeEnrichmentScore(net.up_stage3_4, species = 'human')


#up is stage1 ebriched
#extract LR and run pathways
table(net.up_stage1_2$pathway_name)
#stage1
stage1.gene <- sapply(net.up_stage1_2$interaction_name, function(x){str_split(x,"_")})
stage1.gene <- unlist(stage1.gene)
stage1.gene <- unique(stage1.gene)

#stage2
stage2vs1.gene <- sapply(net.down_stage1_2$interaction_name, function(x){str_split(x,"_")})
stage2vs1.gene <- unlist(stage2vs1.gene)
stage2vs1.gene <- unique(stage2vs1.gene)

stage2vs3.gene <- sapply(net.up_stage2_3$interaction_name, function(x){str_split(x,"_")})
stage2vs3.gene <- unlist(stage2vs3.gene)
stage2vs3.gene <- unique(stage2vs3.gene)

stage2.gene <- intersect(stage2vs1.gene,stage2vs3.gene)

#stage3
stage3vs2.gene <- sapply(net.down_stage2_3$interaction_name, function(x){str_split(x,"_")})
stage3vs2.gene <- unlist(stage3vs2.gene)
stage3vs2.gene <- unique(stage3vs2.gene)

stage3vs4.gene <- sapply(net.up_stage3_4$interaction_name, function(x){str_split(x,"_")})
stage3vs4.gene <- unlist(stage3vs4.gene)
stage3vs4.gene <- unique(stage3vs4.gene)

stage3.gene <- intersect(stage3vs2.gene,stage2vs3.gene)

#stage4
stage4.gene <- sapply(net.down_stage3_4$interaction_name, function(x){str_split(x,"_")})
stage4.gene <- unlist(stage4.gene)
stage4.gene <- unique(stage4.gene)

#GO analysis
symbols_list <- list(
  Stage1=stage1.gene,
  Stage2=stage2.gene,
  Stage3=stage3.gene,
  Stage4=stage4.gene
)

gcSample = lapply(symbols_list, function(y){ 
  y=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                               keys = y,
                                               columns = 'ENTREZID',
                                               keytype = 'SYMBOL')[,2]))
  y
})

formula_res <- compareCluster(
  gcSample, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

lineage1_ego <- clusterProfiler::simplify(
  formula_res, 
  cutoff=0.5, 
  by="p.adjust", 
  select_fun=min
) 

dotplot(lineage1_ego, showCategory=15)  + scale_y_discrete(labels=function(x) str_wrap(x, width=50)) +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.position = "right",
        axis.title = element_text(size=12))

ggsave("/data/khaoxian/project/Stagecrc/output/pathway/stage4_cellchat.de.go_bp.pdf",width = 17.8,height = 22, units="cm")

saveRDS(lineage1_ego, "/data/khaoxian/project/Stagecrc/output/pathway/stage1-4.LR.pathway.enrichgo.simplified.rds")

