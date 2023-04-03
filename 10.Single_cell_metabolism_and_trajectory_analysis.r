#CD8T: metabolic trajectory
library(parallel)
library(Vision)
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(gam)
options(mc.cores = 48)
#use magic imputated data
cd8t.palantir <- anndata::read_h5ad("/data/khaoxian/project/Stagecrc/output/palantir/cd8t/module1/Crc_cd8t.module1.h5ad")
cd8t.lay=cd8t.palantir$layers[["MAGIC_imputed_data"]]
cd8t.lay <- t(cd8t.lay)

#1. run scmetabolism
signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                     package = "scMetabolism")
gmtFile <- signatures_KEGG_metab
vis <- Vision(cd8t.lay, signatures = gmtFile)
vis <- analyze(vis)
signature_exp <- data.frame(t(vis@SigScores))
saveRDS(signature_exp,"/data/khaoxian/project/Stagecrc/save/scmetabolism/cd8t.metabolic.score.rds")

a <- rownames(signature_exp)
a <- gsub(" |\\(|\\)|\\/|,","_",a)
a <- gsub("__|___|_\\-_","_",a)
rownames(signature_exp) <- a
rownames(signature_exp) <- gsub("_"," ",rownames(signature_exp))

#2. analyze metabolic trajectory
CD8T_metabolic_trends <- compute_gene_trends(seur.palantir=cd8t.palantir,
                                             gene_exprs=t(signature_exp),
                                             branch=c("CD8T_proliferative", "CD8T_IEL_CD160"),
                                             markers=rownames(signature_exp))

p <- plot_gene_trends(CD8T_metabolic_trends,
                      legend.text.size=5,
                      legend.key.size=2)
ggsave("/data/khaoxian/project/Stagecrc/output/palantir/cd8t/module1/CD8T_metabolic_trends_palantir.pdf", p, height = 20, width = 28, units = "cm")

#metabolic trend heatmap plot
p1 <- plot_gene_trend_heatmaps(gene_trends=CD8T_metabolic_trends, 
                               genes=NULL, #rownames(signature_exp)[c(1,2,12,15,51,79)]
                               lineages="CD8T_IEL_CD160",
                               Annotategenecluster=NULL,
                               show_row_names=T,
                               right_annotation=NULL,
                               clust_by_Annotategenecluster=F,
                               clust_by_Annotategenecluster_rankded_clust = F,
                               fontsize=5,
                               link_width.size=5,
                               link_height.size=5,
                               padding.size=1)

pdf("/data/khaoxian/project/Stagecrc/output/palantir/cd8t/module1/CD8T_IEL_metabolic_trends.heatmap.pdf", height = 14, width = 10)
p1$CD8T_IEL_CD160
dev.off()

p2 <- plot_gene_trend_heatmaps(gene_trends=CD8T_metabolic_trends, 
                               genes=NULL, #rownames(signature_exp)[c(1,2,12,15,51,79)]
                               lineages="CD8T_proliferative",
                               Annotategenecluster=NULL,
                               show_row_names=T,
                               right_annotation=NULL,
                               clust_by_Annotategenecluster=F,
                               clust_by_Annotategenecluster_rankded_clust = F,
                               fontsize=5,
                               link_width.size=5,
                               link_height.size=5,
                               padding.size=1)

pdf("/data/khaoxian/project/Stagecrc/output/palantir/cd8t/module1/CD8T_proliferative_metabolic_trends.heatmap.pdf", height = 14, width = 10)
p2$CD8T_proliferative
dev.off()
