#import packages
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
from scipy.sparse import issparse
import random, torch
import warnings
import matplotlib
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cv2
import seaborn as sns
from anndata import AnnData
from numpy.random import default_rng
import glob
import scipy.io
import sys
from scipy.spatial import distance_matrix
warnings.filterwarnings("ignore")
import palantir
import scvi
from matplotlib import rcParams

#set environment                                                                                         
working_dir="/data/khaoxian/project/Stagecrc/output/palantir_correct/cd8t/module1"                                                             #<-----------------------                                                          #<-----------------------                                                  #<-----------------------

if not os.path.exists(working_dir):
  os.makedirs(working_dir)
os.chdir(working_dir)

#read dataset
#Representative cluster: CD8+ T cell
scadata=sc.read_h5ad('/data/khaoxian/project/Stagecrc/save/corrected_norm_adata/Crc_cd8t.h5ad')
sc_ref=sc.read_h5ad('/data/khaoxian/project/Stagecrc/save/corrected_norm_adata/Crc_cd8t.h5ad')
sc.pp.highly_variable_genes(scadata, n_top_genes=2000, flavor='seurat')
scadata.var["highly_variable"]=scadata.var["vst.variable"]

sc.pp.pca(scadata)
scadata.obsm['X_pca']=scadata.obsm["X_harmony"] 

# Run diffusion maps
pca_projections = pd.DataFrame(scadata.obsm['X_pca'], index=scadata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=5, knn=50)
ms_data = palantir.utils.determine_multiscale_space(dm_res)

#Run umap
sc.pp.neighbors(scadata)
sc.tl.umap(scadata)

#ploting
# Use scanpy functions to visualize umaps or FDL
sc.pl.embedding(scadata, basis='umap',show=True)
plt.savefig("./umap.png")
#ploting
sc.pl.umap(scadata, color=["Cluster"])
plt.savefig("./umapCluster.png")

#replace umap
umapbackup=scadata.obsm['X_umap'].copy()
scadata.obsm['X_umap']=sc_ref.obsm['X_umap'].copy()

#ploting
sc.pl.umap(scadata, color=["Cluster"])
plt.savefig("./umapCluster_ref.png")

#imputated data by MAGIC
scadata.layers['MAGIC_imputed_data'] = palantir.utils.run_magic_imputation(scadata, dm_res)

#plot diffusion components
umap = pd.DataFrame(scadata.obsm['X_umap'], index=scadata.obs_names)
umap.columns=["x","y"]

#ploting
palantir.plot.plot_diffusion_components(umap, dm_res)
plt.savefig("./diffusion.umap.png")

#check expression of representative markers 
sc.pl.embedding(scadata, basis='umap', layer='MAGIC_imputed_data',
                color=['MKI67', 'OLFM4', 'TFF3', 'EPCAM', "CEACAM7", "PTPRC", "CD3D", "CD8A", "MS4A1", "CD68", "COL1A1", "SSR4", "PECAM1"])
plt.savefig("./markers.umap1.png")

sc.pl.embedding(scadata, basis='umap', layer='MAGIC_imputed_data',
                color=['CD3D', 'CD8A', 'MKI67', 'PCNA', "CCR7", "SELL", "LEF1", "GZMK", "GZMH", "FGFBP2", "TRDC", "PDCD1", "TIGIT", "XCL1", "XCL2", "SELENOK"])
plt.savefig("./markers.umap2.png")

#palantir analysis
#star cell: naive cd8 CTACACCAGAGGGATA.32
pr_res = palantir.core.run_palantir(ms_data, early_cell="CTACACCAGAGGGATA.32",num_waypoints=500)
#P_P_P0410_08701 is proliferative cells

#save file
import joblib
joblib.dump(pr_res, 'cd8t_palantir_result1.pkl')

#plot pseudotime and differential potiential
palantir.plot.plot_palantir_results(pr_res, umap)
plt.savefig("./diffusion.palantir.png")

#check terminal states
pr_res.branch_probs.columns
#['CGTTAGAAGAGCTGCA.39', 'GSM5688709_CATGAGTAGGGCGAAG-1']

#annotate terminal states
terminal_states = pd.Series(['T09_CD8_proliferative', 'T06_CD8_IEL_CD160'],
                            index=['CGTTAGAAGAGCTGCA.39',
                                   'GSM5688709_CATGAGTAGGGCGAAG-1'])
pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]

#check gene trend
import rpy2
import rpy2.rinterface_lib.embedded as embedded
from rpy2.robjects.packages import importr
os.environ["R_HOME"]=r"/usr/lib/R"

genes = ['CD3D', 'CD8A', 'MKI67', 'PCNA', "CCR7", "SELL", "LEF1", "GZMK", "GZMH", "FGFBP2", "TRDC", "PDCD1", "TIGIT", "XCL1", "XCL2", "SELENOK"]
imp_df = pd.DataFrame(scadata[:, genes].layers['MAGIC_imputed_data'],
                      index=scadata.obs_names, columns=genes)
gene_trends = palantir.presults.compute_gene_trends( pr_res, imp_df.loc[:, genes])
gene_trends.keys()

palantir.plot.plot_gene_trends(gene_trends)
plt.savefig("./palantir_genetrnd1.png")

#get gene module
#add interesting genes into hvg

#using hvg to construct gene modules
genes=scadata.var[scadata.var.highly_variable==1].highly_variable.index
#2000 genes

imp_df = pd.DataFrame(scadata[:, genes].layers['MAGIC_imputed_data'],
                      index=scadata.obs_names, columns=genes)

gene_trends = palantir.presults.compute_gene_trends(pr_res,imp_df)

#save file
import joblib
joblib.dump(gene_trends, 'cd8_gene_trends.pkl')

gene_trends.keys()

#module1:CD8T_proliferative_gene_clusters
CD8T_proliferative_trends = gene_trends['T09_CD8_proliferative']['trends']
CD8T_proliferative_gene_clusters = palantir.presults.cluster_gene_trends(CD8T_proliferative_trends,k=200)

#plot clustering of genes
palantir.plot.plot_gene_trend_clusters(CD8T_proliferative_trends, CD8T_proliferative_gene_clusters)
plt.savefig("./CD8T_proliferative_gene_cluster.png")

#save gene module
CD8T_proliferative_gene_clusters=pd.DataFrame(CD8T_proliferative_gene_clusters)
CD8T_proliferative_gene_clusters.to_csv("./CD8T_proliferative_gene_clusters1.csv")

#save the trend of gene expressionand associated standard deviation
CD8T_proliferative_trends.to_csv("./CD8T_proliferative_gene_trend.csv")

CD8T_proliferative_std = gene_trends['T09_CD8_proliferative']['std']
CD8T_proliferative_std.to_csv("./CD8T_proliferative_gene_std.csv")

#module2:CD8T_IEL_CD160
CD8T_IEL_CD160_trends = gene_trends['T06_CD8_IEL_CD160']['trends']
CD8T_IEL_CD160_gene_clusters = palantir.presults.cluster_gene_trends(CD8T_IEL_CD160_trends,k=200)

#plot trend cluster
palantir.plot.plot_gene_trend_clusters(CD8T_IEL_CD160_trends, CD8T_IEL_CD160_gene_clusters)
plt.savefig("./CD8T_IEL_CD160_gene_cluster.png")

#save gene module
CD8T_IEL_CD160_gene_clusters=pd.DataFrame(CD8T_IEL_CD160_gene_clusters)
CD8T_IEL_CD160_gene_clusters.to_csv("./CD8T_IEL_CD160_gene_clusters1.csv")

#save gene trend and std
CD8T_IEL_CD160_trends.to_csv("./CD8T_IEL_CD160_gene_trend.csv")

CD8T_IEL_CD160_std = gene_trends['T06_CD8_IEL_CD160']['std']
CD8T_IEL_CD160_std.to_csv("./CD8T_IEL_CD160_gene_std.csv")

scadata.obs["pseudotime"]=pr_res.pseudotime
scadata.obs["entropy"]=pr_res.entropy


for id in pr_res.branch_probs.columns:
  scadata.obs[id]=pr_res.branch_probs[id]

sc.pl.umap(scadata, color="pseudotime")
plt.savefig("./palantir_pseudotime.png")

#save h5ad
scadata.write_h5ad("./Crc_cd8t.module1.h5ad")
