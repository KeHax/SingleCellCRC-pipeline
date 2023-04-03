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
import squidpy as sq
import seaborn as sns
from anndata import AnnData
from numpy.random import default_rng
import glob
import scipy.io
import sys
from scipy.spatial import distance_matrix
warnings.filterwarnings("ignore")
import cell2location
import scvi
from matplotlib import rcParams

#set environment                                                                                          
working_dir="/data/khaoxian/project/Stagecrc/output/cell2location"                                                          
scRNA_path="/data/khaoxian/project/cytof/data/SC_filtered/sampled_stagecrc"                                                        
output_dir="/data/khaoxian/project/Stagecrc/output/cell2location"
os.chdir(working_dir)

#ScRNA-seq data processing
#read scRNA-seq data
scprojects="stage_crc"                                                
scprefixes=scprojects+"_"
scfeatures_file = scprefixes+"features.tsv"
scbarcodes_file=scprefixes+"barcodes.tsv"
sccounts_file=scprefixes+"matrix.mtx"
scmeta_file=scprefixes+"meta.txt"

features_path = os.path.join(scRNA_path, scfeatures_file)
barcodes_path = os.path.join(scRNA_path, scbarcodes_file)
matrix_path = os.path.join(scRNA_path, sccounts_file)
meta_path = os.path.join(scRNA_path, scmeta_file)

gene_names = [row[0] for row in csv.reader(open(features_path), delimiter="\t")]
barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]
counts = sc.read_mtx(matrix_path)

#construct adata
scadata = counts.copy()
scadata = scadata.transpose()
scadata.obs_names = barcodes
scadata.var_names = gene_names
#read metatdata
scmeta=pd.read_csv(meta_path,sep="\t",header=0,na_filter=False)

for id in scmeta.columns:
  scadata.obs[id]=scmeta[id]

#write h5ad
scadata.write_h5ad("/data/khaoxian/project/Stagecrc/output/cell2location/stagecrc_sampled1000.raw.h5ad")

#if read
scadata=sc.read_h5ad("/data/khaoxian/project/Stagecrc/output/cell2location/stagecrc_sampled1000.raw.h5ad")

#select genes
from cell2location.utils.filtering import filter_genes
selected = filter_genes(scadata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
scadata = scadata[:, selected].copy()

##################################################################################################################
#train scRNA-seq data
#Estimation of reference cell type signatures (NB regression)
import cell2location
import scvi
from matplotlib import rcParams
# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=scadata,
                                                   batch_key='orig.ident',
                                                   labels_key='Cluster'
)

# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(scadata)
# view anndata_setup as a sanity check
mod.view_anndata_setup()
mod.train(max_epochs=250, use_gpu=True)


# create paths and names to results folders for reference regression and cell2location models
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
scadata = mod.export_posterior(
  scadata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
ref_run_name = "/data/khaoxian/project/Stagecrc/output/cell2location"
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/stagecrc_module.h5ad"
scadata.write(adata_file)
adata_file

#make signature
if 'means_per_cluster_mu_fg' in scadata.varm.keys():
  inf_aver = scadata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in scadata.uns['mod']['factor_names']]].copy()
else:
  inf_aver = scadata.var[[f'means_per_cluster_mu_fg_{i}' for i in scadata.uns['mod']['factor_names']]].copy()
inf_aver.columns = scadata.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

#cell2location analysis
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
import squidpy as sq
import seaborn as sns
from anndata import AnnData
from numpy.random import default_rng
import glob
import scipy.io
import sys
from scipy.spatial import distance_matrix
warnings.filterwarnings("ignore")
import cell2location
import scvi
from matplotlib import rcParams

#set environment
working_dir="/data/khaoxian/project/Stagecrc/output/cell2location/border6/stagecrc"                                                             #<-----------------------                                            #<-----------------------
output_dir="/data/khaoxian/project/Stagecrc/output/cell2location/border6/stagecrc" 

if not os.path.exists(working_dir):
    os.makedirs(working_dir)


import cell2location
import scvi
from matplotlib import rcParams
from cell2location.models import RegressionModel

ref_run_name = "/data/khaoxian/project/Stagecrc/output/cell2location"
adata_file = f"{ref_run_name}/stagecrc_module.h5ad"
scadata=sc.read_h5ad(adata_file)
#read in sc mode
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", scadata)

##read in st data------------------
#set input data--------
stdata_dir="/data/khaoxian/project/cytof/data/ST_filtered/border6"
projects="border6"
prefixes=projects+"_"
features_file = prefixes+"features.tsv"
barcodes_file=prefixes+"barcodes.tsv"
counts_file=prefixes+"matrix.mtx"
position_file=prefixes+"position.txt"

features_path = os.path.join(stdata_dir, features_file)
barcodes_path = os.path.join(stdata_dir, barcodes_file)
matrix_path = os.path.join(stdata_dir, counts_file)
spatial_path = os.path.join(stdata_dir, position_file)
img_path = glob.glob(os.path.join(stdata_dir, '*.png'))
img_path=img_path[0]
scalefactor_path=os.path.join(stdata_dir, "scale_factor.txt")
#input data
gene_names = [row[0] for row in csv.reader(open(features_path), delimiter="\t")]
barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]
counts = sc.read_mtx(matrix_path)
scale_factor=pd.read_csv(scalefactor_path,sep=" ",header=None,na_filter=False)
scale_factor=scale_factor.iloc[[0],[0]].values[0][0]
#input spatial information
spatial=pd.read_csv(spatial_path,sep=" ",header=None,na_filter=False,index_col=0)
coord = spatial.iloc[:,[-2,-1]]
coord=np.array(coord)

#input image
img=cv2.imread(img_path)

#construct adata
adata = counts.copy()
adata = adata.transpose()

#add barcode and features
adata.obs_names = barcodes
adata.var_names = gene_names

adata.obsm={"spatial": coord}


#add image
spatial_key = "spatial"
library_id = projects
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": img}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": scale_factor, "spot_diameter_fullres": 90}

#add meta.data
meta_file=prefixes+"meta.txt"
meta_path = os.path.join(stdata_dir, meta_file)
metadata=pd.read_csv(meta_path,sep="\t",header=0,na_filter=False)
adata.obs["Cluster"]=metadata["Cluster"].tolist()

#Calculate MT- features
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

#save raw
adata_backup= adata.copy()

#format for cell2location
adata.obs['sample'] = projects
adata.var['SYMBOL'] = adata.var_names
adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)

# Calculate QC metrics
from scipy.sparse import csr_matrix
adata.X = adata.X.toarray()
sc.pp.calculate_qc_metrics(adata, inplace=True)
adata.X = csr_matrix(adata.X)
adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

# add sample name to obs names
adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
adata.obs_names = adata.obs["sample"] \
                  + '_' + adata.obs_names
adata.obs.index.name = 'spot_id'

# mitochondria-encoded (MT) genes should be removed for spatial mapping
adata.obsm['mt'] = adata[:, adata.var['mt'].values].X.toarray()
adata = adata[:, ~adata.var['mt'].values]

#cell2location
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in scadata.varm.keys():
    inf_aver = scadata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in scadata.uns['mod']['factor_names']]].copy()
else:
    inf_aver = scadata.var[[f'means_per_cluster_mu_fg_{i}' for i in scadata.uns['mod']['factor_names']]].copy()

inf_aver.columns = scadata.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]
#obs['x'].astype(str)
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata.var_names, inf_aver.index)
adata = adata[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()


# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata, batch_key="sample")

# create and train the model
mod = cell2location.models.Cell2location(
    adata, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=10,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()

#cell2location
mod.train(max_epochs=30000,
          batch_size=None,
          train_size=1,
          use_gpu=True)
          
run_name=stdata_dir+"/"+"stagecrc_cell2location_map"
mod.save(f"{run_name}", overwrite=True)

adata = mod.export_posterior(
    adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata.write(adata_file)
adata_file


#fi read
adata=sc.read_h5ad(adata_file)

#plot abundance########################################

adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": scale_factor, "spot_diameter_fullres": 450}
fig = mod.plot_spatial_QC_across_batches()
plt.savefig(os.path.join(output_dir,"qc.pdf"))
plt.show()

adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": scale_factor, "spot_diameter_fullres": 90}

adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
clusters=adata.uns['mod']['factor_names']  #.tolist()

sc.pl.spatial(adata, cmap='magma',
              color=clusters[0:8],
              ncols=4, size=6,
              img_key='hires',
              # limit color scale at 99.2% quantile of cell abundance
              vmin=0, vmax='p99.2',
              show=False
              )
plt.savefig(os.path.join(output_dir,"crc1_c1.pdf"))

sc.pl.spatial(adata, cmap='magma',
              color=clusters[8:16],
              ncols=4, size=6,
              img_key='hires',
              # limit color scale at 99.2% quantile of cell abundance
              vmin=0, vmax='p99.2',
              show=False
              )
plt.savefig(os.path.join(output_dir,"crc1_c2.pdf"))

sc.pl.spatial(adata, cmap='magma',
              color=clusters[16:24],
              ncols=4, size=6,
              img_key='hires',
              # limit color scale at 99.2% quantile of cell abundance
              vmin=0, vmax='p99.2',
              show=False
              )
plt.savefig(os.path.join(output_dir,"crc1_c3.pdf"))

sc.pl.spatial(adata, cmap='magma',
              color=clusters[24:32],
              ncols=4, size=6,
              img_key='hires',
              # limit color scale at 99.2% quantile of cell abundance
              vmin=0, vmax='p99.2',
              show=False
              )
plt.savefig(os.path.join(output_dir,"crc1_c4.pdf"))

sc.pl.spatial(adata, cmap='magma',
              color=clusters[32:40],
              ncols=4, size=6,
              img_key='hires',
              # limit color scale at 99.2% quantile of cell abundance
              vmin=0, vmax='p99.2',
              show=False
              )
plt.savefig(os.path.join(output_dir,"crc1_c5.pdf"))

sc.pl.spatial(adata, cmap='magma',
              color=clusters[40:48],
              ncols=4, size=6,
              img_key='hires',
              # limit color scale at 99.2% quantile of cell abundance
              vmin=0, vmax='p99.2',
              show=False
              )
plt.savefig(os.path.join(output_dir,"crc1_c6.pdf"))

sc.pl.spatial(adata, cmap='magma',
              color=clusters[48:53],
              ncols=4, size=6,
              img_key='hires',
              # limit color scale at 99.2% quantile of cell abundance
              vmin=0, vmax='p99.2',
              show=False
              )
plt.savefig(os.path.join(output_dir,"crc1_c7.pdf"))

#Identifying cellular compartments / tissue zones using matrix factorisation (NMF)
from cell2location import run_colocation
res_dict, adata = run_colocation(
    adata,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
        'n_fact': np.arange(9, 12), # IMPORTANT: use a wider range of the number of factors (5-30)
        'sample_name_col': 'sample', # columns in adata_vis.obs thaoutput_dirt identifies sample
        'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{output_dir}/CoLocatedComb_subset/'}
)

res_dict['n_fact9']['mod'].plot_cell_type_loadings()
plt.savefig(os.path.join(output_dir,"crc1_nmf_factor9.pdf"))
res_dict['n_fact10']['mod'].plot_cell_type_loadings()
plt.savefig(os.path.join(output_dir,"crc1_nmf_factor10.pdf"))
res_dict['n_fact11']['mod'].plot_cell_type_loadings()
plt.savefig(os.path.join(output_dir,"crc1_nmf_factor11.pdf"))


















