#scenic work flow
#1. prefiltered scRNA-seq data (must by seurat)
#2. format data (by R)
#3. run scenic.py
#4. analyze scenic output
#usage
#scenic.py ouput_dir input_h5ad cluster_column
#scenic.py /data/khaoxian/project/cytof/output/scenic /data/khaoxian/project/cytof/save/imm.h5ad Clusters
#require /home/khaoxian/anaconda3/envs/pyscenic/bin/python


# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import sys
from MulticoreTSNE import MulticoreTSNE as TSNE

#set wdir
wdir = sys.argv[1]      #parameter[1]
os.chdir( wdir )

f_mtx_dir = sys.argv[2]               #parameter[2]

adata = sc.read_h5ad( f_mtx_dir )

#1. grn analysis
# create basic row and column attributes for the loom file:
f_loom_path_scenic = "adata_scenic.loom"

row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)        

os.system("/home/khaoxian/anaconda3/envs/pyscenic/bin/pyscenic grn adata_scenic.loom /data/khaoxian/config/scenic/allTFs_hg38.txt -o adj.csv --num_workers 30")

#2. ctx analysis
os.system("/home/khaoxian/anaconda3/envs/pyscenic/bin/pyscenic ctx adj.csv \
             /data/khaoxian/config/scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
             /data/khaoxian/config/scenic/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
             --annotations_fname /data/khaoxian/config/scenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
             --expression_mtx_fname adata_scenic.loom \
             --output reg.csv \
             --mask_dropouts \
             --num_workers 30")

#3. auccell analysis
os.system("/home/khaoxian/anaconda3/envs/pyscenic/bin/pyscenic aucell \
    adata_scenic.loom \
    reg.csv \
    --output pyscenic_output.loom \
    --num_workers 20")

#4. export file
#Visualization of SCENIC's AUC matrix
import json
import zlib
import base64

f_pyscenic_output="pyscenic_output.loom"
# collect SCENIC AUCell output
lf = lp.connect("pyscenic_output.loom", mode='r+', validate=False )#           #trash2 but important
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

import umap
# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_umap.txt", sep='\t')
# tSNE
tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_tsne.txt", sep='\t')

# scenic output
lf = lp.connect("pyscenic_output.loom", mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))

auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
dr_umap = pd.read_csv( 'scenic_umap.txt', sep='\t', header=0, index_col=0 )           #trash3
dr_tsne = pd.read_csv( 'scenic_tsne.txt', sep='\t', header=0, index_col=0 )           #trash4

#Fix regulon objects to display properly in SCope
auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
# regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )


tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])

Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
Embeddings_X = pd.concat( [
    pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
    pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
    dr_tsne['X'] ,
    dr_umap['X']
], sort=False, axis=1, join='outer' )
Embeddings_X.columns = ['1','2','3','4']

Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
Embeddings_Y = pd.concat( [
    pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
    pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
    dr_tsne['Y'] ,
    dr_umap['Y']
], sort=False, axis=1, join='outer' )
Embeddings_Y.columns = ['1','2','3','4']

### metadata
metaJson = {}

metaJson['embeddings'] = [
    {
        "id": -1,
        "name": f"Scanpy t-SNE (highly variable genes)"
    },
    {
        "id": 1,
        "name": f"Scanpy UMAP  (highly variable genes)"
    },
    {
        "id": 2,
        "name": "Scanpy PC1/PC2"
    },
    {
        "id": 3,
        "name": "SCENIC AUC t-SNE"
    },
    {
        "id": 4,
        "name": "SCENIC AUC UMAP"
    },
]

metaJson["clusterings"] = [{
    "id": 0,
    "group": "Scanpy",
    "name": "Scanpy louvain default resolution",
    "clusters": [],
}]

metaJson["metrics"] = [
    {
        "name": "nUMI"
    }, {
        "name": "nGene"
    }, {
        "name": "Percent_mito"
    }
]

#add annotation
cluster_annotated=sys.argv[3]
#cluster_annotated="Clusters"
metaJson["annotations"] = [
    {
        "name": "Louvain_clusters_Scanpy",
        "values": list(set( adata.obs[cluster_annotated].astype(str) )) #np.str -- > str
    },
    #{
    #    "name": "Genotype",
    #    "values": list(set(adata.obs['Genotype'].values))
    #},
    #{
    #    "name": "Timepoint",
    #    "values": list(set(adata.obs['Timepoint'].values))
    #},
    #{
    #    "name": "Sample",
    #    "values": list(set(adata.obs['Sample'].values))
    #}
]

# SCENIC regulon thresholds:
metaJson["regulonThresholds"] = rt

#cluster_annotated=Clusters
for i in range(max(set([int(x) for x in adata.obs[cluster_annotated]])) + 1):
    clustDict = {}
    clustDict['id'] = i
    clustDict['description'] = f'Unannotated Cluster {i + 1}'
    metaJson['clusterings'][0]['clusters'].append(clustDict)

clusterings = pd.DataFrame()
clusterings["0"] = adata.obs[cluster_annotated].values.astype(np.int64)
#Assemble loom file row and column attributes

def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr
#cluster_annotated
col_attrs = {
    "CellID": np.array(adata.obs.index),
    "nUMI": np.array(adata.obs['nCount_RNA'].values),        #'n_counts'
    "nGene": np.array(adata.obs['nFeature_RNA'].values),        #'n_genes'
    "Louvain_clusters_Scanpy": np.array( adata.obs[cluster_annotated].values ),
    #"Genotype": np.array(adata.obs['Genotype'].values),
    #"Timepoint": np.array(adata.obs['Timepoint'].values),
    #"Sample": np.array(adata.obs['Sample'].values),
    "Percent_mito": np.array(adata.obs['pMT'].values),   #percent_mito
    "Embedding": dfToNamedMatrix(tsneDF),
    "Embeddings_X": dfToNamedMatrix(Embeddings_X),
    "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
    "RegulonsAUC": dfToNamedMatrix(auc_mtx),
    "Clusterings": dfToNamedMatrix(clusterings),
    "ClusterID": np.array(adata.obs[cluster_annotated].values)
}

row_attrs = {
    "Gene": lf.ra.Gene,
    "Regulons": regulons,
}

attrs = {
    "title": "sampleTitle",
    "MetaData": json.dumps(metaJson),
    "Genome": 'hg38',
    "SCopeTreeL1": "",
    "SCopeTreeL2": "",
    "SCopeTreeL3": ""
}

# compress the metadata field:
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

#Create a new loom file, copying the expression matrix from the open loom connection:
f_final_loom = 'scenic_integrated-output.loom'

lp.create(
    filename = f_final_loom ,
    layers=lf[:,:],
    row_attrs=row_attrs,
    col_attrs=col_attrs,
    file_attrs=attrs
)
lf.close()













