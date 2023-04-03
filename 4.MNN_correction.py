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
import mnnpy

#analyze trajectory of epi cells
projects="Crc_epi"                                                                                                   
working_dir="/data/khaoxian/project/Stagecrc/output/palantir/epi/module1"                                                          
scRNA_path="/data/khaoxian/project/Stagecrc/save/epi_palantir"                                                                                  

os.chdir(working_dir)

#readin sc.adata and transfer harmony, hvg
scadata=sc.read_h5ad("/data/khaoxian/project/Stagecrc/save/epi_palantir/Crc_epi.h5ad")

from scipy.sparse import csr_matrix

##normalize
sc.pp.normalize_per_cell(scadata)
sc.pp.log1p(scadata,base=2)

# split per batch into new objects.
adata_mnn = scadata.copy()
adata_list = [adata_mnn[adata_mnn.obs['orig.ident'] == i] for i in adata_mnn.obs['orig.ident'].unique()]
cdata = sc.external.pp.mnn_correct(*adata_list,
                                   svd_dim = 50, batch_key = 'orig.ident', save_raw = True, var_subset = scadata.var.index)

corr_data = cdata[0]
corr_data.write_h5ad("/data/khaoxian/project/Stagecrc/save/MNNcorrect/epi.mnncorrect.h5ad")


#analyze trajectory of CD4 cells 
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
import mnnpy
#set environment
#/data/khaoxian/project/Stagecrc/save/epi_palantir
working_dir="/data/khaoxian/project/Stagecrc/output/palantir/cd4t/module1"                                                           
scRNA_path="/data/khaoxian/project/Stagecrc/save/cd4t_palantir"                                                         

if not os.path.exists(working_dir):
  os.makedirs(working_dir)

os.chdir(working_dir)

#readin sc.adata and transfer harmony, hvg
scadata=sc.read_h5ad("/data/khaoxian/project/Stagecrc/save/cd4t_palantir/Crc_cd4t_senic.h5ad")

##normalize
sc.pp.normalize_per_cell(scadata)
sc.pp.log1p(scadata,base=2)

# split datasets per batch into new objects.
adata_mnn = scadata.copy()
adata_list = [adata_mnn[adata_mnn.obs['orig.ident'] == i] for i in adata_mnn.obs['orig.ident'].unique()]
cdata = sc.external.pp.mnn_correct(*adata_list,
                                   svd_dim = 50, batch_key = 'orig.ident', save_raw = True, var_subset = scadata.var.index)

corr_data = cdata[0]
corr_data.write_h5ad("/data/khaoxian/project/Stagecrc/save/MNNcorrect/cd4t.mnncorrect.h5ad")



#analyze trajectory of CD8 cells
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
import mnnpy
#set environment
#/data/khaoxian/project/Stagecrc/save/epi_palantir
working_dir="/data/khaoxian/project/Stagecrc/output/palantir/cd8t/module1"                                                         
scRNA_path="/data/khaoxian/project/Stagecrc/save/cd8t_palantir"                                                        

if not os.path.exists(working_dir):
  os.makedirs(working_dir)

os.chdir(working_dir)
#readin sc.adata and transfer harmony, hvg
scadata=sc.read_h5ad("/data/khaoxian/project/Stagecrc/save/cd8t_palantir/Crc_cd8t.h5ad")

##normalize
sc.pp.normalize_per_cell(scadata)
sc.pp.log1p(scadata,base=2)

# split datasets per batch into new objects.
adata_mnn = scadata.copy()
adata_list = [adata_mnn[adata_mnn.obs['orig.ident'] == i] for i in adata_mnn.obs['orig.ident'].unique()]
cdata = sc.external.pp.mnn_correct(*adata_list,
                                   svd_dim = 50, batch_key = 'orig.ident', save_raw = True, var_subset = scadata.var.index)

corr_data = cdata[0]
corr_data.write_h5ad("/data/khaoxian/project/Stagecrc/save/MNNcorrect/cd8t.mnncorrect.h5ad")



#analyze trajectory of B cells 
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
import mnnpy
#set environment
#/data/khaoxian/project/Stagecrc/save/epi_palantir
working_dir="/data/khaoxian/project/Stagecrc/output/palantir/bcell/module1"                                                            
scRNA_path="/data/khaoxian/project/Stagecrc/save/b_palantir"                                                        

if not os.path.exists(working_dir):
  os.makedirs(working_dir)

os.chdir(working_dir)
#readin sc.adata and transfer harmony, hvg
scadata=sc.read_h5ad("/data/khaoxian/project/Stagecrc/save/b_palantir/Crc_b.h5ad")

##normalize
sc.pp.normalize_per_cell(scadata)
sc.pp.log1p(scadata,base=2)

# split datasets per batch into new objects.
adata_mnn = scadata.copy()
adata_list = [adata_mnn[adata_mnn.obs['orig.ident'] == i] for i in adata_mnn.obs['orig.ident'].unique()]
cdata = sc.external.pp.mnn_correct(*adata_list,
                                   svd_dim = 50, batch_key = 'orig.ident', save_raw = True, var_subset = scadata.var.index)

corr_data = cdata[0]
corr_data.write_h5ad("/data/khaoxian/project/Stagecrc/save/MNNcorrect/b.mnncorrect.h5ad")


#analyze trajectory of myeloid cells
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
import mnnpy
#set environment
working_dir="/data/khaoxian/project/Stagecrc/output/palantir/myeloid/module1"                                                       
scRNA_path="/data/khaoxian/project/Stagecrc/save/myeloid_palantir"                                                      

if not os.path.exists(working_dir):
  os.makedirs(working_dir)

os.chdir(working_dir)
#readin sc.adata and transfer harmony, hvg
scadata=sc.read_h5ad("/data/khaoxian/project/Stagecrc/save/myeloid_palantir/Crc_myeloid.h5ad")

##normalize
sc.pp.normalize_per_cell(scadata)
sc.pp.log1p(scadata,base=2)

# split datasets per batch into new objects.
adata_mnn = scadata.copy()
adata_list = [adata_mnn[adata_mnn.obs['orig.ident'] == i] for i in adata_mnn.obs['orig.ident'].unique()]
cdata = sc.external.pp.mnn_correct(*adata_list,
                                   svd_dim = 50, batch_key = 'orig.ident', save_raw = True, var_subset = scadata.var.index)

corr_data = cdata[0]
corr_data.write_h5ad("/data/khaoxian/project/Stagecrc/save/MNNcorrect/myeloid.mnncorrect.h5ad")



