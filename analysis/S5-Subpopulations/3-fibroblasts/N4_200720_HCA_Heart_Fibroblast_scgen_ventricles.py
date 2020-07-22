##########
#
# Running scGen for VENTRICULAR fibroblasts
# Eric Lindberg, MDC Berlin
# Date: 2020 07 20
# eric.lindberg@mdc-berlin.de
# 
##########

## This was run on a Centos7 machine with NVIDIA-SMI 430.46 and CUDA 10.1
##

# Libraries
import scanpy as sc
import tensorflow as tf
import scgen
import numpy as np
import re
import random

# Global settings
conf = tf.ConfigProto()
conf.gpu_options.per_process_gpu_memory_fraction=0.85
session = tf.Session(config=conf)
random_state=10

# Functions needed for downstream analysis
def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def batch_removal_ct5(network, adata, batch_key="batch", cell_label_key="cell_type"):
    from scipy.sparse import issparse
    import anndata
    import pandas as pd
    if issparse(adata.X):
        latent_all = network.to_latent(adata.X.A)
    else:
        latent_all = network.to_latent(adata.X)
    adata_latent = anndata.AnnData(latent_all)
    adata_latent.obs = adata.obs.copy(deep=True)
    unique_cell_types = np.unique(adata_latent.obs[cell_label_key])
    shared_ct = []
    not_shared_ct = []
    for cell_type in unique_cell_types:
        temp_cell = adata_latent[adata_latent.obs[cell_label_key] == cell_type]
        if len(np.unique(temp_cell.obs[batch_key])) < 2:
            cell_type_ann = adata_latent[adata_latent.obs[cell_label_key] == cell_type]
            not_shared_ct.append(cell_type_ann)
            continue
        temp_cell = adata_latent[adata_latent.obs[cell_label_key] == cell_type]
        batch_list = {}
        batch_ind = {}
        max_batch = 0
        max_batch_ind = ""
        batches = np.unique(temp_cell.obs[batch_key])
        for i in batches:
            temp = temp_cell[temp_cell.obs[batch_key] == i]
            temp_ind = temp_cell.obs[batch_key] == i
            if max_batch < len(temp):
                max_batch = len(temp)
                max_batch_ind = i
            batch_list[i] = temp
            batch_ind[i] = temp_ind
        max_batch_ann = batch_list[max_batch_ind]
        for study in batch_list:
            delta = np.average(max_batch_ann.X, axis=0) - np.average(batch_list[study].X, axis=0)
            batch_list[study].X = delta + batch_list[study].X
            temp_cell[batch_ind[study]].X = batch_list[study].X
        shared_ct.append(temp_cell)
    all_shared_ann = anndata.AnnData.concatenate(*shared_ct, batch_key="concat_batch", index_unique=None)
    if "concat_batch" in all_shared_ann.obs.columns:
        del all_shared_ann.obs["concat_batch"]
    if len(not_shared_ct) < 1:
        corrected = anndata.AnnData(network.reconstruct(all_shared_ann.X, use_data=True))
        corrected.obs = all_shared_ann.obs.copy(deep=True)
        corrected.var_names = adata.var_names.tolist()
        corrected = corrected[adata.obs_names]
        if adata.raw is not None:
            adata_raw = anndata.AnnData(X=adata.raw.X, var=adata.raw.var)
            adata_raw.obs_names = adata.obs_names
            corrected.raw = adata_raw
        return corrected
    else:
        all_not_shared_ann = anndata.AnnData.concatenate(*not_shared_ct, batch_key="concat_batch", index_unique=None)
        all_corrected_data = anndata.AnnData.concatenate(all_shared_ann, all_not_shared_ann, batch_key="concat_batch",  index_unique=None)
        if "concat_batch" in all_shared_ann.obs.columns:
            del all_corrected_data.obs["concat_batch"]
        corrected = anndata.AnnData(network.reconstruct(all_corrected_data.X, use_data=True), )
        corrected.obs = pd.concat([all_shared_ann.obs, all_not_shared_ann.obs])
        corrected.var_names = adata.var_names.tolist()
        corrected = corrected[adata.obs_names]
        if adata.raw is not None:
            adata_raw = anndata.AnnData(X=adata.raw.X, var=adata.raw.var)
            adata_raw.obs_names = adata.obs_names
            corrected.raw = adata_raw
        return corrected


print("SETUP COMPLETE")

# Load files needed for analysis
path_out='/fast/AG_Huebner/huebner3/ANALYSES/20191112_cs_HH_England/scanpy/20200423/FB/'

nuclei=sc.read_h5ad(path_out + "fibroblasts_nucleis_filt_V2_VT_LEIDEN_RAW.h5ad")
cells=sc.read_h5ad(path_out + "fibroblasts_cells_filt_V2_VT_LEIDEN_RAW.h5ad")

# Assign source
nuclei.obs['source']="nuclei"
cells.obs['source']="cells"

# Process nuclei
tmp_nuclei=nuclei.copy()
sc.pp.normalize_per_cell(tmp_nuclei, counts_per_cell_after=1e4)
sc.pp.log1p(tmp_nuclei)
tmp_nuclei.raw=tmp_nuclei
sc.pp.highly_variable_genes(tmp_nuclei, n_top_genes=7000)

# Process cells
tmp_cells=cells.copy()
sc.pp.normalize_per_cell(tmp_cells, counts_per_cell_after=1e4)
sc.pp.log1p(tmp_cells)
tmp_cells.raw=tmp_cells
sc.pp.highly_variable_genes(tmp_cells, n_top_genes=7000)


# Concatenate objects

all=nuclei.concatenate([cells])
all.write(path_out + "ALL_corrected_FB_V_double_RAW.h5ad")
## Process
sc.pp.normalize_per_cell(all, counts_per_cell_after=1e4)
sc.pp.log1p(all)
all.raw=all

# Continue with the intersect of HVGs
all = all[:, ((tmp_nuclei.var['highly_variable'] | tmp_cells.var['highly_variable']))].copy()

# Defining a Trainingsset will decrease computational resources, as we only have few cells, we use all for training
# All cells and nuclei are embedded in the final manifold
nuclei=sc.pp.subsample(nuclei, n_obs=10000, random_state=random_state, copy=True) 
all_train=nuclei.concatenate([cells])
sc.pp.normalize_per_cell(all_train, counts_per_cell_after=1e4)
sc.pp.log1p(all_train)
all_train = all_train[:, ((tmp_nuclei.var['highly_variable'] | tmp_cells.var['highly_variable']))].copy()



print("DATA PROCESSING COMPLETE")

# Round 1: Correct for source
## Train scGen on the trainings dataset and apply batch removal to both trainings and the complete set of nuclei/cells
all.obs['cell_type']=all.obs['leiden_annotated'].values
all.obs['batch']=all.obs['cell_source'].values
all_train.obs['cell_type']=all_train.obs['leiden_annotated'].values
all_train.obs['batch']=all_train.obs['cell_source'].values

network = scgen.VAEArith(x_dimension= all_train.shape[1], model_path= path_out + 'scgen_model/scgen_model_V_source')
print("Network made")
network.train(train_data=all_train, n_epochs=20)
print("Network trained")

corrected_adata =  batch_removal_ct5(network, all, batch_key="batch", cell_label_key="cell_type")
corrected_adata_train =  batch_removal_ct5(network, all_train, batch_key="batch", cell_label_key="cell_type")

corrected_adata.write(path_out + "ALL_corrected_FB_ATV_source.h5ad")


# Round 2: Correct for donor
all=corrected_adata.copy()
all_train=corrected_adata_train.copy()

all.obs['cell_type']=all.obs['leiden_annotated'].values
all.obs['batch']=all.obs['donor'].values
all_train.obs['cell_type']=all_train.obs['leiden_annotated'].values
all_train.obs['batch']=all_train.obs['donor'].values


network = scgen.VAEArith(x_dimension= all_train.shape[1], model_path= path_out + 'scgen_model/scgen_model_V_sourcedonor')
print("Network made")
network.train(train_data=all_train, n_epochs=20)
print("Network trained")
corrected_adata =  batch_removal_ct5(network, all, batch_key="batch", cell_label_key="cell_type")

# Finished with batch alignment

corrected_adata.write(path_out + "ALL_corrected_FB_V_double.h5ad")

sc.pp.neighbors(corrected_adata, random_state = random_state)
sc.tl.umap(corrected_adata, min_dist = 0.3, spread = 2, random_state = random_state)

corrected_adata.write(path_out + "ALL_corrected_FB_V_double.h5ad")

print("FINISHED")
