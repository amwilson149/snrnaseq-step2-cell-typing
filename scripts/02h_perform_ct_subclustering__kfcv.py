import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import os
import sys
import pickle
import json
import argparse
import yaml
from utils.get_rand_seed import *
from utils.sparse_funcs import *

# This script takes expression data with
# validated, generalized cell types (meaning
# that some initial algorithm-generated
# clusters have been assigned cell types
# and possibly recombined), and it performs
# subclustering on a specified subset
# of those cell type clusters.
# These specified clusters may e.g. appear to be
# composed of multiple manifolds in UMAP space.
# This step facilitates more specific delineation
# of nuclei into cell types.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Subclustering of selected major cell type' + \
        ' clusters in expression data.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
parser.add_argument(f'--pass-number',type=str,required=True)
args = parser.parse_args()

# 01. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)

# 01a. Cell typing pass number.
# Convert this string to a 2-digit,
# zero padded integer if not in that
# form already
pn = args.pass_number
pn = pn.rjust(2,'0')

# 01b. Dataset information
# 01b.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01b.ii. Input (original) dataset parent directory
dataset_parent_dir = qc_root_dir + f'/' + cfg.get('data_parent_dir')
# 01b.iii. Output (intermediate) dataset parent directory
output_parent_dir = qc_root_dir + f'/' + cfg.get('out_parent_dir')
# 01b.iv. Input expression data directory
input_data_dir = qc_root_dir + f'/' + cfg.get(f'input_data_dir__{pn}')
# 01b.v. Input expression data file name
input_data_fn = cfg.get(f'input_data_fn__{pn}')
# 01b.vi. Output expression data directory
output_data_dir = qc_root_dir + f'/' + cfg.get(f'output_data_dir__{pn}')
# 01b.vii. Output expression data file name
output_data_fn = cfg.get(f'output_data_fn__{pn}')

# 01c. Subclustering information
# 01c.i. Input cluster label
ct_label_col = cfg.get(f'cluster_label__{pn}')
# 01c.ii. List of clusters to subcluster 
cts_to_subcluster = cfg.get(f'ct_groups_to_subcluster__{pn}')
# 01c.iii. Post-subclustering cluster label
post_subclust_ct_label = cfg.get(f'output_cluster_label__{pn}')

# 01d. Random state information
input_random_state_dir = qc_root_dir + f'/' + cfg.get(f'input_rs_dir__{pn}')
input_random_state_fn = cfg.get(f'input_rs_fn__{pn}')
output_random_state_dir = qc_root_dir + f'/' + cfg.get(f'output_rs_dir__{pn}')
output_random_state_fn = cfg.get(f'output_rs_fn__{pn}')

# Run subclustering

# 02. Import expression data
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 03. Import random number generator state and
# assign it to the current random number generator
working_random_state = None
random_seed_full_fn = f'{input_random_state_dir}/{input_random_state_fn}'
with open(random_seed_full_fn,'rb') as f:
    working_random_state = pickle.load(f)
numpy_random_state = np.random.RandomState()
numpy_random_state.set_state(working_random_state)

# 04. Define a partial reproduction of the original QC workflow 
# applied to expression data, but without QC filtering
# steps, to allow for subclustering
def qc_and_typical_workflow_partial(adata: ad.AnnData,
        hvg_flavor: str = 'seurat_v3',
        hvg_scaling_cap: float = None,
        n_hvgs_to_keep: int = 2000,
        n_pcs: int = 25,
        leiden_resolution: float = 0.5,
        working_dtype_int: str = 'int64',
        working_dtype_float: str = 'float64',
        scaling_accumulator_dtype: str = 'float64'):
    """Repeat HVG detection, network mapping, and cluster detection\nof a typical scanpy workflow"""
    # Typical scanpy workflow
    # Dataset-wide standardization for highly variable gene (hvg) detection
    # We default to using 'seurat_v3'-flavored hvg detection, since it
    # incorporates some harmonization into its approach.
    # This flavor expects raw counts data (not lognormalized data as the other
    # flavors do).
    adata_qc = adata.copy()
    # 04a. Clean observation data to remove
    # old leiden cluster information
    obs_cols_to_remove = ['leiden']
    for obs_col in obs_cols_to_remove:
        if obs_col in adata_qc.obs.columns:
            del adata_qc.obs[obs_col]
    # 04b. Clean gene data to remove old
    # highly variable gene information
    var_cols_to_remove = ['highly_variable',
            'highly_variable_rank',
            'means',
            'variances',
            'variances_norm']
    for var_col in var_cols_to_remove:
        if var_col in adata_qc.var.columns:
            del adata_qc.var[var_col]
    # 04c. Clean old scaling and network related data
    # from this copy of the expression data
    del adata_qc.obsm, adata_qc.obsp, adata_qc.uns
    sys.stderr.write(f'Detecting highly variable genes with {hvg_flavor}...')
    if hvg_flavor != 'seurat_v3':
        sys.stderr.write(f'Lognormalizing aggregated count data')
        sc.pp.normalize_total(adata_qc)
        sc.pp.log1p(adata_qc)
    # 04d. Detect HVGs
    print(f'Finding top {n_hvgs_to_keep} highly variable genes...')
    sc.pp.highly_variable_genes(adata_qc,
                                n_top_genes=n_hvgs_to_keep,
                                flavor=hvg_flavor)
    adata_qc.obsm['hvgX'] = adata_qc.X[:, adata_qc.var['highly_variable']].copy()
    # 04e. Do population-wide expression level scaling per HVG
    print(f'Converting a copy of hvgX matrix to {working_dtype_float}...')
    hvg_cp = adata_qc.obsm['hvgX'].copy().astype(working_dtype_float)
    adata_qc.obsm['scaledX'], means, stds = scale_array__dtype_mod(hvg_cp,
            zero_center=True,
            max_value=hvg_scaling_cap,
            copy=True,
            return_mean_std=True,
            working_dtype_float=working_dtype_float)
    # 04f. Convert scaled count matrix to numpy array so it can be written
    # to hd5 format
    adata_qc.obsm['scaledX'] = np.asarray(adata_qc.obsm['scaledX'])
    # 04g. Because we know that the scaled counts should be much closer
    # to zero, we reduce this matrix to float16 to reduce memory loading.
    adata_qc.uns['hvg_scale_means'] = means.copy()
    adata_qc.uns['hvg_scale_stds'] = stds.copy()
    del means, stds
    adata_qc.obsm['X_pca'] = sc.tl.pca(adata_qc.obsm['scaledX'],
            random_state=numpy_random_state)
    sc.pp.neighbors(adata_qc,
                    use_rep='X_pca',
                    n_pcs=n_pcs,
                    method='umap',
                    metric='cosine',
                    n_neighbors=20,
                    random_state=numpy_random_state)
    # 04h. Perform knn-clustering-based manifold detection, UMAP
    # generation, and Leiden clustering
    sc.tl.umap(adata_qc,
            random_state=numpy_random_state)
    sc.tl.leiden(adata_qc,
        resolution=leiden_resolution,
        random_state=get_rand_seed(numpy_random_state))
    return adata_qc

# 05. Get parameters from configuration file for running the partial
# workflow, using network parameters that encourage more splitting
# of clusters.
sc_hyperparameters = cfg.get(f'sc_hyperparams__{pn}')
N_HVGS = sc_hyperparameters['n_hvgs']
max_val_hvg_scaling = sc_hyperparameters['hvg_expression_scale_cap_val']
if max_val_hvg_scaling == 'None':
    max_val_hvg_scaling = None
working_dtype_int = 'int32'
working_dtype_float = 'float32'
ACCUMULATOR_DTYPE = 'float32'
N_PCS = sc_hyperparameters['n_pcs']
LEIDEN_RES = sc_hyperparameters['leiden_res']

# 06. For each specified cluster, isolate and perform subclustering on the data,
# and then assign cluster labels in a new column
# 06a. Get list of cluster names
ct_label_col_vals = list(
        np.unique(
            adatas_human_qc.obs[ct_label_col].values.tolist()
            )
        )
if len(cts_to_subcluster) == 0:
    cts_to_subcluster = ct_label_col_vals
new_ct_list = []
# 06b. Work with a copy of adata.obs and then save it to 
# the anndata object later
ahc_obs_cp = adatas_human_qc.obs.copy()
for ct in ct_label_col_vals:
    idxs_ct = ahc_obs_cp[
            ahc_obs_cp[ct_label_col] == ct
            ].index.values.tolist()
    if ct in cts_to_subcluster:
        print('\n\n')
        print(f'Running subclustering on {ct}.')
        print(f'There are {len(idxs_ct)} nuclei in this cluster.')
        # 06c. Run a new round of HVG detection,
        # scaling, network mapping, and cluster detection
        # for expression data for the current cell type 
        ahc_ct_new_clust = qc_and_typical_workflow_partial(
                adata=adatas_human_qc[
                    idxs_ct,:
                    ].copy(),
                hvg_scaling_cap=max_val_hvg_scaling,
                n_hvgs_to_keep=N_HVGS,
                n_pcs=N_PCS,
                leiden_resolution=LEIDEN_RES,
                working_dtype_int=working_dtype_int,
                working_dtype_float=working_dtype_float,
                scaling_accumulator_dtype=ACCUMULATOR_DTYPE
                )
        # 06d. Pull the updated Leiden clusters for this
        # isolated cluster's expression data and use it 
        # to produce new names for subclusters
        ahc_ct_new_clust.obs[post_subclust_ct_label] = [
                f'{ct}__sc_{_}'
                for _ in
                ahc_ct_new_clust.obs['leiden'].copy()
                ]
        # 06e. Join the subclustering results for this cluster to
        # the copy of observational values for the full anndata
        # object
        if post_subclust_ct_label not in ahc_obs_cp.columns.values.tolist():
            ahc_obs_cp = ahc_obs_cp.join(
                    ahc_ct_new_clust.obs[post_subclust_ct_label].copy()
                    )
        else:
            for aidx in ahc_ct_new_clust.obs.index.values.tolist():
                ahc_obs_cp.at[aidx,post_subclust_ct_label] = ahc_ct_new_clust.obs.loc[
                        aidx][post_subclust_ct_label] 
        del ahc_ct_new_clust
        print(f'\n')
    else:
        # 06f. Assign the same cluster name to the new column
        print(f'Preserving current cluster label for {ct}.')
        df_vals_to_add = ahc_obs_cp.loc[idxs_ct].copy()
        if post_subclust_ct_label not in ahc_obs_cp.columns.values.tolist():
            df_vals_to_add.rename(
                    columns={ct_label_col:post_subclust_ct_label},
                    inplace=True
                    )
            # 06f.i. Set the column with values to add to have
            # a dtype 'object' to facilitate addition
            df_vals_to_add[post_subclust_ct_label] = df_vals_to_add[post_subclust_ct_label].astype('object')
            ahc_obs_cp = ahc_obs_cp.join(
                    df_vals_to_add[post_subclust_ct_label].copy()
                    )
        else:
            for aidx in idxs_ct:
                ahc_obs_cp.at[aidx,post_subclust_ct_label] = df_vals_to_add.loc[aidx
                        ][ct_label_col]
        del df_vals_to_add
        print(f'\n')
# 06g. Add the observation data with subcluster information back into
# the expression dataset
adatas_human_qc.obs = ahc_obs_cp
# 06h. Delete the working copy of the data to save space
del ahc_obs_cp

# 07. Add the working state of the random number generator to
# the updated AnnData object

# 08. Save the subclustered data using the specified file name tag
subclust_data_fn = f'{output_data_dir}/{output_data_fn}'
adatas_human_qc.write(subclust_data_fn)

# 06. Save the working state of the random number generator
numpy_random_state_dict = numpy_random_state.get_state(
        legacy=False
        )
numpy_random_state_output_full_fn = f'{output_random_state_dir}/{output_random_state_fn}'
with open(numpy_random_state_output_full_fn,'wb') as outfile:
    pickle.dump(numpy_random_state_dict,
            outfile)

sys.exit()

