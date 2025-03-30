import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import os
import sys
import argparse
import yaml

# This script produces a plot of the
# clusters for the current cell typing
# pass that were assigned putative
# cell types, as an initial inpsection
# of this process.

# 00. Get color names for plotting
color_names = mcolors.cnames
# 00a. get a sequence of visible, distinct colors
# to plot on UMAPs
color_list_1 = ['forestgreen','darkturquoise','dodgerblue','orangered',
        'lightpink','deeppink','yellowgreen','darkblue',
        'crimson','goldenrod','darkorange','red',
        'slategrey','slateblue','saddlebrown','tan',
        'teal','darkmagenta','midnightblue','black',
        'limegreen','mediumpurple','gold','magenta',
        'mediumblue','darkgreen','darkred','olive',
        'lightsalmon','rebeccapurple','lightskyblue','darkgray',
        'lightcoral']

# Script setup

# 01. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Visualization of cell type likelihood for initial' + \
        ' expression data clusters.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
parser.add_argument(f'--pass-number',type=str,required=True)
args = parser.parse_args()

# 01a. Get the cell typing pass number
pn = args.pass_number
pn = pn.rjust(2,'0')

# 02. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01a.ii. Original (input) dataset parent directory
dataset_parent_dir = qc_root_dir + f'/' + cfg.get('data_parent_dir')
# 01a.iii. Output (intermediate) dataset parent directory
output_parent_dir = qc_root_dir + f'/' + cfg.get('out_parent_dir')
# 01a.iii. Input expression data directory
input_data_dir = qc_root_dir + f'/' + cfg.get(f'input_data_dir__{pn}')
# 01a.iv. Input expression data file name
input_data_fn = cfg.get(f'input_data_fn__{pn}')
# 01a.v. Cell typing data directory
cell_typing_data_dir = output_parent_dir + f'/' + cfg.get('cell_typing_data_dir')

# 01b. Cell typing parameters
# 01b.i. Cluster label used for specified cell typing pass
cluster_label = cfg.get(f'cluster_label__{pn}')
# 01b.ii. Directory with file containing n folds
# used for distinguishing HVG detection
n_fold_log_dir = qc_root_dir + f'/' + cfg.get('kfcv_n_fold_log_dir')
n_fold_log_fn = cfg.get(f'n_folds_used_fn__02_ct__{pn}')

# 01c. Read the number of folds used in from file
n_folds_full_fn = f'{n_fold_log_dir}/{n_fold_log_fn}'
n_folds = None
with open(n_folds_full_fn,'r') as f:
    lines = f.readlines()
    if len(lines) > 0:
        n_folds = int(lines[0])

# Inspect clusters with putative cell types 

# 02. Import expression data
print(f'Importing QCed, batch-corrected, preprocessed expression data...')
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 03. Read in the putative cell type file for the current grouping
put_ct_fn = f'{cell_typing_data_dir}'
put_ct_fn_name_start = f'putative_cell_types_per_clust__{cluster_label}'
put_ct_fn_name_end = f'__kfcv_{n_folds}'
ct_files = os.listdir(cell_typing_data_dir)
overexpression_threshold = None
# 03a. Add file name with kfcv tag to file path
for ct_file in ct_files:
    if (
            (put_ct_fn_name_start in ct_file)
            and
            (put_ct_fn_name_end in ct_file)
            ):
        put_ct_fn += f'/{ct_file}'
        # 03b. Pull this file's
        # overexpression threshold
        overexpression_threshold = float(
                ct_file[:
                    ct_file.find('l2fc')
                    ].split('__')[-1].split('_')[0]
                )
cell_type_all_df = pd.read_csv(
        put_ct_fn,
        index_col=0
        )
    
# 04. Plot a 2D UMAP of expression data showing
# which clusters have putative cell type likelihood vectors
# 04a. Add a label to the anndata object that is equivalent
# to the cluster label, but with value 0 for clusters
# not assigned putative cell types
clust_w_ct_info = []
for idx in cell_type_all_df.index:
    row = cell_type_all_df.loc[idx]
    llks = row['log_likelihoods']
    n_llks = len(llks.split(','))
    if n_llks > 1:
        clust_w_ct_info.append(str(idx))
print(f'Clusters with putative cell types: {clust_w_ct_info}')
# 04b. Plot the UMAP with the cell type masked labels
fig, ax = plt.subplots()
sc.pl.embedding(adatas_human_qc,
        basis='umap',
        color=cluster_label,
        palette=color_list_1,
        groups=clust_w_ct_info,
        show=False,
        ax=ax)
plt.tight_layout()
umap_fig_fn = f'{cell_typing_data_dir}/2D_UMAP_clusters_w_cell_type_info' + \
        f'__{cluster_label}_pass_{pn}' + \
        f'__{overexpression_threshold}_l2fc_min.png'
plt.savefig(umap_fig_fn,
        dpi=300)
plt.close('all')

sys.exit()


