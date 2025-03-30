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

# This script produces visualizations
# to facilitate handling of subclusters
# (i.e. determining how to isolate vs.
# recombine them).

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
parser = argparse.ArgumentParser(description='Visualization of marker genes for' + \
        ' subclustered expression data.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
parser.add_argument(f'--pass-number',type=str,required=True)
args = parser.parse_args()

# 02. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)

# 02a. Cell typing pass number.
# Convert this string to a 2-digit,
# zero padded integer if not in that
# form already
pn = args.pass_number
pn = pn.rjust(2,'0')

# 02b. Dataset information
# 02b.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 02b.ii. Input (original) dataset parent directory
dataset_parent_dir = qc_root_dir + f'/' + cfg.get('data_parent_dir')
# 02b.iii. Output (intermediate) dataset parent directory
output_parent_dir = qc_root_dir + f'/' + cfg.get('out_parent_dir')
# 02b.iii. Output (post-subclustering) expression data directory for current pass
input_data_dir = qc_root_dir + f'/' + cfg.get(f'output_data_dir__{pn}')
# 02b.iv. Output expression data file name
input_data_fn = cfg.get(f'output_data_fn__{pn}')
# 02b.v. Cell typing data directory
cell_typing_data_dir = output_parent_dir + f'/' + cfg.get('cell_typing_data_dir')

# 02b. Cluster label information
# 02b.i. Pre-subclustering cell type labels for current pass
input_cluster_lbl = cfg.get(f'cluster_label__{pn}')
# 02b.ii. Post-subclustering cell type labels for current pass
output_cluster_lbl = cfg.get(f'output_cluster_label__{pn}')

# Generate visualizations of subclustering results

# 03. Import expression data
print(f'Importing QCed, batch-corrected, preprocessed expression data...')
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 04. Make output directory
output_dir = f'{cell_typing_data_dir}/2D_UMAP_subclusters'
if not os.path.isdir(output_dir):
    os.system(f'mkdir -p {output_dir}')

# 05. Plot a 2D UMAP showing each input cluster's subclusters
input_clusts = list(
        np.unique(
            adatas_human_qc.obs[
                input_cluster_lbl
                ].copy().values.tolist()
            )
        )
# 05a. Iterate over input cell types
for i_clust in input_clusts:
    # 05b. Get the current cell types' subclusters
    ao = adatas_human_qc.obs.loc[
            adatas_human_qc.obs[input_cluster_lbl] == i_clust
            ].copy()
    subclusters_curr = list(
            np.unique(
                ao[output_cluster_lbl].values.tolist()
                )
            )
    # 05c. Delete observation data copy to save space
    del ao
    plt.figure()
    ax = plt.subplot(1,1,1)
    sc.pl.embedding(adatas_human_qc,
            basis='umap',
            color=output_cluster_lbl,
            palette=color_list_1,
            groups=subclusters_curr,
            show=False,
            ax=ax)
    i_clust_str_pieces = [_[0].upper()+_[1:].lower()
            if _.lower() not in ['odc','opc','smc']
            else _.upper()
            for _ in i_clust.split('_')
            ]
    i_clust_str = ' '.join(i_clust_str_pieces)
    title_str = f'Subclusters, {i_clust_str}'
    ax.set_title(title_str)
    plt.tight_layout()
    umap_fig_fn = f'{output_dir}/2D_UMAP__{i_clust}_subclusters.png'
    plt.savefig(umap_fig_fn,
            dpi=300)
    plt.close('all')

sys.exit()
