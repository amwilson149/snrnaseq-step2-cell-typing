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

# This script updates expression data
# with refined cluster labels following a
# subclustering pass.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Mapping of initial expression data clusters' + \
        ' to major cell types.')
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
dataset_parent_dir = cfg.get('root_dir') + f'/' + cfg.get('data_parent_dir')
# 01b.iii. Output (intermediate) dataset parent directory
output_parent_dir = cfg.get('root_dir') + f'/' + cfg.get('out_parent_dir')
# 01b.iv. Output (subclustered) expression data directory
input_data_dir = qc_root_dir + f'/' + cfg.get(f'output_data_dir__{pn}')
# 01b.v. Output (subclustered) expression data file name
input_data_fn = cfg.get(f'output_data_fn__{pn}')
# 01b.vi. Cell typing output data directory
cell_typing_data_dir = output_parent_dir + f'/' + cfg.get('cell_typing_data_dir')
# 01b.vii. Cell typing expression data output directory (in case different)
ct_output_data_dir = qc_root_dir + f'/' + cfg.get(f'output_data_dir__{pn}')
# 01b.vii. Cell typing output data file name
ct_output_data_fn = cfg.get(f'output_data_fn__{pn}')

# 01c. Cell typing information
# 01c.i. Annotated cluster-to-cell-type directory
scc_map_dir = qc_root_dir + f'/' + cfg.get(f'ann_cluster_to_ct_map_dir__{pn}')
# 01c.ii. Cluster-to-cell-type map file
scc_map_fn = cfg.get(f'ann_cluster_to_ct_map_fn__{pn}')
# 01c.iii. Post-subcluster label (pre-recombination)
clust_label = cfg.get(f'output_cluster_label__{pn}')
# 01c.iv. Specify updated cluster label (based on manually validated file format)
new_clust_label = f'new_{clust_label}'

# Apply refined subcluster labels to expression data

# 02. Import expression data
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)
print(adatas_human_qc)
print('\n')

# 03. Make a copy of the observation
# annotations to update
ao = adatas_human_qc.obs.copy()

# 04. Import cluster remapping file
cluster_map_file_full_fn = f'{scc_map_dir}/{scc_map_fn}'
clust_mapping_df = pd.read_csv(cluster_map_file_full_fn)

# 05. Remap subcluster names in expression data,
# making changes in-place to the current pass
# subcluster column
# 05a. First, make sure the cluster column is
# an object data type to ensure that any categories
# removed during cluster recombination are also
# removed from the category list
ao[clust_label] = ao[clust_label].astype('object')
for clust_curr,refined_clust_curr in zip(
        clust_mapping_df[clust_label].values.tolist(),
        clust_mapping_df[new_clust_label].values.tolist()
        ):
    aidxs = ao.loc[
            ao[clust_label] == clust_curr
            ].index.values.tolist()
    for aidx in aidxs:
        ao.at[aidx,clust_label] = refined_clust_curr
ao[clust_label] = ao[clust_label].astype('category')

print('updated names:')
print(ao.head(10)[clust_label])

# 06. Add the observation data with the updated
# subcluster column to the expression data
adatas_human_qc.obs = ao.copy()

# 07. Delete the working copy of the observation
# data to clear memory
del ao

# 07. Save the updated expression data
output_fn = f'{ct_output_data_dir}/{ct_output_data_fn}'
adatas_human_qc.write(output_fn)

sys.exit()


