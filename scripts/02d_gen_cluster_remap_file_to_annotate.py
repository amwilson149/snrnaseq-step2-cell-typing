import numpy as np
import pandas as pd
import anndata as ad
import datetime
import os
import sys
import argparse
import yaml

# This script generates a skeleton file with
# putative cell types per cluster to be used
# for manual validation of the specified
# cell typing pass.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Subcluster-to-refined-subcluster-cell-type map' + \
        ' file skeleton generation.')
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
# 01b.iii. Input epression data directory
input_data_dir = qc_root_dir + f'/' + cfg.get(f'input_data_dir__{pn}')
# 01b.iv. Input expression data file name
input_data_fn = cfg.get(f'input_data_fn__{pn}')
# 01b.v. Cell typing data directory
cell_typing_data_dir = output_parent_dir + f'/' + cfg.get('cell_typing_data_dir')

# 01b. Cluster label information
# 01b.i. Cluster label name to remap from (source)
clust_label = cfg.get(f'cluster_label__{pn}')
# 01b.ii. Cluster label name to remap to (target)
clust_label_out = cfg.get(f'output_cluster_label__{pn}')
# 01b.iii. Name of subcluster-to-refined-cell-type map file
cluster_mapping_fn = cfg.get(f'cluster_to_ct_map_fn__{pn}')

# 02. Import expression data
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 03. Pull specified column from observation data
ao = adatas_human_qc.obs.copy()
clust_list = []
if clust_label in ao.columns.values.tolist():
    clust_list = list(
            np.unique(
                ao[clust_label].values.tolist()
                )
            )
else:
    err_message = f'The column\n\t{clust_label}\nis not in the input expression data.' + \
            f'\nPlease check your input argmuments.'
    sys.exit(err_message)
del ao

# 04. Populate remap file
cluster_mapping_df = pd.DataFrame(
        data={
            f'{clust_label}':clust_list,
            f'{clust_label_out}':None,
            f'top_cell_types':None # to be filled in using the text file output from 02c_gen_est_cell_types__kfcv.py
            }
        )

# 05. Save map file skeleton
output_fn = f'{cell_typing_data_dir}/{cluster_mapping_fn}'
cluster_mapping_df.to_csv(output_fn,
        index=False)

sys.exit()
