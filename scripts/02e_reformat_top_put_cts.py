import numpy as np
import pandas as pd
import os
import sys
import argparse
import yaml

# This script takes a cluster-to-cell-type map
# with a column labeled 'top_cell_types' and,
# using the cell-type-to-collapsed-cell-type
# map, converts these putative cell types to
# collapsed cell types, puts them in the column
# 'top_collapsed_cell_types', and saves the
# cluster-to-cell-type map again.

# Script setup

# 00. Create an argparse object for easier reading
# of input arguments
parser = argparse.ArgumentParser(description='Collapsing of top cell types to standardized' + \
        'major cell type values in initial-expression-data-cluster-to-cell-type map.')
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
# 01b.ii. Output (intermediate) dataset parent directory
output_parent_dir = qc_root_dir + f'/' + cfg.get('out_parent_dir')
# 01b.v. Cell typing data directory
cell_typing_data_dir = output_parent_dir + f'/' + cfg.get('cell_typing_data_dir')

# 01b. Cluster label information
# 01b.i. Cluster label name to remap from (source)
clust_label = cfg.get(f'cluster_label__{pn}')
# 01b.ii. Cluster label name to remap to (target)
clust_label_out = cfg.get(f'output_cluster_label__{pn}')
# 01b.iii. Annotated cluster-to-cell-type map directory
ann_clust_ct_map_dir = qc_root_dir + f'/' + cfg.get(f'ann_cluster_to_ct_map_dir__{pn}')
# 01b.iii. Annotated cluster-to-cell-type map file name
ann_clust_map_fn = cfg.get(f'ann_cluster_to_ct_map_fn__{pn}')
# 01b.iv. Cell-type-to-collapsed-cell-type directory
ct_coll_ct_dir = cfg.get('mg_db_dir')
# 01b.v. Cell-type-to-collapsed-cell-type file name
ct_coll_ct_fn = cfg.get('ct_col_ct_fn')

# 02. Read in cluster-to-cell-type map
scc_map_full_fn = f'{ann_clust_ct_map_dir}/{ann_clust_map_fn}'
scc_map_df = pd.read_csv(scc_map_full_fn)

# 03. Read in the annotated cell-type-to-collapsed-cell-type-map
ct_coll_ct_full_fn = f'{ct_coll_ct_dir}/{ct_coll_ct_fn}' 
ct_coll_ct_df = pd.read_csv(ct_coll_ct_full_fn,index_col=0)
# 03a. Make sure all columns (which should correspond to cluster names)
# have an 'object' dtype
for col_curr in ct_coll_ct_df.columns.values.tolist():
    ct_coll_ct_df[col_curr] = ct_coll_ct_df[col_curr].astype('object')

# 04. Translate each cluster's list of top putative cell
# types into a list of unique, collapsed cell types
for idx in scc_map_df.index:
    print(f'Processing cluster {idx}....')
    row = scc_map_df.loc[idx]
    print(row)
    top_cts = [_.strip()
            for _ in
            row['top_cell_types'].split(',')]
    coll_cts = []
    for tct in top_cts:
        if tct != '':
            coll_ct_curr = ct_coll_ct_df.loc[
                        ct_coll_ct_df['cell_type'] == tct
                        ]['collapsed_cell_type'].values.tolist()[0]
            # 04a. Only add if the collapsed cell type is not already
            # in the list; this approach results in a unique list
            # of collapsed cell types with preserved order
            if coll_ct_curr not in coll_cts:
                coll_cts.append(coll_ct_curr)
    scc_map_df.at[idx,'top_collapsed_cell_types'] = ','.join(coll_cts)

scc_map_df.to_csv(scc_map_full_fn,
        index=False)

sys.exit()
