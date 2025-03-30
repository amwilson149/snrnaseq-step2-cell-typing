import numpy as np
import pandas as pd
import anndata as ad
import os
import sys
import argparse
import yaml

# This script assigns putative cell types
# from the current cell typing pass to
# expression data.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Mapping of initial expression data clusters' + \
        ' to major cell types.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
parser.add_argument(f'--pass-number',type=str,required=True)
parser.add_argument(f'--finalize',type=bool,required=False)
args = parser.parse_args()

# 00a. Determine whether this run will be finalizing cell types
to_finalize=False
if args.finalize:
    to_finalize=args.finalize

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
# 01b.iv. Input expression data directory
input_data_dir = qc_root_dir + f'/' + cfg.get(f'input_data_dir__{pn}')
# 01b.v. Input expression data file name
input_data_fn = cfg.get(f'input_data_fn__{pn}')
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
# 01c.iii. Input cluster label
input_ct_lbl = cfg.get(f'cluster_label__{pn}')
# 01c.iv. Output cluster label
output_ct_lbl = cfg.get(f'output_cluster_label__{pn}')

# 02. Import expression data
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 03. Import annotated cluster-to-cell-type file
scc_map_full_fn = f'{scc_map_dir}/{scc_map_fn}'
scc_map_df = pd.read_csv(scc_map_full_fn)

# 04. Add new cell type labels to expression data
ahc_obs_cp = adatas_human_qc.obs.copy()
if output_ct_lbl in ahc_obs_cp.columns:
    ahc_obs_cp[output_ct_lbl] = ahc_obs_cp[output_ct_lbl].astype('object')
for in_ctl, out_ctl in zip(
        scc_map_df[input_ct_lbl].values.tolist(),
        scc_map_df[output_ct_lbl].values.tolist()
        ):
    # 04a. Find barcodes in expression data with
    # input cluster label
    ahc_clust_idxs = ahc_obs_cp.loc[
            ahc_obs_cp[input_ct_lbl]==f'{in_ctl}'
            ].index.values.tolist()
    # 04b. Assign new output cluster label to those barcodes
    # 04b.i. Reformat output cluster label
    out_ctl_reform = '_'.join(
            [
                _.strip().lower()
                for _ in f'{out_ctl}'.split(' ')
                ]
            )
    ahc_obs_cp.loc[ahc_clust_idxs,output_ct_lbl] = out_ctl_reform

# 04c. Add updated annotations to expression data
adatas_human_qc.obs = ahc_obs_cp.copy()
del ahc_obs_cp

# 05. If it is specified that this update should lead to finalized
# cell types, run extra steps to retain only barcodes with finalized
# cell types, and output the finalized dataset with the appropriate
# finalized cell type name
if to_finalize==True:
    # 05a. Filter to keep only finalized cell types
    adatas_human_qc = adatas_human_qc[
            (
                (~adatas_human_qc.obs[output_ct_lbl].isna())
                &
                (adatas_human_qc.obs[output_ct_lbl]!='')
                )
            ]
    finalized_ct_label = cfg.get('finalized_ct_label')
    adatas_human_qc.obs.rename(
            columns={output_ct_lbl:finalized_ct_label},
            inplace=True
            )
    # 05b. Pull final data output file name
    ct_output_data_fn = cfg.get('output_data_fn__final')

# 05. Save updated expression data
ahc_ct_fn = f'{ct_output_data_dir}/{ct_output_data_fn}'
adatas_human_qc.write(ahc_ct_fn)

sys.exit()

