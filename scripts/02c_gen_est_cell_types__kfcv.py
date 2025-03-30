import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import sys
import json
import argparse
import yaml

# This script takes k-fold cross-validated output
# from a cell typing pass, and uses it to
# generate a list of putative cell types
# per cluster for validation, along with
# visualization of putative cell type
# log-likelihoods.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Marker-gene-likelihood-based cell typing' + \
        'for initial expression data clusters.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
parser.add_argument(f'--pass-number',type=str,required=True)
args = parser.parse_args()
# 00a. Get color names for plotting
color_names = mcolors.cnames
color_name_keys = [_ for _ in color_names.keys()]

# 01. Get analysis parameters
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Cell typing pass number.
# Convert this string to a 2-digit,
# zero padded integer if not in that
# form already
pn = args.pass_number
pn = pn.rjust(2,'0')

# 02. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01a.ii. Input (original) dataset parent directory
dataset_parent_dir = qc_root_dir + f'/' + cfg.get('data_parent_dir')
# 01a.iii. Output (intermediate) dataset parent directory
output_parent_dir = qc_root_dir + f'/' + cfg.get('out_parent_dir')
# 01a.iii. Analysis expression data directory
input_data_dir = qc_root_dir + f'/' + cfg.get(f'input_data_dir__{pn}')
# 01a.iv. Analysis expression data file name
input_data_fn = cfg.get(f'input_data_fn__{pn}')
# 01a.v. Cell typing data directory
cell_typing_data_dir = output_parent_dir + f'/' + cfg.get('cell_typing_data_dir')
# 01a.vi. Marker gene database directory
marker_db_dir = qc_root_dir + f'/' + cfg.get('mg_db_dir')
marker_db_fn = cfg.get('mg_db_fn')

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

# Assign putative cell types by log likelihood
# to each cluster

# 02. Get cell typing results DataFrame
# (lowest powered log2 fold change 
# overexpression threshold)
put_ct_fn = f'{cell_typing_data_dir}'
put_ct_fn_name_start = f'putative_cell_types_per_clust__{cluster_label}'
put_ct_fn_name_end = f'__kfcv_{n_folds}'
ct_files = os.listdir(cell_typing_data_dir)
overexpression_threshold = None
# 02a. Add file name with kfcv tag to file path
for ct_file in ct_files:
    if (
            (put_ct_fn_name_start in ct_file)
            and
            (put_ct_fn_name_end in ct_file)
            ):
        put_ct_fn += f'/{ct_file}'
        # 02b. Pull the file's
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
# 02a. Sort DataFrame by index to print putative cell type
# information in order of ascending cluster label
cell_type_all_df.sort_index(inplace=True)

# 03. For each cluster in the DataFrame, convert the log likelihoods
# into probabilities and inspect the dropoff in cell type probability.
# Note: since putative cell types in the database are each unique, 
# the cell type probability is calculated by normalizing the likelihoods.
color_idx = 10
color_idx_increment = 2
plt.style.use('dark_background')
fig,axs = plt.subplots(nrows=1,
        ncols=3)
threshold_pct = 95
threshold_frac = threshold_pct*1.0/100
slope_threshold = -0.005
threshold_lk_rel_to_max = 0.0001
pct_dict_no_thresh = {}
pct_dict_toppct_thresh = {}
pct_dict_slope_thresh = {}
pct_dict_lkrm_thresh = {}
for idx in cell_type_all_df.index:
    # 03a. Convert log likelihoods into likelihoods,
    # then into probabilities for visualization
    # Note: log likelihoods were computed by adding
    # by adding log2 fold changes are thus in log base 2
    log_lks = [float(_) for _ in 
            cell_type_all_df.loc[idx][
                'log_likelihoods'
                ].split('[')[-1].split(']')[0].split(',')
            ]
    lks = [2**(_) for _ in log_lks]
    probs = [_/np.sum(lks) for _ in lks]
    cell_type_all_df.at[idx,'probs'] = str(probs)
    # 03b. Sort probabilities and their corresponding
    # putative cell types in decreasing order
    idxs_prob_sort = np.argsort(probs)[::-1]
    probs_decreasing = [probs[_]
            for _ in idxs_prob_sort]
    cts = [_.strip() for _ in
            cell_type_all_df.loc[idx][
                'cell_types'
                ].replace(
                    '\'',''
                    ).split(
                        '['
                        )[-1].split(
                            ']'
                            )[0].split(
                                ','
                                )
                            ]
    cts_decreasing = [cts[_]
            for _ in idxs_prob_sort]
    # 03c. Plot probabilities in decreasing order
    axs[0].plot(np.arange(len(probs)),
            probs_decreasing,
            color=color_name_keys[color_idx],
            alpha=0.7,
            label=f'{idx}'
            )
    # 03d. Plot the CDF of the probabilities in decreasing order
    cdf = [np.sum(probs_decreasing[:(i+1)]) 
            for i in 
            range(len(probs_decreasing))]
    axs[1].plot(np.arange(len(cdf)),
            cdf,
            color=color_name_keys[color_idx],
            alpha=0.7,
            label=f'{idx}'
            )
    # 03e. Plot the approximate centered finite-difference-based slopes of the
    # probabilities in decreasing order
    cfd_slopes = [(probs_decreasing[_+1]-probs_decreasing[_-1])/2.0
            for _ in range(1,len(probs_decreasing)-1)]
    axs[2].plot(np.arange(len(cfd_slopes)),
            cfd_slopes,
            color=color_name_keys[color_idx],
            alpha=0.7,
            label=f'{idx}'
            )
    color_idx += color_idx_increment
    # 03f. Explore thresholds for more vs. less likely
    # putative cell types
    # 03f.i. Find the first point in the CDF that is above pre-specified
    # cutoff thresholds (see above). This and greater-likelihood putative
    # cell types would be included in thresholded putative cell type lists
    idx_cdf_thresh = next(
            (
                i+1 for i,_ in enumerate(cdf)
                if (
                    (_ > threshold_frac)
                    )
                ),
            None
            )
    # 03f.ii. Find the first point at which two consecutive slopes
    # are below the defined threshold (meaning these and less-likely
    # cell types no longer have very different likelihoods)
    idx_slope_thresh = next(
            (
                i for i,_ in enumerate(cfd_slopes[:-1])
                if (
                    (_ > slope_threshold)
                    and
                    (cfd_slopes[i+1] > slope_threshold)
                    )
                ),
            None
            )
    # 03f.ii. Compute a threshold based on likelihood relative to the
    # most likely option, in which anything less than the threshold
    # fraction of the maximum value is removed.
    rel_lks_decreasing = [_/np.max(probs_decreasing)
            for _ in probs_decreasing]
    idx_lk_rel_max_thresh = next(
            (
                i for i,_ in enumerate(rel_lks_decreasing)
                if _ < threshold_lk_rel_to_max
                ),
            None
            )

    probs_decreasing_to_thresh = probs_decreasing[:idx_cdf_thresh]
    cts_decreasing_to_thresh = cts_decreasing[:idx_cdf_thresh]
    probs_descreasing_to_slope_thresh = probs_decreasing[:idx_slope_thresh]
    cts_decreasing_to_slope_thresh = cts_decreasing[:idx_slope_thresh]#
    probs_decreasing_to_lkrm_thresh = probs_decreasing[:idx_lk_rel_max_thresh]
    cts_decreasing_to_lkrm_thresh = cts_decreasing[:idx_lk_rel_max_thresh]
    
    # 03g. Add putative cell types and their probabilites
    # for the current cluster to the appropriate dictionaries
    pct_dict_no_thresh[idx] = [
            [cts_decreasing[_],probs_decreasing[_]]
            for _ in range(len(cts_decreasing))
            ]
    pct_dict_toppct_thresh[idx] = [
            [cts_decreasing_to_thresh[_],probs_decreasing_to_thresh[_]]
            for _ in range(len(cts_decreasing_to_thresh))
            ]
    pct_dict_slope_thresh[idx] = [
            [cts_decreasing_to_slope_thresh[_],probs_descreasing_to_slope_thresh[_]]
            for _ in range(len(cts_decreasing_to_slope_thresh))
            ]
    pct_dict_lkrm_thresh[idx] = [
            [cts_decreasing_to_lkrm_thresh[_],probs_decreasing_to_lkrm_thresh[_]]
            for _ in range(len(cts_decreasing_to_lkrm_thresh))
            ]

# 03h. Plot a guideline on the CDF plot
# to show the top-pct threshold line
axs[0].set_ylim(bottom=0,
        top=1.1)
axs[1].set_ylim(bottom=0,
        top=1.1)
axs[0].set_yticks(np.arange(0,1.1,0.2))
axs[0].set_yticklabels([f'{_:0.1}'
    for _ in
    np.arange(0,1.1,0.2)],
    fontsize=6,
    rotation=90.0)
axs[1].set_yticks(np.arange(0,1.1,0.2))
axs[1].set_yticklabels([f'{_:0.1}'
    for _ in 
    np.arange(0,1.1,0.2)],
    fontsize=6,
    rotation=90.0)
axs[2].set_yticks(np.arange(-0.5,0.2,0.1))
axs[2].set_yticklabels([f'{_:0.1}'
    for _ in np.arange(-0.5,0.2,0.1)],
    fontsize=6,
    rotation=90.0)
axs[1].plot([0,len(cell_type_all_df)],
        [threshold_frac]*2,
        '--',
        color='w')
axs[2].plot([0,len(cell_type_all_df)],
        [slope_threshold]*2,
        '--',
        color='w')
axs[0].set_xlabel('Cell type ordinal index\n(unique per group; descending likelihood)',
        fontsize=6)
axs[1].set_xlabel('Cell type ordinal index',
        fontsize=6)
axs[2].set_xlabel('Cell type ordinal index',
        fontsize=6)
axs[0].set_ylabel('Probability',
        fontsize=6)
axs[0].legend(fontsize=4)
axs[1].legend(fontsize=4)
axs[0].set_title('Cell type probabilities',
        fontsize=8)
axs[1].set_title('CDF',
        fontsize=8)
axs[2].set_title('Slopes',
        fontsize=8)
prob_plot_fn = f'{cell_typing_data_dir}/put_ct_probs__{cluster_label}_pass_{pn}__{overexpression_threshold}_l2fc__{n_folds}_kfcv.png'
plt.savefig(prob_plot_fn,
        dpi=300)
plt.close()

# 04. Write the sorted putative cell types per cluster to file,
# both with and without CDF-based thresholding
pct_no_thresh_summary_fn = f'{cell_typing_data_dir}/putative_cell_types_per_clust__no_threshold__{cluster_label}_pass_{pn}__{overexpression_threshold}_l2fc__{n_folds}_kfcv.txt'
with open(pct_no_thresh_summary_fn,'w') as f:
    for key,values in pct_dict_no_thresh.items():
        f.write(f'Cluster {key}:\n')
        for ct in values:
            f.write(f'\t{ct[0]}:\t\t{ct[1]:0.04}\n')
        f.write(f'\n')

# 04a. Writing of files useful for inspection only is commented out to avoid confusion
#pct_cdf_thresh_summary_fn = f'{cell_typing_data_dir}/putative_cell_types_per_clust__top_{threshold_pct}_pct_probs__{cluster_label}_pass_{pn}__{overexpression_threshold}_l2fc__{n_folds}_kfcv__INSPECTION_ONLY.txt'
#with open(pct_cdf_thresh_summary_fn,'w') as f:
#    for key,values in pct_dict_toppct_thresh.items():
#        f.write(f'Cluster {key}:\n')
#        for ct in values:
#            f.write(f'\t{ct[0]}:\t\t{ct[1]:0.04}\n')
#        f.write(f'\n')
#
#pct_slope_thresh_summary_fn = f'{cell_typing_data_dir}/putative_cell_types_per_clust__slope_gt_n{abs(slope_threshold)}__{cluster_label}_pass_{pn}__{overexpression_threshold}_l2fc__{n_folds}_kfcv__INSPECTION_ONLY.txt'
#with open(pct_slope_thresh_summary_fn,'w') as f:
#    for key,values in pct_dict_slope_thresh.items():
#        f.write(f'Cluster {key}:\n')
#        for ct in values:
#            f.write(f'\t{ct[0]}:\t\t{ct[1]:0.04}\n')
#        f.write(f'\n')
#
#pct_lkrm_thresh_summary_fn = f'{cell_typing_data_dir}/putative_cell_types_per_clust__lk_rel_max_lt{threshold_lk_rel_to_max}__{cluster_label}_pass_{pn}__{overexpression_threshold}_l2fc__{n_folds}_kfcv__INSPECTION_ONLY.txt'
#with open(pct_lkrm_thresh_summary_fn,'w') as f:
#    for key,values in pct_dict_lkrm_thresh.items():
#        f.write(f'Cluster {key}:\n')
#        for ct in values:
#            f.write(f'\t{ct[0]}:\t\t{ct[1]:0.04}\n')
#        f.write(f'\n')


sys.exit()

    






































