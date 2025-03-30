import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import scipy.stats as st
import matplotlib.pyplot as plt
import pickle
import os
import sys
import argparse
import yaml
from utils.power_analysis import get_sample_size
from utils.get_rand_seed import get_rand_seed

# This script performs reads in expresion data
# and assigns to each input cluster a list
# of likelihoods, one for each of the
# putative cell types listed in the input
# marker gene database.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Marker gene detection and cell type likelihood' + \
        ' estimation for initial expression data clusters.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
parser.add_argument(f'--pass-number',type=str,required=True)
args = parser.parse_args()

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

# 01b. Dataset information
# 01b.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01b.ii. Input (original) data parent directory
dataset_parent_dir = qc_root_dir + f'/' + cfg.get('data_parent_dir')
# 01b.iii. Output (processed) data parent directory
output_parent_dir = qc_root_dir + f'/' + cfg.get('out_parent_dir')
# 01b.iv. Validated expression data directory
input_data_dir = qc_root_dir + f'/' + cfg.get(f'input_data_dir__{pn}')
# 01b.v. Input expression data file name
input_data_fn = cfg.get(f'input_data_fn__{pn}')
# 01b.vi. Cell typing data output directory
cell_typing_data_dir = output_parent_dir + f'/' + cfg.get('cell_typing_data_dir')
if not os.path.isdir(cell_typing_data_dir):
    os.system(f'mkdir -p {cell_typing_data_dir}')
# 01b.v. Marker gene database directory
marker_db_dir = qc_root_dir + f'/' + cfg.get('mg_db_dir')
# 01b.vi. Marker gene database file name
marker_db_fn = cfg.get('mg_db_fn')
# 01b.vii. Root file name for list
# of each cluster's distinguishing HVGs
clust_dist_HVG_fn_tag = cfg.get(f'dist_HVG_list__{pn}')

# 01c. Parameters to use for cluster distinguishing
# HVG identification in the current cell typing pass
# 01c.i. Cluster label
cluster_lbl = cfg.get(f'cluster_label__{pn}')
# 01c.ii. Adjusted p-value
padj_dist_HVGs = cfg.get(f'padj_distinguishing_HVGs')
# 01c.iii. Whether to use k+1 for a more
# conservative k-fold cross-validation
go_conservative_kfcv = cfg.get('do_more_conservative_k_fold_cv')

# 01d. Parameters for cell type log-likelihood estimation
# 01d.i. Whether to use negative marker genes
use_neg_marker_genes = cfg.get('use_neg_marker_genes_in_llks')

# 01e. Random state file setup
# 01e.i. Input random state directory
input_random_state_dir = qc_root_dir + f'/' + cfg.get(f'input_rs_dir__{pn}')
# 01e.ii. Input random state file name
input_random_state_fn = cfg.get(f'input_rs_fn__{pn}')
# 01e.ii.A. If no random state is specified, read in a random seed
input_random_seed = None
if input_random_state_fn is None:
    input_random_seed = cfg.get(f'input_rs__{pn}')
# 01e.iii. Output random state file directory
output_random_state_dir = qc_root_dir + f'/' + cfg.get(f'output_rs_dir__{pn}')
if not os.path.isdir(output_random_state_dir):
    os.system(f'mkdir -p {output_random_state_dir}')
# 01e.iv. Output random state file name
output_random_state_fn = cfg.get(f'output_rs_fn__{pn}')

# 01f. k-fold cross-validation file setup
# 01f.i. Directory for log files containing the number
# of folds k used in the current pass
kfcv_log_dir = qc_root_dir + f'/' + cfg.get('kfcv_n_fold_log_dir')
if not os.path.isdir(kfcv_log_dir):
    os.system(f'mkdir -p {kfcv_log_dir}')
n_fold_output_log_fn = cfg.get(f'n_folds_used_fn__02_ct__{pn}')

# 02. Set up the analysis random state from the
# most recent checkpoint
numpy_random_state = None
if input_random_state_fn is not None:
    numpy_random_state = np.random.RandomState()
    working_np_state = None
    full_random_state_input_fn = f'{input_random_state_dir}/{input_random_state_fn}'
    with open(full_random_state_input_fn,'rb') as f:
        working_np_state = pickle.load(f)
    numpy_random_state.set_state(working_np_state)
    del working_np_state
else:
    numpy_random_state = np.random.RandomState(
            seed=input_random_seed
            )

# Identify cluster distinguishing marker genes
# and compute putative-cell-type likelihoods

# 03. Import marker gene database
marker_db_input_fn = f'{marker_db_dir}/{marker_db_fn}'
marker_db = pd.read_csv(
        marker_db_input_fn,
        index_col=0)

# 04. Build dictionary mapping marker genes (as keys)
# to cell types
marker_dict = {}
negative_marker_dict = {}
# 04a. Define the reference genome prefix so that dictionary
# marker gene names correspond with dataset gene names
genome_prefix = 'GRCh38_______________'
marker_gene_names_list = []
negative_marker_gene_names_list = []
for idx in marker_db.index:
    row = marker_db.loc[idx]
    cell_type = row['cell_type']
    marker_genes = row['marker_genes'].split(',')
    negative_marker_genes_raw = row['negative_markers']
    negative_marker_genes = []
    # 04b. If there are negative markers (such that the
    # DataFrame entry is not nan), add them to the list
    if negative_marker_genes_raw==negative_marker_genes_raw:
        negative_marker_genes = negative_marker_genes_raw.split(',')
    marker_gene_names_list.extend(marker_genes)
    negative_marker_gene_names_list.extend(negative_marker_genes)
    for mg in marker_genes:
        mg_str = f'{genome_prefix}{mg}'
        if mg_str in marker_dict.keys():
            marker_dict[mg_str].extend([cell_type])
        else:
            marker_dict[mg_str] = [cell_type]
    for nm in negative_marker_genes:
        nm_str = f'{genome_prefix}{nm}'
        if nm_str in negative_marker_dict.keys():
            negative_marker_dict[nm_str].extend([cell_type])
        else:
            negative_marker_dict[nm_str] = [cell_type]
del marker_gene_names_list, negative_marker_gene_names_list

# 05. Import expression data
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 06. Define helper functions to run distinguishing
# HVG detection with k-fold cross-validation
# 06a. Function to perform distinguishing HVG
# detection for one k-fold data subset
def get_dist_HVGs(adata,
        grouping='leiden',
        log2_min: float = None,
        pval_max: float = None):
    if 'log1p' in adatas_human_qc.uns.keys():
        del adatas_human_qc.uns['log1p']
    # 06a.i. Normalize expression data (in-place)
    print('Computing log1p values...')
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # 06a.ii. Rank expression by cluster
    print('Running rank_genes_groups...')
    sc.tl.rank_genes_groups(adata,groupby=grouping,method='wilcoxon')
    # 06a.iii. Pull significantly up- or down-regulated genes for the
    # current cluster vs. the rest (significant padj,
    # log2 fold change magnitude > specified log2_min)
    print('Pulling DataFrame with gene ranks...')
    genes_df = sc.get.rank_genes_groups_df(adata,
            group=adata.obs[grouping].unique(),
            log2fc_min=log2_min,
            pval_cutoff=pval_max)
    # 06a.iv. Add gene counts to DataFrame for context
    print('Adding count data to DataFrame...')
    counts_df = (adata.var[['n_counts_all_cells']]
            .reset_index()
            .rename(columns={'gene_name': 'names'}))
    df = pd.merge(left=genes_df,
            right=counts_df,
            how='left',
            on='names')
    return df

# 06b. Wrapper function for performing
# k-fold cross-validated distinguishing
# HVG detection
def k_fold_get_dist_HVGs(adata,
        grouping = 'leiden',
        log2_min: float = None,
        pval_max: float = None,
        fold_col: str = None,
        k: int = 5,
        testing = False):
    ahc_df_all = pd.DataFrame()
    # 06b.i. Generate k-fold
    # data subsets and run
    # HVG detection
    for fold in range(k):
        # 06b.ii. Generate a data subset
        # consisting of all but the current
        # fold of barcodes
        adata_for_mg = adata[
                adata.obs[fold_col]!=fold
                ].copy()
        # 06b.iii. Run distinguishing HVG detection
        ahc_df_curr = get_dist_HVGs(adata_for_mg,
                grouping=grouping,
                log2_min=log2_min,
                pval_max=pval_max)
        # 06b.iv. Add the resulting DataFrame to
        # the all-subset DataFrame
        if len(ahc_df_all) == 0:
            ahc_df_all = ahc_df_curr.copy()
        else:
            ahc_df_all = pd.concat(
                    [ahc_df_all,
                        ahc_df_curr.copy()],
                    ignore_index=True)

    ahc_df_final = pd.DataFrame()
    if testing:
        per_group_mg_redundancy = {}
        for g,g_df in ahc_df_all.groupby('group'):
            mg_names, mg_counts = np.unique(
                    g_df['names'].values.tolist(),
                    return_counts = True)
            per_group_mg_redundancy[g] = [_ for _ in mg_counts]
        return ahc_df_all, per_group_mg_redundancy
    else:
        # 06b.v. Include only distinguishing HVGs present
        # for all data subsets
        for g,g_df in ahc_df_all.groupby('group'):
            mg_names, mg_counts = np.unique(
                    g_df['names'].values.tolist(),
                    return_counts = True)
            # 06b.v.A. Pull HVGs appearing in all data subsets
            mgs_to_keep = [_ for i,_ in enumerate(mg_names)
                    if mg_counts[i] == k]
            ahc_df_f_curr = g_df.loc[
                    g_df['names'].isin(mgs_to_keep)]
            # 06b.v.B. Compute the fold-averages of all scores, log2 fold
            # changes, and p-values for distinguishing
            # HVGs appearing across all data subsets
            ahc_df_avg = pd.DataFrame()
            df_avg_idx = 0
            for mg,mg_df in ahc_df_f_curr.groupby('names'):
                ahc_df_avg.at[df_avg_idx,'group'] = g
                ahc_df_avg.at[df_avg_idx,'names'] = mg
                ahc_df_avg.at[df_avg_idx,'n_counts_all_cells'] = list(
                        np.unique(
                            mg_df['n_counts_all_cells'].values.tolist()
                            )
                        )[0]
                ahc_df_avg.at[df_avg_idx,'scores'] = np.mean(
                        mg_df['scores'].values.tolist()
                        )
                ahc_df_avg.at[df_avg_idx,'logfoldchanges'] = np.mean(
                        mg_df['logfoldchanges'].values.tolist()
                        )
                ahc_df_avg.at[df_avg_idx,'pvals'] = np.mean(
                        mg_df['pvals'].values.tolist()
                        )
                ahc_df_avg.at[df_avg_idx,'pvals_adj'] = np.mean(
                        mg_df['pvals_adj'].values.tolist()
                        )
                df_avg_idx += 1
            # 06b.vi. Add the resulting DataFrame to the combined one
            if len(ahc_df_final) == 0:
                ahc_df_final = ahc_df_avg.copy()
            else:
                ahc_df_final = pd.concat(
                        [ahc_df_final,
                            ahc_df_avg.copy()],
                        ignore_index = True)
        return ahc_df_final, None

# 06c. Helper functions for saving and reading distinguishing
# HVG information for the current dataset.
# To save run time, by default this DataFrame is created and
# a copy saved on the first script run, and it is read in
# on subsequent runs (unless the file is deleted).
# 06c.i. Save distinguishing HVG information 
def save_dist_HVG_data(df,
        df_fname = './adata_marker_gene_df'):
    df.to_csv(df_fname,index=True)
# 06c.ii. Read previously generated distinguishing
# information
def read_dist_HVG_data(df_fname):
    df = pd.read_csv(df_fname,index_col=0)
    return df

# 06d. Function to compute the number k folds that would
# produce large enough data subsets to power detection
# of each cluster's distinguishing HVGs
# This function assesses power for data at a given
# k value, at the specified adjusted p-value, and
# over a range of log2 fold change thresholds
over_expression_thresholds = np.arange(0.5,3.5,0.5)
def compute_n_folds(anndata,
        grouping = 'leiden',
        min_sample_size = 0):
    # 06d.i. Pull cluster sizes and identify the minimum
    sizes_df = anndata.obs.groupby(grouping,
            as_index = False
            ).size()
    sizes_var = sizes_df['size'].values.tolist()
    min_size = np.min(sizes_var)
    clusters_too_small = []
    sizes_var = [_ for _ in sizes_var if _ != min_size]
    # 06d.ii. Compute the minimum integer number of folds, k,
    # that can be used to create data subsets (by leaving one
    # fold out) that are sufficiently large to power analysis.
    # If a cluster has n nuclei, each of its data subsets contains
    # (1 - 1/k)*n nuclei.
    # This analysis is thus attempting to satisfy the inequality
    # (1 - 1/k)*n >= min_powered_sample_size
    # which is equivalent to
    # k >= 1/( 1 - (min_powered_sample_size/n) ).
    # (min_powered_sample_size is min_sample_size below)
    # Note: this approach fails for clusters where
    # min_powered_sample_size >= n; any failed clusters
    # in this search are reported.
    n_folds = 0
    while n_folds <= 0:
        if min_size > min_sample_size:
            n_folds = int(
                    np.ceil(
                        1.0/( 1.0 - (min_sample_size/min_size) )
                        )
                    )
        else:
            # 06d.iii. Record the cluster with the minimum size
            # as being too small for the current min_sample_size
            # and check the next-smallest cluster
            clusters_too_small.extend(
                    sizes_df.loc[
                        sizes_df['size'] == min_size
                        ][grouping].values.tolist()
                    )
            min_size = np.min(sizes_var)
            sizes_var = [_ for _ in sizes_var if _ != min_size]
    return n_folds, clusters_too_small

# 06e. Function to split a cluster into data
# subsets for k-fold cross-validation
def generate_folds(anndata,k,rand_state):
    # 06e.i. Randomly assign fold indices
    # of 1 through k to cluster barcodes,
    # saved in-place as an anndata .obs
    # annotation (in the column 'k{k}_fold')
    fold_col_name = f'k{k}_fold'
    anndata.obs[fold_col_name] = rand_state.randint(k,
            size=len(anndata.obs))
    return fold_col_name

# 07. Define a function to build a list of marker gene database
# hits for each cluster's distinguishing HVGs, which will be
# used to construct log likelihoods per putative cell type
# for that cluster.
def build_cell_type_marker_gene_vector(mg_df,
        mg_dict,
        group='group'):
    # 07a. For each cluster with distinguishing HVGs,
    # compute the likelihood per database cell type
    cluster_mg_cell_types_df = pd.DataFrame(index = np.unique(mg_df[group].values).tolist().sort(),
            columns = ['marker_genes_hits','cell_types_hits','l2fc_weights_hits'],
            dtype='object')
    for g,g_df in mg_df.groupby(group):
        if len(g_df) > 0:
            print(f'Processing cluster {g}...')
            mgs_hits = []
            cts_hits = []
            wts_hits = []
            rows_mgs_in_dict = g_df.loc[g_df['names'].isin(mg_dict.keys())]
            for mg_idx in rows_mgs_in_dict.index:
                mg_row = rows_mgs_in_dict.loc[mg_idx]
                mg_name = mg_row['names']
                mg_l2fc = mg_row['logfoldchanges']
                mg_cells = mg_dict[mg_name]
                # 07b. Record information about any overlaps
                # between the current distinguishing HVG and
                # any putative cell type marker genes
                if len(mg_cells) > 0:
                    for mg_cell in mg_cells:
                        # 07b.i. Record the gene name
                        mgs_hits.append(mg_name)
                        # 07b.ii. Record the current cell type with
                        # the overlapping marker gene 
                        cts_hits.append(mg_cell)
                        # 07b.iii. Record the cluster's overexpression
                        # level of the current distinguishing HVG vs.
                        # all other clusters (the log2 fold change)
                        # as a measure of the likelihood it is a true
                        # cluster marker gene
                        wts_hits.append(mg_l2fc)
            # 07c. Add the marker gene, putative cell type, and log2 fold
            # change data for all marker gene hits to the current cluster's
            # putative cell type information
            cluster_mg_cell_types_df.at[g,'marker_genes_hits'] = mgs_hits
            cluster_mg_cell_types_df.at[g,'cell_types_hits'] = cts_hits
            cluster_mg_cell_types_df.at[g,'l2fc_weights_hits'] = wts_hits
    return cluster_mg_cell_types_df

# 08. Define a function to compute each cluster's log likelihood
# for each putative cell types with marker genes that overlapp
# with its distinguishing HVGs.
def compute_cell_type_log_likelihoods(adata,
        cell_type_df,
        neg_marker_cell_type_df,
        adata_grouping = 'leiden',
        use_neg_marker_genes = False):
    cell_type_llk_df = pd.DataFrame(index = cell_type_df.index,
            columns = ['cell_types',
                'log_likelihoods',
                'max_likelihood_cell_type',
                'n_barcodes',
                'fraction_of_all_barcodes'],
            dtype = 'object')
    cell_type_llk_df['n_barcodes'] = cell_type_llk_df['n_barcodes'].astype('float32')
    cell_type_llk_df['fraction_of_all_barcodes'] = cell_type_llk_df['fraction_of_all_barcodes'].astype('float32')
    for idx in cell_type_df.index:
        print(f'Computing cell type log likelihoods for cluster {idx}...')
        row = cell_type_df.loc[idx]
        mg_list = row['marker_genes_hits']
        ct_list = row['cell_types_hits']
        l2fc_list = row['l2fc_weights_hits']
        putative_cell_types = np.unique(ct_list).tolist()
        cell_type_log_likelihoods = []
        # 08.a. Get negative marker gene hits
        nm_row = neg_marker_cell_type_df.loc[idx]
        nm_list = nm_row['marker_genes_hits']
        nm_ct_list = nm_row['cell_types_hits']
        nm_l2fc_list = nm_row['l2fc_weights_hits']
        # 08b. Compute the log likelihood that the cluster is each
        # putative cell type for which it over-expresses marker
        # genes
        for pct in putative_cell_types:
            # 08b.i. Compute positive part of log likelihood for marker genes
            ct_locs = [i for i,q in enumerate(ct_list) if q == pct]
            ct_llk_curr = np.sum([
                l2fc_list[_]
                for _ in ct_locs])
            # 08b.ii. If negative marker gene usage is specified,
            # compute the negative part of the log likelihood for any
            # negative marker genes that overlap the cluster's distinguishing
            # distinguishing HVGs
            if use_neg_marker_genes:
                nm_ct_locs = [i for i,q in enumerate(nm_ct_list) if q == pct]
                ct_llk_curr -= np.sum([
                    nm_l2fc_list[_]
                    for _ in nm_ct_locs])
            # 08c. Record the log likelihood for the current putative cell type
            cell_type_log_likelihoods.append(
                    ct_llk_curr)
        cell_type_llk_df.at[idx,'cell_types'] = putative_cell_types
        cell_type_llk_df.at[idx,'log_likelihoods'] = cell_type_log_likelihoods
        if len(cell_type_log_likelihoods) > 0:
            cell_type_llk_df.at[idx,'max_likelihood_cell_type'] = putative_cell_types[np.argmax(cell_type_log_likelihoods)]
        # 08d. Record the number of barcodes in the current cluster
        n_barcodes = len(adata.obs.loc[adata.obs[adata_grouping] == str(idx)])
        frac_all_barcodes = n_barcodes*1.0/adata.obs.shape[0]
        cell_type_llk_df.at[idx,'n_barcodes'] = n_barcodes
        cell_type_llk_df.at[idx,'fraction_of_all_barcodes'] = frac_all_barcodes
    return cell_type_llk_df

# 09. Run k-fold cross-validated distinguishing HVG detection and cell type
# log likelihood estimation for each cluster
for over_expression_threshold in over_expression_thresholds:
    print(f'Processing overexpression threshold of {over_expression_threshold}...')
    # 09a. Get the sample size needed to power one-vs.-rest distinguishing
    # HVG detection
    sample_size_needed,rough_power = get_sample_size(test_type='wilcoxon',
        min_l2fc=over_expression_threshold,
        adj_p_val=padj_dist_HVGs,
        rand_seed=get_rand_seed(numpy_random_state))
    print(f'Approximate minimum sample size needed to power marker gene' + \
            f'\ndetection at alpha={padj_dist_HVGs} for ' + \
            f'l2fc={over_expression_threshold}: {sample_size_needed}')
    print(f'Approximate power achieved with this sample size: {rough_power}')
    
    # 09b. Compute the fold size needed to achieve the minimum
    # powered sample size for all clusters
    print('Computing largest powered fold size for k-fold cross-validation...')
    n_folds, clusts_too_small = compute_n_folds(
            adatas_human_qc,
            grouping=cluster_lbl,
            min_sample_size=sample_size_needed)
    print(f'\tNumber of folds returned: {n_folds}')
    print(f'\tNumber of groups that are too small to undergo k-fold cross-validation: {len(clusts_too_small)}')
    print(f'\t\tGroups: {clusts_too_small}')

    if len(clusts_too_small) > 0:
        print(f'Not all clusters are powered for marker gene detection at this level. Skipping.')

    else:
        # 09c. Perform marker gene detection, and further, because we want the most substantial data
        # fold removal (or close to it) that powers distinguishing HVG detection, move to the next
        # overexpression threshold once a k powering all clusters is identified.
        if go_conservative_kfcv:
            n_folds += 1
            print(f'Performing more conservative k-fold clustering (leaving more data per computation)\n' + \
                    f'with {n_folds} folds.')

        # 09d. Generate or read in cluster distinguishing HVGs
        ahc_df_fn = f'{cell_typing_data_dir}/{clust_dist_HVG_fn_tag}' + \
                f'__{cluster_lbl}_pass_{pn}' + \
                f'__{over_expression_threshold}_l2fc_min__kfcv_{n_folds}.csv'
        if not os.path.isfile(ahc_df_fn):
            # 09d.i. Make a copy of the expression data for distinguishing HVG
            # detection
            print('Generating marker gene data and summary DataFrame...')
            adatas_human_qc_for_mg = adatas_human_qc.copy()
            # 09d.ii. Generate k data folds 
            fold_col_name = generate_folds(adatas_human_qc_for_mg,
                    k = n_folds,
                    rand_state = numpy_random_state)
            # 09d.iii. Perform k-fold distinguishing HVG detection
            TO_TEST = False
            ahc_df, redundancy_dict = k_fold_get_dist_HVGs(adatas_human_qc_for_mg,
                    grouping = cluster_lbl,
                    log2_min = over_expression_threshold,
                    pval_max = padj_dist_HVGs,
                    fold_col = fold_col_name,
                    k = n_folds,
                    testing = TO_TEST)
            # 09d.iv. If a redundancy dictionary is returned
            # as a result of running distinguishing HVG detection
            # in test mode (setting TO_TEST=True), plot the distribution
            # of redundancies across data subsets
            if redundancy_dict:
                test_redundancy_data_dir = f'{cell_typing_data_dir}/mg_redundancies_{cluster_lbl}_pass_{pn}__{over_expression_threshold}_l2fc_min__kfcv_{n_folds}'
                if not os.path.isdir(test_redundancy_data_dir):
                    os.system(f'mkdir -p {test_redundancy_data_dir}')
                print(f'Producing statistics on k-fold redundancy groups....')
                for key,value in redundancy_dict.items():
                    print(f'Current key: {key}: Length of marker gene repeat counts: {len(value)}')
                    plt.figure()
                    plt.hist(value)
                    plt.xlabel(f'Repeats in {n_folds}-fold testing per marker gene')
                    plt.ylabel(f'Frequency')
                    plt.title(f'{cluster_lbl} cluster {key}')
                    plt.savefig(f'{test_redundancy_data_dir}/mg_redundancies__clust_{key}.png',
                            dpi=300)
                    plt.close('all')
            # 09d.v. Save the resulting distinguishing HVG data
            save_dist_HVG_data(ahc_df,
                    ahc_df_fn)
        else:
            # 09e. If a distinguishing HVG .csv already exists, read it in
            print('Reading marker gene summary DataFrame from file...')
            ahc_df = read_dist_HVG_data(ahc_df_fn)
        
        # 09a. Generate a DataFrame showing the names, associated cell types, and 
        # relative over-expression (in the form of log2 fold change values) of
        # identified marker genes that overlap with marker genes or negative marker
        # genes of cell types in the marker gene dictionary
        cluster_cell_type_expression_df = build_cell_type_marker_gene_vector(ahc_df,
                marker_dict)
        cluster_neg_marker_expression_df = build_cell_type_marker_gene_vector(ahc_df,
                negative_marker_dict)
        
        # 09b. Compute the log likelihood that each cluster is of each of the cell types
        # for which its marker genes overlap
        clust_put_cell_type_df = compute_cell_type_log_likelihoods(adata=adatas_human_qc,
                cell_type_df=cluster_cell_type_expression_df,
                neg_marker_cell_type_df=cluster_neg_marker_expression_df,
                adata_grouping = cluster_lbl,
                use_neg_marker_genes=use_neg_marker_genes)

        # 09c. Save this DataFrame with putative cell types to file
        put_ct_fn = f'{cell_typing_data_dir}/putative_cell_types_per_clust' + \
                f'__{cluster_lbl}_pass_{pn}' + \
                f'__{over_expression_threshold}_l2fc_min__kfcv_{n_folds}.csv'
        clust_put_cell_type_df.to_csv(put_ct_fn,index=True)

        # 09d. Save the current state of the expression dataset random
        # number generator for further processing
        numpy_random_state_dict = numpy_random_state.get_state(
                legacy=False
                )
        numpy_random_state_fn = f'{output_random_state_dir}/{output_random_state_fn}'
        with open(numpy_random_state_fn,'wb') as outfile:
            pickle.dump(numpy_random_state_dict,
                    outfile)

        # 09e. Save the number of folds used in this calculation to a pickle file
        full_n_folds_output_fn = f'{kfcv_log_dir}/{n_fold_output_log_fn}'
        with open(full_n_folds_output_fn,'w') as outfile:
            outfile.write(f'{n_folds}')


        sys.exit()



