import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import sys
import json
import argparse
import yaml

# This script generates plots to help
# validate manual cluster cell type assignment,
# by plotting putative cell type marker gene
# expression for input clusters that are
# potentially of the same cell type.

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

# 01. Create an argparse object for easier reading
# of input arguments
parser = argparse.ArgumentParser(description='Visualization of marker genes for most likely' + \
        ' cell types in major cell type clusters of expression data.')
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
# 02b.iv. Output expression data directory from current cell typing pass
input_data_dir = qc_root_dir + f'/' + cfg.get(f'output_data_dir__{pn}')
# 02b.v. Output expression data file name
input_data_fn = cfg.get(f'output_data_fn__{pn}')
# 02b.vi. Cell typing data directory
cell_typing_data_dir = output_parent_dir + f'/' + cfg.get('cell_typing_data_dir')

# 02c. Cell type annotation information
# 02c.i. Annotated cluster-to-cell-type map file directory
scc_map_dir = qc_root_dir + f'/' + cfg.get(f'ann_cluster_to_ct_map_dir__{pn}')
# 02c.ii. Annotated cluster-to-cell-type map file
scc_map_fn = cfg.get(f'ann_cluster_to_ct_map_fn__{pn}')
# 02c.iii. Cell-type-to-collapsed-cell-type directory
ct_coll_ct_dir = qc_root_dir + f'/' + cfg.get('mg_db_dir')
# 02c.iv. Cell-type-to-collapsed-cell-type file name
ct_coll_ct_fn = cfg.get('ct_col_ct_fn')
# 02c.v. Input cluster label
input_ct_lbl = cfg.get(f'cluster_label__{pn}')
# 02c.vi. Output cluster label
output_ct_lbl = cfg.get(f'output_cluster_label__{pn}')

# 02d. K-fold cross-validation info (for reading
# in cluster distinguishing HVGs)
# 02d.i. N folds log file directory
n_folds_dir = qc_root_dir + f'/' + cfg.get('kfcv_n_fold_log_dir')
# 02d.ii. N folds log file name
n_folds_fn = cfg.get(f'n_folds_used_fn__02_ct__{pn}')

# 02e. Visualization directory information
# 02e.i. Visualization output directory
subclust_mg_data_dir = cell_typing_data_dir + f'/' + cfg.get(f'ct_subclust_vis_dir__{pn}')
if not os.path.isdir(subclust_mg_data_dir):
    os.system(f'mkdir -p {subclust_mg_data_dir}')

# 02b.xiv. Marker gene file name root
dist_HVG_fn_root = cfg.get(f'dist_HVG_list__{pn}')
# 02b.xv. Marker gene database directory
mg_db_dir = qc_root_dir + f'/' + cfg.get('mg_db_dir')
# 02b.xvi. Marker gene database file name
mg_db_fn = cfg.get('mg_db_fn')

# Inspect output putative cell type marker genes
# for input clusters potentially needing recombining

# 03. Import expression data
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 04. Pull n folds used in current cell typing
# pass to facilitate reading results files
n_folds_full_fn = f'{n_folds_dir}/{n_folds_fn}'
n_folds = None
with open(n_folds_full_fn,'r') as f:
    lines = f.readlines()
    for line in lines:
        n_folds = int(line)

# 05. Specify thresholds to filter input cluster distinguishing
# HVGs, to enable more conservative inspection of clusters
# to determine if they should be recombined
# 05a. Mimimum counts required for any potential distinguishing HVG
n_tot_counts_threshold = 1000
# 05b. Number of distinguishing HVGs to use per input cluster
# pulled in order of descending over-expression vs. other clusters)
n_top_genes = 20

# 06. Import distinguishing HVG data
kfcv_root = f'__kfcv_{n_folds}'
dist_HVG_fn = [
        _ for _ in os.listdir(cell_typing_data_dir)
        if (
            (dist_HVG_fn_root in _)
            and
            (kfcv_root in _)
            )
        ][0]
overexpression_threshold = float(
        dist_HVG_fn[
            :dist_HVG_fn.find('l2fc')
            ].split('__')[-1].split('_')[0]
        )
dist_HVG_fn = f'{cell_typing_data_dir}/{dist_HVG_fn}'
dist_HVG_df = pd.read_csv(dist_HVG_fn,index_col=0)

# 07. Import annotated cluster-to-cell-type map
scc_map_full_fn = f'{scc_map_dir}/{scc_map_fn}'
scc_map_df = pd.read_csv(scc_map_full_fn)
# 07a. Build a map of output clusters onto input
# clusters
output_to_input_clust_map = {}
group_to_top_ct_map = {}
for g,g_df in scc_map_df.groupby(output_ct_lbl):
    output_to_input_clust_map[g] = g_df[input_ct_lbl].values.tolist()
    all_top_cts = [
            _.strip()
            for item in
            g_df['top_collapsed_cell_types'].values.tolist()
            for _ in item.split(',')
            ]
    group_to_top_ct_map[g] = list(np.unique(all_top_cts))
    
# 08. Remove genome prefixes from dataset gene names for visualization
adatas_human_qc.var['gene_symbol'] = [
        _.split('GRCh38_______________')[-1]
        for _ in
        adatas_human_qc.var_names.copy().values.tolist()
        ]
adatas_human_qc.var.set_index('gene_symbol',
        drop=True,
        inplace=True)
adatas_human_qc.var.index.rename('gene_name',
        inplace=True)

# 09. Get marker genes for the union of collapsed
# top cell types that correspond to each output
# cluster
# 09a. Import marker gene database
marker_db_full_fn = f'{mg_db_dir}/{mg_db_fn}'
marker_db = pd.read_csv(marker_db_full_fn,index_col=0)

# 09b. Break degeneracies in database putative cell
# type labels by prepending the tissue type as needed
mdb_ct_u, mdb_ct_c = np.unique(marker_db['cell_type'],
        return_counts=True)
mdb_ct_repeat = mdb_ct_u[
        [i for i,_ in enumerate(mdb_ct_c) if _ > 1]
        ]
for c_repeat in mdb_ct_repeat:
    r_idxs_repeat = marker_db.loc[
            marker_db['cell_type'] == c_repeat
            ].index.values.tolist()
    for r_idx in r_idxs_repeat:
        tissue_row = marker_db.loc[r_idx]['tissue']
        ct = marker_db.loc[r_idx]['cell_type']
        marker_db.at[r_idx,'cell_type'] = f'{tissue_row} {ct}'

# 09c. Import cell-type-to-collapsed-cell-type map to generate
# a collapsed-cell-type marker gene database
ct_collapsed_ct_map_full_fn = f'{ct_coll_ct_dir}/{ct_coll_ct_fn}'
ct_col_map = pd.read_csv(ct_collapsed_ct_map_full_fn,index_col=0)
top_collapsed_ct_mg_df = pd.DataFrame(columns=['tissue','collapsed_cell_type','marker_genes'])
mdc_idx_curr = 0
for col_ct,ct_df in ct_col_map.groupby('collapsed_cell_type'):
    cts = ct_df['cell_type'].values.tolist()
    marker_db_mgs_raw = marker_db.loc[marker_db['cell_type'].isin(cts)]['marker_genes'].values.tolist()
    marker_db_mgs = []
    for _ in marker_db_mgs_raw:
        marker_db_mgs.extend(_.split(','))
    # 09c.i. Store union of marker genes for each collapsed cell type
    # in this new version of the marker gene database, giving each marker
    # gene its own row for easier querying
    for _ in marker_db_mgs:
        top_collapsed_ct_mg_df.at[mdc_idx_curr,'tissue'] = ct_df['tissue'].values.tolist()[0]
        top_collapsed_ct_mg_df.at[mdc_idx_curr,'collapsed_cell_type'] = col_ct
        top_collapsed_ct_mg_df.at[mdc_idx_curr,'marker_genes'] = str(_)
        mdc_idx_curr += 1

# 10. Compare top distinguishing HVGs for the input
# clusters that form each output cluster, and compare
# input cluster expression for marker genes of all
# the top putative cell types they were identified
# as during cell typing
# 10a. Helper function to generate a colorized UMAP
# for each output cluster
def plot_UMAP_mge_by_expression_island(
        adatas,
        output_clust,
        input_clust_list,
        input_clust_label,
        dh_df,
        marker_gene_df,
        most_likely_cell_types,
        n_tot_counts_threshold=1000,
        n_top_genes=n_top_genes):
    # 10a. Plot UMAP showing individual input
    # clusters in each output cluster (to show
    # any cluster that may be recombined after
    # the current cell typing pass)
    fig, ax = plt.subplots()
    sc.pl.embedding(
            adatas,
            basis='umap',
            color=input_clust_label,
            palette=color_list_1,
            groups=input_clust_list,
            title=f'2D UMAP, {output_clust}',
            show=False,
            ax=ax
            )
    plt.tight_layout()
    umap_fig_fn = f'{subclust_mg_data_dir}/2D_UMAP__{output_clust}_input_clusters' + \
            f'__{overexpression_threshold}_l2fc_min' + \
            f'__kfcv_{n_folds}.png'
    plt.savefig(
            umap_fig_fn,
            dpi=300,
            bbox_inches='tight')
    plt.close('all')

    # 10b. Get the top distinguishing HVGs for input clusters
    print(f'Collecting top {n_top_genes} distinguishing HVGs per input cluster....')
    all_top_mgs = []
    all_mgs = []
    for g_num,g_df in dh_df.groupby('group'):
        g = f'{g_num}'
        if g in input_clust_list:
            # 10c. Filter distinguishing HVGs to keep only those
            # with greater than the specified minimum number of counts
            g_df_filt = g_df.loc[g_df['n_counts_all_cells']
                    > n_tot_counts_threshold
                    ].copy()
            # 10d. Pull the top n of the resulting distinguishing HVGs
            g_df_top_n = g_df_filt.sort_values(by=['logfoldchanges'],
                    ascending=False).copy()
            g_mgs = g_df_top_n['names'].values.tolist()
            g_top_n_mgs = g_mgs[:n_top_genes]
            all_mgs.extend(g_mgs)
            all_top_mgs.extend(g_top_n_mgs)
    # 10e. Print statistics about the numbers of unique and shared
    # distinguishing HVGs among the input clusters
    print(f'\n\tDistinguishing HVGs statistics for current output cluster inputs:')
    all_mgs_u, n_repeats_all = list(np.unique(all_mgs,return_counts=True))
    print(f'\tN unique distinguishing HVGs across input clusters: {len(all_mgs_u)}')
    all_mgs_shared_some = [all_mgs_u[i] for i,_ in enumerate(n_repeats_all) if _ > 1]
    print(f'\tN distinguishing HVGs shared across >1 input clusters: {len(all_mgs_shared_some)}')
    all_mgs_shared_all = [all_mgs_u[i] for i,_ in enumerate(n_repeats_all) if _ == len(input_clust_list)]
    print(f'\tN distinguishing HVGs share across all input clusters: {len(all_mgs_shared_all)}\n')
    # 10f. Print the same statistics for only the top n distinguishing HVGs
    all_top_mgs_u, n_repeats = list(np.unique(all_top_mgs,return_counts=True))
    print(f'\tUnique top-{n_top_genes} among input clusters:\n\t{all_top_mgs_u}\n\t({len(all_top_mgs_u)} total)')
    all_top_mgs_shared_some = [all_top_mgs_u[i] for i,_ in enumerate(n_repeats) if _ > 1]
    print(f'\tTop-{n_top_genes} distinguishing HVGs shared across > 1 input cluster:' + \
            f'\n\t{all_top_mgs_shared_some}\n\t({len(all_top_mgs_shared_some)} total)')
    all_top_mgs_shared_all = [all_top_mgs_u[i] for i,_ in enumerate(n_repeats) if _ == len(input_clust_list)]
    print(f'\tTop-{n_top_genes} distinguishing HVGs shared across all input clusters:' + \
            f'\n\t{all_top_mgs_shared_all}\n\t({len(all_top_mgs_shared_all)} total)\n')
    # 10g. Remove the genome prefix from distinguishing HVG lists for consistency
    # with expression data for this visualization
    all_mgs_u_no_genome = [_.split('GRCh38_______________')[-1] for _ in all_mgs_u]
    all_top_mgs_u_no_genome = [_.split('GRCh38_______________')[-1] for _ in all_top_mgs_u]
    # 10h. Get expression data for input clusters
    adatas_dp = adatas[adatas.obs[input_clust_label].isin(input_clust_list)].copy()
    # 10i. Make a collapsed-cell-type marker gene dictionary for all input cluster distinguishing HVGs
    marker_dict_curr_all = {}
    md_rows_all = marker_gene_df.loc[marker_gene_df['marker_genes'].isin(
        all_mgs_u_no_genome)]
    if len(md_rows_all) > 0:
        for ct_dict,rows_ct_dict in md_rows_all.groupby('collapsed_cell_type'):
            marker_dict_curr_all[ct_dict] = rows_ct_dict['marker_genes'].values.tolist()
    # 10j. Make a similar dictionary for only the top n distinguishing HVGs
    marker_dict_curr_top_n = {}
    # 10j.i. Top n distinguishing HVGs that are in the marker gene dictionary
    md_rows = marker_gene_df.loc[marker_gene_df['marker_genes'].isin(all_top_mgs_u_no_genome)]
    # 10j.ii. Top n distinguishing HVGs not in the marker gene dictionary
    non_md_top_mgs = list(set(all_top_mgs_u_no_genome) - set(md_rows['marker_genes'].values.tolist()))
    if len(md_rows) > 0:
        for ct_dict,rows_ct_dict in md_rows.groupby('collapsed_cell_type'):
            marker_dict_curr_top_n[ct_dict] = rows_ct_dict['marker_genes'].values.tolist()
    marker_dict_curr_top_n['other'] = non_md_top_mgs
    
    # 10k. Make a dotplot comparing all overexpressed genes for the current output cluster
    # to marker genes for collapsed cell types
    plt.figure()
    plt.rcParams.update({'font.size': 6})
    ax = plt.subplot(1,1,1)
    sc.pl.dotplot(adatas_dp,
            var_names = marker_dict_curr_all,
            groupby = input_clust_label,
            dendrogram = True,
            dot_min = 0.0,
            dot_max = 1.0,
            smallest_dot = 0.0,
            figsize = (20,10),
            var_group_rotation = 0.0,
            title = f'Cell Type Marker Genes, All DEGs, {output_clust}\n',
            show = False,
            ax = ax)
    dp_fig_fn = f'{subclust_mg_data_dir}/diff_exp_by_{input_clust_label}_cluster__all_de_genes__{output_clust}' + \
            f'__kfcv_{n_folds}.png'
    plt.savefig(dp_fig_fn,
            dpi=200)
    plt.close('all')
    # 10k. Make a dotplot showing top-n gene overexpressed genes that correspond
    # to collapsed cell types for clusters in the island
    plt.figure()
    plt.rcParams.update({'font.size': 6})
    ax = plt.subplot(1,1,1)
    sc.pl.dotplot(adatas_dp,
            var_names = marker_dict_curr_top_n,
            groupby = input_clust_label,
            dendrogram = True,
            dot_min = 0.0,
            dot_max = 1.0,
            smallest_dot = 0.0,
            figsize = (20,10),
            var_group_rotation = 0.0,
            title = f'Cell Type Marker Genes, Top-{n_top_genes} DEGs, {output_clust}\n',
            show = False,
            ax = ax)
    dp_fig_fn = f'{subclust_mg_data_dir}/diff_exp_by_{input_clust_label}_cluster__top_{n_top_genes}_de_genes__{output_clust}' + \
            f'__kfcv_{n_folds}.png'
    plt.savefig(dp_fig_fn,
            dpi=200)
    plt.close('all')
    # 10l. Clear the copy of the anndata object to avoid memory issues
    del adatas_dp

# 11. Generate UMAPs showing the constituents of each output cluster
# and the relationships between their distinguishing HVGs and putative
# cell type marker genes
for output_clust_lbl_curr, ahc_obs_curr in adatas_human_qc.obs.groupby(output_ct_lbl):
    print(f'Output cluster: {output_clust_lbl_curr}')
    input_clust_lbl_list_curr = list(np.unique(ahc_obs_curr[input_ct_lbl].values.tolist()))
    print(f'Related input clusters: {input_clust_lbl_list_curr}')
    plot_UMAP_mge_by_expression_island(
            adatas=adatas_human_qc,
            output_clust=output_clust_lbl_curr,
            input_clust_list=input_clust_lbl_list_curr,
            input_clust_label=input_ct_lbl,
            dh_df=dist_HVG_df,
            marker_gene_df=top_collapsed_ct_mg_df,
            most_likely_cell_types = group_to_top_ct_map[output_clust_lbl_curr],
            n_tot_counts_threshold=n_tot_counts_threshold,
            n_top_genes=n_top_genes)

sys.exit()

