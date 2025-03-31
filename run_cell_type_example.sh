#!/bin/bash

set -ev

# 01. Set up environment
exec_dir=$( pwd )
cd "${exec_dir}"
ct_dir="${exec_dir}/scripts"

# 02. Set up config files
# 02a. Specify config file path
cfg="${exec_dir}/configs/config_cell_type.yaml"
echo "${cfg}"
# 02b. Add root directory to config
# file if it does not exist in there
# yet (meaning the script hasn't been
# run before)
if ! $( grep -q 'root_dir' ${cfg} ); then
	echo "Initializing config file with current directory as root directory"
	echo "root_dir: '${exec_dir}'" >> ${cfg}
else
	echo "Config file already contains root directory; proceeding"
fi

# 03. Run Pass 1: cell typing on full dataset
# 03a. Pass 1 cell typing
echo -e "Pass 1...."
python "${ct_dir}/02a_assign_cell_types__kfcv.py" --config-yaml-path ${cfg} --pass-number "01"
# 03b. Pass 1 cell typing inspection
echo -e "Checking Pass 1 output...."
python "${ct_dir}/02b_check_clusts_w_put_cts.py" --config-yaml-path ${cfg} --pass-number "01"
# 03c. Pass 1 cluster putative cell type generation
echo -e "Generating cluster putative cell types from Pass 1...."
python "${ct_dir}/02c_gen_est_cell_types__kfcv.py" --config-yaml-path ${cfg} --pass-number "01"
# 03d. Pass 1 skeleton cluster putative cell type file for manual annotation
echo -e "Generating cluster-to-cell-type map for manual validation...."
python "${ct_dir}/02d_gen_cluster_remap_file_to_annotate.py" --config-yaml-path ${cfg} --pass-number "01"
# 03e. Print statement telling user to annotate file then run next script
#echo -e "Annotate skeleton file with top cell types if not already done before continuing (for peer review file is provided)"
#echo -e "Steps to annotate: 1. fill in top_cell_type from likelihood .txt output file; 2. run 02e_reformat_top_put_cts.py; 
#3. fill in column with Pass 1 putative cluster cell types." 
# **** 03f. Provide commented-out run statement for 02e_reformat_top_put_cts.py
# for as-needed future use
#echo -e "Standardizing top_cell_type entries to facilitate Pass 1 cluster putative cell type annotation...."
#python "${ct_dir}/02e_reformat_top_put_cts.py" --config-yaml-path ${cfg} --pass-number "01"
# ****

# 03g. Update expression data with Pass 1 cell types
echo -e "Updating expression data with Pass 1 cell types...."
python "${ct_dir}/02f_assign_new_ct_labels.py" --config-yaml-path ${cfg} --pass-number "01"
# 03h. Visualize expression of putative cell type marker genes for updated clusters, for validation purposes
echo -e "Producing visualizations of putative cell type marker gene expression for validation of annotated output clusters...."
python "${ct_dir}/02g_inspect_clust_put_ct_mgs_to_validate_clust_ct_ann.py" --config-yaml-path ${cfg} --pass-number "01"

# Print note
#echo -e "For future users: if a cell typing pass indicates that two or more clusters should be combined, or that an 
#existing cluster should be split, it is highly recommended to perform an additional cell typing on the resulting 
#clusters to validate their putative types."

# 04. Run Pass 2: subclustering
# 04a. Run subclustering on specified Pass 1 cell type clusters
echo -e "Running Pass 2 subclustering on specified Pass 1 cell type clusters...."
python "${ct_dir}/02h_perform_ct_subclustering__kfcv.py" --config-yaml-path ${cfg} --pass-number "02"
# 04b. Inspect expression of new subclusters to help determine how to isolate or recombine them
echo -e "Inspecting expression of Pass 2 new subclusters...."
python "${ct_dir}/02i_inspect_ct_subcluster_expression_labels.py" --config-yaml-path ${cfg} --pass-number "02"
# 04c. Generate Pass 2 skeleton cell-type-to-refined-cell-type map file for annotation
echo -e "Generating Pass 2 cell-type-to-refined-cell-type map for manual validation...."
python "${ct_dir}/02j_gen_cluster_remap_file_subcluster_refinement.py" --config-yaml-path ${cfg} --pass-number "02"
# 04d. Assign new Pass 2 cell type labels to expression data
echo -e "Updating expression data with Pass 2 cell types...."
python "${ct_dir}/02k_assign_new_ct_labels_subcluster_refinement.py" --config-yaml-path ${cfg} --pass-number "02"

# 05. Run Pass 3: cell typing of refined clusters
# 05a. Pass 3 cell typing
echo -e "Pass 3...."
python "${ct_dir}/02a_assign_cell_types__kfcv.py" --config-yaml-path ${cfg} --pass-number "03"
# 05b. Pass 3 cell typing inspection
echo -e "Checking Pass 3 output...."
python "${ct_dir}/02b_check_clusts_w_put_cts.py" --config-yaml-path ${cfg} --pass-number "03"
# 05c. Pass 3 cluster putative cell type generation
echo -e "Generating cluster putative cell types from Pass 3...."
python "${ct_dir}/02c_gen_est_cell_types__kfcv.py" --config-yaml-path ${cfg} --pass-number "03"
# 05d. Pass 3 skeleton cluster putative cell type file for manual annotation
echo -e "Generating cluster-to-cell-type map for manual validation...."
python "${ct_dir}/02d_gen_cluster_remap_file_to_annotate.py" --config-yaml-path ${cfg} --pass-number "03"
# 05e. Print statement telling user to annotate file then run next script
#echo -e "Annotate skeleton file with top cell types if not already done before continuing (for peer review file is provided)"
#echo -e "Steps to annotate: 1. fill in top_cell_type from likelihood .txt output file; 2. run 02e_reformat_top_put_cts.py;
#3. fill in column with Pass 3 putative cluster cell types."
# **** 05f. Provide commented-out run statement for 02e_reformat_top_put_cts.py
# for as-needed future use
# echo -e "Standardizing top_cell_type entries to facilitate Pass 3 cluster putative cell type annotation...."
# python "${ct_dir}/02e_reformat_top_put_cts.py" --config-yaml-path ${cfg} --pass-number "03"
# ****

# 05g. Update expression data with Pass 3 refined cell types and finalize these
echo -e "Updating expression data with Pass 3 cell types...."
python "${ct_dir}/02f_assign_new_ct_labels.py" --config-yaml-path ${cfg} --pass-number "03" --finalize True



