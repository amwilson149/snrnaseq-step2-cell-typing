
# 01. Set up environment
# 01a. Activate conda environment
ml anaconda3/2020.11
module  load 'anaconda3/2020.11'
Shell debugging temporarily silenced: export LMOD_SH_DBG_ON=1 for Lmod's output

The following have been reloaded with a version change:
  1) gcc/14.2.0 => gcc/8.3.0

Shell debugging restarted
ml -python
module  load '-python'
Shell debugging temporarily silenced: export LMOD_SH_DBG_ON=1 for Lmod's output
Shell debugging restarted
source /hpc/packages/minerva-centos7/anaconda3/2020.11/etc/profile.d/conda.sh
export CONDA_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/conda'
export _CE_M=''
export _CE_CONDA=''
export CONDA_PYTHON_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/python'

# Copyright (C) 2012 Anaconda, Inc
# SPDX-License-Identifier: BSD-3-Clause

__add_sys_prefix_to_path() {
    # In dev-mode CONDA_EXE is python.exe and on Windows
    # it is in a different relative location to condabin.
    if [ -n "${_CE_CONDA}" ] && [ -n "${WINDIR+x}" ]; then
        SYSP=$(\dirname "${CONDA_EXE}")
    else
        SYSP=$(\dirname "${CONDA_EXE}")
        SYSP=$(\dirname "${SYSP}")
    fi

    if [ -n "${WINDIR+x}" ]; then
        PATH="${SYSP}/bin:${PATH}"
        PATH="${SYSP}/Scripts:${PATH}"
        PATH="${SYSP}/Library/bin:${PATH}"
        PATH="${SYSP}/Library/usr/bin:${PATH}"
        PATH="${SYSP}/Library/mingw-w64/bin:${PATH}"
        PATH="${SYSP}:${PATH}"
    else
        PATH="${SYSP}/bin:${PATH}"
    fi
    \export PATH
}

__conda_hashr() {
    if [ -n "${ZSH_VERSION:+x}" ]; then
        \rehash
    elif [ -n "${POSH_VERSION:+x}" ]; then
        :  # pass
    else
        \hash -r
    fi
}

__conda_activate() {
    if [ -n "${CONDA_PS1_BACKUP:+x}" ]; then
        # Handle transition from shell activated with conda <= 4.3 to a subsequent activation
        # after conda updated to >= 4.4. See issue #6173.
        PS1="$CONDA_PS1_BACKUP"
        \unset CONDA_PS1_BACKUP
    fi

    \local cmd="$1"
    shift
    \local ask_conda
    CONDA_INTERNAL_OLDPATH="${PATH}"
    __add_sys_prefix_to_path
    ask_conda="$(PS1="$PS1" "$CONDA_EXE" $_CE_M $_CE_CONDA shell.posix "$cmd" "$@")" || \return $?
    rc=$?
    PATH="${CONDA_INTERNAL_OLDPATH}"
    \eval "$ask_conda"
    if [ $rc != 0 ]; then
        \export PATH
    fi
    __conda_hashr
}

__conda_reactivate() {
    \local ask_conda
    CONDA_INTERNAL_OLDPATH="${PATH}"
    __add_sys_prefix_to_path
    ask_conda="$(PS1="$PS1" "$CONDA_EXE" $_CE_M $_CE_CONDA shell.posix reactivate)" || \return $?
    PATH="${CONDA_INTERNAL_OLDPATH}"
    \eval "$ask_conda"
    __conda_hashr
}

conda() {
    if [ "$#" -lt 1 ]; then
        "$CONDA_EXE" $_CE_M $_CE_CONDA
    else
        \local cmd="$1"
        shift
        case "$cmd" in
            activate|deactivate)
                __conda_activate "$cmd" "$@"
                ;;
            install|update|upgrade|remove|uninstall)
                CONDA_INTERNAL_OLDPATH="${PATH}"
                __add_sys_prefix_to_path
                "$CONDA_EXE" $_CE_M $_CE_CONDA "$cmd" "$@"
                \local t1=$?
                PATH="${CONDA_INTERNAL_OLDPATH}"
                if [ $t1 = 0 ]; then
                    __conda_reactivate
                else
                    return $t1
                fi
                ;;
            *)
                CONDA_INTERNAL_OLDPATH="${PATH}"
                __add_sys_prefix_to_path
                "$CONDA_EXE" $_CE_M $_CE_CONDA "$cmd" "$@"
                \local t1=$?
                PATH="${CONDA_INTERNAL_OLDPATH}"
                return $t1
                ;;
        esac
    fi
}

if [ -z "${CONDA_SHLVL+x}" ]; then
    \export CONDA_SHLVL=0
    # In dev-mode CONDA_EXE is python.exe and on Windows
    # it is in a different relative location to condabin.
    if [ -n "${_CE_CONDA+x}" ] && [ -n "${WINDIR+x}" ]; then
        PATH="$(\dirname "$CONDA_EXE")/condabin${PATH:+":${PATH}"}"
    else
        PATH="$(\dirname "$(\dirname "$CONDA_EXE")")/condabin${PATH:+":${PATH}"}"
    fi
    \export PATH

    # We're not allowing PS1 to be unbound. It must at least be set.
    # However, we're not exporting it, which can cause problems when starting a second shell
    # via a first shell (i.e. starting zsh from bash).
    if [ -z "${PS1+x}" ]; then
        PS1=
    fi
fi
conda activate CO-cell-typing-env
PS1='(CO-cell-typing-env) '
export PATH='/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/bin:/hpc/packages/minerva-centos7/anaconda3/2020.11/bin:/hpc/packages/minerva-centos7/gcc/8.3.0_32b/bin:/hpc/packages/minerva-centos7/anaconda3/2020.11/condabin:/hpc/users/wilsoa28/google-cloud-sdk/bin:/hpc/users/wilsoa28/git-filter-repo:/hpc/packages/minerva-rocky9/git/2.46.0/bin:/hpc/packages/minerva-common/vim/8.0/bin:/hpc/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/hpc/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/bin:/usr/bin:/usr/mbin:/local/bin:/usr/local:/usr/ucb:/usr/local/sbin:/usr/sbin:/usr/lpp/mmfs/bin:/hpc/users/wilsoa28/.local/bin:/hpc/users/wilsoa28/bin'
export CONDA_PREFIX='/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env'
export CONDA_SHLVL='1'
export CONDA_DEFAULT_ENV='CO-cell-typing-env'
export CONDA_PROMPT_MODIFIER='(CO-cell-typing-env) '
export CONDA_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/conda'
export _CE_M=''
export _CE_CONDA=''
export CONDA_PYTHON_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/python'
. "/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/etc/conda/activate.d/libglib_activate.sh"
export GSETTINGS_SCHEMA_DIR_CONDA_BACKUP="${GSETTINGS_SCHEMA_DIR:-}"
export GSETTINGS_SCHEMA_DIR="$CONDA_PREFIX/share/glib-2.0/schemas"
. "/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/etc/conda/activate.d/libxml2_activate.sh"
#!/bin/sh

if test -n "${XML_CATALOG_FILES:-}"; then
    xml_catalog_files_libxml2="${XML_CATALOG_FILES}"
    XML_CATALOG_FILES="${XML_CATALOG_FILES} "
else
    xml_catalog_files_libxml2=""
    XML_CATALOG_FILES=""
fi


# Replace space with '%20'; equivalent to
# conda_catalog_files=${CONDA_PREFIX// /%20}, except trailing space is
# ignored.
conda_catalog_files=""
ifs_libxml2="${IFS}"
IFS=" "
rem="${CONDA_PREFIX}"
for pre in ${rem}; do
    while test "${rem#"${pre}"}" = "${rem}"; do
	conda_catalog_files="${conda_catalog_files}%20"
	rem=${rem#" "}
    done
    conda_catalog_files="${conda_catalog_files}${pre}"
    rem=${rem#"${pre}"}
done
IFS="${ifs_libxml2}"

conda_catalog_files="file://${conda_catalog_files}/etc/xml/catalog file:///etc/xml/catalog"
export XML_CATALOG_FILES="${XML_CATALOG_FILES}${conda_catalog_files}"
unset conda_catalog_files ifs_libxml2 rem

# 01b. Get root directory
exec_dir=$( pwd )
cd "${exec_dir}"
# 01c. Define code path
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
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
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
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
WARNING: dendrogram data not found (using key=dendrogram_leiden). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
WARNING: Groups are not reordered because the `groupby` categories and the `var_group_labels` are different.
categories: 6, 7, 21
var_group_labels: Astrocytes, B cells, Brain cancer cells, etc.
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
  dot_ax.scatter(x, y, **kwds)
WARNING: Groups are not reordered because the `groupby` categories and the `var_group_labels` are different.
categories: 6, 7, 21
var_group_labels: Smooth muscle cells, other
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
  dot_ax.scatter(x, y, **kwds)
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
WARNING: dendrogram data not found (using key=dendrogram_leiden). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
WARNING: Groups are not reordered because the `groupby` categories and the `var_group_labels` are different.
categories: 4, 9, 16
var_group_labels: Astrocytes, B cells, Brain Endothelial cells, etc.
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
  dot_ax.scatter(x, y, **kwds)
WARNING: Groups are not reordered because the `groupby` categories and the `var_group_labels` are different.
categories: 4, 9, 16
var_group_labels: B cells, Brain Endothelial cells, Erythroid-like and erythroid precursor cells, etc.
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
  dot_ax.scatter(x, y, **kwds)
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
WARNING: Dendrogram not added. Dendrogram is added only when the number of categories to plot > 2
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
  dot_ax.scatter(x, y, **kwds)
WARNING: Dendrogram not added. Dendrogram is added only when the number of categories to plot > 2
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_dotplot.py:749: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
  dot_ax.scatter(x, y, **kwds)

# Print note
#echo -e "For future users: if a cell typing pass indicates that two or more clusters should be combined, or that an 
#existing cluster should be split, it is highly recommended to perform an additional cell typing on the resulting 
#clusters to validate their putative types."

# 04. Run Pass 2: subclustering
# 04a. Run subclustering on specified Pass 1 cell type clusters
echo -e "Running Pass 2 subclustering on specified Pass 1 cell type clusters...."
python "${ct_dir}/02h_perform_ct_subclustering__kfcv.py" --config-yaml-path ${cfg} --pass-number "02"
Detecting highly variable genes with seurat_v3...# 04b. Inspect expression of new subclusters to help determine how to isolate or recombine them
echo -e "Inspecting expression of Pass 2 new subclusters...."
python "${ct_dir}/02i_inspect_ct_subcluster_expression_labels.py" --config-yaml-path ${cfg} --pass-number "02"
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
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
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
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
/sc/arion/work/wilsoa28/.conda/envs/CO-cell-typing-env/lib/python3.8/site-packages/anndata/_core/anndata.py:1235: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
  df[key] = c




