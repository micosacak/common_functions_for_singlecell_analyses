Repository Summary: common_functions_for_singlecell_analyses
This is a personal utility library for single-cell multiome (scRNA-seq + scATAC-seq) analyses, primarily working with zebrafish (Danio rerio) data using Seurat/Signac in R, with a Python bridge for Python-based tools. Here's what each file does:

`77_CommonFunctionS5.R` — The core utility file
This is the heart of the repo. It contains many helper functions grouped by purpose:
Environment & setup

ldFnx() — Master loader: sets working directories, loads all required R libraries (Seurat, Signac, CellChat, FigR, etc.), and registers parallel cores.
ldData() — Loads and caches heavy reference data: gene annotations for human/mouse/zebrafish, JASPAR transcription factor motif databases, and cross-species ortholog tables. Saves .rds files to avoid re-downloading.
setDir(), crtD(), crtDfull() — Helper functions for creating and navigating output directories.
setRAM(), setCPU() — Shortcuts to configure memory and CPU usage for parallel computing.

Gene/ortholog utilities

getOrthos() — Fetches ortholog mappings between any two species (human, mouse, zebrafish) via Ensembl BioMart.
find_match() / find_match_improved() — Fuzzy gene name matching: finds the closest matching gene name in a dataset even with typos, capitalisation differences, or partial names.

Plotting

getVln() — Generates violin plots with jitter and quartile lines from a data frame.
getQuantiles() — Computes per-group quantile thresholds and returns which values fall above them (useful for filtering outlier cells).
DoHeatmap_custom() — Custom heatmap using ComplexHeatmap: takes a Seurat object and a named list of gene groups, draws a grouped heatmap with colour-coded row annotations.
getDensPlot(), getDensPlot1(), getDensPlot2() — UMAP density plots per sample or group, with options for 2D contour or heatmap-style density colouring.
get_DP1() — Normalises density plots across groups to share a common colour scale.
ggplotColours(), get_vCols(), oo() — Colour palette helpers and plot size shortcuts.
my.ggvenn() / prepare_venn_data() — Custom Venn diagram that overlays two different count sets (e.g. gene counts from two conditions) in the same diagram with different coloured labels.
get_factors() — Creates ordered factor combinations from multiple grouping vectors (useful for multi-condition UMAP colouring).

Gene regulatory network (FigR/DORC)

runGenePeakcorr_dre() — Extended version of FigR's peak–gene correlation function, adapted for zebrafish (GRCz11 genome). Finds which chromatin-accessible peaks are correlated with nearby gene expression.
run_FigGRN_dre_final() — Runs the full FigR gene regulatory network pipeline for zebrafish: matches TF motifs to peaks, computes TF–DORC (domain of regulatory chromatin) associations.
get_network_data_old() — Prepares network graph data (igraph object) from FigR output for visualisation.
get_figR_rds_files(), get_figR_folders() — Utility functions to load batches of FigR result files from a directory.


cellchat_analyses.R — Cell–cell communication pipeline

getCellChat() — Wrapper that builds a complete CellChat object from a Seurat object: creates the CellChat object, identifies overexpressed ligand–receptor interactions, computes communication probabilities, and calculates network centrality scores. Returns a fully analysed CellChat object ready for plotting.
The rest of the file is a worked example script showing how to run CellChat in parallel across two conditions (control vs treatment), select communication pattern numbers (k), identify incoming/outgoing patterns, compute functional/structural network similarity, merge conditions for comparison, and generate bubble/chord plots.


cellChatPlottings.R — CellChat visualisation helpers
Contains custom modifications to CellChat's built-in plotting functions (e.g. netVisual_bubble_modified), allowing finer control over bubble plot aesthetics for comparing ligand–receptor interactions between conditions.

custom_plottings.R — Custom visualisation functions

get_fold_RectPlot() — Creates a bidirectional horizontal bar (lollipop/waterfall) chart for comparing fold-changes across cell clusters and samples. Includes threshold lines to highlight biologically meaningful fold-changes. Useful for showing cell-type composition shifts between conditions.


export_seu5_signac_for_anndata.R — Seurat → flat files exporter (R side)
A robust exporter that converts a Seurat/Signac object into flat files readable by Python. Works with ATAC-only, RNA-only, or multiome objects, and handles Seurat v4/v5 compatibility. Exports: raw counts (MTX format), normalised data layers, dimensional reductions (UMAP, PCA, LSI), peak metadata, pre-computed marker genes, colour palettes, cell–cell graphs, and a manifest.json that describes everything exported. Designed to pair with convert_to_anndata.py.

convert_to_anndata.py — Flat files → AnnData converter (Python side)
The Python counterpart of the exporter above. Takes the directory of flat files produced by export_seurat() and assembles them into AnnData objects ready for tools like CREsted, scanpy, or scrublet.

files_to_anndata() — Main function: reads the manifest, loads ATAC and/or RNA count matrices, normalises peak names to chrN:start-end format (handling multiple naming conventions), attaches cell metadata, and optionally creates pseudobulk profiles (summing raw counts per cell type).
_make_pseudobulk() — Aggregates per-cell ATAC counts into one row per cell type, skipping types with too few cells.
_normalise_peak_names() — Standardises peak names from various formats (Signac dashes, underscores, missing chr prefix) into a uniform chrN:start-end format.
save_anndatas() — Saves the resulting AnnData objects as .h5ad files with type-safe columns (prevents common h5py write errors).


generate_custom_cellchatDB.R — Custom CellChat ligand–receptor database
Script for building a custom CellChat ligand–receptor interaction database, likely extended with zebrafish-specific or manually curated interactions not present in the default CellChatDB.

installPackages.R — Package installer
One-time script to install all required R packages (Bioconductor + CRAN + GitHub packages) needed for the analyses in this repository.

create_conda_env.sh — Environment setup
Shell script to create a conda environment with the required Python packages (anndata, scanpy, scipy, etc.) for the Python components of the workflow.
