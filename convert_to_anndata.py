#
#Note: this code was partially written or optimized by Claude.ai
#
# =============================================================================
# files_to_anndata.py
# Convert files exported by the R function `export_seurat()` into AnnData
# objects ready for CREsted model training.
#
# Usage:
#   from files_to_anndata import files_to_anndata, save_anndatas
#   adata_atac, adata_rna = files_to_anndata("crested_export")
#
# Requirements:
#   pip install anndata scipy pandas numpy
#
# Fixes applied vs v1:
#   - MTX orientation auto-detected (fixes all-zero pseudobulk matrix)
#   - Peak names: chr prefix added if missing (e.g. "1-100-200" → "chr1:100-200")
#   - Peak name parser handles dash/underscore/colon AND bare chromosome numbers
#   - _build_peak_var: no more None values (h5py write fix)
#   - save_anndatas: var columns cast to correct dtypes before writing
# =============================================================================

import json
import warnings
from pathlib import Path
from typing import Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import scanpy as sc

# -----------------------------------------------------------------------------
# Helper: load a sparse MTX matrix + features + barcodes → AnnData
# Auto-detects whether MTX is (features × barcodes) or (barcodes × features)
# -----------------------------------------------------------------------------
def _load_matrix_as_anndata(
    mtx_path: str,
    features_path: str,
    barcodes_path: str,
    feature_col: str,
) -> ad.AnnData:
    """Load a MatrixMarket file + feature/barcode CSVs into an AnnData.

    Seurat's WriteMM writes features × barcodes (peaks × cells), so the
    standard approach is to transpose (.T) to get cells × features.
    However we verify dimensions against the CSV files and correct if needed.
    """
    mat      = scipy.io.mmread(mtx_path).tocsr()
    features = pd.read_csv(features_path)[feature_col].astype(str).values
    barcodes = pd.read_csv(barcodes_path)["barcode"].astype(str).values

    n_features = len(features)
    n_barcodes = len(barcodes)

    # Auto-detect orientation
    if mat.shape == (n_features, n_barcodes):
        # Standard Seurat output: features x barcodes -> transpose to cells x features
        mat = mat.T.tocsr()
    elif mat.shape == (n_barcodes, n_features):
        # Already cells x features - no transpose needed
        mat = mat.tocsr()
    else:
        raise ValueError(
            f"MTX shape {mat.shape} does not match "
            f"features ({n_features}) x barcodes ({n_barcodes}) "
            f"in either orientation. Check your export files."
        )

    assert mat.shape == (n_barcodes, n_features), \
        f"After orientation fix, expected ({n_barcodes}, {n_features}), got {mat.shape}"

    adata           = ad.AnnData(X=mat)
    adata.obs_names = barcodes
    adata.var_names = features
    return adata


# -----------------------------------------------------------------------------
# Helper: normalise peak names to "chrN:start-end" format
# Handles ALL of:
#   "chr1-1000-2000"   (Signac default)
#   "1-1000-2000"      (missing chr prefix)
#   "chr1:1000-2000"   (already correct)
#   "chr1_1000_2000"   (underscore-separated)
#   "1_1000_2000"      (underscore, no chr)
# -----------------------------------------------------------------------------
def _normalise_peak_names(names: np.ndarray) -> np.ndarray:
    import re
    # Matches optional "chr" prefix, chromosome name, any separator, start, end
    pattern = re.compile(r"^(chr)?([^\W_\-:]+)[:\-_](\d+)[:\-_](\d+)$")
    out = []
    for n in names:
        m = pattern.match(str(n))
        if m:
            chr_prefix = m.group(1) or "chr"   # add "chr" if missing
            chrom      = chr_prefix + m.group(2)
            start      = m.group(3)
            end        = m.group(4)
            out.append(f"{chrom}:{start}-{end}")
        else:
            out.append(str(n))  # leave unchanged
    return np.array(out)


# -----------------------------------------------------------------------------
# Helper: build peak-coordinate DataFrame for adata.var
# All values are non-null (h5py safe)
# -----------------------------------------------------------------------------
def _build_peak_var(peak_names: np.ndarray) -> pd.DataFrame:
    chroms, starts, ends = [], [], []
    for name in peak_names:
        try:
            chrom, coords = name.split(":")
            start, end    = coords.split("-")
            chroms.append(str(chrom))
            starts.append(int(start))
            ends.append(int(end))
        except Exception:
            # Non-parseable peak - use safe defaults (never None)
            chroms.append("unknown")
            starts.append(-1)
            ends.append(-1)

    return pd.DataFrame(
        {
            "chrom":      pd.array(chroms, dtype="str"),
            "chromStart": pd.array(starts, dtype="int64"),
            "chromEnd":   pd.array(ends,   dtype="int64"),
        },
        index=peak_names,
    )


# -----------------------------------------------------------------------------
# Pseudobulk helper: sum raw counts per cell type -> (cell_types x peaks)
# -----------------------------------------------------------------------------
def _make_pseudobulk(
    adata: ad.AnnData,
    celltype_col: str,
    min_cells_per_type: int,
    verbose: bool,
) -> ad.AnnData:
    """Aggregate per-cell counts into pseudobulk per cell type."""

    cell_types = adata.obs[celltype_col].astype(str)
    unique_cts = cell_types.unique()
    is_sparse  = scipy.sparse.issparse(adata.X)

    ct_list, mat_rows, ncells_list = [], [], []

    for ct in sorted(unique_cts):
        mask = (cell_types == ct).values
        n    = mask.sum()
        if n < min_cells_per_type:
            warnings.warn(
                f"Cell type '{ct}' has only {n} cells "
                f"(< {min_cells_per_type}); skipping.",
                stacklevel=2,
            )
            continue

        if is_sparse:
            row = np.asarray(adata.X[mask].sum(axis=0)).flatten()
        else:
            row = adata.X[mask].sum(axis=0).flatten()

        mat_rows.append(row)
        ct_list.append(ct)
        ncells_list.append(int(n))

    if len(ct_list) == 0:
        raise ValueError(
            "No cell types passed the min_cells_per_type filter. "
            f"Try lowering min_cells_per_type (current: {min_cells_per_type})."
        )

    pb_matrix = scipy.sparse.csr_matrix(np.vstack(mat_rows))  # cell_types x peaks

    # Sanity check: pseudobulk must not be all zeros
    if pb_matrix.nnz == 0:
        raise ValueError(
            "Pseudobulk matrix is all zeros!\n"
            "This usually means the MTX was read in the wrong orientation.\n"
            "Check that atac_counts.mtx, atac_features.csv, and atac_barcodes.csv "
            "were exported correctly from Seurat."
        )

    adata_pb           = ad.AnnData(X=pb_matrix)
    adata_pb.obs_names = ct_list
    adata_pb.var_names = adata.var_names.copy()
    adata_pb.var       = adata.var.copy()
    adata_pb.obs["cell_type"] = ct_list
    adata_pb.obs["n_cells"]   = ncells_list

    if verbose:
        total_counts = pb_matrix.sum()
        density      = 100 * pb_matrix.nnz / (pb_matrix.shape[0] * pb_matrix.shape[1])
        print(f"  Pseudobulk : {adata_pb.n_obs} cell types x {adata_pb.n_vars} peaks")
        print(f"  Total counts       : {total_counts:,.0f}")
        print(f"  Non-zero entries   : {pb_matrix.nnz:,} ({density:.1f}% dense)")

    return adata_pb


# -----------------------------------------------------------------------------
# Main conversion function
#
# Parameters
# ----------
# export_dir : str
#     Path to the directory produced by the R function `export_seurat()`.
# celltype_col : str
#     Column in cell_metadata.csv with cell type labels. Default: "cell_type".
# pseudobulk : bool
#     If True (default), return pseudobulk AnnData (cell types x peaks).
#     CREsted training requires pseudobulk=True.
# min_cells_per_type : int
#     Cell types with fewer cells than this are excluded. Default: 10.
# verbose : bool
#     Print progress messages. Default: True.
#
# Returns
# -------
# (adata_atac, adata_rna)
#   adata_atac : AnnData  - pseudobulk or per-cell ATAC (always returned)
#   adata_rna  : AnnData  - per-cell RNA (multiome only, else None)
# -----------------------------------------------------------------------------
def files_to_anndata(
    export_dir: str,
    celltype_col: str = "cell_type",
    pseudobulk: bool = True,
    min_cells_per_type: int = 10,
    verbose: bool = True,
) -> Tuple[ad.AnnData, Optional[ad.AnnData]]:

    export_dir = Path(export_dir)
    if not export_dir.exists():
        raise FileNotFoundError(f"Export directory not found: {export_dir}")

    # --- 1. Load manifest ---------------------------------------------------
    manifest_path = export_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            "manifest.json not found. Did you run export_seurat() in R first?"
        )
    with open(manifest_path) as f:
        manifest = json.load(f)

    modality   = manifest["modality"]
    export_rna = manifest.get("export_rna", False)

    if verbose:
        print(f"Modality          : {modality}")
        print(f"Cells             : {manifest['n_cells']:,}")
        print(f"Cell types        : {manifest['n_cell_types']}  "
              f"({', '.join(manifest['cell_types'][:5])}"
              f"{'...' if len(manifest['cell_types']) > 5 else ''})")

    # --- 2. Load cell metadata ----------------------------------------------
    meta = pd.read_csv(export_dir / "cell_metadata.csv")
    meta.index = meta["barcode"].astype(str)

    if celltype_col not in meta.columns:
        raise ValueError(
            f"Column '{celltype_col}' not found in cell_metadata.csv.\n"
            f"Available columns: {list(meta.columns)}"
        )
    cell_types = meta[celltype_col].astype(str)

    # --- 3. Load ATAC -------------------------------------------------------
    adata_atac = None
    if modality in ("atac_only", "multiome"):
        if verbose:
            print("\nLoading ATAC counts ...")

        adata_atac = _load_matrix_as_anndata(
            mtx_path      = export_dir / "atac_counts.mtx",
            features_path = export_dir / "atac_features.csv",
            barcodes_path = export_dir / "atac_barcodes.csv",
            feature_col   = "peak",
        )

        if verbose:
            print(f"  Raw matrix     : {adata_atac.n_obs:,} cells x {adata_atac.n_vars:,} peaks")
            nnz = adata_atac.X.nnz if scipy.sparse.issparse(adata_atac.X) else np.count_nonzero(adata_atac.X)
            print(f"  Non-zero counts: {nnz:,}")

        # Normalise peak names to "chrN:start-end"
        norm_peaks           = _normalise_peak_names(adata_atac.var_names.values)
        adata_atac.var_names = norm_peaks
        adata_atac.var       = _build_peak_var(norm_peaks)

        # Report unparseable peaks
        n_unknown = (adata_atac.var["chrom"] == "unknown").sum()
        if n_unknown > 0:
            warnings.warn(f"{n_unknown} peak names could not be parsed -> chrom='unknown'")
        elif verbose:
            print(f"  Peak names     : all {adata_atac.n_vars:,} parsed successfully "
                  f"(e.g. {norm_peaks[0]})")

        # Attach cell metadata
        adata_atac.obs               = meta.reindex(adata_atac.obs_names).copy()
        adata_atac.obs[celltype_col] = cell_types.reindex(adata_atac.obs_names).values

        # Load optional Signac peak-level metadata
        peak_meta_path = export_dir / "peak_metadata.csv"
        if peak_meta_path.exists():
            peak_meta = pd.read_csv(peak_meta_path).set_index("peak_name")
            shared = adata_atac.var_names.intersection(peak_meta.index)
            if len(shared) > 0:
                for col in peak_meta.columns:
                    adata_atac.var[col] = peak_meta.reindex(adata_atac.var_names)[col].values

        # --- 4. Pseudobulk --------------------------------------------------
        if pseudobulk:
            adata_atac = _make_pseudobulk(
                adata_atac,
                celltype_col       = celltype_col,
                min_cells_per_type = min_cells_per_type,
                verbose            = verbose,
            )

    # --- 5. Load RNA (multiome only) ----------------------------------------
    adata_rna = None
    if modality == "multiome" and export_rna:
        if verbose:
            print("\nLoading RNA counts ...")
        adata_rna = _load_matrix_as_anndata(
            mtx_path      = export_dir / "rna_counts.mtx",
            features_path = export_dir / "rna_features.csv",
            barcodes_path = export_dir / "rna_barcodes.csv",
            feature_col   = "gene",
        )
        adata_rna.obs               = meta.reindex(adata_rna.obs_names).copy()
        adata_rna.obs[celltype_col] = cell_types.reindex(adata_rna.obs_names).values
        if verbose:
            print(f"  RNA : {adata_rna.n_obs:,} cells x {adata_rna.n_vars:,} genes")

    # --- 6. Summary ---------------------------------------------------------
    if verbose:
        print("\nConversion complete.")
        if adata_atac is not None:
            label = "cell types" if pseudobulk else "cells"
            print(f"  adata_atac : {adata_atac.n_obs} {label} x {adata_atac.n_vars} peaks")
        if adata_rna is not None:
            print(f"  adata_rna  : {adata_rna.n_obs:,} cells x {adata_rna.n_vars:,} genes")
        print("\n  Next step -> crested.pp.normalize_peaks(adata_atac)")

    return adata_atac, adata_rna


# -----------------------------------------------------------------------------
# Save both AnnData objects to disk (h5py-safe dtypes enforced)
# -----------------------------------------------------------------------------
def save_anndatas(
    adata_atac: ad.AnnData,
    adata_rna: Optional[ad.AnnData] = None,
    out_dir: str = "anndata_output",
) -> None:
    """Save converted AnnData objects to .h5ad files with h5py-safe dtypes."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Ensure var columns have h5py-safe dtypes (no mixed types, no None/NaN)
    for col in adata_atac.var.columns:
        if adata_atac.var[col].dtype == object:
            adata_atac.var[col] = adata_atac.var[col].fillna("unknown").astype(str)
        elif pd.api.types.is_float_dtype(adata_atac.var[col]):
            adata_atac.var[col] = adata_atac.var[col].fillna(-1.0)
        elif pd.api.types.is_integer_dtype(adata_atac.var[col]):
            adata_atac.var[col] = adata_atac.var[col].fillna(-1).astype("int64")

    atac_path = out / "atac_pseudobulk.h5ad"
    adata_atac.write_h5ad(atac_path)
    print(f"Saved: {atac_path}")

    if adata_rna is not None:
        rna_path = out / "rna_percell.h5ad"
        adata_rna.write_h5ad(rna_path)
        print(f"Saved: {rna_path}")

# -----------------------------------------------------------------------------
# Example usage (run as script)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    import sys
    import os
    os.chdir("post_analyses_by_scanpy")
    export_dir = sys.argv[1] if len(sys.argv) > 1 else "mymultiome_seu_exported"
    export_dir = "elavl3_multiome_seu_exported"
    
    adata_atac, adata_rna = files_to_anndata(
        export_dir         = f"{export_dir}",
        celltype_col       = "main2",
        pseudobulk         = True,
        min_cells_per_type = 3,
        verbose            = True,
    )
    # Load Seurat UMAP
    umap = pd.read_csv(f"{export_dir}/reduction_umap_gex.csv", index_col=0)
    print(umap.columns.tolist())
    print(umap.shape)
    
    # Inject into adata
    adata_rna.obsm["X_umap"] = umap.loc[adata_rna.obs_names].values
    import scipy.io, pandas as pd
    
    raw = scipy.io.mmread(f"{export_dir}/rna_raw_counts.mtx").T.tocsr()
    barcodes = pd.read_csv(f"{export_dir}/rna_raw_barcodes.csv")["barcode"].values
    cell_order = pd.Index(barcodes).get_indexer(adata_rna.obs_names)
    adata_rna.layers["counts"] = raw[cell_order, :]
    
    # Test UMAP
    import pandas as pd
    import scanpy as sc
    
    # Plot
    sc.pl.umap(adata_rna, color="main2", title="RNA – Seurat UMAP")
  
    save_anndatas(adata_atac, adata_rna, out_dir="mymultiome_anndata_output")
