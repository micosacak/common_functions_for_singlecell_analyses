# Note this code has been partly optimised by Claude.ai!
# The standart export sometimes raises errors
#
# Standart export!
#library(Seurat)
#library(SeuratDisk)
# Convert and save with embeddings intact
#SaveH5Seurat(seu, filename = "mymultiome_multiome.h5seurat", overwrite = TRUE)
#Convert("mymultiome_multiome.h5seurat", dest = "h5ad", overwrite = TRUE)
#
# PLEASE USE BELOW CODE FOR SAFER EXPORT!
#
# =============================================================================
# seurat_to_files.R  (improved)
# Export a Seurat object (ATAC-only, RNA-only, or Multiome) to flat files
# readable by the companion Python function `files_to_anndata()`.
#
# Improvements over v1:
#   - Exports all dimensional reductions (UMAP, PCA, LSI, …)
#   - Exports log-normalised data layer (so Python needs zero preprocessing)
#   - Exports raw counts separately → adata.layers["counts"] (needed for scrublet/DESeq2)
#   - Exports pre-computed marker genes (if available)
#   - Exports cluster-colour palettes stored by Seurat
#   - Safer Seurat v4 / v5 compatibility throughout
#   - Richer manifest + per-step timing
#
# Usage:
#   source("seurat_to_files.R")
#   export_seurat(seu, out_dir = "my_export", celltype_col = "seurat_cluster")
# =============================================================================

library(Matrix)
library(jsonlite)

# ── Internal helpers ──────────────────────────────────────────────────────────

.detect_modality <- function(seu) {
  nms <- tolower(names(seu@assays))
  has_rna  <- any(nms %in% c("rna", "scrnaseq", "gex", "gene expression"))
  has_atac <- any(nms %in% c("peaks", "atac", "chromatin", "atacseq"))
  if (has_rna && has_atac) "multiome"
  else if (has_atac)        "atac_only"
  else if (has_rna)         "rna_only"
  else                      "unknown"
}

.find_assay <- function(seu, patterns) {
  nms <- names(seu@assays)
  hit <- nms[tolower(nms) %in% patterns]
  if (length(hit) == 0) return(NULL)
  if (length(hit) >  1) message("  Multiple matches (", paste(hit, collapse=", "), "); using: ", hit[1])
  hit[1]
}

# Seurat v4 / v5 safe layer getter ───────────────────────────────────────────
.get_layer <- function(seu, assay, layer = "counts") {
  tryCatch(
    LayerData(seu, assay = assay, layer = layer),   # v5
    error = function(e)
      tryCatch(
        GetAssayData(seu, assay = assay, slot = layer),  # v4
        error = function(e2) NULL
      )
  )
}

# Peak-name → BED parser ──────────────────────────────────────────────────────
.peaks_to_bed <- function(peak_names) {
  norm <- gsub("_(?=\\d)", "-", peak_names, perl = TRUE)
  norm <- gsub(":", "-", norm)
  parts <- strsplit(norm, "-")
  ok    <- lengths(parts) == 3
  if (!all(ok)) warning(sum(!ok), " peak name(s) could not be parsed and are skipped.")
  parts <- parts[ok]
  data.frame(
    chr       = sapply(parts, `[`, 1),
    start     = as.integer(sapply(parts, `[`, 2)),
    end       = as.integer(sapply(parts, `[`, 3)),
    peak_name = peak_names[ok],
    stringsAsFactors = FALSE
  )
}

# ── NEW: export all dimensional reductions ────────────────────────────────────
.export_reductions <- function(seu, out_dir, verbose) {
  reds <- names(seu@reductions)
  if (length(reds) == 0) {
    if (verbose) message("  No reductions found in Seurat object.")
    return(character(0))
  }
  exported <- character(0)
  for (red in reds) {
    emb <- tryCatch(Embeddings(seu, reduction = red), error = function(e) NULL)
    if (is.null(emb) || nrow(emb) == 0) next
    fname <- paste0("reduction_", tolower(red), ".csv")
    write.csv(as.data.frame(emb), file.path(out_dir, fname))
    if (verbose) message("  Wrote ", fname, "  [", nrow(emb), " cells x ", ncol(emb), " dims]")
    exported <- c(exported, setNames(fname, red))
  }
  exported
}

# ── NEW: export normalised / scaled data layer ────────────────────────────────
.export_data_layer <- function(seu, assay, out_dir, prefix, verbose) {
  exported <- list(data = NULL, scale_data = NULL)

  # ── data layer (log-normalised) ───────────────────────────────────────────
  mat <- .get_layer(seu, assay, "data")
  if (is.null(mat)) {
    if (verbose) message("    No 'data' layer found for assay ", assay, " – skipping.")
  } else {
    counts <- .get_layer(seu, assay, "counts")
    if (!is.null(counts) && identical(dim(mat), dim(counts)) &&
        isTRUE(all.equal(mat, counts))) {
      if (verbose) message("    'data' layer identical to counts for ", assay, " – skipping.")
    } else {
      fname <- paste0(prefix, "_data.mtx")
      writeMM(mat, file.path(out_dir, fname))
      if (verbose) message("    Wrote ", fname,
                           "  [log-normalised, ", nrow(mat), " x ", ncol(mat), "]")
      exported$data <- fname
    }
  }

  # ── scale.data layer (z-scored dense — GEX typically has this) ───────────
  scaled <- tryCatch(.get_layer(seu, assay, "scale.data"), error = function(e) NULL)
  if (!is.null(scaled) && length(scaled) > 0) {
    fname_rds <- paste0(prefix, "_scale_data.rds")
    saveRDS(scaled, file.path(out_dir, fname_rds))
    if (verbose) message("    Wrote ", fname_rds,
                         "  [scaled dense, ", nrow(scaled), " x ", ncol(scaled), "]")
    exported$scale_data <- fname_rds
  }

  exported
}

# ── NEW: export pre-computed markers ─────────────────────────────────────────
.export_markers <- function(seu, out_dir, verbose) {
  # Look in common slots where analysts stash markers
  markers <- NULL
  for (slot_name in c("markers", "all.markers", "deg", "FindAllMarkers")) {
    markers <- tryCatch(seu@misc[[slot_name]], error = function(e) NULL)
    if (!is.null(markers) && is.data.frame(markers) && nrow(markers) > 0) {
      if (verbose) message("  Found pre-computed markers in seu@misc$", slot_name)
      break
    }
    markers <- NULL
  }
  if (is.null(markers)) {
    if (verbose) message("  No pre-computed markers found in seu@misc – skipping.")
    return(NULL)
  }
  fname <- "markers.csv"
  write.csv(markers, file.path(out_dir, fname), row.names = FALSE)
  if (verbose) message("  Wrote ", fname, "  [", nrow(markers), " rows]")
  fname
}

# ── NEW: export Seurat colour palettes ───────────────────────────────────────
.export_colours <- function(seu, out_dir, verbose) {
  # Seurat stores palettes in seu@misc as <col>_colors
  colour_keys <- grep("_colors$|_colours$", names(seu@misc), value = TRUE)
  if (length(colour_keys) == 0) return(NULL)
  palettes <- lapply(seu@misc[colour_keys], function(x) as.list(x))
  fname <- "colour_palettes.json"
  writeLines(toJSON(palettes, pretty = TRUE, auto_unbox = TRUE),
             file.path(out_dir, fname))
  if (verbose) message("  Wrote ", fname, "  [", length(colour_keys), " palette(s)]")
  fname
}

# ── NEW: export graph / neighbour info ───────────────────────────────────────
.export_graphs <- function(seu, out_dir, verbose) {
  graphs <- names(seu@graphs)
  if (length(graphs) == 0) return(character(0))
  exported <- character(0)
  for (g in graphs) {
    snn <- tryCatch(seu@graphs[[g]], error = function(e) NULL)
    if (is.null(snn)) next
    fname <- paste0("graph_", g, ".mtx")
    writeMM(snn, file.path(out_dir, fname))
    if (verbose) message("  Wrote ", fname)
    exported <- c(exported, setNames(fname, g))
  }
  exported
}

# =============================================================================
# Main export function
#
# Parameters:
#   seu              Seurat object
#   out_dir          Output directory (created if absent)
#   celltype_col     Metadata column with cell-type labels (NULL → active Idents)
#   export_rna       Export RNA modality in multiome mode (default TRUE)
#   export_reductions  Export UMAP / PCA / LSI etc. (default TRUE)         ← NEW
#   export_data        Export log-normalised data layer (default TRUE)      ← NEW
#   export_raw_counts  Export raw counts as separate layer files (TRUE)     ← NEW
#                      → rna_raw_counts.mtx loaded into adata.layers["counts"]
#                      → required for scrublet, DESeq2, and pydeseq2
#   export_markers     Export pre-computed markers from seu@misc (default TRUE) ← NEW
#   export_graphs      Export SNN graphs (default FALSE – can be large)     ← NEW
#   assay_atac       Override ATAC assay name (NULL = auto)
#   assay_rna        Override RNA assay name  (NULL = auto)
#   verbose          Print progress (default TRUE)
# =============================================================================
export_seurat <- function(seu,
                          out_dir            = "crested_export",
                          celltype_col       = NULL,
                          export_rna         = TRUE,
                          export_reductions  = TRUE,
                          export_data        = TRUE,
                          export_raw_counts  = TRUE,
                          export_markers     = TRUE,
                          export_graphs      = FALSE,
                          assay_atac         = NULL,
                          assay_rna          = NULL,
                          verbose            = TRUE) {

  t0 <- proc.time()
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ── 1. Modality ─────────────────────────────────────────────────────────────
  modality <- .detect_modality(seu)
  if (verbose) message("[1/9] Modality: ", modality)
  if (modality == "unknown") stop("Cannot detect any RNA or ATAC assay.")

  assay_atac <- assay_atac %||%
    .find_assay(seu, c("peaks", "atac", "chromatin", "atacseq"))
  assay_rna  <- assay_rna  %||%
    .find_assay(seu, c("rna", "scrnaseq", "gex", "gene expression"))

  if (verbose && !is.null(assay_atac)) message("    ATAC assay : ", assay_atac)
  if (verbose && !is.null(assay_rna))  message("    RNA assay  : ", assay_rna)

  # ── 2. Cell metadata ────────────────────────────────────────────────────────
  if (verbose) message("[2/9] Exporting cell metadata …")
  barcodes <- colnames(seu)
  cell_types <- if (is.null(celltype_col)) {
    if (verbose) message("    Using active Idents.")
    as.character(Idents(seu))
  } else {
    if (!celltype_col %in% colnames(seu@meta.data))
      stop("Column '", celltype_col, "' not found. Available: ",
           paste(colnames(seu@meta.data), collapse = ", "))
    as.character(seu@meta.data[[celltype_col]])
  }
  meta <- seu@meta.data
  meta$barcode   <- barcodes
  meta$cell_type <- cell_types
  write.csv(meta, file.path(out_dir, "cell_metadata.csv"), row.names = FALSE)
  if (verbose) message("    Wrote cell_metadata.csv  [", nrow(meta), " cells x ",
                       ncol(meta), " cols]")

  # ── 3. ATAC counts ──────────────────────────────────────────────────────────
  atac_data_file <- NULL
  if (modality %in% c("atac_only", "multiome")) {
    if (verbose) message("[3/9] Exporting ATAC counts …")
    atac_counts <- .get_layer(seu, assay_atac, "counts")
    writeMM(atac_counts, file.path(out_dir, "atac_counts.mtx"))
    write.csv(data.frame(peak    = rownames(atac_counts)),
              file.path(out_dir, "atac_features.csv"), row.names = FALSE)
    write.csv(data.frame(barcode = colnames(atac_counts)),
              file.path(out_dir, "atac_barcodes.csv"), row.names = FALSE)
    bed <- .peaks_to_bed(rownames(atac_counts))
    write.table(bed[, c("chr","start","end","peak_name")],
                file.path(out_dir, "peaks.bed"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    if (verbose) message("    Wrote atac_counts.mtx  [", nrow(atac_counts),
                         " peaks x ", ncol(atac_counts), " cells]")

    # peak-level metadata (Signac)
    tryCatch({
      pm <- seu@assays[[assay_atac]]@meta.data
      if (nrow(pm) > 0) {
        pm$peak_name <- rownames(pm)
        write.csv(pm, file.path(out_dir, "peak_metadata.csv"), row.names = FALSE)
        if (verbose) message("    Wrote peak_metadata.csv")
      }
    }, error = function(e) NULL)

    # normalised ATAC data layer
    if (export_data) {
      atac_layers    <- .export_data_layer(seu, assay_atac, out_dir, "atac", verbose)
      atac_data_file <- atac_layers$data
    }
  } else {
    if (verbose) message("[3/9] Skipping ATAC (modality = ", modality, ")")
  }

  # ── 4. RNA counts ───────────────────────────────────────────────────────────
  rna_data_file  <- NULL
  rna_raw_file   <- NULL
  rna_scale_file <- NULL
  if (modality %in% c("multiome", "rna_only") && export_rna && !is.null(assay_rna)) {
    if (verbose) message("[4/9] Exporting RNA counts …")
    rna_counts <- .get_layer(seu, assay_rna, "counts")
    writeMM(rna_counts, file.path(out_dir, "rna_counts.mtx"))
    write.csv(data.frame(gene    = rownames(rna_counts)),
              file.path(out_dir, "rna_features.csv"), row.names = FALSE)
    write.csv(data.frame(barcode = colnames(rna_counts)),
              file.path(out_dir, "rna_barcodes.csv"), row.names = FALSE)
    if (verbose) message("    Wrote rna_counts.mtx  [", nrow(rna_counts),
                         " genes x ", ncol(rna_counts), " cells]")

    # normalised RNA data layer
    if (export_data) {
      rna_layers       <- .export_data_layer(seu, assay_rna, out_dir, "rna", verbose)
      rna_data_file    <- rna_layers$data
      rna_scale_file   <- rna_layers$scale_data
    }

    # ── raw counts as dedicated layer files (for scrublet / DESeq2) ──────────
    rna_raw_file <- NULL
    if (export_raw_counts) {
      raw <- .get_layer(seu, assay_rna, "counts")
      if (!is.null(raw)) {
        writeMM(raw, file.path(out_dir, "rna_raw_counts.mtx"))
        write.csv(data.frame(gene    = rownames(raw)),
                  file.path(out_dir, "rna_raw_features.csv"), row.names = FALSE)
        write.csv(data.frame(barcode = colnames(raw)),
                  file.path(out_dir, "rna_raw_barcodes.csv"), row.names = FALSE)
        rna_raw_file <- "rna_raw_counts.mtx"
        if (verbose) message("    Wrote rna_raw_counts.mtx  [raw integer counts → adata.layers['counts']]")
      } else {
        if (verbose) message("    Raw counts not found for assay ", assay_rna, " – skipping.")
      }
    }
  } else {
    if (verbose) message("[4/9] Skipping RNA export.")
  }

  # ── 5. Dimensional reductions ────────────────────────────────────────────────
  reduction_files <- character(0)
  if (export_reductions) {
    if (verbose) message("[5/9] Exporting dimensional reductions …")
    reduction_files <- .export_reductions(seu, out_dir, verbose)
    if (length(reduction_files) == 0 && verbose)
      message("    (none exported)")
  } else {
    if (verbose) message("[5/9] Skipping reductions (export_reductions = FALSE)")
  }

  # ── 6. Markers ───────────────────────────────────────────────────────────────
  markers_file <- NULL
  if (export_markers) {
    if (verbose) message("[6/9] Exporting markers …")
    markers_file <- .export_markers(seu, out_dir, verbose)
  } else {
    if (verbose) message("[6/9] Skipping markers (export_markers = FALSE)")
  }

  # ── 7. Colour palettes & graphs ─────────────────────────────────────────────
  if (verbose) message("[7/9] Exporting colour palettes …")
  colour_file <- .export_colours(seu, out_dir, verbose)

  # ── 8. Graphs ────────────────────────────────────────────────────────────────
  if (verbose) message("[8/9] Exporting graphs …")
  graph_files <- if (export_graphs) .export_graphs(seu, out_dir, verbose) else {
    if (verbose) message("    Skipping graphs (export_graphs = FALSE)") ; character(0)
  }

  # ── 9. Manifest ─────────────────────────────────────────────────────────────
  if (verbose) message("[9/9] Writing manifest …")
  elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

  manifest <- list(
    modality         = modality,
    n_cells          = length(barcodes),
    n_cell_types     = length(unique(cell_types)),
    cell_types       = sort(unique(cell_types)),
    assay_atac       = assay_atac,
    assay_rna        = assay_rna,
    seurat_version   = tryCatch(as.character(packageVersion("Seurat")),
                                error = function(e) "unknown"),
    export_rna       = export_rna && modality %in% c("multiome", "rna_only"),
    elapsed_sec      = elapsed,
    reductions       = as.list(reduction_files),          # NEW
    files = list(
      cell_metadata  = "cell_metadata.csv",
      # ATAC
      atac_counts    = if (modality != "rna_only") "atac_counts.mtx"   else NULL,
      atac_features  = if (modality != "rna_only") "atac_features.csv" else NULL,
      atac_barcodes  = if (modality != "rna_only") "atac_barcodes.csv" else NULL,
      peaks_bed      = if (modality != "rna_only") "peaks.bed"         else NULL,
      atac_data      = atac_data_file,                                  # NEW
      # RNA
      rna_counts     = if (modality %in% c("multiome","rna_only") && export_rna) "rna_counts.mtx"   else NULL,
      rna_features   = if (modality %in% c("multiome","rna_only") && export_rna) "rna_features.csv" else NULL,
      rna_barcodes   = if (modality %in% c("multiome","rna_only") && export_rna) "rna_barcodes.csv" else NULL,
      rna_data       = rna_data_file,                                   # NEW
      rna_scale_data = rna_scale_file,                                  # NEW scale.data
      # Raw counts — loaded into adata.layers["counts"] in Python
      # Use for: scrublet doublet detection, DESeq2, pydeseq2
      rna_raw_counts   = rna_raw_file,                                  # NEW
      rna_raw_features = if (!is.null(rna_raw_file)) "rna_raw_features.csv" else NULL,
      rna_raw_barcodes = if (!is.null(rna_raw_file)) "rna_raw_barcodes.csv" else NULL,
      # Extras
      markers        = markers_file,                                    # NEW
      colour_palettes = colour_file,                                    # NEW
      graphs         = as.list(graph_files)                             # NEW
    )
  )

  writeLines(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE),
             file.path(out_dir, "manifest.json"))

  if (verbose) {
    message("\n✓ Export complete in ", elapsed, "s → ", normalizePath(out_dir))
    message("  Reductions exported : ",
            if (length(reduction_files) > 0) paste(names(reduction_files), collapse=", ") else "none")
    message("  Pass this directory to: files_to_anndata('", out_dir, "')")
  }

  invisible(manifest)
}

# Null-coalescing operator (base R ≥ 4.4 has it; define here for safety)
`%||%` <- function(a, b) if (!is.null(a)) a else b
