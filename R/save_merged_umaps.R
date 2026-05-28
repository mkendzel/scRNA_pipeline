#' Save UMAPs across available clustering resolutions for a single Seurat object
#'
#' Detects all clustering resolution columns (`RNA_snn_res.*`) in the Seurat object,
#' sets identities, generates UMAP plots, and saves one image per resolution inside `out_dir`.
#' If a UMAP embedding is missing, it is computed using the provided `dims`.
#'
#' @param obj A Seurat object.
#' @param dims Numeric vector of PCs (e.g., `1:20`) used by `RunUMAP()` if UMAP is absent.
#' @param out_dir Character. Output folder where UMAP plots will be saved.
#' @param pt.size Numeric. Point size passed to `DimPlot()`.
#' @param width Numeric. Plot width (in inches) for `ggsave()`.
#' @param height Numeric. Plot height (in inches) for `ggsave()`.
#' @param dpi Integer. Resolution (dots per inch) for `ggsave()`.
#'
#' @return Invisibly returns the updated Seurat object (if UMAP was computed).
#'
#' @details
#' Resolution columns are detected via `grep("^RNA_snn_res\\.", colnames(obj@meta.data))`.
#' Files are written as `<out_dir>/UMAP_res-<res>.png`.
#'
#' @examples
#' # save_umaps_single(seurat_obj, dims = 1:20, out_dir = "figures/umap_by_res", pt.size = 1.2)
#'


save_merged_umaps <- function(obj, dims = 1:20,
                                  out_dir = "figures/merged/umap_by_res",
                                  pt.size = 1.2, width = 6, height = 5, dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Run UMAP if missing
  if (!"umap" %in% names(obj@reductions)) {
    obj <- RunUMAP(obj, dims = dims)
  }
  
  # Find all available resolutions
  res_cols <- grep("^RNA_snn_res\\.", colnames(obj@meta.data), value = TRUE)
  res_vals <- sub("^RNA_snn_res\\.", "", res_cols)
  
  # Plot and save
  for (res in res_vals) {
    Idents(obj) <- obj@meta.data[[paste0("RNA_snn_res.", res)]]
    p <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = pt.size) +
      ggtitle(sprintf("UMAP (res=%s)", res))
    ggsave(file.path(out_dir, sprintf("UMAP_res-%s.png", res)),
           plot = p, width = width, height = height, dpi = dpi)
  }
  
  invisible(obj)
}