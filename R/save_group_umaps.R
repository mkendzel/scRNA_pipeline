#' Save UMAPs for all groups across available resolutions
#'
#' For each object in `singlet_list`, detect all clustering resolution columns
#' (`RNA_snn_res.*`), set identities, generate a UMAP plot, and save one image
#' per resolution inside a group-specific subfolder of `out_dir`.
#' If a UMAP embedding is missing, it is computed using `dims_use[[group]]`.
#'
#' @param singlet_list Named list of Seurat objects, one per group.
#' @param dims_use Named list mapping group name -> numeric vector of PCs (e.g., `1:20`)
#'   used by `RunUMAP()` when the UMAP reduction is absent.
#' @param out_dir Character. Top-level output folder. One subfolder per group is created.
#' @param pt.size Numeric. Point size passed to `DimPlot()`.
#' @param width Numeric. Plot width (in inches) for `ggsave()`.
#' @param height Numeric. Plot height (in inches) for `ggsave()`.
#' @param dpi Integer. Resolution (dots per inch) for `ggsave()`.
#'
#' @return Invisibly returns `singlet_list` (updated in memory if UMAP was computed).
#'
#' @details
#' Resolution columns are detected via `grep("^RNA_snn_res\\.", colnames(obj@meta.data))`.
#' Files are written as `<out_dir>/<group>/UMAP_res-<res>.png`.
#'
#' @examples
#' # save_group_umaps(singlet_list, dims_use, out_dir = "figures/umap_by_res", pt.size = 1.2)
#'
#' @export

save_group_umaps <- function(singlet_list, dims_use,
                             out_dir = "figures/umap_by_res",
                             pt.size = 1.2, width = 6, height = 5, dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (group in names(singlet_list)) {
    obj <- singlet_list[[group]]
    
    # Run UMAP if missing
    if (!"umap" %in% names(obj@reductions)) {
      obj <- RunUMAP(obj, dims = dims_use[[group]])
      singlet_list[[group]] <- obj
    }
    
    # Find all available resolutions
    res_cols <- grep("^RNA_snn_res\\.", colnames(obj@meta.data), value = TRUE)
    res_vals <- sub("^RNA_snn_res\\.", "", res_cols)
    
    group_dir <- file.path(out_dir, group)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Plot and save
    for (res in res_vals) {
      Idents(obj) <- obj@meta.data[[paste0("RNA_snn_res.", res)]]
      p <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = pt.size) +
        ggtitle(sprintf("%s (res=%s)", group, res))
      ggsave(file.path(group_dir, sprintf("UMAP_res-%s.png", res)),
             plot = p, width = width, height = height, dpi = dpi)
    }
  }
  
  invisible(singlet_list)
}
