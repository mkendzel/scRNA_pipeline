#' Preview UMAP for a specific sample and resolution
#'
#' Quickly visualize clustering results for a single sample at a specified resolution.
#' Updates active identities based on the given resolution column
#' (e.g., `"RNA_snn_res.0.6"`) and generates a UMAP plot without saving it.
#'
#' @param group Character. Name of the sample in `singlet_list` to visualize.
#' @param res Numeric or character. Resolution value used in clustering (e.g., `0.6`).
#' @param dims Optional numeric vector. Principal components to use when running UMAP.
#'   If `NULL`, uses `dims_use[[group]]`.
#'
#' @return A `ggplot` object showing the UMAP embedding for the selected sample and resolution.
#' @details
#' Checks for `"RNA_snn_res.<res>"` in `obj@meta.data`.
#' If UMAP has not yet been computed, it will run `RunUMAP()` automatically.
#'
#' @examples
#' preview_umap("ZIKV-IC High", 0.7)
#' preview_umap("ZIKV, MOI 1", 0.6)
#'
#' @export

preview_umap <- function(group, res, pt.size, dims = NULL) {
  obj <- singlet_list[[group]]
  res_col <- paste0("RNA_snn_res.", res)
  if (!(res_col %in% colnames(obj@meta.data))) {
    stop(sprintf("Resolution %s not found for %s", res, group))
  }
  Idents(obj) <- obj@meta.data[[res_col]]
  
  if (!"umap" %in% names(obj@reductions)) {
    if (is.null(dims)) dims <- dims_use[[group]]
    obj <- RunUMAP(obj, dims = dims)
  }
  
  DimPlot(obj, reduction = "umap", label = TRUE, pt.size = pt.size) +
    ggtitle(sprintf("%s (res=%.2f)", group, as.numeric(res)))
}

# Example:
# preview_umap("ZIKV-IC High", 0.7)
# preview_umap("ZIKV, MOI 1", 0.6)




