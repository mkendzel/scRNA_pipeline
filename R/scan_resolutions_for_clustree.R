#' Scan clustering resolutions and keep all for clustree
#' @param obj Seurat object already PCA'd
#' @param dims integer vector of PCs to use, e.g. 1:20
#' @param res sequence of resolutions to test
#' @return list(obj=Seurat with all RNA_snn_res.* columns, summary=data.frame)
scan_resolutions_for_clustree <- function(obj,
                                          dims,
                                          res = seq(0.2, 1.2, by = 0.1)) {
  obj <- FindNeighbors(obj, dims = dims)
  nclust <- numeric(length(res))
  for (i in seq_along(res)) {
    r <- res[i]
    obj <- FindClusters(obj, resolution = r)
    nclust[i] <- length(levels(Idents(obj)))
  }
  list(obj = obj, summary = data.frame(resolution = res, n_clusters = nclust))
}