#' Add SoupX grouping clusters to a Seurat object
#'
#' Runs a lightweight Seurat workflow (normalization, variable feature detection,
#' scaling, PCA, neighbor finding, and clustering) to assign cluster identities
#' and store them directly as `soup_group` for SoupX contamination estimation.
#'
#' @param sobj Seurat object containing raw or filtered counts. Works also for Lists
#'
#' @return The same Seurat object with a metadata column `soup_group`.
#' @examples
#' \dontrun{
#' seurat <- add_soup_groups(seurat)
#' head(seurat$soup_group)
#' }
#' @export
add_soup_groups <- function(sobj) {
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(sobj, nfeatures = 2000,
                               verbose = FALSE, selection.method = "vst")
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  
  sobj@meta.data$soup_group <- sobj@meta.data$seurat_clusters
  sobj@meta.data$seurat_clusters <- NULL
  sobj
}