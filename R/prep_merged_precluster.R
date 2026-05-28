#' Preprocess and pick PCs (stop before clustering)
#' @param obj Seurat object
#' @param assay assay to use
#' @param max_pcs PCs to compute in RunPCA
#' @param var_cutoff cumulative variance target (e.g., 0.90)
#' @param nfeatures number of variable features to select
#' @param vars.to.regress variables to regress out in ScaleData
#' @param exclude_features character vector of features to exclude from PCA (e.g., viral)
#' @return list(obj=Seurat with PCA, dims_use=integer vector of PCs)
prep_merged_precluster <- function(
    obj,
    assay = "RNA",
    max_pcs = 50,
    var_cutoff = 0.90,
    nfeatures,
    vars.to.regress = c("nCount_RNA","percent.mt"),
    exclude_features = NULL
) {
  DefaultAssay(obj) <- assay
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj, vars.to.regress = vars.to.regress)
  
  pca_features <- VariableFeatures(obj)
  if (!is.null(exclude_features)) {
    pca_features <- setdiff(pca_features, exclude_features)
  }
  
  obj <- RunPCA(obj, features = pca_features, npcs = max_pcs)
  
  sdev <- obj[["pca"]]@stdev
  var  <- sdev^2 / sum(sdev^2)
  k    <- which(cumsum(var) >= var_cutoff)[1]
  
  list(obj = obj, dims_use = 1:min(k, max_pcs))
}
