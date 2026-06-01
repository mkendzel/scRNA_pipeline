#' Preprocess and pick PCs (stop before clustering)
#' @param obj Seurat object
#' @param assay assay to use
#' @param max_pcs PCs to compute in RunPCA
#' @param var_cutoff cumulative variance target (e.g., 0.90)
#' @return list(obj=Seurat with PCA, dims_use=integer vector of PCs)
prep_precluster <- function(obj, assay = "RNA", max_pcs = 50, 
                            var_cutoff = 0.90, nfeatures) {
  DefaultAssay(obj) <- assay
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = max_pcs)
  
  sdev <- obj[["pca"]]@stdev
  var  <- sdev^2 / sum(sdev^2)
  k    <- which(cumsum(var) >= var_cutoff)[1]
  
  list(obj = obj, dims_use = 1:min(k, max_pcs))
}