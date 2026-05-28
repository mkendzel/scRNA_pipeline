#' Run additive linear models for a subset of genes within a cluster
#'
#' Fits: expression ~ viral_load + treat_group
#'
#' @param obj Seurat object
#' @param cluster_id Cluster identifier (matching obj$seurat_clusters)
#' @param gene_list Character vector of gene names
#' @param assay Assay to use (default "RNA")
#' @param slot Data slot to use (default "data")
#'
#' @return Data frame with coefficients and p-values for each gene
run_additive_subset <- function(obj, cluster_id, gene_list,
                                assay = "RNA", slot = "data") {
  
  # Extract normalized expression matrix
  expr <- GetAssayData(obj, assay = assay, slot = slot)
  
  # Subset cells in the chosen cluster
  cells <- which(obj$seurat_clusters == cluster_id)
  meta  <- obj@meta.data[cells, ]
  
  # Subset expression for selected genes
  expr_mat <- expr[gene_list, cells, drop = FALSE]
  
  # Fit model for each gene
  results <- lapply(gene_list, function(g) {
    y   <- as.numeric(expr_mat[g, ])
    fit <- lm(y ~ meta$viral_load + meta$treat_group)
    cf  <- summary(fit)$coefficients
    
    data.frame(
      gene          = g,
      cluster       = cluster_id,
      coef_viral_load = cf["meta$viral_load", "Estimate"],
      p_viral_load    = cf["meta$viral_load", "Pr(>|t|)"],
      coef_treat_IC   = cf["meta$treat_groupIC", "Estimate"],
      p_treat_IC      = cf["meta$treat_groupIC", "Pr(>|t|)"]
    )
  })
  
  do.call(rbind, results)
}