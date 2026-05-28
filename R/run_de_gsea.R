#' Differential Expression and GSEA Analysis
#'
#' Performs differential expression (DE) analysis between two groups in a Seurat object
#' and runs Gene Set Enrichment Analysis (GSEA) using `fgsea`.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param group_col Column name in `seurat_obj@meta.data` used to define identities.
#' @param ident_treat Identity for the treatment group.
#' @param ident_base Identity for the baseline/control group.
#' @param pathways A list of gene sets (e.g., Hallmark pathways) for GSEA.
#' @param assay Assay to use for DE analysis. Default is `"RNA"`.
#' @param test.use Statistical test for DE (passed to `FindMarkers`). Default is `"wilcox"`.
#' @param logfc.threshold Minimum log fold-change threshold for DE. Default is `0`.
#' @param min.pct Minimum fraction of cells expressing the gene. Default is `0.01`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{de}{Data frame of differential expression results from `FindMarkers`.}
#'   \item{gsea}{Data frame of GSEA results from `fgsea`.}
#' }
#'
#' @details
#' Sets the default assay, assigns identities based on `group_col`,
#' computes DE between `ident_treat` and `ident_base`, ranks genes by log2 fold-change,
#' and runs GSEA using the provided pathways.
#'
#' @examples
#' \dontrun{
#' gsea_results <- run_de_gsea(
#'   seurat_obj = my_seurat,
#'   group_col = "condition",
#'   ident_treat = "treated",
#'   ident_base = "control",
#'   pathways = hallmark_pathways
#' )
#' }
#'

run_de_gsea <- function(seurat_obj,
                        group_col,
                        ident_treat,
                        ident_base,
                        pathways,
                        assay = "RNA",
                        test.use = "wilcox",
                        logfc.threshold = 0,
                        min.pct = 0.01,
                        exclude_features = NULL) {
  
  DefaultAssay(seurat_obj) <- assay
  Idents(seurat_obj) <- group_col
  
  # DE: treatment vs baseline
  features_use <- rownames(seurat_obj[[assay]])
  if (!is.null(exclude_features) && length(exclude_features) > 0) {
    exclude_features <- intersect(exclude_features, features_use)
    features_use <- setdiff(features_use, exclude_features)
  }
  
  de <- FindMarkers(
    seurat_obj,
    ident.1         = ident_treat,
    ident.2         = ident_base,
    logfc.threshold = logfc.threshold,
    min.pct         = min.pct,
    assay           = assay,
    test.use        = test.use,
    features        = features_use
  )
  
  de$gene <- rownames(de)
  
  # Replace p-values of zero with the smallest positive double -------------- 
  smallest_pos <- .Machine$double.xmin * .Machine$double.eps               
  de$p_val_fixed <- ifelse(de$p_val == 0, smallest_pos, de$p_val)          
  
  ### ============================
  ### THREE RANKING METHODS
  ### ============================
  
  # 1) log2FC --------------------------------------------------------------- 
  ranks_logfc <- de$avg_log2FC                                              
  names(ranks_logfc) <- de$gene                                             
  ranks_logfc <- ranks_logfc[is.finite(ranks_logfc)]                        
  ranks_logfc <- sort(ranks_logfc, decreasing = TRUE)                       
  
  # 2) sign(log2FC) * -log10(pval_fixed) ----------------------------------- 
  ranks_sign_p <- sign(de$avg_log2FC) * -log10(de$p_val_fixed)              
  names(ranks_sign_p) <- de$gene                                            
  ranks_sign_p <- ranks_sign_p[is.finite(ranks_sign_p)]                     
  ranks_sign_p <- sort(ranks_sign_p, decreasing = TRUE)                     
  
  # 3) log2FC * -log10(pval_fixed) ------------------------------------------ 
  ranks_fc_p <- de$avg_log2FC * -log10(de$p_val_fixed)                      
  names(ranks_fc_p) <- de$gene                                              
  ranks_fc_p <- ranks_fc_p[is.finite(ranks_fc_p)]                           
  ranks_fc_p <- sort(ranks_fc_p, decreasing = TRUE)                         
  
  ### ============================
  ### Run GSEA
  ### ============================
  
  gsea_logfc <- fgsea(pathways = pathways, stats = ranks_logfc)             
  gsea_sign_p <- fgsea(pathways = pathways, stats = ranks_sign_p)           
  gsea_fc_p   <- fgsea(pathways = pathways, stats = ranks_fc_p)             
  
  # Add comparison labels --------------------------------------------------- 
  gsea_logfc$baseline  <- ident_base
  gsea_logfc$treatment <- ident_treat
  
  gsea_sign_p$baseline  <- ident_base
  gsea_sign_p$treatment <- ident_treat
  
  gsea_fc_p$baseline  <- ident_base
  gsea_fc_p$treatment <- ident_treat
  
  ### Return ---------------------------------------------------------------- 
  list(
    de           = de,
    gsea_logfc   = gsea_logfc,
    gsea_sign_p  = gsea_sign_p,
    gsea_fc_p    = gsea_fc_p
  )
}


