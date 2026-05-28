plot_hallmark_umap <- function(seurat_obj, hallmark_sets, hallmark_name, sample, infection,
                               method = c("zscore", "module")) {
  
  method <- match.arg(method)
  
  # Get gene set and subset cells
  genes <- hallmark_sets[[hallmark_name]]
  obj <- subset(seurat_obj, subset = sample_id == sample & infection_status == infection)
  genes <- genes[genes %in% rownames(obj)]
  
  # Clean up hallmark name for plot labels
  label <- gsub("^HALLMARK_", "", hallmark_name)
  label <- gsub("_", " ", label)
  label <- tools::toTitleCase(tolower(label))
  
  if (method == "zscore") {
    expr <- GetAssayData(obj, slot = "data")[genes, ]
    z_scores <- t(scale(t(as.matrix(expr))))
    obj$hallmark_score <- colMeans(z_scores, na.rm = TRUE)
    score_label <- "Z-score"
  } else {
    obj <- AddModuleScore(obj, features = list(genes), name = "hallmark_score")
    obj$hallmark_score <- obj$hallmark_score1
    score_label <- "Module Score"
  }
  
  # Plot on UMAP
  FeaturePlot(obj, features = "hallmark_score", order = TRUE, pt.size = 1) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = score_label) +
    ggtitle(NULL)
}