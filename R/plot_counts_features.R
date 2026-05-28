#' Plot nCount_RNA vs nFeature_RNA for Seurat objects
#' @param data Named list of Seurat objects or a single Seurat object
#' @return Patchwork object (or single ggplot) showing LM fit and Pearson R
#' @export
plot_counts_features <- function(data, fig_dir = "figures/counts_features") {
  if (inherits(data, "Seurat")) data <- list(Object = data)
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  
  safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  
  plots <- lapply(names(data), function(nm) {
    obj <- data[[nm]]
    df  <- obj[[]]
    r   <- cor(df$nCount_RNA, df$nFeature_RNA, method = "pearson", use = "complete.obs")
    lab <- paste0("R = ", round(r, 3))
    
    p <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red") +
      ggplot2::annotate("text", x = Inf, y = Inf, label = lab,
                        hjust = 1.1, vjust = 2, size = 4) +
      ggplot2::ggtitle(nm)
    
    ggplot2::ggsave(
      filename = paste0("counts_vs_features_", safe_name(nm), ".png"),
      path = fig_dir,
      plot = p,
      width = 6, height = 5, dpi = 150
    )
    p
  })
  
  if (length(plots) == 1) plots[[1]] else patchwork::wrap_plots(plots)
}