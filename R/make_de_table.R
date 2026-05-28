#' @title Differential expression for infection comparison
#' @description
#' Runs Seurat::FindMarkers() to identify genes differentially expressed
#' between infected and bystander cells within specified treatment groups.
#' Filters for significance (p_val_adj < 0.05), ranks by log2 fold change,
#' and saves the top 20 upregulated and downregulated genes as both .rds
#' and markdown (.md) tables for documentation.
#'
#' @param obj A Seurat object containing integrated and annotated data.
#' @param samples Character vector of sample IDs (e.g. c("WNV-TX-1","WNV-TX-2","WNV-TX-3")).
#' @param group_label Short label for output file naming (e.g. "WNV-TX").
#' @return A data.frame of the top significantly up- and down-regulated genes.
#' @examples
#' tx_samples  <- c("WNV-TX-1","WNV-TX-2","WNV-TX-3")
#' mad_samples <- c("WNV-MAD-1","WNV-MAD-2","WNV-MAD-3")
#' de_tx  <- make_de_table(combined, tx_samples,  "WNV-TX")
#' de_mad <- make_de_table(combined, mad_samples, "WNV-MAD")
#' @export


make_de_table <- function(obj, samples, group_label) {
  subset_obj <- subset(obj, cells = WhichCells(obj, expression = sample %in% samples))
  Idents(subset_obj) <- "infection_status"
  
  markers <- FindMarkers(subset_obj,
                         ident.1 = "infected",
                         ident.2 = "bystander")
  
  markers <- subset(markers, p_val_adj < 0.05)
  markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
  
  top_up   <- head(markers, 20)
  top_down <- tail(markers, 20)
  
  top_combined <- rbind(
    data.frame(Direction = "Upregulated", Gene = rownames(top_up), top_up),
    data.frame(Direction = "Downregulated", Gene = rownames(top_down), top_down)
  )
  
  # save outputs
  sink(paste0("DEgenes(", group_label, "_infection_sig).md"))
  print(kable(top_combined, digits = 3, align = "c",
              caption = paste0("Top significantly up- and down-regulated genes in ",
                               group_label, " infected vs. bystander cells (p_val_adj < 0.05)")))
  sink()
  
  return(top_combined)
}