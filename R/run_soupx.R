#' Run SoupX contamination estimation across multiple Seurat objects
#'
#' For each sample in a named list of Seurat objects, constructs a
#' `SoupChannel` from raw and filtered 10x matrices, aligns gene sets,
#' assigns clusters, and attempts contamination estimation (`autoEstCont`)
#' using three increasingly permissive parameter sets:
#' \enumerate{
#'   \item Default parameters.
#'   \item `forceAccept = TRUE`.
#'   \item Relaxed thresholds: `tfidfMin = 0`, `soupQuantile = 0.3`, `forceAccept = TRUE`.
#' }
#' Records which attempt succeeded and the resulting `rho` value
#' for each sample. If all attempts fail, it stores the unadjusted
#' `SoupChannel` for that sample.
#'
#' @param data_list Named list of Seurat objects. Each must contain
#'   a `"RNA"` assay with a `"counts"` layer and a column `soup_group`
#'   defining clustering for SoupX.
#' @param dir_map Named character vector mapping each sample name in
#'   `data_list` to its corresponding 10x output directory.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `soup_results`: named list of `SoupChannel` objects (one per sample)
#'   \item `soup_summary`: data.frame with columns
#'     \code{sample}, \code{rho}, \code{method}, \code{status}
#' }
#'
#' @details
#' Aligns genes between raw and filtered matrices to prevent mismatch errors
#' (`tod` and `toc` must have the same genes). Each estimation attempt is
#' wrapped in `try()` to continue processing even if one sample fails.
#'
#' @examples
#' \dontrun{
#' res <- run_soupx(data_list, dir_map)
#' soup_results <- res$soup_results
#' soup_summary <- res$soup_summary
#' }
#'
#' @export
run_soupx <- function(data_list, dir_map) {
  soup_results <- setNames(vector("list", length(data_list)), names(data_list))
  soup_summary <- data.frame(sample = names(data_list),
                             rho = NA_real_, method = NA_character_,
                             status = "failed", stringsAsFactors = FALSE)
  
  for (nm in names(data_list)) {
    raw_mat <- Read10X(file.path(dir_map[[nm]], "outs", "raw_feature_bc_matrix"))
    counts  <- GetAssayData(data_list[[nm]], assay = "RNA", layer = "counts")
    genes_all <- intersect(rownames(raw_mat), rownames(counts))
    
    sc0 <- SoupChannel(raw_mat[genes_all, , drop = FALSE],
                       counts[genes_all, , drop = FALSE])
    sc0 <- setClusters(sc0, data_list[[nm]]$soup_group)
    
    # try default
    sc1 <- try(autoEstCont(sc0, doPlot = FALSE), silent = TRUE)
    if (!inherits(sc1, "try-error")) {
      soup_results[[nm]] <- sc1
      soup_summary[soup_summary$sample == nm, c("rho","method","status")] <- list(sc1$meta$rho, "default", "ok")
      next
    }
    
    # try forceAccept
    sc2 <- try(autoEstCont(sc0, doPlot = FALSE, forceAccept = TRUE), silent = TRUE)
    if (!inherits(sc2, "try-error")) {
      soup_results[[nm]] <- sc2
      soup_summary[soup_summary$sample == nm, c("rho","method","status")] <- list(sc2$meta$rho, "forceAccept", "ok")
      next
    }
    
    # try relaxed + forceAccept
    sc3 <- try(autoEstCont(sc0, doPlot = FALSE, tfidfMin = 0, soupQuantile = 0.3, forceAccept = TRUE),
               silent = TRUE)
    if (!inherits(sc3, "try-error")) {
      soup_results[[nm]] <- sc3
      soup_summary[soup_summary$sample == nm, c("rho","method","status")] <- list(sc3$meta$rho, "relaxed", "ok")
      next
    }
    
    # all attempts failed: keep pre-estimation object
    soup_results[[nm]] <- sc0
  }
  
  list(soup_results = soup_results, soup_summary = soup_summary)
}
