#' Summarize total counts before and after SoupX adjustment
#'
#' For each Seurat object in a list, compute total RNA counts in the
#' `"original.counts"` assay and in the adjusted `"RNA"` assay. Returns a
#' markdown table summarizing raw totals and the adjustment ratio.
#'
#' @param adjusted_list Named list of Seurat objects that include both
#'   `"original.counts"` and adjusted `"RNA"` assays.
#'
#' @return A markdown-formatted table (via `knitr::kable`) showing:
#' \itemize{
#'   \item `sample`: name of the sample
#'   \item `original_sum`: total counts before adjustment
#'   \item `adjusted_sum`: total counts after adjustment
#'   \item `ratio`: adjusted_sum / original_sum
#' }
#'
#' @examples
#' \dontrun{
#' adjusted_count_summary(adjusted_list)
#' }
#'
#' @export
adjusted_count_summary <- function(adjusted_list) {
  results <- lapply(names(adjusted_list), function(nm) {
    sobj <- adjusted_list[[nm]]
    orig_sum <- sum(GetAssayData(sobj[["original.counts"]], layer = "counts"))
    adj_sum  <- sum(GetAssayData(sobj[["RNA"]], layer = "counts"))
    data.frame(
      sample = nm,
      original_sum = orig_sum,
      adjusted_sum = adj_sum,
      ratio = adj_sum / orig_sum,
      row.names = NULL
    )
  })
  
  counts_summary <- do.call(rbind, results)
  
  knitr::kable(
    counts_summary,
    format = "markdown",
    align = "c",
    caption = "SoupX adjusted RNA counts per sample"
  )
}