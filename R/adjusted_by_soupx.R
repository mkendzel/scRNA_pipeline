#' Adjust Seurat objects using estimated SoupX contamination
#'
#' For each sample in a list, apply `adjustCounts()` using corresponding
#' `SoupChannel` objects from `soup_results`. The original counts are preserved
#' in a new assay layer called `"original.counts"`.
#'
#' @param data_list Named list of Seurat objects that have been QC and filtered
#' @param soup_results Named list of `SoupChannel` objects with contamination
#'   estimates produced by `run_soupx()` or equivalent.
#'
#' @return A named list of Seurat objects where:
#' \itemize{
#'   \item The `"RNA"` assay’s `"counts"` layer contains the now decontaminated counts.
#'   \item A new layer`"original.counts"` stores the unadjusted matrix.
#' }
#'
#' @examples
#' \dontrun{
#' adjusted_list <- adjust_by_soupx(data_list, soup_results)
#' }
#'
#' @export
adjust_by_soupx <- function(data_list, soup_results) {
  adjusted_list <- lapply(names(data_list), function(nm) {
    sc   <- soup_results[[nm]]
    sobj <- data_list[[nm]]
    counts <- GetAssayData(sobj, assay = "RNA", layer = "counts")
    # 1) decontaminate counts
    out <- adjustCounts(sc, roundToInt = TRUE)
    # 2) keep original counts as a layer
    sobj[["original.counts"]] <- CreateAssayObject(counts = counts)
    # 3) replace counts layer with decontaminated matrix
    sobj <- SetAssayData(
      object = sobj, assay = "RNA", layer = "counts",
      new.data = out
    )
    sobj
  })
  names(adjusted_list) <- names(data_list)
  adjusted_list
}