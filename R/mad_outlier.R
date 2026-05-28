#' Identify MAD outliers in Seurat metadata
#'
#' Flags cells whose given metric is more than `nmads` MADs from the median.
#'
#' @param sobj Seurat object
#' @param metric Character. Column name in `sobj@meta.data`.
#' @param nmads Numeric. Number of MADs from the median used as the cutoff.
#'
#' @return Logical vector marking outlier cells.
#' @export

mad_outlier <- function(sobj, metric, nmads){
  M <- sobj@meta.data[[metric]]
  median_M <- median(M, na.rm = TRUE)
  mad_M <- mad(M, na.rm = TRUE)
  outlier <- (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))
  return(outlier)
}