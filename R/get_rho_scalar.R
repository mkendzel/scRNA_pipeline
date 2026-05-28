#' Extract a single rho; average when multiple values exist
#'
#' Returns a scalar contamination fraction (`rho`) from a `SoupChannel`.
#' Lookup order:
#' \enumerate{
#'   \item \code{sc$meta$rho}: if length 1, return that scalar.
#'   \item \code{sc$metaData$rho}: if vector of per-cell rhos, return the
#'         arithmetic mean with \code{na.rm = TRUE}.
#'   \item \code{sc$meta$rho}: if a vector, return the arithmetic mean with
#'         \code{na.rm = TRUE}.
#' }
#' If none yield a finite number, returns \code{NA_real_}.
#'
#' @param sc A \code{SoupChannel} object.
#'
#' @return Numeric scalar rho, or \code{NA_real_}.
#'
#' @examples
#' \dontrun{
#' rho_value <- get_rho_scalar(soup_results[["Mock-1"]])
#' }
#' @export
get_rho_scalar <- function(sc) {
  # 1) global rho stored in meta list
  r1 <- try(sc$meta$rho, silent = TRUE)
  if (!inherits(r1, "try-error") && length(r1) == 1 && is.finite(r1)) return(as.numeric(r1))
  
  # 2) per-cell rho in metaData -> arithmetic mean
  r2 <- try(sc$metaData$rho, silent = TRUE)
  if (!inherits(r2, "try-error") && length(r2) > 0) {
    r2m <- suppressWarnings(mean(as.numeric(r2), na.rm = TRUE))
    if (is.finite(r2m)) return(r2m)
  }
  
  # 3) vector in meta$rho -> arithmetic mean
  if (!inherits(r1, "try-error") && length(r1) > 1) {
    r1m <- suppressWarnings(mean(as.numeric(r1), na.rm = TRUE))
    if (is.finite(r1m)) return(r1m)
  }
  
  NA_real_
}
