#' Override rho for selected samples using the average from others
#' re-adjusts counts by said average rho
#'
#' Computes the mean contamination fraction (`rho`) from all samples **not** in
#' `targets`, then applies that mean to the `targets`. For each target:
#' restores the `"RNA"` counts from `"original.counts"`, sets the new rho on the
#' corresponding `SoupChannel`, re-runs `adjustCounts()`, and writes the result
#' back. Optionally updates `soup_summary`.
#'
#' @param adjusted_list Named list of Seurat objects that contain an
#'   `"original.counts"` assay and adjusted `"RNA"` counts.
#' @param soup_results Named list of `SoupChannel` objects aligned to `adjusted_list`.
#' @param targets Character vector of sample names to override, e.g. `c("Mock-1","Mock-2")`.
#' @param soup_summary Optional data.frame with columns `sample`, `rho`, `method`, `status`.
#'
#' @return A list with updated objects:
#' \itemize{
#'   \item `adjusted_list`: after re-adjustment for `targets`
#'   \item `soup_results`: with overridden rho for `targets`
#'   \item `soup_summary`: updated if supplied
#'   \item `avg_rho`: scalar mean rho used for override
#' }
#'
#' @examples
#' \dontrun{
#' res <- override_targets_by_avg_rho(
#'   adjusted_list = adjusted_list,
#'   soup_results  = soup_results,
#'   targets       = c("Mock-1","Mock-2"),
#'   soup_summary  = soup_summary
#' )
#' adjusted_list <- res$adjusted_list
#' soup_results  <- res$soup_results
#' soup_summary  <- res$soup_summary
#' res$avg_rho
#' }
#' @export
override_targets_by_avg_rho <- function(adjusted_list, soup_results, targets, soup_summary = NULL) {
  others  <- setdiff(names(soup_results), targets)
  avg_rho <- mean(vapply(others, function(nm) get_rho_scalar(soup_results[[nm]]), numeric(1)), na.rm = TRUE)
  stopifnot(is.finite(avg_rho))
  
  for (nm in intersect(targets, names(adjusted_list))) {
    sobj <- adjusted_list[[nm]]
    sc   <- soup_results[[nm]]
    
    # reset from saved original counts
    orig <- GetAssayData(sobj[["original.counts"]], layer = "counts")
    sobj <- SetAssayData(sobj, assay = "RNA", layer = "counts", new.data = orig)
    
    # override rho and re-adjust
    sc_fix <- setContaminationFraction(sc, contFrac = avg_rho, forceAccept = TRUE)
    out    <- adjustCounts(sc_fix, roundToInt = TRUE)
    sobj   <- SetAssayData(sobj, assay = "RNA", layer = "counts", new.data = out)
    
    adjusted_list[[nm]] <- sobj
    soup_results[[nm]]  <- sc_fix
  }
  
  if (!is.null(soup_summary)) {
    idx <- soup_summary$sample %in% targets
    soup_summary$rho[idx]    <- avg_rho
    soup_summary$method[idx] <- "avg_rho_override"
    soup_summary$status[idx] <- "ok"
  }
  
  list(
    adjusted_list = adjusted_list,
    soup_results  = soup_results,
    soup_summary  = soup_summary,
    avg_rho       = avg_rho
  )
}
