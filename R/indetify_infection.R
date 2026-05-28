#' Identify infected vs bystander vs unidentified cells in a list of Seurat objects
#'
#' @description
#' Assigns infection status ("infected", "bystander", or "unidentified") to each
#' cell in a list of Seurat objects based on expression of a viral gene pattern.
#' Optionally compares results between adjusted and original counts layers if the
#' latter exists in all objects.
#'
#' @param obj_list A named list of Seurat objects.
#' @param viral_pattern Character string (regex) to identify viral gene(s).
#' @param bystander_threshold Numeric. Cells with summed viral counts
#'   **strictly less than** this value are labeled "bystander" (default = 1).
#' @param infected_threshold Numeric. Cells with summed viral counts
#'   **strictly greater than** this value are labeled "infected" (default = 5).
#'   Must be >= bystander_threshold. Cells between the two thresholds are "unidentified".
#'
#' @details
#' Classification logic per cell:
#' \itemize{
#'   \item `viral counts > infected_threshold` → "infected"
#'   \item `viral counts < bystander_threshold` → "bystander"
#'   \item otherwise → "unidentified"
#' }
#' If an assay named `"original.counts"` exists in all objects, infection status
#' is computed for both adjusted and original counts, and summary/compare tables
#' are returned.
#'
#' @return A list containing:
#' \itemize{
#'   \item `updated_list`: list of Seurat objects with infection status columns.
#'   \item `viral_feat`: detected viral feature(s) per sample.
#'   \item `infection_summary`: table of infected/bystander/unidentified counts from adjusted data.
#'   \item `viral_feat_original`, `infection_summary_original`,
#'         `compare_infection`: included only if `original.counts` exists.
#' }
#'
#' @examples
#' \dontrun{
#' out <- identify_infection(
#'   singlet_list,
#'   viral_pattern = "ABB01533",
#'   bystander_threshold = 1,
#'   infected_threshold = 5
#' )
#' saveRDS(out, "infection(ABB01533_t0).rds")
#' }
#'
#' @export

identify_infection <- function(obj_list,
                               viral_pattern = "ABB01533",
                               bystander_threshold = 1,
                               infected_threshold = 5) {
  
  # Validate thresholds
  if (bystander_threshold > infected_threshold) {
    stop("`bystander_threshold` must be <= `infected_threshold`.")
  }
  
  # Helper: classify counts into three categories
  classify_infection <- function(total_v) {
    dplyr::case_when(
      total_v > infected_threshold  ~ "infected",
      total_v < bystander_threshold ~ "bystander",
      TRUE                          ~ "unidentified"
    )
  }
  
  # Helper: summarise a status column into a one-row data.frame
  summarise_status <- function(status_vec) {
    tbl <- table(status_vec)
    get_n <- function(k) if (k %in% names(tbl)) tbl[[k]] else 0
    data.frame(
      infected     = get_n("infected"),
      bystander    = get_n("bystander"),
      unidentified = get_n("unidentified")
    )
  }
  
  # 1) Viral features from adjusted counts
  viral_feat <- lapply(obj_list, function(obj) {
    rn   <- rownames(GetAssayData(obj, layer = "counts"))
    hits <- unlist(lapply(viral_pattern, function(p) grep(p, rn, value = TRUE)))
    unique(hits)
  })
  
  # 2) Assign infection status (adjusted counts)
  obj_list <- mapply(function(obj, feats) {
    if (length(feats) == 0) {
      obj$infection_status <- rep(NA_character_, ncol(obj))
    } else {
      v       <- GetAssayData(obj, layer = "counts")[feats, , drop = FALSE]
      total_v <- colSums(v)
      obj$infection_status <- classify_infection(total_v)
    }
    obj
  }, obj_list, viral_feat, SIMPLIFY = FALSE)
  
  # 3) Summary (adjusted)
  infection_summary <- do.call(rbind, lapply(obj_list, function(obj) {
    summarise_status(obj$infection_status)
  }))
  rownames(infection_summary) <- names(obj_list)
  
  # 4) Check for original counts assay
  has_original <- all(vapply(
    obj_list, function(obj) "original.counts" %in% names(obj@assays), logical(1)
  ))
  
  if (!has_original) {
    return(list(
      updated_list      = obj_list,
      viral_feat        = viral_feat,
      infection_summary = infection_summary
    ))
  }
  
  # 5) Viral features from original counts
  viral_feat_original <- lapply(obj_list, function(obj) {
    rn   <- rownames(GetAssayData(obj, assay = "original.counts", layer = "counts"))
    hits <- unlist(lapply(viral_pattern, function(p) grep(p, rn, value = TRUE)))
    unique(hits)
  })
  
  # 6) Assign infection status (original counts)
  obj_list <- mapply(function(obj, feats) {
    if (length(feats) == 0) {
      obj$infection_status_original <- rep(NA_character_, ncol(obj))
    } else {
      v       <- GetAssayData(obj, assay = "original.counts", layer = "counts")[feats, , drop = FALSE]
      total_v <- colSums(v)
      obj$infection_status_original <- classify_infection(total_v)
    }
    obj
  }, obj_list, viral_feat_original, SIMPLIFY = FALSE)
  
  # 7) Summary (original)
  infection_summary_original <- do.call(rbind, lapply(obj_list, function(obj) {
    summarise_status(obj$infection_status_original)
  }))
  rownames(infection_summary_original) <- names(obj_list)
  
  # 8) Compare table
  compare_infection <- do.call(rbind, lapply(names(obj_list), function(nm) {
    x <- obj_list[[nm]]
    data.frame(
      sample                = nm,
      adjusted_infected     = sum(x$infection_status          == "infected",     na.rm = TRUE),
      adjusted_bystander    = sum(x$infection_status          == "bystander",    na.rm = TRUE),
      adjusted_unidentified = sum(x$infection_status          == "unidentified", na.rm = TRUE),
      original_infected     = sum(x$infection_status_original == "infected",     na.rm = TRUE),
      original_bystander    = sum(x$infection_status_original == "bystander",    na.rm = TRUE),
      original_unidentified = sum(x$infection_status_original == "unidentified", na.rm = TRUE)
    )
  }))
  
  list(
    updated_list               = obj_list,
    viral_feat                 = viral_feat,
    infection_summary          = infection_summary,
    viral_feat_original        = viral_feat_original,
    infection_summary_original = infection_summary_original,
    compare_infection          = compare_infection
  )
}