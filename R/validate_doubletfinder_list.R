# --- Detect the DF classification column ------------------------------------
.df_col <- function(obj) {
  hits <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
  if (!length(hits)) return(NA_character_)
  hits[which.max(nchar(hits))]
}

# --- Main: summarize DF across a list, with optional expected rates ----------
# expected_rates: named numeric vector of expected % (0-100) per sample
#                 OR a data.frame with columns: sample, expected_pct
# tol_pp: absolute tolerance in percentage points for "within_range" check
validate_doubletfinder_list <- function(doublet_list,
                                        cluster_col = "seurat_clusters",
                                        expected_rates = NULL,
                                        tol_pp = 2) {
  stopifnot(is.list(doublet_list), length(doublet_list) > 0)
  
  # Normalize expected_rates into a named numeric vector of % (0–100)
  exp_vec <- NULL
  if (!is.null(expected_rates)) {
    if (is.data.frame(expected_rates)) {
      stopifnot(all(c("sample","expected_pct") %in% names(expected_rates)))
      exp_vec <- setNames(expected_rates$expected_pct, expected_rates$sample)
    } else if (is.numeric(expected_rates)) {
      stopifnot(!is.null(names(expected_rates)))
      exp_vec <- expected_rates
    } else {
      stop("expected_rates must be a named numeric vector or a data.frame with sample and expected_pct.")
    }
  }
  
  per_sample  <- list()
  per_cluster <- list()
  qc_medians  <- list()
  
  for (nm in names(doublet_list)) {
    obj <- doublet_list[[nm]]
    dfcol <- .df_col(obj)
    if (is.na(dfcol)) next
    
    md <- obj@meta.data
    if (!cluster_col %in% names(md)) md[[cluster_col]] <- Idents(obj)
    
    n_cells <- nrow(md)
    n_dbl   <- sum(md[[dfcol]] == "Doublet", na.rm = TRUE)
    pct_dbl <- n_dbl / n_cells  # percent
    
    exp_pct <- if (!is.null(exp_vec) && nm %in% names(exp_vec)) exp_vec[[nm]] else NA_real_
    within  <- if (is.finite(exp_pct)) abs(pct_dbl - exp_pct) <= tol_pp else NA
    
    per_sample[[nm]] <- data.frame(
      sample        = nm,
      df_col        = dfcol,
      n_cells       = n_cells,
      n_doublets    = n_dbl,
      pct_doublets  = pct_dbl,
      expected_pct  = exp_pct,
      diff_pp       = pct_dbl - exp_pct,
      within_range  = within,
      row.names = NULL
    )
    
    # per-cluster
    tab_all <- as.data.frame(table(md[[cluster_col]]), stringsAsFactors = FALSE)
    tab_dbl <- as.data.frame(table(md[[cluster_col]][md[[dfcol]] == "Doublet"]),
                             stringsAsFactors = FALSE)
    names(tab_all) <- c("cluster","n_cells")
    names(tab_dbl) <- c("cluster","n_doublets")
    m <- merge(tab_all, tab_dbl, by = "cluster", all.x = TRUE)
    m$n_doublets[is.na(m$n_doublets)] <- 0L
    m$pct_doublets <- 100 * m$n_doublets / m$n_cells
    m$sample <- nm
    m$df_col <- dfcol
    per_cluster[[nm]] <- m[, c("sample","cluster","n_cells","n_doublets","pct_doublets","df_col")]
    
    # QC medians by DF call
    if (all(c("nCount_RNA","nFeature_RNA") %in% names(md))) {
      med <- aggregate(cbind(nCount_RNA, nFeature_RNA) ~ DF_call + sample, data = transform(
        md, DF_call = md[[dfcol]], sample = nm
      ), FUN = median, na.rm = TRUE)
      qc_medians[[nm]] <- med
    }
  }
  
  list(
    per_sample  = do.call(rbind, per_sample),
    per_cluster = do.call(rbind, per_cluster),
    qc_medians  = if (length(qc_medians)) do.call(rbind, qc_medians) else NULL
  )
}
