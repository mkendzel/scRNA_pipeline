# Format Seurat's resolution suffix: 0.70 -> "0.7", 1.0 -> "1"
.fmt_res <- function(x) sub("\\.?0+$", "", as.character(x))

# Set Idents() for each object using your chosen_res mapping
set_idents_by_res <- function(obj_list, chosen_res) {
  mapply(function(obj, nm) {
    if (!nm %in% names(chosen_res)) return(obj)
    col <- paste0("RNA_snn_res.", .fmt_res(chosen_res[[nm]]))
    if (col %in% colnames(obj@meta.data)) Idents(obj) <- obj[[col]][, 1]
    obj
  }, obj_list, names(obj_list), SIMPLIFY = FALSE)
}

# Quick confirmation summary
confirm_idents <- function(obj_list, chosen_res) {
  do.call(rbind, lapply(names(obj_list), function(nm) {
    obj <- obj_list[[nm]]
    lev <- levels(Idents(obj))
    data.frame(
      sample    = nm,
      used_res  = if (nm %in% names(chosen_res)) chosen_res[[nm]] else NA_real_,
      n_levels  = length(lev),
      levels_preview = paste(utils::head(lev, 6), collapse = ", "),
      stringsAsFactors = FALSE
    )
  }))
}