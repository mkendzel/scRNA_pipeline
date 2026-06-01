# ---- Load Library, Functions, and objects ----
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(clustree)

#All project functions
invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))

#QC and SoupX adjusted data
adjusted_list <- load_checkpoint("allData_MADqc_SoupX", dir = "projects/YOUR_PROJECT_NAME/data/processed")

# ---- Predict Multiplet (doublet) rate for DoubletFinder ----
# Load Samples from experiment
new_samples <- data.frame(
  sample = names(adjusted_list),
  CellsRecovered = vapply(adjusted_list, ncol, numeric(1))
)

# helper function to predict multiplet rate using pre-trained model
results <- predict_multiplet_rate(new_samples)

# View results
knitr::kable(results$detailed_predictions, format = "markdown")
print(results$detailed_predictions)
print(results$quick_lookup_matrix)

# Assign predicted multiplet % per sample for later
expected_df <- data.frame(sample = results$detailed_predictions$sample, expected_pct = results$detailed_predictions$upr)

# ---- Re-normalize, scale, PCA, and Find Variable features on new count matrix ----
pre_out <- lapply(adjusted_list, prep_precluster, nfeatures = 2000)

# Access results
adjusted_list <- lapply(pre_out, `[[`, "obj")
dims_use     <- lapply(pre_out, `[[`, "dims_use")

# ================================#
# Code to run DoubletFinder and adjust counts
# ================================#
# DoubletFinder requires clustering at a fine resolution to capture heterogeneity, but not so high that clusters are mostly doublets.
# Scan resolution for stable higher resolution clusters 
# run across samples
res_out <- mapply(
  function(obj, dims) scan_resolutions_for_clustree(obj, dims, res = seq(0.2, 1.2, 0.1)),
  adjusted_list, dims_use, SIMPLIFY = FALSE
)

# update objects (contain all res columns for clustree)
adjusted_list <- lapply(res_out, `[[`, "obj")


# save clustree plots, one per sample
mapply(function(obj, nm) {
  p <- clustree::clustree(obj)   # returns a ggplot
  ggplot2::ggsave(
    filename = file.path("YOUR_Project_NAME/diagnostics/DoubletFinder/clustertree", paste0(nm, "_clustree.png")),
    plot = p, width = 8, height = 8, dpi = 300
  )
}, adjusted_list, names(adjusted_list))


#Save umaps per resolution
save_group_umaps(adjusted_list, dims_use, out_dir = "YOUR_Project_NAME/diagnostics/DoubletFinder/umaps", pt.size = 0.3)


# ADJUST: Select resolution per sample after reviewing clustree and UMAPs above
chosen_res <- c(
  "Condition_1" = 0.8,  # ADJUST
  "Condition_2" = 0.8   # ADJUST
  # Add an entry for each sample using its condition name
)

adjusted_list <- mapply(function(obj, dims) {
  if (!"umap" %in% names(obj@reductions))
    obj <- RunUMAP(obj, dims = dims)
  obj
}, adjusted_list, dims_use, SIMPLIFY = FALSE)

# ---- DoubletFinder Parameter sweep ----
sweep_results <- mapply(function(obj, group) {
  res <- chosen_res[[group]]
  res_col <- paste0("RNA_snn_res.", res)

  Idents(obj) <- obj@meta.data[[res_col]]
  pcs <- dims_use[[group]]

  paramSweep(obj, PCs = pcs, sct = FALSE)
}, adjusted_list, names(adjusted_list), SIMPLIFY = FALSE)
beepr::beep(1)
sweep_stats <- lapply(sweep_results, function(obj){
  summarizeSweep(obj, GT = FALSE)
})

bcmvn_list <- lapply(sweep_stats, find.pK)

#Find best pK value
pK_list <- lapply(bcmvn_list, function(df) {
  pk_val <- df %>%
    filter(BCmetric == max(BCmetric)) %>%
    slice(1) %>%                 # handle ties
    pull(pK)
  as.numeric(as.character(pk_val))  # convert factor/char to numeric
})


annotations_list <- mapply(function(obj, group) {
  res <- chosen_res[[group]]
  res_col <- paste0("RNA_snn_res.", res)

  if (res_col %in% colnames(obj@meta.data)) {
    obj@meta.data[[res_col]]
  } else {
    stop(sprintf("Resolution %s not found for %s", res, group))
  }
}, adjusted_list, names(adjusted_list), SIMPLIFY = FALSE)


homotypic_list <- lapply(annotations_list, modelHomotypic)


# ---- compute expected doublets per object ----
nExp_poi <- setNames(
  vapply(names(adjusted_list), function(nm) {
    upr_val <- results$detailed_predictions$upr[
      match(nm, results$detailed_predictions$sample)
    ]
    round(upr_val * nrow(adjusted_list[[nm]]@meta.data))
  }, numeric(1)),
  names(adjusted_list)
)

nExp_poi_adj_list <- setNames(
  lapply(names(nExp_poi), function(nm) {
    nExp <- nExp_poi[[nm]]
    homo <- homotypic_list[[nm]]
    round(nExp * (1 - homo))
  }),
  names(nExp_poi)
)


# ---- Run Doublet finder ----
doublet_list <- setNames(lapply(names(adjusted_list), function(nm) {
  obj  <- adjusted_list[[nm]]
  pK   <- pK_list[[nm]]
  nExp <- nExp_poi_adj_list[[nm]]
  pcs  <- dims_use[[nm]]          # object-specific PCs

  doubletFinder(
    obj,
    PCs  = pcs,
    pN   = 0.25,
    pK   = pK,
    nExp = nExp,
    sct  = FALSE
  )
}), names(adjusted_list))


#Save
save_checkpoint(doublet_list, "allData_MADqc_SoupX_Doublet", dir = "projects/YOUR_PROJECT_NAME/data/processed")

# Optional load
# doublet_list <- load_checkpoint("allData_MADqc_SoupX_Doublet", dir = "projects/YOUR_PROJECT_NAME/data/processed")


doublet_umap_plots <- lapply(names(doublet_list), function(nm) {
  obj <- doublet_list[[nm]]
  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)[1]

  dims_val <- length(dims_use[[nm]])
  res_val  <- chosen_res[[nm]]

  title_text <- paste0(nm, " | ", df_col, " | dims=", dims_val, " | res=", res_val)

  p <- DimPlot(obj, reduction = "umap", group.by = df_col) + ggtitle(title_text)

  filename <- paste0("UMAP_", nm, "_doublets_dims", dims_val, "_res", res_val, ".png")
  ggsave(file.path("projects/YOUR_PROJECT_NAME/diagnostics/DoubletFinder/doublets", filename), plot = p, width = 6, height = 5, dpi = 300)
  p
})



# ---- Validate DoubletFinder ----
doublet_list <- set_idents_by_res(doublet_list, chosen_res)
confirm_idents(doublet_list, chosen_res)

#Validation table for Doublet finder
res <- validate_doubletfinder_list(
  doublet_list,
  expected_rates = expected_df,
  tol_pp = 2
)


# Examples
writeLines(knitr::kable(res$per_sample[order(-res$per_sample$pct_doublets), ],
                        format = "pipe", row.names = FALSE, digits = 4))

df <- res$per_sample[order(-res$per_sample$pct_doublets), ]
cat(knitr::kable(df, format = "pipe"))
knitr::kable(subset(res$per_cluster, pct_doublets >= 10)[order(-subset(res$per_cluster, pct_doublets >= 10)$pct_doublets), ], format = "pipe", digits = 2)

knitr::kable(res$qc_medians, format = "pipe", digits = 1)

# ---- Subset out singlets -----
singlet_list <- lapply(doublet_list, function(obj) {
  df_col <- tail(grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE), 1)
  keep   <- rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Singlet"]
  subset(obj, cells = keep)
})
names(singlet_list) <- names(doublet_list)

#Save
save_checkpoint(singlet_list, "singletData", dir = "projects/YOUR_PROJECT_NAME/data/processed")

#Plots
umap_singlet_plots <- lapply(names(singlet_list), function(nm) {
  obj <- singlet_list[[nm]]

  dims_val <- length(dims_use[[nm]])
  res_val  <- chosen_res[[nm]]

  title_text <- paste0(nm, " | dims=", dims_val, " | res=", res_val)

  p <- DimPlot(obj, reduction = "umap") + ggtitle(title_text)

  filename <- paste0("UMAP_", nm, "_singlets_dims", dims_val, "_res", res_val, ".png")
  ggsave(
    filename = file.path("projects/YOUR_PROJECT_NAME/diagnostics/DoubletFinder/singlets", filename),
    plot = p,
    width = 6, height = 5, dpi = 300
  )
  p
})



#See how many doublets calculated
removed_summary <- data.frame(
  sample = names(doublet_list),
  before = sapply(doublet_list, ncol),
  after  = sapply(singlet_list, ncol)
)

removed_summary$removed <- removed_summary$before - removed_summary$after

knitr::kable(removed_summary, format = "markdown")



res_dims_tbl <- data.frame(
  sample = names(singlet_list),
  dims_used = vapply(names(singlet_list), function(nm) length(dims_use[[nm]]), numeric(1)),
  res_used  = vapply(names(singlet_list), function(nm) chosen_res[[nm]], numeric(1))
)

# Print as markdown table
knitr::kable(
  res_dims_tbl,
  caption = "Dimensions and Resolution Used per Sample",
  format = "markdown"
)
