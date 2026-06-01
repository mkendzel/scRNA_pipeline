# ==============================================#
# Clustering and UMAP Visualization of Singlet Data per sample#
# Includes optional regression out of viral/pathogen features if applicable #
# Final per sample qc checkpoint before merger and integration: 
# check umaps across res and samples for batch effects, failed library, contamination, bad doublet removal, etc. before proceeding to final clustering and downstream analyses #
# ==============================================#

# ---- Load Library, Functions, and objects ----
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(clustree)

#All project functions
invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))

#Singlet Data
singlet_list <- load_checkpoint("singletData", dir = "projects/YOUR_PROJECT_NAME/data/processed")

#---- Prepare Singlet data for clustering ----

# OPTIONAL (viral/pathogen experiments only): Calculate the percentage of reads
# mapping to viral or pathogen features. Update the pattern to match your feature
# naming convention (e.g., "^HIV-", "^SARS-"). Remove this entire block if your
# experiment does not include viral/pathogen features.
singlet_list <- lapply(singlet_list, function(obj) {
  obj[["percent.viral"]] <- PercentageFeatureSet(
    obj,
    pattern = "^YOUR_VIRAL_FEATURE_PATTERN"  # ADJUST: match your viral feature names
  )
  return(obj)
})
# END OPTIONAL viral block

# prep for clustering
# OPTIONAL (viral/pathogen experiments only): Remove "percent.viral" from
# vars.to.regress and remove the exclude_features argument entirely if you did
# not run the percent.viral block above.
pre_out <- lapply(singlet_list, prep_merged_precluster,
                  nfeatures = 2000,
                  vars.to.regress = c("nCount_RNA", "percent.mt", "percent.viral"),  # ADJUST: remove percent.viral if not applicable
                  exclude_features = c("YOUR_FEATURE_PATTERN"))                       # OPTIONAL: remove if no viral/pathogen features to exclude

# Access results
singlet_list <- lapply(pre_out, `[[`, "obj")
dims_use     <- lapply(pre_out, `[[`, "dims_use")

# ---- Scan Resolutions for proper clustering ----

# run across samples
res_out <- mapply(
  function(obj, dims) scan_resolutions_for_clustree(obj, dims, res = seq(0.2, 1.2, 0.1)),
  singlet_list, dims_use, SIMPLIFY = FALSE
)

# update objects (contain all res columns for clustree)
singlet_list <- lapply(res_out, `[[`, "obj")


# save clustree plots, one per sample
mapply(function(obj, nm) {
  p <- clustree::clustree(obj)   # returns a ggplot
  ggplot2::ggsave(
    filename = file.path("projects/YOUR_PROJECT_NAME/diagnostics/sample_clustering/clustertree", paste0(nm, "_clustree.png")),
    plot = p, width = 8, height = 8, dpi = 300
  )
}, singlet_list, names(singlet_list))


# ---- Cluster based on chosen resolution ----

#Save all resolution umaps for selection
save_group_umaps(singlet_list, out_dir = "projects/YOUR_PROJECT_NAME/diagnostics/sample_clustering/umaps", dims_use, pt.size = 0.3)

# ADJUST: Replace "Condition_1" with a sample name from your data to preview a specific UMAP
preview_umap("Condition_1", res = 0.4, pt.size = 1)


# ADJUST: Select resolution per sample after reviewing clustree and UMAPs above
chosen_res <- c(
  "Condition_1" = 0.2,  # ADJUST
  "Condition_2" = 0.2   # ADJUST
  # Add an entry for each sample
)

knitr::kable(as.data.frame(chosen_res), col.names = c("Sample", "Chosen Resolution"))

singlet_list <- mapply(function(obj, res) {
  res_col <- paste0("RNA_snn_res.", res)
  if (res_col %in% colnames(obj@meta.data)) {
    Idents(obj) <- obj@meta.data[[res_col]]
  } else {
    stop(paste("Resolution", res, "not found for", unique(obj$sample_id)))
  }
  obj
}, singlet_list, chosen_res[names(singlet_list)], SIMPLIFY = FALSE)

# ---- Run final UMaps ----

singlet_list <- mapply(function(obj, dims) {
  if (!"umap" %in% names(obj@reductions))
    obj <- RunUMAP(obj, dims = dims)
  obj
}, singlet_list, dims_use, SIMPLIFY = FALSE)

plots <- mapply(function(obj, nm) {
  res_val <- chosen_res[[nm]]
  DimPlot(obj, reduction = "umap", label = TRUE) +
    ggtitle(sprintf("%s (res=%.2f)", nm, res_val))
}, singlet_list, names(singlet_list), SIMPLIFY = FALSE)

# save as before
mapply(function(p, nm)
  ggsave(file.path("projects/YOUR_PROJECT_NAME/results/graphs/sample_umaps", paste0(nm, "_UMAP.png")),
         plot = p, width = 6, height = 5, dpi = 300),
  plots, names(plots))

####
# ---- Save Final edited objects ----
####
save_checkpoint(singlet_list, "singletData_clustered", dir = "projects/YOUR_PROJECT_NAME/data/processed")
