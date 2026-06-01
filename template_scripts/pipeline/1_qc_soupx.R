# ---- Library and functions ----
library(Seurat)
library(SoupX)
library(ggplot2)

#All project functions
invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))

# ---- Discover samples ----
# ADJUST: Set the path to your Cell Ranger output directory
data_dir <- "path/to/your/cellranger/output"
dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)

# ADJUST: Update pattern to match your Cell Ranger output directory naming convention
sample_dirs <- dirs[grepl("YOUR_SAMPLE_DIR_PATTERN", basename(dirs))]
sample_dirs <- unique(sample_dirs)

# IDs
bn  <- basename(sample_dirs)
# ADJUST: Update regex to extract sample ID from directory name
sid <- sub("^.*-(s\\d+).*$", "\\1", bn)

# ADJUST: Map extracted sample IDs to condition/treatment names
sample_lookup <- c(
  s001 = "Condition_1",
  s002 = "Condition_2"
  # Add entries for all samples
)

sample_ids <- unname(sample_lookup[sid])   # character vector of conditions
dir_map    <- setNames(sample_dirs, sample_ids)


# ---- Import and pre-process ----

data_list <- Map(function(d, id) {
  pp_dir(
    sample_dir   = d,
    sample_id    = id,
    min_counts   = 800,    # ADJUST: minimum UMI counts per cell
    min_features = 500,    # ADJUST: minimum genes detected per cell
    mad_counts   = 5,      # ADJUST: MAD multiplier for count filtering
    mad_features = 5,      # ADJUST: MAD multiplier for feature filtering
    mad_mito     = 3       # ADJUST: MAD multiplier for mitochondrial filtering
  )
}, dir_map, names(dir_map))


# ---- Graph Data to confirm Data Filtering and QC worked ----

#Full datalist
plot_counts_features(data_list)

# ADJUST: Set upper count threshold based on your data distribution
data_list <- lapply(data_list, function(obj) subset(obj, subset = nCount_RNA <= 50000))

#Table of remaining cells
cell_table <- data.frame(
  Cells  = vapply(data_list, ncol, numeric(1))
)

knitr::kable(cell_table, format = "markdown")

#Save checkpoint of pre-processed data before SoupX
save_checkpoint(data_list, "allData_MADqc", dir = "YOUR_PROJECT_NAME/data/processed")

# Optional load
# data_list <- load_checkpoint("allData_MADqc", dir = "YOUR_PROJECT_NAME/data/processed")

# ================================#
# Code to run SoupX and adjust counts
# ================================#
# Without hashing or known empty droplets, we will use SoupX for estimating ambient RNA contamination. 
# This is a common approach when no additional information is available to define empty droplets.
# ---- Prepare data set for SoupX ----
pre_out <- lapply(data_list, prep_precluster, nfeatures = 2000)

# Access results
data_list <- lapply(pre_out, `[[`, "obj")
dims_use     <- lapply(pre_out, `[[`, "dims_use")

# ---- Scan resolution for proper clustering used in SoupX ----
# SoupX runs better on lower resolution clusters that capture broad cell types, so select the lowest stable resolution
# run across samples
# Adjust resolution range as needed based on your data and clustree results
res_out <- mapply(
  function(obj, dims) scan_resolutions_for_clustree(obj, dims, res = seq(0.2, 1.2, 0.1)),
  data_list, dims_use, SIMPLIFY = FALSE
)

# update objects (contain all res columns for clustree)
data_list <- lapply(res_out, `[[`, "obj")

# save clustree plots, one per sample
mapply(function(obj, nm) {
  p <- clustree::clustree(obj)   # returns a ggplot
  ggplot2::ggsave(
    filename = file.path("YOUR_Project_NAME/diagnostics/SoupX/clustertree", paste0(nm, "_clustree.png")),
    plot = p, width = 8, height = 8, dpi = 300
  )
}, data_list, names(data_list))

#Save umaps per resolution
save_group_umaps(data_list, dims_use, out_dir = "YOUR_Project_NAME/diagnostics/SoupX/umaps", pt.size = 0.3)

# View clustree and UMAPs to select resolution for each sample, then update chosen_res below
# ADJUST: Select resolution per sample after reviewing clustree and UMAPs above
chosen_res <- c(
  "Condition_1" = 0.3,  # ADJUST
  "Condition_2" = 0.3   # ADJUST
  # Add an entry for each sample using its condition name from sample_lookup
)

# ---- Set Cluster groups for SoupX ----
data_list <- mapply(function(obj, res) {
  res_col <- paste0("RNA_snn_res.", res)
  if (res_col %in% colnames(obj@meta.data)) {
    obj$soup_group <- obj@meta.data[[res_col]]   # save chosen resolution clusters
    Idents(obj) <- obj$soup_group
  } else {
    stop(paste("Resolution", res, "not found for", unique(obj$sample_id)))
  }
  obj
}, data_list, chosen_res[names(data_list)], SIMPLIFY = FALSE)


# ---- Explore RHO per group ----
# dir_map is created in "import and pre-process" and "Discover samples" sections
set.seed(1234)
res <- run_soupx(data_list, dir_map)

soup_results <- res$soup_results
soup_summary <- res$soup_summary
knitr::kable(soup_summary, format = "markdown")

# ---- Adjust Counts based on estimated Ambient RNA ----

adjusted_list <- adjust_by_soupx(data_list, soup_results)

#Summary information
adjusted_count_summary(adjusted_list)

#Save
save_checkpoint(adjusted_list, "allData_MADqc_SoupX",  dir = "YOUR_PROJECT_NAME/data/processed")
# See all versions saved
list_checkpoints("allData_MADqc_SoupX", dir = "YOUR_PROJECT_NAME/data/processed")

# Optional load
# adjusted_list <- load_checkpoint("allData_MADqc_SoupX", dir = "YOUR_PROJECT_NAME/data/processed")
