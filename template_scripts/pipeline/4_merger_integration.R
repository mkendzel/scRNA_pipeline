# ---- Library and original data import ----
library(Seurat)
library(tidyverse)
library(pbapply)
library(harmony)
library(knitr)
library(cowplot)


# load
singlet_list <- load_checkpoint("singletData_clustered", dir = "projects/YOUR_PROJECT_NAME/data/processed")

# ---- Project functions ----

invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))
####

# ---- Identify Infected and Bystander Cells ----
# ADJUST: Update viral_pattern to match your pathogen feature naming convention
# Remove this section entirely if your experiment has no viral/pathogen features
out <- identify_infection_normalized(singlet_list,
                                     viral_pattern = c("YOUR_VIRAL_FEATURE_PATTERN"),  # ADJUST
                                     bystander_threshold = 0,
                                     infected_threshold = 1)

#Data summary and exploration
out$infection_summary
out$compare_infection
singlet_list <- out$updated_list
knitr::kable(out$compare_infection)

save_checkpoint(singlet_list, "singletData_infectionLabeled", dir = "data")

# ---- Infection load histograms ----
# ADJUST: Update sample factor levels to match your condition names and desired order
viral_feat <- out$viral_feat
viral_expr <- do.call(rbind, lapply(names(singlet_list), function(nm) {
  obj <- singlet_list[[nm]]
  idx <- obj$infection_status == "infected"
  if (!any(idx)) return(NULL)
  feat <- viral_feat[[nm]]
  v <- GetAssayData(obj, layer = "counts")[feat, idx, drop = TRUE]
  data.frame(sample = nm, expression = as.numeric(v))
}))

# ADJUST: Set factor levels to control facet order
viral_expr$sample <- factor(viral_expr$sample, levels = c(
  "Condition_1",
  "Condition_2"
  # Add all infected conditions in desired display order
))


p <- ggplot(viral_expr, aes(x = expression)) +
  geom_histogram(bins = 50, fill = "grey70", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Viral gene counts (log1p)", y = "Cells") +
  facet_wrap(~ sample, scales = "free_y")

p

ggsave("projects/YOUR_Project_name/results/graphs/viral_expression_histogram.png", plot = p, width = 10, height = 6, dpi = 300)


# ---- Merge data for standalone and integration ----

# Merge
merged.seurat <- merge(singlet_list[[1]], singlet_list[-1], add.cell.ids = names(singlet_list))

# ADJUST: Rename sample_id values if needed to shorter display names (remove if not needed)
# merged.seurat$sample_id <- dplyr::recode(
#   merged.seurat$sample_id,
#   "Original Long Name 1" = "Short Name 1",
#   "Original Long Name 2" = "Short Name 2"
# )

# Precluster prep
# ADJUST: Update vars.to.regress and exclude_features to match your experiment
pre_out <- prep_merged_precluster(merged.seurat,
                                  nfeatures = 2000,
                                  vars.to.regress = c("nCount_RNA", "percent.mt", "percent.viral"),  # ADJUST
                                  exclude_features = c("YOUR_FEATURE_PATTERN"))                       # ADJUST or remove: example could exclude viral/pathogen features if they are not relevant for clustering
dims_use <- pre_out$dims_use

#Scan resolutions for cluster tree
res_out <- scan_resolutions_for_clustree(pre_out$obj, dims_use, res = seq(0.2, 1.2, 0.1))

# Save umaps to select resolution
save_merged_umaps(res_out$obj, out_dir = "projects/YOUR_Project_name/diagnostics/merging",
                  dims = dims_use, pt.size = 0.10)

#Final clustering
merged.seurat <- res_out$obj
merged.seurat <- RunUMAP(merged.seurat, dims = dims_use)
merged.seurat <- FindNeighbors(merged.seurat, dims = dims_use)
merged.seurat <- FindClusters(merged.seurat, resolution = 0.13)  # ADJUST: set after reviewing UMAPs above
merged.seurat <- JoinLayers(merged.seurat)

# ADJUST: Update viral_feats to match your pathogen feature naming convention
# Set viral_feats <- c() and merged.seurat$viral_load <- NA if no viral features
viral_feats <- c("YOUR_VIRAL_FEATURE_PATTERN")  # ADJUST or set to c()
if (length(viral_feats) > 0) {
  vmat <- GetAssayData(merged.seurat, layer = "data")[viral_feats, , drop = FALSE]
  merged.seurat$viral_load <- colSums(vmat)
} else {
  merged.seurat$viral_load <- NA
}

# ---- Merged Graphs ----
sample   <- DimPlot(merged.seurat, reduction = "umap", group.by = "sample_id")
clusters <- DimPlot(merged.seurat, reduction = "umap", group.by = "seurat_clusters")
infected <- DimPlot(merged.seurat, reduction = "umap", group.by = "infection_status")

sample | clusters
infected | sample
clusters | infected
clusters | infected | sample

#Save graphs

# Save individual plots
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/sample.png",   sample,   width = 6, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/clusters.png", clusters, width = 6, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/infected.png", infected, width = 6, height = 5)

# Save combined plots
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/sample_clusters.png",        sample | clusters,           width = 12, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/infected_sample.png",        infected | sample,           width = 12, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/clusters_infected.png",      clusters | infected,         width = 12, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/clusters_sample_infected.png", clusters | sample | infected, width = 12, height = 5)



# ---- Save Merged Object ----

save_checkpoint(merged.seurat, "merged_clustered_normData", dir = "projects/YOUR_PROJECT_NAME/data/processed")
# Optional load
# merged.seurat <- load_checkpoint("merged_clustered_normData", dir = "projects/YOUR_PROJECT_NAME/data/processed")

# ================================================================= #
# ---- Integrate experimental data for downstream comparisons ----
# ================================================================= #
#      Optional Step to correct for batch effects if present         #

# Harmony batch correction
# use the metadata column you want to correct by: here "sample" or whatever column contains your sample/condition labels. Adjust if needed.
integrated.harmony <- RunHarmony(merged.seurat, group.by.vars = "sample_id")

# Scan resolutions for clustering
res_out <- scan_resolutions_harmony(integrated.harmony, dims_use, res = seq(0.2, 1.2, 0.1))

save_integrated_umaps(res_out$obj, out_dir = "[projects/YOUR_PROJECT_NAME/diagnostics/integration/umap_by_res",
                      dims = dims_use, pt.size = 0.15)



# Select resolution and save to new object
integrated.harmony <- res_out$obj
integrated.harmony <- RunUMAP(integrated.harmony, reduction = "harmony", dims = dims_use)
integrated.harmony <- FindNeighbors(integrated.harmony, reduction = "harmony", dims = dims_use)
integrated.harmony <- FindClusters(integrated.harmony, resolution = 0.2)  # ADJUST


# Example plots
sample   <- DimPlot(integrated.harmony, reduction = "umap", group.by = "sample_id")
clusters <- DimPlot(integrated.harmony, reduction = "umap", group.by = "seurat_clusters")
infected <- DimPlot(integrated.harmony, reduction = "umap", group.by = "infection_status")

sample | clusters
sample | infected
clusters | infected
clusters | infected | sample


#Save graphs

ggsave("projects/YOUR_PROJECT_NAME/results/graphs/integrated/sample.png",          sample,          width = 6, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/integrated/clusters.png",        clusters,        width = 6, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/integrated/infected.png",        infected,        width = 6, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/integrated/sample_clusters.png", sample | clusters,  width = 12, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/integrated/infected_sample.png", infected | sample,  width = 12, height = 5)
ggsave("projects/YOUR_PROJECT_NAME/results/graphs/integrated/clusters_infected.png", clusters | infected, width = 12, height = 5)


save_checkpoint(integrated.harmony, "integrated_harmony_normData", dir = "projects/YOUR_PROJECT_NAME/data/processed")
