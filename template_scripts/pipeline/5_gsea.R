# ---- Load Libraries -----
library(Seurat)
library(fgsea)

# ---- Project functions ----
invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))
####
# ---- Load merged data set ----
combined <- load_checkpoint("merged_clustered_normData", dir = "data")

# Check Sample Id names
table(combined$sample_id)

# ADJUST: Set condition levels to match your sample names and desired factor order
condition_levels <- c(
  "Condition_1",
  "Condition_2"
  # Add all condition names in desired display order
)

# Set factor levels
combined$sample_id <- factor(
  combined$sample_id,
  levels = condition_levels
)

# ADJUST: Update path to your MSigDB gene set file
# Download from https://www.gsea-msigdb.org/gsea/msigdb/ and put in pathways/gsea/ directory
gmt_path <- "pathways/gsea/h.all.vYYYY.N.Hs.symbols.gmt"   # ADJUST: update version string
hallmark_sets <- gmtPathways(gmt_path)


# =========================================================
# ---- GSEA: Single treatment vs baseline ----
# =========================================================

# ADJUST: Set your comparison levels
treat_level    <- "Treatment_1"   # ADJUST
baseline_level <- "Baseline"      # ADJUST (e.g., "Mock", "Control", "Vehicle")

# Subset to comparison samples
obj <- subset(
  combined,
  subset = sample_id %in% c(baseline_level, treat_level)
)

obj$gsea_group <- ifelse(obj$sample_id == treat_level, treat_level, baseline_level)
table(obj$gsea_group)

# Run GSEA
gsea_result <- run_de_gsea(
  seurat_obj       = obj,
  group_col        = "gsea_group",
  ident_treat      = treat_level,
  ident_base       = baseline_level,
  pathways         = hallmark_sets,
  exclude_features = c("YOUR_FEATURE_PATTERN")   # ADJUST or remove argument
)

# Save results
saveRDS(
  gsea_result,
  file.path("results", paste0("gsea_", gsub("[^A-Za-z0-9]", "_", treat_level),
                              "_vs_", gsub("[^A-Za-z0-9]", "_", baseline_level), ".rds"))
)


# =========================================================
# ---- GSEA: Multiple treatments vs baseline (infected and bystander) ----
# =========================================================

# ADJUST: Set comparison levels
comparison_levels <- c("Treatment_1", "Treatment_2", "Treatment_3")  # ADJUST
baseline_level    <- "Baseline"                                       # ADJUST

# Build subsetted Seurat objects
obj_list <- lapply(comparison_levels, function(lvl) {
  obj <- subset(
    combined,
    subset = sample_id %in% c(baseline_level, lvl)
  )
  obj$gsea_group <- baseline_level

  obj$gsea_group[
    obj$sample_id == lvl & obj$infection_status == "infected"
  ] <- paste0(lvl, "_infected")

  obj$gsea_group[
    obj$sample_id == lvl & obj$infection_status == "bystander"
  ] <- paste0(lvl, "_bystander")

  obj
})
names(obj_list) <- comparison_levels

sapply(obj_list, ncol)

# Run GSEA for each level
gsea_results <- list()

for (lvl in names(obj_list)) {
  message("Starting ", lvl, " ...")

  obj        <- obj_list[[lvl]]
  treat_inf  <- paste0(lvl, "_infected")
  treat_byst <- paste0(lvl, "_bystander")
  base       <- baseline_level

  message("   Running infected vs ", base)
  res_inf <- run_de_gsea(
    seurat_obj       = obj,
    group_col        = "gsea_group",
    ident_treat      = treat_inf,
    ident_base       = base,
    pathways         = hallmark_sets,
    exclude_features = c("YOUR_FEATURE_PATTERN")   # ADJUST or remove argument
  )

  message("   Running bystander vs ", base)
  res_byst <- run_de_gsea(
    seurat_obj       = obj,
    group_col        = "gsea_group",
    ident_treat      = treat_byst,
    ident_base       = base,
    pathways         = hallmark_sets,
    exclude_features = c("YOUR_FEATURE_PATTERN")   # ADJUST or remove argument
  )

  gsea_results[[paste0(lvl, "_infected_vs_",  gsub(" ", "", baseline_level))]]  <- res_inf
  gsea_results[[paste0(lvl, "_bystander_vs_", gsub(" ", "", baseline_level))]]  <- res_byst

  message("Finished ", lvl)
}

# Save
saveRDS(
  gsea_results,
  file.path("results", paste0("gsea_", gsub("[^A-Za-z0-9]", "_", baseline_level),
                              "_comparisons.rds"))
)


# =========================
# ---- Optional Comps  ----
# =========================

# ---- GSEA: Treatment groups (infected) vs a non-baseline infected group ----
# Example: compare infected cells across treatment levels against a reference treatment
#
# comparison_levels <- c("Treatment_1", "Treatment_2", "Treatment_3")  # ADJUST
# baseline_level    <- "Reference_Treatment"                            # ADJUST
#
# obj_list <- list()
# for (lvl in comparison_levels) {
#   obj <- subset(combined, subset = sample_id %in% c(baseline_level, lvl))
#
#   obj$gsea_group <- NA
#   obj$gsea_group[obj$sample_id == baseline_level & obj$infection_status == "infected"] <-
#     paste0("Baseline_", gsub(" ", "_", baseline_level), "_infected")
#   obj$gsea_group[obj$sample_id == lvl & obj$infection_status == "infected"] <-
#     paste0(lvl, "_infected")
#
#   obj <- subset(obj, subset = !is.na(gsea_group))
#   obj_list[[lvl]] <- obj
# }
#
# gsea_results <- list()
# baseline_group <- paste0("Baseline_", gsub(" ", "_", baseline_level), "_infected")
#
# for (lvl in comparison_levels) {
#   obj       <- obj_list[[lvl]]
#   treat_inf <- paste0(lvl, "_infected")
#
#   res_inf <- run_de_gsea(
#     seurat_obj  = obj,
#     group_col   = "gsea_group",
#     ident_treat = treat_inf,
#     ident_base  = baseline_group,
#     pathways    = hallmark_sets
#   )
#   gsea_results[[paste0(lvl, "_infected_vs_ref_infected")]] <- res_inf
# }
#
# saveRDS(gsea_results, file.path("results", "gsea_vs_ref_infected.rds"))


# ---- GSEA: Within-condition infected vs bystander ----
# Example: compare infected and bystander cells within each treatment
#
# within_levels <- c("Treatment_1", "Treatment_2")  # ADJUST
# obj_list <- vector("list", length(within_levels))
# names(obj_list) <- within_levels
#
# for (lvl in within_levels) {
#   obj <- subset(
#     combined,
#     subset = sample_id == lvl & infection_status %in% c("infected", "bystander")
#   )
#   obj$gsea_group <- ifelse(obj$infection_status == "infected",
#                            paste0(lvl, "_infected"),
#                            paste0(lvl, "_bystander"))
#   obj_list[[lvl]] <- obj
# }
#
# gsea_results <- list()
# for (lvl in within_levels) {
#   obj       <- obj_list[[lvl]]
#   treat_inf <- paste0(lvl, "_infected")
#   base      <- paste0(lvl, "_bystander")
#
#   res_inf <- run_de_gsea(
#     seurat_obj  = obj,
#     group_col   = "gsea_group",
#     ident_treat = treat_inf,
#     ident_base  = base,
#     pathways    = hallmark_sets
#   )
#   gsea_results[[paste0(lvl, "_infected_vs_bystander")]] <- res_inf
# }
#
# saveRDS(gsea_results, file.path("results", "gsea_infected_vs_bystander.rds"))
