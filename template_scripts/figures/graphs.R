# ---- Libraries ----
library(Seurat)
library(tidyverse)
library(fgsea)
library(cowplot)

invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))

# ---- Load pipeline outputs ----
# Objects produced by pipeline scripts 4 and 5.
# ADJUST: update project directory path to match your project
project_dir <- "projects/YOUR_PROJECT_NAME"

# Load merged and integrated data sets
merged.seurat     <- load_checkpoint("merged_clustered_normData",    dir = file.path(project_dir, "data/processed"))
integrated.seurat <- load_checkpoint("integrated_harmony_normData",  dir = file.path(project_dir, "data/processed"))

# ADJUST: point combined at whichever object you want for final plots
combined <- merged.seurat
# combined <- integrated.seurat

# ADJUST: set condition levels to match your sample_id values and desired display order
condition_levels <- c(
  "Condition_1",
  "Condition_2"
  # Add all conditions in desired order
)

merged.seurat$sample_id     <- factor(merged.seurat$sample_id,     levels = condition_levels)
integrated.seurat$sample_id <- factor(integrated.seurat$sample_id, levels = condition_levels)
combined$sample_id          <- factor(combined$sample_id,          levels = condition_levels)

# ADJUST: update to your GMT file path and version string
gmt_path      <- "pathways/gsea/h.all.vYYYY.N.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)

# Load all GSEA result RDS files produced by pipeline script 5
gsea_files        <- list.files(file.path(project_dir, "results"), pattern = "^gsea_.*\\.rds$", full.names = TRUE)
gsea_results_full <- lapply(gsea_files, readRDS)
names(gsea_results_full) <- sub("^gsea_", "", sub("\\.rds$", "", basename(gsea_files)))

# ---- Cluster UMAP: merged vs integrated split ----
# Side-by-side comparison to highlight the effect of integration on sample mixing.
p_merged <- DimPlot(merged.seurat, reduction = "umap", group.by = "sample_id") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Merged") +
  theme_minimal()

p_integrated <- DimPlot(
  integrated.seurat,
  reduction = "umap",
  group.by  = "sample_id",
  split.by  = "sample_id",
  pt.size   = 0.1
) + labs(x = "UMAP 1", y = "UMAP 2", title = "Integrated")

p_merged | p_integrated


# Adjust combined above to point at either merged.seurat or integrated.seurat before running the rest of the code.
# ---- Sample UMAP ----
sample_umap <- DimPlot(combined, reduction = "umap", group.by = "sample_id") +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme_minimal()

sample_umap

clusters | sample_umap

# ---- Infection status UMAP ----
# Remove this section if your experiment does not have infection_status metadata.
infected <- DimPlot(combined, reduction = "umap", group.by = "infection_status") +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme_minimal()

infected

clusters | infected | sample_umap

# ---- Per-sample split UMAP ----
p_split <- DimPlot(
  combined,
  reduction = "umap",
  group.by  = "sample_id",
  split.by  = "sample_id",
  pt.size   = 0.1
) + labs(x = "UMAP 1", y = "UMAP 2")

p_split

# ---- Per-sample highlight UMAP loop ----
# Highlights one sample at a time against a grey background.
# Useful for spotting batch effects or uneven distributions before grouping.
meta <- combined@meta.data
for (smp in levels(meta$sample_id)) {
  cells_in <- rownames(meta)[meta$sample_id == smp]
  p <- DimPlot(
    combined,
    reduction       = "umap",
    cells.highlight = list(cells_in),
    cols            = c("gray80", "red")
  ) + ggtitle(smp) + theme(legend.position = "none")
  print(p)
}

# ---- Group UMAP: collapse infection + sample into display groups ----
# ADJUST: define how your conditions map to display groups.
# Remove the infection_status cases if your experiment has no infection metadata.
combined$group <- dplyr::case_when(
  # ADJUST: add a case for each condition that should show as a single group
  combined$sample_id == "Condition_1"              ~ "Control",
  combined$infection_status == "infected"          ~ "Infected",
  combined$infection_status == "bystander"         ~ "Bystander",
  TRUE                                             ~ "Other"
)

# ADJUST: levels and colors
group_levels <- c("Control", "Bystander", "Infected", "Other")
group_colors <- c(
  "Control"   = "darkgreen",
  "Bystander" = "orange",
  "Infected"  = "red",
  "Other"     = "grey50"
)

combined$group <- factor(combined$group, levels = group_levels)

p_group <- DimPlot(combined, reduction = "umap", group.by = "group") +
  scale_color_manual(values = group_colors) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme_minimal()

clusters | p_group

# ---- Feature expression UMAP ----
# ADJUST: any gene or numeric metadata column in the Seurat object
FeaturePlot(
  combined,
  features  = "MX1",
  reduction = "umap",
  cols      = c("lightgray", "red")
) + labs(x = "UMAP 1", y = "UMAP 2") + theme_minimal()

# ---- OPTIONAL: viral / pathogen feature heatmap ----
# Remove this section if your experiment has no viral or pathogen features.
# viral_load was added to merged.seurat in pipeline script 4 if viral features are present.
if ("viral_load" %in% colnames(combined@meta.data) && !all(is.na(combined$viral_load))) {
  FeaturePlot(
    combined,
    features  = "viral_load",
    reduction = "umap",
    cols      = c("lightgray", "red")
  ) + labs(x = "UMAP 1", y = "UMAP 2", title = "Log-normalized Viral Load") + theme_minimal()
}

# ---- GSEA bubble plot ----
# ADJUST: which GSEA ranking method to use
# Options: "gsea_logfc" | "gsea_fc_p" | "gsea_sign_p"
gsea_slot <- "gsea_logfc"

# ADJUST: baseline control label as it appears in contrast names (from pipeline script 5)
base_ctrl <- "Baseline"

# ADJUST: experimental groups that have bystander/infected sub-contrasts.
# Must match the treatment levels used in pipeline script 5.
include_groups <- c("Condition_2", "Condition_3")

# ADJUST: human-readable labels for the x-axis
pretty_group_labels <- c(
  "Condition_2" = "Treatment 2",
  "Condition_3" = "Treatment 3"
)

# Builds contrast strings matching the naming convention from pipeline script 5
make_contrast_names <- function(groups, base_ctrl) {
  unlist(lapply(groups, function(g) {
    c(paste0(g, "_bystander_vs_", base_ctrl),
      paste0(g, "_infected_vs_",  base_ctrl))
  }))
}

# Parses group name and infection status back out of a contrast string
parse_group_status <- function(comparison, base_ctrl) {
  if (grepl("_(infected|bystander)_vs_", comparison)) {
    status <- sub(paste0("^.*_(infected|bystander)_vs_", base_ctrl, "$"), "\\1", comparison)
    group  <- sub(paste0("_(infected|bystander)_vs_", base_ctrl, "$"), "", comparison)
  } else {
    status <- "none"
    group  <- sub(paste0("_vs_", base_ctrl, "$"), "", comparison)
  }
  list(group = group, status = status)
}

wanted <- make_contrast_names(include_groups, base_ctrl)
pieces <- list()
k      <- 0

for (family in names(gsea_results_full)) {
  fam  <- gsea_results_full[[family]]
  if (!is.list(fam)) next
  hits <- intersect(names(fam), wanted)
  if (!length(hits)) next
  for (nm in hits) {
    g <- fam[[nm]][[gsea_slot]]
    if (is.null(g)) next
    df <- as.data.frame(g, stringsAsFactors = FALSE)
    df$comparison <- nm
    k <- k + 1
    pieces[[k]] <- df
  }
}

df_plot <- do.call(rbind, pieces)
gs <- lapply(df_plot$comparison, parse_group_status, base_ctrl = base_ctrl)
df_plot$group  <- vapply(gs, `[[`, character(1), "group")
df_plot$status <- vapply(gs, `[[`, character(1), "status")

df_plot$pathway      <- sub("^HALLMARK_", "", df_plot$pathway)
df_plot$group_pretty <- unname(pretty_group_labels[df_plot$group])
df_plot$x_group      <- paste(df_plot$group_pretty, df_plot$status)

# Bystander columns left, infected columns right
x_levels <- c(
  paste(unname(pretty_group_labels[include_groups]), "bystander"),
  paste(unname(pretty_group_labels[include_groups]), "infected")
)
df_plot$x_group <- factor(df_plot$x_group, levels = x_levels)

df_plot <- df_plot[!is.na(df_plot$x_group) & !is.na(df_plot$padj) & df_plot$padj <= 0.05, , drop = FALSE]
df_plot$pathway <- with(df_plot, reorder(pathway, NES, mean))

split_line <- length(include_groups) + 0.5

GSEA_bubble <- ggplot(
  df_plot,
  aes(x = x_group, y = pathway, size = -log10(padj), fill = NES)
) +
  geom_point(shape = 21, color = "black") +
  geom_vline(xintercept = split_line, linetype = "dotted", linewidth = 0.6, color = "black") +
  scale_size(name = "-log10(padj)", range = c(1, 10)) +
  ylab("Gene Set") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_text(colour = "black"),
    axis.text.x  = element_text(colour = "black", angle = 45, hjust = 0),
    axis.title.x = element_blank()
  ) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES")

GSEA_bubble

# ---- Hallmark pathway module score UMAPs ----
# Requires plot_hallmark_umap() from R/. Subsets to a specific sample and
# infection status before scoring to avoid diluting the signal.
# ADJUST: pathway name, sample, and infection status
plot_hallmark_umap(
  combined, hallmark_sets,
  hallmark_name = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  sample        = "Condition_2",
  infection     = "bystander",
  method        = "module"
)
