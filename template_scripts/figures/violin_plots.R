# ---- Libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)

invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))

# ---- Load pipeline output ----
# ADJUST: update project directory path to match your project
project_dir <- "projects/YOUR_PROJECT_NAME"

combined <- load_checkpoint("merged_clustered_normData", dir = file.path(project_dir, "data/processed"))

# ---- Condition setup ----
# ADJUST: ordered levels matching your sample_id values
condition_levels <- c(
  "Condition_1",
  "Condition_2",
  "Condition_3"
)

combined$sample_id <- factor(combined$sample_id, levels = condition_levels)

# Verify names before building groups
table(combined$sample_id)

# ---- Simple split violin: one plot per gene, split by infection_status ----
# Useful as a quick sanity check before building grouped violins.
# Remove split.by if your experiment has no infection_status metadata.
VlnPlot(
  combined,
  features = c("MX1"),   # ADJUST: gene to preview
  group.by = "sample_id",
  split.by = "infection_status",  # remove if no infection metadata
  pt.size  = 0,
  ncol     = 1
) + RotatedAxis()

# ---- Grouped violin setup ----
# Builds one violin per sample x infection_status combination, arranged as:
# [control block] | [bystander block] | [infected block]
#
# ADJUST: conditions shown as a single group (no infection split)
simple_groups <- c("Condition_1")

# ADJUST: conditions shown split into bystander and infected sub-violins
split_groups <- c("Condition_2", "Condition_3")

# ADJUST: infection status values to include in the split groups
include_states <- c("bystander", "infected")

# The raw group key is "<sample_id> <infection_status>" or just "<sample_id>"
x_levels_vln <- c(
  simple_groups,
  paste(split_groups, "bystander"),
  paste(split_groups, "infected")
)

# ADJUST: display labels matching the order above
x_labels_vln <- c(
  simple_groups,
  paste(split_groups, "(Bystander)"),  # ADJUST: customize if needed
  paste(split_groups, "(Infected)")
)

combined$violin_group <- dplyr::case_when(
  combined$sample_id %in% simple_groups ~ as.character(combined$sample_id),
  combined$sample_id %in% split_groups &
    combined$infection_status %in% include_states ~
    paste(as.character(combined$sample_id), as.character(combined$infection_status)),
  TRUE ~ NA_character_
)

combined$violin_group <- factor(
  combined$violin_group,
  levels = x_levels_vln,
  labels = x_labels_vln
)

combined_vln <- subset(combined, subset = !is.na(violin_group))
Idents(combined_vln) <- "violin_group"

# Vertical divider positions: after the simple_groups block, after the bystander block
split_line_1 <- length(simple_groups) + 0.5
split_line_2 <- length(simple_groups) + length(split_groups) + 0.5

# Centered x positions for block header annotations
header_x <- c(
  mean(seq_len(length(simple_groups))),
  length(simple_groups) + mean(seq_len(length(split_groups))),
  length(simple_groups) + length(split_groups) + mean(seq_len(length(split_groups)))
)

# ---- Output directory ----
# ADJUST: update to your project output path
out_dir <- file.path(project_dir, "results/graphs/violin_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Feature violin loop ----
# ADJUST: genes to plot
features <- c("GENE1", "GENE2", "GENE3")

for (feature in features) {

  if (!feature %in% rownames(combined_vln)) {
    warning(paste("Skipping", feature, "- not found in combined_vln"))
    next
  }

  safe_feature <- gsub("[^A-Za-z0-9_.-]+", "_", feature)

  # Grouped violin: control / bystander / infected blocks with block headers
  p <- VlnPlot(
    combined_vln,
    features = feature,
    group.by = "violin_group",
    pt.size  = 0,
    ncol     = 1
  ) +
    RotatedAxis() +
    theme(legend.position = "none") +
    geom_vline(
      xintercept = c(split_line_1, split_line_2),
      linetype   = "dashed",
      linewidth  = 0.5,
      color      = "grey40"
    ) +
    annotate(
      "text",
      x        = header_x,
      y        = Inf,
      label    = c("CTRL", "Bystander", "Infected"),
      vjust    = 1.3,
      fontface = "bold",
      size     = 4
    ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 14, r = 5, b = 5, l = 5)) +
    # Strip the infection status suffix from x tick labels; the block header carries that info
    scale_x_discrete(labels = function(x) sub("\\s+(bystander|infected)\\s*$", "", x))

  # Full split violin: sample_id on x-axis, split by infection_status
  b <- VlnPlot(
    combined,
    features = feature,
    group.by = "sample_id",
    split.by = "infection_status",
    pt.size  = 0,
    ncol     = 1
  ) +
    RotatedAxis() +
    theme(legend.position = "none")

  ggsave(
    filename = file.path(out_dir, paste0(safe_feature, "_grouped_violin.png")),
    plot     = p,
    width    = 10, height = 5, units = "in",
    dpi      = 300, bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, paste0(safe_feature, "_split_violin.png")),
    plot     = b,
    width    = 10, height = 5, units = "in",
    dpi      = 300, bg = "white"
  )
}
