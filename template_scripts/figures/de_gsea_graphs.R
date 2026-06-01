# ---- Libraries ----
library(Seurat)
library(fgsea)
library(dplyr)
library(tidyr)
library(tibble)
library(knitr)
library(pheatmap)
library(grid)
library(openxlsx)

invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))

# ---- Load pipeline outputs ----
# ADJUST: update project directory path to match your project
project_dir <- "projects/YOUR_PROJECT_NAME"

# Load merged and integrated data sets
merged.seurat     <- load_checkpoint("merged_clustered_normData",   dir = file.path(project_dir, "data/processed"))
integrated.seurat <- load_checkpoint("integrated_harmony_normData", dir = file.path(project_dir, "data/processed"))

# ADJUST: point combined at whichever object you want for the analyses below
combined <- merged.seurat
# combined <- integrated.seurat

# ADJUST: set condition levels to match your sample_id values and desired display order
condition_levels <- c(
  "Condition_1",
  "Condition_2"
)

merged.seurat$sample_id     <- factor(merged.seurat$sample_id,     levels = condition_levels)
integrated.seurat$sample_id <- factor(integrated.seurat$sample_id, levels = condition_levels)
combined$sample_id          <- factor(combined$sample_id,           levels = condition_levels)

# ADJUST: update to your GMT file path and version string
gmt_path      <- "pathways/gsea/h.all.vYYYY.N.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)

# Load all GSEA result RDS files produced by pipeline script 5
gsea_files        <- list.files(file.path(project_dir, "results"), pattern = "^gsea_.*\\.rds$", full.names = TRUE)
gsea_results_full <- lapply(gsea_files, readRDS)
names(gsea_results_full) <- sub("^gsea_", "", sub("\\.rds$", "", basename(gsea_files)))

# ADJUST: significance threshold for DE filtering
sig_threshold <- 0.05

# ---- Top DE table ----
# Pulls top 10 up- and down-regulated significant genes from every sub-comparison
# and writes them to an Excel file for quick review.
top_de_df <- do.call(rbind, lapply(names(gsea_results_full), function(comp_nm) {
  comp_obj <- gsea_results_full[[comp_nm]]
  do.call(rbind, lapply(names(comp_obj), function(sub_nm) {
    de_tbl <- comp_obj[[sub_nm]]$de
    de_tbl <- de_tbl[!is.na(de_tbl$avg_log2FC) & !is.na(de_tbl$p_val_adj), , drop = FALSE]
    de_tbl <- de_tbl[de_tbl$p_val_adj < sig_threshold, , drop = FALSE]
    if (nrow(de_tbl) == 0) return(NULL)

    up_tbl   <- head(de_tbl[order(-de_tbl$avg_log2FC), , drop = FALSE], 10)
    down_tbl <- head(de_tbl[order( de_tbl$avg_log2FC), , drop = FALSE], 10)

    rbind(
      cbind(comparison = comp_nm, subcomp = sub_nm, gene = rownames(up_tbl),
            direction = "up",   as.data.frame(up_tbl)),
      cbind(comparison = comp_nm, subcomp = sub_nm, gene = rownames(down_tbl),
            direction = "down", as.data.frame(down_tbl))
    )
  }))
}))

top_de_df <- as.data.frame(top_de_df)
desired_cols <- c("comparison", "subcomp", "gene", "direction",
                  "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
top_de_df <- top_de_df[, intersect(desired_cols, colnames(top_de_df))]

openxlsx::write.xlsx(
  top_de_df,
  file      = file.path(project_dir, "results/df/top_de_up_down.xlsx"),
  overwrite = TRUE
)

# ---- GSEA bubble plot ----
# ADJUST: which GSEA ranking method to use
# Options: "gsea_logfc" | "gsea_fc_p" | "gsea_sign_p"
gsea_slot <- "gsea_logfc"

# ADJUST: baseline control label as it appears in contrast names (from pipeline script 5)
base_ctrl <- "Baseline"

# ADJUST: experimental groups with bystander/infected sub-contrasts
include_groups <- c("Condition_2", "Condition_3")

# ADJUST: human-readable x-axis labels
pretty_group_labels <- c(
  "Condition_2" = "Treatment 2",
  "Condition_3" = "Treatment 3"
)

make_contrast_names <- function(groups, base_ctrl) {
  unlist(lapply(groups, function(g) {
    c(paste0(g, "_bystander_vs_", base_ctrl),
      paste0(g, "_infected_vs_",  base_ctrl))
  }))
}

parse_group_status <- function(comparison, base_ctrl) {
  status <- sub(paste0("^.*_(infected|bystander)_vs_", base_ctrl, "$"), "\\1", comparison)
  group  <- sub(paste0("_(infected|bystander)_vs_", base_ctrl, "$"), "", comparison)
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
  scale_size(name = "-log10(padj)", range = c(4, 10)) +
  ylab("Gene Set") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_text(colour = "black"),
    axis.text.x  = element_text(colour = "black", angle = 45, hjust = 0),
    axis.title.x = element_blank()
  ) +
  scale_fill_distiller(palette = "Spectral")

GSEA_bubble

# ---- GSEA NES heatmap (all contrasts x all pathways) ----
# NES values are masked to NA for non-significant results, then clustered.
# Gives an overview of enrichment across all comparisons simultaneously.

tests <- c("gsea_logfc", "gsea_fc_p", "gsea_sign_p")

heatmap_dir <- file.path(project_dir, "results/graphs")

for (top_nm in names(gsea_results_full)) {
  top_obj <- gsea_results_full[[top_nm]]
  if (!is.list(top_obj) || is.data.frame(top_obj)) next

  for (test_nm in tests) {

    gsea_long <- do.call(rbind, lapply(names(top_obj), function(comp_name) {
      df <- top_obj[[comp_name]][[test_nm]]
      if (is.null(df)) return(NULL)
      df$comparison <- comp_name
      df
    }))

    if (is.null(gsea_long) || nrow(gsea_long) == 0) {
      message("Skipping ", top_nm, " / ", test_nm, " (no data).")
      next
    }

    gsea_long$NES_masked <- ifelse(gsea_long$padj < 0.05, gsea_long$NES, NA)

    gsea_wide <- reshape(
      gsea_long[, c("pathway", "comparison", "NES_masked")],
      timevar   = "comparison",
      idvar     = "pathway",
      direction = "wide"
    )
    colnames(gsea_wide) <- sub("^NES_masked\\.", "", colnames(gsea_wide))

    mat_plot <- as.matrix(gsea_wide[, -1])
    rownames(mat_plot) <- gsea_wide$pathway

    mat_plot <- mat_plot[
      rowSums(!is.na(mat_plot)) > 0,
      colSums(!is.na(mat_plot)) > 0,
      drop = FALSE
    ]

    if (nrow(mat_plot) == 0 || ncol(mat_plot) == 0) {
      message("Skipping ", top_nm, " / ", test_nm, " (no significant pathways).")
      next
    }

    rownames(mat_plot) <- sub("^HALLMARK_", "", rownames(mat_plot))

    mat_cluster <- mat_plot
    mat_cluster[is.na(mat_cluster)] <- 0
    hc_rows <- hclust(dist(mat_cluster))
    hc_cols <- hclust(dist(t(mat_cluster)))

    pal     <- colorRampPalette(c("blue", "white", "red"))(100)
    max_nes <- max(abs(mat_cluster))
    breaks  <- seq(-max_nes, max_nes, length.out = 101)

    safe_top  <- gsub("[/:*?\"<>|]", "_", top_nm)
    safe_test <- gsub("[/:*?\"<>|]", "_", test_nm)
    out_sub   <- file.path(heatmap_dir, safe_test)
    dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)
    fname <- file.path(out_sub, paste0(safe_top, "_", safe_test, "_heatmap.png"))

    message("Plotting: ", top_nm, " / ", test_nm, " -> ", fname)

    png(fname, width = 1200, height = 1600, res = 150)
    pheatmap(
      mat_plot,
      color         = pal,
      na_col        = "darkgrey",
      cluster_rows  = hc_rows,
      cluster_cols  = hc_cols,
      fontsize_row  = 10,
      fontsize_col  = 10,
      angle_col     = 45,
      breaks        = breaks,
      legend_breaks = c(-max_nes, 0, max_nes),
      legend_labels = round(c(-max_nes, 0, max_nes), 1),
      main          = paste0(top_nm, " (", test_nm, ")")
    )
    grid.rect(x = 0.8, y = 0.1, width = 0.06, height = 0.06,
              gp = gpar(fill = "darkgrey"))
    grid.text("Non-significant", x = 0.9, y = 0.1,
              gp = gpar(fontsize = 10))
    dev.off()
  }
}

# ---- Hallmark pathway leading-edge heatmap ----
# Shows avg_log2FC for genes in the leading edge of a specific pathway across
# selected contrasts. Useful for pinpointing which genes drive enrichment.

# ADJUST: pathway and contrast family to use
hm_name     <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
comp_family <- "comparison_family_name"   # top-level key in gsea_results_full
score_col   <- "avg_log2FC"

# ADJUST: contrast names within that family
keep_contrasts <- c(
  "Condition_2_bystander_vs_Baseline",
  "Condition_2_infected_vs_Baseline",
  "Condition_3_bystander_vs_Baseline",
  "Condition_3_infected_vs_Baseline"
)

hm_genes <- hallmark_sets[[hm_name]]
res_list  <- gsea_results_full[[comp_family]]
res_list  <- res_list[names(res_list) %in% keep_contrasts]

de_long <- do.call(rbind, lapply(names(res_list), function(contrast_name) {
  de <- as.data.frame(res_list[[contrast_name]][["de"]])
  de <- de[de$gene %in% hm_genes, , drop = FALSE]
  data.frame(
    gene     = de$gene,
    col_name = contrast_name,
    score    = de[[score_col]],
    stringsAsFactors = FALSE
  )
}))

mat <- de_long %>%
  dplyr::select(gene, col_name, score) %>%
  tidyr::pivot_wider(names_from = col_name, values_from = score) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

max_abs <- max(abs(mat), na.rm = TRUE)
breaks  <- seq(-max_abs, max_abs, length.out = 101)
cols    <- colorRampPalette(c("blue", "white", "red"))(100)

pheatmap::pheatmap(
  mat,
  color        = cols,
  breaks       = breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main         = sprintf("%s (%d genes)", hm_name, nrow(mat))
)

# ---- Cluster composition table ----
# Checks whether any cluster is dominated by a single sample or condition.
# Run this when you notice a cluster that looks biologically unusual.
cluster_tables <- combined@meta.data %>%
  dplyr::count(seurat_clusters, sample_id) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  dplyr::ungroup() %>%
  split(.$seurat_clusters)

for (cl in names(cluster_tables)) {
  cat("\n\n### Cluster", cl, "\n\n")
  print(knitr::kable(
    cluster_tables[[cl]][, c("sample_id", "n", "percent")],
    digits = 2
  ))
}

# ---- Cluster-level DE ----
Idents(combined) <- "seurat_clusters"

# ADJUST: set ident.1 to the cluster number you want to characterize
de_cluster <- FindMarkers(combined, ident.1 = 0, only.pos = FALSE)
de_sig     <- de_cluster[de_cluster$p_val_adj < sig_threshold, ]

top_up   <- head(de_sig[order(-de_sig$avg_log2FC), ], 10)
top_down <- head(de_sig[order( de_sig$avg_log2FC), ], 10)

list(top_up = top_up, top_down = top_down)

# ---- Gene-level DE explorer ----
# Pulls one gene's DE statistics across every sub-comparison. Useful for
# quickly checking whether a gene of interest is consistently regulated.

# ADJUST: gene symbol to look up
gene <- "GENE1"

gene_de_df <- do.call(rbind, lapply(names(gsea_results_full), function(comp_nm) {
  comp_obj <- gsea_results_full[[comp_nm]]
  if (!is.list(comp_obj) || is.data.frame(comp_obj)) return(NULL)
  do.call(rbind, lapply(names(comp_obj), function(sub_nm) {
    de_tbl <- comp_obj[[sub_nm]]$de
    if (!gene %in% rownames(de_tbl)) return(NULL)
    row <- de_tbl[gene, , drop = FALSE]
    cbind(comparison = comp_nm, subcomp = sub_nm, gene = rownames(row), as.data.frame(row))
  }))
}))

gene_de_df <- as.data.frame(gene_de_df)
gene_de_df <- gene_de_df[, intersect(
  c("comparison", "subcomp", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"),
  colnames(gene_de_df)
)]
gene_de_df$p_val_adj <- as.numeric(gene_de_df$p_val_adj)

knitr::kable(gene_de_df, digits = 3)

# Export with conditional formatting: grey out rows where p_val_adj > sig_threshold
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "DE")
openxlsx::writeData(wb, "DE", gene_de_df)
openxlsx::conditionalFormatting(
  wb, "DE",
  cols  = 4:7,
  rows  = 2:(nrow(gene_de_df) + 1),
  type  = "expression",
  rule  = paste0("$", openxlsx::int2col(match("p_val_adj", names(gene_de_df))), "2>",
                 sig_threshold),
  style = openxlsx::createStyle(bgFill = "#CCCCCC", fontColour = "#999999")
)

openxlsx::saveWorkbook(
  wb,
  file.path(project_dir, "results/df", paste0(gene, "_DE_summary.xlsx")),
  overwrite = TRUE
)
