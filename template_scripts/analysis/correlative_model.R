# ============================================ #
# ---- scRNA Correlative Model ----
# ============================================ #
# Scores a gene module on cells and correlates it with a continuous
# variable (e.g., viral load) using linear regression.

# ---- Libraries ----
library(Seurat)
library(ggplot2)
library(UCell)
library(fgsea)
library(tidyverse)

# ---- Project functions ----
invisible(lapply(list.files(here::here("R"), full.names = TRUE, pattern = "\\.R$"), source))

# ---- Load Seurat Object ----
# ADJUST: Pick merged or harmony-integrated object
combined <- load_checkpoint("merged_clustered_normData",    dir = "projects/YOUR_PROJECT_NAME/data/processed")
# combined <- load_checkpoint("integrated_harmony_normData", dir = "projects/YOUR_PROJECT_NAME/data/processed")

table(combined$sample_id)

# ADJUST: Set condition levels to match your sample names and desired factor order
condition_levels <- c(
  "Condition_1",
  "Condition_2"
  # Add all condition names in desired display order
)

combined$sample_id <- factor(combined$sample_id, levels = condition_levels)

# ---- Load GSEA Results ----
# ADJUST: Update the filename to match what 5_gsea.R saved
gsea_results <- readRDS("projects/YOUR_PROJECT_NAME/results/gsea_Baseline_comparisons.rds")  # ADJUST

# ADJUST: Set the comparison key used to define your module genes
# Keys follow the pattern produced by 5_gsea.R: "{Treatment}_infected_vs_{Baseline}"
module_comparison <- "Treatment_1_infected_vs_Baseline"  # ADJUST
de <- gsea_results[[module_comparison]][["de"]]

# ---- Load Hallmark Gene Sets ----
# ADJUST: Update path and version string to your MSigDB GMT file
gmt_path <- "pathways/gsea/h.all.vYYYY.N.Hs.symbols.gmt"  # ADJUST
hallmark_sets <- gmtPathways(gmt_path)

# ---- Define Module Genes ----
# ADJUST: Select hallmark sets relevant to your module
module_hallmarks <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",  # ADJUST
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",  # ADJUST
  "HALLMARK_INFLAMMATORY_RESPONSE",      # ADJUST
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"    # ADJUST
)

module_universe <- unique(unlist(hallmark_sets[module_hallmarks]))

# ADJUST: Number of top DE genes to include in the module
n_module_genes <- 50  # ADJUST

module_genes <- de %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0, gene %in% module_universe) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = n_module_genes) %>%
  pull(gene)

message("Module gene count: ", length(module_genes))

# ---- Gene Table ----
# ADJUST: Update category labels to match module_hallmarks above
hallmark_lookup <- bind_rows(
  tibble(gene = hallmark_sets[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], category = "ISG"),           # ADJUST
  tibble(gene = hallmark_sets[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], category = "ISG"),           # ADJUST
  tibble(gene = hallmark_sets[["HALLMARK_INFLAMMATORY_RESPONSE"]],     category = "Inflammatory"),  # ADJUST
  tibble(gene = hallmark_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],  category = "Inflammatory")   # ADJUST
) %>%
  distinct() %>%
  group_by(gene) %>%
  summarise(category = paste(unique(category), collapse = " / "), .groups = "drop")

gene_table <- de %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0, gene %in% module_universe) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = n_module_genes) %>%
  dplyr::select(gene, avg_log2FC, p_val_adj) %>%
  left_join(hallmark_lookup, by = "gene") %>%
  mutate(
    category   = replace_na(category, "Other"),
    avg_log2FC = round(avg_log2FC, 3),
    p_val_adj  = signif(p_val_adj, 3)
  ) %>%
  arrange(category, desc(avg_log2FC))

gene_table

# ---- Score Module ----
# ADJUST: Subset to conditions and infection statuses relevant for scoring
scoring_conditions <- c("Condition_1", "Condition_2")  # ADJUST

score_cells <- subset(
  combined,
  subset = sample_id %in% scoring_conditions &
    infection_status %in% c("bystander", "infected", "unidentified")  # ADJUST or remove if no infection labels
)

# ADJUST: Short descriptive label for the module (no spaces)
module_name <- "ModuleScore"  # ADJUST

score_cells <- AddModuleScore(
  score_cells,
  features = list(module_genes),
  name     = module_name
)

# Transfer module score back to combined (AddModuleScore appends "1" to the name)
module_col <- paste0(module_name, "1")
combined[[module_col]] <- NA
combined[[module_col]][colnames(score_cells)] <- score_cells[[module_col]]

# ---- Build Regression Data Frame ----
# ADJUST: Conditions that have the continuous predictor variable (e.g., viral load)
infected_conditions <- c("Condition_1", "Condition_2")  # ADJUST

regression_df <- combined@meta.data %>%
  filter(
    infection_status == "infected",     # ADJUST if needed
    sample_id %in% infected_conditions,
    !is.na(.data[[module_col]])
  ) %>%
  mutate(
    condition = factor(sample_id, levels = infected_conditions)
  )

# ---- Fit Linear Model ----
# ADJUST: Replace viral_load with your continuous predictor if different; update formula as needed
fit <- lm(as.formula(paste0(module_col, " ~ viral_load * condition")), data = regression_df)  # ADJUST

summary(fit)

# ---- Save Checkpoints ----
save_checkpoint(score_cells,   paste0("scored_cells_", module_name), dir = "projects/YOUR_PROJECT_NAME/data/processed")
save_checkpoint(regression_df, "regression_df",                      dir = "projects/YOUR_PROJECT_NAME/results/correlative_model")
save_checkpoint(fit,           "lm_fit",                             dir = "projects/YOUR_PROJECT_NAME/results/correlative_model")
