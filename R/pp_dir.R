#' Preprocess and QC-filter 10x samples from dir.
#'
#' Reads a 10x Genomics sample directory, creates a Seurat object,
#' computes quality-control metrics, applies MAD-based outlier removal,
#' and filters cells based on minimum count and feature thresholds.
#'
#' @param sample_dir Character. Path to the sample directory containing `outs/filtered_feature_bc_matrix`.
#' @param sample_id Character. Unique identifier for the sample (stored in metadata).
#' @param min_counts Numeric. Minimum total counts (nCount_RNA) required to retain a cell.
#' @param min_features Numeric. Minimum number of detected features (nFeature_RNA) required to retain a cell.
#'#' @param mad_counts Numeric. MAD cutoff for `log1p_total_counts`.
#' @param mad_features Numeric. MAD cutoff for `log1p_n_genes_by_counts`.
#' @param mad_mito Numeric. MAD cutoff for `percent.mt`.
#'
#' @details
#' Adds log1p-transformed total counts and feature counts, as well as mitochondrial percentage (`percent.mt`),
#' then removes cells identified as outliers using the `mad_outlier()` function with fixed thresholds:
#'
#' @return A filtered Seurat object containing only high-quality cells.
#' @examples
#' \dontrun{
#' sobj <- preprocess_sample("data/sample1", "Mock-1", min_counts = 800, min_features = 500)
#' }
#' @export

pp_dir <- function(sample_dir, sample_id, 
                   min_counts, min_features, 
                   mad_counts, mad_features, mad_mito) {
  path <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
  sobj <- Read10X(data.dir = path)
  sobj <- CreateSeuratObject(counts = sobj, min.cells = 0, min.features = 200)
  sobj$sample_id <- sample_id
  
  # auto-detect organism by feature prefix
  feat_names <- rownames(sobj[["RNA"]])
  mt_pattern <- if (any(grepl("^mt-", feat_names))) "^mt-" else "^MT-"
  
  sobj$log1p_total_counts      <- log1p(sobj$nCount_RNA)
  sobj$log1p_n_genes_by_counts <- log1p(sobj$nFeature_RNA)
  sobj[["percent.mt"]]         <- PercentageFeatureSet(sobj, pattern = mt_pattern)
  
  n_before <- ncol(sobj)
  qc <- sobj@meta.data
  
  keep <- !mad_outlier(sobj, "log1p_total_counts", mad_counts) &
    !mad_outlier(sobj, "log1p_n_genes_by_counts", mad_features) &
    !mad_outlier(sobj, "percent.mt", mad_mito) &
    qc$nCount_RNA   > min_counts &
    qc$nFeature_RNA > min_features
  
  sobj <- subset(sobj, cells = which(keep))
  cat(sample_id, ": kept", ncol(sobj), "of", n_before,
      "(", round(ncol(sobj) / n_before * 100, 1), "%)\n")
  sobj
}