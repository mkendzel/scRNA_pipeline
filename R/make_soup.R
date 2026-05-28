#' Decontaminate a Seurat object with SoupX
#'
#' Creates a SoupX `SoupChannel` from raw and filtered matrices,
#' estimates the contamination fraction using non-mitochondrial and
#' non-viral genes, and replaces the RNA counts layer with adjusted values.
#'
#' @param sobj Seurat object containing metadata column `sample_id`
#'   and `soup_group` cluster assignments.
#' @param dir_map Named character vector mapping sample IDs to 10x output directories.
#'
#' @details
#' 1. Reads the raw (TOD) matrix from `outs/raw_feature_bc_matrix/`.  
#' 2. Builds a full SoupX model on all genes.  
#' 3. Builds an estimation model excluding mitochondrial (`^MT-` or `^mt-`)
#'    and viral (`^cds-`) genes to estimate the contamination fraction (`rho`).  
#' 4. Applies the estimated `rho` to the full model and adjusts all genes.  
#' 5. Replaces the Seurat object's RNA `counts` layer with decontaminated data.
#'
#' @return Seurat object with decontaminated RNA counts.
#' @examples
#' \dontrun{
#' data_list <- lapply(data_list, make_soup, dir_map = dir_map)
#' }
#' @export
make_soup <- function(sobj, dir_map){
  sample_id <- as.character(sobj$sample_id[1])
  raw <- Read10X(data.dir = file.path(dir_map[[sample_id]], "outs", "raw_feature_bc_matrix"))
  
  toc_all <- GetAssayData(sobj, assay = "RNA", layer = "counts")
  genes_all <- intersect(rownames(raw), rownames(toc_all))
  
  # clusters
  cl <- setNames(as.character(sobj$soup_group), colnames(sobj))
  
  # 1) Full model (used for final adjustment, keeps all genes)
  sc_full <- SoupChannel(tod = raw[genes_all, , drop = FALSE],
                         toc = toc_all[genes_all, , drop = FALSE])
  sc_full <- setClusters(sc_full, cl)
  
  # 2) Estimation model (exclude MT and viral for rho estimation only)
  mt     <- grep("^(MT-|mt-)", genes_all,  value = TRUE)
  viral  <- grep("^cds-",      genes_all,  value = TRUE)   # adjust pattern if needed
  genes_use <- setdiff(genes_all, c(mt, viral))
  
  sc_est <- SoupChannel(tod = raw[genes_use, , drop = FALSE],
                        toc = toc_all[genes_use, , drop = FALSE])
  sc_est <- setClusters(sc_est, cl)
  sc_est <- autoEstCont(sc_est, doPlot = FALSE)            # estimate rho on filtered genes
  rho <- sc_est$meta$rho
  
  # 3) Apply rho to full model and adjust all genes
  sc_full <- setContaminationFraction(sc_full, rho, forceAccept = TRUE)
  out <- adjustCounts(sc_full, roundToInt = TRUE)
  
  # write back counts
  sobj <- SetAssayData(sobj, assay = "RNA", layer = "counts", new.data = out)
  sobj
}