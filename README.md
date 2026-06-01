# scRNA-seq Pipeline

This is a generalized scRNA-seq pipeline. It handles everything from raw Cell Ranger output through QC, contamination removal, doublet removal, clustering, integration, and GSEA. Each experiment gets its own folder under `projects/` with a copy of the analysis scripts.

---

## Setting up a new project

### 1. Copy the template folder

The `projects/_template` folder has the directory structure every project needs. Copy it and rename it after your experiment.

**Mac / Linux:**
```bash
cp -r projects/_template projects/YOUR_PROJECT_NAME
```

**Windows (PowerShell):**
```powershell
Copy-Item -Recurse projects\_template projects\YOUR_PROJECT_NAME
```

Leave `_template` as-is. It's just there to keep the folder layout for future projects.

### 2. Copy the pipeline scripts

Copy everything from `template_scripts/` into your project's `scripts/` folder. All three folders

**Mac / Linux:**
```bash
cp -r template_scripts/. projects/YOUR_PROJECT_NAME/scripts/
```

**Windows (PowerShell):**
```powershell
Copy-Item -Recurse template_scripts\* projects\YOUR_PROJECT_NAME\scripts\
```

Once copied, do a find-and-replace of `YOUR_PROJECT_NAME` with your actual folder name across all the scripts. Then search for `# ADJUST:` to find every line that needs to be changed for your experiment (sample names, conditions, thresholds, file paths, etc.).

Everything was written to be customizable to any specific scRNA design, but some edits may be necessary

### 3. Restore the R environment

Open `scRNA_pipeline.Rproj` in RStudio or Positron first, then run:

```r
renv::restore()
```

This installs the exact package versions the pipeline was built with.

### 4. Download gene set files

The GSEA scripts need an MSigDB Hallmark GMT file. Download it from the MSigDB website and put it here:

```
pathways/gsea/h.all.vYYYY.N.Hs.symbols.gmt
```

Update `gmt_path` in any script that uses it (marked with `# ADJUST:`).

---

## What's in a project folder

```
projects/YOUR_PROJECT_NAME/
├── data/
│   ├── raw/                      Cell Ranger output or other raw input files
│   └── processed/                RDS checkpoints saved as the pipeline runs
├── diagnostics/
│   ├── SoupX/
│   │   ├── clustertree/          Clustree plots for picking the SoupX resolution
│   │   └── umaps/                UMAPs across resolutions for SoupX QC
│   ├── DoubletFinder/
│   │   ├── clustertree/          Clustree plots for picking the DoubletFinder resolution
│   │   ├── umaps/                UMAPs across resolutions for DoubletFinder QC
│   │   ├── doublets/             UMAPs showing predicted doublets
│   │   └── singlets/             UMAPs after doublets are removed
│   ├── sample_clustering/
│   │   ├── clustertree/          Clustree plots for final per-sample clustering
│   │   └── umaps/                UMAPs across resolutions for final per-sample QC
│   ├── merging/                  UMAPs from the merged pre-integration object
│   └── integration/              UMAPs from the Harmony-integrated object
├── results/
│   ├── df/                       summary tables
│   ├── graphs/
│   │   ├── sample_umaps/         Final per-sample UMAP figures
│   │   ├── merged/               Figures from the merged object
│   │   └── integrated/           Figures from the integrated object
│   └── correlative_model/        Regression outputs and model checkpoints
└── scripts/
    ├── pipeline/                 Run these in order (1 through 5)
    ├── figures/                  Figure scripts, run after the pipeline finishes
    └── analysis/                 Downstream analysis (correlative model, etc.)
```

---

## Running the pipeline

Run scripts in `scripts/pipeline/` in order. Each one saves a checkpoint that the next script picks up.

| Script | What it does | Checkpoint saved |
|--------|-------------|-----------------|
| `1_qc_soupx.R` | MAD-based QC filtering and SoupX ambient RNA correction | `data/processed/allData_MADqc_SoupX` |
| `2_doubletfinder.R` | Doublet detection and removal | `data/processed/singletData` |
| `3_clustering_umap.R` | Per-sample clustering and UMAPs | `data/processed/singletData_clustered` |
| `4_merger_integration.R` | Merge samples, optional Harmony integration | `data/processed/merged_clustered_normData`, `integrated_harmony_normData` |
| `5_gsea.R` | Differential expression and GSEA | `results/gsea_*.rds` |

After the pipeline finishes, run the scripts in `scripts/figures/` for exploratory figures, and `scripts/analysis/` for any downstream modeling.

---

## Repo layout

```
scRNA_pipeline/
├── R/                    Helper functions used across all scripts
├── pathways/
│   ├── gsea/             MSigDB GMT files (not included, download separately)
│   └── kegg/             KEGG pathway files (download separately if needed)
├── projects/
│   ├── _template/        The folder skeleton for new projects
│   └── YOUR_PROJECT_NAME/ Your working project
├── template_scripts/     Master copies of all pipeline scripts
├── renv/                 Package library managed by renv
├── renv.lock             Locked package versions
├── install_packages.R    One-time setup script
└── scRNA_pipeline.Rproj  Open this in RStudio before running anything
```
