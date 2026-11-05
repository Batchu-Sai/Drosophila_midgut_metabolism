# Drosophila_midgut_metabolism

Single-cell RNA-seq clustering, metabolic trajectory analysis, and intestinal progenitor dissection for:

**Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022).** *The microRNA miR-277 controls physiology and pathology of the adult Drosophila midgut by regulating the expression of fatty acid beta-oxidation-related genes in intestinal stem cells.* **eLife** 11:e78093. https://doi.org/10.7554/eLife.78093

Also includes the 2025 follow-up single-cell co-expression and guidance receptor analysis scripts (`Scripts/Publication_Figures/Plots_2025.R`).

- **Interactive Web Explorer:** https://drosophila-midgut-metabolism.pages.dev
- **Pre-Annotated & Original `.rds` Objects (Release v1.0.0):** https://github.com/Batchu-Sai/Drosophila_midgut_metabolism/releases/tag/v1.0.0

---

## Clustering Annotations vs. Hung et al. (2020, GSE120537)

Starting from the adult *Drosophila* midgut 10x and inDrop datasets (Hung et al., *PNAS* 2020, `GSE120537`), Seurat integration, `scImpute` dropout imputation (0.4 threshold), `scran` deconvolution normalization (`computeSumFactors` + `calculateTPM`), and `slingshot` trajectory inference were used to define two refinements over the baseline Hung et al. clustering:

1. **Differentiating Enterocytes (`dEC`, Cluster `"2"`)**:
   Cluster `"2"` (806 cells in `integrated.rds`; 653 cells in `selected_impute_sce.rds`) was unannotated in the original dataset. Metabolic gene clustering and Slingshot trajectory reconstruction (`Scripts/Publication_Figures/Metabolic_Genes_Trajectory/Plots.R`) resolved Cluster `"2"` as differentiating enterocytes (`dEC`) bridging progenitors to differentiated lineages along four trajectories:
   - `ISC/EB -> dEC -> pEC`
   - `ISC/EB -> dEC -> mEC`
   - `ISC/EB -> dEC -> aEC`
   - `ISC/EB -> dEC -> EE`

2. **Progenitor Subtypes (`qISC`, `pISC`, `EB`)**:
   Within the 651 imputed `ISC/EB` cells (`selected_impute_ISC_EB_sce.rds`, `tsne_final_ISC_EB.csv`), progenitors were resolved into three distinct states (`Data/qISC_ISC_EB_metadata.csv`, `Data/trajectory-sce_ISC_EB.rds`):
   - `qISC` (quiescent ISCs, $n = 104$): `esg > 0, Dl < 1, klu < 1`
   - `pISC` (proliferating ISCs, $n = 102$): `Dl >= 1, klu < 1`
   - `EB` (enteroblasts, $n = 274$): `Dl < 1, klu >= 1`
   - Transitioning / other `ISC/EB` ($n = 171$; total = 651 cells)

---

## `.rds` Objects (Available in Release v1.0.0)

In the original analysis pipeline, relabeling Cluster `"2"` to `"dEC"` and merging the `qISC`/`pISC`/`EB` annotations occurred at script runtime. For ease of downstream reuse, both **pre-annotated drop-in `.rds` objects** and the **original unmodified intermediate `.rds` objects** are hosted under [Releases (v1.0.0)](https://github.com/Batchu-Sai/Drosophila_midgut_metabolism/releases/tag/v1.0.0):

| File | Class | Exact Dimensions (Cells x Genes) | Description |
| :--- | :--- | :--- | :--- |
| `Zipper2022_Midgut_Seurat_Annotated.rds` | `Seurat` | 10,605 cells x 16,986 genes (`RNA`) / 10,141 (`SCT`, `integrated`) | Full midgut Seurat object (`pca`, `umap`) with `zipper2022_broad_celltype` (`"2"` relabeled to `"dEC"`) and `zipper2022_celltype` (`qISC`, `pISC`, `EB`, `dEC`, `aEC`, `mEC`, `pEC`, `EE`). |
| `Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds` | `Seurat` | 651 cells x 16,827 genes | Imputed `ISC/EB` Seurat object with publication `tsne` reduction, `slingshot_Identity` (`qISC` [104], `pISC` [102], `EB` [274]), and `marker_subtype`. Tracked directly in the repo root. |
| `Zipper2022_Midgut_Imputed_SCE_Annotated.rds` | `SingleCellExperiment` | 10,605 cells x 16,827 genes | Imputed (`scImpute`) and normalized (`scran`) SCE object with `zipper2022_celltype` plus `all_genes_tSNE` and `metabolic_tSNE` coordinates. |
| `trajectory-sce_ISC_EB.rds` | `SingleCellExperiment` | 480 cells x 10,141 genes | Dissected progenitor Slingshot trajectory (`qISC` [104] -> `pISC` [102] -> `EB` [274]). Tracked in `./Data/`. |
| `slingshot-trajectory_ISC_EB.rds` | `SingleCellExperiment` | 480 cells x 10,141 genes | Companion Slingshot pseudotime and principal curve object for `qISC -> pISC -> EB`. Tracked in `./Data/`. |
| `trajectory-sce_diff.rds` | `SingleCellExperiment` | 5,913 cells x 10,141 genes | Whole-midgut 4-lineage differentiation trajectory object (`799 ISC/EB, 806 dEC, 1,828 aEC, 449 mEC, 1,182 pEC, 849 EE`). |
| `slingshot-trajectory_diff.rds` | `SingleCellExperiment` | 5,913 cells x 10,141 genes | Companion Slingshot pseudotime (`slingPseudotime_1..4`) and 4 lineage curves for whole-midgut differentiation. |
| `Deconvolution_tpm.rds` | `matrix` | 16,827 genes x 10,605 cells | `scran` deconvolution-normalized TPM matrix used for KEGG pathway activity scoring (`Scripts/Dissection_Heatmaps.R`). |
| `selected_impute_ISC_EB_sce.rds` | `SingleCellExperiment` | 651 cells x 16,827 genes | Original imputed `ISC/EB` SCE object prior to subtype merge. |
| `integrated.rds` | `Seurat` | 10,605 cells x 16,986 genes (`RNA`) / 10,141 (`integrated`) | Original Seurat v3 integrated object prior to `dEC` relabeling. |
| `selected_sce.rds` | `SingleCellExperiment` | 10,605 cells x 16,827 genes | Original whole-midgut SCE object prior to dropout imputation. |
| `tsne_final_ISC_EB.csv` | `CSV` | 651 cells x 16,829 cols | Imputed `ISC/EB` expression matrix and `tSNE_1, tSNE_2` coordinates. |

---

## Step-by-Step Local Replication Guide

To clone this repository and reproduce all publication figures or run the interactive analyses on your own machine:

### 1. Clone the Repository
```bash
git clone https://github.com/Batchu-Sai/Drosophila_midgut_metabolism.git
cd Drosophila_midgut_metabolism
```

### 2. Install R Dependencies (R >= 4.2 recommended)
Open R or RStudio and install the required CRAN and Bioconductor packages:
```r
install.packages(c(
  "Seurat", "tidyverse", "viridis", "cowplot", "ggpubr",
  "shiny", "pheatmap", "reshape2", "scales", "RColorBrewer", "matrixStats"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("SingleCellExperiment", "scater"))
```

### 3. Reproduce All Publication Figures (One-Command Pipeline)
From the repository root, run:
```bash
Rscript reproduce_figures.R
```
- **Automatic Large-File Retrieval:** Files $\le 50\text{ MB}$ (`Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds`, `Data/trajectory-sce_ISC_EB.rds`, `Data/slingshot-trajectory_ISC_EB.rds`, `Data/all_tsne.csv`, `Data/metabolic_tsne.csv`) are included directly in the Git clone. On first run, `reproduce_figures.R` automatically downloads the two $>50\text{ MB}$ whole-midgut trajectory objects (`trajectory-sce_diff.rds` and `slingshot-trajectory_diff.rds`) from GitHub Release `v1.0.0` into `./Data/`.
- **Output:** Writes all 7 SVG figures to `./Reproduced_Figures/`:
  1. `All_cells_tSNE_plot.svg` (All genes vs. Metabolic genes tSNE)
  2. `Cell_type_tSNE_plot.svg` (Cell-type-specific metabolic sub-cluster tSNEs)
  3. `ISC_EB_expression_plot.svg` (12-gene imputed expression tSNE grid across 651 ISC/EB cells)
  4. `Gradient_plot.svg` (Seurat two-color blended `klu`/`Dl` and `esg`/`pros` gradients)
  5. `ISC_trajectory_plot.svg` (480-cell `qISC -> pISC -> EB` Slingshot trajectory & pseudotime curves)
  6. `trajectory_plot.svg` (5,913-cell 4-lineage whole-midgut differentiation trajectories)
  7. `SingleAndCoexpr_Grid.svg` (Single-gene and pairwise co-expression tSNE grid)

### 4. Run Individual Analysis Scripts
Each figure directory under `Scripts/` can also be executed standalone from the repository root:
```bash
Rscript Scripts/Publication_Figures/TSNE_plots/Plots.R
Rscript Scripts/Publication_Figures/ISC_EB_clusters/Plots.R
Rscript Scripts/Publication_Figures/Gradient_Plots/Plots.R
Rscript Scripts/Publication_Figures/ISC_Trajectory/Plots.R
Rscript Scripts/Publication_Figures/Metabolic_Genes_Trajectory/Plots.R
Rscript Scripts/Publication_Figures/Plots_2025.R
Rscript Scripts/Dissection_Heatmaps.R
```

### 5. Launch the Interactive R Shiny Application
Explore any of the 16,827 genes interactively across single-gene tSNEs, pairwise co-expression grids, Slingshot pseudotime trajectories, and subtype mean expression tables:
```r
# Run directly from GitHub without cloning:
shiny::runGitHub("Drosophila_midgut_metabolism", "Batchu-Sai")

# Or launch from your local clone:
shiny::runApp("app.R")
```

### 6. Load Pre-Annotated Seurat Objects Directly in R
```r
library(Seurat)

# 1. Whole midgut (10,605 cells; download from Release v1.0.0 if not in working directory)
midgut <- readRDS("Zipper2022_Midgut_Seurat_Annotated.rds")
table(midgut$zipper2022_broad_celltype)
table(midgut$zipper2022_celltype)
DimPlot(midgut, group.by = "zipper2022_celltype", label = TRUE)

# 2. ISC/EB progenitors (651 cells; included in repository root)
isc_eb <- readRDS("Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds")
FeaturePlot(isc_eb, features = c("esg", "Dl", "klu", "Eip75B", "stg", "whd"), reduction = "tsne")
```

---

## Repository Structure

```text
.
├── README.md                                       # Documentation and replication guide
├── app.R                                           # Interactive R Shiny application
├── reproduce_figures.R                             # Master figure reproduction script
├── Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds       # Annotated 651-cell ISC/EB Seurat object
├── Zipper2022_ISC_EB_All_Genes_Mean_Expression.csv # Mean expression across qISC/pISC/EB (16,827 genes)
├── Zipper2022_ISC_EB_Key_Genes_Mean_Expression.csv # Curated marker & metabolic mean expression table
├── Data/                                           # Coordinates, metadata, KEGG pathways, and trajectory RDS
├── Scripts/                                        # Modular figure and KEGG heatmap scripts
├── Figures/                                # Pre-rendered reference SVGs and JPEGs
├── Reproduced_Figures/                             # Output directory populated by reproduce_figures.R
└── web_explorer/                                   # Static web explorer (deployed to Cloudflare Pages)
```
