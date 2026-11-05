# Drosophila Adult Midgut Single-Cell RNA-seq Atlas & Progenitor Dissection
**Reference Publication:**  
Zipper L, Batchu S, Kaya N, Antonello ZA, Reiff T (2022). *The microRNA miR-277 controls physiology and pathology of the adult Drosophila midgut by regulating the expression of fatty acid β-oxidation-related genes in intestinal stem cells.* **eLife** 11:e78093.  
**Follow-up Analyses Included:** April 2025 Netrin / Guidance & Co-Expression Analyses (`Plots_2025.R`).

---

## 1. Overview & Why This Clustering Improves on Hung et al. (Perrimon 2020, GSE120537)
Starting from the raw 10x + inDrop adult *Drosophila* midgut dataset (Hung et al., *PNAS* 2020, GEO `GSE120537`), our workflow (`Methods.docx`) performed Seurat CCA/SCT integration, `scImpute` dropout imputation (0.4 threshold), `scran` deconvolution size-factor normalization (`computeSumFactors` + `calculateTPM`), and `Slingshot` pseudotime trajectory inference to achieve two critical biological resolutions:

1. **Identification of Differentiating Enterocytes (`dEC`, Cluster `"2"`)**:  
   In the original Perrimon clustering (`integrated.rds` / `selected_sce.rds`), Cluster `"2"` (806 cells in the full integrated object; 653 cells in the quality-filtered imputed SCE object) was left as an unannotated numeric cluster `"2"`. Using metabolic gene clustering and Slingshot lineage reconstruction (`Metabolic_Genes_Trajectory/Plots.R`), we identified Cluster `"2"` as **differentiating enterocytes (`dEC`)** bridging progenitors (`ISC/EB`) to mature regionally specialized enterocytes (`aEC`, `mEC`, `pEC`) across 4 trajectories:
   - Lineage 1: `ISC/EB -> dEC -> pEC`
   - Lineage 2: `ISC/EB -> dEC -> mEC`
   - Lineage 3: `ISC/EB -> dEC -> aEC`
   - Lineage 4: `ISC/EB -> EE`

2. **High-Resolution Dissection of Progenitors (`ISC/EB` $\rightarrow$ `qISC`, `pISC`, `EB`)**:  
   Using `scImpute` + `scran` deconvolution normalization on the 651 `ISC/EB` progenitor cells (`selected_impute_ISC_EB_sce.rds`, `tsne_final_ISC_EB.csv`), we resolved the lumped `ISC/EB` cluster into three distinct progenitor states along a single Slingshot differentiation trajectory (`qISC -> pISC -> EB`, `ISC_Trajectory/Plots.R`):
   - **`qISC` (Quiescent ISCs, $n = 104$)**: `esg+ Dl- klu-`
   - **`pISC` (Proliferating / Active ISCs, $n = 102$)**: `Dl+ klu-`
   - **`EB` (Enteroblasts, $n = 274$)**: `Dl- klu+`
   - Plus 171 intermediate/transitioning `ISC/EB` cells (`total = 651` cells).

---

## 2. Quick Start: Pre-Annotated Drop-In `.rds` Objects
In the original 2022 scripts, relabeling `"2" -> "dEC"` and merging the `qISC / pISC / EB` metadata took place at script runtime. For immediate convenience, we provide **three pre-annotated `.rds` objects** at the root of this package with all Zipper et al. (2022) labels and coordinates already embedded:

### A. `Zipper2022_Midgut_Seurat_Annotated.rds` (Seurat v5 Object — 10,605 cells $\times$ 16,827 genes)
- **Assays**: `RNA`, `SCT`, `integrated`
- **Reductions**: `pca`, `umap`
- **Pre-baked Metadata (`@meta.data`)**:
  - `zipper2022_broad_celltype` / `celltype`: Broad midgut cell types with Cluster `"2"` annotated as `"dEC"` (`ISC/EB` [799], `dEC` [806], `aEC` [1,828], `mEC` [449], `pEC` [1,182], `EE` [849], `cardia` [225], `LFC` [549], `iron/copper` [559], plus unannotated clusters).
  - `zipper2022_celltype` (also set as `Idents(obj)`): High-resolution annotation with `qISC` (104), `pISC` (102), and `EB` (274) mapped directly into the whole-midgut atlas.
  - `zipper2022_progenitor_subtype`: `"qISC"`, `"pISC"`, `"EB"`, or `NA`.

```r
library(Seurat)
midgut <- readRDS("Zipper2022_Midgut_Seurat_Annotated.rds")
table(midgut$zipper2022_broad_celltype)
table(midgut$zipper2022_celltype)
DimPlot(midgut, group.by = "zipper2022_celltype", label = TRUE)
```

### B. `Zipper2022_Midgut_Imputed_SCE_Annotated.rds` (`SingleCellExperiment` — 10,605 cells $\times$ 16,827 genes)
- Contains `scImpute` + `scran`-normalized expression (`counts`, `logcounts`) used for all publication heatmaps and tSNE comparisons.
- **Pre-baked `colData(sce)`**:
  - `zipper2022_broad_celltype` (`"2"` relabeled to `"dEC"`)
  - `zipper2022_celltype` (`qISC` [104], `pISC` [102], `EB` [274], `ISC/EB` [171], `dEC` [653], `aEC` [1,706], `mEC` [456], `pEC` [1,251], `EE` [914])
  - `all_genes_tSNE1`, `all_genes_tSNE2` & `metabolic_tSNE1`, `metabolic_tSNE2` (for the 5,631 core midgut cells).

### C. `Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds` (Seurat v5 Object — 651 Progenitor Cells $\times$ 16,827 genes)
- Exact Seurat object (`obj_seur`) used in `Gradient_Plots/Plots.R` (2022) and `Plots_2025.R` (2025) with the publication `tsne` reduction (`tSNE_1`, `tSNE_2`), `slingshot_Identity` (`qISC`, `pISC`, `EB`), and `marker_subtype` gates (`esg` / `Dl` / `klu`) pre-loaded.

```r
library(Seurat)
isc_eb <- readRDS("Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds")
FeaturePlot(isc_eb, features = c("esg", "Dl", "klu", "Eip75B", "stg", "whd"), reduction = "tsne")
DimPlot(isc_eb, group.by = "slingshot_Identity", reduction = "tsne")
```

---

## 3. Interactive Shiny App (`app.R`)
You can launch the included interactive R Shiny application locally in R or RStudio:

```r
shiny::runApp("app.R")
```

**Features inside `app.R`:**
1. **Tab 1 — ISC/EB Imputed tSNE & Co-Expression (651 cells)**: Plot any combination of the 16,827 *Drosophila* genes as single-gene expression grids, pairwise co-expression grids (from `Plots_2025.R` with adjustable quantile threshold), or `qISC / pISC / EB` classifications, and export publication-ready `.svg` figures.
2. **Tab 2 — `qISC -> pISC -> EB` Slingshot Pseudotime Trajectory (480 cells)**: Plot any set of genes along the `qISC -> pISC -> EB` differentiation trajectory (`trajectory-sce_ISC_EB.rds`).
3. **Tab 3 — Whole-Midgut All-Genes vs. Metabolic-Genes tSNE (5,631 cells)**: Interactive comparison of all-genes vs. metabolic-genes tSNE (`ISC/EB`, `dEC`, `aEC`, `mEC`, `pEC`, `EE`).
4. **Tab 4 — Subtype Mean Expression Table & CSV Export**: Instant calculation and CSV download of mean expression across `qISC`, `pISC`, `EB`, and all `ISC/EB` cells.

---

## 4. Directory Contents

- **Root Pre-Annotated Files & Tables**:
  - `Zipper2022_Midgut_Seurat_Annotated.rds` — Full 10,605-cell annotated Seurat object.
  - `Zipper2022_Midgut_Imputed_SCE_Annotated.rds` — Full 10,605-cell `scImpute` annotated `SingleCellExperiment` object.
  - `Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds` — 651-cell `ISC/EB` Seurat object with publication tSNE & subtype labels.
  - `Zipper2022_ISC_EB_All_Genes_Mean_Expression.csv` — Mean expression and % expressing across `qISC` ($n=104$), `pISC` ($n=102$), `EB` ($n=274$), and all `ISC/EB` ($n=651$) for all 16,827 genes.
  - `Zipper2022_ISC_EB_Key_Genes_Mean_Expression.csv` — Curated table of key Notch, `E(spl)`, `mir-277` target, β-oxidation, glycolytic/TCA, and 2025 Netrin/guidance genes (`esg`, `Dl`, `klu`, `stg`, `Eip75B`, `whd`, `yip2`, `Mtpalpha`, `CG4389`, `NetA`, `NetB`, `fra`, `unc-5`, etc.).
  - `app.R` — Standalone interactive R Shiny application.
  - `Methods.docx` — Complete single-cell bioinformatics methods from Zipper et al. (*eLife* 2022).
- **`Data/` (Original Unmodified 2021–2022 `.rds` and `.csv` Files)**:
  - `selected_sce.rds` — Original 10,605-cell `SingleCellExperiment` object before imputation.
  - `selected_impute_ISC_EB_sce.rds` — Original 651-cell `ISC/EB` imputed `SingleCellExperiment` object.
  - `trajectory-sce_ISC_EB.rds`, `slingshot-trajectory_ISC_EB.rds`, `trajectory-PCA_ISC_EB.rds` — `qISC -> pISC -> EB` Slingshot trajectory objects (480 cells).
  - `trajectory-sce_diff.rds`, `slingshot-trajectory_diff.rds`, `trajectory-PCA_diff.rds` — Whole-midgut 4-lineage differentiation Slingshot objects (`ISC/EB -> dEC -> pEC/mEC/aEC` & `EE`).
  - `tsne_final_ISC_EB.csv` — Imputed 651-cell $\times$ 16,827-gene expression + tSNE coordinate matrix.
  - `qISC_ISC_EB_metadata.csv`, `ISC_EB_EC_metadata.csv` — Barcode-to-subtype lookup tables.
  - `all_tsne.csv`, `metabolic_tsne.csv` — 5,631-cell tSNE coordinate tables.
- **`Scripts/` (All Original R Scripts from 2022 *eLife* & 2025 Follow-Up)**:
  - `Publication_Figures/Plots_2025.R` — April 2025 single & pairwise co-expression tSNE grid script (`unc-5`, `trio`, `elav`, `SoxN`, `Sox102F`, `numb`, `fra`, `Abl`, `NetA`, `NetB`, `Hr3`, `Eip75B`).
  - `Publication_Figures/Gradient_Plots/Plots.R` — 2022 *eLife* ISC/EB tSNE gradient plots (`esg`, `Dl`, `klu`, `Eip75B`, `stg`, `pros`).
  - `Publication_Figures/ISC_EB_clusters/Plots.R` — 2022 *eLife* `mir-277` target & β-oxidation gene tSNE plots (`CG4389`, `CG31075`, `CG3902`, `CG4860`, `CG9547`, `CG5599`, `yip2`, `whd`, `Mtpalpha`).
  - `Publication_Figures/ISC_Trajectory/Plots.R` — 2022 *eLife* `qISC -> pISC -> EB` Slingshot pseudotime plots.
  - `Publication_Figures/Metabolic_Genes_Trajectory/Plots.R` — 2022 *eLife* whole-midgut 4-lineage Slingshot trajectory plots.
  - `Publication_Figures/TSNE_plots/Plots.R` — 2022 *eLife* all-genes vs. metabolic-genes tSNE comparison.
  - `Dissection_Heatmaps.R` — 2022 *eLife* KEGG metabolic pathway heatmaps (`qISC`, `pISC`, `EB`, `dEC`, `aEC`, `mEC`, `pEC`).
- **`Figures/` (All Rendered `.svg`, `.jpeg`, and `.pdf` Publication Figures & Pathway Heatmaps)**.
