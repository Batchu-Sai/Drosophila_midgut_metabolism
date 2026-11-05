# Drosophila Midgut `.rds` & Dataset Verification Report
**Date:** November 4, 2026  
**Source Host:** `saibatchu@100.84.8.123` (`~/Research/Drosophila_midgut/`)  
**Mirrored on Mac mini (`Sais-Mac-mini.local`):** `/Users/sai_mini/Documents/Research/Drosophila_midgut/` and `/tmp/rds_details.txt`

---

## 1. Summary of Verified Objects (for Tobias Reiff & Zeus Antonello)

All target genes requested across both sets of analyses (`esg`, `klu`, `Dl`, `pros`, `stg`, `Eip75B`, `unc-5` / `unc.5`, `fra`, `NetA`, `NetB`, `Hr3`, `trio`, `elav`, `SoxN`, `Sox102F`, `numb`, `Abl`) were verified across the datasets below:

| File (M1 & Mac mini Mirror) | Object Class | Dimensions | Assays / Coordinates | Cell Populations / Clusters | Role in Figures |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `Publication_Figures/Gradient_Plots/Data/tsne_final_ISC_EB.csv` | CSV (`CreateSeuratObject(t(tsne_final))`) | **651 cells × 16,829 cols** (16,827 genes + `X`, `Y`) | Imputed expression + `X`, `Y` (`tSNE_1`, `tSNE_2`) | **651 ISC/EB cells** (`quiescent_ISC` subset via `Dl < 1 & klu < 1 & esg > 0`) | **Primary input** for `Plots_2025.R` (April 2025 co-expression grids for `unc.5`, `fra`, `NetA`, `NetB`, `Hr3`, `trio`, `elav`, `SoxN`, `Sox102F`, `numb`, `Abl`), `Gradient_Plots/Plots.R`, and `ISC_EB_clusters/Plots.R` |
| `metabolism/selected_impute_ISC_EB_sce.rds` | `SingleCellExperiment` | **16,829 genes × 651 cells** | `tpm`, `exprs` | `integrated_celltype`: **EB (330)**, **ISC (321)** | Imputed ISC/EB `SingleCellExperiment` object corresponding directly to the 651 cells in `tsne_final_ISC_EB.csv` (uses `unc.5` naming) |
| `Publication_Figures/ISC_Trajectory/trajectory-sce_ISC_EB.rds` | `SingleCellExperiment` | **10,141 genes × 480 cells** | `logcounts`; ReducedDims: `PCA`, `UMAP` | `Identity`: **qISC (104)**, **pISC (102)**, **EB (274)** | Slingshot trajectory object (`qISC -> pISC -> EB`) used in `ISC_Trajectory/Plots.R` (Zipper et al., *eLife* 2022 trajectory figures) |
| `Publication_Figures/Metabolic_Genes_Trajectory/Data/integrated.rds` (same as `metabolism/integrated.rds`) | `Seurat` (v3, update with `UpdateSeuratObject()`) | **16,986 genes × 11,804 cells** (26 metadata cols) | Assays: `RNA`, `SCT`, `integrated`; Reductions: `pca`, `umap` | `integrated_celltype`: `ISC/EB` (799), `EE` (849), `aEC` (1,828), `mEC` (449), `pEC` (1,182), `LFC` (549), `iron/copper` (559), `cardia` (225), plus unannotated clusters | Full midgut integrated Seurat object used as the parent dataset for `Metabolic_Genes_Trajectory/Plots.R` and `ISC_Trajectory/Plots.R` |
| `metabolism/Figures/selected_sce.rds` | `SingleCellExperiment` | **16,827 genes × 10,605 cells** | `tpm`, `exprs` | `integrated_celltype`: `ISC/EB` (651), `EE` (914), `aEC` (1,706), `mEC` (456), `pEC` (1,251), `LFC` (471), `iron/copper` (472), `cardia` (267), etc. | Full midgut TPM/expression `SingleCellExperiment` object |
| `zeus/clustering/integrated.rds` | `Seurat` (v3, update with `UpdateSeuratObject()`) | **16,986 genes × 11,805 cells** (26 metadata cols) | Assays: `RNA`, `SCT`, `integrated`; Reductions: `pca`, `umap` | `integrated_celltype`: `ISC/EB` (824), `EE` (805), `aEC` (1,794), `mEC` (467), `pEC` (1,185), `LFC` (552), `iron/copper` (565), `cardia` (225) | Earlier clustering snapshot in `zeus/clustering/` |

---

## 2. Important Technical Notes for Re-Running Plots

1. **Gene Symbol Formatting (`unc-5` vs `unc.5`)**:
   - In `integrated.rds`, `selected_sce.rds`, and `trajectory-sce_ISC_EB.rds`, the gene is stored with a hyphen as **`unc-5`**.
   - In `Plots_2025.R`, `read.csv("./Data/tsne_final_ISC_EB.csv", row.names = 1)` uses R's default `check.names = TRUE`, which converts `unc-5` to **`unc.5`**.
2. **Seurat Version Compatibility**:
   - `integrated.rds` was saved under Seurat v3. When loading in Seurat v4/v5, run `obj <- UpdateSeuratObject(readRDS("integrated.rds"))` before calling `DimPlot` or `FeaturePlot`.
3. **Quiescent ISC (`qISC`) Definition**:
   - In `Gradient_Plots/Plots.R`, putative quiescent ISCs are subsetted from the 651 ISC/EB cells via:
     ```r
     quiescent_ISC <- subset(obj_seur, subset = Dl < 1 & klu < 1 & esg > 0)
     ```
   - In `ISC_Trajectory/Plots.R`, explicit `qISC` (104 cells), `pISC` (102 cells), and `EB` (274 cells) labels come from `qISC_ISC_EB_metadata.csv` and are stored in `colData(sce)$Identity` inside `trajectory-sce_ISC_EB.rds`.
