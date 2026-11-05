#!/usr/bin/env Rscript
# ==============================================================================
# Script:       Scripts/Publication_Figures/Gradient_Plots/Plots.R
# Reference:    Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022)
#               eLife 11:e78093. https://doi.org/10.7554/eLife.78093
# Description:
#   Generates Seurat two-color blended FeaturePlot gradients on the 651 ISC/EB
#   cells (klu vs. Dl) and on the marker-gated quiescent ISC subset (esg vs.
#   pros; gated by Dl < 1 & klu < 1 & esg > 0).
# Usage:
#   Rscript Scripts/Publication_Figures/Gradient_Plots/Plots.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

out_dir <- if (dir.exists("./Scripts")) "./Reproduced_Figures" else "./Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load pre-annotated ISC/EB Seurat object (or construct from CSV if needed)
rds_candidates <- c(
  "./Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds",
  "../../../Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds"
)
rds_path <- rds_candidates[file.exists(rds_candidates)][1]

if (!is.na(rds_path)) {
  obj_seur <- readRDS(rds_path)
} else {
  csv_candidates <- c("./Data/tsne_final_ISC_EB.csv", "../../../Data/tsne_final_ISC_EB.csv")
  tsne_final <- read.csv(csv_candidates[file.exists(csv_candidates)][1], row.names = 1, check.names = FALSE)
  expr_mat <- t(as.matrix(tsne_final[, setdiff(colnames(tsne_final), c("X", "Y"))]))
  obj_seur <- CreateSeuratObject(counts = expr_mat)
  obj_seur <- SetAssayData(obj_seur, layer = "data", new.data = as(expr_mat, "dgCMatrix"))
  emb <- as.matrix(data.frame(tSNE_1 = tsne_final$X, tSNE_2 = tsne_final$Y, row.names = rownames(tsne_final)))
  obj_seur[["tsne"]] <- CreateDimReducObject(embeddings = emb, key = "tSNE_", assay = DefaultAssay(obj_seur))
}

# Marker-based progenitor classification gates:
#   Quiescent ISC (qISC): esg > 0 & Dl < 1 & klu < 1
#   Proliferating ISC:    Dl > 0 & klu < 1
#   Enteroblast (EB):     Dl < 1 & klu > 1
quiescent_ISC <- subset(obj_seur, subset = Dl < 1 & klu < 1 & esg > 0)
ISC           <- subset(obj_seur, subset = Dl > 0 & klu < 1)
EB            <- subset(obj_seur, subset = Dl < 1 & klu > 1)

clean_theme <- theme(
  panel.border     = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line        = element_line(colour = "black"),
  axis.text        = element_blank(),
  axis.ticks       = element_blank()
)

# Seurat v5 FeaturePlot(blend = TRUE) accepts 3 colors: c(low, high_gene1, high_gene2)
p1 <- FeaturePlot(
  obj_seur,
  features        = c("klu", "Dl"),
  blend           = TRUE,
  blend.threshold = 0.0001,
  cols            = c("lightgrey", "#cc5500", "#46cb18"),
  reduction       = "tsne",
  pt.size         = 1.25,
  min.cutoff      = 0,
  combine         = FALSE
)
for (k in 1:3) p1[[k]] <- p1[[k]] + clean_theme
p1[[4]] <- p1[[4]] + ggtitle("")

p2 <- FeaturePlot(
  quiescent_ISC,
  features        = c("esg", "pros"),
  blend           = TRUE,
  blend.threshold = 0.0001,
  cols            = c("lightgrey", "red", "blue"),
  reduction       = "tsne",
  pt.size         = 1.25,
  min.cutoff      = 0,
  combine         = FALSE
)
for (k in 1:3) p2[[k]] <- p2[[k]] + xlim(c(-60, 60)) + ylim(c(-20, 20)) + clean_theme
p2[[4]] <- p2[[4]] + ggtitle("")

svg(file.path(out_dir, "Gradient_plot.svg"), width = 12, height = 8)
print(cowplot::plot_grid(plotlist = c(p1, p2), ncol = 4))
dev.off()
