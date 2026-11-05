#!/usr/bin/env Rscript
# ==============================================================================
# Script:       Scripts/Publication_Figures/ISC_EB_clusters/Plots.R
# Reference:    Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022)
#               eLife 11:e78093. https://doi.org/10.7554/eLife.78093
# Description:
#   Plots single-cell imputed expression (log2 TPM) across the 651 ISC/EB
#   progenitor cells for 12 key stem-cell and fatty-acid oxidation markers:
#   esg, Dl, klu, Eip75B, CG31075, CG3902, CG4860, CG5599, CG9547, Mtpalpha,
#   whd, and yip2.
# Usage:
#   Rscript Scripts/Publication_Figures/ISC_EB_clusters/Plots.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(viridis)
  library(cowplot)
})

set.seed(10)

out_dir <- if (dir.exists("./Scripts")) "./Reproduced_Figures" else "./Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load either the pre-annotated ISC/EB Seurat object or the raw CSV table
rds_candidates <- c(
  "./Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds",
  "../../../Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds"
)
rds_path <- rds_candidates[file.exists(rds_candidates)][1]

if (!is.na(rds_path)) {
  obj_seur <- readRDS(rds_path)
  emb_df   <- as.data.frame(Embeddings(obj_seur, "tsne"))
  colnames(emb_df) <- c("X", "Y")
  expr_mat <- as.matrix(GetAssayData(obj_seur, layer = "data"))
} else {
  csv_candidates <- c("./Data/tsne_final_ISC_EB.csv", "../../../Data/tsne_final_ISC_EB.csv")
  tsne_final <- read.csv(csv_candidates[file.exists(csv_candidates)][1], row.names = 1, check.names = FALSE)
  emb_df   <- data.frame(X = tsne_final$X, Y = tsne_final$Y, row.names = rownames(tsne_final))
  expr_mat <- t(as.matrix(tsne_final[, setdiff(colnames(tsne_final), c("X", "Y"))]))
}

genes_interest <- c(
  "esg", "Dl", "klu", "Eip75B",
  "CG31075", "CG3902", "CG4860", "CG5599",
  "CG9547", "Mtpalpha", "whd", "yip2"
)

isc_plots <- lapply(genes_interest, function(g) {
  df_g <- cbind(emb_df, val = expr_mat[g, ])
  ggplot(df_g) +
    geom_point(aes(as.numeric(X), as.numeric(Y), colour = val), size = 0.5) +
    labs(x = "tSNE1", y = "tSNE2", title = g, color = "log2(TPM)") +
    theme_bw() +
    scale_color_viridis(option = "inferno") +
    theme(
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      text             = element_text(face = "bold", size = 7)
    )
})

svg(file.path(out_dir, "ISC_EB_expression_plot.svg"), width = 12, height = 8)
print(cowplot::plot_grid(plotlist = isc_plots, ncol = 4))
dev.off()
