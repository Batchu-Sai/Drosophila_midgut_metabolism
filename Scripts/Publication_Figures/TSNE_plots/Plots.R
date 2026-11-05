#!/usr/bin/env Rscript
# ==============================================================================
# Script:       Scripts/Publication_Figures/TSNE_plots/Plots.R
# Reference:    Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022)
#               eLife 11:e78093. https://doi.org/10.7554/eLife.78093
# Description:
#   Generates whole-midgut tSNE embeddings comparing clustering on all genes vs.
#   metabolic genes alone, plus cell-type-specific metabolic sub-cluster tSNEs
#   (ISC/EB, dEC, aEC, mEC, pEC, EE).
# Usage:
#   Can be run from the repository root or from inside this script's directory:
#     Rscript Scripts/Publication_Figures/TSNE_plots/Plots.R
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

set.seed(10)

# Resolve paths whether run from repo root or from inside TSNE_plots/
data_root <- if (dir.exists("./Data/Cell_Type_Specific")) {
  "./Data"
} else if (dir.exists("../../../Data/Cell_Type_Specific")) {
  "../../../Data"
} else {
  stop("Could not locate Data/ directory. Run from repository root or script folder.")
}

out_dir <- if (dir.exists("./Scripts")) "./Reproduced_Figures" else "./Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Standard colorblind-safe palette across midgut cell types
cbp2 <- c("#009E73", "#D55E00", "#293352", "#4E84C4", "#56B4E9", "#CC79A7")
cell_levels <- c("ISC/EB", "dEC", "aEC", "mEC", "pEC", "EE")

# ------------------------------------------------------------------------------
# 1. Whole-midgut tSNE: All genes vs. Metabolic genes (5,631 cells)
# ------------------------------------------------------------------------------
all_tsne_path <- if (file.exists(file.path(data_root, "all_tsne.csv"))) {
  file.path(data_root, "all_tsne.csv")
} else {
  file.path(data_root, "All_Cells/all_tsne.csv")
}

metab_tsne_path <- if (file.exists(file.path(data_root, "metabolic_tsne.csv"))) {
  file.path(data_root, "metabolic_tsne.csv")
} else {
  file.path(data_root, "All_Cells/metabolic_tsne.csv")
}

all_tsne <- read.csv(all_tsne_path, row.names = 1)
all_tsne$group <- factor(all_tsne$group, levels = cell_levels, labels = cell_levels)

g1 <- ggplot(all_tsne) +
  geom_point(aes(x, y, colour = group), size = 0.5) +
  labs(x = "tSNE1", y = "tSNE2", title = "All genes") +
  theme_bw() +
  scale_color_manual(values = cbp2) +
  theme(
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    text             = element_text(face = "bold", size = 8)
  )

metab_tsne <- read.csv(metab_tsne_path, row.names = 1)
metab_tsne$group <- factor(metab_tsne$group, levels = cell_levels, labels = cell_levels)

g2 <- ggplot(metab_tsne) +
  geom_point(aes(x, y, colour = group), size = 0.5) +
  labs(x = "tSNE1", y = "tSNE2", title = "Metabolic genes") +
  theme_bw() +
  scale_color_manual(values = cbp2) +
  theme(
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    text             = element_text(face = "bold", size = 8)
  )

svg(file.path(out_dir, "All_cells_tSNE_plot.svg"), width = 8, height = 4)
print(cowplot::plot_grid(g1, g2, ncol = 2))
dev.off()

# ------------------------------------------------------------------------------
# 2. Cell-type-specific metabolic sub-cluster tSNE plots
# ------------------------------------------------------------------------------
types <- c("ISC_EB", "dEC", "aEC", "mEC", "pEC", "EE")
names(cbp2) <- types
plots_list <- list()

for (ii in types) {
  coord_file <- file.path(data_root, "Cell_Type_Specific", paste0(ii, "_tSNE_output.csv"))
  coordinates <- read.csv(coord_file, row.names = 1)
  plots_list[[ii]] <- ggplot(coordinates) +
    geom_point(aes(x, y), color = cbp2[ii], size = 0.5) +
    labs(x = "tSNE1", y = "tSNE2", title = ii) +
    theme_bw() +
    theme(
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      text             = element_text(face = "bold", size = 8)
    )
}

svg(file.path(out_dir, "Cell_type_tSNE_plot.svg"), width = 12, height = 8)
print(cowplot::plot_grid(plotlist = plots_list, ncol = 3))
dev.off()
