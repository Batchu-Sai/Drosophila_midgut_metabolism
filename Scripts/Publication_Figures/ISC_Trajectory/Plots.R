#!/usr/bin/env Rscript
# ==============================================================================
# Script:       Scripts/Publication_Figures/ISC_Trajectory/Plots.R
# Reference:    Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022)
#               eLife 11:e78093. https://doi.org/10.7554/eLife.78093
# Description:
#   Reconstructs and plots the Slingshot pseudotime trajectory across the 480
#   dissected progenitor cells (qISC [104] -> pISC [102] -> EB [274]) along with
#   smoothed pseudotime expression profiles for 12 key stem-cell and fatty-acid
#   beta-oxidation genes.
# Usage:
#   Rscript Scripts/Publication_Figures/ISC_Trajectory/Plots.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(viridis)
  library(cowplot)
  library(ggpubr)
})

set.seed(10)

data_root <- if (dir.exists("./Data")) {
  "./Data"
} else if (dir.exists("../../../Data")) {
  "../../../Data"
} else {
  "."
}

out_dir <- if (dir.exists("./Scripts")) "./Reproduced_Figures" else "./Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load SingleCellExperiment and Slingshot trajectory objects (480 cells x 10,141 genes)
sce_path <- if (file.exists(file.path(data_root, "trajectory-sce_ISC_EB.rds"))) {
  file.path(data_root, "trajectory-sce_ISC_EB.rds")
} else {
  "./trajectory-sce_ISC_EB.rds"
}

sl_path <- if (file.exists(file.path(data_root, "slingshot-trajectory_ISC_EB.rds"))) {
  file.path(data_root, "slingshot-trajectory_ISC_EB.rds")
} else {
  "./slingshot-trajectory_ISC_EB.rds"
}

sce <- readRDS(sce_path)
sl1 <- readRDS(sl_path)

# Extract Slingshot principal curve directly from int_metadata (no slingshot binary dependency required)
sshot <- sl1@int_metadata$slingshot
crv1  <- attr(sshot, "curves")[[1]]

pseudotimes <- data.frame(
  cell        = colnames(sl1),
  pseudotime1 = sl1$slingPseudotime_1,
  PC_1        = reducedDims(sl1)$PCA[, "PC_1"],
  PC_2        = reducedDims(sl1)$PCA[, "PC_2"],
  Identity    = sl1$Identity
)

curve1 <- data.frame(crv1$s[crv1$ord, ])

# 12 stem-cell and lipid/fatty-acid oxidation genes profiled along pseudotime
interest_genes <- c(
  "esg", "Dl", "klu", "CG4389", "CG31075", "CG3902",
  "CG4860", "CG9547", "CG5599", "yip2", "whd", "Mtpalpha"
)

keep_genes <- intersect(interest_genes, rownames(sce))
fattydat <- logcounts(sce)[keep_genes, ] %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "logcount")

fattymarkers <- pseudotimes %>%
  left_join(fattydat, by = "cell") %>%
  mutate(lineagestring = "qISC->pISC->EB")

p1 <- ggplot(pseudotimes, aes(PC_1, PC_2, color = Identity)) +
  geom_point(size = 2.0, alpha = 0.9) +
  scale_color_viridis(option = "viridis", discrete = TRUE) +
  theme_classic() +
  theme(text = element_text(face = "bold", size = 12)) +
  geom_path(
    data        = curve1,
    aes(x = PC_1, y = PC_2),
    inherit.aes = FALSE,
    col         = "black",
    linewidth   = 1.5
  ) +
  ggtitle("qISC -> pISC -> EB") +
  coord_equal()

p2 <- ggplot(fattymarkers, aes(pseudotime1, logcount, color = lineagestring)) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.5) +
  facet_wrap(~gene, ncol = 1, strip.position = "top", scales = "free_y") +
  theme_pubr() +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.spacing = unit(0.2, "lines"),
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.text.x   = element_text(size = 8),
    text          = element_text(face = "bold", size = 10)
  ) +
  ylab("log count") +
  xlab("pseudotime")

svg(file.path(out_dir, "ISC_trajectory_plot.svg"), width = 12, height = 18)
print(cowplot::plot_grid(p1, p2, ncol = 1))
dev.off()
