#!/usr/bin/env Rscript
# ==============================================================================
# Script:       Scripts/Publication_Figures/Metabolic_Genes_Trajectory/Plots.R
# Reference:    Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022)
#               eLife 11:e78093. https://doi.org/10.7554/eLife.78093
# Description:
#   Plots the 4 whole-midgut Slingshot differentiation lineages across 5,913
#   cells (ISC/EB -> dEC -> aEC, pEC, mEC, EE) and smoothed pseudotime
#   expression profiles for key metabolic and progenitor marker genes.
# Usage:
#   Rscript Scripts/Publication_Figures/Metabolic_Genes_Trajectory/Plots.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(viridis)
  library(cowplot)
  library(ggpubr)
})

set.seed(10)
options(timeout = max(600, getOption("timeout")))

data_root <- if (dir.exists("./Data")) {
  "./Data"
} else if (dir.exists("../../../Data")) {
  "../../../Data"
} else {
  "."
}

out_dir <- if (dir.exists("./Scripts")) "./Reproduced_Figures" else "./Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

release_base_url <- "https://github.com/Batchu-Sai/Drosophila_midgut_metabolism/releases/download/v1.0.0"

ensure_file <- function(fname) {
  local_cand <- file.path(data_root, fname)
  if (!file.exists(local_cand)) {
    message(sprintf("Downloading %s from GitHub Release v1.0.0...", fname))
    download.file(paste0(release_base_url, "/", fname), destfile = local_cand, mode = "wb")
  }
  return(local_cand)
}

sce_path <- ensure_file("trajectory-sce_diff.rds")
sl_path  <- ensure_file("slingshot-trajectory_diff.rds")

sce <- readRDS(sce_path)
sl1 <- readRDS(sl_path)

# Relabel Cluster "2" to differentiating enterocytes ("dEC") as in Zipper et al.
sl1@colData$integrated_celltype[sl1@colData$integrated_celltype == "2"] <- "dEC"
sshot      <- sl1@int_metadata$slingshot
raw_curves <- attr(sshot, "curves")

pseudotimes <- data.frame(
  cell        = colnames(sl1),
  pseudotime1 = sl1$slingPseudotime_1,
  pseudotime2 = sl1$slingPseudotime_2,
  pseudotime3 = sl1$slingPseudotime_3,
  pseudotime4 = sl1$slingPseudotime_4,
  PC_1        = reducedDims(sl1)$PCA[, "PC_1"],
  PC_2        = reducedDims(sl1)$PCA[, "PC_2"],
  Celltype    = factor(
    sl1$integrated_celltype,
    levels = c("ISC/EB", "dEC", "aEC", "mEC", "pEC", "EE")
  )
)

curves <- lapply(1:4, function(i) {
  crv <- raw_curves[[i]]
  data.frame(crv$s[crv$ord, ])
})

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
  pivot_longer(
    cols      = c(pseudotime1, pseudotime2, pseudotime3, pseudotime4),
    names_to  = "lineage",
    values_to = "pseudotime"
  ) %>%
  filter(!is.na(pseudotime)) %>%
  mutate(lineagestring = case_when(
    lineage == "pseudotime1" ~ "ISC/EB->dEC->pEC",
    lineage == "pseudotime2" ~ "ISC/EB->dEC->mEC",
    lineage == "pseudotime3" ~ "ISC/EB->dEC->aEC",
    lineage == "pseudotime4" ~ "ISC/EB->dEC->EE"
  ))

cbp2 <- c("#009E73", "#D55E00", "#293352", "#4E84C4", "#56B4E9", "#CC79A7")

build_lineage_panel <- function(curve_idx, lin_str, title_str) {
  ptop <- ggplot(pseudotimes, aes(PC_1, PC_2, color = Celltype)) +
    geom_point(size = 0.5, alpha = 0.9) +
    scale_color_manual(values = cbp2) +
    theme_classic() +
    geom_path(
      data        = curves[[curve_idx]],
      aes(x = PC_1, y = PC_2),
      inherit.aes = FALSE,
      col         = "black",
      linewidth   = 1.2
    ) +
    ggtitle(title_str) +
    coord_equal()

  pbot <- ggplot(
    fattymarkers[fattymarkers$lineagestring == lin_str, ],
    aes(pseudotime, logcount)
  ) +
    geom_smooth(se = FALSE, color = "red", linewidth = 1.2) +
    facet_wrap(~gene, ncol = 1, strip.position = "top", scales = "free_y") +
    theme_pubr() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  cowplot::plot_grid(ptop, pbot, ncol = 1)
}

aEC <- build_lineage_panel(3, "ISC/EB->dEC->aEC", "ISC/EB -> dEC -> aEC")
pEC <- build_lineage_panel(1, "ISC/EB->dEC->pEC", "ISC/EB -> dEC -> pEC")
mEC <- build_lineage_panel(2, "ISC/EB->dEC->mEC", "ISC/EB -> dEC -> mEC")
EE  <- build_lineage_panel(4, "ISC/EB->dEC->EE",  "ISC/EB -> dEC -> EE")

svg(file.path(out_dir, "trajectory_plot.svg"), width = 25, height = 15)
print(cowplot::plot_grid(aEC, pEC, mEC, EE, ncol = 4))
dev.off()
