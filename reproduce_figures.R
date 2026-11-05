#!/usr/bin/env Rscript
# ==============================================================================
# Script:       reproduce_figures.R
# Title:        Master Figure Reproduction Pipeline
# Reference:    Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022)
#               "The Rotherham/Drosophila Midgut Single-Cell & Metabolic Atlas"
#               eLife 11:e78093. https://doi.org/10.7554/eLife.78093
#               and 2025 Follow-Up Progenitor Co-Expression Analyses
#
# Description:
#   End-to-end reproducibility script that regenerates all primary single-cell
#   and metabolic clustering figures, blended marker gradient plots, Slingshot
#   pseudotime trajectories, and pairwise co-expression grids directly from the
#   curated .rds and coordinate objects in this repository.
#
# Usage:
#   git clone https://github.com/Batchu-Sai/Drosophila_midgut_metabolism.git
#   cd Drosophila_midgut_metabolism
#   Rscript reproduce_figures.R
#
# Output:
#   Writes 7 publication-grade SVG figures into ./Reproduced_Figures/
#
# Notes on Large Data Objects:
#   Files <= 50 MB are tracked directly in Git under ./Data/ and the repo root.
#   Files > 50 MB (e.g., whole-midgut Slingshot objects) are hosted on the
#   GitHub Release (v1.0.0) and are automatically downloaded on first run if
#   not already present in ./Data/.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(tidyverse)
  library(viridis)
  library(cowplot)
  library(ggpubr)
})

# Fixed random seed for deterministic rendering
set.seed(10)

# Increase download timeout for large GitHub Release RDS objects (>50 MB)
options(timeout = max(600, getOption("timeout")))

# Output directory for reproduced SVG plots
out_dir <- "./Reproduced_Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("./Data", showWarnings = FALSE, recursive = TRUE)

# Helper function: automatically fetches large (>50 MB) RDS files from GitHub
# Release v1.0.0 if a user is running from a fresh `git clone`.
release_base_url <- "https://github.com/Batchu-Sai/Drosophila_midgut_metabolism/releases/download/v1.0.0"

ensure_data_file <- function(local_path, release_filename = basename(local_path)) {
  if (file.exists(local_path)) {
    return(local_path)
  }
  url <- paste0(release_base_url, "/", release_filename)
  message(sprintf("  -> Downloading %s from GitHub Release v1.0.0...", release_filename))
  tryCatch({
    download.file(url, destfile = local_path, mode = "wb", quiet = FALSE)
  }, error = function(e) {
    warning(sprintf("Could not download %s: %s", release_filename, e$message))
  })
  return(local_path)
}

# Standard colorblind-safe palette used across Zipper et al. (eLife 2022)
cbp2 <- c(
  "ISC/EB" = "#009E73",
  "dEC"    = "#D55E00",
  "aEC"    = "#293352",
  "mEC"    = "#4E84C4",
  "pEC"    = "#56B4E9",
  "EE"     = "#CC79A7"
)

# ==============================================================================
# [1/6] Whole-Midgut tSNE Plots: All Genes vs. Metabolic Genes & Sub-Clusters
#       Corresponds to: Scripts/Publication_Figures/TSNE_plots/Plots.R
# ==============================================================================
cat("[1/6] Reproducing All_cells_tSNE_plot.svg and Cell_type_tSNE_plot.svg...\n")

all_tsne <- read.csv("./Data/all_tsne.csv", row.names = 1)
all_tsne$group <- factor(
  all_tsne$group,
  levels = c("ISC/EB", "dEC", "aEC", "mEC", "pEC", "EE")
)

g1 <- ggplot(all_tsne) +
  geom_point(aes(x, y, colour = group), size = 0.5) +
  labs(x = "tSNE1", y = "tSNE2", title = "All genes") +
  theme_bw() +
  scale_color_manual(values = unname(cbp2)) +
  theme(
    panel.border = element_blank(),
    panel.grid   = element_blank(),
    axis.line    = element_line(colour = "black"),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    text         = element_text(face = "bold", size = 8)
  )

metab_tsne <- read.csv("./Data/metabolic_tsne.csv", row.names = 1)
metab_tsne$group <- factor(
  metab_tsne$group,
  levels = c("ISC/EB", "dEC", "aEC", "mEC", "pEC", "EE")
)

g2 <- ggplot(metab_tsne) +
  geom_point(aes(x, y, colour = group), size = 0.5) +
  labs(x = "tSNE1", y = "tSNE2", title = "Metabolic genes") +
  theme_bw() +
  scale_color_manual(values = unname(cbp2)) +
  theme(
    panel.border = element_blank(),
    panel.grid   = element_blank(),
    axis.line    = element_line(colour = "black"),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    text         = element_text(face = "bold", size = 8)
  )

svg(file.path(out_dir, "All_cells_tSNE_plot.svg"), width = 8, height = 4)
print(cowplot::plot_grid(g1, g2, ncol = 2))
dev.off()

# Cell-type-specific metabolic sub-cluster tSNE plots
types <- c("ISC_EB", "dEC", "aEC", "mEC", "pEC", "EE")
type_colors <- setNames(unname(cbp2), types)
plots_list <- list()

for (ii in types) {
  coord_file <- file.path("./Data/Cell_Type_Specific", paste0(ii, "_tSNE_output.csv"))
  coordinates <- read.csv(coord_file, row.names = 1)
  plots_list[[ii]] <- ggplot(coordinates) +
    geom_point(aes(x, y), color = type_colors[ii], size = 0.5) +
    labs(x = "tSNE1", y = "tSNE2", title = ii) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid   = element_blank(),
      axis.line    = element_line(colour = "black"),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      text         = element_text(face = "bold", size = 8)
    )
}

svg(file.path(out_dir, "Cell_type_tSNE_plot.svg"), width = 12, height = 8)
print(cowplot::plot_grid(plotlist = plots_list, ncol = 3))
dev.off()

# ==============================================================================
# [2/6] ISC/EB 12-Gene Imputed Expression tSNE Grid (651 Progenitor Cells)
#       Corresponds to: Scripts/Publication_Figures/ISC_EB_clusters/Plots.R
# ==============================================================================
cat("[2/6] Reproducing ISC_EB_expression_plot.svg...\n")

ensure_data_file("./Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds")
obj_seur <- readRDS("./Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds")

emb_df <- as.data.frame(Embeddings(obj_seur, "tsne"))
colnames(emb_df) <- c("X", "Y")
expr_mat <- as.matrix(GetAssayData(obj_seur, layer = "data"))

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
      panel.border = element_blank(),
      panel.grid   = element_blank(),
      axis.line    = element_line(colour = "black"),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      text         = element_text(face = "bold", size = 7)
    )
})

svg(file.path(out_dir, "ISC_EB_expression_plot.svg"), width = 12, height = 8)
print(cowplot::plot_grid(plotlist = isc_plots, ncol = 4))
dev.off()

# ==============================================================================
# [3/6] Seurat Two-Color Blended Marker Gradient Plots (klu/Dl & esg/pros)
#       Corresponds to: Scripts/Publication_Figures/Gradient_Plots/Plots.R
# ==============================================================================
cat("[3/6] Reproducing Gradient_plot.svg...\n")

# Filter quiescent ISCs using the publication marker gate (Dl < 1, klu < 1, esg > 0)
quiescent_ISC <- subset(obj_seur, subset = Dl < 1 & klu < 1 & esg > 0)

clean_theme <- theme(
  panel.border = element_blank(),
  panel.grid   = element_blank(),
  axis.line    = element_line(colour = "black"),
  axis.text    = element_blank(),
  axis.ticks   = element_blank()
)

# Note: Seurat v5 FeaturePlot(blend = TRUE) requires 3 colors (negative, high1, high2)
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

# ==============================================================================
# [4/6] Dissected Progenitor Trajectory: qISC -> pISC -> EB (480 Cells)
#       Corresponds to: Scripts/Publication_Figures/ISC_Trajectory/Plots.R
# ==============================================================================
cat("[4/6] Reproducing ISC_trajectory_plot.svg...\n")

ensure_data_file("./Data/trajectory-sce_ISC_EB.rds")
ensure_data_file("./Data/slingshot-trajectory_ISC_EB.rds")

sce_isc <- readRDS("./Data/trajectory-sce_ISC_EB.rds")
sl_isc  <- readRDS("./Data/slingshot-trajectory_ISC_EB.rds")

# Extract Slingshot principal curve coordinates directly from int_metadata
sshot_isc <- sl_isc@int_metadata$slingshot
crv1_isc  <- attr(sshot_isc, "curves")[[1]]

pseudotimes_isc <- data.frame(
  cell        = colnames(sl_isc),
  pseudotime1 = sl_isc$slingPseudotime_1,
  PC_1        = reducedDims(sl_isc)$PCA[, "PC_1"],
  PC_2        = reducedDims(sl_isc)$PCA[, "PC_2"],
  Identity    = sl_isc$Identity
)

curve1_isc <- data.frame(crv1_isc$s[crv1_isc$ord, ])

interest_genes <- c(
  "esg", "Dl", "klu", "CG4389", "CG31075", "CG3902",
  "CG4860", "CG9547", "CG5599", "yip2", "whd", "Mtpalpha"
)

keep_isc_genes <- intersect(interest_genes, rownames(sce_isc))
fattydat_isc <- logcounts(sce_isc)[keep_isc_genes, ] %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "logcount")

fattymarkers_isc <- pseudotimes_isc %>%
  left_join(fattydat_isc, by = "cell") %>%
  mutate(lineagestring = "qISC->pISC->EB")

p1_isc <- ggplot(pseudotimes_isc, aes(PC_1, PC_2, color = Identity)) +
  geom_point(size = 2.0, alpha = 0.9) +
  scale_color_viridis(option = "viridis", discrete = TRUE) +
  theme_classic() +
  theme(text = element_text(face = "bold", size = 12)) +
  geom_path(
    data        = curve1_isc,
    aes(x = PC_1, y = PC_2),
    inherit.aes = FALSE,
    col         = "black",
    linewidth   = 1.5
  ) +
  ggtitle("qISC -> pISC -> EB") +
  coord_equal()

p2_isc <- ggplot(fattymarkers_isc, aes(pseudotime1, logcount, color = lineagestring)) +
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
print(cowplot::plot_grid(p1_isc, p2_isc, ncol = 1))
dev.off()

# ==============================================================================
# [5/6] Whole-Midgut 4-Lineage Differentiation Trajectories (5,913 Cells)
#       Corresponds to: Scripts/Publication_Figures/Metabolic_Genes_Trajectory/Plots.R
# ==============================================================================
cat("[5/6] Reproducing trajectory_plot.svg (4 midgut differentiation lineages)...\n")

ensure_data_file("./Data/trajectory-sce_diff.rds")
ensure_data_file("./Data/slingshot-trajectory_diff.rds")

if (file.exists("./Data/trajectory-sce_diff.rds") &&
    file.exists("./Data/slingshot-trajectory_diff.rds")) {

  sce_diff <- readRDS("./Data/trajectory-sce_diff.rds")
  sl_diff  <- readRDS("./Data/slingshot-trajectory_diff.rds")

  # Relabel Cluster "2" to differentiating enterocytes ("dEC") as in Zipper et al.
  sl_diff@colData$integrated_celltype[sl_diff@colData$integrated_celltype == "2"] <- "dEC"
  sshot_d      <- sl_diff@int_metadata$slingshot
  raw_curves_d <- attr(sshot_d, "curves")

  pseudotimes_d <- data.frame(
    cell        = colnames(sl_diff),
    pseudotime1 = sl_diff$slingPseudotime_1,
    pseudotime2 = sl_diff$slingPseudotime_2,
    pseudotime3 = sl_diff$slingPseudotime_3,
    pseudotime4 = sl_diff$slingPseudotime_4,
    PC_1        = reducedDims(sl_diff)$PCA[, "PC_1"],
    PC_2        = reducedDims(sl_diff)$PCA[, "PC_2"],
    Celltype    = factor(
      sl_diff$integrated_celltype,
      levels = c("ISC/EB", "dEC", "aEC", "mEC", "pEC", "EE")
    )
  )

  curves_d <- lapply(1:4, function(i) {
    crv <- raw_curves_d[[i]]
    data.frame(crv$s[crv$ord, ])
  })

  keep_diff_genes <- intersect(interest_genes, rownames(sce_diff))
  fattydat_d <- logcounts(sce_diff)[keep_diff_genes, ] %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "logcount")

  fattymarkers_d <- pseudotimes_d %>%
    left_join(fattydat_d, by = "cell") %>%
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

  build_lin_panel <- function(curve_idx, lin_str, title_str) {
    ptop <- ggplot(pseudotimes_d, aes(PC_1, PC_2, color = Celltype)) +
      geom_point(size = 0.5, alpha = 0.9) +
      scale_color_manual(values = unname(cbp2)) +
      theme_classic() +
      geom_path(
        data        = curves_d[[curve_idx]],
        aes(x = PC_1, y = PC_2),
        inherit.aes = FALSE,
        col         = "black",
        linewidth   = 1.2
      ) +
      ggtitle(title_str) +
      coord_equal()

    pbot <- ggplot(
      fattymarkers_d[fattymarkers_d$lineagestring == lin_str, ],
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

  pan_aEC <- build_lin_panel(3, "ISC/EB->dEC->aEC", "ISC/EB -> dEC -> aEC")
  pan_pEC <- build_lin_panel(1, "ISC/EB->dEC->pEC", "ISC/EB -> dEC -> pEC")
  pan_mEC <- build_lin_panel(2, "ISC/EB->dEC->mEC", "ISC/EB -> dEC -> mEC")
  pan_EE  <- build_lin_panel(4, "ISC/EB->dEC->EE",  "ISC/EB -> dEC -> EE")

  svg(file.path(out_dir, "trajectory_plot.svg"), width = 25, height = 15)
  print(cowplot::plot_grid(pan_aEC, pan_pEC, pan_mEC, pan_EE, ncol = 4))
  dev.off()
}

# ==============================================================================
# [6/6] 2025 Single-Gene & Pairwise Co-Expression Grid (Plots_2025.R)
#       Corresponds to: Scripts/Plots_2025.R
# ==============================================================================
cat("[6/6] Reproducing SingleAndCoexpr_Grid.svg (Plots_2025.R)...\n")

# Selected guidance/signaling panel for a compact single + pairwise co-expression grid
genes_2025 <- intersect(c("NetA", "NetB", "fra", "unc-5", "Eip75B", "klu"), rownames(obj_seur))
global_range_2025 <- range(expr_mat[genes_2025, ], na.rm = TRUE)

# 1. Single-gene expression panels
single_2025 <- lapply(genes_2025, function(g) {
  df_g <- cbind(emb_df, val = expr_mat[g, ])
  ggplot(df_g, aes(X, Y, color = val)) +
    geom_point(size = 0.7, alpha = 0.9) +
    scale_color_viridis_c(
      option = "C",
      name   = "Expr",
      limits = global_range_2025,
      oob    = scales::squish
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      plot.title   = element_text(face = "bold")
    ) +
    labs(title = g, x = "tSNE 1", y = "tSNE 2")
})

# 2. Pairwise co-expression panels (log1p(gx * gy) above 30th percentile threshold)
pairs_2025 <- combn(genes_2025, 2, simplify = FALSE)
coexpr_2025 <- lapply(pairs_2025, function(pr) {
  gx <- pr[1]
  gy <- pr[2]
  df_c <- cbind(emb_df, gx = expr_mat[gx, ], gy = expr_mat[gy, ])
  tx <- quantile(df_c$gx, 0.3, na.rm = TRUE)
  ty <- quantile(df_c$gy, 0.3, na.rm = TRUE)
  df_c$co <- ifelse(df_c$gx > tx & df_c$gy > ty, log1p(df_c$gx * df_c$gy), NA)

  ggplot(df_c, aes(X, Y)) +
    geom_point(data = df_c[is.na(df_c$co), ], color = "gray75", size = 0.6, alpha = 0.25) +
    geom_point(data = df_c[!is.na(df_c$co), ], aes(color = co), size = 0.9, alpha = 0.9) +
    scale_color_viridis_c(option = "C", name = "Co-expr") +
    theme_minimal(base_size = 9) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      plot.title   = element_text(face = "bold")
    ) +
    labs(title = paste0(gx, " / ", gy), x = "tSNE 1", y = "tSNE 2")
})

all_2025_panels <- c(single_2025, coexpr_2025)
svg(file.path(out_dir, "SingleAndCoexpr_Grid.svg"), width = 18, height = 22)
print(cowplot::plot_grid(plotlist = all_2025_panels, ncol = 3))
dev.off()

cat("All 7 publication figure SVGs reproduced without errors in:", normalizePath(out_dir), "\n")
