#!/usr/bin/env Rscript
# ==============================================================================
# Script:       Scripts/Publication_Figures/Plots_2025.R
# Reference:    Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022)
#               eLife 11:e78093 & 2025 Follow-Up Co-Expression Analyses
# Description:
#   Provides reusable functions for plotting single-gene imputed expression and
#   pairwise co-expression grids (log1p(gene_x * gene_y) above a configurable
#   quantile threshold) on the 651 ISC/EB tSNE embedding.
# Usage:
#   Rscript Scripts/Publication_Figures/Plots_2025.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(viridis)
  library(cowplot)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

out_dir <- if (dir.exists("./Scripts")) "./Reproduced_Figures" else "./Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load pre-annotated ISC/EB Seurat object (651 cells x 16,827 genes)
rds_candidates <- c(
  "./Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds",
  "../../Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds"
)
rds_path <- rds_candidates[file.exists(rds_candidates)][1]

if (!is.na(rds_path)) {
  obj_seur <- readRDS(rds_path)
} else {
  csv_candidates <- c("./Data/tsne_final_ISC_EB.csv", "../../Data/tsne_final_ISC_EB.csv")
  tsne_final <- read.csv(csv_candidates[file.exists(csv_candidates)][1], row.names = 1, check.names = FALSE)
  expr_mat <- t(as.matrix(tsne_final[, setdiff(colnames(tsne_final), c("X", "Y"))]))
  obj_seur <- CreateSeuratObject(counts = expr_mat)
  obj_seur <- NormalizeData(obj_seur, verbose = FALSE)
  obj_seur <- ScaleData(obj_seur, verbose = FALSE)
  emb <- as.matrix(data.frame(tSNE_1 = tsne_final$X, tSNE_2 = tsne_final$Y, row.names = rownames(tsne_final)))
  obj_seur[["tsne"]] <- CreateDimReducObject(embeddings = emb, key = "tSNE_", assay = DefaultAssay(obj_seur))
}

# ------------------------------------------------------------------------------
# Function: plot_gene_coexpression_pairs_with_singletons()
#   Plots single-gene expression panels with a shared color scale followed by
#   all pairwise co-expression panels (log1p(gx * gy) where both exceed their
#   respective quantile threshold).
# ------------------------------------------------------------------------------
plot_gene_coexpression_pairs_with_singletons <- function(
  seurat_obj,
  genes,
  reduction          = "tsne",
  threshold_quantile = 0.3,
  pt_size            = 0.5,
  pt_alpha           = 0.9,
  output_file        = NULL,
  ncol               = 3
) {
  DefaultAssay(seurat_obj) <- "RNA"
  coord_cols  <- colnames(Embeddings(seurat_obj, reduction))
  valid_genes <- genes[genes %in% rownames(seurat_obj)]

  if (length(valid_genes) == 0) {
    stop("None of the requested genes were found in the Seurat object.")
  }

  expression_values <- unlist(lapply(valid_genes, function(gene) {
    FetchData(seurat_obj, vars = gene)[[gene]]
  }))
  global_expr_range <- range(expression_values, na.rm = TRUE)

  # 1. Single-gene expression panels
  single_gene_plots <- lapply(valid_genes, function(gene) {
    df <- FetchData(seurat_obj, vars = c(coord_cols, gene))
    ggplot(df, aes(x = .data[[coord_cols[1]]], y = .data[[coord_cols[2]]], color = .data[[gene]])) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      scale_color_viridis_c(
        option = "C",
        name   = "Expression",
        limits = global_expr_range,
        oob    = scales::squish
      ) +
      theme(
        panel.background = element_blank(),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.2),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10, face = "bold"),
        axis.title       = element_text(size = 9),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7)
      ) +
      xlab("tSNE 1") +
      ylab("tSNE 2") +
      labs(title = gene)
  })

  # 2. Pairwise co-expression panels
  gene_pairs <- combn(valid_genes, 2, simplify = FALSE)
  coexpr_plots <- lapply(gene_pairs, function(pair) {
    gene_x <- pair[1]
    gene_y <- pair[2]
    df <- FetchData(seurat_obj, vars = c(coord_cols, gene_x, gene_y))
    threshold_x <- quantile(df[[gene_x]], threshold_quantile, na.rm = TRUE)
    threshold_y <- quantile(df[[gene_y]], threshold_quantile, na.rm = TRUE)
    df$coexpr <- df[[gene_x]] > threshold_x & df[[gene_y]] > threshold_y
    df$coexpr_level <- ifelse(df$coexpr, log1p(df[[gene_x]] * df[[gene_y]]), NA)

    ggplot(df, aes(x = .data[[coord_cols[1]]], y = .data[[coord_cols[2]]])) +
      geom_point(
        data  = df[is.na(df$coexpr_level), ],
        color = "gray60",
        size  = pt_size,
        alpha = 0.2
      ) +
      geom_point(
        data  = df[!is.na(df$coexpr_level), ],
        aes(color = coexpr_level),
        size  = pt_size + 0.3,
        alpha = pt_alpha
      ) +
      scale_color_viridis_c(option = "C", name = "Co-expression") +
      theme(
        panel.background = element_blank(),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.2),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10, face = "bold"),
        axis.title       = element_text(size = 9),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7)
      ) +
      xlab("tSNE 1") +
      ylab("tSNE 2") +
      labs(title = paste0(gene_x, " / ", gene_y))
  })

  all_plots <- c(single_gene_plots, coexpr_plots)
  plot_grid_clean <- cowplot::plot_grid(
    plotlist = lapply(all_plots, function(p) p + theme(legend.position = "none")),
    ncol     = ncol
  )
  legend <- cowplot::get_legend(all_plots[[1]] + theme(legend.position = "right"))
  final_plot <- cowplot::plot_grid(plot_grid_clean, legend, rel_widths = c(1, 0.25))

  if (!is.null(output_file)) {
    ggsave(
      output_file,
      plot   = final_plot,
      width  = 5 * ncol + 2,
      height = 5 * ceiling(length(all_plots) / ncol),
      units  = "in"
    )
  }
  return(final_plot)
}

# ------------------------------------------------------------------------------
# Function: plot_single_gene_expression()
#   Plots single-gene expression panels with a shared plasma ("C") color scale.
# ------------------------------------------------------------------------------
plot_single_gene_expression <- function(
  seurat_obj,
  genes,
  reduction   = "tsne",
  pt_size     = 0.5,
  pt_alpha    = 0.9,
  output_file = NULL,
  ncol        = 3
) {
  DefaultAssay(seurat_obj) <- "RNA"
  coord_cols  <- colnames(Embeddings(seurat_obj, reduction))
  valid_genes <- genes[genes %in% rownames(seurat_obj)]

  if (length(valid_genes) == 0) {
    stop("None of the specified genes are present in the Seurat object.")
  }

  expression_values <- unlist(lapply(valid_genes, function(gene) {
    FetchData(seurat_obj, vars = gene)[[gene]]
  }))
  global_expr_range <- range(expression_values, na.rm = TRUE)

  gene_plots <- lapply(valid_genes, function(gene) {
    df <- FetchData(seurat_obj, vars = c(coord_cols, gene))
    ggplot(df, aes(x = .data[[coord_cols[1]]], y = .data[[coord_cols[2]]], color = .data[[gene]])) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      scale_color_viridis_c(
        option = "C",
        name   = "Expression",
        limits = global_expr_range,
        oob    = scales::squish
      ) +
      theme_minimal(base_size = 9) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
        plot.title   = element_text(size = 10, face = "bold"),
        axis.title   = element_text(size = 9),
        axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7)
      ) +
      xlab(paste0(reduction, " 1")) +
      ylab(paste0(reduction, " 2")) +
      labs(title = gene)
  })

  plot_grid_clean <- cowplot::plot_grid(
    plotlist = lapply(gene_plots, function(p) p + theme(legend.position = "none")),
    ncol     = ncol
  )
  legend <- cowplot::get_legend(gene_plots[[1]] + theme(legend.position = "right"))
  final_plot <- cowplot::plot_grid(plot_grid_clean, legend, rel_widths = c(1, 0.25))

  if (!is.null(output_file)) {
    ggsave(
      output_file,
      plot   = final_plot,
      width  = 5 * ncol + 2,
      height = 5 * ceiling(length(gene_plots) / ncol),
      units  = "in"
    )
  }
  return(final_plot)
}

# Resolve unc-5 / unc.5 naming automatically depending on whether loaded from RDS or CSV
unc5_sym <- if ("unc-5" %in% rownames(obj_seur)) "unc-5" else "unc.5"

# Generate the 2025 single + pairwise co-expression grid
plot_gene_coexpression_pairs_with_singletons(
  seurat_obj         = obj_seur,
  genes              = c(unc5_sym, "trio", "elav", "SoxN", "Sox102F", "numb"),
  reduction          = "tsne",
  threshold_quantile = 0.3,
  output_file        = file.path(out_dir, "SingleAndCoexpr_Grid.svg")
)

# Generate the Netrin/Fra/Unc-5/Hr3 single-gene expression grid
plot_single_gene_expression(
  seurat_obj  = obj_seur,
  genes       = c("fra", unc5_sym, "Abl", "NetB", "NetA", "Hr3"),
  reduction   = "tsne",
  pt_size     = 0.3,
  pt_alpha    = 0.9,
  ncol        = 3,
  output_file = file.path(out_dir, "Netrin_Guidance_SingleGene_Grid.svg")
)
