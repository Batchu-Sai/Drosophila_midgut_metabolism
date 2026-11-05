#!/usr/bin/env Rscript
# ==============================================================================
# Script:       Scripts/Dissection_Heatmaps.R
# Reference:    Zipper L, Batchu S, Kaya NH, Antonello ZA, Reiff T (2022)
#               eLife 11:e78093. https://doi.org/10.7554/eLife.78093
# Description:
#   Computes cell-type-specific KEGG metabolic pathway activity scores using
#   deconvolved TPM values and 5,000-permutation cell-label shuffling (Xiao et
#   al., Nat Commun 2019 methodology adapted for the Drosophila midgut).
#   Generates pathway activity heatmaps and violin plots across:
#     1. Dissected progenitors: qISC, pISC (ISC), EB, dEC, aEC, mEC, pEC, EE
#     2. Combined progenitor + enterocyte hierarchy (qISC, pISC, ISC_all, EB,
#        dEC, EC_all, aEC, mEC, pEC, EE)
# Usage:
#   Rscript Scripts/Dissection_Heatmaps.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(stringr)
  library(reshape2)
  library(scales)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(RColorBrewer)
  library(matrixStats)
})

options(timeout = max(600, getOption("timeout")))

# Resolve repository Data/ directory
data_root <- if (dir.exists("./Data")) {
  "./Data"
} else if (dir.exists("../Data")) {
  "../Data"
} else {
  stop("Cannot find Data/ directory. Run from repository root or Scripts/ directory.")
}

out_base <- if (dir.exists("./Scripts")) "./Reproduced_Figures" else "./Plots"
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

release_base_url <- "https://github.com/Batchu-Sai/Drosophila_midgut_metabolism/releases/download/v1.0.0"

ensure_release_rds <- function(filename, local_dir = data_root) {
  target_path <- file.path(local_dir, filename)
  if (file.exists(target_path)) return(target_path)
  root_cand <- file.path(".", filename)
  if (file.exists(root_cand)) return(root_cand)
  message(sprintf("Downloading %s from GitHub Release v1.0.0...", filename))
  download.file(paste0(release_base_url, "/", filename), destfile = target_path, mode = "wb")
  return(target_path)
}

pathway_file <- file.path(data_root, "KEGG_dros_pathways.gmt")
source(file.path(data_root, "utils.R"))

sce_path      <- ensure_release_rds("Zipper2022_Midgut_Imputed_SCE_Annotated.rds")
deconv_path   <- ensure_release_rds("Deconvolution_tpm.rds")
qisc_meta_csv <- file.path(data_root, "qISC_ISC_EB_metadata.csv")

# ------------------------------------------------------------------------------
# Helper function: compute pathway activity matrix + permutation p-values
# ------------------------------------------------------------------------------
compute_kegg_pathway_activity <- function(norm_tpm, all_cell_types, pathways, gene_pathway_number, n_perm = 1000) {
  pathway_names <- names(pathways)
  cell_types    <- unique(all_cell_types)

  mean_expression_shuffle   <- matrix(NA, nrow = length(pathway_names), ncol = length(cell_types), dimnames = list(pathway_names, cell_types))
  mean_expression_noshuffle <- matrix(NA, nrow = length(pathway_names), ncol = length(cell_types), dimnames = list(pathway_names, cell_types))
  pvalues_mat               <- matrix(NA, nrow = length(pathway_names), ncol = length(cell_types), dimnames = list(pathway_names, cell_types))

  set.seed(123)
  times <- seq_len(n_perm)
  shuffle_cell_types_list <- lapply(times, function(x) sample(all_cell_types))
  names(shuffle_cell_types_list) <- times

  for (p in pathway_names) {
    genes      <- pathways[[p]]
    genes_comm <- intersect(genes, rownames(norm_tpm))
    if (length(genes_comm) < 5) next

    pathway_metabolic_tpm <- norm_tpm[genes_comm, , drop = FALSE]
    pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm) > 0, , drop = FALSE]
    if (nrow(pathway_metabolic_tpm) < 3) next

    mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x) by(x, all_cell_types, mean))
    keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType > 0.001)]
    if (length(keep) < 3) next

    pathway_metabolic_tpm <- pathway_metabolic_tpm[keep, , drop = FALSE]
    pathway_metabolic_tpm <- t(apply(pathway_metabolic_tpm, 1, function(x) {
      x[x <= 0] <- min(x[x > 0])
      x
    }))

    pathway_number_weight  <- 1 / gene_pathway_number[keep, ]
    mean_exp_eachCellType  <- apply(pathway_metabolic_tpm, 1, function(x) by(x, all_cell_types, mean))
    ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)

    col_quantile <- apply(ratio_exp_eachCellType, 2, function(x) quantile(x, na.rm = TRUE))
    col_q1       <- col_quantile["25%", ]
    col_q3       <- col_quantile["75%", ]
    col_upper    <- col_q3 * 3
    col_lower    <- col_q1 / 3
    outliers     <- apply(ratio_exp_eachCellType, 1, function(x) any((x > col_upper) | (x < col_lower)))

    if (sum(!outliers) < 3) next

    keep <- names(outliers)[!outliers]
    pathway_metabolic_tpm  <- pathway_metabolic_tpm[keep, , drop = FALSE]
    pathway_number_weight  <- 1 / gene_pathway_number[keep, ]
    mean_exp_eachCellType  <- apply(pathway_metabolic_tpm, 1, function(x) by(x, all_cell_types, mean))
    ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
    weight_values          <- pathway_number_weight / sum(pathway_number_weight)
    mean_exp_pathway       <- apply(ratio_exp_eachCellType, 2, function(x) weighted.mean(x, weight_values))

    mean_expression_shuffle[p, ]   <- mean_exp_pathway[cell_types]
    mean_expression_noshuffle[p, ] <- mean_exp_pathway[cell_types]

    group_mean <- function(idx) {
      sapply(cell_types, function(y) rowMeans(pathway_metabolic_tpm[, shuffle_cell_types_list[[idx]] == y, drop = FALSE]))
    }
    column_weight_mean <- function(idx) {
      apply(ratio_exp_eachCellType_list[[idx]], 2, function(y) weighted.mean(y, weight_values))
    }

    mean_exp_eachCellType_list  <- lapply(times, group_mean)
    ratio_exp_eachCellType_list <- lapply(times, function(idx) mean_exp_eachCellType_list[[idx]] / rowMeans(mean_exp_eachCellType_list[[idx]]))
    mean_exp_pathway_list       <- lapply(times, column_weight_mean)

    shuffle_results <- matrix(unlist(mean_exp_pathway_list), ncol = length(cell_types), byrow = TRUE)
    colnames(shuffle_results) <- cell_types

    for (c_type in cell_types) {
      if (is.na(mean_expression_shuffle[p, c_type])) next
      if (mean_expression_shuffle[p, c_type] > 1) {
        pval <- sum(shuffle_results[, c_type] > mean_expression_shuffle[p, c_type]) / n_perm
      } else if (mean_expression_shuffle[p, c_type] < 1) {
        pval <- sum(shuffle_results[, c_type] < mean_expression_shuffle[p, c_type]) / n_perm
      } else {
        pval <- 1
      }
      if (pval > 0.01) mean_expression_shuffle[p, c_type] <- NA
      pvalues_mat[p, c_type] <- pval
    }
  }

  all_NA <- rowAlls(is.na(mean_expression_shuffle))
  mean_expression_shuffle <- mean_expression_shuffle[!all_NA, , drop = FALSE]

  list(
    shuffle   = mean_expression_shuffle,
    noshuffle = mean_expression_noshuffle,
    pvalues   = pvalues_mat
  )
}

# ------------------------------------------------------------------------------
# Run Dissected Midgut KEGG Pathway Activity Analysis (qISC, ISC, EB, EC, EE)
# ------------------------------------------------------------------------------
selected_impute_sce <- readRDS(sce_path)
colData(selected_impute_sce)$integrated_celltype[colData(selected_impute_sce)$integrated_celltype == "2"] <- "dEC"

qISC_metadata <- read.csv(qisc_meta_csv, row.names = 1)
new_coldata   <- colData(selected_impute_sce)
new_coldata$Identity <- new_coldata$integrated_celltype
new_coldata[rownames(new_coldata) %in% rownames(qISC_metadata)[qISC_metadata$Identity == "qISC"], ]$Identity <- "qISC"
new_coldata[rownames(new_coldata) %in% rownames(qISC_metadata)[qISC_metadata$Identity == "EB"], ]$Identity   <- "EB"
new_coldata[rownames(new_coldata) %in% rownames(qISC_metadata)[qISC_metadata$Identity == "ISC"], ]$Identity  <- "ISC"
colData(selected_impute_sce) <- new_coldata

keep_types <- c("dEC", "EE", "aEC", "mEC", "pEC", "qISC", "ISC", "EB")
selected_impute_sce <- selected_impute_sce[, colData(selected_impute_sce)$Identity %in% keep_types]

pathways       <- gmtPathways(pathway_file)
all_cell_types <- as.vector(selected_impute_sce$Identity)
metab_genes    <- rownames(selected_impute_sce)[rowData(selected_impute_sce)$metabolic]
gene_pathway_number <- num_of_pathways(pathway_file, metab_genes)

norm_tpm <- readRDS(deconv_path)
norm_tpm <- norm_tpm[, colnames(norm_tpm) %in% rownames(colData(selected_impute_sce))]

res <- compute_kegg_pathway_activity(norm_tpm, all_cell_types, pathways, gene_pathway_number, n_perm = 1000)
dat <- res$shuffle

sort_row <- c()
for (i in colnames(dat)) {
  select_row <- which(rowMaxs(dat, na.rm = TRUE) == dat[, i])
  tmp <- rownames(dat)[select_row][order(dat[select_row, i], decreasing = TRUE)]
  sort_row <- c(sort_row, tmp)
}
sort_column <- names(apply(dat[sort_row, ], 2, function(x) order(x)[nrow(dat)]))
dat[is.na(dat)] <- 1

mybreaks <- c(
  seq(0, 0.5, length.out = 33),
  seq(0.51, 1.5, length.out = 33),
  seq(1.51, max(dat), length.out = 34)
)
color <- colorRampPalette(c("#0072B2", "white", "#D55E00"))(100)

svg(file.path(out_base, "heatmap_qISC.svg"), width = 8, height = 14)
pheatmap(
  dat[sort_row, sort_column],
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color        = color,
  fontsize     = 9,
  breaks       = mybreaks,
  labels_row   = make_bold_names(dat[sort_row, sort_column], rownames, rownames(dat[sort_row, sort_column])),
  labels_col   = make_bold_names(dat[sort_row, sort_column], colnames, colnames(dat[sort_row, sort_column]))
)
dev.off()
