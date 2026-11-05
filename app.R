# ==============================================================================
# Interactive Shiny Explorer: Drosophila Midgut Single-Cell & Metabolic Atlas
# Associated Publication: Zipper, Batchu, Kaya, Antonello, Reiff (eLife 2022)
#                         & 2025 Follow-Up Co-Expression Analyses
# ==============================================================================
# Run locally in R / RStudio:
#   shiny::runApp("~/Research/Drosophila_midgut/Zipper2022_Share_Package")
# ==============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(SingleCellExperiment)
  library(ggplot2)
  library(viridis)
  library(cowplot)
  library(dplyr)
  library(tidyr)
})

# Resolve base directory (works whether run inside Zipper2022_Share_Package or parent)
base_dir <- if (file.exists("Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds")) {
  normalizePath(".")
} else if (file.exists("../Publication_Figures/Gradient_Plots/Data/tsne_final_ISC_EB.csv")) {
  normalizePath("..")
} else {
  normalizePath("~/Research/Drosophila_midgut")
}

find_path <- function(rel_candidates) {
  for (cand in rel_candidates) {
    full <- file.path(base_dir, cand)
    if (file.exists(full)) return(full)
  }
  return(NULL)
}

# Load lightweight / primary datasets on startup
isc_eb_rds <- find_path(c(
  "Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds",
  "Zipper2022_Share_Package/Zipper2022_ISC_EB_Imputed_Seurat_tSNE.rds"
))

if (!is.null(isc_eb_rds)) {
  obj_isc_eb <- readRDS(isc_eb_rds)
} else {
  csv_path <- find_path(c(
    "Data/tsne_final_ISC_EB.csv",
    "Publication_Figures/Gradient_Plots/Data/tsne_final_ISC_EB.csv"
  ))
  tsne_final <- read.csv(csv_path, row.names = 1, check.names = FALSE)
  expr_mat <- t(as.matrix(tsne_final[, setdiff(colnames(tsne_final), c("X", "Y"))]))
  obj_isc_eb <- CreateSeuratObject(counts = expr_mat)
  obj_isc_eb <- NormalizeData(obj_isc_eb, verbose = FALSE)
  obj_isc_eb <- ScaleData(obj_isc_eb, verbose = FALSE)
  emb <- as.matrix(data.frame(tSNE_1 = tsne_final$X, tSNE_2 = tsne_final$Y, row.names = rownames(tsne_final)))
  obj_isc_eb[["tsne"]] <- CreateDimReducObject(embeddings = emb, key = "tSNE_", assay = DefaultAssay(obj_isc_eb))
}

isc_traj_path <- find_path(c(
  "Data/trajectory-sce_ISC_EB.rds",
  "Publication_Figures/ISC_Trajectory/trajectory-sce_ISC_EB.rds"
))
isc_slingshot_path <- find_path(c(
  "Data/slingshot-trajectory_ISC_EB.rds",
  "Publication_Figures/ISC_Trajectory/slingshot-trajectory_ISC_EB.rds"
))
all_tsne_path <- find_path(c(
  "Data/all_tsne.csv",
  "Publication_Figures/TSNE_plots/Data/All_Cells/all_tsne.csv"
))
metab_tsne_path <- find_path(c(
  "Data/metabolic_tsne.csv",
  "Publication_Figures/TSNE_plots/Data/All_Cells/metabolic_tsne.csv"
))

sce_isc_traj <- if (!is.null(isc_traj_path)) readRDS(isc_traj_path) else NULL
sl_isc_traj  <- if (!is.null(isc_slingshot_path)) readRDS(isc_slingshot_path) else NULL
all_tsne_df  <- if (!is.null(all_tsne_path)) read.csv(all_tsne_path, row.names = 1) else NULL
met_tsne_df  <- if (!is.null(metab_tsne_path)) read.csv(metab_tsne_path, row.names = 1) else NULL

available_genes <- sort(rownames(obj_isc_eb))
default_2025_genes <- intersect(c("unc-5", "trio", "elav", "SoxN", "Sox102F", "numb", "fra", "Abl", "NetA", "NetB", "Hr3"), available_genes)
default_2022_genes <- intersect(c("esg", "Dl", "klu", "Eip75B", "stg", "pros", "CG4389", "CG31075", "CG3902", "CG4860", "CG9547", "CG5599", "yip2", "whd", "Mtpalpha"), available_genes)

cbp2 <- c(
  "ISC/EB" = "#009E73",
  "dEC"    = "#D55E00",
  "aEC"    = "#293352",
  "mEC"    = "#4E84C4",
  "pEC"    = "#56B4E9",
  "EE"     = "#CC79A7"
)

ui <- navbarPage(
  title = "Drosophila Midgut Atlas (Zipper et al., eLife 2022)",
  tabPanel(
    "1. ISC/EB tSNE & Co-Expression (651 cells)",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("ISC/EB Imputed tSNE Explorer"),
        p("651 progenitor cells (330 EB, 321 ISC; resolved into qISC, pISC, EB) across 16,827 genes."),
        selectizeInput(
          "isc_genes",
          "Select Gene(s):",
          choices = NULL,
          multiple = TRUE
        ),
        radioButtons(
          "isc_mode",
          "Plot Mode:",
          choices = c(
            "Single-Gene Expression Grid" = "single",
            "Single + Pairwise Co-Expression Grid (Plots_2025.R)" = "coexpr",
            "Progenitor Subtype Classification (qISC / pISC / EB)" = "subtype"
          ),
          selected = "single"
        ),
        sliderInput("thresh_q", "Co-expression Quantile Threshold:", min = 0, max = 0.9, value = 0.3, step = 0.05),
        sliderInput("pt_size", "Point Size:", min = 0.3, max = 2.5, value = 0.8, step = 0.1),
        hr(),
        actionButton("preset_2022", "Load 2022 Markers (esg, Dl, klu, Eip75B, stg)", class = "btn-sm btn-info"),
        br(), br(),
        actionButton("preset_2025", "Load 2025 Netrin/Guidance Panel", class = "btn-sm btn-primary"),
        hr(),
        downloadButton("dl_isc_svg", "Download Plot (SVG)")
      ),
      mainPanel(
        width = 9,
        plotOutput("isc_plot", height = "780px")
      )
    )
  ),
  tabPanel(
    "2. qISC -> pISC -> EB Trajectory (Slingshot)",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Slingshot Pseudotime Trajectory"),
        p("480 dissected progenitor cells: qISC (104) -> pISC (102) -> EB (274) from Zipper et al., eLife 2022."),
        selectizeInput(
          "traj_genes",
          "Genes Along Pseudotime:",
          choices = NULL,
          multiple = TRUE
        ),
        downloadButton("dl_traj_svg", "Download Trajectory Figure (SVG)")
      ),
      mainPanel(
        width = 9,
        plotOutput("traj_plot", height = "820px")
      )
    )
  ),
  tabPanel(
    "3. Whole-Midgut Clustering (All vs Metabolic Genes)",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Whole-Midgut Atlas (5,631 / 11,804 cells)"),
        p("Comparison of All-Genes tSNE vs Metabolic-Genes tSNE showing superior separation of ISC/EB, dEC (Cluster 2), aEC, mEC, pEC, and EE.")
      ),
      mainPanel(
        width = 9,
        plotOutput("whole_tsne_plot", height = "520px")
      )
    )
  ),
  tabPanel(
    "4. Mean Expression Table & Export",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Subtype Mean Expression Table"),
        p("Computes mean imputed expression across qISC, pISC, EB, and overall ISC/EB for selected genes."),
        downloadButton("dl_table_csv", "Download Table (CSV)")
      ),
      mainPanel(
        width = 9,
        tableOutput("expr_table")
      )
    )
  )
)

server <- function(input, output, session) {
  updateSelectizeInput(session, "isc_genes", choices = available_genes, selected = c("esg", "Dl", "klu", "Eip75B", "stg", "pros"), server = TRUE)
  if (!is.null(sce_isc_traj)) {
    traj_gene_choices <- sort(rownames(sce_isc_traj))
    updateSelectizeInput(session, "traj_genes", choices = traj_gene_choices, selected = intersect(c("esg", "Dl", "klu", "CG4389", "CG31075", "CG3902", "yip2", "whd", "Mtpalpha"), traj_gene_choices), server = TRUE)
  }

  observeEvent(input$preset_2022, {
    updateSelectizeInput(session, "isc_genes", choices = available_genes, selected = c("esg", "Dl", "klu", "Eip75B", "stg", "pros"), server = TRUE)
  })
  observeEvent(input$preset_2025, {
    updateSelectizeInput(session, "isc_genes", choices = available_genes, selected = default_2025_genes[1:min(6, length(default_2025_genes))], server = TRUE)
  })

  build_isc_plot <- reactive({
    emb <- as.data.frame(Embeddings(obj_isc_eb, "tsne"))
    colnames(emb) <- c("tSNE_1", "tSNE_2")

    if (input$isc_mode == "subtype") {
      emb$Slingshot_Identity <- obj_isc_eb$slingshot_Identity
      emb$Marker_Subtype <- obj_isc_eb$marker_subtype
      p1 <- ggplot(emb, aes(tSNE_1, tSNE_2, color = Slingshot_Identity)) +
        geom_point(size = input$pt_size) +
        theme_bw(base_size = 12) +
        ggtitle("Slingshot Dissected Identity (qISC / pISC / EB)")
      p2 <- ggplot(emb, aes(tSNE_1, tSNE_2, color = Marker_Subtype)) +
        geom_point(size = input$pt_size) +
        theme_bw(base_size = 12) +
        ggtitle("Marker Gate (qISC: esg+ Dl<1 klu<1 | pISC: Dl+ klu<1 | EB: klu+)")
      return(cowplot::plot_grid(p1, p2, ncol = 2))
    }

    genes <- intersect(input$isc_genes, rownames(obj_isc_eb))
    validate(need(length(genes) >= 1, "Please select at least one gene."))

    expr_vals <- unlist(lapply(genes, function(g) FetchData(obj_isc_eb, vars = g)[[g]]))
    global_range <- range(expr_vals, na.rm = TRUE)

    single_plots <- lapply(genes, function(g) {
      df <- cbind(emb, expr = FetchData(obj_isc_eb, vars = g)[[g]])
      ggplot(df, aes(tSNE_1, tSNE_2, color = expr)) +
        geom_point(size = input$pt_size, alpha = 0.9) +
        scale_color_viridis_c(option = "C", name = "Expr", limits = global_range, oob = scales::squish) +
        theme_minimal(base_size = 10) +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
              axis.text = element_blank(), axis.ticks = element_blank(),
              plot.title = element_text(face = "bold")) +
        labs(title = g, x = "tSNE 1", y = "tSNE 2")
    })

    if (input$isc_mode == "coexpr" && length(genes) >= 2) {
      pairs <- combn(genes, 2, simplify = FALSE)
      co_plots <- lapply(pairs, function(pr) {
        gx <- pr[1]; gy <- pr[2]
        df <- cbind(emb, gx = FetchData(obj_isc_eb, vars = gx)[[gx]], gy = FetchData(obj_isc_eb, vars = gy)[[gy]])
        tx <- quantile(df$gx, input$thresh_q, na.rm = TRUE)
        ty <- quantile(df$gy, input$thresh_q, na.rm = TRUE)
        df$co <- ifelse(df$gx > tx & df$gy > ty, log1p(df$gx * df$gy), NA)
        ggplot(df, aes(tSNE_1, tSNE_2)) +
          geom_point(data = df[is.na(df$co), ], color = "gray70", size = input$pt_size, alpha = 0.25) +
          geom_point(data = df[!is.na(df$co), ], aes(color = co), size = input$pt_size + 0.3, alpha = 0.9) +
          scale_color_viridis_c(option = "C", name = "Co-expr") +
          theme_minimal(base_size = 10) +
          theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
                axis.text = element_blank(), axis.ticks = element_blank(),
                plot.title = element_text(face = "bold")) +
          labs(title = paste0(gx, " / ", gy), x = "tSNE 1", y = "tSNE 2")
      })
      all_p <- c(single_plots, co_plots)
    } else {
      all_p <- single_plots
    }

    cowplot::plot_grid(plotlist = all_p, ncol = min(3, length(all_p)))
  })

  output$isc_plot <- renderPlot({ build_isc_plot() })

  output$dl_isc_svg <- downloadHandler(
    filename = function() { "Zipper2022_ISC_EB_tSNE_Plot.svg" },
    content = function(file) {
      ggsave(file, plot = build_isc_plot(), width = 14, height = 10)
    }
  )

  build_traj_plot <- reactive({
    validate(need(!is.null(sce_isc_traj), "trajectory-sce_ISC_EB.rds not found."))
    pca_df <- as.data.frame(reducedDims(sce_isc_traj)$PCA[, 1:2])
    colnames(pca_df) <- c("PC_1", "PC_2")
    pca_df$Identity <- factor(colData(sce_isc_traj)$Identity, levels = c("qISC", "pISC", "EB"))
    pca_df$cell <- rownames(pca_df)

    p1 <- ggplot(pca_df, aes(PC_1, PC_2, color = Identity)) +
      geom_point(size = 2.2, alpha = 0.9) +
      scale_color_viridis_d(option = "viridis") +
      theme_classic(base_size = 12) +
      ggtitle("Slingshot PCA Trajectory: qISC -> pISC -> EB (480 cells)")

    if (!is.null(sl_isc_traj)) {
      pca_df$pseudotime <- sl_isc_traj$slingPseudotime_1
    } else {
      pca_df$pseudotime <- rank(pca_df$PC_1)
    }

    tgenes <- intersect(input$traj_genes, rownames(sce_isc_traj))
    if (length(tgenes) == 0) return(p1)

    mat <- as.matrix(logcounts(sce_isc_traj)[tgenes, , drop = FALSE])
    long_df <- as.data.frame(mat) %>%
      tibble::rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "cell", values_to = "logcount") %>%
      left_join(pca_df[, c("cell", "pseudotime", "Identity")], by = "cell")

    p2 <- ggplot(long_df, aes(pseudotime, logcount)) +
      geom_point(aes(color = Identity), size = 0.8, alpha = 0.4) +
      geom_smooth(se = FALSE, color = "#D55E00", linewidth = 1.4) +
      facet_wrap(~gene, ncol = 2, scales = "free_y") +
      scale_color_viridis_d(option = "viridis") +
      theme_bw(base_size = 11) +
      labs(x = "Slingshot Pseudotime (qISC -> pISC -> EB)", y = "logcount")

    cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.38, 0.62))
  })

  output$traj_plot <- renderPlot({ build_traj_plot() })

  output$dl_traj_svg <- downloadHandler(
    filename = function() { "Zipper2022_ISC_Trajectory_Plot.svg" },
    content = function(file) {
      ggsave(file, plot = build_traj_plot(), width = 12, height = 14)
    }
  )

  output$whole_tsne_plot <- renderPlot({
    validate(need(!is.null(all_tsne_df) && !is.null(met_tsne_df), "all_tsne.csv / metabolic_tsne.csv not found."))
    all_tsne_df$group <- factor(all_tsne_df$group, levels = names(cbp2))
    met_tsne_df$group <- factor(met_tsne_df$group, levels = names(cbp2))

    g1 <- ggplot(all_tsne_df, aes(x, y, color = group)) +
      geom_point(size = 0.6) +
      scale_color_manual(values = cbp2) +
      theme_classic(base_size = 12) +
      labs(title = "All Genes tSNE (5,631 cells)", x = "tSNE 1", y = "tSNE 2")

    g2 <- ggplot(met_tsne_df, aes(x, y, color = group)) +
      geom_point(size = 0.6) +
      scale_color_manual(values = cbp2) +
      theme_classic(base_size = 12) +
      labs(title = "Metabolic Genes tSNE (5,631 cells)", x = "tSNE 1", y = "tSNE 2")

    cowplot::plot_grid(g1, g2, ncol = 2)
  })

  calc_table <- reactive({
    genes <- intersect(input$isc_genes, rownames(obj_isc_eb))
    validate(need(length(genes) >= 1, "Select at least one gene in Tab 1."))
    mat <- as.data.frame(FetchData(obj_isc_eb, vars = genes, layer = "counts"))
    mat$Subtype <- obj_isc_eb$marker_subtype
    by_sub <- mat %>%
      group_by(Subtype) %>%
      summarise(across(all_of(genes), ~ round(mean(.x, na.rm = TRUE), 3)), .groups = "drop")
    overall <- mat %>%
      summarise(across(all_of(genes), ~ round(mean(.x, na.rm = TRUE), 3))) %>%
      mutate(Subtype = "All ISC/EB (651 cells)") %>%
      select(Subtype, everything())
    bind_rows(overall, by_sub)
  })

  output$expr_table <- renderTable({ calc_table() })

  output$dl_table_csv <- downloadHandler(
    filename = function() { "Zipper2022_ISC_EB_Mean_Expression.csv" },
    content = function(file) {
      write.csv(calc_table(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
