setwd("~/Desktop/Drosophila_midgut/Publication_Figures/TSNE_plots/")
rm(list = ls())
library(mclust)
library(Seurat)
library(Matrix)
library(tidyverse)
library(viridis)
library(broom)
library(SingleCellExperiment)
library(slingshot)
library(RColorBrewer)
library("extrafont")
library(cowplot)
library(ggpubr)
set.seed(10)

#########################################################################################################################################
## Script for TSNE plots of all cells for only metabolism genes, and all genes of all cells
#########################################################################################################################################

outDir <- file.path("./Plots/")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

# Color palette
cbp2 <- c("#009E73", "#D55E00","#293352", "#4E84C4", "#56B4E9", "#CC79A7")

# load in tSNE coordinates
all_tsne <- read.csv("./Data/All_Cells/all_tsne.csv", row.names = 1)
all_tsne$group <- factor(all_tsne$group, levels=c("ISC/EB", "dEC", "aEC", 'mEC', 'pEC', 'EE'), labels = c("ISC/EB", "dEC", "aEC", 'mEC', 'pEC', 'EE'))

g1 <- ggplot(all_tsne) + geom_point(aes(x, y, colour = group), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  scale_color_manual(values = cbp2)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("All genes")+
  theme(text=element_text(family="Arial", face="bold", size=8))
g1


# load in tSNE coordinates
metab_tsne <- read.csv("./Data/All_Cells/metabolic_tsne.csv", row.names = 1)
metab_tsne$group <- factor(metab_tsne$group, levels=c("ISC/EB", "dEC", "aEC", 'mEC', 'pEC', 'EE'), labels = c("ISC/EB", "dEC", "aEC", 'mEC', 'pEC', 'EE'))


g2 <- ggplot(metab_tsne) + geom_point(aes(x, y, colour = group), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  scale_color_manual(values = cbp2)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Metabolic genes")+
  theme(text=element_text(family="Arial", face="bold", size=8))
g2


svg(paste0(outDir,"All_cells_tSNE_plot.svg"), width = 8, height = 4)
cowplot::plot_grid(g1, g2, ncol = 2)
dev.off()

jpeg(paste0(outDir,"All_cells_tSNE_plot.jpeg"), width = 8000, height = 4000, res = 600)
cowplot::plot_grid(g1, g2, ncol = 2)
dev.off()



#########################################################################################################################################
## Script for TSNE plots of each cell type for only metabolism genes
#########################################################################################################################################

types <- c("ISC_EB", "dEC", "aEC", 'mEC', 'pEC', 'EE')
names(cbp2) <- c("ISC_EB", "dEC", "aEC", 'mEC', 'pEC', 'EE')
plots_list <- list()

for (ii in types){
  coordinates <- read.csv(paste0("./Data/Cell_Type_Specific/",ii,"_tSNE_output.csv"), row.names = 1)
  plots_list[[ii]] <- ggplot(coordinates) + 
    geom_point(aes(x, y), color = cbp2[ii], size = .5) +
    labs(x = "tSNE1",y = "tSNE2") +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank())+
    ggtitle(paste0(ii))+
    theme(text=element_text(family="Arial", face="bold", size=8))
}

cowplot::plot_grid(plotlist = plots_list, ncol = 3)


svg(paste0(outDir,"Cell_type_tSNE_plot.svg"), width = 12, height = 8)
cowplot::plot_grid(plotlist = plots_list, ncol = 3)
dev.off()

jpeg(paste0(outDir,"Cell_type__tSNE_plot.jpeg"), width = 8000, height = 4000, res = 600)
cowplot::plot_grid(plotlist = plots_list, ncol = 3)
dev.off()

