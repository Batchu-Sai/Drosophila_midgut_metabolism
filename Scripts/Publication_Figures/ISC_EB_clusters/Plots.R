setwd("~/Desktop/Drosophila_midgut/Publication_Figures/ISC_EB_clusters/")
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
## Script for ISC/EB clusters and gene expression plots
#########################################################################################################################################

outDir <- file.path("./Plots/")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

tsne_final <- read.csv("./Data/tsne_final_ISC_EB.csv", row.names = 1)

genes_interest <- c("esg", "Dl", "klu", "Eip75B", 
                    "CG31075", "CG3902", "CG4860", "CG5599", 
                    "CG9547", "Mtpalpha", "whd", "yip2"
                    )

p_1.1 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'esg')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('esg')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_1.2 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'Dl')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('Dl')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_1.3 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'klu')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('klu')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_1.4 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'Eip75B')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('Eip75B')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_1.5 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'CG31075')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  ggtitle('CG31075')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_1.6 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'CG3902')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('CG3902')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_1.7 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'CG4860')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('CG4860')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_1.8 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'CG5599')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('CG5599')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_1.9 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'CG9547')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('CG9547')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_2.0 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'Mtpalpha')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('Mtpalpha')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_2.1 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'whd')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('whd')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

p_2.2 <- ggplot(tsne_final) + 
  geom_point(aes(as.numeric(X), as.numeric(Y), colour =tsne_final[[which(colnames(tsne_final) == 'yip2')]]), size = .5) +
  labs(x = "tSNE1",y = "tSNE2") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle('yip2')+
  scale_color_viridis(option = 'inferno')+
  labs(color='log2(TPM)')+
  theme(text=element_text(family="Arial", face="bold", size=7))

svg(paste0(outDir,"ISC_EB_expression_plot.svg"), width = 12, height = 8)
cowplot::plot_grid(p_1.1, p_1.2, p_1.3, p_1.4, p_1.5, p_1.6, p_1.7, p_1.8, p_1.9, p_2.0, p_2.1, p_2.2, ncol = 4)
dev.off()


jpeg(paste0(outDir,"ISC_EB_expression_plot.jpeg"), width = 9000, height = 7000, res = 600)
cowplot::plot_grid(p_1.1, p_1.2, p_1.3, p_1.4, p_1.5, p_1.6, p_1.7, p_1.8, p_1.9, p_2.0, p_2.1, p_2.2, ncol = 4)
dev.off()
