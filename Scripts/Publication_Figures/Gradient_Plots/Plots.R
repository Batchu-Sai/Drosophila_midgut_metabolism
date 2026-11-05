setwd("~/Desktop/Drosophila_midgut/Publication_Figures/Gradient_Plots/")
rm(list = ls())
library(scater)
library(reshape2)
library(viridis)
library(umap)
library(cowplot)
library(Seurat)
library("extrafont")
library(dplyr)
options(stringsAsFactors=FALSE)
outDir <- file.path("./Plots/")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

### Make Seurat object using ISC and EB data

tsne_final <- read.csv("./Data/tsne_final_ISC_EB.csv", row.names = 1)
obj_seur <- CreateSeuratObject(t(tsne_final))

### Add tSNE reductions from metabolism to Seurat object

df <- data.frame(tSNE_1 = tsne_final$X,tSNE_2 = tsne_final$Y)
rownames(df) <- rownames(tsne_final)
df <- as.matrix(df)
obj_seur[["tsne"]] <- CreateDimReducObject(embeddings = df , key = "tSNE_", assay = DefaultAssay(obj_seur))

### Check out expression of klu, Dl, and esg

# esg+, klu-, Dl- === quiescent ISC
# Dl+, klu- === ISC
# Dl-, klu+ === EB

# subset putative quiescent ISC
quiescent_ISC <- subset(obj_seur, subset = Dl < 1 & klu < 1 & esg > 0)
# subset putative ISC
ISC <- subset(obj_seur, subset = Dl > 0 & klu < 1)
# subset putative EB
EB <- subset(obj_seur, subset = Dl < 1 & klu > 1)

### Plot and change axes to match

p1 <- FeaturePlot(obj_seur, features = c("klu", "Dl"), blend=TRUE, blend.threshold=0.0001, cols=c("#cc5500", "#46cb18"), reduction = 'tsne', pt.size = 1.25, min.cutoff = 0)
p1[[1]] <- p1[[1]]+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
p1[[2]] <-p1[[2]]+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                 axis.text.y=element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x = element_blank())
p1[[3]] <-p1[[3]]+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                 axis.text.y=element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x = element_blank())
p1[[4]] <- p1[[4]]+ggtitle("")

p2 <- FeaturePlot(quiescent_ISC, features = c("esg", "pros"), blend=TRUE, blend.threshold=0.0001, cols=c("red", "blue"), reduction = 'tsne', pt.size = 1.25, min.cutoff = 0)
p2[[1]] <- p2[[1]] + xlim(c(-60,60)) + ylim(c(-20,20))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
p2[[2]] <- p2[[2]] + xlim(c(-60,60)) + ylim(c(-20,20))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
p2[[3]] <- p2[[3]] + xlim(c(-60,60)) + ylim(c(-20,20))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
p2[[4]] <- p2[[4]]+ggtitle("")
cowplot::plot_grid(p1, p2, ncol =1)


svg(paste0(outDir,"Gradient_plot.svg"), width = 12, height = 8)
cowplot::plot_grid(p1, p2, ncol =1)
dev.off()

jpeg(paste0(outDir,"Gradient_plot.jpeg"), width = 8000, height = 4000, res = 600)
cowplot::plot_grid(p1, p2, ncol =1)
dev.off()
