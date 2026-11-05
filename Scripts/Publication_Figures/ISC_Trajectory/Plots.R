

setwd("~/Desktop/Drosophila_midgut/Publication_Figures/ISC_Trajectory/")
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
library(ggpubr)
set.seed(10)

#########################################################################################################################################
## Script for Slingshot data analysis
#########################################################################################################################################

outDir <- file.path("./Plots/")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
# Load the integrated data from Seurat and make a SummarizedSinglecellExperiment object using the integrated data

integrated = readRDS("~/Desktop/Drosophila_midgut/Publication_Figures/Metabolic_Genes_Trajectory/Data/integrated.rds")
integrated_celltype = integrated@meta.data$integrated_celltype
simple_celltype = case_when(
  integrated_celltype %in% c("ISC/EB", "ISC") ~ "ISC/EB",
  integrated_celltype %in% c("EE") ~ "EE",
  integrated_celltype %in% c("aEC", "pEC", "mEC", "EC") ~ "EC",
  integrated_celltype %in% c("2") ~ "dEC",
  TRUE ~ "other")
integrated@meta.data$simple_celltype = simple_celltype


# subset for only qISC, ISC, EB
qISC_metadata <- read.csv('~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/qISC_ISC_EB_metadata.csv', row.names = 1)
metadata1 <- integrated@meta.data
rownames(metadata1) <- gsub(":", "\\.", rownames(metadata1))
rownames(metadata1) <- gsub("-", "\\.", rownames(metadata1))
metadata2 <- merge(qISC_metadata, metadata1, by = 0)
rownames(metadata2) <- metadata2$Row.names

#include dEC
#subsetted <- integrated[,colnames(integrated) %in% c(metadata2$cell, integrated@meta.data[integrated@meta.data$simple_celltype == 'dEC',]$cell)]

# don't include dEC
subsetted <- integrated[,colnames(integrated) %in% metadata2$cell]

rownames(metadata2) <- metadata2$cell
metadata2[metadata2$Identity == 'ISC',]$Identity <- 'pISC'
subsetted <- AddMetaData(subsetted, metadata2, col.name = NULL)
table(subsetted@meta.data$Identity)


# slingshot with different end clusters
sce = as.SingleCellExperiment(subsetted)

# load and use metabolism genes

metab_genes <- read.table("~/Desktop/Drosophila_midgut/Publication_Figures/Metabolic_Genes_Trajectory/Data/Metabolic_Genes.tsv", header = F)[,1]
metabGenes <- metab_genes[metab_genes %in% rownames(assays(sce)$logcounts)]
length(metabGenes)

pca1 = prcomp(t(assays(sce)$logcounts[metabGenes,]), scale=FALSE)
pcadat = pca1$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(colData(sce) %>% as.data.frame(), by="cell")

# plot by simple cell type
ggplot(pcadat, aes(PC1, PC2, color=Identity)) + geom_point() +
  scale_color_viridis(discrete=TRUE)


sl1 = slingshot(sce, clusterLabels="Identity", reducedDim="PCA", start.clus="qISC", end.clus='EB')
pal = brewer.pal(9, 'Set1')
pointcolors = pal[as.factor(sce$Identity)]
plot(reducedDims(sce)$PCA, col = pointcolors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sl1), lwd=2, type = 'lineages', col = 'black')
legendlabels = unique(as.factor(sce$Identity))
legendcolors = unique(pointcolors)
legend("topleft", pch=16, legend=legendlabels, col=legendcolors)

ggplot(pcadat, aes(PC1, PC2, color=simple_celltype)) + geom_point() +
  scale_color_viridis(discrete=TRUE)

ggplot(pcadat, aes(PC1, PC2, color=Identity)) + geom_point() +
  scale_color_viridis(discrete=TRUE)
saveRDS(sce, file.path("./trajectory-PCA_ISC_EB.rds"))
saveRDS(sce, file.path("./trajectory-sce_ISC_EB.rds"))
saveRDS(sl1, file.path("./slingshot-trajectory_ISC_EB.rds"))


#########################################################################################################################################
## Script to make Publication Plots
#########################################################################################################################################

# Load data from previous script
pcadat = readRDS("trajectory-PCA_diff.rds")
sce = readRDS("trajectory-sce_diff.rds")
sl1 = readRDS("slingshot-trajectory_diff.rds")

### Make slingshot data set class to hold data

sshot = SlingshotDataSet(sl1)

pseudotimes = data.frame(cell=colnames(sl1),
                         pseudotime1=sl1$slingPseudotime_1,
                         PC_1=reducedDims(sl1)$PCA[, "PC_1"],
                         PC_2=reducedDims(sl1)$PCA[, "PC_2"],
                         Identity=sl1$Identity
)

pcamarkers = data.frame(pseudotime1=sl1$slingPseudotime_1,
                        PC_1=reducedDims(sl1)$PCA[, "PC_1"],
                        PC_2=reducedDims(sl1)$PCA[, "PC_2"],
                        cell=colnames(sl1)) %>%
  left_join(markerdat, by="cell")

plotdata = data.frame(
  cell = rownames(sshot@clusterLabels),
  curve1x=sshot@curves[["curve1"]]$s[, "PC_1"],
  curve1y=sshot@curves[["curve1"]]$s[, "PC_2"],
  celltype=sshot@clusterLabels
) %>% left_join(pseudotimes, by="cell")

curve1 = data.frame(sshot@curves[["curve1"]]$s[sshot@curves[["curve1"]]$ord,])



p = ggplot(plotdata, aes(PC_1, PC_2, color=Identity)) +
  geom_point(size=2.5, alpha=3.8) +
  scale_color_viridis(option = 'viridis', discrete=TRUE) + 
  theme_classic()+
  theme(text=element_text(family="Arial", face="bold", size=12))
p = p + geom_path(data=curve1, col="black", lwd=1, linetype="solid")

p + scale_fill_discrete(limits = c("qISC", "pISC", "EB"))


###### Plot genes of interest

interest_genes <- c("esg", 'Dl', 'klu', "CG4389", "CG31075", "CG3902", "CG4860", "CG9547", "CG5599", "yip2", "whd", "Mtpalpha")

fattydat = logcounts(sce) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  gather("cell", "logcount", -gene) %>%
  filter(gene %in% interest_genes)

fattymarkers = data.frame(pseudotime1=sl1$slingPseudotime_1,
                          PC_1=reducedDims(sl1)$PCA[, "PC_1"],
                          PC_2=reducedDims(sl1)$PCA[, "PC_2"],
                          cell=colnames(sl1)) %>%
  left_join(fattydat, by="cell")

fattymarkers = fattymarkers %>%
  gather(c(pseudotime1), key="lineage", value="pseudotime")
fattymarkers = fattymarkers %>%
  mutate(lineagestring=ifelse(lineage == "pseudotime1", paste0(sshot@lineages$Lineage1, collapse="->"),
                              ifelse(lineage == "pseudotime2", paste0(sshot@lineages$Lineage2, collapse="->"),
                                     ifelse(lineage == "pseudotime3", paste0(sshot@lineages$Lineage3, collapse="->"),
                                            paste0(sshot@lineages$Lineage4, collapse="->")))))




################################




p1 = ggplot(plotdata, aes(PC_1, PC_2, color=Identity)) +
  geom_point(size=2.0, alpha=1.4) +
  scale_color_viridis(option = 'viridis', discrete=TRUE) + 
  #scale_color_manual(values = cbp2)+
  theme_classic()+
  theme(text=element_text(family="Arial", face="bold", size=12))
p1 = p1 + geom_path(data=curve1, col="black", lwd=1.5, linetype="solid")
p1 = p1 + ggtitle('qISC -> pISC -> EB')+coord_equal()

p2 <- ggplot(fattymarkers[fattymarkers$lineagestring == 'qISC->pISC->EB',], aes(pseudotime, logcount, color=lineagestring)) +
  geom_smooth(se=FALSE, color= 'red', size = 2) +
  facet_wrap(~gene, ncol=1, strip.position="top", scale='free_y') +
  theme_pubr() +
  theme(panel.margin=unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing=unit(0.2, "lines"),
        strip.text.y = element_text(angle = 180),
        #axis.text.y=element_text(size=8),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(size=8)) +
  scale_color_viridis(discrete=TRUE, begin=0) +
  ylab("log count") +
  labs(color="lineage")+
  theme(text=element_text(family="Arial", face="bold", size=12))
cowplot::plot_grid(p1, p2, ncol = 1)


svg(paste0(outDir,"ISC_trajectory_plot.svg"), width = 12, height = 18)
cowplot::plot_grid(p1, p2, ncol = 1)
dev.off()

jpeg(paste0(outDir,"ISC_trajectory_plot.jpeg"), width = 8000, height = 12000, res = 600)
cowplot::plot_grid(p1, p2, ncol = 1)
dev.off()
