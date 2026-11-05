setwd("~/Desktop/Drosophila_midgut/Publication_Figures/Metabolic_Genes_Trajectory/")
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

integrated = readRDS("./Data/integrated.rds")
integrated_celltype = integrated@meta.data$integrated_celltype
simple_celltype = case_when(
  integrated_celltype %in% c("ISC/EB", "ISC") ~ "ISC/EB",
  integrated_celltype %in% c("EE") ~ "EE",
  integrated_celltype %in% c("aEC", "pEC", "mEC", "EC") ~ "EC",
  integrated_celltype %in% c("2") ~ "2",
  TRUE ~ "other")
integrated@meta.data$simple_celltype = simple_celltype
sce = as.SingleCellExperiment(integrated)

# Select cells in the aEC, mEC, pEC, dEC ('2'), ISC/EB, EE

as.data.frame.table(table(colData(sce)$integrated_celltype)) %>%
  dplyr::rename(celltype=Var1, ncells=Freq) %>%
  knitr::kable()

keep_celltypes = c("aEC", "EE", "ISC", "ISC/EB","mEC", "pEC", "2")
keep_cells = colData(sce) %>%
  as.data.frame() %>%
  dplyr::filter(integrated_celltype %in% keep_celltypes) %>%
  pull(cell)
sce = sce[, keep_cells]
as.data.frame.table(table(colData(sce)$integrated_celltype)) %>%
  dplyr::rename(celltype=Var1, ncells=Freq) %>%
  knitr::kable()

# load and use metabolism genes

metab_genes <- read.table("./Data/Metabolic_Genes.tsv", header = F)[,1]
metabGenes <- metab_genes[metab_genes %in% rownames(assays(sce)$logcounts)]
length(metabGenes)

pca1 = prcomp(t(assays(sce)$logcounts[metabGenes,]), scale=FALSE)
pcadat = pca1$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(colData(sce) %>% as.data.frame(), by="cell")

# plot by simple cell type
ggplot(pcadat, aes(PC1, PC2, color=simple_celltype)) + geom_point() +
  scale_color_viridis(discrete=TRUE)

# plot by integrated cell type

ggplot(pcadat, aes(PC1, PC2, color=integrated_celltype)) + geom_point() +
  scale_color_viridis(discrete=TRUE)

# slingshot with different end clusters

sl1 = slingshot(sce, clusterLabels="integrated_celltype", reducedDim="PCA", start.clus="ISC/EB",
                end.clus=c("EE", 'mEC', 'aEC', 'pEC'))
pal = brewer.pal(9, 'Set1')
pointcolors = pal[as.factor(sce$simple_celltype)]
plot(reducedDims(sce)$PCA, col = pointcolors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sl1), lwd=2, type = 'lineages', col = 'black')
legendlabels = unique(as.factor(sce$simple_celltype))
legendcolors = unique(pointcolors)
legend("topleft", pch=16, legend=legendlabels, col=legendcolors)

ggplot(pcadat, aes(PC1, PC2, color=simple_celltype)) + geom_point() +
  scale_color_viridis(discrete=TRUE)

ggplot(pcadat, aes(PC1, PC2, color=integrated_celltype)) + geom_point() +
  scale_color_viridis(discrete=TRUE)
saveRDS(sce, file.path("./trajectory-PCA_diff.rds"))
saveRDS(sce, file.path("./trajectory-sce_diff.rds"))
saveRDS(sl1, file.path("./slingshot-trajectory_diff.rds"))


#########################################################################################################################################
## Script to make Publication Plots
#########################################################################################################################################

# Load data from previous script
pcadat = readRDS("trajectory-PCA_diff.rds")
sce = readRDS("trajectory-sce_diff.rds")
sl1 = readRDS("slingshot-trajectory_diff.rds")

# Name Lineages
colnames(sl1@int_metadata$slingshot@clusterLabels)[2] <- 'dEC'
sl1@int_metadata$slingshot@lineages[[1]] <- c("ISC/EB", "dEC", "pEC")
sl1@int_metadata$slingshot@lineages[[2]] <- c("ISC/EB", "dEC", "mEC")
sl1@int_metadata$slingshot@lineages[[3]] <- c("ISC/EB", "dEC", "aEC")
sl1@int_metadata$slingshot@lineages[[4]] <- c("ISC/EB", "dEC", "EE")
sl1@colData$integrated_celltype[sl1@colData$integrated_celltype == '2'] <- 'dEC'

# Plot Pseudotime gene expression plots
toplot = c("Dl", "esg", "klu", "E(spl)m3-HLH", "E(spl)malpha-BFM", "E(spl)mbeta-HLH", "sna", "stg", "Piezo",
           "mesh", "ssk", "pros", "betaTry", "lambdaTry", "lab")

markerdat = logcounts(sce) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  gather("cell", "logcount", -gene) %>%
  filter(gene %in% toplot)

pseudotimes = data.frame(cell=colnames(sl1),
                         pseudotime1=sl1$slingPseudotime_1,
                         pseudotime2=sl1$slingPseudotime_2,
                         pseudotime3=sl1$slingPseudotime_3,
                         PC_1=reducedDims(sl1)$PCA[, "PC_1"],
                         PC_2=reducedDims(sl1)$PCA[, "PC_2"],
                         integrated_celltype=sl1$integrated_celltype
)

pcamarkers = data.frame(pseudotime1=sl1$slingPseudotime_1,
                        pseudotime2=sl1$slingPseudotime_2,
                        pseudotime3=sl1$slingPseudotime_3,
                        PC_1=reducedDims(sl1)$PCA[, "PC_1"],
                        PC_2=reducedDims(sl1)$PCA[, "PC_2"],
                        cell=colnames(sl1)) %>%
  left_join(markerdat, by="cell")

### Make slingshot data set class to hold data

sshot = SlingshotDataSet(sl1)

plotdata = data.frame(
  cell = rownames(sshot@clusterLabels),
  curve1x=sshot@curves[["curve1"]]$s[, "PC_1"],
  curve1y=sshot@curves[["curve1"]]$s[, "PC_2"],
  celltype=sshot@clusterLabels
) %>% left_join(pseudotimes, by="cell")

curve1 = data.frame(sshot@curves[["curve1"]]$s[sshot@curves[["curve1"]]$ord,])
curve2 = data.frame(sshot@curves[["curve2"]]$s[sshot@curves[["curve2"]]$ord,])
curve3 = data.frame(sshot@curves[["curve3"]]$s[sshot@curves[["curve3"]]$ord,])
curve4 = data.frame(sshot@curves[["curve4"]]$s[sshot@curves[["curve4"]]$ord,])


plotdata$integrated_celltype <- factor(plotdata$integrated_celltype, levels=c("ISC/EB", "dEC", "aEC", 'mEC', 'pEC', 'EE'), labels = c("ISC/EB", "dEC", "aEC", 'mEC', 'pEC', 'EE'))

cbp2 <- c("#009E73", "#D55E00","#293352", "#4E84C4", "#56B4E9", "#CC79A7")


p = ggplot(plotdata, aes(PC_1, PC_2, color=integrated_celltype)) +
  geom_point(size=0.5, alpha=0.8) +
  scale_color_viridis(option = 'viridis', discrete=TRUE) + 
  theme_classic()+
  theme(text=element_text(family="Arial", face="bold", size=12))
p = p + geom_path(data=curve1, col="black", lwd=1, linetype="solid")
p = p + geom_path(data=curve2, col="black", lwd=1, linetype="solid")
p = p + geom_path(data=curve3, col="black", lwd=1, linetype="solid")
p = p + geom_path(data=curve4, col="black", lwd=1, linetype="solid")
p + scale_fill_discrete(limits = c("ISC/EB", "dEC", "aEC", 'mEC', 'pEC'))

target_genes = c("CG4389", "CG31075", "CG3902", "CG4860", "CG9547", "CG5599", "yip2", "whd", "Mtpalpha") #CG4389 is same as Mtpalpha
fatty = c("CG4389", "CG31075", "CG3902", "CG4860", "CG9547", "CG5599", "yip2", "whd", "Mtpalpha", "esg", 'Dl', 'klu', 'E(spl)malpha-BFM', 'E(spl)mbeta-HLH', 'E(spl)m3-HLH')

interest_genes <- c("esg", 'Dl', 'klu', "CG4389", "CG31075", "CG3902", "CG4860", "CG9547", "CG5599", "yip2", "whd", "Mtpalpha")

fattydat = logcounts(sce) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  gather("cell", "logcount", -gene) %>%
  filter(gene %in% interest_genes)

fattymarkers = data.frame(pseudotime1=sl1$slingPseudotime_1,
                          pseudotime2=sl1$slingPseudotime_2,
                          pseudotime3=sl1$slingPseudotime_3,
                          pseudotime4=sl1$slingPseudotime_4,
                          PC_1=reducedDims(sl1)$PCA[, "PC_1"],
                          PC_2=reducedDims(sl1)$PCA[, "PC_2"],
                          cell=colnames(sl1)) %>%
  left_join(fattydat, by="cell")

ggplot(fattymarkers %>%
         filter(gene == "Dl"), aes(PC_1, PC_2, color=logcount)) +
  geom_point(size=0.5, alpha=0.1) +
  scale_color_viridis() + 
  geom_path(data=curve1, col="black", lwd=1, linetype="solid") +
  geom_path(data=curve2, col="black", lwd=1, linetype="solid") +
  geom_path(data=curve3, col="black", lwd=1, linetype="solid") +
  geom_path(data=curve4, col="black", lwd=1, linetype="solid")


fattymarkers = fattymarkers %>%
  gather(c(pseudotime1, pseudotime2, pseudotime3, pseudotime4), key="lineage", value="pseudotime")
fattymarkers = fattymarkers %>%
  mutate(lineagestring=ifelse(lineage == "pseudotime1", paste0(sshot@lineages$Lineage1, collapse="->"),
                              ifelse(lineage == "pseudotime2", paste0(sshot@lineages$Lineage2, collapse="->"),
                                     ifelse(lineage == "pseudotime3", paste0(sshot@lineages$Lineage3, collapse="->"),
                                            paste0(sshot@lineages$Lineage4, collapse="->")))))
ggplot(fattymarkers, aes(pseudotime, logcount, group=gene, color=gene)) +
  geom_smooth(se=FALSE, lwd=0.5) +
  facet_wrap(~lineagestring, ncol=1)


p2 <-ggplot(fattymarkers[fattymarkers$lineagestring == 'ISC/EB->EE' | fattymarkers$lineagestring == 'ISC/EB->dEC->EE' ,], aes(pseudotime, logcount, color=lineagestring)) +
  geom_smooth(se=FALSE) +
  facet_wrap(~gene, ncol=1, strip.position="left", scale='free_y') +
  theme_pubr() +
  theme(panel.margin=unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing=unit(0.2, "lines"),
        strip.text.y = element_text(angle = 180),
        axis.text.y=element_text(size=8)) +
  scale_color_viridis(discrete=TRUE, begin=0) +
  ylab("log count") +
  labs(color="lineage")+
  theme(text=element_text(family="Arial", face="bold", size=8))

colnames(plotdata)[[15]] <- 'Celltype'

#### ISC/EB -> dEC -> aEC

p1 = ggplot(plotdata, aes(PC_1, PC_2, color=Celltype)) +
  geom_point(size=0.5, alpha=1.4) +
  #scale_color_viridis(option = 'viridis', discrete=TRUE) + 
  scale_color_manual(values = cbp2)+
  theme_classic()+
  theme(text=element_text(family="Arial", face="bold", size=12))
p1 = p1 + geom_path(data=curve3, col="black", lwd=1.5, linetype="solid")
p1 = p1 + ggtitle('ISC/EB -> dEC -> aEC')+coord_equal()

p2 <- ggplot(fattymarkers[fattymarkers$lineagestring == 'ISC/EB->dEC->aEC',], aes(pseudotime, logcount, color=lineagestring)) +
  geom_smooth(se=FALSE, color= 'red', size = 1.5) +
  facet_wrap(~gene, ncol=1, strip.position="top", scale='free_y') +
  theme_pubr() +
  theme(panel.margin=unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing=unit(0.2, "lines"),
        strip.text.y = element_text(angle = 180),
        #axis.text.y=element_text(size=8),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(size=12)) +
  scale_color_viridis(discrete=TRUE, begin=0) +
  ylab("log count") +
  labs(color="lineage")+
  theme(text=element_text(family="Arial", face="bold", size=15))
p2
aEC <- cowplot::plot_grid(p1, p2, ncol = 1)

#### ISC/EB -> dEC -> mEC

p3 = ggplot(plotdata, aes(PC_1, PC_2, color=Celltype)) +
  geom_point(size=0.5, alpha=1.4) +
  #scale_color_viridis(option = 'viridis', discrete=TRUE) + 
  scale_color_manual(values = cbp2)+
  theme_classic()+
  theme(text=element_text(family="Arial", face="bold", size=12))
p3 = p3 + geom_path(data=curve2, col="black", lwd=1.5, linetype="solid")
p3 = p3 + ggtitle('ISC/EB -> dEC -> mEC')+coord_equal()

p4 <- ggplot(fattymarkers[fattymarkers$lineagestring == 'ISC/EB->dEC->mEC',], aes(pseudotime, logcount, color=lineagestring)) +
  geom_smooth(se=FALSE, color= 'red', size = 1.5) +
  facet_wrap(~gene, ncol=1, strip.position="top", scale='free_y') +
  theme_pubr() +
  theme(panel.margin=unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing=unit(0.2, "lines"),
        strip.text.y = element_text(angle = 180),
        #axis.text.y=element_text(size=8),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(size=12)) +
  scale_color_viridis(discrete=TRUE, begin=0) +
  ylab("log count") +
  labs(color="lineage")+
  theme(text=element_text(family="Arial", face="bold", size=15))

mEC <- cowplot::plot_grid(p3, p4, ncol = 1)

#### ISC/EB -> dEC -> pEC

p5 = ggplot(plotdata, aes(PC_1, PC_2, color=Celltype)) +
  geom_point(size=0.5, alpha=1.4) +
  #scale_color_viridis(option = 'viridis', discrete=TRUE) + 
  scale_color_manual(values = cbp2)+ 
  theme_classic()+
  theme(text=element_text(family="Arial", face="bold", size=12))+
  geom_path(data=curve1, col="black", lwd=1.5, linetype="solid")+
  ggtitle('ISC/EB -> dEC -> pEC')+coord_equal()

p6 <- ggplot(fattymarkers[fattymarkers$lineagestring == 'ISC/EB->dEC->pEC',], aes(pseudotime, logcount, color=lineagestring)) +
  geom_smooth(se=FALSE, color= 'red', size = 1.5) +
  facet_wrap(~gene, ncol=1, strip.position="top", scale='free_y') +
  theme_pubr() +
  theme(panel.margin=unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing=unit(0.2, "lines"),
        strip.text.y = element_text(angle = 180),
        #axis.text.y=element_text(size=8),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(size=12)) +
  scale_color_viridis(discrete=TRUE, begin=0) +
  ylab("log count") +
  labs(color="lineage")+
  theme(text=element_text(family="Arial", face="bold", size=15))

pEC <- cowplot::plot_grid(p5, p6, ncol = 1)


#### ISC/EB -> dEC -> EE

p7 = ggplot(plotdata, aes(PC_1, PC_2, color=Celltype)) +
  geom_point(size=0.5, alpha=1.4) +
  #scale_color_viridis(option = 'viridis', discrete=TRUE) + 
  scale_color_manual(values = cbp2)+ 
  theme_classic()+
  theme(text=element_text(family="Arial", face="bold", size=12))+
  geom_path(data=curve4, col="black", lwd=1.5, linetype="solid")+
  ggtitle('ISC/EB -> dEC -> EE')+coord_equal()

p8 <- ggplot(fattymarkers[fattymarkers$lineagestring == 'ISC/EB->dEC->EE',], aes(pseudotime, logcount, color=lineagestring)) +
  geom_smooth(se=FALSE, color= 'red', size = 1.5) +
  facet_wrap(~gene, ncol=1, strip.position="top", scale='free_y') +
  theme_pubr() +
  theme(panel.margin=unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing=unit(0.2, "lines"),
        strip.text.y = element_text(angle = 180),
        #axis.text.y=element_text(size=8),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(size=12)) +
  scale_color_viridis(discrete=TRUE, begin=0) +
  ylab("log count") +
  labs(color="lineage")+
  theme(text=element_text(family="Arial", face="bold", size=15))

EE <- cowplot::plot_grid(p7, p8, ncol = 1)


svg(paste0(outDir,"trajectory_plot_changed_fonts.svg"), width = 25, height = 15)
cowplot::plot_grid(aEC, pEC,mEC, EE, ncol = 4)
dev.off()
dev.off()

jpeg(paste0(outDir,"trajectory_plot.jpeg"), width = 12000, height = 9000, res = 600)
cowplot::plot_grid(aEC, pEC,mEC, EE, ncol = 4)
dev.off()


