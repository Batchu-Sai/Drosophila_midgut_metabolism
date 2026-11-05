rm(list=ls())
library(stringr)
library(reshape2)
library(scales)
library(scater)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyverse)
library(RColorBrewer)
library("extrafont")
setwd("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps")
pathway_file <- "./KEGG_dros_pathways.gmt"
source("./utils.R")


################################################################################################################################################ 
###### First heatmap using quiescent ISC, and ISC, and EB
################################################################################################################################################ 
outDir <- file.path("./Heatmap_qISC_ISC_EB")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

#1. Loading the data, replace the metadata with one from dissection analysis, and keep only cell types of interest
selected_impute_sce <- readRDS("~/Desktop/Drosophila_midgut/metabolism/Figures/selected_impute_sce.rds")
colData(selected_impute_sce)$integrated_celltype[colData(selected_impute_sce)$integrated_celltype == '2'] <- 'dEC'
qISC_metadata <- read.csv('~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/qISC_ISC_EB_metadata.csv', row.names = 1)
new_coldata <- colData(selected_impute_sce)
new_coldata$Identity <- new_coldata$integrated_celltype
new_coldata[rownames(new_coldata) %in% (rownames(qISC_metadata)[qISC_metadata$Identity == 'qISC']),]$Identity <- 'qISC'
new_coldata[rownames(new_coldata) %in% (rownames(qISC_metadata)[qISC_metadata$Identity == 'EB']),]$Identity <- 'EB'
new_coldata[rownames(new_coldata) %in% (rownames(qISC_metadata)[qISC_metadata$Identity == 'ISC']),]$Identity <- 'ISC'
table(new_coldata$Identity)
colData(selected_impute_sce) <- new_coldata
selected_impute_sce <- selected_impute_sce[,colData(selected_impute_sce)$Identity %in% c('dEC', 'EE', 'aEC', 'mEC', 'pEC', 'qISC', 'ISC', 'EB')]
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
all_cell_types <- as.vector(selected_impute_sce$Identity)
cell_types <- unique(all_cell_types)

#some genes occur in multiple pathways.
gene_pathway_number <- num_of_pathways(pathway_file,rownames(selected_impute_sce)[rowData(selected_impute_sce)$metabolic])

set.seed(123)
norm_tpm <- readRDS("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/Deconvolution_TPM/Deconvolution_tpm.rds")

norm_tpm <- norm_tpm[,colnames(norm_tpm) %in% rownames(colData(selected_impute_sce))]

##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next
  
  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  
  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]
  
  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
  
  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
  
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
  names(shuffle_cell_types_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- cell_types
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 5000 
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 5000
    }
    if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}
all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]


######### Draw heatmap

dat <- mean_expression_shuffle
saveRDS(dat, file = file.path(outDir, 'heatmap_qISC.rds'))

sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1
mybreaks <- c(
  seq(0, 0.5, length.out=33),
  seq(0.51, 1.5, length.out=33),
  seq(1.51, max(dat),length.out=34)
) 
color <- colorRampPalette(c("#0072B2","white","#D55E00"))(100)
jpeg(file = file.path(outDir,"heatmap_qISC.jpeg"), width = 5000, height = 6000, res = 600)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,color=color,fontsize = 9, breaks = mybreaks,
         labels_row = make_bold_names(dat[sort_row,sort_column], rownames, c(rownames(dat[sort_row,sort_column]))),
         labels_col = make_bold_names(dat[sort_row,sort_column], colnames, c(colnames(dat[sort_row,sort_column]))))
dev.off()

#write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
#write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
#write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue.txt"),row.names=T,col.names=T,quote=F,sep="\t")

#boxplot show the distribution of pathway activity
scRNA_dat <- as.data.frame(mean_expression_noshuffle)
scRNA_dat$X <- NULL

scRNA_df <- melt(scRNA_dat)
scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]
g <- ggplot(scRNA_df,aes(x=variable,y=value,fill=variable)) +
  scale_y_continuous(limits=c(0,4),breaks=0:4,labels=0:4)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=1,color="blue")+
  scale_fill_brewer(palette="Set2") +
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 6),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))

jpeg(file.path(outDir,"pathway_activity_violinplot_qISC.jpeg"),width = 1000, height = 1200, res = 600)
g
dev.off()





################################################################################################################################################ 
###### Second heatmap using ISC, and EC, and EB, and dEC
################################################################################################################################################ 
rm(list = ls())
outDir <- file.path("./Heatmap_ISC_EC_EB")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
source("./utils.R")
pathway_file <- "./KEGG_dros_pathways.gmt"

#1. Loading the data, replace the metadata with one from dissection analysis, and keep only cell types of interest
selected_impute_sce <- readRDS("~/Desktop/Drosophila_midgut/metabolism/Figures/selected_impute_sce.rds")
colData(selected_impute_sce)$integrated_celltype[colData(selected_impute_sce)$integrated_celltype == '2'] <- 'dEC'
ISC_metadata <- read.csv('~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/ISC_EB_EC_metadata.csv', row.names = 1)
new_coldata <- colData(selected_impute_sce)
new_coldata$Identity <- new_coldata$integrated_celltype
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity == 'ISC']),]$Identity <- 'ISC'
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity == 'EB']),]$Identity <- 'EB'
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity == 'dEC']),]$Identity <- 'dEC'
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity == 'EC']),]$Identity <- 'EC'
table(new_coldata$Identity)
colData(selected_impute_sce) <- new_coldata
selected_impute_sce <- selected_impute_sce[,colData(selected_impute_sce)$Identity %in% c('dEC', 'EE', 'EC', 'ISC', 'EB')]
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
all_cell_types <- as.vector(selected_impute_sce$Identity)
cell_types <- unique(all_cell_types)

#some genes occur in multiple pathways.
gene_pathway_number <- num_of_pathways(pathway_file,rownames(selected_impute_sce)[rowData(selected_impute_sce)$metabolic])

set.seed(123)
norm_tpm <- readRDS("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/Deconvolution_TPM/Deconvolution_tpm.rds")

norm_tpm <- norm_tpm[,colnames(norm_tpm) %in% rownames(colData(selected_impute_sce))]

##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next
  
  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  
  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]
  
  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
  
  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
  
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
  names(shuffle_cell_types_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- cell_types
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 5000 
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 5000
    }
    if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}
all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]


######### Draw heatmap

dat <- mean_expression_shuffle
saveRDS(dat, file = file.path(outDir, 'heatmap_ISC_EB_allEC.rds'))

sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1
mybreaks <- c(
  seq(0, 0.5, length.out=33),
  seq(0.51, 1.5, length.out=33),
  seq(1.51, max(dat),length.out=34)
) 
color <- colorRampPalette(c("#0072B2","white","#D55E00"))(100)
jpeg(file = file.path(outDir,"heatmap_ISC_EB_allEC.jpeg"), width = 5000, height = 6000, res = 600)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,color=color,fontsize = 9, breaks = mybreaks,
         labels_row = make_bold_names(dat[sort_row,sort_column], rownames, c(rownames(dat[sort_row,sort_column]))),
         labels_col = make_bold_names(dat[sort_row,sort_column], colnames, c(colnames(dat[sort_row,sort_column]))))
dev.off()

#write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
#write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
#write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue.txt"),row.names=T,col.names=T,quote=F,sep="\t")

#boxplot show the distribution of pathway activity
scRNA_dat <- as.data.frame(mean_expression_noshuffle)
scRNA_dat$X <- NULL

scRNA_df <- melt(scRNA_dat)
scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]
g <- ggplot(scRNA_df,aes(x=variable,y=value,fill=variable)) +
  scale_y_continuous(limits=c(0,4),breaks=0:4,labels=0:4)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=1,color="blue")+
  scale_fill_brewer(palette="Set2") +
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 6),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))

jpeg(file.path(outDir,"pathway_activity_violinplot_ISC_EB_allEC.jpeg"),width = 1000, height = 1200, res = 600)
g
dev.off()


################################################################################################################################################ 
###### Third heatmap using ISC, and each EC, and EB, and dEC
################################################################################################################################################ 
rm(list = ls())
setwd("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps")
pathway_file <- "./KEGG_dros_pathways.gmt"
source("./utils.R")
outDir <- file.path("./Heatmap_ISC_eachEC_EB")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

#1. Loading the data, replace the metadata with one from dissection analysis, and keep only cell types of interest
selected_impute_sce <- readRDS("~/Desktop/Drosophila_midgut/metabolism/Figures/selected_impute_sce.rds")
colData(selected_impute_sce)$integrated_celltype[colData(selected_impute_sce)$integrated_celltype == '2'] <- 'dEC'
ISC_metadata <- read.csv('~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/ISC_EB_EC_metadata.csv', row.names = 1)
new_coldata <- colData(selected_impute_sce)
new_coldata$Identity <- new_coldata$integrated_celltype
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity2 == 'ISC']),]$Identity <- 'ISC'
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity2 == 'EB']),]$Identity <- 'EB'
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity2 == 'dEC']),]$Identity <- 'dEC'
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity2 == 'aEC']),]$Identity <- 'aEC'
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity2 == 'mEC']),]$Identity <- 'mEC'
new_coldata[rownames(new_coldata) %in% (rownames(ISC_metadata)[ISC_metadata$Identity2 == 'pEC']),]$Identity <- 'pEC'
table(new_coldata$Identity)
colData(selected_impute_sce) <- new_coldata
selected_impute_sce <- selected_impute_sce[,colData(selected_impute_sce)$Identity %in% c('dEC', 'EE', 'aEC','mEC','pEC', 'ISC', 'EB')]
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
all_cell_types <- as.vector(selected_impute_sce$Identity)
cell_types <- unique(all_cell_types)

#some genes occur in multiple pathways.
gene_pathway_number <- num_of_pathways(pathway_file,rownames(selected_impute_sce)[rowData(selected_impute_sce)$metabolic])

set.seed(123)
norm_tpm <- readRDS("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/Deconvolution_TPM/Deconvolution_tpm.rds")

norm_tpm <- norm_tpm[,colnames(norm_tpm) %in% rownames(colData(selected_impute_sce))]

##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next
  
  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  
  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]
  
  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
  
  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
  
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
  names(shuffle_cell_types_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- cell_types
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 5000
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 5000
    }
    if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}
all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]


######### Draw heatmap

dat <- mean_expression_shuffle
saveRDS(dat, file = file.path(outDir, 'heatmap_ISC_EB_eachEC.rds'))
dat <- readRDS('./Heatmap_ISC_eachEC_EB/heatmap_ISC_EB_eachEC.rds')
sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1

### Sort by pathway scores for qISC and reorder cell types
rownames(dat) <- gsub("_", " ", rownames(dat))
sort_row <- rownames(dat)[order(dat[,'ISC'],decreasing = T)]
sort_column <- c('ISC', 'EB', 'dEC', 'aEC', 'mEC', 'pEC', 'EE')

mybreaks <- c(
  seq(0, 0.5, length.out=33),
  seq(0.51, 1.5, length.out=33),
  seq(1.51, max(dat),length.out=34)
) 
color <- colorRampPalette(c("#0072B2","white","#D55E00"))(100)
svg("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/Heatmap_ISC_eachEC_EB/heatmap_ISC_EB_eachEC.svg", width = 8, height = 14)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,color=color,fontsize = 9, breaks = mybreaks,
         labels_row = make_bold_names(dat[sort_row,sort_column], rownames, c(rownames(dat[sort_row,sort_column]))),
         labels_col = make_bold_names(dat[sort_row,sort_column], colnames, c(colnames(dat[sort_row,sort_column]))))
dev.off()

write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue.txt"),row.names=T,col.names=T,quote=F,sep="\t")

#boxplot show the distribution of pathway activity
scRNA_dat <- as.data.frame(mean_expression_noshuffle)
scRNA_dat$X <- NULL

scRNA_df <- melt(scRNA_dat)
scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]
scRNA_df$variable <- factor(scRNA_df$variable, levels=c('ISC', 'EB', 'dEC', 'aEC', 'mEC', 'pEC', 'EE'), labels = c('ISC', 'EB', 'dEC', 'aEC', 'mEC', 'pEC', 'EE'))
# Color palette
cbp2 <- c("#009E73", "#D55E00","#293352", "#4E84C4", "#56B4E9", "#CC79A7")


g <- ggplot(scRNA_df,aes(x=variable,y=value,fill=variable)) +
  scale_y_continuous(limits=c(0,4),breaks=0:4,labels=0:4)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=1,color="black")+
  theme_classic() + 
  scale_color_manual(values = cbp2)+
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 6),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))

jpeg(file.path(outDir,"pathway_activity_violinplot_ISC_EB_eachEC.jpeg"),width = 1000, height = 1200, res = 600)
g
dev.off()


################################################################################################################################################ 
###### Heatmap using qISC, ISC, all ISC, EB, dEC, all EC, aEC, mEC, pEC, EE
################################################################################################################################################ 
setwd("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps")
pathway_file <- "./KEGG_dros_pathways.gmt"
source("./utils.R")
rm(list = ls())
outDir <- file.path("./Heatmap_combined_heatmap")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

#1. Loading the data, replace the metadata with one from dissection analysis, and keep only cell types of interest
selected_impute_sce <- readRDS("~/Desktop/Drosophila_midgut/metabolism/Figures/selected_impute_sce.rds")
colData(selected_impute_sce)$integrated_celltype[colData(selected_impute_sce)$integrated_celltype == '2'] <- 'dEC'
qISC_metadata <- read.csv('~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/qISC_ISC_EB_metadata.csv', row.names = 1)
new_coldata <- colData(selected_impute_sce)
new_coldata$Identity <- new_coldata$integrated_celltype
new_coldata[rownames(new_coldata) %in% (rownames(qISC_metadata)[qISC_metadata$Identity == 'qISC']),]$Identity <- 'qISC'
new_coldata[rownames(new_coldata) %in% (rownames(qISC_metadata)[qISC_metadata$Identity == 'EB']),]$Identity <- 'EB'
new_coldata[rownames(new_coldata) %in% (rownames(qISC_metadata)[qISC_metadata$Identity == 'ISC']),]$Identity <- 'prolif_ISC'

# Make new metadata with qISC and ISC cells, add different '_1' cell ending and add to metadata, 

colData_copy <- new_coldata
rownames(colData_copy) <- paste0(rownames(colData_copy), "_1")
head(rownames(colData_copy))
colData_copy[colData_copy$Identity %in% c("qISC", "prolif_ISC"),]$Identity <- 'ISC_all'
colData_copy[colData_copy$Identity %in% c("aEC", "mEC", "pEC"),]$Identity <- 'EC_all'
colData_copy[!(colData_copy$Identity %in% c('ISC_all', 'EC_all')),]$Identity <- 'Duplicate'
new_coldata <- rbind(new_coldata, colData_copy)
table(new_coldata$Identity)

# Make duplicate of selected_impute_sce data, and append '_1' and make new scater object, so new metadata can match
tpm_original <- tpm(selected_impute_sce)
tpm_copy <- tpm_original
head(colnames(tpm_copy))
colnames(tpm_copy) <- paste0(colnames(tpm_copy), "_1")
head(colnames(tpm_copy))
all_data <- cbind(tpm_original, tpm_copy)

# mark the metabolic genes 
metabolics2 <- read.delim("metabolic_genes.csv", sep = ",")
metabolics <- as.character(metabolics2[,2])
row_data <- data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE


# Build scater object 

all_data <- data.matrix(all_data)
rm(selected_impute_sce)
selected_impute_sce <- SingleCellExperiment(
  assays = list(tpm=all_data,exprs=all_data),
  colData = new_coldata,
  rowData = row_data
)

table(colData(selected_impute_sce)$Identity)
selected_impute_sce <- selected_impute_sce[,colData(selected_impute_sce)$Identity %in% c('qISC', 'prolif_ISC', 'ISC_all', 'EB', 'dEC', 'EC_all', 'aEC', 'mEC', 'pEC', 'EE')]
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
all_cell_types <- as.vector(selected_impute_sce$Identity)
cell_types <- unique(all_cell_types)

#some genes occur in multiple pathways.
gene_pathway_number <- num_of_pathways(pathway_file,rownames(selected_impute_sce)[rowData(selected_impute_sce)$metabolic])

# Make duplicate of norm_tpm, and append '_1' and make new scater object, so new metadata can match

norm_tpm <- readRDS("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/Deconvolution_TPM/Deconvolution_tpm.rds")
norm_tpm_dup <- norm_tpm
colnames(norm_tpm_dup) <- paste0(colnames(norm_tpm_dup), "_1")
norm_tpm <- cbind(norm_tpm, norm_tpm_dup)
norm_tpm <- norm_tpm[,colnames(norm_tpm) %in% rownames(colData(selected_impute_sce))]


set.seed(123)
##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next
  
  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  
  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]
  
  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
  
  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
  
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:1000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
  names(shuffle_cell_types_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- cell_types
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 1000 
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 1000
    }
    if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}
all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]


######### Draw heatmap

dat <- mean_expression_shuffle
saveRDS(dat, file = file.path(outDir, 'heatmap_combined.rds'))
dat <- readRDS('~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/Heatmap_combined_heatmap/heatmap_combined.rds')
sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1

### Sort by pathway scores for qISC and reorder cell types
sort_row <- rownames(dat)[order(dat[,'qISC'],decreasing = T)]
sort_column <- c('qISC', 'ISC', 'ISC_all', 'EB', 'dEC', 'EC_all', 'aEC', 'mEC', 'pEC', 'EE')

mybreaks <- c(
  seq(0, 0.5, length.out=33),
  seq(0.51, 1.5, length.out=33),
  seq(1.51, max(dat),length.out=34)
) 
color <- colorRampPalette(c("#0072B2","white","#D55E00"))(100)
jpeg(file = file.path("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/Heatmap_combined_heatmap/heatmap_combined.jpeg"), width = 5000, height = 6000, res = 600)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,color=color,fontsize = 9, breaks = mybreaks,
         labels_row = make_bold_names(dat[sort_row,sort_column], rownames, c(rownames(dat[sort_row,sort_column]))),
         labels_col = make_bold_names(dat[sort_row,sort_column], colnames, c(colnames(dat[sort_row,sort_column]))))
dev.off()

svg("~/Desktop/Drosophila_midgut/metabolism/Dissection_Figures/More_Pathway_Heatmaps/Heatmap_combined_heatmap/heatmap_combined.svg", width = 8, height = 14)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,color=color,fontsize = 9, breaks = mybreaks,
         labels_row = make_bold_names(dat[sort_row,sort_column], rownames, c(rownames(dat[sort_row,sort_column]))),
         labels_col = make_bold_names(dat[sort_row,sort_column], colnames, c(colnames(dat[sort_row,sort_column]))))
dev.off()

write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue.txt"),row.names=T,col.names=T,quote=F,sep="\t")

#boxplot show the distribution of pathway activity
scRNA_dat <- as.data.frame(mean_expression_noshuffle)
scRNA_dat$X <- NULL

scRNA_df <- melt(scRNA_dat)
scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]
g <- ggplot(scRNA_df,aes(x=variable,y=value,fill=variable)) +
  scale_y_continuous(limits=c(0,4),breaks=0:4,labels=0:4)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=1,color="blue")+
  scale_fill_brewer(palette="Set2") +
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 6),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))

jpeg(file.path(outDir,"pathway_activity_violinplot_ISC_EB_eachEC.jpeg"),width = 1000, height = 1200, res = 600)
g
dev.off()
