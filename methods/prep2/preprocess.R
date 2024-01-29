library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(Signac)

S = readRDS("../../data/ex2/kidney.rds")
S = S[,S@meta.data$cell_type %in% c('Podocyte','Fibroblast')]
S = NormalizeData(S)
S = FindVariableFeatures(S, selection.method = "vst", nfeatures = 6000)
write.table(VariableFeatures(S),"../../data/ex22/FS.csv",row.names=T,col.names=TRUE,sep=",")
# S = ScaleData(S)
# S = RunPCA(S)
# S = RunUMAP(S, dims = 1:10)
# colnames(S@meta.data)[11] = 'cell_type'
# DimPlot(S, reduction = 'umap', group.by = 'cell_type', raster = F)+ggtitle('S_umap')
# ggsave("../../data/ex1/S_umap.pdf", width = 7, height = 7, units = "in", dpi = 300)
lS = S@meta.data$cell_type
write.table(lS,"../../data/ex22/lS.csv",row.names=T,col.names=TRUE,sep=",")
Ta = readRDS("../../data/ex2/T.rds")
U = readRDS("../../data/ex2/U.rds")
Ta = Ta[,Ta@meta.data$cell_type %in% c('Connecting Tubule Cell','Podocyte','Ascending Vasa Recta Endothelial Cell','Fibroblast')]
U = U[,U@meta.data$cell_type %in% c('Connecting Tubule Cell','Podocyte','Ascending Vasa Recta Endothelial Cell','Fibroblast')]
Ta = NormalizeData(Ta)
Ta = FindVariableFeatures(Ta, selection.method = "vst", nfeatures = 6000)
write.table(VariableFeatures(Ta),"../../data/ex22/FT.csv",row.names=T,col.names=TRUE,sep=",")
# TU = ScaleData(TU)
# TU = RunPCA(TU)
# TU = RunUMAP(TU, dims = 1:10)
# DimPlot(TU, reduction = 'umap', group.by = 'cell_type')+ggtitle('T_umap')
# ggsave("../../data/ex1/T_umap.pdf", width = 7, height = 7, units = "in",dpi = 300)
lT = Ta@meta.data$cell_type
write.table(lT,"../../data/ex22/lT.csv",row.names=T,col.names=TRUE,sep=",")
SF = read.csv("../../data/ex22/FS.csv")
SF = SF[[1]]
TF = read.csv("../../data/ex22/FT.csv")
TF = TF[[1]]
sum(SF %in% TF)
Feature = SF[SF %in% TF]
write.table(Feature,"../../data/ex22/F.csv",row.names=T,col.names=TRUE,sep=",")
mS = S@assays$RNA@data[Feature,]
write.csv(mS,"../../data/ex22/mS.csv")
mT = Ta@assays$RNA@data[Feature,]
write.csv(mT,"../../data/ex22/mT.csv")
mU = U@assays$peaks@counts
ct = rowSums(mU)
UF = names(sort(ct, decreasing = T)[1:2000])
mU = U@assays$peaks@counts[UF,]
Matrix::writeMM(mU,"../../data/ex22/mU.csv")






