library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(Signac)

S = readRDS("../../data/ex1/S.rds")
colnames(S@meta.data)[11] = 'cell_type'
S = S[,!(S@meta.data$cell_type %in% c('Acinar'))]
S = NormalizeData(S)
S = FindVariableFeatures(S, selection.method = "vst", nfeatures = 6000)
write.table(VariableFeatures(S),"../../data/ex19/FS.csv",row.names=T,col.names=TRUE,sep=",")
# S = ScaleData(S)
# S = RunPCA(S)
# S = RunUMAP(S, dims = 1:10)
# colnames(S@meta.data)[11] = 'cell_type'
# DimPlot(S, reduction = 'umap', group.by = 'cell_type', raster = F)+ggtitle('S_umap')
# ggsave("../../data/ex1/S_umap.pdf", width = 7, height = 7, units = "in", dpi = 300)
lS = S@meta.data$cell_type
write.table(lS,"../../data/ex19/lS.csv",row.names=T,col.names=TRUE,sep=",")
TU = readRDS("../../data/ex1/TU.rds")
TU = TU[,!(TU@meta.data$cell_type %in% c('Acinar'))]
TU = NormalizeData(TU)
TU = FindVariableFeatures(TU, selection.method = "vst", nfeatures = 6000)
write.table(VariableFeatures(TU),"../../data/ex19/FT.csv",row.names=T,col.names=TRUE,sep=",")
# TU = ScaleData(TU)
# TU = RunPCA(TU)
# TU = RunUMAP(TU, dims = 1:10)
# DimPlot(TU, reduction = 'umap', group.by = 'cell_type')+ggtitle('T_umap')
# ggsave("../../data/ex1/T_umap.pdf", width = 7, height = 7, units = "in",dpi = 300)
lT = TU@meta.data$cell_type
write.table(lT,"../../data/ex19/lT.csv",row.names=T,col.names=TRUE,sep=",")
SF = read.csv("../../data/ex19/FS.csv")
SF = SF[[1]]
TF = read.csv("../../data/ex19/FT.csv")
TF = TF[[1]]
sum(SF %in% TF)
Feature = SF[SF %in% TF]
write.table(Feature,"../../data/ex19/F.csv",row.names=T,col.names=TRUE,sep=",")
mS = S@assays$RNA@data[Feature,]
write.csv(mS,"../../data/ex19/mS.csv")
mT = TU@assays$RNA@data[Feature,]
write.csv(mT,"../../data/ex19/mT.csv")
mU = TU@assays$ATAC@counts
ct = rowSums(mU)
UF = names(sort(ct, decreasing = T)[1:2000])
mU = TU@assays$ATAC@counts[UF,]
Matrix::writeMM(mU,"../../data/ex19/mU.csv")






