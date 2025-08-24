library(tidyverse)
library(Seurat)
library(harmony)

library(RColorBrewer)
getPalette.1 <- colorRampPalette(brewer.pal(n = 12,name = "Set3"))
pal <- c(brewer.pal(n = 12,name = "Set3"),brewer.pal(n = 12,name = "Paired"),brewer.pal(n = 8,name = "Accent"))
pal[2] <- "darkred"
pal[9] <- "grey25"
pal[23] <- "darkblue"
pal[19] <- colors()[60]
pal[13] <- "turquoise4"
pal[28] <- colors()[16]
pal[31] <- colors()[54] 
pal <- c(pal,colors()[c(85,84,83,89,45,31,41,42,26,29,10,139,107,108,120,109,119,121,143)])

pal2 <- pal[c(1:8,10,13:14,16:20, 22:25,27:30, 33:35 )] # use pal2


RDSdir <- "RDS/"


lymphocytes <- readRDS(paste0(RDSdir,"lymphocytes_V1p4.RDS"))
myeloid <- readRDS(paste0(RDSdir,"myeloid_V1p3-2.RDS"))

# in case of previously added metadata
lymphocytes$umap_1 <- NULL
lymphocytes$umap_2 <- NULL
lymphocytes$annotation_V2 <- NULL
lymphocytes$barcode <- NULL

whole <- merge(x= myeloid, y = lymphocytes)
whole <- JoinLayers(whole)

whole <- NormalizeData(object = whole, normalization.method = "LogNormalize", scale.factor = 10000)
whole <- FindVariableFeatures(object = whole, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
VariableFeatures(whole) <-  VariableFeatures(whole)[!grepl("^Tra|^Trb|^Igh|^Igk|Mamu-a", VariableFeatures(whole))]
VariableFeaturePlot(whole)
whole <- ScaleData(object = whole, features = VariableFeatures(object = whole), vars.to.regress = c("nCount_RNA", "percent.mt"))
whole <- RunPCA(object = whole,
                features =  VariableFeatures(object = whole),
                dims = 1:50)
gc()
ElbowPlot(whole,ndims = 50)

whole <- RunHarmony(object = whole, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
whole <- RunUMAP(whole,dims = 1:30,reduction = "harmony")
whole <- RunTSNE(whole,dims = 1:30,reduction = "harmony")

whole <- FindNeighbors(whole, dims = 1:30,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  whole <- FindClusters(object = whole, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
whole <- JoinLayers(whole)

saveRDS(whole, file = paste0(RDSdir, "whole_V1p4-2.RDS"))