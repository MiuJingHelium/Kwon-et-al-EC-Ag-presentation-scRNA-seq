library(tidyverse)
library(Seurat)
library(scater)

library(RColorBrewer)
getPalette.1 <- colorRampPalette(brewer.pal(n = 12,name = "Set3"))
pal <- c(brewer.pal(n = 12,name = "Set3"),brewer.pal(n = 12,name = "Paired"),brewer.pal(n = 8,name = "Accent"))
pal[2] <- "darkred"
pal[9] <- "grey25"
pal[23] <- "turquoise4"
pal[19] <- colors()[60]
pal[13] <- "darkblue"
pal[28] <- colors()[16]
pal[31] <- colors()[54] 
pal <- c(pal,colors()[c(85,84,83,89,45,31,41,42,26,29,10,139,107,108,120,109,119,121,143)])

RDSdir <- "RDS/"

lymphocytes <- readRDS(file = paste0(RDSdir,"lymphocytes_V1p3.RDS"))
T_cells <- readRDS(paste0(RDSdir,"T_cells_V1p4.Rds"))

T_cells <- readRDS(paste0(RDSdir,"T_cells_V1p4.Rds"))
lymphocytes <- readRDS(file = paste0(RDSdir,"lymphocytes_V1p3.RDS"))

non_T <- subset(lymphocytes, subset = annotation_V1 %in% c("L0: NK","L7: B cells","L10: ILC3"))
lymphocytes <- merge(T_cells, non_T)
lymphocytes <- JoinLayers(lymphocytes)

lymphocytes <- FindVariableFeatures(object = lymphocytes, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
VariableFeatures(lymphocytes) <-  VariableFeatures(lymphocytes)[!grepl("^Tra|^Trb|^Igh|^Igk|Mamu-a", VariableFeatures(lymphocytes))]
VariableFeaturePlot(lymphocytes)
lymphocytes <- ScaleData(object = lymphocytes, features = VariableFeatures(object = lymphocytes), vars.to.regress = c("nCount_RNA", "percent.mt"))
lymphocytes <- RunPCA(object = lymphocytes,
                      features =  VariableFeatures(object = lymphocytes),
                      dims = 1:50)
gc()
ElbowPlot(lymphocytes,ndims = 50)

lymphocytes <- RunHarmony(object = lymphocytes, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
lymphocytes <- RunUMAP(lymphocytes,dims = 1:30,reduction = "harmony")
lymphocytes <- RunTSNE(lymphocytes,dims = 1:30,reduction = "harmony")

lymphocytes <- FindNeighbors(lymphocytes, dims = 1:30,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  lymphocytes <- FindClusters(object = lymphocytes, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
lymphocytes <- JoinLayers(lymphocytes)

saveRDS(lymphocytes, file = paste0(RDSdir,"lymphocytes_V1p4.RDS"))
