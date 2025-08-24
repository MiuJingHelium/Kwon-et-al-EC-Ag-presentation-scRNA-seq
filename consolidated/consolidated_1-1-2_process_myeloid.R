library(tidyverse)
library(Seurat)
library(scater)
library(harmony)
library(SeuratWrappers)

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

indir <- "RDS/"
outdir <- "RDS/"
if (!dir.exists(outdir)) dir.create(outdir)

whole <- readRDS(paste0(outdir,"whole_V1p2.RDS"))

## ---merge myeloid with doublet-like cluster

myeloid <- subset(whole, subset = RNA_snn_res.0.2 %in% c(0, 4, 5, 6, 7, 10))
myeloid_s <- readRDS(paste0(outdir,"myeloid-s_V1p2.RDS"))

myeloid <- merge(myeloid, myeloid_s)
myeloid <- JoinLayers(myeloid)

myeloid <- FindVariableFeatures(object = myeloid, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
VariableFeaturePlot(myeloid)
myeloid <- ScaleData(object = myeloid, features = VariableFeatures(object = myeloid), vars.to.regress = c("nCount_RNA", "percent.mt"))
myeloid <- RunPCA(object = myeloid,
                  features =  VariableFeatures(object = myeloid),
                  dims = 1:50)
gc()
ElbowPlot(myeloid,ndims = 50)

myeloid <- RunHarmony(object = myeloid, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
myeloid <- RunUMAP(myeloid,dims = 1:30,reduction = "harmony")
myeloid <- RunTSNE(myeloid,dims = 1:30,reduction = "harmony")

myeloid <- FindNeighbors(myeloid, dims = 1:30,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  myeloid <- FindClusters(object = myeloid, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
myeloid <- JoinLayers(myeloid)

saveRDS(myeloid, file = paste0(outdir,"myeloid_V1p2.RDS"))

## ---second iter without c2

Idents(myeloid) <- "RNA_snn_res.0.2"
myeloid <- subset(myeloid, idents=c(0:1,3:11))

myeloid <- FindVariableFeatures(object = myeloid, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
VariableFeaturePlot(myeloid)
myeloid <- ScaleData(object = myeloid, features = VariableFeatures(object = myeloid), vars.to.regress = c("nCount_RNA", "percent.mt"))
myeloid <- RunPCA(object = myeloid,
                  features =  VariableFeatures(object = myeloid),
                  dims = 1:50)
gc()
ElbowPlot(myeloid,ndims = 50)

myeloid <- RunHarmony(object = myeloid, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
myeloid <- RunUMAP(myeloid,dims = 1:30,reduction = "harmony")
myeloid <- RunTSNE(myeloid,dims = 1:30,reduction = "harmony")

myeloid <- FindNeighbors(myeloid, dims = 1:30,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  myeloid <- FindClusters(object = myeloid, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
myeloid <- JoinLayers(myeloid)

myeloid$compartment <- "Myeloid"

myeloid@meta.data <- myeloid@meta.data %>% mutate(
  annotation_V1 = case_match(
    as.character(RNA_snn_res.0.2),
    "0" ~ "M0: mo-Mac1",
    "1" ~ "M1: mo-Mac2",
    "2" ~ "M2: Mac-like/artifact-like",
    "3" ~ "M3: mo-DC/monocytes",
    "4" ~ "M4: Neutrophil-like",
    "5" ~ "M5: migDC",
    "6" ~ "M6: proliferating mo-Mac",
    "7" ~ "M7: cDC1/DC2",
    "8" ~ "M8: mo-Mac3",
    "9" ~ "M9: mo-Mac4",
    "10" ~ "M10: granulocytes"
  )
)

saveRDS(myeloid, file = paste0(outdir,"myeloid_V1p3.RDS"))


## originally the M2-artifact like cluster was removed. It was later determined to be basophil and merged back. 
## The V1p3 myeloid was used for the final analysis. Alternatively, there is one more cluster removed from V1.3 myeloid

## ---- further update:

myeloid <- subset(myeloid,subset = (!RNA_snn_res.0.6 %in% c(7) ))

myeloid <- FindVariableFeatures(object = myeloid, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
VariableFeaturePlot(myeloid)
myeloid <- ScaleData(object = myeloid, features = VariableFeatures(object = myeloid), vars.to.regress = c("nCount_RNA", "percent.mt"))
myeloid <- RunPCA(object = myeloid,
                  features =  VariableFeatures(object = myeloid),
                  dims = 1:50)
gc()
ElbowPlot(myeloid,ndims = 50)

myeloid <- RunHarmony(object = myeloid, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
myeloid <- RunUMAP(myeloid,dims = 1:20,reduction = "harmony")
myeloid <- RunTSNE(myeloid,dims = 1:20,reduction = "harmony")

myeloid <- FindNeighbors(myeloid, dims = 1:20,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  myeloid <- FindClusters(object = myeloid, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
myeloid <- JoinLayers(myeloid)

myeloid@meta.data <- myeloid@meta.data %>% mutate(
  annotation_V1 = case_match(
    as.character(RNA_snn_res.0.6),
    "0" ~ "M0: Mo-Mac1",
    "1" ~ "M1: Mo-Mac2",
    "2" ~ "M2: Mo-Mac3",
    "3" ~ "M3: Basophil/granulocytes",
    "4" ~ "M4: Mo-Mac4",
    "5" ~ "M5: Mo-DC/monocytes",
    "6" ~ "M6: Neutrophil",
    "7" ~ "M7: Mo-Mac5",
    "8" ~ "M8: Proliferating Mo-Mac",
    "9" ~ "M9: migDC",
    "10" ~ "M10: cDC1/DC2",
    "11" ~ "M11: Mo-Mac6",
    "12" ~ "M12: Mo-Mac7",
    "13" ~ "M13: granulocytes/mast cells"
  )
)

saveRDS(myeloid, file = paste0(RDSdir,"myeloid_V1p3-2.RDS"))



