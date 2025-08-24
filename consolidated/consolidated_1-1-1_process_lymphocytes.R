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

indir <- "RDS/"
outdir <- "RDS/"
if (!dir.exists(outdir)) dir.create(outdir)

whole <- readRDS(paste0(outdir,"whole_V1p2.RDS"))

## ---------------isolate (pan)-lymphocytes

lymphocytes <- subset(whole, subset = RNA_snn_res.0.2 %in% c(1,3,9))

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

saveRDS(lymphocytes, file = paste0(outdir,"lymphocytes_V1p2.RDS"))

## ------ remove doublet

myeloid_s <- subset(lymphocytes , subset= RNA_snn_res.0.2 %in% c(0) )
saveRDS(myeloid_s, file = paste0(outdir,"myeloid-s_V1p2.RDS"))



lymphocytes <- subset(lymphocytes, subset = RNA_snn_res.0.2 %in% c(1:8))

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
lymphocytes <- RunUMAP(lymphocytes,dims = 1:20,reduction = "harmony")
lymphocytes <- RunTSNE(lymphocytes,dims = 1:20,reduction = "harmony")

lymphocytes <- FindNeighbors(lymphocytes, dims = 1:20,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  lymphocytes <- FindClusters(object = lymphocytes, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
lymphocytes <- JoinLayers(lymphocytes)

saveRDS(lymphocytes, file = paste0(outdir,"lymphocytes_V1p3.RDS"))

## -------------prelim annotation



lymphocytes$compartment <- "Lymphocytes"

lymphocytes@meta.data <- lymphocytes@meta.data %>% mutate(
  annotation_V1 = case_match(
    as.character(RNA_snn_res.0.7),
    "0" ~ "L0: NK",
    "1" ~ "L1: Stem-like CD8 Tex",
    "2" ~ "L2: NKT",
    "3" ~ "L3: Treg/Th2",
    "4" ~ "L4: Tem/ex",
    "5" ~ "L5: Tcm/em",
    "6" ~ "L6: Th2",
    "7" ~ "L7: B cells",
    "8" ~ "L8: Th",
    "9" ~ "L9: gdT",
    "10" ~ "L10: ILC3",
    "11" ~ "L11: proliferating cells"
  )
)

saveRDS(lymphocytes, file = paste0(outdir,"lymphocytes_V1p3.RDS"))

