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
outdir <- "T_cells_out/"
if (!dir.exists(outdir)) dir.create(outdir)

lymphocytes <- readRDS(file = paste0(indir,"lymphocytes_V1p3.RDS"))

### ------- isolate all Cd3e+ clusters

T_cells <- subset(lymphocytes, idents = c(1,2,3,4,5,6,8,9,11))

T_cells <- FindVariableFeatures(object = T_cells, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
VariableFeatures(T_cells) <-  VariableFeatures(T_cells)[!grepl("^Tra|^Trb|^Igh|^Igk|Mamu-a", VariableFeatures(T_cells))]
VariableFeaturePlot(T_cells)
T_cells <- ScaleData(object = T_cells, features = VariableFeatures(object = T_cells), vars.to.regress = c("nCount_RNA", "percent.mt"))
T_cells <- RunPCA(object = T_cells,
                  features =  VariableFeatures(object = T_cells),
                  dims = 1:50)
gc()
ElbowPlot(T_cells,ndims = 50)

T_cells <- RunHarmony(object = T_cells, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
T_cells <- RunUMAP(T_cells,dims = 1:20,reduction = "harmony")
T_cells <- RunTSNE(T_cells,dims = 1:20,reduction = "harmony")

T_cells <- FindNeighbors(T_cells, dims = 1:20,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  T_cells <- FindClusters(object = T_cells, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
T_cells <- JoinLayers(T_cells)

### ----- remove additional artifact

Idents(T_cells) <- "RNA_snn_res.0.6"
T_cells <- subset(T_cells,idents = 6, invert = T)

T_cells <- FindVariableFeatures(object = T_cells, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
VariableFeatures(T_cells) <-  VariableFeatures(T_cells)[!grepl("^Tra|^Trb|^Igh|^Igk|Mamu-a", VariableFeatures(T_cells))]
VariableFeaturePlot(T_cells)
T_cells <- ScaleData(object = T_cells, features = VariableFeatures(object = T_cells), vars.to.regress = c("nCount_RNA", "percent.mt"))
T_cells <- RunPCA(object = T_cells,
                  features =  VariableFeatures(object = T_cells),
                  dims = 1:50)
gc()
ElbowPlot(T_cells,ndims = 50)

T_cells <- RunHarmony(object = T_cells, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
T_cells <- RunUMAP(T_cells,dims = 1:20,reduction = "harmony")
T_cells <- RunTSNE(T_cells,dims = 1:20,reduction = "harmony")

T_cells <- FindNeighbors(T_cells, dims = 1:20,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  T_cells <- FindClusters(object = T_cells, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
T_cells <- JoinLayers(T_cells)

saveRDS(T_cells, file = paste0(indir,"T_cells_V1p4.Rds"))

T_cells@meta.data <- T_cells@meta.data %>% mutate(annotation_V1 = case_match(
  as.character(RNA_snn_res.0.6),
  "0" ~ "c0: Naive CD8 T",
  "1" ~ "c1: NKT-like",
  "2" ~ "c2: Treg",
  "3" ~ "c3: Th2",
  "4" ~ "c4: Tem/ex",
  "5" ~ "c5: Th",
  "6" ~ "c6: Tex",
  "7" ~ "c7: gdT",
  "8" ~ "c8: Proliferating"
))

saveRDS(T_cells, file = paste0(indir,"T_cells_V1p4.Rds"))

