# processed data and metadata were downloaded from GEO.

library(Seurat)
library(tidyverse)
library(harmony)
library(miQC)
library(SeuratWrappers)
library(scater)

indir <- "GSE284013_RAW/"
outdir <- "GSE284013_RDS/"
if(!dir.exists(outdir))dir.create(outdir)

## ------ set up object
files <- list.files(indir)
files <- files[grepl("day_5|day_45", files)]

day5 <- Read10X_h5(filename = paste0(indir,files[grep("day_5",files)]))
day45 <- Read10X_h5(filename = paste0(indir,files[grep("day_45",files)]))

seu_list <- list(day5, day45)
names(seu_list) <- c("day5","day45")

seu_list <- lapply(seu_list, function(x){
  x <- CreateSeuratObject(x)
  x[["percent.mt"]] <- PercentageFeatureSet(object = x, pattern = "^mt|^MT|^Mt")
  x$log10_nCount_RNA <- log10(x$nCount_RNA+1) # add a pseudcount offset to prevent log10(0)
  x$log10_nFeature_RNA <- log10(x$nFeature_RNA+1)
  x$percent.ribo <- PercentageFeatureSet(x,pattern="^Rp[ls]")
  x[["high_mito"]] <- isOutlier(x$percent.mt, type="higher", min.diff=0.5)
  
  x@meta.data <- x@meta.data %>% mutate(
    mitoQC = case_when(
      high_mito ~ "high_mito",
      .default = "normal"
    ))
  
  x <- RunMiQC(x, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.8,
               model.slot = "flexmix_model")
  x
})

seu_list$day5$orig.ident <- "Day 5"
seu_list$day45$orig.ident <- "Day 45"

whole <- merge(seu_list$day5, seu_list$day45)

whole <- subset(whole, subset = percent.mt <= 8 )

whole <- subset(whole, subset = nFeature_RNA <= 5000 & nFeature_RNA >= 1000)

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
whole <- RunUMAP(whole,dims = 1:20,reduction = "harmony")
whole <- RunTSNE(whole,dims = 1:20,reduction = "harmony")

whole <- FindNeighbors(whole, dims = 1:20,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  whole <- FindClusters(object = whole, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
whole <- JoinLayers(whole)

## ------ additional clean-up

whole <- subset(whole, subset = RNA_snn_res.0.6 != 10)

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
whole <- RunUMAP(whole,dims = 1:20,reduction = "harmony")
whole <- RunTSNE(whole,dims = 1:20,reduction = "harmony")

whole <- FindNeighbors(whole, dims = 1:20,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  whole <- FindClusters(object = whole, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
whole <- JoinLayers(whole)

whole <- subset(whole, subset = RNA_snn_res.0.6 != 9)

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
whole <- RunUMAP(whole,dims = 1:20,reduction = "harmony")
whole <- RunTSNE(whole,dims = 1:20,reduction = "harmony")

whole <- FindNeighbors(whole, dims = 1:20,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  whole <- FindClusters(object = whole, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
whole <- JoinLayers(whole)

whole <- subset(whole, subset = percent.mt <= 5)

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
whole <- RunUMAP(whole,dims = 1:15,reduction = "harmony")
whole <- RunTSNE(whole,dims = 1:15,reduction = "harmony")

whole <- FindNeighbors(whole, dims = 1:15,reduction = "harmony")
for(res in seq(0.1, 1, 0.1))  {   
  whole <- FindClusters(object = whole, resolution = res, print.output = 0, save.SNN = TRUE)
} 
gc()
whole <- JoinLayers(whole)

saveRDS(whole, file = paste0(outdir,"LCMV_day5-45_processed_V2.RDS"))

## -------------------extract signature

comparison <- subset(whole, subset = RNA_snn_res.0.4 %in% c(0,1,2,5))

comparison@meta.data <- comparison@meta.data %>% mutate(
  annotation = case_match(
    as.character(RNA_snn_res.0.4),
    "0" ~ "late TD",
    "1" ~ "early stem-like",
    "2" ~ "early stem-like",
    "5" ~ "late TD"
  )
)

early_stem <- subset(comparison, subset = annotation == "early stem-like" & orig.ident == "Day 5")
late_TD <- subset(comparison, subset = annotation == "late TD" & orig.ident == "Day 45")

comparison <- merge(early_stem,late_TD)
comparison <- JoinLayers(comparison) 
# no need to go through workflow since DEG are calculated using data slot aka the normalized reads.

Idents(comparison) <- "annotation"
late_TD_DEG <- FindMarkers(comparison, ident.1 = "late TD", ident.2 = "early stem-like",test.use = "MAST", min.pct = 0.2)
late_TD_DEG$gene <- rownames(late_TD_DEG)
late_TD_sig <- unlist(late_TD_DEG %>% slice_max(order_by = avg_log2FC, n = 150) %>% select(gene))

write.table(late_TD_sig,file  = "signature_mapping_LCMV/late_TD_vs_early_stem_signature.txt",sep = "\t",col.names = F, row.names = F, quote = F)

early_stem_DEG <- FindMarkers(comparison, ident.2 = "late TD", ident.1 = "early stem-like",test.use = "MAST", min.pct = 0.2)


early_stem_DEG$gene <- rownames(early_stem_DEG)

early_stem_sig <- unlist(early_stem_DEG %>% slice_max(order_by = avg_log2FC, n = 150) %>% select(gene))

write.table(early_stem_sig,file  = "signature_mapping_LCMV/early_stem_vs_late_TD_signature.txt",sep = "\t",col.names = F, row.names = F, quote = F)






