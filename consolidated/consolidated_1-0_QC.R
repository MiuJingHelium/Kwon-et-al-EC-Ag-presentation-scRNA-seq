library(tidyverse)
library(Seurat)
library(scater)
library(DoubletFinder)
library(miQC)
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

indir <- "SoupX/"
outdir <- "RDS/"
if (!dir.exists(outdir)) dir.create(outdir)

# ------------------------------------------------------------------------------

# code credit: 
# http://github.com/gatelabNW/csf_aging/blob/main/code/0_preprocessing/04_preprocessing_doubletfinder_seurat.R

run_doubletfinder <- function(s) {
  ## Process normally
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
  
  # Run TSNE clustering
  s <- RunPCA(s, features = VariableFeatures(object = s))
  s <- FindNeighbors(s, dims = 1:12)
  s <- FindClusters(s, resolution = 0.3)
  s <- RunTSNE(s, dims = 1:12)
  
  ## pK Identification (no ground-truth) ---------------------------------------
  sweep.res.list <- paramSweep(s, PCs = 1:12, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  print(s)
  
  ## Run DoubletFinder with varying classification stringencies ----------------
  s <- doubletFinder(s, PCs = 1:12, pN = 0.25, pK = optimal_pK,
                     nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  ## Rename column name for consistency
  colnames(s@meta.data)[ grep("DF.classifications*",
                              colnames(s@meta.data)) ] <- "DF.classifications"
  print(table(s[["DF.classifications"]]))
  
  return(s)
}

# ------------------------------------------------------------------------------


samples <- list.files(indir)

seu_list <- lapply(samples, function(x){
  seu_temp <- Read10X(paste0(indir,x,"/")) # a temporary object
  seu <- CreateSeuratObject(seu_temp)
  seu$orig.ident <- x
  seu
})
names(seu_list) <- samples # name the elements by the sample name

seu_list <- lapply(seu_list, function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(object = x, pattern = "^mt|^Mt")
  x$log10_nCount_RNA <- log10(x$nCount_RNA+1) # add a pseudcount offset to prevent log10(0)
  x$log10_nFeature_RNA <- log10(x$nFeature_RNA+1)
  x[["high_mito"]] <- isOutlier(x$percent.mt, type="higher", min.diff=0.5)
  
  x$percent.ribo <- PercentageFeatureSet(x,pattern="^Rp[ls]")
  x@meta.data <- x@meta.data %>% mutate(
    mitoQC = case_when(
      high_mito ~ "remove",
      .default = "keep"
    )
  )
  
  x
})

doublet_formation_rate <- 0.054
seu_list <- lapply(seu_list, run_doubletfinder)

seu_list <- lapply(seu_list, function(x){
  RunMiQC(x, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75,
          model.slot = "flexmix_model")
})

whole <- Reduce(function(x, y) merge(x, y), seu_list)

whole$percent.ribo <- PercentageFeatureSet(whole,pattern="^Rp[ls]")

whole@meta.data[, grep("pANN",colnames(whole@meta.data))] <- NULL

##################### sanity check plots #######################################
VlnPlot(whole,features = c("log10_nFeature_RNA","log10_nCount_RNA"),group.by = "orig.ident",alpha = 0.1,
        split.by = "miQC.keep") #+ scale_fill_manual(values = c("Singlet" = "darkblue","Doublet" = "red3"))#the sample quality are very similar
VlnPlot(whole,features = c("log10_nFeature_RNA","log10_nCount_RNA"),group.by = "orig.ident",alpha = 0.1,
        split.by = "DF.classifications")

ggplot(whole@meta.data,aes(x = log10_nCount_RNA, y = percent.mt, color = miQC.keep))+
  facet_wrap(~orig.ident)+
  geom_point()+
  scale_color_manual(values = c("keep" = "darkblue","discard" = "red3"))+
  theme_classic()
ggplot(whole@meta.data,aes(x = log10_nCount_RNA, y = percent.mt, color = mitoQC))+
  facet_wrap(~orig.ident)+
  geom_point()+
  scale_color_manual(values = c("keep" = "darkblue","remove" = "red3"))+
  theme_classic()

VlnPlot(whole,features = c("percent.mt","percent.ribo"),group.by = "orig.ident",alpha = 0.1) #the sample quality are very similar

ggplot(whole@meta.data,aes(x = log10_nCount_RNA, fill = miQC.keep))+
  facet_wrap(~orig.ident+miQC.keep, scales = "free")+
  geom_histogram()+
  scale_fill_manual(values = c("keep" = "darkblue","discard" = "red3"))+
  theme_classic()
ggplot(whole@meta.data,aes(x = log10_nFeature_RNA, fill = miQC.keep))+
  facet_wrap(~orig.ident+miQC.keep, scales = "free")+
  geom_histogram()+
  scale_fill_manual(values = c("keep" = "darkblue","discard" = "red3"))+
  theme_classic()
ggplot(whole@meta.data,aes(x = percent.mt, fill = miQC.keep))+
  facet_wrap(~orig.ident+miQC.keep, scales = "free")+
  scale_fill_manual(values = c("keep" = "darkblue","discard" = "red3"))+
  geom_histogram()+
  theme_classic()

ggplot(whole@meta.data,aes(x = log10_nCount_RNA, fill = mitoQC))+
  facet_wrap(~orig.ident+mitoQC, scales = "free")+
  geom_histogram()+
  scale_fill_manual(values = c("keep" = "darkblue","remove" = "red3"))+
  theme_classic()
ggplot(whole@meta.data,aes(x = log10_nFeature_RNA, fill = mitoQC))+
  facet_wrap(~orig.ident+mitoQC, scales = "free")+
  geom_histogram()+
  scale_fill_manual(values = c("keep" = "darkblue","remove" = "red3"))+
  theme_classic()
ggplot(whole@meta.data,aes(x = percent.mt, fill = mitoQC))+
  facet_wrap(~orig.ident+mitoQC, scales = "free")+
  scale_fill_manual(values = c("keep" = "darkblue","remove" = "red3"))+
  geom_histogram()+
  theme_classic()

################################################################################


saveRDS(whole, file = paste0(outdir,"CR8_whole_prefilter.RDS"))

table(whole$DF.classifications, whole$miQC.keep)

### Round 1
whole <- subset(whole, subset =  (miQC.keep == "keep"))

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

saveRDS(whole, file = paste0(outdir,"whole_V1.RDS"))

### Round 2

whole <- subset(whole, subset = DF.classifications == "Singlet")

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

saveRDS(whole, file = paste0(outdir,"whole_V1p1.RDS"))

### Round 3

# remove c10,12 at res=0.2
whole <- subset(whole, subset = RNA_snn_res.0.2 %in% c(0:9,11))

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

Idents(whole) <- "RNA_snn_res.0.2"

DimPlot(whole, group.by = "RNA_snn_res.0.2",label = T)+scale_color_manual(values = pal)
DimPlot(whole, group.by = "orig.ident",shuffle = T)+scale_color_manual(values = pal)
FeaturePlot(whole, features = "Ptprc",label = T)+scale_color_gradientn(colors = c("grey90","red2","red3"))

VlnPlot(whole,features = c("log10_nFeature_RNA"),alpha = 0.1) +scale_fill_manual(values = pal)
VlnPlot(whole,features = c("log10_nCount_RNA"),alpha = 0.1) +scale_fill_manual(values = pal)
VlnPlot(whole,features = c("percent.mt"),alpha = 0.1) +scale_fill_manual(values = pal)

DotPlot(whole,features = c("Ptprc","Cd3e","Cd8a","Cd4","Cd79a","Cd19","Ccr2","Adgre1","Lyz2","Klrb1c","Tyrobp","Ncam1","Xcr1","Ncr1","Siglech","Mpo","S100a8","Cd200r3","Fcer1a","Fcer2a","Flt3","Pecam1","Epcam","Krt1","Myct1","Hbb-bs","Pf4","Col1a1","Mki67","Jchain"),group.by = "RNA_snn_res.0.2")+RotatedAxis()+scale_color_gradientn(colors = c("grey90","red2","red3"))

DimPlot(whole,group.by = "DF.classifications")+scale_color_manual(values = pal[c(19,21)])

DotPlot(whole,features = c("Ptprc","Cd3e","Cd8a","Cd8b1","Cd4","Ccr7","Tcf7","Lef1","Foxp3","Klrb1c","Tyrobp","Prf1","Gzma","Gzmk","Lag3","Cd44","Pdcd1","Tox","Eomes","Ikzf2","Mki67","H2-Ab1","Cd74"),group.by = "RNA_snn_res.0.2")+RotatedAxis()+scale_color_gradientn(colors = c("grey90","red2","red3"))

saveRDS(whole, file = paste0(outdir,"whole_V1p2.RDS"))



