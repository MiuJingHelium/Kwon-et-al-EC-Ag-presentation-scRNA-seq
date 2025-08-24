library(tidyverse)
library(Seurat)
library(scater)
library(mascarade)
library(data.table)
library(ComplexHeatmap)
library(circlize)

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
pal2 <- pal[c(1:8,10,13:14,16:30)] # use pal2

RDSdir <- "RDS/"

outdir <- "V1p4_lymphocyte_plots/"
if (!dir.exists(outdir)) dir.create(outdir)

T_cells <- readRDS(paste0(RDSdir,"T_cells_V1p4.Rds"))
lymphocytes <- readRDS(file = paste0(RDSdir,"lymphocytes_V1p4.RDS"))
whole <- readRDS(file = paste0(RDSdir,"whole_V1p4-2.RDS"))

### UMAP
whole@meta.data <- whole@meta.data %>%
  mutate(highlight = case_match(
    compartment,
    "Lymphocytes" ~ "Lymphocytes",
    .default = "NA"
  ))

DimPlot(whole,group.by = "highlight") + 
  scale_color_manual(values = c("Lymphocytes" = colorRampPalette(brewer.pal(11, "RdYlBu"))(8)[3],na.value = "grey90"))
g <- DimPlot(whole,group.by = "highlight") + 
  scale_color_manual(values = c("Lymphocytes" = colorRampPalette(brewer.pal(11, "RdYlBu"))(8)[3],na.value = "grey90"))
ggsave(filename = paste0(outdir,"Supp_Whole_Lymphocytes_unsplit.pdf"), plot = g, units = "cm",height = 10,width = 12)


### FeaturePlot
# lymphocytes@meta.data <- lymphocytes@meta.data %>% select(-umap_1, -umap_2)

maskTable <- generateMask(dims=lymphocytes@reductions$umap@cell.embeddings, 
                          clusters=lymphocytes$annotation_V1)
data <- data.table(lymphocytes@reductions$umap@cell.embeddings, 
                   cluster=lymphocytes$annotation_V1,
                   lymphocytes@meta.data)
centers <- data %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::summarize(umap_1 = median(umap_1), umap_2 = median(umap_2)) 


data$orig.ident <- factor(data$orig.ident, levels = c("WT_HC","KO_HC"))

g <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=cluster), size = 0.5) + 
  geom_path(data=maskTable, aes(group=group)) +
  #geom_text(data = centers, aes(label=cluster), size = 3.5, color = "black") +
  scale_color_manual(values = pal2)+
  coord_fixed() + 
  theme_classic()+
  facet_wrap(~orig.ident)+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
g
ggsave(filename = paste0(outdir,"lymphocytes_V1p4_ann-V1_umap_split_mascarade.pdf"), plot = g, width = 18, height = 7)

g <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=cluster), size = 0.5) + 
  geom_path(data=maskTable, aes(group=cluster)) +
  ggrepel::geom_text_repel(data = centers, aes(label=cluster), size = 5, color = "black",min.segment.length = unit(0, 'lines'),force = 5, label.padding = 5,max.overlaps = Inf,direction = "both") +
  scale_color_manual(values = pal2)+
  coord_fixed() + 
  theme_classic()+theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
g
ggsave(filename = paste0(outdir,"lymphocytes_V1p4_ann-V1_umap_mascarade.pdf"), plot = g, width = 18, height = 7)

p1 <- DotPlot(lymphocytes, group.by = "annotation_V1",features = c("Cd19","Kit","Gata3","Ncr1","Tyrobp","Cd3e","Rorc","Klrb1c","Cd4","Foxp3","Il4","Il5","Il13","Cd8a","Lef1","Ccr7","Slamf6","Ccr4","Cd69","Il2ra","Cd44","Prf1","Gzmb","Gzmk","Lag3","Pdcd1","Tox","Havcr2","Mki67"))+scale_color_gradientn(colors = c("grey90","red2","red3"))+ RotatedAxis()+coord_flip()
ggsave(filename = paste0(outdir,"Supp_lymphocytes_Dotplot_annotation.pdf"), plot = p1, units = "cm",height = 24,width = 20)

dat <- GetAssayData(lymphocytes, slot = "data")

markers <- c("Cd19","Kit","Gata3","Ncr1","Tyrobp","Cd3e","Rorc","Klrb1c","Cd4","Foxp3","Il4","Il5","Il13","Cd8a","Lef1","Ccr7","Slamf6","Ccr4","Cd69","Il2ra","Cd44","Prf1","Gzmb","Gzmk","Lag3","Pdcd1","Tox","Havcr2","Mki67")

marker_exp <- as.data.frame(t(as.matrix(dat[rownames(dat) %in% markers,])))

marker_exp$barcode <- rownames(marker_exp)

data <- left_join(data, marker_exp, by = "barcode")

# underlay cluster color
g <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_polygon(data=maskTable, aes(group=group,fill = cluster), alpha = 1) +
  geom_point(aes(color=Foxp3), size = 0.5) + 
  #geom_text(data = centers, aes(label=cluster), size = 3.5, color = "black") +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  scale_fill_manual(values = pal2)+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Foxp3")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
g

g1 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Foxp3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group), color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c2: Treg", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  scale_fill_manual(values = pal2)+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Foxp3")+
  theme(legend.position = "right")
g1


g2 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd19), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "L7: B cells", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd19")+
  theme(legend.position = "right")
g2

g3 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Mki67), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c8: Proliferating", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Mki67")+
  theme(legend.position = "right")
g3

g4 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd3e), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=cluster), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd3e")+
  theme(legend.position = "right")
g4

g5 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd4), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c2: Treg","c3: Th2","c5: Th","c6: Tex"), cluster, "")), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd4")+
  theme(legend.position = "right")
g5

g6 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd8a), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c0: Naive CD8 T","c4: Tem/ex"), cluster, "")), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd8a")+
  theme(legend.position = "right")
g6

g7 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Rorc), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c7: gdT", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Rorc")+
  theme(legend.position = "right")
g7

g8 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Ncr1), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c1: NKT-like","L0: NK"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Ncr1")+
  theme(legend.position = "right")
g8

g9 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Kit), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "L10: ILC2", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Kit")+
  theme(legend.position = "right")
g9

g10 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Gata3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("L10: ILC2","c3: Th2","c5: Th"), cluster, "")), size = 3.5, color = "black", nudge_x = 2, nudge_y = -1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Gata3")+
  theme(legend.position = "right")
g10


g11 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Tyrobp), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c1: NKT-like","L0: NK"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Tyrobp")+
  theme(legend.position = "right")
g11

g12 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Il13), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c3: Th2"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Il13")+
  theme(legend.position = "right")
g12

g <- gridExtra::grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6, g7, g8,g11, g9, g10, g12), ncol = 4)
ggsave(paste0(outdir, "lymphocyte_marker_mascarade_ann-V1_featureplot.pdf"),plot = g, width = 15, height = 12 )

### ann-V1 proportion plots

sample_imm_size = data.frame(table(whole$orig.ident))
sample_imm_size <- sample_imm_size %>% dplyr::rename(n_immune = Freq, orig.ident = Var1)

clust_prop <- lymphocytes@meta.data %>% group_by(annotation_V1,compartment,orig.ident) %>% summarise(n_subpop = n()) %>% ungroup()
clust_prop <- left_join(clust_prop,sample_imm_size, by = "orig.ident")
clust_prop <- clust_prop %>% group_by(orig.ident) %>% mutate(percentage = 100*n_subpop/n_immune)

ratio <- clust_prop %>% group_by(annotation_V1) %>% reframe(FC = log2(percentage[orig.ident == "KO_HC"]/percentage[orig.ident == "WT_HC"]))

write.table(clust_prop, file = paste0(outdir,"lymphocyte_ann-V1_proportion.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(ratio, file = paste0(outdir,"lymphocyte_ann-V1_log2-proportion.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

g3 <- ggplot(clust_prop,aes(x = factor(orig.ident, levels = c("WT_HC", "KO_HC")), y = percentage, fill = annotation_V1))+
  geom_bar(aes(color = orig.ident), fill = "white",stat = "identity", position = "dodge")+
  scale_color_manual("condition",values = c("KO_HC" = "darkblue", "WT_HC" = "red3"))+theme_classic()+
  ylab("% in immune cells")+
  xlab("conditions")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  facet_wrap(~annotation_V1, scales = "free")
g3  
ggsave(filename = paste0(outdir,"lymphocyte_percentage_bar_fill_by_condition_ann-V1.pdf"), plot = g3, width = 8, height = 5)



g4 <- ggplot(ratio, aes( x = reorder(annotation_V1, -FC), y = FC))+
  geom_bar(aes(fill = annotation_V1),stat = "identity", position = "dodge")+
  scale_fill_manual(values = pal2)+theme_classic()+
  ylab("log2 % in KO v.s. WT immune cells")+
  xlab("subpopulations")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  geom_hline(yintercept = c(-1,1), color = "red3", linetype = 3)
g4
ggsave(filename = paste0(outdir,"lymphocyte_log2_percentage_bar_ann-V1.pdf"), plot = g4, width = 6, height = 4)

### T cell plot

outdir <- "V1p4_T-cell_plots/"
if (!dir.exists(outdir)) dir.create(outdir)

markers <- c("Tbx21","Gata3","Bcl6","Tyrobp","Cd3e","Rorc","Klrb1c","Cd4","Foxp3","Il4","Il5","Il13","Cd8a","Tcf7","Lef1","Ccr7","Id2","Id3","Prdm1","Slamf6","Ccr4","Cd69","Il2ra","Cd44","Prf1","Gzmb","Gzmk","Lag3","Pdcd1","Tox","Havcr2","Mki67","Itgal")

glist <- lapply(markers, function(x){
  FeaturePlot(T_cells, features = x,label = F)+scale_color_gradientn(colors = c("grey90","red2","red3"))
})
g <- gridExtra::grid.arrange(grobs = glist, ncol = 5)
ggsave(paste0(outdir, "T_markers_featurePlot.pdf"), plot = g, width = 30, height = 30)

#T_cells@meta.data <- T_cells@meta.data %>% select(-umap_1,-umap_2)

maskTable <- generateMask(dims=T_cells@reductions$umap@cell.embeddings, 
                          clusters=T_cells$annotation_V1)
data <- data.table(T_cells@reductions$umap@cell.embeddings, 
                   cluster=T_cells$annotation_V1,
                   T_cells@meta.data)
centers <- data %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::summarize(umap_1 = median(umap_1), umap_2 = median(umap_2)) 


data$orig.ident <- factor(data$orig.ident, levels = c("WT_HC","KO_HC"))

g <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=cluster), size = 0.5) + 
  geom_path(data=maskTable, aes(group=group)) +
  #geom_text(data = centers, aes(label=cluster), size = 3.5, color = "black") +
  ggrepel::geom_text_repel(data = centers, aes(label=cluster), size = 5, color = "black",min.segment.length = unit(0, 'lines'),force = 5, label.padding = 5,max.overlaps = Inf,direction = "both") +
  scale_color_manual(values = pal2)+
  coord_fixed() + 
  theme_classic()+
  facet_wrap(~orig.ident)+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
g
ggsave(filename = paste0(outdir,"T_V1p4_ann-V1_umap_split_mascarade.pdf"), plot = g, width = 18, height = 7)

g <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=cluster), size = 0.5) + 
  geom_path(data=maskTable, aes(group=cluster)) +
  ggrepel::geom_text_repel(data = centers, aes(label=cluster), size = 5, color = "black",min.segment.length = unit(0, 'lines'),force = 5, label.padding = 5,max.overlaps = Inf,direction = "both") +
  scale_color_manual(values = pal2)+
  coord_fixed() + 
  theme_classic()+theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
g
ggsave(filename = paste0(outdir,"T_V1p4_ann-V1_umap_mascarade.pdf"), plot = g, width = 18, height = 7)

markers <- c("Tbx21","Gata3","Bcl6","Tyrobp","Cd3e","Rorc","Klrb1c","Cd4","Foxp3","Il4","Il5","Il13","Cd8a","Tcf7","Lef1","Ccr7","Id2","Id3","Prdm1","Slamf6","Ccr4","Cd69","Il2ra","Cd44","Prf1","Gzmb","Gzmk","Lag3","Pdcd1","Tox","Havcr2","Mki67","Itgal")

dat <- GetAssayData(T_cells, slot = "data")

marker_exp <- as.data.frame(t(as.matrix(dat[rownames(dat) %in% markers,])))

marker_exp$barcode <- rownames(marker_exp)

data <- left_join(data, marker_exp, by = "barcode")

g1 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Foxp3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group), color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c2: Treg", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  scale_fill_manual(values = pal2)+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Foxp3")+
  theme(legend.position = "right")
g1


g2 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Klrb1c), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c1: NKT-like", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Klrb1c")+
  theme(legend.position = "right")
g2

g13 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Tyrobp), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c1: NKT-like", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Tyrobp")+
  theme(legend.position = "right")
g13

g3 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Mki67), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c8: Proliferating", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Mki67")+
  theme(legend.position = "right")
g3

```

```{r}
g4 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd3e), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=cluster), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd3e")+
  theme(legend.position = "right")
g4

g5 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd4), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c2: Treg","c3: Th2","c5: Th","c6: Tex"), cluster, "")), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd4")+
  theme(legend.position = "right")
g5

g6 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd8a), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c0: Naive CD8 T","c4: Tem/ex"), cluster, "")), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd8a")+
  theme(legend.position = "right")
g6

g7 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Rorc), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c7: gdT", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Rorc")+
  theme(legend.position = "right")
g7

g8 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Gzmb), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c4: Tem/ex","c2: Treg"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Gzmb")+
  theme(legend.position = "right")
g8

g19 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Gzmk), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c4: Tem/ex"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Gzmk")+
  theme(legend.position = "right")
g19

g9 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Slamf6), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster == "c0: Naive CD8 T", cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Slamf6")+
  theme(legend.position = "right")
g9


g10 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Gata3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("L10: ILC2","c3: Th2","c5: Th"), cluster, "")), size = 3.5, color = "black", nudge_x = 2, nudge_y = -1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Gata3")+
  theme(legend.position = "right")
g10


g11 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Ccr7), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c0: Naive CD8 T"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Tyrobp")+
  theme(legend.position = "right")
g11

g14 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Tcf7), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c1: NKT-like","c5: Th","c3: Th2","c0: Naive CD8 T"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Tcf7")+
  theme(legend.position = "right")
g14

g15 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Lag3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c4: Tem/ex","c6: Tex"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Lag3")+
  theme(legend.position = "right")
g15

g16 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Havcr2), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c4: Tem/ex","c6: Tex"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Havcr2")+
  theme(legend.position = "right")
g16

g17 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Pdcd1), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c4: Tem/ex","c6: Tex"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Pdcd1")+
  theme(legend.position = "right")
g17

g18 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Id2), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  # ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c3: Th2","c0: Naive CD8 T"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Id2")+
  theme(legend.position = "right")
g18

g12 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Il13), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(cluster %in% c("c3: Th2"), cluster, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = 2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Il13")+
  theme(legend.position = "right")
g12

g20 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Itgal), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Itgal (LFA-1)")+
  theme(legend.position = "right")
g20

g <- gridExtra::grid.arrange(grobs = list(g1, g2, g13, g3, g4, g5, g6, g7, g8, g19,g11, g9, g10,g14, g15, g16,g17,g18, g12,g20), ncol = 4)
ggsave(paste0(outdir, "T_marker_mascarade_ann-V1_featureplot.pdf"),plot = g, width = 15, height = 12 )

### T cell subset violin plot

T_focus <- subset(T_cells, subset = annotation_V1 %in% c("c0: Naive CD8 T","c2: Treg","c4: Tem/ex","c6: Tex"))

T_focus@meta.data <- T_focus@meta.data %>% mutate(
  label = case_match(
    annotation_V1,
    "c0: Naive CD8 T" ~ "Naive CD8T",
    "c2: Treg" ~ "Treg",
    "c4: Tem/ex" ~ "Tem/ex",
    "c6: Tex" ~ "Tex",
    .default = annotation_V1
  )
)
T_focus$label = factor(T_focus$label, levels = c("Treg", "Tex", "Naive CD8T", "Tem/ex"))

g <- VlnPlot(T_focus, features = c("Cd4", "Cd8a", "Foxp3", "Sell","Pdcd1"), group.by = "label", cols = pal2[c(3,7,1,5)])
g
ggsave(filename = paste0(outdir,"June23_T-marker_VlnPlot.pdf"), plot = g, width = 8 , height = 6)

### cell state heatmap

# markers from LCMV
# https://www.pnas.org/doi/10.1073/pnas.2221985120

inhibitory_receptor <- c("Lag3","Pdcd1","Cd160","Ctla4","Cd244","Havcr2")
costimu <- c("Tnfsf14","Cd28","Tnfsf4","Icos","Tnfsf9")
self_renewel <- c("Prickle1","Kit","Axin2","Dvl2","Sox4","Lef1")
chemokine_and_receptor <- c("Cxcl11","Xcl1","Cxcl10","Ccl5","Ccl4","Ccl3","Csf1","Ccr5","Cxcr4","Cxcr5")
Cytokine_and_effector <- c("Fasl","Ifng","Tnfsf10","Prf1","Gzma","Gzmb","Il10","Il2","Tnf")
Tnx_factor <- c("Plagl1","Bcl6","Id3", "Foxo1","Id2","Prdm1","Eomes","Tbx21","Tox")
Mem_precursors <- c("Sell", "Il7r","Il2rb", "Klrg1")

avg_T <- AggregateExpression(T_cells, assays = "RNA",features = c(inhibitory_receptor,costimu,self_renewel,chemokine_and_receptor,Cytokine_and_effector,Tnx_factor, Mem_precursors), return.seurat = T,group.by = c("annotation_V1","orig.ident"))

avg_T$annotation_V1 <- factor(avg_T$annotation_V1, levels = c("g Naive CD8 T",
                                                              "g NKT-like",
                                                              "g Treg",
                                                              "g Th2",
                                                              "g Tem/ex",
                                                              "g Th",
                                                              "g Tex",
                                                              "g gdT",
                                                              "g Proliferating"))
avg_T$sample <- rep(c("WT_HC","KO_HC"),9)
avg_T$sample <- factor(avg_T$sample, levels = c("WT_HC","KO_HC"))
T_pal <- pal2[1:9]
names(T_pal) <- unique(unname(avg_T$annotation_V1))
col = list(
  Cell_type = T_pal,
  Sample = c("WT_HC"="red3", "KO_HC" = "darkblue")
)
ha <- HeatmapAnnotation(
  Cell_type = unname(avg_T$annotation_V1), Sample = unname(avg_T$sample),
  col = col
)
mat <-  GetAssayData(avg_T,layer = "scale.data")
Heatmap(mat,name = "T cells",top_annotation = ha,cluster_rows = F,cluster_columns = F,show_row_names = TRUE,show_column_names = F,row_names_side = "left",row_names_gp = grid::gpar(fontsize = 6),col = colorRamp2(c(-2,0,2),c("darkblue","white","red")))


pdf(paste0(outdir,"T_marker_Heatmap_LCMV_PNAS-panels.pdf"))
Heatmap(mat,name = "T cells",top_annotation = ha,cluster_rows = F,cluster_columns = F,show_row_names = TRUE,show_column_names = F,row_names_side = "left",row_names_gp = grid::gpar(fontsize = 6),col = colorRamp2(c(-2,0,2),c("darkblue","white","red")))
dev.off()

avg_T <- AggregateExpression(subset(T_cells, subset = annotation_V1 %in% c("c0: Naive CD8 T", "c4: Tem/ex") ), assays = "RNA",features = c(inhibitory_receptor,costimu,self_renewel,chemokine_and_receptor,Cytokine_and_effector,Tnx_factor, Mem_precursors), return.seurat = T,group.by = c("annotation_V1","orig.ident"))

avg_T$annotation_V1 <- factor(avg_T$annotation_V1, levels = c("g Naive CD8 T",
                                                              
                                                              "g Tem/ex"
))
avg_T$sample <- rep(c("WT_HC","KO_HC"),2)
avg_T$sample <- factor(avg_T$sample, levels = c("WT_HC","KO_HC"))
T_pal <- pal2[c(1,5)]
names(T_pal) <- unique(unname(avg_T$annotation_V1))
col = list(
  Cell_type = T_pal,
  Sample = c("WT_HC"="red3", "KO_HC" = "darkblue")
)
ha <- HeatmapAnnotation(
  Cell_type = unname(avg_T$annotation_V1), Sample = unname(avg_T$sample),
  col = col
)
mat <-  GetAssayData(avg_T,layer = "scale.data")
Heatmap(mat,name = "T cells",top_annotation = ha,cluster_rows = F,cluster_columns = F,show_row_names = TRUE,show_column_names = F,row_names_side = "left",row_names_gp = grid::gpar(fontsize = 6),col = colorRamp2(c(-2,0,2),c("darkblue","white","red")))


pdf(paste0(outdir,"T_marker_Heatmap_LCMV_PNAS-panels_CD8-only.pdf"), width = 4, height = 10)
Heatmap(mat,name = "CD8 T cells",top_annotation = ha,cluster_rows = F,cluster_columns = F,show_row_names = TRUE,show_column_names = F,row_names_side = "left",row_names_gp = grid::gpar(fontsize = 6),col = colorRamp2(c(-2,0,2),c("darkblue","white","red")))
dev.off()

avg_T <- AggregateExpression(subset(T_cells, subset = annotation_V1 %in% c("c0: Naive CD8 T", "c4: Tem/ex") ), assays = "RNA",features = c(inhibitory_receptor,costimu,self_renewel,chemokine_and_receptor,Cytokine_and_effector,Tnx_factor, Mem_precursors), return.seurat = T,group.by = c("orig.ident"))

LCMV_like_exp <- GetAssayData(avg_T, slot = "scale.data")
LCMV_like_exp_gene <- rownames(LCMV_like_exp)
LCMV_like_exp <- cbind(LCMV_like_exp_gene, LCMV_like_exp)
colnames(LCMV_like_exp)[1] <- "gene"

write.table(LCMV_like_exp,file = paste0(outdir,"LCMV_WT-vs-KO_scaled_exp_for_heatmap.tsv"), sep = "\t", row.names = F,col.names = T,quote = F )

avg_T <- AggregateExpression(subset(T_cells, subset = annotation_V1 %in% c("c0: Naive CD8 T", "c4: Tem/ex") ), assays = "RNA",features = c(inhibitory_receptor,costimu,self_renewel,chemokine_and_receptor,Cytokine_and_effector,Tnx_factor, Mem_precursors), return.seurat = T,group.by = c("orig.ident"))

avg_T$sample <- rep(c("WT_HC","KO_HC"),2)
avg_T$sample <- factor(avg_T$sample, levels = c("WT_HC","KO_HC"))

col = list(
  Sample = c("WT_HC"="red3", "KO_HC" = "darkblue")
)
ha <- HeatmapAnnotation(
  Sample = unname(avg_T$sample),
  col = col
)
mat <-  GetAssayData(avg_T,layer = "scale.data")
Heatmap(mat,name = "CD8 T cells",top_annotation = ha,cluster_rows = F,cluster_columns = F,show_row_names = TRUE,show_column_names = F,row_names_side = "left",row_names_gp = grid::gpar(fontsize = 6),col = colorRamp2(c(-1,0,1),c("darkblue","white","red3")))


pdf(paste0(outdir,"T_marker_Heatmap_LCMV_PNAS-panels_WT-vs_KO_CD8.pdf"), width = 3, height = 10)
Heatmap(mat,name = "CD8 T cells",top_annotation = ha,cluster_rows = F,cluster_columns = F,show_row_names = TRUE,show_column_names = F,row_names_side = "left",row_names_gp = grid::gpar(fontsize = 6),col = colorRamp2(c(-1,0,1),c("darkblue","white","red3")))
dev.off()

### stem-like gene expression comparison 

T_cells$orig.ident <- factor(T_cells$orig.ident, levels = c("WT_HC","KO_HC"))

markers <- c("Pdcd1","Tcf7","Lef1","Xcl1","Id3","Slamf6")
glist <- lapply(markers, function(x){
  VlnPlot(T_cells, group.by = "annotation_V1", features = x, split.by = "orig.ident")+scale_fill_manual(values = c("WT_HC"="red3", "KO_HC" = "darkblue"))+xlab("")+theme(axis.text.x = element_text(size = 8))
})

g <- gridExtra::grid.arrange(grobs = glist, ncol = 3)
ggsave(filename = paste0(outdir,"stemness-exhaustion-marker-comparison_all-T-cells.pdf"), plot = g, width = 12, height = 8)

glist <- lapply(markers, function(x){
  VlnPlot(subset(T_cells, subset = annotation_V1 %in% c("c0: Naive CD8 T", "c4: Tem/ex")), group.by = "annotation_V1", features = x, split.by = "orig.ident")+scale_fill_manual(values = c("WT_HC"="red3", "KO_HC" = "darkblue"))+theme(axis.text.x = element_text(size = 8))+xlab("")
})

g <- gridExtra::grid.arrange(grobs = glist, ncol = 3)
ggsave(filename = paste0(outdir,"stemness-exhaustion-marker-comparison_CD8-only.pdf"), plot = g, width = 8 , height = 6)



