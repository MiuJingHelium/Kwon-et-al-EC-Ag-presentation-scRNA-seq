library(tidyverse)
library(Seurat)
library(scater)
library(harmony)
library(SeuratWrappers)
library(mascarade)
library(data.table)
library(ComplexHeatmap)
library(circlize)

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

pal2 <- pal[c(1:8,10,13:14,16:20, 22:25,27:30, 33:35 )] # use pal2

RDSdir <- "RDS/"
outdir <- "V1p3_myeloid_plots/"
if (!dir.exists(outdir)) dir.create(outdir)

myeloid <- readRDS(paste0(RDSdir,"myeloid_V1p3-2.RDS"))
whole <- readRDS(paste0(RDSdir, "whole_V1p4-2.RDS"))


outdir <- "V1p3-2_myeloid_plots/"
if (!dir.exists(outdir)) dir.create(outdir)


whole@meta.data <- whole@meta.data %>%
  mutate(highlight = case_match(
    compartment,
    "Myeloid" ~ "Myeloid",
    .default = "NA"
  ))

### overview from whole

DimPlot(whole,group.by = "highlight") + 
  scale_color_manual(values = c("Myeloid" = colorRampPalette(brewer.pal(11, "BrBG"))(8)[3],na.value = "grey90"))
g <- DimPlot(whole,group.by = "highlight") + 
  scale_color_manual(values = c("Myeloid" = colorRampPalette(brewer.pal(11, "BrBG"))(8)[3],na.value = "grey90"))
ggsave(filename = paste0(outdir,"Supp_Whole_Myeloid_unsplit.pdf"), plot = g, units = "cm",height = 10,width = 14)

### proportion plot
sample_imm_size = data.frame(table(whole$orig.ident))
sample_imm_size <- sample_imm_size %>% dplyr::rename(n_immune = Freq, orig.ident = Var1)

clust_prop <- myeloid@meta.data %>% group_by(annotation_V1,compartment,orig.ident) %>% summarise(n_subpop = n()) %>% ungroup()
clust_prop <- left_join(clust_prop,sample_imm_size, by = "orig.ident")
clust_prop <- clust_prop %>% group_by(orig.ident) %>% mutate(percentage = 100*n_subpop/n_immune)

ratio <- clust_prop %>% group_by(annotation_V1) %>% reframe(FC = log2(percentage[orig.ident == "KO_HC"]/percentage[orig.ident == "WT_HC"]))

write.table(clust_prop, file = paste0(outdir,"myeloid_ann-V1_proportion.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(ratio, file = paste0(outdir,"myeloid_ann-V1_log2-proportion.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

g3 <- ggplot(clust_prop,aes(x = factor(orig.ident, levels = c("WT_HC", "KO_HC")), y = percentage, fill = annotation_V1))+
  geom_bar(aes(color = orig.ident), fill = "white",stat = "identity", position = "dodge")+
  scale_color_manual("condition",values = c("KO_HC" = "darkblue", "WT_HC" = "red3"))+theme_classic()+
  ylab("% in immune cells")+
  xlab("conditions")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  facet_wrap(~annotation_V1, scales = "free")
g3  
ggsave(filename = paste0(outdir,"myeloid_percentage_bar_fill_by_condition_ann-V1.pdf"), plot = g3, width = 10, height = 8)



g4 <- ggplot(ratio, aes( x = reorder(annotation_V1, -FC), y = FC))+
  geom_bar(aes(fill = annotation_V1),stat = "identity", position = "dodge")+
  scale_fill_manual(values = pal2[13:26])+theme_classic()+
  ylab("log2 % in KO v.s. WT immune cells")+
  xlab("subpopulations")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  geom_hline(yintercept = c(-1,1), color = "red3", linetype = 3)
g4
ggsave(filename = paste0(outdir,"myeloid_log2_percentage_bar_ann-V1.pdf"), plot = g4, width = 6, height = 4)

### myeloid UMAP and featureplots

myeloid$annotation_V1 <- factor(myeloid$annotation_V1, levels = c("M0: Mo-Mac1",
                                                                  "M1: Mo-Mac2",
                                                                  "M2: Mo-Mac3",
                                                                  "M4: Mo-Mac4",
                                                                  "M7: Mo-Mac5",
                                                                  "M8: Proliferating Mo-Mac",
                                                                  "M11: Mo-Mac6",
                                                                  "M12: Mo-Mac7",
                                                                  "M6: Neutrophil",
                                                                  "M5: Mo-DC/monocytes",
                                                                  "M9: migDC",
                                                                  "M10: cDC1/DC2",
                                                                  "M3: Basophil/granulocytes",
                                                                  "M13: granulocytes/mast cells"
                                                                  
))

maskTable <- generateMask(dims=myeloid@reductions$umap@cell.embeddings, 
                          clusters=myeloid$annotation_V1)
data <- data.table(myeloid@reductions$umap@cell.embeddings, 
                   cluster=myeloid$annotation_V1,
                   myeloid@meta.data)
centers <- data %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::summarize(umap_1 = median(umap_1), umap_2 = median(umap_2)) 


data$orig.ident <- factor(data$orig.ident, levels = c("WT_HC","KO_HC"))

g <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=cluster), size = 0.5) + 
  geom_path(data=maskTable, aes(group=group)) +
  #geom_text(data = centers, aes(label=cluster), size = 3.5, color = "black") +
  scale_color_manual(values = pal2[13:26])+
  coord_fixed() + 
  theme_classic()+
  facet_wrap(~orig.ident)+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
g
ggsave(filename = paste0(outdir,"myeloid_V1p3-2_ann-V1_umap_split_mascarade.pdf"), plot = g, width = 18, height = 7)

g <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=cluster), size = 0.5) + 
  geom_path(data=maskTable, aes(group=cluster)) +
  ggrepel::geom_text_repel(data = centers, aes(label=cluster), size = 5, color = "black",min.segment.length = unit(0, 'lines'),force = 5, label.padding = 5,max.overlaps = Inf,direction = "both") +
  scale_color_manual(values = pal2[13:26])+
  coord_fixed() + 
  theme_classic()+theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
g
ggsave(filename = paste0(outdir,"myeloid_V1p3-2_ann-V1_umap_mascarade.pdf"), plot = g, width = 18, height = 7)

dat <- GetAssayData(myeloid, slot = "data")

data$barcode <- rownames(myeloid@meta.data)

markers <- c("Cd24a","Ccr3","S100a8","Cd200r3",
             "Flt3","Sirpa","Xcr1","Siglech","Zbtb46","Ccr2","Clec9a","Clec10a","Ccr7",
             "Mki67",
             "Adgre1","Mertk","Ifit3","Rsad2","Hopx",
             "Arg1","Mrc1","Chil3","Nos2","Cxcl3","Trem2")
marker_exp <- as.data.frame(t(as.matrix(dat[rownames(dat) %in% markers,])))

marker_exp$barcode <- rownames(marker_exp)

data <- left_join(data, marker_exp, by = "barcode")

table(data$cluster)
centers$annotation_V1 <- as.character(centers$cluster)

g1 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd24a), size = 0.5)+
  geom_path(data=maskTable, aes(group=group), color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(annotation_V1 == "M3: Basophil/granulocytes", annotation_V1, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  scale_fill_manual(values = pal2)+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd24a")+
  theme(legend.position = "right")
g1

g2 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=S100a8), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(annotation_V1 == "M6: Neutrophil", annotation_V1, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("S100a8")+
  theme(legend.position = "right")
g2

g13 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cd200r3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(annotation_V1 == "M13: granulocytes/mast cells", annotation_V1, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cd200r3")+
  theme(legend.position = "right")
g13

##### DC ####

g3 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Flt3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(annotation_V1 %in% c("M5: Mo-DC/monocytes","M9: migDC","M10: cDC1/DC2"), annotation_V1, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Flt3")+
  theme(legend.position = "right")
g3

g4 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Zbtb46), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(annotation_V1 %in% c("M5: Mo-DC/monocytes","M9: migDC","M10: cDC1/DC2"), annotation_V1, "")), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Zbtb46")+
  theme(legend.position = "right")
g4

g5 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Xcr1), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(annotation_V1 %in% c("M10: cDC1/DC2"), annotation_V1, "")), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Xcr1")+
  theme(legend.position = "right")
g5

g6 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Ccr7), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(annotation_V1 %in% c("M9: migDC"), annotation_V1, "")), size = 3.5, color = "black",nudge_x = 0.5, nudge_y = 1) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Ccr7")+
  theme(legend.position = "right")
g6


#### Mo-Mac ####


g7 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Mki67), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(annotation_V1 == "M8: Proliferating Mo-Mac", annotation_V1, "")), size = 3.5, color = "black", nudge_x = 1, nudge_y = -2) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Mki67")+
  theme(legend.position = "right")
g7

g8 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Adgre1), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(grepl("Mo-Mac", annotation_V1), annotation_V1, "")), size = 2, color = "black",box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Adgre1")+
  theme(legend.position = "right")
g8

g19 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Mertk), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(!grepl("Mo-Mac3|Mo-Mac5", annotation_V1) & grepl("Mo-Mac", annotation_V1), annotation_V1, "")), size = 2, color = "black",,box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Mertk")+
  theme(legend.position = "right")
g19

g9 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Ifit3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(grepl("Mo-Mac6", annotation_V1), annotation_V1, "")), size = 3.5, color = "black",box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Ifit3")+
  theme(legend.position = "right")
g9


g10 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Arg1), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(!grepl("Mo-Mac3", annotation_V1) & grepl("Mo-Mac", annotation_V1), annotation_V1, "")), size = 2, color = "black",box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Arg1")+
  theme(legend.position = "right")
g10


g11 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Mrc1), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(!grepl("Mo-Mac3", annotation_V1) & grepl("Mo-Mac", annotation_V1), annotation_V1, "")), size = 2, color = "black",box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Mrc1")+
  theme(legend.position = "right")
g11

g14 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Chil3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(grepl("Mo-Mac", annotation_V1), annotation_V1, "")), size = 2, color = "black",box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Chil3")+
  theme(legend.position = "right")
g14

g15 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Nos2), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(grepl("Mo-Mac5", annotation_V1), annotation_V1, "")), size = 3.5, color = "black",box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Nos2")+
  theme(legend.position = "right")
g15

g16 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Cxcl3), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(grepl("Mo-Mac5", annotation_V1), annotation_V1, "")), size = 3.5, color = "black",box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Cxcl3")+
  theme(legend.position = "right")
g16

g17 <- ggplot(data, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=Trem2), size = 0.5)+
  geom_path(data=maskTable, aes(group=group),color = "grey60") +
  ggrepel::geom_text_repel(data = centers, aes(label=ifelse(!grepl("Mo-Mac5|Mo-Mac3", annotation_V1) & grepl("Mo-Mac", annotation_V1), annotation_V1, "")), size = 2, color = "black",box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_gradientn(colors = c("grey90", "red2","red3"))+
  coord_fixed() + 
  theme_classic()+
  ggtitle("Trem2")+
  theme(legend.position = "right")
g17


g <- gridExtra::grid.arrange(grobs = list(g1, g2, g13, g3, g4, g5, g6, g7, g8, g19,g11, g9, g10,g14, g15, g16,g17), ncol = 5)
ggsave(paste0(outdir, "Myeloid_marker_mascarade_ann-V1_featureplot.pdf"),plot = g, width = 25, height = 20 )

g <- DotPlot(myeloid, features = markers, group.by = "annotation_V1")+scale_color_gradientn(colors = c("grey90", "red2","red3"))+RotatedAxis()
ggsave(paste0(outdir, "Myeloid_marker_ann-V1_Dotplot.pdf"),plot = g, width = 10, height = 6)



### mo-Mac heatmap


mo_Mac <- subset(myeloid, subset = annotation_V1 %in% c("M0: Mo-Mac1","M1: Mo-Mac2","M2: Mo-Mac3","M4: Mo-Mac4","M7: Mo-Mac5","M8: Proliferating Mo-Mac","M11: Mo-Mac6","M12: Mo-Mac7"))

g <- DotPlot(mo_Mac, features = c("Mrc1", "Arg1", "Trem2", "Cx3cr1", "Mki67", "Nos2", "Ccr2","Chil3"), group.by = "annotation_V1",col.min = 0)+scale_color_gradientn(colors = c("grey90", "red2","red3"))+RotatedAxis()
g
ggsave(paste0(outdir, "Mo-Mac_markers_ann-V1_Dotplot_min0.pdf"),plot = g, width = 10, height = 6)

mo_Mac$orig.ident <- factor(mo_Mac$orig.ident, levels = c("WT_HC","KO_HC"))

g <- VlnPlot(mo_Mac, features = c("Mrc1", "Arg1", "Trem2", "Cx3cr1", "Mki67", "Nos2", "Ccr2","Chil3"), group.by = "annotation_V1",alpha = 0.2, split.by = "orig.ident",cols = c("darkblue","red3"))
g
ggsave(paste0(outdir, "Mo-Mac_markers_ann-V1_Vlnplot.pdf"),plot = g, width = 10, height = 12)

g <- VlnPlot(mo_Mac, features = c("Mrc1", "Arg1", "Trem2", "Cx3cr1", "Mki67", "Nos2", "Ccr2","Chil3"), group.by = "annotation_V1", alpha = 0.2,cols = pal2[13:26][c(1:3,5,8,9,12,13)])
g
ggsave(paste0(outdir, "Mo-Mac_markers_ann-V1_Vlnplot_unsplit.pdf"),plot = g, width = 10, height = 12)

g <- VlnPlot(mo_Mac, features = c("Mrc1", "Arg1", "Trem2", "Cx3cr1", "Mki67", "Nos2", "Ccr2","Chil3"), group.by = "orig.ident",alpha = 0.2, cols = c("darkblue","red3"))
g
ggsave(paste0(outdir, "Mo-Mac_markers_WT-vs_KO_ann-V1_Vlnplot.pdf"),plot = g, width = 8, height = 8)

g <- VlnPlot(mo_Mac, features = c("Mrc1"), group.by = "annotation_V1",alpha = 0.2, split.by = "orig.ident",cols = c("darkblue","red3"))
g
ggsave(paste0(outdir, "Mo-Mac_markers_ann-V1_Vlnplot_example-legend.pdf"),plot = g, width = 4, height = 4)



