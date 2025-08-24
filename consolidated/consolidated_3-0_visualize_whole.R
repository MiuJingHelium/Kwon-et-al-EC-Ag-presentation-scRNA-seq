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
outdir <- "V1p4-2_whole_plots/"
if(!dir.exists(outdir))dir.create(outdir)


whole <- readRDS(paste0(RDSdir, "whole_V1p4-2.RDS"))

### ------ subcompartment plots


whole@meta.data <- whole@meta.data %>%
  mutate(highlight = case_match(
    compartment,
    "Myeloid" ~ "Myeloid",
    .default = "NA"
  ))

DimPlot(whole,group.by = "highlight") + 
  scale_color_manual(values = c("Myeloid" = colorRampPalette(brewer.pal(11, "BrBG"))(8)[3],na.value = "grey90"))
g <- DimPlot(whole,group.by = "highlight") + 
  scale_color_manual(values = c("Myeloid" = colorRampPalette(brewer.pal(11, "BrBG"))(8)[3],na.value = "grey90"))
ggsave(filename = paste0(outdir,"Supp_Whole_Myeloid_unsplit.pdf"), plot = g, units = "cm",height = 10,width = 14)

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
ggsave(filename = paste0(outdir,"Supp_Whole_Lymphocytes_unsplit.pdf"), plot = g, units = "cm",height = 10,width = 14)


whole$annotation_V1 <- factor(whole$annotation_V1, levels = c("c0: Naive CD8 T",
                                                              "c1: NKT-like",
                                                              "c2: Treg",
                                                              "c3: Th2",
                                                              "c4: Tem/ex",
                                                              "c5: Th",
                                                              "c6: Tex",
                                                              "c7: gdT",
                                                              "c8: Proliferating",
                                                              "L0: NK",
                                                              "L10: ILC2",
                                                              "L7: B cells",
                                                              "M0: Mo-Mac1",
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
                                                              "M13: granulocytes/mast cells"))

g <- DimPlot(whole, group.by = "annotation_V1",label = T, repel = T)+
  scale_color_manual(values = pal2)# + ggrepel::geom_text_repel(aes( label = ifelse(grepl("^c[0-9]*|^M7|^M0", annotation_V1), annotation_V1, "")))
g
ggsave(filename = paste0(outdir,"whole_V1p4-2_annotated_umap.pdf"), plot = g, width = 10, height = 7)

g <- DimPlot(whole, group.by = "annotation_V1",label = F, repel = T, split.by = "orig.ident")+
  scale_color_manual(values = pal2)
g
ggsave(filename = paste0(outdir,"whole_V1p4-2_annotated_umap_split.pdf"), plot = g, width = 15, height = 7)

