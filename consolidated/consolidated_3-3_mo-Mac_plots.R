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
outdir <- "V1p3-2_moMac_plots/"
if (!dir.exists(outdir)) dir.create(outdir)

myeloid <- readRDS(paste0(RDSdir,"myeloid_V1p3-2.RDS"))
whole <- readRDS(paste0(RDSdir, "whole_V1p4-2.RDS"))

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


saveRDS(myeloid, file = paste0(RDSdir,"myeloid_V1p3-2.RDS"))

mo_Mac <- subset(myeloid, subset = annotation_V1 %in% c("M0: Mo-Mac1","M1: Mo-Mac2","M2: Mo-Mac3","M4: Mo-Mac4","M7: Mo-Mac5","M8: Proliferating Mo-Mac","M11: Mo-Mac6","M12: Mo-Mac7"))

mo_Mac@meta.data <- mo_Mac@meta.data %>% mutate(
  condition = factor(
    case_match(
      orig.ident,
      "WT_HC" ~ "WT",
      "KO_HC" ~ "MKO"
    ), levels = c("WT", "MKO")
  )
)

g <- VlnPlot(mo_Mac, features = c("Nos2", "Ccr2","Chil3","Mrc1", "Arg1", "Trem2", "Cx3cr1", "Mki67"),pt.size = 0.1, ,group.by = "condition",alpha = 0.05, cols = c("grey90","red3"), ncol = 4)
g
ggsave(paste0(outdir, "Mo-Mac_markers_WT-vs_KO_ann-V1_Vlnplot.pdf"),plot = g, width = 5, height = 5)

g <- VlnPlot(subset(mo_Mac,subset = annotation_V1 %in% c("M7: Mo-Mac5","M11: Mo-Mac6")), features = c("Nos2"),pt.size = 0.1, ,group.by = "annotation_V1", split.by = "condition",alpha = 0.05, cols = c("grey90","red3"), ncol = 4)
g
ggsave(paste0(outdir, "Mo-Mac_Nos_WT-vs_KO_ann-V1_Vlnplot.pdf"),plot = g, width =6, height = 4)



