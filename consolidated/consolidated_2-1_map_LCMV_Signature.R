library(tidyverse)
library(Seurat)
library(fgsea)
library(data.table)

RDSdir <- "RDS/"
outdir <- "signature_mapping_LCMV/"

late_TD_sig <- unlist(read.delim(file  = "signature_mapping_LCMV/late_TD_vs_early_stem_signature.txt",sep = "\t",header = F))
early_stem_sig <- unlist(read.delim(file  = "signature_mapping_LCMV/early_stem_vs_late_TD_signature.txt",sep = "\t",header = F))

## ---------- create geneset 
fgsea_LCMV <- list(late_TD = late_TD_sig, early_stem = early_stem_sig)

## --------- perform signature mapping over annotation V1
T_cells <- readRDS(paste0(RDSdir,"T_cells_V1p4.Rds"))

Idents(T_cells) <- "annotation_V1"
Tex_markers <- FindMarkers(T_cells, ident.1 = "c4: Tem/ex", ident.2 = "c0: Naive CD8 T",test.use = "MAST", min.pct = 0.2)
Tex_markers$gene <- rownames(Tex_markers)

DE <- Tex_markers %>% mutate(metric = avg_log2FC) %>% arrange(desc(metric)) %>% dplyr::select(gene,metric)
rankings <- DE$metric
names(rankings) <- DE$gene
#print(rankings)
plot(rankings)
title(main = "Tem/ex v.s naive")
fgseaRes <- fgsea(fgsea_LCMV, stats = rankings,nPermSimple = 100000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>% filter(padj < 0.05 ) %>%
  arrange(desc(NES))

fwrite(fgseaResTidy, file=paste0(outdir,"LCMV_fgsea_V1p4_annotation.tsv"), sep="\t", sep2=c("", " ", ""))

pdf(paste0(outdir,"V1p4_Tem-ex_vs_Tn_late-TD.pdf"), height = 3, width = 5)
plotEnrichment(fgsea_LCMV[["late_TD"]],rankings
) + labs(title="late_TD")
dev.off()
pdf(paste0(outdir,"V1p4_Tem-ex_vs_Tn_early-stem.pdf"), height = 3, width = 5)
plotEnrichment(fgsea_LCMV[["early_stem"]],rankings
) + labs(title="early_stem")
dev.off()

## --------- perform signature mapping over condition DE

Idents(T_cells) <- "orig.ident"
cond_markers <- FindMarkers(T_cells, ident.1 = "WT_HC", ident.2 = "KO_HC",test.use = "MAST", min.pct = 0.2)

cond_markers$gene <- rownames(cond_markers)

DE <- cond_markers %>% mutate(metric = avg_log2FC) %>% arrange(desc(metric)) %>% dplyr::select(gene,metric)
rankings <- DE$metric
names(rankings) <- DE$gene
#print(rankings)
plot(rankings)
title(main = "WT v.s Myct1-KO")
fgseaRes <- fgsea(fgsea_LCMV, stats = rankings,nPermSimple = 100000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>% filter(padj < 0.05 ) %>%
  arrange(desc(NES))

fwrite(fgseaResTidy, file=paste0(outdir,"LCMV_fgsea_WT-vs-KO.tsv"), sep="\t", sep2=c("", " ", ""))

pdf(paste0(outdir,"V1p4_WT_vs_Myct1-KO_late-TD.pdf"), height = 3, width = 5)
plotEnrichment(fgsea_LCMV[["late_TD"]],rankings
) + labs(title="late_TD")
dev.off()
pdf(paste0(outdir,"V1p4_WT_vs_Myct1-KO_early-stem.pdf"), height = 3, width = 5)
plotEnrichment(fgsea_LCMV[["early_stem"]],rankings
) + labs(title="early_stem")
dev.off()


