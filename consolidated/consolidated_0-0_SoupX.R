library(SoupX)#Soup removal
library(Seurat)
library(tidyverse)
library(Matrix)
library(DropletUtils)

######### define palette #######
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
#################################

################# Run SoupX ###################

# code credit: https://github.com/gatelabNW/csf_aging/blob/main/code/0_preprocessing/02_preprocessing_soupx.R

# ------------------------------------------------------------------------------
# Define inputs
gex_dir <- "align_V8/"
soupx_dir <- "SoupX/"

# Create output directory
dir.create(soupx_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Run SoupX

# Create lists of directories to load as SoupX objects and Seurat objects
sample_dirs <- list.dirs(gex_dir, recursive = FALSE)

# Initialize gene list to estimate contamination fraction
# Selecting for monocyte/dendritic markers which are highly specific and abundant
M.genes <- c("Cd14", "Cd68","Adgre1", "Ms4a7")

# Initialize contamination fraction dataframe 
contamination_frac <- setNames(data.frame(matrix(ncol = 2, nrow = 0)),
                               c("sample", "contamination_frac"))

# For each sample...
for (dir in sample_dirs) {
  # Isolate sample name
  sample <- unlist(strsplit(dir, "/")) %>%
    tail(1)
  
  # Print what sample is being processed
  print(paste0("Processing sample ", sample))
  
  # Load in SoupX and Seurat data
  soupx <- load10X(paste0(dir, "/outs"))
  seurat <- Read10X(paste0(dir, "/outs/filtered_feature_bc_matrix")) %>%
    CreateSeuratObject()
  seurat
  
  # Normalize and run PCA on Seurat object
  seurat <-  seurat  %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
  
  # Run TSNE clustering
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
  seurat <- FindNeighbors(seurat, dims = 1:30)
  seurat <- FindClusters(seurat, resolution = 0.3)
  seurat <- RunTSNE(seurat, dims = 1:30)
  
  # Add cluster info to SoupX object
  soupx <- setClusters(soupx, setNames(seurat$seurat_clusters, colnames(seurat)))
  
  # Estimate contamination fraction
  useToEst <- estimateNonExpressingCells(soupx, nonExpressedGeneList = list(M.genes = M.genes))
  soupx <- calculateContaminationFraction(soupx,
                                          list(M.genes = M.genes),
                                          useToEst = useToEst)
  
  # Save contamination fraction
  contamination_frac[nrow(contamination_frac) + 1, ] <- c(sample, soupx$metaData$rho[1])
  
  # Create adjusted counts
  adj_counts <- adjustCounts(soupx)
  
  # Save Corrected Counts
  dir.create(paste0(soupx_dir, "/", sample), showWarnings = FALSE)
  DropletUtils:::write10xCounts(paste0(soupx_dir, "/", sample), adj_counts, overwrite = TRUE)
}

# Print contamination fractions
contamination_frac