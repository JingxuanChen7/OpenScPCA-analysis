library(dplyr)
library(Seurat)
library(enrichR)
library(ggplot2)
# connect to human database
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
source("exploratory/functions.R")

path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)

sample <- "SCPCS000517"; library <- "SCPCL000849"
obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))
obj <- SeuratObject::LoadSeuratRds(paste0("results/SCPCP000014_merged.h5Seurat"))
# calculate marker genes for compartments
Idents(obj) <- obj@meta.data[["seurat_clusters"]]
markers_all <- FindAllMarkers(object = obj, 
                              only.pos = TRUE,
                              logfc.threshold = 0.25) 
# readr::write_rds(markers_all, file = paste0("results/SCPCP000014_merged_markers.rds"))
markers <- markers_all %>%
  filter(p_val_adj < 0.05 & pct.1 > 0.3 & avg_log2FC > 0.5)
markers <- by(markers$gene, markers$cluster, head, n=20)

markers_plot <- clean_gslist(markers, obj)
Seurat::DotPlot(obj, features = markers_plot, cols = c("blue","red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
DimPlot(obj, reduction = "umap", label = T)

dbslist <- listEnrichrDbs() %>% filter(categoryId == 5)
dbs <- c("Azimuth_2023","CellMarker_2024","PanglaoDB_Augmented_2021")

enriched <- enrichr(markers[['0']], dbs)


