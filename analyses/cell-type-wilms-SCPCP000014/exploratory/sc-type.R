library("openxlsx")
library("dplyr")
library("Seurat")
library("HGNChelper")
library("ggpubr")

path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)

sample <- "SCPCS000517"; library <- "SCPCL000849"
obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Kidney" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
gs_list <- gs_list$gs_positive

###### marker genes from atlas-compart
markers_file <- read.csv(paste0("./results/degenes_compart.csv"))
markers <- markers_file %>%
  filter(p_val_adj < 0.05 & pct.1 > 0.5 & pct.2 < 0.1 )
gs_list <- by(markers$gene, markers$cluster, head, n=10)
gs_list$fetal_nephron <- gs_list$fetal_nephron[gs_list$fetal_nephron != "EPCAM"]
# marker genes from atlas-celltype
markers_file <- read.csv(paste0("./results/degenes_celltype.csv"))
markers <- markers_file %>%
  filter(p_val_adj < 0.05 & pct.1 > 0.5 & pct.2 < 0.1 )
gs_list <- by(markers$gene, markers$cluster, head, n=10)


# remove genes that not present in object and duplicate hits
dupgenes <- unique(unlist(gs_list)[duplicated(unlist(gs_list))])
availgenes <- annoobj@assays[["SCT"]]@data@Dimnames[[1]]
for (i in 1:length(gs_list)) {
  #name <- names(gs_list[i])
  element <- gs_list[[i]][gs_list[[i]] %in% availgenes]
  element <- element[!element %in% dupgenes]
  gs_list[[i]] <- element
  #names(gs_list[i]) <- name
}



# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- obj@assays[["SCT"]]$scale.data %>% as.matrix()

########### run ScType code ########### 
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. pbmc[["RNA"]]$counts or pbmc[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
#Seurat::DimPlot(obj, reduction = "umap", group.by = "seurat_clusters")
cL_resutls <- do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"

# add annotation
annoobj <- obj
annoobj@meta.data <- obj@meta.data %>%
  left_join(sctype_scores[,1:2], by = c("seurat_clusters" = "cluster")) %>%
  mutate(scType = factor(type)) %>%
  data.frame(row.names = .$barcodes)
p1 <- Seurat::DimPlot(annoobj, reduction = "umap", group.by = "scType", label = T) + 
  theme(legend.position = "bottom") + 
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
p2 <- Seurat::DimPlot(annoobj, reduction = "umap", group.by = "seurat_clusters") + 
  theme(legend.position = "bottom")
ggarrange(p1, p2, ncol = 2)

# # sample level prediction pretty messy
# per_sample <- apply(es.max, MARGIN = 2, which.max) %>%
#   data.frame(scType_sample = rownames(es.max)[.]) %>%
#   tibble::rownames_to_column(var = "barcodes")
# annoobj@meta.data <- annoobj@meta.data %>%
#   left_join(per_sample[,c(1,3)], by = c("barcodes")) %>%
#   mutate(scType_sample = factor(scType_sample)) %>%
#   data.frame(row.names = .$barcodes)
# Seurat::DimPlot(annoobj, reduction = "umap", group.by = "scType_sample") +
#   theme(legend.position = "right") +
#   guides(color = guide_legend(ncol=1, override.aes = list(size = 3)))

## dot plot of marker genes
markers <- gs_list[unique(sctype_scores$type)]
DotPlot(annoobj, features = markers, group.by = "scType", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))



