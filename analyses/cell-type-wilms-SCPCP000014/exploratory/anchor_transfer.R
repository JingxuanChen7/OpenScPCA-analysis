library(Seurat)
library(dplyr)
library(ggpubr)

######## atlas

ref_obj <- SeuratObject::LoadSeuratRds(paste0("results/kidneyatlas.h5Seurat"))

# ######### wilms tumor
# path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
# sce <- readRDS(paste0(path_ref, "/aat1699-young.rds"))
# # sce <- sce[,grepl( "Wilms", sce$Source)]
# sce <- sce[,!sce$Cell_type1 == "Junk"]
# sce <- sce[,!sce$Cell_type1 == "Normal_cell"]
# seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(sce),
#                                                assay = "RNA",
#                                                project = "wilms")
# # convert colData and rowData to data.frame for use in the Seurat object
# cell_metadata <- as.data.frame(SingleCellExperiment::colData(sce))
# row_metadata <- as.data.frame(SingleCellExperiment::rowData(sce))
# # add cell metadata (colData) from SingleCellExperiment to Seurat
# seurat_obj@meta.data <- cell_metadata
# # add row metadata (rowData) from SingleCellExperiment to Seurat
# seurat_obj[["RNA"]]@meta.data <- row_metadata
# # add metadata from SingleCellExperiment to Seurat
# seurat_obj@misc <- S4Vectors::metadata(sce)
# # make a copy for processing
# wilms_obj <- seurat_obj
# 
# # log transform counts
# wilms_obj <- Seurat::NormalizeData(wilms_obj, normalization.method = "LogNormalize")
# wilms_obj <- Seurat::FindVariableFeatures(wilms_obj, selection.method = "vst", nfeatures = 2000)
# wilms_obj <- Seurat::ScaleData(wilms_obj, features = Seurat::VariableFeatures(object = wilms_obj))
# wilms_obj <- Seurat::RunPCA(wilms_obj, features = Seurat::VariableFeatures(object = wilms_obj))
# 
# # ndims <- 50
# # wilms_obj <- Seurat::FindNeighbors(wilms_obj, dims = 1:ndims)
# # wilms_obj <- Seurat::FindClusters(wilms_obj, resolution = 0.8, algorithm = 1)
# # 
# # wilms_obj <- Seurat::RunUMAP(wilms_obj, dims = 1:ndims)
# # Seurat::DimPlot(wilms_obj, reduction = "umap")


######### dataset per sample ######### 

path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"

sample <- "SCPCS000517"; library <- "SCPCL000849"
sample_obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))

# # find anchors
# anchors <- FindTransferAnchors(reference = wilms_obj, query = sample_obj)
# 
# # transfer labels
# predictions <- TransferData(
#   anchorset = anchors,
#   refdata = wilms_obj$Cell_type1
# )
# sample_obj <- AddMetaData(object = sample_obj, metadata = predictions)
# Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", label = F)
# Seurat::DimPlot(sample_obj, reduction = "umap", label = T)


# find anchors
anchors <- FindTransferAnchors(reference = ref_obj, query = sample_obj)

# clean annotation
ref_obj@meta.data <- ref_obj@meta.data %>%
  mutate(annot = case_when(compartment == "stroma" ~ "stroma",
                           compartment == "immune" ~ "immune",
                           # grepl("S shaped body", celltype) ~ "S shaped body",
                           TRUE ~ celltype))
ref_obj@meta.data$annot <- factor(ref_obj@meta.data$annot)
# transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = ref_obj$annot
)
predictions <- mutate(predictions, predicted.id = case_when(prediction.score.max < 0.5 ~ "Unknown",
                                                            TRUE ~ predicted.id))
sample_obj <- AddMetaData(object = sample_obj, metadata = predictions)

# plot
color <- Polychrome::glasbey.colors( length( unique(predictions$predicted.id) ) +1  )
color <- color[-1]
names(color) <- unique(predictions$predicted.id)
color['Unknown'] <- "gray90"

p1 <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", label = F, cols = color)
p2 <- Seurat::DimPlot(sample_obj, reduction = "umap", label = T)
df <- sample_obj@meta.data %>%
  dplyr::group_by(seurat_clusters, predicted.id) %>%
  dplyr::count(name = "sum")
p3 <- ggplot(df, aes(x = seurat_clusters, y = sum, fill =  predicted.id)) +
  geom_bar(width = 0.5, stat = "identity", position = "fill") +
  scale_fill_manual(values = color)

toprow <- ggpubr::ggarrange(p1, p2, widths = c(1.5,1))
ggpubr::ggarrange(toprow, p3, ncol = 1)
Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", 
                label = F, cols = color, split.by = "predicted.id", ncol = 3, alpha = 0.1)

######### dataset merged ######### 

path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"

sample_obj <- SeuratObject::LoadSeuratRds(paste0("results/SCPCP000014_merged.h5Seurat"))


# find anchors
anchors <- FindTransferAnchors(reference = ref_obj, query = sample_obj)

# clean annotation
ref_obj@meta.data <- ref_obj@meta.data %>%
  mutate(annot = case_when(compartment == "stroma" ~ "stroma",
                           compartment == "immune" ~ "immune",
                           #grepl("S shaped body", celltype) ~ "S shaped body",
                           TRUE ~ celltype))
ref_obj@meta.data$annot <- factor(ref_obj@meta.data$annot)
# transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = ref_obj$annot
)
predictions <- mutate(predictions, predicted.id = case_when(prediction.score.max < 0.5 ~ "Unknown",
                                                            TRUE ~ predicted.id))
sample_obj <- AddMetaData(object = sample_obj, metadata = predictions)

# plot
color <- Polychrome::glasbey.colors( length( unique(predictions$predicted.id) ) +1  )
color <- color[-1]
names(color) <- unique(predictions$predicted.id)
color['Unknown'] <- "gray90"

p1 <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", 
                label = F, cols = color, alpha = 0.1)
Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", 
                      label = F, cols = color, ncol = 3, alpha = 0.1, split.by = "predicted.id")
Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", 
                label = F, cols = color, ncol = 3, alpha = 0.1, split.by = "library_id")
#Seurat::FeaturePlot(sample_obj, reduction = "umap", features = "prediction.score.max", cols = c("blue","red"))
p2 <- Seurat::DimPlot(sample_obj, reduction = "umap", label = T)
df <- sample_obj@meta.data %>%
  dplyr::group_by(seurat_clusters, predicted.id) %>%
  dplyr::count(name = "sum")
p3 <- ggplot(df, aes(x = seurat_clusters, y = sum, fill =  predicted.id)) +
  geom_bar(width = 0.5, stat = "identity", position = "fill") +
  scale_fill_manual(values = color)

toprow <- ggpubr::ggarrange(p1, p2, widths = c(1.5,1))
ggpubr::ggarrange(toprow, p3, ncol = 1)

Seurat::DimPlot(sample_obj, reduction = "umap", label = F)
# var_genes <- VariableFeatures(sample_obj)
var_genes <- list(
  list("PECAM1","PLVAP","TIMP3"),
  list("SIX1","CITED1","PAX2","ALDOB","GLYAT","GPX3","SLC12A1","CLCNKA","ATP6V1B1","KRT7","S100P","UPK1A"),
  list("PTPRC","NKG7","CD3D","MS4A1","CD14","FCGR3A","CPA3","TPSAB1","ITGB3","GP6","HBM","HBZ"),
  list("LUM","PDGFRB","PRELP","TNC"),
  list("WT1","CTNNB1","AMER1","IGF2","NCAM1")
)
var_genes <- setNames(object = var_genes, c("Vasc","DevNephron","Immune","Stroma","Tumor"))
Seurat::DimPlot(sample_obj, reduction = "umap", label = T)
Seurat::DotPlot(sample_obj, features = var_genes, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


