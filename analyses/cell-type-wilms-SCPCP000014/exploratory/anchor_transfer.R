library(Seurat)
library(dplyr)

# ######## atlas
# path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
# sce <- zellkonverter::readH5AD(paste0(path_ref, "/Fetal_full_v3.h5ad"))
# seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(sce),
#                                                assay = "RNA",
#                                                project = "kidneyatlas")
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
# obj <- seurat_obj
# # log transform counts
# obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
# obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
# obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(object = obj))
# obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))
# SeuratObject::SaveSeuratRds(obj, file = paste0("results/kidneyatlas.h5Seurat"))
obj <- SeuratObject::LoadSeuratRds(paste0("results/kidneyatlas.h5Seurat"))

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


######### dataset

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
anchors <- FindTransferAnchors(reference = obj, query = sample_obj)

# clean annotation
obj@meta.data <- obj@meta.data %>%
  mutate(annot = case_when(compartment == "stroma" ~ "stroma",
                           TRUE ~ celltype))
obj@meta.data$annot <- factor(obj@meta.data$annot)
# transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = obj$annot
)
sample_obj <- AddMetaData(object = sample_obj, metadata = predictions)
Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", label = F)
Seurat::DimPlot(sample_obj, reduction = "umap", label = T)