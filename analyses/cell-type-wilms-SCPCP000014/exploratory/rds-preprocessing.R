library(dplyr)
library(Seurat)
path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)
db_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/results/doublet-detection/SCPCP000014"

########### per sample rds preprocessing ##########

sample <- "SCPCS000517"; library <- "SCPCL000849"
rds <- readRDS(paste0(path_proj,"/",sample,"/",library,"_processed.rds"))
# out_rds <- paste0(path_proj,"/",sample,"/",library,"_processed_genesym.rds")
db <- read.table(paste0(db_proj,"/",sample,"/",library,"_processed_scdblfinder.tsv"), header = T)

######## convert ensembl IDs to symbols and dedup
rownames(rds) <- rds@rowRanges@elementMetadata@listData[["gene_symbol"]]
rds <- singleCellTK::dedupRowNames(rds, as.rowData = F, return.list = F)
rds <- rds[!is.na(rownames(rds)),]
rds <- rds[!duplicated(rownames(rds)),]
# readr::write_rds(rds, out_rds, compress = "bz2")

######## remove doublets
rds <- rds[,which(db$class == "singlet")]

######## create seurat object from the SCE counts matrix
seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(rds),
                                               assay = "RNA",
                                               project = library)
# convert colData and rowData to data.frame for use in the Seurat object
cell_metadata <- as.data.frame(SingleCellExperiment::colData(rds))
row_metadata <- as.data.frame(SingleCellExperiment::rowData(rds))
# add cell metadata (colData) from SingleCellExperiment to Seurat
seurat_obj@meta.data <- cell_metadata
# add row metadata (rowData) from SingleCellExperiment to Seurat
seurat_obj[["RNA"]]@meta.data <- row_metadata
# add metadata from SingleCellExperiment to Seurat
seurat_obj@misc <- S4Vectors::metadata(rds)

source(file = "exploratory/functions.R")
seurat_obj <- pre_seuratobj(seurat_obj, nfeatures = 3000, run_harmony = F, reduction = "pca")

# Seurat::DimPlot(seurat_obj, reduction = "umap", label = T)
# obj <- Seurat::RunTSNE(obj, dims = 1:ndims)
# Seurat::DimPlot(obj, reduction = "tsne")

SeuratObject::SaveSeuratRds(seurat_obj, file = paste0("results/",sample,".h5Seurat"))
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
SeuratDisk::SaveH5Seurat(seurat_obj, file = paste0("results/",sample,".h5seurat"), overwrite = T )

########### merged rds preprocessing ##########

rds <- readRDS(paste0(path_proj,"/SCPCP000014_merged.rds"))
db <- data.frame()
for (i in 1:10) {
  sample <- meta[i,2]
  library <- meta[i,3]
  tmp <- read.table(paste0(db_proj,"/",sample,"/",library,"_processed_scdblfinder.tsv"), header = T) %>%
    mutate(barcodes = paste0(library,"-",barcodes))
  db <- rbind(db, tmp)
}

######## convert ensembl IDs to symbols and dedup
rownames(rds) <- rds@rowRanges@elementMetadata@listData[["gene_symbol"]]
rds <- singleCellTK::dedupRowNames(rds, as.rowData = F, return.list = F)
rds <- rds[!is.na(rownames(rds)),]
rds <- rds[!duplicated(rownames(rds)),]
# readr::write_rds(rds, out_rds, compress = "bz2")

######## remove doublets

rds <- rds[,colnames(rds) %in% db$barcodes[db$class == "singlet"]]
counts <- as(SingleCellExperiment::counts(rds), "dgCMatrix")
# counts <- counts[-which(Matrix::rowSums(counts) == 0),]

######## create seurat object from the SCE counts matrix
seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts,
                                               assay = "RNA",
                                               project = "SCPCP000014_merged")
# convert colData and rowData to data.frame for use in the Seurat object
cell_metadata <- as.data.frame(SingleCellExperiment::colData(rds))
row_metadata <- as.data.frame(SingleCellExperiment::rowData(rds))
# add cell metadata (colData) from SingleCellExperiment to Seurat
seurat_obj@meta.data <- cell_metadata
# add row metadata (rowData) from SingleCellExperiment to Seurat
seurat_obj[["RNA"]]@meta.data <- row_metadata
# add metadata from SingleCellExperiment to Seurat
seurat_obj@misc <- S4Vectors::metadata(rds)

SeuratObject::SaveSeuratRds(seurat_obj, file = paste0("results/SCPCP000014_merged_clean.h5Seurat"))

# preprocess
source(file = "exploratory/functions.R")
obj <- SeuratObject::LoadSeuratRds(paste0("results/SCPCP000014_merged_clean.h5Seurat"))
obj <- pre_seuratobj(obj)
obj_3000 <- pre_seuratobj(obj, nfeatures = 3000)

Seurat::DimPlot(obj_3000, reduction = "umap", label = T)
Seurat::DimPlot(obj_3000, reduction = "umap", group.by = "library_id")
Seurat::DimPlot(obj_3000, reduction = "umap", group.by = "library_id", split.by= "library_id", ncol = 3)
# obj <- Seurat::RunTSNE(obj, dims = 1:ndims)
# Seurat::DimPlot(obj, reduction = "tsne")

SeuratObject::SaveSeuratRds(obj_3000, file = paste0("results/SCPCP000014_merged.h5Seurat"))
obj <- SeuratObject::LoadSeuratRds(paste0("results/SCPCP000014_merged.h5Seurat"))
obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay")
SeuratDisk::SaveH5Seurat(obj, file = paste0("results/SCPCP000014_merged.h5seurat") )

########### merged by patient (one tumor + one PDX) ##########
# sample_obj <- SeuratObject::LoadSeuratRds(paste0("results/SCPCP000014_merged.h5Seurat"))
# # Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "library_id", split.by= "participant_id", ncol = 3)
# sample_obj <- sample_obj[,sample_obj$participant_id == "SJWLM046146" ]
# sample_obj <- Seurat::NormalizeData(sample_obj, normalization.method = "LogNormalize")
# sample_obj <- Seurat::FindVariableFeatures(sample_obj, selection.method = "vst", nfeatures = 2000)
# sample_obj <- Seurat::ScaleData(sample_obj, features = Seurat::VariableFeatures(object = sample_obj))
# # obj <- Seurat::SCTransform(obj)
# sample_obj <- Seurat::RunPCA(sample_obj, features = Seurat::VariableFeatures(object = sample_obj))
# sample_obj <- harmony::RunHarmony(sample_obj, group.by.vars = "library_id")
# ndims <- 50
# sample_obj <- Seurat::FindNeighbors(sample_obj, dims = 1:ndims, reduction = "harmony")
# sample_obj <- Seurat::FindClusters(sample_obj, resolution = 0.8, algorithm = 1)
# sample_obj <- Seurat::RunUMAP(sample_obj, dims = 1:ndims, reduction = "harmony")
# Seurat::DimPlot(sample_obj, reduction = "umap", label = T, group.by = "library_id")

######## atlas ######## 
path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
sce <- zellkonverter::readH5AD(paste0(path_ref, "/Fetal_full_v3.h5ad"))
seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(sce),
                                               assay = "RNA",
                                               project = "kidneyatlas")
# convert colData and rowData to data.frame for use in the Seurat object
cell_metadata <- as.data.frame(SingleCellExperiment::colData(sce))
row_metadata <- as.data.frame(SingleCellExperiment::rowData(sce))
# add cell metadata (colData) from SingleCellExperiment to Seurat
seurat_obj@meta.data <- cell_metadata
# add row metadata (rowData) from SingleCellExperiment to Seurat
seurat_obj[["RNA"]]@meta.data <- row_metadata
# add metadata from SingleCellExperiment to Seurat
seurat_obj@misc <- S4Vectors::metadata(sce)
# make a copy for processing
obj <- seurat_obj
# log transform counts
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(object = obj))
obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))
# ndims <- 50
# obj <- Seurat::FindNeighbors(obj, dims = 1:ndims)
# obj <- Seurat::FindClusters(obj, resolution = 0.8, algorithm = 1)
# obj <- Seurat::RunUMAP(obj, dims = 1:ndims)
# Seurat::DimPlot(obj, reduction = "umap", label = T, group.by = "celltype")
SeuratObject::SaveSeuratRds(obj, file = paste0("results/kidneyatlas.h5Seurat"))