library(dplyr)

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
# make a copy for processing
obj <- seurat_obj

######## Normalize, scale, feature selection
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(object = obj))
# obj <- Seurat::SCTransform(obj)
obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))

# Seurat::ElbowPlot(obj, ndims = 50)

######## Clustering and dimentional reduction
ndims <- 50
obj <- Seurat::FindNeighbors(obj, dims = 1:ndims)
obj <- Seurat::FindClusters(obj, resolution = 0.8, algorithm = 1)

obj <- Seurat::RunUMAP(obj, dims = 1:ndims)
Seurat::DimPlot(obj, reduction = "umap")
obj <- Seurat::RunTSNE(obj, dims = 1:ndims)
Seurat::DimPlot(obj, reduction = "tsne")

SeuratObject::SaveSeuratRds(obj, file = paste0("results/",sample,".h5Seurat"))

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
# make a copy for processing
obj <- seurat_obj


######## Normalize, scale, feature selection
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(object = obj))
# obj <- Seurat::SCTransform(obj)
obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))
# Seurat::ElbowPlot(obj, ndims = 50)

######## batch effect correction
obj <- SeuratObject::LoadSeuratRds(paste0("results/SCPCP000014_merged.h5Seurat"))

obj <- harmony::RunHarmony(obj, group.by.vars = "library_id")

######## Clustering and dimentional reduction
ndims <- 50
obj <- Seurat::FindNeighbors(obj, dims = 1:ndims, reduction = "harmony")
obj <- Seurat::FindClusters(obj, resolution = 0.8, algorithm = 1)

obj <- Seurat::RunUMAP(obj, dims = 1:ndims, reduction = "harmony")
Seurat::DimPlot(obj, reduction = "umap", label = T)
Seurat::DimPlot(obj, reduction = "umap", group.by = "library_id")
Seurat::DimPlot(obj, reduction = "umap", group.by = "library_id", split.by= "library_id", ncol = 3)
# obj <- Seurat::RunTSNE(obj, dims = 1:ndims)
# Seurat::DimPlot(obj, reduction = "tsne")

SeuratObject::SaveSeuratRds(obj, file = paste0("results/SCPCP000014_merged.h5Seurat"))

######## atlas
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
SeuratObject::SaveSeuratRds(obj, file = paste0("results/kidneyatlas.h5Seurat"))