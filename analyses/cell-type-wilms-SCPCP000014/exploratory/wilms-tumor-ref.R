library(dplyr)
library(Seurat)
library(SingleCellExperiment)

######## sce object preparation ######## 
path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
path_supp_table <- paste0(path_ref, "/aat1699-young-tabless1-s12-revision2.xlsx")
cell_manifest <- readxl::read_xlsx(path_supp_table, sheet = 11, skip = 1) %>%
  select(c(barcode, SangerID, DropletID, ClusterID, QCpass, Source))
cluster_info <- readxl::read_xlsx(path_supp_table, sheet = 2, skip = 1) %>%
  select(c(Cluster_ID,Category,Cell_type1,Cell_type2,Cell_type3))
# supp object
path_supp_mtx <- paste0(path_ref, "/tableOfCounts.mtx")
path_supp_mtx_col <- paste0(path_ref, "/tableOfCounts_colLabels.tsv")
path_supp_mtx_row <- paste0(path_ref, "/tableOfCounts_rowLabels.tsv")
# =======================
supp_mtx <- Matrix::readMM(path_supp_mtx)
supp_mtx <- as(supp_mtx, "CsparseMatrix")
# =======================
coldata <- read.table(path_supp_mtx_col, header = T) %>% 
  left_join(cell_manifest, by = c("SangerID","DropletID","Barcode"="barcode")) %>%
  left_join(cluster_info, by = c("ClusterID" = "Cluster_ID"))
rownames(coldata) <- coldata$DropletID
rowdata <- read.table(path_supp_mtx_row, header = T)
rownames(rowdata) <- rowdata$GeneLabel
dimnames(supp_mtx) <- list(rowdata$GeneLabel, coldata$DropletID)

sce <- SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = supp_mtx), 
                                                  colData = coldata,
                                                  rowData = rowdata)
# remove qc fail & un-annotated cells
sce <- sce[, sce$QCpass == "TRUE"]
sce <- sce[, !is.na(sce$Cell_type1)]
sce <- sce[, !is.na(sce$Category)]
sce <- as(sce, "SingleCellExperiment")
######## convert ensembl IDs to symbols and dedup
rownames(sce) <- sce@rowRanges@elementMetadata@listData[["Symbol"]]
sce <- singleCellTK::dedupRowNames(sce, as.rowData = F, return.list = F)
sce <- sce[!is.na(rownames(sce)),]
sce <- sce[!duplicated(rownames(sce)),]

readr::write_rds(sce, paste0(path_ref, "/aat1699-young.rds"))

######## create seurat obj for wilms ######## 
path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"

obj <- readr::read_rds( paste0(path_ref, "/aat1699-young.rds"))
seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(obj),
                                               assay = "RNA",
                                               project = "aat1699-young")
# convert colData and rowData to data.frame for use in the Seurat object
cell_metadata <- as.data.frame(SingleCellExperiment::colData(obj))
row_metadata <- as.data.frame(SingleCellExperiment::rowData(obj))
# add cell metadata (colData) from SingleCellExperiment to Seurat
seurat_obj@meta.data <- cell_metadata
# add row metadata (rowData) from SingleCellExperiment to Seurat
seurat_obj[["RNA"]]@meta.data <- row_metadata
# add metadata from SingleCellExperiment to Seurat
seurat_obj@misc <- S4Vectors::metadata(obj)
obj <- seurat_obj[,grepl("Wilms",seurat_obj$Source)]
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(object = obj))
obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))

SeuratObject::SaveSeuratRds(obj, file = paste0("results/aat1699-young_wilms.h5Seurat"))

######## marker gene ######## 
obj <- SeuratObject::LoadSeuratRds(paste0("results/aat1699-young_wilms.h5Seurat"))

coldata <- obj@meta.data %>% as.data.frame()
# calculate marker genes for compartments
Idents(obj) <- obj@meta.data[["Category"]]
markers_all <- FindMarkers(object = obj, 
                           only.pos = TRUE,
                           assay = "RNA",
                           ident.1 = c("Kidney_tumour","Kidney_tumour_immune"),
                           logfc.threshold = 0.25) 
write.csv(markers_all, file = paste0("./results/aat1699-young_wilms_degenes_tumor.csv"))

