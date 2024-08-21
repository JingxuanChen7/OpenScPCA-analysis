library(dplyr)
library(Seurat)
####################### wilms tumor markers - Young et al. ####################### 
path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
path_supp_table <- paste0(path_ref, "/aat1699-young-tabless1-s12-revision2.xlsx")
markers <- readxl::read_xlsx(path_supp_table, sheet = 4, skip = 1) %>%
  filter(isKeyGene == 1)
cell_manifest <- readxl::read_xlsx(path_supp_table, sheet = 11, skip = 1) %>%
  select(c(barcode, SangerID, DropletID, ClusterID, QCpass, Source))
cluster_info <- readxl::read_xlsx(path_supp_table, sheet = 2, skip = 1) %>%
  select(c(Cluster_ID,Category,Cell_type1,Cell_type2,Cell_type3))

filter_markers <- left_join(markers, cluster_info, by = c("Cluster" = "Cluster_ID") ) %>%
  filter(!is.na(Category) & !is.na(Cell_type1)) %>%
  filter(grepl("tumour", Category))

####################### kidney atlas ####################### 

path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
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
obj <- SeuratObject::LoadSeuratRds(paste0("results/kidneyatlas.h5Seurat"))
# calculate marker genes for compartments
Idents(obj) <- obj@meta.data[["compartment"]]
markers_all <- FindAllMarkers(object = obj, 
                       only.pos = TRUE,
                       logfc.threshold = 0.25) 
write.csv(markers_all, file = paste0("./results/degenes_compart.csv"))
# calculate marker genes for all celltypes
Idents(obj) <- obj@meta.data[["celltype"]]
markers_celltype <- FindAllMarkers(object = obj, 
                              only.pos = TRUE,
                              logfc.threshold = 0.25) 
write.csv(markers_celltype, file = paste0("./results/degenes_celltype.csv"))

####################### manual gene collection ####################### 
stroma <- c("COL1A2","COL6A3","PDGFRA")
Seurat::FeaturePlot(annoobj, features = stroma)

Seurat::FeaturePlot(annoobj, features = c("REN"))
