library(dplyr)

path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)

sample <- "SCPCS000517"; library <- "SCPCL000849"
rds <- readRDS(paste0(path_proj,"/",sample,"/",library,"_processed.rds"))

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
obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))

# Seurat::ElbowPlot(obj, ndims = 50)

######## Clustering and dimentional reduction
ndims <- 20
obj <- Seurat::FindNeighbors(obj, dims = 1:ndims)
obj <- Seurat::FindClusters(obj, resolution = 0.8, algorithm = 1)

obj <- Seurat::RunUMAP(obj, dims = 1:ndims)
Seurat::DimPlot(obj, reduction = "umap")
obj <- Seurat::RunTSNE(obj, dims = 1:ndims)
Seurat::DimPlot(obj, reduction = "tsne")

SeuratObject::SaveSeuratRds(obj, file = paste0("results/",sample,".h5Seurat"))

# ######## automated cell annotation - singleR
# celldex::surveyReferences()
# blueprint_ref <- celldex::BlueprintEncodeData(ensembl = TRUE)
# hpca_ref <- celldex::HumanPrimaryCellAtlasData(ensembl = TRUE)
# count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
# #count_mat@Dimnames[[1]] <- seurat_obj[["RNA"]]@meta.data$gene_symbol
# pred.obj <- SingleR::SingleR(test = count_mat, 
#                     ref = blueprint_ref, 
#                     assay.type.test=1,
#                     labels = blueprint_ref$label.main)
# SingleR::plotScoreHeatmap(pred.obj)
# SingleR::plotDeltaDistribution(pred.obj, ncol = 3)
# # add singleR results to seurat obj
# singler_out <- pred.obj$labels
# names(singler_out) <- rownames(pred.obj)
# obj@meta.data$singler <- singler_out
# Seurat::DimPlot(obj, reduction = "umap", group.by = "singler")
# Seurat::DimPlot(obj, reduction = "umap", group.by = "singler", split.by = "singler", ncol = 3)


######## automated cell annotation - singleR - literature ref
obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))
path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
sce <- readRDS(paste0(path_ref, "/aat1699-young.rds"))
sce@NAMES <- sce@elementMetadata@listData[["EnsemblID"]]

# explore reference
ref_coldata <- as(SummarizedExperiment::colData(sce), "DataFrame") %>% 
  as.data.frame() %>%
  filter(startsWith(Source, "Wilms"))
freq <- table(ref_coldata$Category, ref_coldata$Cell_type1) %>% 
  as.data.frame() %>%
  filter(Freq > 0)

# cell exclusion
sce <- sce[, !sce$Cell_type1 == "Junk"]
sce <- sce[, startsWith(sce$Source, "Wilms")]
# cell_types <- c("Wilms_tumour_and_fibroblast","Wilms_tumour","Ureter_epithelium","Plasmacytoid dendritic cell",
#                 "Renal_cell_carcinoma","Plasma cell","Nephron_others")
# sce <- sce[, sce$Cell_type1 %in% cell_types]
#sce <- sce[, sce$Category == "Foetal_kidney" | sce$Category == "Foetal_kidney_immune" | sce$Category == "Kidney_tumour"]

sce <- scuttle::logNormCounts(sce) 
count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
# label category
pred.obj <- SingleR::SingleR(test = count_mat, # normalized count_mat
                             ref = sce, 
                             assay.type.test=1,
                             labels = sce$Cell_type1,
                             num.threads = 6)
SingleR::plotScoreHeatmap(pred.obj)
SingleR::plotDeltaDistribution(pred.obj, ncol = 3)
# all.markers <- metadata(pred.obj)$de.genes
# rds$labels <- pred.obj$labels
# scater::plotHeatmap(rds, order_columns_by="labels", features = unique(unlist(all.markers)) )

singler_out <- pred.obj$labels
names(singler_out) <- rownames(pred.obj)
obj@meta.data$singler <- singler_out
Seurat::DimPlot(obj, reduction = "umap", group.by = "singler")
Seurat::DimPlot(obj, reduction = "umap", group.by = "singler", split.by = "singler", ncol = 3)
Seurat::DimPlot(obj, reduction = "umap", group.by = "seurat_clusters")
Seurat::DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", split.by = "seurat_clusters", ncol = 3)


