path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)

sample <- "SCPCS000514"; library <- "SCPCL000846"
rds <- readRDS(paste0(path_proj,"/",sample,"/",library,"_processed.rds"))

######## create seurat object from the SCE counts matrix
seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts(rds),
                                 assay = "RNA",
                                 project = library)
# convert colData and rowData to data.frame for use in the Seurat object
cell_metadata <- as.data.frame(colData(rds))
row_metadata <- as.data.frame(rowData(rds))
# add cell metadata (colData) from SingleCellExperiment to Seurat
seurat_obj@meta.data <- cell_metadata
# add row metadata (rowData) from SingleCellExperiment to Seurat
seurat_obj[["RNA"]]@meta.data <- row_metadata
# add metadata from SingleCellExperiment to Seurat
seurat_obj@misc <- metadata(rds)
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

######## automated cell annotation - singleR
celldex::surveyReferences()
blueprint_ref <- celldex::BlueprintEncodeData(ensembl = TRUE)
hpca_ref <- celldex::HumanPrimaryCellAtlasData(ensembl = TRUE)
count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
#count_mat@Dimnames[[1]] <- seurat_obj[["RNA"]]@meta.data$gene_symbol
pred.obj <- SingleR::SingleR(test = count_mat, 
                    ref = blueprint_ref, 
                    assay.type.test=1,
                    labels = blueprint_ref$label.main)
SingleR::plotScoreHeatmap(pred.obj)
SingleR::plotDeltaDistribution(pred.obj, ncol = 3)
# add singleR results to seurat obj
singler_out <- pred.obj$labels
names(singler_out) <- rownames(pred.obj)
obj@meta.data$singler <- singler_out
Seurat::DimPlot(obj, reduction = "umap", group.by = "singler")
Seurat::DimPlot(obj, reduction = "umap", group.by = "singler", split.by = "singler", ncol = 3)
