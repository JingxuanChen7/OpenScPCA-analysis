library(dplyr)
library(Seurat)

path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
path_meta <- file.path(path_repo,"data","current","SCPCP000014","single_cell_metadata.tsv") # keep for debug
# path_meta <- file.path(opt$metadata)
meta <- read.table(path_meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

source(file = file.path(path_anal,"scripts","utils","02_clustering_functions.R"))
source(file = file.path(path_anal,"scripts","utils","01_anchor_transfer_seurat_functions.R"))

library <- "SCPCL000850"; sample <- "SCPCS000518"

obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library,".rdsSeurat")) )
# DimPlot(obj, group.by = "cluster")

## output files
scratch_out_dir <- file.path(path_anal, "scratch", "02_clustering")
dir.create(scratch_out_dir, showWarnings = FALSE, recursive = TRUE)
results_out_dir <- file.path(path_anal, "results", "02_clustering")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)
plots_out_dir <- file.path(path_anal, "plots", "02_clustering")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)

# cluster_df <- calculate_clusters(obj,
#                                  algorithm = "leiden",
#                                  resolution = 0.006, objective_function = "CPM",
#                                  seed = 233)
cluster_df <- calculate_clusters(obj, 
                                 algorithm = "leiden", 
                                 resolution = 1, objective_function = "modularity",
                                 seed = 233)


obj <- AddMetaData(object = obj, metadata = cluster_df)
DimPlot(obj, group.by = "cluster")

obj <- Seurat::SetIdent(obj, value = "cluster")
ndims <- 50; reduction <- "pca"
obj <- Seurat::RunUMAP(obj, dims = 1:ndims, reduction = reduction)
DimPlot(obj, group.by = "cluster")

# SCTransform
obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library,".rdsSeurat")) )
obj <- Seurat::SCTransform(obj)
obj <- Seurat::RunPCA(obj)
obj <- Seurat::FindNeighbors(obj, dims = 1:ndims, reduction = reduction)
cluster_df <- calculate_clusters(obj, 
                                 algorithm = "leiden", 
                                 resolution = 1, objective_function = "modularity",
                                 seed = 233)
obj <- AddMetaData(object = obj, metadata = cluster_df)
obj <- Seurat::SetIdent(obj, value = "cluster")
obj <- Seurat::RunUMAP(obj, dims = 1:ndims, reduction = reduction)
DimPlot(obj, group.by = "cluster")

# SAM algorithm
library(reticulate)
use_condaenv("openscpca-cell-type-wilms-tumor-14")
samalg <- import("samalg") #https://github.com/atarashansky/self-assembling-manifold/tree/master

rds <- readRDS( file.path(path_repo,"data","current","SCPCP000014", sample, paste0(library,"_processed.rds")) )
rds <- rds[!is.na(SingleCellExperiment::rowData(rds)$gene_symbol),]

sam = samalg$SAM(counts = c(r_to_py(t(counts(rds))),
                            r_to_py(as.array(rownames(rds))),
                            r_to_py(as.array(colnames(rds)))))

# count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "count")
# count_mat <- t(count_mat)
# sam = samalg$SAM(counts = c(r_to_py(count_mat),
#                             r_to_py(as.array(rownames(obj))),
#                             r_to_py(as.array(colnames(obj)))))
sam$preprocess_data()
sam$run(distance = 'correlation')
sam$clustering(method = "leiden") #leiden clustering is the default algorithm in SAM
scipy <- import("scipy")
sam$adata = scipy$sparse$csr_matrix(sam$adata)
sam$save_anndata( file.path(scratch_out_dir, paste0("sam_",library,".h5ad")) )
library(sceasy)
loompy <- reticulate::import('loompy')
sceasy::convertFormat(file.path(scratch_out_dir, paste0("sam_",library,".h5ad")), 
                      from = "anndata", to = "seurat",
                      outFile = file.path(scratch_out_dir, paste0("sam_",library,".rds")) )
# plt <- import("matplotlib.pyplot")
# sam$scatter()
# plt$show()
obj <- readRDS( file.path(scratch_out_dir, paste0("sam_",library,".rds")) )
DimPlot(obj, group.by = "leiden_clusters")

cluster_df <- calculate_clusters(obj, 
                                 algorithm = "leiden", 
                                 resolution = 1, objective_function = "modularity",
                                 seed = 233)


obj <- AddMetaData(object = obj, metadata = cluster_df)
DimPlot(obj, group.by = "cluster")
# SeuratDisk::Convert(file.path(scratch_out_dir, paste0("sam_",library,".h5ad")) , dest = "h5seurat", overwrite = T)
# sam_obj <- SeuratDisk::LoadH5Seurat( file.path(scratch_out_dir, paste0("sam_",library,".h5seurat")) )
sam_obj <- schard::h5ad2seurat( filename = file.path(scratch_out_dir, paste0("sam_",library,".h5ad")) )
# need this because of schard error (even convert to h5Seurat) 
# potential related to this issue https://github.com/mojaveazure/seurat-disk/issues/141
sam_sce <- zellkonverter::readH5AD( file.path(scratch_out_dir, paste0("sam_",library,".h5ad")) )




# anchor transfer results
level = "celltype"
predictions <- read.csv(file.path(path_anal, "results", "01_anchor_transfer_seurat", paste0(library, "_", level,".csv")))
obj <- AddMetaData(object = obj, metadata = predictions)
plot_anchorTrans(path_anal = path_anal, 
                 sample_obj = obj,
                 library = library, 
                 level = level,
                 plot_cluster = "cluster",
                 internal = TRUE)

