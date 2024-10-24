library(dplyr)
library(Seurat)

path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
path_meta <- file.path(path_repo,"data","current","SCPCP000014","single_cell_metadata.tsv") # keep for debug
# path_meta <- file.path(opt$metadata)
meta <- read.table(path_meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

# import clustering functions from main repo
source(file = file.path(path_repo,"packages","rOpenScPCA","R","calculate-clusters.R"))
source(file = file.path(path_anal,"scripts","utils","02_clustering_functions.R"))
source(file = file.path(path_anal,"scripts","utils","01_anchor_transfer_seurat_functions.R"))

# library <- "SCPCL000850"; sample <- "SCPCS000518"

## output files
scratch_out_dir <- file.path(path_anal, "scratch", "02_clustering")
dir.create(scratch_out_dir, showWarnings = FALSE, recursive = TRUE)
results_out_dir <- file.path(path_anal, "results", "02_clustering")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)
plots_out_dir <- file.path(path_anal, "plots", "02_clustering")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)

purrr::walk2(
  meta$scpca_sample_id,
  meta$scpca_library_id,
  \(sample, library) run_SAM( path_repo = path_repo,
                              sample = sample, library = library,
                              scratch_out_dir = scratch_out_dir
                              )
)


# obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library,".rdsSeurat")) )
# # DimPlot(obj, group.by = "cluster")
# # cluster_df <- calculate_clusters(obj,
# #                                  algorithm = "leiden",
# #                                  resolution = 0.006, objective_function = "CPM",
# #                                  seed = 233)
# cluster_df <- calculate_clusters(obj, 
#                                  algorithm = "leiden", 
#                                  resolution = 1, objective_function = "modularity",
#                                  seed = 233)
# 
# 
# obj <- AddMetaData(object = obj, metadata = cluster_df)
# DimPlot(obj, group.by = "cluster")


# # SCTransform
# obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library,".rdsSeurat")) )
# obj <- Seurat::SCTransform(obj)
# obj <- Seurat::RunPCA(obj)
# obj <- Seurat::FindNeighbors(obj, dims = 1:ndims, reduction = reduction)
# cluster_df <- calculate_clusters(obj, 
#                                  algorithm = "leiden", 
#                                  resolution = 1, objective_function = "modularity",
#                                  seed = 233)
# obj <- AddMetaData(object = obj, metadata = cluster_df)
# obj <- Seurat::SetIdent(obj, value = "cluster")
# obj <- Seurat::RunUMAP(obj, dims = 1:ndims, reduction = reduction)
# DimPlot(obj, group.by = "cluster")

# SAM algorithm

# plt <- import("matplotlib.pyplot")
# sam$scatter()
# plt$show()
# obj <- readRDS( file.path(scratch_out_dir, paste0("sam_",library,".rds")) )
# DimPlot(obj, group.by = "leiden_clusters")
# 
# cluster_df <- calculate_clusters(obj, 
#                                  algorithm = "leiden", 
#                                  resolution = 1, objective_function = "modularity",
#                                  seed = 233)
# 
# 
# obj <- AddMetaData(object = obj, metadata = cluster_df)
# DimPlot(obj, group.by = "cluster")
# SeuratDisk::Convert(file.path(scratch_out_dir, paste0("sam_",library,".h5ad")) , dest = "h5seurat", overwrite = T)
# sam_obj <- SeuratDisk::LoadH5Seurat( file.path(scratch_out_dir, paste0("sam_",library,".h5seurat")) )
#sam_obj <- schard::h5ad2seurat( filename = file.path(scratch_out_dir, paste0("sam_",library,".h5ad")) )
# need this because of schard error (even convert to h5Seurat) 
# potential related to this issue https://github.com/mojaveazure/seurat-disk/issues/141
#sam_sce <- zellkonverter::readH5AD( file.path(scratch_out_dir, paste0("sam_",library,".h5ad")) )




# anchor transfer results
level = "celltype"
predictions <- read.csv(file.path(path_anal, "results", "01_anchor_transfer_seurat", "RNA", paste0(library, "_", level,".csv")))
obj <- AddMetaData(object = obj, metadata = predictions)
plots <- plot_anchorTrans(path_anal = path_anal, 
                 sample_obj = obj,
                 library = library, 
                 level = level,
                 plot_cluster = "leiden_clusters",
                 internal = TRUE)

