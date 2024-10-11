library(dplyr)
library(Seurat)

path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
path_meta <- file.path(path_repo,"data","current","SCPCP000014","single_cell_metadata.tsv") # keep for debug
# path_meta <- file.path(opt$metadata)
meta <- read.table(path_meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

source(file = file.path(path_anal,"scripts","utils","02_clustering_functions.R"))
source(file = file.path(path_anal,"scripts","utils","01_anchor_transfer_seurat_functions.R"))

library <- "SCPCL000850"

obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library,".rdsSeurat")) )
DimPlot(obj, group.by = "cluster")

cluster_df <- calculate_clusters(obj,
                                 algorithm = "leiden",
                                 resolution = 0.006, objective_function = "CPM",
                                 seed = 233)
cluster_df <- calculate_clusters(obj, 
                                 algorithm = "leiden", 
                                 resolution = 1, objective_function = "modularity",
                                 seed = 233)
cluster_df <- calculate_clusters(obj, 
                                 algorithm = "leiden", 
                                 resolution = 1, objective_function = "modularity",
                                 seed = 233)

obj <- AddMetaData(object = obj, metadata = cluster_df)
DimPlot(obj, group.by = "cluster")

# anchor transfer results
level = "compartment"
predictions <- read.csv(file.path(path_anal, "results", "01_anchor_transfer_seurat", paste0(library, "_", level,".csv")))
obj <- AddMetaData(object = obj, metadata = predictions)
plot_anchorTrans(path_anal = path_anal, 
                 sample_obj = obj,
                 library = library, 
                 level = level,
                 plot_cluster = "cluster",
                 internal = TRUE)

