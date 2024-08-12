library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(openxlsx)
library(HGNChelper)

reticulate::use_condaenv("tf")
library(cellassign)
tensorflow::tf_config()

path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)

sample <- "SCPCS000517"; library <- "SCPCL000849"
rds <- readRDS(paste0(path_proj,"/",sample,"/",library,"_processed_genesym.rds"))

##### prepare markergenes from scType db
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Kidney" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

s <- sizeFactors(rds)

marker_mat <- marker_list_to_mat(gs_list$gs_positive)
marker_mat <- marker_mat[rownames(marker_mat) %in% rownames(rds),]
input_sce <- rds[rownames(rds) %in% rownames(marker_mat),]
s <- sizeFactors(rds)
## ModuleNotFoundError: No module named 'tf_keras'
fit <- cellassign(exprs_obj = input_sce, 
                  marker_gene_info = marker_mat, 
                  s = s, 
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)
