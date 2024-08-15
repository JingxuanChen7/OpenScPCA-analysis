library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(openxlsx)
library(HGNChelper)
library(ggplot2)

reticulate::use_condaenv("tf")
library(cellassign)
tensorflow::tf_config()

path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)

sample <- "SCPCS000517"; library <- "SCPCL000849"
rds <- readRDS(paste0(path_proj,"/",sample,"/",library,"_processed_genesym.rds"))
obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))

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
gs_list$gs_positive[[2]]
marker_mat <- marker_list_to_mat(gs_list$gs_positive)

##### prepare markergenes from kidney atlas


marker_mat <- marker_mat[rownames(marker_mat) %in% rownames(rds),]
input_sce <- rds[rownames(rds) %in% rownames(marker_mat),]

counts <- counts(input_sce)

# remove low expression genes
input_sce <- input_sce[-which(rowSums(counts) <= 100)]
marker_mat <- marker_mat[rownames(marker_mat) %in% rownames(input_sce),]
marker_mat <- marker_mat[,-which(colSums(marker_mat)==0)]
marker_mat <- marker_mat[rownames(input_sce),]

s <- sizeFactors(input_sce)
fit <- cellassign(exprs_obj = input_sce, 
                  marker_gene_info = marker_mat, 
                  s = s, 
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)

table(celltypes(fit))
pheatmap::pheatmap(cellprobs(fit))
# add annotation
annoobj <- obj
annoobj@meta.data$cellassign <- factor(celltypes(fit))
color <- Polychrome::glasbey.colors( length(unique(celltypes(fit))) +1  )
color <- color[-1]
names(color) <- unique(celltypes(fit))
color['unassigned'] <- "gray90"

                       
p1 <- Seurat::DimPlot(annoobj, reduction = "umap", group.by = "cellassign", cols=color) + 
  theme(legend.position = "bottom") + 
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3))) 
p2 <- Seurat::DimPlot(annoobj, reduction = "umap", group.by = "seurat_clusters") + 
  theme(legend.position = "bottom")
ggpubr::ggarrange(p1, p2, ncol = 2)
