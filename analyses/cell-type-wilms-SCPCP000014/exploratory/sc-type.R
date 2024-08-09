library("openxlsx")
library("dplyr")
library("Seurat")
library("HGNChelper")

path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
path_meta <- paste0(path_proj,"/single_cell_metadata.tsv")
meta <- read.table(path_meta, sep = "\t", header = T, stringsAsFactors = F)

sample <- "SCPCS000517"; library <- "SCPCL000849"
obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Kidney" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

# extract scaled scRNA-seq matrix
gene_meta <- obj@assays[["RNA"]]@meta.data %>%
  select(c(gene_ids, gene_symbol))
scRNAseqData_scaled <- obj@assays[["RNA"]]$scale.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  left_join(gene_meta, by = c("rowname" = "gene_ids")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), rowname, gene_symbol)) %>%
  select(-c(rowname)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "gene_symbol") %>%
  as.matrix()

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. pbmc[["RNA"]]$counts or pbmc[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])