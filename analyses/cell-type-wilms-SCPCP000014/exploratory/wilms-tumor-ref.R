library(dplyr)

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
rowdata <- read.table(path_supp_mtx_row, header = T)


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


########################################################################
path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
sce <- zellkonverter::readH5AD(paste0(path_ref, "/Fetal_full_v3.h5ad"))

