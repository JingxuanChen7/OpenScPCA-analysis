library(dplyr)

path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
path_supp_table <- paste0(path_ref, "/aat1699-young-tabless1-s12-revision2.xlsx")
cell_manifest <- readxl::read_xlsx(path_supp_table, sheet = 11, skip = 1) %>%
  select(c(barcode, SangerID, DropletID, ClusterID, QCpass))
cluster_info <- readxl::read_xlsx(path_supp_table, sheet = 2, skip = 1) %>%
  select(c(Cluster_ID,Category,Cell_type1,Cell_type2,Cell_type3))
# supp object
path_supp_mtx <- paste0(path_ref, "/tableOfCounts.mtx")
path_supp_mtx_col <- paste0(path_ref, "/tableOfCounts_colLabels.tsv")
path_supp_mtx_row <- paste0(path_ref, "/tableOfCounts_rowLabels.tsv")
supp_mtx <- Matrix::readMM(path_supp_mtx)
coldata <- read.table(path_supp_mtx_col, header = T) %>% 
  left_join(cell_manifest, by = c("SangerID","DropletID","Barcode"="barcode")) %>%
  left_join(cluster_info, by = c("ClusterID" = "Cluster_ID"))
rowdata <- read.table(path_supp_mtx_row, header = T)


sce <- SummarizedExperiment::SummarizedExperiment(assays = supp_mtx, 
                                                  colData = coldata,
                                                  rowData = rowdata)

sce <- sce[, sce$QCpass == "TRUE"]
readr::write_rds(sce, paste0(path_ref, "/aat1699-young.rds"))


sceM <- scRNAseq::MuraroPancreasData()
