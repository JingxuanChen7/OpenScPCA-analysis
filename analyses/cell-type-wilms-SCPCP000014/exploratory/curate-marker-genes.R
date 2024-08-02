library(dplyr)

####################### wilms tumor markers - Young et al. ####################### 
path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
path_supp_table <- paste0(path_ref, "/aat1699-young-tabless1-s12-revision2.xlsx")
markers <- readxl::read_xlsx(path_supp_table, sheet = 4, skip = 1) %>%
  filter(isKeyGene == 1)
cell_manifest <- readxl::read_xlsx(path_supp_table, sheet = 11, skip = 1) %>%
  select(c(barcode, SangerID, DropletID, ClusterID, QCpass, Source))
cluster_info <- readxl::read_xlsx(path_supp_table, sheet = 2, skip = 1) %>%
  select(c(Cluster_ID,Category,Cell_type1,Cell_type2,Cell_type3))

filter_markers <- left_join(markers, cluster_info, by = c("Cluster" = "Cluster_ID") ) %>%
  filter(!is.na(Category) & !is.na(Cell_type1)) %>%
  filter(grepl("tumour", Category))
