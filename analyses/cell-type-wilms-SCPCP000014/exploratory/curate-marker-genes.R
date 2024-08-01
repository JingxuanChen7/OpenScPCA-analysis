library(dplyr)

path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
path_supp_table <- paste0(path_ref, "/aat1699-young-tabless1-s12-revision2.xlsx")
markers <- readxl::read_xlsx(path_supp_table, sheet = 4, skip = 1) %>%
  filter(isKeyGene == 1)