library(dplyr)
library(Seurat)
library(copykat)
library(ggplot2)

path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
path_meta <- file.path(path_repo,"data","current","SCPCP000014","single_cell_metadata.tsv") # keep for debug
# path_meta <- file.path(opt$metadata)
meta <- read.table(path_meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

# create output dirs
scratch_out_dir <- file.path(path_anal, "scratch", "03_cnv")
dir.create(scratch_out_dir, showWarnings = FALSE, recursive = TRUE)
results_out_dir <- file.path(path_anal, "results", "03_cnv")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)
plots_out_dir <- file.path(path_anal, "plots", "03_cnv")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)

library <- "SCPCL000850"
obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library,".rdsSeurat")) )
level = "compartment"
predictions <- read.csv(file.path(path_anal, "results", "01_anchor_transfer_seurat", paste0(library, "_", level,".csv"))) 
obj <- AddMetaData(object = obj, metadata = predictions)

# save copykat results
library_out_dir <- file.path(scratch_out_dir, library)
setwd(library_out_dir)

# run copykat
count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "count")
normal_cells <- predictions %>%
  tibble::column_to_rownames(var = "X") %>%
  filter(predicted.id == "stroma")
copykat_result <- copykat(
  rawmat = count_mat,
  id.type = "E",
  sam.name = library,
  norm.cell.names = rownames(normal_cells),
  ngene.chr = 2,
  plot.genes = FALSE,
  output.seg = FALSE,
  n.cores = 8
)

copykat_df <- copykat_result$prediction %>%
  select(copykat.pred) %>% 
  rename(c("stroma_ref.copykat" = "copykat.pred"))
table(copykat_df)
obj <- AddMetaData(object = obj, metadata = copykat_df)
DimPlot(obj, group.by = "stroma_ref.copykat", alpha = 0.3)

copykat_result_noref <- copykat(
  rawmat = count_mat,
  id.type = "E",
  sam.name = paste0(library,"_noref"),
  norm.cell.names = NULL,
  ngene.chr = 2,
  plot.genes = FALSE,
  output.seg = FALSE,
  n.cores = 8
)
copykat_noref_df <- copykat_result_noref$prediction %>%
  select(copykat.pred) %>% 
  rename(c("noref.copykat" = "copykat.pred"))
table(copykat_noref_df)
obj <- AddMetaData(object = obj, metadata = copykat_noref_df)
DimPlot(obj, group.by = "noref.copykat", alpha = 0.3)


