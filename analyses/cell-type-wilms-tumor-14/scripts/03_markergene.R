library(dplyr)
library(Seurat)
library(copykat)
library(ggplot2)
library("org.Hs.eg.db")

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

### dot plot
var_genes <- list(
  list("COL1A1","COL1A2","ENOX1","PRICKLE1","RERG","FBN2","ROBO2","PLEKHH2"),
  list("VEGFA"),
  list("CNTNAP2","DPF3","BNC2","CSMD1","ITGB8","NLGN1"),
  list("ERBB4","PAX2","GFRA1","LTBP1","IGF2BP3","EDIL3","CNTN5","PTPRG","GPC3","DMD")
  
)
var_genes <- setNames(object = var_genes, c("Stroma","EPI","UB","MC"))


var_genes <- lapply(var_genes,function(x){
  AnnotationDbi::mapIds(org.Hs.eg.db,
                        keys=as.character(x), 
                        column="ENSEMBL",
                        keytype="SYMBOL",
                        multiVals="first")
})

# plotting
plot_obj <- obj
colors <- c("blue","gray","red")
Seurat::DotPlot(plot_obj, features = var_genes, group.by = "predicted.id") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_gradientn(colours = colors)
# plots <- lapply(var_genes, function(x){
#   Seurat::VlnPlot(plot_obj, features = x, 
#                   split.by = "predicted.id", group.by = "predicted.id",
#                   stack = T, flip = T) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# })
plots <- lapply(var_genes, function(x){
  scCustomize::Stacked_VlnPlot(plot_obj, features = x, 
                               split.by = "predicted.id", group.by = "predicted.id",
                               plot_legend = T, x_lab_rotate = T) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})
ggpubr::ggarrange(plotlist = plots, ncol = 2, nrow = 3, common.legend = TRUE)

Seurat::VlnPlot(plot_obj, features = var_genes$Vasc, 
                split.by = "predicted.id", group.by = "predicted.id",
                stack = T, flip = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))