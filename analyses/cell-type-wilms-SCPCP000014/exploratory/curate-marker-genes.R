library(dplyr)
library(Seurat)
library(ggpubr)
library(RColorBrewer)

####################### kidney atlas ####################### 

path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
# sce <- zellkonverter::readH5AD(paste0(path_ref, "/Fetal_full_v3.h5ad"))
# seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(sce),
#                                                assay = "RNA",
#                                                project = "kidneyatlas")
# # convert colData and rowData to data.frame for use in the Seurat object
# cell_metadata <- as.data.frame(SingleCellExperiment::colData(sce))
# row_metadata <- as.data.frame(SingleCellExperiment::rowData(sce))
# # add cell metadata (colData) from SingleCellExperiment to Seurat
# seurat_obj@meta.data <- cell_metadata
# # add row metadata (rowData) from SingleCellExperiment to Seurat
# seurat_obj[["RNA"]]@meta.data <- row_metadata
# # add metadata from SingleCellExperiment to Seurat
# seurat_obj@misc <- S4Vectors::metadata(sce)
# # make a copy for processing
# obj <- seurat_obj
# # log transform counts
# obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
obj <- SeuratObject::LoadSeuratRds(paste0("results/kidneyatlas.h5Seurat"))
# calculate marker genes for compartments
Idents(obj) <- obj@meta.data[["compartment"]]
markers_all <- FindAllMarkers(object = obj, 
                       only.pos = TRUE,
                       logfc.threshold = 0.25) 
write.csv(markers_all, file = paste0("./results/degenes_compart.csv"))
# calculate marker genes for all celltypes
Idents(obj) <- obj@meta.data[["celltype"]]
markers_celltype <- FindAllMarkers(object = obj, 
                              only.pos = TRUE,
                              logfc.threshold = 0.25) 
write.csv(markers_celltype, file = paste0("./results/degenes_celltype.csv"))

####################### cell markers ####################### 
# get human marker list from http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html
path_ref <- "/home/lightsail-user/wilms_tumor/ref_data"
path_supp_table <- paste0(path_ref, "/Cell_marker_Human.xlsx")
markers <- readxl::read_xlsx(path_supp_table) %>%
  filter(tissue_class == "Kidney")
  

####################### manual gene collection ####################### 
path_proj <- "/home/lightsail-user/wilms_tumor/OpenScPCA-analysis/data/current/SCPCP000014"
sample <- "SCPCS000517"; library <- "SCPCL000849"
sample_obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))

# from azimuth sample-only
var_genes <- list(
  list("SLC8A1","NFIB","SLIT3","FOXO1","PCDH7"),
  list("COL1A2","COL6A3","PDGFRA","ZEB2"),
  list("PAPPA2","KCTD8","TMEM182","NRG3","CUX2"),
  list("LSAMP","DCC"),
  list("NAV3")
)
var_genes <- setNames(object = var_genes, c("CNT","FIB","MD","MFIB","OMCD-PC"))
# # from merged cohort (not working)
# var_genes <- list(
#   list("MECOM","NFIA","FMN1","SLC26A7","PAX8"),
#   list("RBPMS","C7","ZFPM2-AS1","COL4A1","THBS1")
# 
# )
# var_genes <- setNames(object = var_genes, c("CNT","FIB"))
# from literature
var_genes <- c(
               "BCL6", "CCNA1", "CTHRC1", "DGKD", "EPB41L4B", "ERRFI1", "LRRC40", "NCEH1", "NEBL", "PDSS1", "ROR1", "RTKN2", # 35480093
               "TRIM28","FBXW7","NYNRIN","KDM3B", # 30885698
               "EMCN", # "CCNA1" 38937666
               "TCF3" # 34278464
               #"COL4A3","COL4A4","KCNJ1","MME","SLC12A1" # 33564352, low in tumor
               )


Seurat::DotPlot(sample_obj, features = var_genes, cols = c("blue","red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
scaled_genes <- GetAssayData(sample_obj, slot = "scale.data") %>% rownames()
var_genes <- unlist(var_genes)[unlist(var_genes) %in% scaled_genes ]
Seurat::DoHeatmap(sample_obj, features = var_genes, draw.lines = T ) 
Seurat::FeaturePlot(sample_obj, features = var_genes)
# # from cell marker
# gs_list <- by(markers$Symbol, markers$cell_name, head, n=10) 
# var_genes <-  unlist(var_genes)[unlist(var_genes) %in% scaled_genes ]

# from science 2019 supp figure
var_genes <- list(
  list("PECAM1","PLVAP","TIMP3"),
  list("SIX1","CITED1","PAX2","ALDOB","GLYAT","GPX3","SLC12A1","CLCNKA","ATP6V1B1","KRT7","S100P","UPK1A"),
  list("PTPRC","NKG7","CD3D","MS4A1","CD14","FCGR3A","CPA3","TPSAB1","ITGB3","GP6","HBM","HBZ"),
  list("LUM","PDGFRB","PRELP","TNC"),
  list("WT1","CTNNB1","AMER1","IGF2","NCAM1")
)
var_genes <- setNames(object = var_genes, c("Vasc","DevNephron","Immune","Stroma","Tumor"))

# from science 2018, findmarkers tumor genes
markers <- read.csv(paste0("./results/aat1699-young_wilms_degenes_tumor.csv")) %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05 & pct.1 > 0.4)
var_genes <- markers$X

colors <- c("blue","gray","red")

sample_obj_plot <- Seurat::AddModuleScore(sample_obj, features = list(var_genes),name = "tumor")
Seurat::FeaturePlot(sample_obj_plot, features = "tumor1") +
  scale_colour_gradientn(colours = colors)
p1 <- Seurat::DotPlot(sample_obj_plot, features = "tumor1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_gradientn(colours = colors)# + ggtitle("avg_modulescore")
p2 <- Seurat::DotPlot(sample_obj_plot, features = var_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradientn(colours = colors) 
ggarrange(p1, p2, widths = c(1,5), common.legend = T, legend = "right")
# Seurat::DoHeatmap(sample_obj, features = var_genes, draw.lines = T, slot = "data" ) 

scaled_genes <- GetAssayData(sample_obj, slot = "scale.data") %>% rownames()
var_genes <- unlist(var_genes)[unlist(var_genes) %in% scaled_genes ]
# Seurat::DoHeatmap(sample_obj, features = var_genes, draw.lines = T ) 
sample_obj_plot <- Seurat::AddModuleScore(sample_obj, features = list(var_genes),name = "tumor")
Seurat::FeaturePlot(sample_obj_plot, features = "tumor1") +
  scale_colour_gradientn(colours = colors)
p1 <- Seurat::DotPlot(sample_obj_plot, features = "tumor1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_gradientn(colours = colors)# + ggtitle("avg_modulescore")
p2 <- Seurat::DotPlot(sample_obj_plot, features = var_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradientn(colours = colors) 
ggarrange(p1, p2, widths = c(1,5), common.legend = T, legend = "right")
