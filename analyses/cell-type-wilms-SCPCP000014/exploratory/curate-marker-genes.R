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

# obj for merged samples
obj <- SeuratObject::LoadSeuratRds(paste0("results/SCPCP000014_merged.h5Seurat"))

# from azimuth (human kidney) run a sample-only
var_genes <- list(
  list("SLC8A1","NFIB","SLIT3","FOXO1","PCDH7"),
  list("COL1A2","COL6A3","PDGFRA","ZEB2"),
  list("PAPPA2","KCTD8","TMEM182","NRG3","CUX2"),
  list("LSAMP","DCC"),
  list("NAV3")
)
var_genes <- setNames(object = var_genes, c("CNT","FIB","MD","MFIB","OMCD-PC"))

# from science 2019 supp figure
var_genes <- list(
  list("PECAM1","PLVAP","TIMP3"),
  list("SIX1","CITED1","PAX2","ALDOB","GLYAT","GPX3","SLC12A1","CLCNKA","ATP6V1B1","KRT7","S100P","UPK1A"),
  list("PTPRC","NKG7","CD3D","MS4A1","CD14","FCGR3A","CPA3","TPSAB1","ITGB3","GP6","HBM","HBZ"),
  list("LUM","PDGFRB","PRELP","TNC"),
  list("WT1","CTNNB1","AMER1","IGF2","NCAM1")
)
var_genes <- setNames(object = var_genes, c("Vasc","DevNephron","Immune","Stroma","Tumor"))

# from science 2018, table S3 # mostly not expressed
var_genes <- list(
  list("SIX2","CITED1","PAX2","SIX1"),
  list("ATP6V0D2","CLCNKB","SLC26A4","SLC4A1"),
  list("SFRP2","EMILIN1","MMP2"),
  list("PTPRO","PODXL"),
  list("CLDN16","SLC12A1"),
  list("PDGFRB","ACTA2"),
  list("HNF1B","RET","GATA3","ELF3","POU3F3","TFCP2L1","CDH16"),
  list("PLVAP","SLC14A1","VCAM1","KDR","PTPRB","PECAM1")
)
var_genes <- setNames(object = var_genes, c("CM","CD","FIB","Glomerulus","LOH","MFIB","UB","Vasc"))

# from running azimuth (fetal development) with one sample
var_genes <- list(
  list("COL1A1","COL1A2","ENOX1","PRICKLE1","RERG","FBN2","ROBO2","PLEKHH2"),
  list("VEGFA"),
  list("CNTNAP2","DPF3","BNC2","CSMD1","ITGB8","NLGN1"),
  list("ERBB4","PAX2","GFRA1","LTBP1","IGF2BP3","EDIL3","CNTN5","PTPRG","GPC3","DMD")

)
var_genes <- setNames(object = var_genes, c("Stroma","EPI","UB","MC"))

# from literature https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2020.00183/full
var_genes <- list(
  list("HOX11","OSR1",'EYA1',"PAX2","SIX1","SIX2","GDNF"),
  list("WNT4","FGF8","PAX8","LHX1","BRN1"),
  list("NOTCH1","NOTCH2"), # WNT4, LHX1
  list("KDR") # LHX1
  
)
var_genes <- setNames(object = var_genes, c("CM","PTA/RV","CSB","SSB"))

# from azimuth database (fetal)
var_genes <- list(
  list('GCSAML', 'RARB', 'ITGA2B', 'TM4SF18', 'HPSE', 'PPBP', 'MYO1B', 'RP11-556I14.2', 'CD86'), #CHRM3
  list('PRRX1', 'SCARA5', 'CTNNA2', 'SFRP2', 'EBF2', 'COL12A1', 'PPARG', 'POSTN', 'MTND4P12', 'IGFBP5'), 
  list('SRGN', 'CD74', 'NR4A3', 'SAMHD1', 'LGMN', 'FOS', 'CSF1R', 'PLAUR', 'HLA-DRA'), # MTND4P12
  list('MEIS2', 'CHRM3', 'LDB2', 'KCNIP4', 'PLVAP', 'ASIC2', 'PPAP2B', 'LINC00478', 'HSPG2'), # IGFBP5
  list('MECOM', 'NAALADL2', 'PKHD1', 'INADL', 'PAX8', 'GATA3', 'COBLL1', 'BMPR1B', 'DCDC2', 'PAX2'),
  list('KIAA1217', 'PALLD', 'MEIS1', 'FLRT2', 'TNC', 'NR2F2', 'KIF26B', 'CCBE1'), # 'MEIS2' GPC3
  list('DACH1', 'EYA1', 'HMGA2', 'GPC3', 'TRABD2B', 'ITGA8', 'PTPN14', 'WT1', 'BMPER') # PAX2
)
var_genes <- setNames(object = var_genes, c("Megakaryocytes","stroma","Myeloid","Vasc","UB","Mesangial","MetaNep"))

# from lifemap (early stage)
var_genes <- list(
  list('EYA1','FOXC1','FOXC2','GDNF','GFRA1','NCAM1','OSR1','PAX2','SALL1','SIX1',"SIX2",'WT1','HOXD11'),
  list('RET','WNT11')
)
var_genes <- setNames(object = var_genes, c("MM","WD"))

# from lifemap ( stage of UB, CM)
var_genes <- list(
  list('CITED1','SIX2','MUC1'),
  list('FOXD1','RARA','RARB'),
  list('GFRA1','RET','WNT9B') # GFRA1 both in MM and UB
)
var_genes <- setNames(object = var_genes, c("CM","InS","UB"))

# from lifemap (following stages)
var_genes <- list(
  list('AQP2','CLDN3','CLDN4','CLDN8','SCNN1B','SCNN1G'),
  list('BMP2','DKK1','DLL1','GREB1','LHX1','PAPSS2','PCSK9','POU3F3'),
  list('TMEM100','WT1'),
  list('WNT4'),
  list('FOXC2','ANO1'),
  list('CDH6','IRX3'),
  list('MAFB','NPHS1','NPHS2')
  
)
var_genes <- setNames(object = var_genes, c("CD","DRV","PRV","PAC","CBB","SBBM","SBBP"))

# interesting genes https://www.nature.com/articles/s41581-022-00598-5

# plotting
plot_obj <- sample_obj
colors <- c("blue","gray","red")
Seurat::DotPlot(plot_obj, features = var_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_gradientn(colours = colors)

scaled_genes <- GetAssayData(plot_obj, slot = "scale.data") %>% rownames()
var_genes <- unlist(var_genes)[unlist(var_genes) %in% scaled_genes ]
Seurat::DoHeatmap(plot_obj, features = var_genes, draw.lines = T ) 
Seurat::FeaturePlot(plot_obj, features = var_genes, cols = colors)
# # from cell marker
# gs_list <- by(markers$Symbol, markers$cell_name, head, n=10) 
# var_genes <-  unlist(var_genes)[unlist(var_genes) %in% scaled_genes ]



# wilms tumor markers from literature
var_genes <- c(
  "BCL6", "CCNA1", "CTHRC1", "DGKD", "EPB41L4B", "ERRFI1", "LRRC40", "NCEH1", "NEBL", "PDSS1", "ROR1", "RTKN2", # 35480093
  "TRIM28","FBXW7","NYNRIN","KDM3B", # 30885698
  "EMCN", # "CCNA1" 38937666
  "TCF3" # 34278464
  #"COL4A3","COL4A4","KCNJ1","MME","SLC12A1" # 33564352, low in tumor
)

# from science 2018, findmarkers tumor genes
markers <- read.csv(paste0("./results/aat1699-young_wilms_degenes_tumor.csv")) %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05 & pct.1 > 0.4)
var_genes <- markers$X

# plotting
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
