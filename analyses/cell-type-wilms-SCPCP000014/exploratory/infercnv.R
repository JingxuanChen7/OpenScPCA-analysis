library(Seurat)
library(dplyr)
library(infercnv)

sample <- "SCPCS000517"; library <- "SCPCL000849"
obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))

# create annotation file for infercnv
annotation_file <- paste0("results/",sample,"_annotations_file_infercnv.txt")
annotation_df <- data.frame(cell = obj$barcodes, cluster = obj$seurat_clusters)
write.table(annotation_df, file = annotation_file, sep = "\t", quote = F, row.names = F, col.names = F)
# check gene order file, from infercnv ftp
gene_order_file <- "/home/lightsail-user/wilms_tumor/ref_data/hg38_gencode_v27.txt"
gene_order <- read.table(file = gene_order_file)

# sub set genes in order file
#length(intersect(gene_order$V1, rownames(count_mat))) # 35498, quite a lot
obj <- obj[rownames(obj) %in% gene_order$V1]
count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "count")

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_mat,
                                    annotations_file=annotation_file,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=NULL) 

# took half hour?
out_dir <- paste0("results/",sample,"infercnv/")
dir.create(out_dir, showWarnings = F)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

readr::write_rds(infercnv_obj, file = paste0(out_dir,sample,"_infercnv.rds"), compress = "bz2")

var_genes <- list(
  list("PECAM1","PLVAP","TIMP3"),
  list("SIX1","CITED1","PAX2","ALDOB","GLYAT","GPX3","SLC12A1","CLCNKA","ATP6V1B1","KRT7","S100P","UPK1A"),
  list("PTPRC","NKG7","CD3D","MS4A1","CD14","FCGR3A","CPA3","TPSAB1","ITGB3","GP6","HBM","HBZ"),
  list("LUM","PDGFRB","PRELP","TNC"),
  list("WT1","CA9")
)
var_genes <- setNames(object = var_genes, c("Vasc","DevNephron","Immune","Stroma","Tumor"))
Seurat::DimPlot(obj, reduction = "umap", label = T)
Seurat::DotPlot(obj, features = var_genes, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
