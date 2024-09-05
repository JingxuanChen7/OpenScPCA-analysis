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



