obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))
count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")