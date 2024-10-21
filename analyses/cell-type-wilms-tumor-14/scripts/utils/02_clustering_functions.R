library(reticulate)
library(sceasy)
library(SingleCellExperiment)

####################### my functions ####################### 

run_SAM <- function( path_repo, 
                     library, sample, 
                     scratch_out_dir
                     ) {
  # setup env
  use_condaenv("openscpca-cell-type-wilms-tumor-14")
  samalg <- import("samalg") #https://github.com/atarashansky/self-assembling-manifold/tree/master
  loompy <- reticulate::import('loompy')
  
  # prepare sce as 00_preprocessing_rds
  # db_proj <- file.path(path_repo,"data","current","results","doublet-detection","SCPCP000014")
  # db <- read.table(file.path(db_proj, sample, paste0(library,"_processed_scdblfinder.tsv")), header = T)
  rds <- readRDS( file.path(path_repo,"data","current","SCPCP000014", sample, paste0(library,"_processed.rds")) )
  # rds$doublet_class <- db$class
  rds <- rds[!is.na(SingleCellExperiment::rowData(rds)$gene_symbol),]
  
  sam = samalg$SAM(counts = c(r_to_py(t(counts(rds))),
                              r_to_py(as.array(rownames(rds))),
                              r_to_py(as.array(colnames(rds)))))
  sam$preprocess_data()
  # more params https://github.com/atarashansky/self-assembling-manifold/blob/master/samalg/__init__.py
  # use default here
  sam$run(distance = "correlation") 
  sam$clustering(method = "leiden") 
  sam$save_anndata( file.path(scratch_out_dir, paste0("sam_",library,".h5ad")) )
  
  # convert sam output to seurat
  # schard::h5ad2seurat doesn't work
  # potential related to this issue https://github.com/mojaveazure/seurat-disk/issues/14
  sceasy::convertFormat(file.path(scratch_out_dir, paste0("sam_",library,".h5ad")), 
                        from = "anndata", to = "seurat",
                        outFile = file.path(scratch_out_dir, paste0("sam_",library,".rds")) )
}

