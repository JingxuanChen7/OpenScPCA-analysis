#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(reticulate)
library(sceasy)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--metadata"),
    type = "character",
    default = NULL,
    help = "Path to cohort metadata"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## general paths
path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
# path_meta <- file.path(path_repo,"data","current","SCPCP000014","single_cell_metadata.tsv") # keep for debug
path_meta <- file.path(opt$metadata)
meta <- read.table(path_meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

## output paths
scratch_out_dir <- file.path(path_anal, "scratch", "02_clustering")
dir.create(scratch_out_dir, showWarnings = FALSE, recursive = TRUE)

################### Functions ################### 

run_SAM <- function( path_repo, 
                     library, sample, 
                     scratch_out_dir
) {
  # output intermediate files
  out_h5ad <- file.path(scratch_out_dir, paste0("sam_",library,".h5ad"))
  out_rds <- file.path(scratch_out_dir, paste0("sam_",library,".rds"))
  
  # setup env
  use_condaenv("openscpca-cell-type-wilms-tumor-14")
  samalg <- import("samalg") #https://github.com/atarashansky/self-assembling-manifold/tree/master
  loompy <- reticulate::import('loompy')
  
  # prepare sce as in 00_preprocessing_rds.R
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
  sam$clustering(method = "louvain") 
  sam$save_anndata( out_h5ad )
  
  # convert sam output to seurat
  # schard::h5ad2seurat doesn't work
  # potential related to this issue https://github.com/mojaveazure/seurat-disk/issues/14
  sceasy::convertFormat( out_h5ad, 
                        from = "anndata", to = "seurat",
                        outFile = out_rds )
}

################### Run SAM with all samples ################### 

purrr::walk2(
  meta$scpca_sample_id,
  meta$scpca_library_id,
  \(sample, library) run_SAM( path_repo = path_repo,
                              sample = sample, library = library,
                              scratch_out_dir = scratch_out_dir
  )
)
