#sudo apt install libglpk40
install.packages("igraph")
install.packages("Seurat")
install.packages("tidyverse")

BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SingleR")
#sudo apt-get install libcurl4-openssl-dev
BiocManager::install("ensembldb")
BiocManager::install("scRNAseq")
BiocManager::install("scuttle")
BiocManager::install("scater")
#sudo apt install jags r-cran-rjags
BiocManager::install("infercnv")

install.packages("viridis")
install.packages("pheatmap")

# install cellassign
install.packages("tensorflow")
install.packages("devtools")
reticulate::install_python("3.10") 
reticulate::virtualenv_create("r-tensorflow", version = "3.10")
reticulate::use_virtualenv("r-tensorflow")
reticulate::use_python("/home/lightsail-user/.virtualenvs/r-tensorflow/bin/python")
tensorflow::install_tensorflow(method = "virtualenv",
                               envname = "r-tensorflow", 
                               extra_packages = "tensorflow-probability",
                               python_version = "3.10")
devtools::install_github("Irrationone/cellassign")
Sys.setenv(RETICULATE_PYTHON="/home/lightsail-user/.virtualenvs/r-tensorflow/bin/python")
library(tensorflow)
tensorflow::tf_config()
install_tensorflow()
