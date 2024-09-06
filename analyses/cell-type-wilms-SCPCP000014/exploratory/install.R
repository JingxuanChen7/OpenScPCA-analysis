#sudo apt install libglpk40
install.packages("igraph")
install.packages("Seurat")
install.packages("tidyverse")
install.packages("HGNChelper")
install.packages("openxlsx")
install.packages("ggpubr")
install.packages("Polychrome")

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
BiocManager::install("zellkonverter")
BiocManager::install("scMerge")
BiocManager::install("harmony")
#sudo apt-get install -y libmagick++-dev
BiocManager::install("SpatialExperiment")
BiocManager::install("singleCellTK")
BiocManager::install('glmGamPoi')

install.packages("viridis")
install.packages("pheatmap")


# install tensorflow & cellassign
install.packages("tensorflow")
install.packages("devtools")
library(reticulate)
py_version <- "3.10"
env_name <- "tf"
reticulate::conda_create(envname = env_name, python_version = py_version)
tensorflow::install_tensorflow(envname = env_name, 
                               method = "conda",
                               extra_packages = "tensorflow-probability",
                               new_env = F)
reticulate::use_condaenv(env_name)
library(tensorflow)
tensorflow::tf_config()
reticulate::py_config()
devtools::install_github("JingxuanChen7/cellassign")
# pip install tf-keras


devtools::install_github('immunogenomics/presto')

devtools::install_github("wjawaid/enrichR")
