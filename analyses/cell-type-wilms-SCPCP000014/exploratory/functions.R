pre_seuratobj <- function(obj, nfeatures = 500, run_harmony = TRUE, reduction = "harmony", ndims = 50){
  ######## Normalize, scale, feature selection
  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(object = obj))
  # obj <- Seurat::SCTransform(obj)
  obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))
  # Seurat::ElbowPlot(obj, ndims = 50)
  
  ######## batch effect correction
  if (run_harmony){
    obj <- harmony::RunHarmony(obj, group.by.vars = "library_id")
  }
  
  ######## Clustering and dimentional reduction
  obj <- Seurat::FindNeighbors(obj, dims = 1:ndims, reduction = reduction)
  obj <- Seurat::FindClusters(obj, resolution = 0.8, algorithm = 1)
  obj <- Seurat::RunUMAP(obj, dims = 1:ndims, reduction = reduction)
  
  return(obj)
}