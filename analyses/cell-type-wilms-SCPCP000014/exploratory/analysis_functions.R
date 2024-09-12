library(Seurat)
library(dplyr)
library(ggpubr)

run_anchorTrans <- function(ref_obj, sample,
                            unknown_cutoff = 0.5) {
  # perspective output files
  filename1 <- paste0("results/plots/",sample,"_anchorTrans.pdf")
  filename2 <- paste0("results/plots/",sample,"_anchorTrans_split.pdf")
  if (file.exists(filename1) & file.exists(filename2)) {
    print("results already exist")
    return(0)
  }
  
  sample_obj <- SeuratObject::LoadSeuratRds(paste0("results/",sample,".h5Seurat"))
  # find anchors
  anchors <- FindTransferAnchors(reference = ref_obj, query = sample_obj)
  
  # clean annotation
  ref_obj@meta.data <- ref_obj@meta.data %>%
    mutate(annot = case_when(compartment == "stroma" ~ "stroma",
                              compartment == "immune" ~ "immune",
                              # grepl("S shaped body", celltype) ~ "S shaped body",
                              TRUE ~ celltype))
  ref_obj@meta.data$annot <- factor(ref_obj@meta.data$annot)
  # transfer labels
  predictions <- TransferData(
    anchorset = anchors,
    refdata = ref_obj$annot
  )
  predictions <- mutate(predictions, predicted.id = case_when(prediction.score.max < unknown_cutoff ~ "Unknown",
                                                              TRUE ~ predicted.id))
  sample_obj <- AddMetaData(object = sample_obj, metadata = predictions)
  
  # plot
  color <- Polychrome::glasbey.colors( length( unique(predictions$predicted.id) ) +1  )
  color <- color[-1]
  names(color) <- unique(predictions$predicted.id)
  color['Unknown'] <- "gray90"
  
  p1 <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", label = F, cols = color) +
    ggtitle(sample)
  p2 <- Seurat::DimPlot(sample_obj, reduction = "umap", label = T)
  df <- sample_obj@meta.data %>%
    dplyr::group_by(seurat_clusters, predicted.id) %>%
    dplyr::count(name = "sum")
  p3 <- ggplot(df, aes(x = seurat_clusters, y = sum, fill =  predicted.id)) +
    geom_bar(width = 0.5, stat = "identity", position = "fill") +
    scale_fill_manual(values = color)
  toprow <- ggpubr::ggarrange(p1, p2, widths = c(1.5,1))
  p <- ggpubr::ggarrange(toprow, p3, ncol = 1)
  p_split <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", 
                  label = F, cols = color, split.by = "predicted.id", ncol = 3, alpha = 0.1) +
    ggtitle(sample)
  
  # save plots

    ggsave(plot = p, filename = filename1, device = "pdf", 
           width = ifelse(length( unique(predictions$predicted.id) ) < 15,9,15), height = 6)
  
    ggsave(plot = p_split, filename = filename2, device = "pdf", 
           width = ifelse(length( unique(predictions$predicted.id) ) < 15,9,15), height = 9)

  
}