library(dplyr)
library(Seurat)
library(infercnv)
library(ggplot2)

path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
path_meta <- file.path(path_repo,"data","current","SCPCP000014","single_cell_metadata.tsv") # keep for debug
# path_meta <- file.path(opt$metadata)
meta <- read.table(path_meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

# create output dirs
scratch_out_dir <- file.path(path_anal, "scratch", "03_cnv")
dir.create(scratch_out_dir, showWarnings = FALSE, recursive = TRUE)
results_out_dir <- file.path(path_anal, "results", "03_cnv")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)
plots_out_dir <- file.path(path_anal, "plots", "03_cnv")
dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)

# load pre-processed object and add anchor transfer results
library <- "SCPCL000850"
obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library,".rdsSeurat")) )
level = "compartment"
predictions <- read.csv(file.path(path_anal, "results", "01_anchor_transfer_seurat", paste0(library, "_", level,".csv")))
obj <- AddMetaData(object = obj, metadata = predictions)

out_dir <- file.path(scratch_out_dir, library)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# create annotation file for infercnv
annotation_file <- file.path(scratch_out_dir, paste0(library,"_annotations_file_infercnv.txt"))
# annotation_df <- data.frame(cell = obj$barcodes, cluster = obj$cluster)
annotation_df <- data.frame(cell = obj$barcodes, cluster = obj$predicted.id)
write.table(annotation_df, file = annotation_file, sep = "\t", quote = F, row.names = F, col.names = F)
# check gene order file, from infercnv ftp
gene_order_file <- file.path(scratch_out_dir,"Homo_sapiens.GRCh38.104.gene_order.txt")
gene_order <- read.table(file = gene_order_file) %>%
  arrange(V2, V3)

obj <- obj[rownames(obj) %in% gene_order$V1]
count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "count")

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_mat,
                                    annotations_file=annotation_file,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names="stroma") 


infercnv_out = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             denoise=TRUE,
                             HMM=TRUE,
                             save_rds = F,
                             num_threads = 8)

readr::write_rds(infercnv_out, file = file.path(out_dir,"infercnv_out.rds"))



infercnv_out <- readr::read_rds(file = file.path(out_dir,"run.final.infercnv_obj"))



## summary strategy 1 https://www.biostars.org/p/9573777/
scores <- apply(infercnv_out@expr.data, 2 ,function(x){ sum(x < 0.95 | x > 1.05)/length(x) }) %>% 
  as.data.frame() %>%
  mutate(pred = obj$predicted.id) %>%
  mutate(cnv_score = ifelse(
    pred == "stroma",
    NA,
    .
  ))
# not bimodal
hist(scores$cnv_score, breaks = 100)
obj <- AddMetaData(object = obj, metadata = scores)
p1 <- FeaturePlot(obj, features = "cnv_score", alpha = 0.3) +
  scale_color_viridis_c()

obj_cnv <- infercnv::add_to_seurat(
  seurat_obj = obj,
  infercnv_output_path = out_dir
)

## summary strategy 2 ewings analysis, based on number of chr that have cnv
cnv_df <- obj_cnv@meta.data %>%
  select(matches("predicted.id") | starts_with("has_cnv_chr") & !matches("has_cnv_chrMT")) %>%
  mutate(count_cnv_chr = ifelse(predicted.id == "stroma",
                            NA,
                            rowSums(across(starts_with("has_cnv_chr")))
                            ) )

# not bimodal
hist(cnv_df$count_cnv_chr)

obj <- AddMetaData(object = obj, metadata = cnv_df)
p2 <- FeaturePlot(obj, features = "count_cnv_chr", alpha = 0.3) +
  scale_color_viridis_c()

## summary strategy 3 ewings analysis, based on number of chr that have cn
chr_weights <- infercnv_out@gene_order |> 
  as.data.frame() |> 
  dplyr::count(chr) |> 
  # only keep chr 1-22, drops MT and GL chr
  dplyr::filter(chr %in% glue::glue("chr{seq(1,22)}")) |> 
  dplyr::pull(n)

cnv_df <- obj_cnv@meta.data %>%
  select(matches("predicted.id") | starts_with("proportion_scaled_cnv_chr") & !ends_with("chrMT")) %>%
  rowwise() %>%
  mutate(weight_mean = ifelse(
    predicted.id == "stroma",
    NA,
    weighted.mean(across(starts_with("proportion_scaled_cnv_chr")), chr_weights/sum(chr_weights))
  ) )

obj <- AddMetaData(object = obj, metadata = cnv_df)
p3 <- FeaturePlot(obj, features = "weight_mean", alpha = 0.3) +
  scale_color_viridis_c()

hist(cnv_df$weight_mean)

## summarizing three evaluation methods and look into correlations
sum_df <- obj@meta.data %>%
  select(c("cnv_score","count_cnv_chr","weight_mean"))
cors <- cor(sum_df, use = "pairwise.complete.obs")
p4 <- ggcorrplot::ggcorrplot(cors)

ggpubr::ggarrange(p1, p2, p3, p4,
                  ncol = 2, nrow = 2)
