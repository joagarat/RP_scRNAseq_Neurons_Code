library(scCustomize)
library(scMerge)
load("data/Smartseq2dataset_seurat_filtered.RData")
seurat_split <- SplitObject(seurat_merged, split.by = "region_label")
regions = names(seurat_split)
seurat_split <- lapply(seurat_split, function(x) {
  DefaultAssay(x) <- "RNA"
  x
})
exprs_mat = list()
segIdx = list()
param = BiocParallel::MulticoreParam(workers = 5, progressbar = TRUE)
for (i in seq(1,11)){
  exprs_mat[[i]] = as.data.frame(seurat_split[[i]]@assays$RNA$data)
  cell_type = seurat_split[[i]]@meta.data$subclass_label
  segIdx[[i]] = scSEGIndex(exprs_mat = exprs_mat[[i]], cell_type = cell_type, return_all = TRUE, BPPARAM = param)
}
for (z in seq(1,11)) {
  dir_to_create <- paste("data/SCSegIdx_smartseq_All_Regions/", as.character(regions[z]), sep = "")
  if (!dir.exists(dir_to_create)) {
    dir.create(dir_to_create, recursive = TRUE)
  }
  write.csv(segIdx[[z]], paste("data/SCSegIdx_smartseq_All_Regions/", as.character(regions[z]),"/","scSEG",".csv", sep = ""))
}