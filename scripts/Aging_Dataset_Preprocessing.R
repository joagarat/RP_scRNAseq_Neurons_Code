library(rhdf5)
library(Seurat)
library(reticulate)
library(anndata)
library(Matrix)
library(zellkonverter)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

X_data <- h5read("data/ArchivoFiltrado_Isocortex.h5ad", "/X/data")
X_data = as.numeric(X_data)
X_indices <- h5read("data/ArchivoFiltrado_Isocortex.h5ad", "/X/indices")
X_indptr <- h5read("data/ArchivoFiltrado_Isocortex.h5ad", "/X/indptr")
X_atts <- h5readAttributes("data/ArchivoFiltrado_Isocortex.h5ad", "/X")
n_cells <- X_atts$shape[1]
n_genes <- X_atts$shape[2]
X_sparse <- sparseMatrix(
  i = rep(seq_len(n_cells), diff(X_indptr)),
  j = X_indices + 1, 
  x = X_data,
  dims = c(n_cells, n_genes)
)

metadata <- list()
genes <- h5read("data/ArchivoFiltrado_Isocortex.h5ad", "/var/_index")
cell_barcodes <- h5read("data/ArchivoFiltrado_Isocortex.h5ad", "/obs/sample_id")
obs_items <- h5ls("data/ArchivoFiltrado_Isocortex.h5ad", recursive = TRUE)
obs_items <- obs_items[grepl("^/obs/", obs_items$group), ]
categories_rows <- obs_items$name == "categories"
codes_rows <- obs_items$name == "codes"
categories_map <- setNames(obs_items$group[categories_rows], obs_items$group[categories_rows])
codes_map <- setNames(obs_items$group[codes_rows], obs_items$group[codes_rows])
metadata <- list()

for (var_path in unique(codes_map)) {
  cat_path <- paste0(var_path, "/categories")
  codes_path <- paste0(var_path, "/codes")
  var_name <- sub("^/obs/", "", var_path) 
  categories <- h5read("data/ArchivoFiltrado_Isocortex.h5ad", cat_path)
  codes <- h5read("data/ArchivoFiltrado_Isocortex.h5ad", codes_path)
  metadata[[var_name]] <- factor(categories[codes + 1])
}

metadata_df <- as.data.frame(metadata)
head(metadata_df)

cell_barcodes = as.character(cell_barcodes)
genes =as.character(genes)
rownames(metadata_df) <- cell_barcodes
rownames(X_sparse) = cell_barcodes
colnames(X_sparse) = genes
seurat <- CreateSeuratObject(counts = t(X_sparse), meta.data = metadata_df)
remove(X_sparse, metadata, genes, cell_barcodes, codes, X_data, X_atts, X_indices, X_indptr)
seurat = NormalizeData(seurat)
seurat= FindVariableFeatures(seurat)
seurat= ScaleData(seurat)
seurat = RunPCA(seurat)
seurat = RunUMAP(seurat, dims = 1:50)
DimPlot(seurat, group.by = "subclass_label")
seurat_adult <- subset(seurat, subset = age_cat == "adult")
seurat_aged <- subset(seurat, subset = age_cat == "aged")
saveRDS(seurat_adult, "data/Normalized_Seurat_adult")
saveRDS(seurat_aged, "data/Normalized_Seurat_aged")
remove(seurat_aged, seurat)
seurat = seurat_adult
remove(seurat_adult)
excluir <- c("MOL NN_4", "CA2-FC-IG Glut_2", "OPC NN_1", "COP NN_1", "MFOL NN_3", "NFOL NN_2")
metadata_cols <- colnames(seurat@meta.data)
columna_a_usar <- NULL
columnas_candidatas <- c("supertype_label") 
for (col in columnas_candidatas) {
  if (col %in% metadata_cols && any(seurat@meta.data[[col]] %in% excluir)) {
    columna_a_usar <- col
    break 
  }
}
cells_to_keep <- rownames(seurat@meta.data)[!(seurat@meta.data[[columna_a_usar]] %in% excluir)]
seurat_filtrado <- subset(seurat, cells = cells_to_keep)
seurat = seurat_filtrado
remove(seurat_filtrado)
table(seurat$roi)
samples = h5read("data/Chromium_expression_matrix.hdf5", "/data/samples")
genes = h5read("data/Chromium_expression_matrix.hdf5", "/data/gene")
metadata = read.csv("data/CTX_Hip_anno_10x.csv")
rownames(metadata) = metadata$exp_component_name
total_regions = unique(as.character(metadata$region_label))
metadata = subset(metadata, subset = region_label %in% c("PL-ILA-ORB", "ACA", "AI", "RSP"))
metadata = subset(metadata, subset = class_label %in% c("Glutamatergic", "GABAergic"))
table(metadata$region_label)
target_counts <- c(
  "ACA" = 12447,
  "AI" = 11976,
  "PL-ILA-ORB" = 10535,
  "RSP" = 12570
)

metadata_subset <- purrr::map_dfr(names(target_counts), function(region_name) {
  n_to_sample <- target_counts[region_name]
  region_subset <- metadata %>%
    filter(region_label == region_name)
  if (nrow(region_subset) < n_to_sample) {
    warning(paste("¡Aviso! Para la región", region_name, "se pidieron", n_to_sample, 
                  "células, pero solo hay", nrow(region_subset), ". Se tomarán todas las disponibles."))
    n_to_sample <- nrow(region_subset)
  }
  region_subset %>%
    slice_sample(n = n_to_sample, replace = FALSE)
})

metadata = metadata_subset
remove(metadata_subset)
index = which(samples%in%metadata$exp_component_name)
mtx = h5read("data/Chromium_expression_matrix.hdf5", name = "data/counts", index = list(index, NULL))
mtx = as.data.frame(mtx)
mtx = t(mtx)
samplenames = samples[index]
colnames(mtx) = samplenames
rownames(mtx) = genes
mtx = as.data.frame(mtx)
seurat2 = CreateSeuratObject(counts = mtx, meta.data = metadata)
remove(mtx, index, samples, genes)
gc()
seurat2 = NormalizeData(seurat2)
set.seed(42)  
split_cells <- sample(Cells(seurat), size = floor(length(Cells(seurat)) / 2))
seurat_part1 <- subset(seurat, cells = split_cells)
seurat_part2 <- subset(seurat, cells = setdiff(Cells(seurat), split_cells))
common_features <- SelectIntegrationFeatures(object.list = list(seurat2, seurat))
length(common_features)
remove(seurat)
gc()
anchors1 <- FindTransferAnchors(
  reference = seurat2,   
  query = seurat_part1,    
  dims = 1:50,    
  features = common_features
)
anchors2 <- FindTransferAnchors(
  reference = seurat2,   
  query = seurat_part2,    
  dims = 1:50,                  
  features = common_features
)
predictions1 <- TransferData(anchors = anchors1, refdata = seurat2$subclass_label)
predictions2 <- TransferData(anchors = anchors2, refdata = seurat2$subclass_label)

predictions_class_1 <- TransferData(anchorset = anchors1, refdata = seurat2$class_label,  
                                    dims = 1:50
)
predictions_class_2 <- TransferData(anchorset = anchors2, refdata = seurat2$class_label,  
                                    dims = 1:50
)

seurat_part1$predicted_cell_class <- predictions_class_1$predicted.id
seurat_part1$predicted_cell_class_score_max = predictions_class_1$prediction.score.max
seurat_part2$predicted_cell_class <- predictions_class_2$predicted.id
seurat_part2$predicted_cell_class_score_max = predictions_class_2$prediction.score.max
seurat_part1 <- AddMetaData(seurat_part1, metadata = predictions1)
seurat_part2 <- AddMetaData(seurat_part2, metadata = predictions2)
seurat_combined <- merge(seurat_part1, seurat_part2)
VlnPlot(seurat_combined, features = "predicted_cell_class_score_max")
seurat_0_95_confidence_filtered <- subset(seurat_combined, subset = predicted_cell_class_score_max > 0.95)
dim(seurat_combined)
dim(seurat_0_95_confidence_filtered)
remove(seurat_part1, seurat_part2, seurat_combined,seurat)
gc()
DefaultAssay(seurat_0_95_confidence_filtered) <- "RNA"
layer1 <- seurat_0_95_confidence_filtered@assays$RNA@layers$counts.1
layer2 <- seurat_0_95_confidence_filtered@assays$RNA@layers$counts.2
layer_combined = cbind(layer1,layer2)
colnames(layer_combined) = colnames(seurat_0_95_confidence_filtered)
rownames(layer_combined) = rownames(seurat_0_95_confidence_filtered)
metadata = as.data.frame(seurat_0_95_confidence_filtered@meta.data)

seurat_combined = CreateSeuratObject(counts = layer_combined, meta.data = metadata)
remove(seurat_0_95_confidence_filtered, layer_combined, layer1, layer2)
gc()
seurat_combined = NormalizeData(seurat_combined)
seurat_combined= FindVariableFeatures(seurat_combined)
seurat_combined= ScaleData(seurat_combined)
seurat_combined = RunPCA(seurat_combined)
seurat_combined = RunUMAP(seurat_combined, dims = 1:50)

saveRDS(seurat_combined, "data/Normalized_Classified_Seurat_Adult")
remove(anchors1, anchors2, seurat_combined)
remove(seurat_combined)
gc()
seurat = readRDS("data/Normalized_Seurat_Aged")
excluir <- c("MOL NN_4", "CA2-FC-IG Glut_2", "OPC NN_1", "COP NN_1", "MFOL NN_3", "NFOL NN_2")
metadata_cols <- colnames(seurat@meta.data)
columna_a_usar <- NULL
columnas_candidatas <- c("supertype_label") 
for (col in columnas_candidatas) {
  if (col %in% metadata_cols && any(seurat@meta.data[[col]] %in% excluir)) {
    columna_a_usar <- col
    break 
  }
}
if (is.null(columna_a_usar)) {
  stop("ERROR: Ninguna de las columnas de metadatos candidatas (predicted.id, subclass_label, etc.) contiene los tipos celulares a excluir. Por favor, verifica los nombres en tu lista 'excluir'.")
}
cells_to_keep <- rownames(seurat@meta.data)[!(seurat@meta.data[[columna_a_usar]] %in% excluir)]

seurat_filtrado <- subset(seurat, cells = cells_to_keep)
seurat = seurat_filtrado
remove(seurat_filtrado)
table(seurat$roi)
set.seed(42)  #
split_cells <- sample(Cells(seurat), size = floor(length(Cells(seurat)) / 2))

seurat_part1 <- subset(seurat, cells = split_cells)
seurat_part2 <- subset(seurat, cells = setdiff(Cells(seurat), split_cells))

gc()
common_features <- SelectIntegrationFeatures(object.list = list(seurat2, seurat))
length(common_features)
remove(seurat)
gc()
anchors1 <- FindTransferAnchors(
  reference = seurat2,   
  query = seurat_part1,   
  dims = 1:50,                  
  features = common_features
)
anchors2 <- FindTransferAnchors(
  reference = seurat2,   
  query = seurat_part2,  
  dims = 1:50,        
  features = common_features
)
predictions1 <- TransferData(anchors = anchors1, refdata = seurat2$subclass_label)
predictions2 <- TransferData(anchors = anchors2, refdata = seurat2$subclass_label)

predictions_class_1 <- TransferData(anchorset = anchors1, refdata = seurat2$class_label,  
                                    dims = 1:50
)
predictions_class_2 <- TransferData(anchorset = anchors2, refdata = seurat2$class_label, 
                                    dims = 1:50
)

seurat_part1$predicted_cell_class <- predictions_class_1$predicted.id
seurat_part1$predicted_cell_class_score_max = predictions_class_1$prediction.score.max
seurat_part2$predicted_cell_class <- predictions_class_2$predicted.id
seurat_part2$predicted_cell_class_score_max = predictions_class_2$prediction.score.max
seurat_part1 <- AddMetaData(seurat_part1, metadata = predictions1)
seurat_part2 <- AddMetaData(seurat_part2, metadata = predictions2)
seurat_combined <- merge(seurat_part1, seurat_part2)
VlnPlot(seurat_combined, features = "predicted_cell_class_score_max")
seurat_0_95_confidence_filtered <- subset(seurat_combined, subset = predicted_cell_class_score_max > 0.95)
dim(seurat_combined)
dim(seurat_0_95_confidence_filtered)
dim(seurat)
remove(seurat_part1, seurat_part2, seurat_combined,seurat)
gc()
DefaultAssay(seurat_0_95_confidence_filtered) <- "RNA"
layer1 <- seurat_0_95_confidence_filtered@assays$RNA@layers$counts.1
layer2 <- seurat_0_95_confidence_filtered@assays$RNA@layers$counts.2
layer_combined = cbind(layer1,layer2)
colnames(layer_combined) = colnames(seurat_0_95_confidence_filtered)
rownames(layer_combined) = rownames(seurat_0_95_confidence_filtered)
metadata = as.data.frame(seurat_0_95_confidence_filtered@meta.data)
seurat_combined = CreateSeuratObject(counts = layer_combined, meta.data = metadata)
remove(seurat_0_95_confidence_filtered, layer_combined, layer1, layer2)
gc()
seurat_combined = NormalizeData(seurat_combined)
seurat_combined= FindVariableFeatures(seurat_combined)
seurat_combined= ScaleData(seurat_combined)
seurat_combined = RunPCA(seurat_combined)
seurat_combined = RunUMAP(seurat_combined, dims = 1:50)
DimPlot(seurat_combined, group.by = "predicted.id")
remove(anchors1, anchors2)
saveRDS(seurat_combined, "data/Normalized_Classified_Seurat_Aged")
saveRDS(seurat2, "data/Normalized_Reference_10x_PlILAORB_ACA_AI_RSP")
gc()

