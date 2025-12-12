library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(RColorBrewer)
library(svglite)
library(dittoSeq)
library(rhdf5)  

load("data/Seurat.ss.rda")
meta_nueva <- read.csv("data/CTX_Hip_anno_SSv4.csv")

seurat_2 <- subset(ss.seurat, cells = meta_nueva$exp_component_name)
rm(ss.seurat)

rownames(meta_nueva) <- meta_nueva$exp_component_name
meta_nueva <- meta_nueva[colnames(seurat_2), ]
seurat_2 <- AddMetaData(seurat_2, metadata = meta_nueva)

seurat_split <- SplitObject(seurat_2, split.by = "region_label")
rm(seurat_2); gc()

for (i in seq_along(seurat_split)) {
  seurat_split[[i]] <- subset(seurat_split[[i]], class_label %in% c("Glutamatergic", "GABAergic"))
}

for (i in seq_along(seurat_split)) {
  df <- as.data.frame(table(seurat_split[[i]]$cluster_label))
  df <- df[df$Freq > 10, ]
  seurat_split[[i]] <- subset(seurat_split[[i]], cluster_label %in% df$Var1)
}

regions_pass_DESeq2 <- c("VISp","VIS","SSs","SSp","RSPv","RSP","MOp","HIP","ALM","AI","ACA")

for (i in seq_along(seurat_split)) {
  names(seurat_split)[i] <- seurat_split[[i]]$region_label[1]
}

seurat_split <- seurat_split[regions_pass_DESeq2]

subclass_pass <- list()
for (z in seq_along(seurat_split)) {
  metadata <- seurat_split[[z]]@meta.data
  subclass_3 <- unique(metadata$subclass_label)
  keep_subs <- character(0)
  
  for (sub in subclass_3) {
    rats1 <- metadata[metadata$subclass_label == sub, "donor_label"]
    rats1_stats <- as.data.frame(table(rats1))
    rats1_stats_mayor_a5 <- rats1_stats[rats1_stats$Freq > 5, ]
    if (nrow(rats1_stats_mayor_a5) > 2) {
      keep_subs <- append(keep_subs, sub)
    }
  }
  
  subclass_pass[[z]] <- keep_subs
  if (length(keep_subs) > 0) {
    seurat_split[[z]] <- subset(seurat_split[[z]], subclass_label %in% keep_subs)
  } else {
    seurat_split[[z]] <- NULL
  }
}

seurat_split <- seurat_split[!sapply(seurat_split, is.null)]
subclass_pass <- unique(unlist(subclass_pass))
regions <- names(seurat_split)

for (i in names(seurat_split)) {
  seurat_split[[i]] <- RenameCells(seurat_split[[i]], add.cell.id = i)
}
seurat_merged <- merge(x = seurat_split[[1]], y = seurat_split[2:length(seurat_split)], add.cell.ids = regions)

seurat_merged <- NormalizeData(seurat_merged)
RP <- read.csv("data/RP.csv")
RP <- gsub("Rack1", "Gnb2l1", as.character(RP$x))

referende_data <- as.data.frame(seurat_merged@assays$RNA$data)
ribosomal_data <- referende_data[RP, ]

cell_ribosomal_sum <- Matrix::colSums(ribosomal_data)
seurat_merged <- AddMetaData(seurat_merged, metadata = cell_ribosomal_sum, col.name = "Global RP expression")

datos_modificados <- apply(ribosomal_data, 2, function(x) ifelse(x > 0, 1, 0))
suma_columnas <- colSums(datos_modificados)
seurat_merged <- AddMetaData(seurat_merged, metadata = suma_columnas, col.name = "Detected_RP")

data_counts <- FetchData(seurat_merged, vars = c("Detected_RP", "subclass_label")) %>%
  group_by(subclass_label, Detected_RP) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(subclass_label) %>%
  mutate(total = sum(count), freq = count / total)

seurat_merged <- SCTransform(seurat_merged)
seurat_merged <- RunPCA(seurat_merged)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:50)

save(seurat_merged, file = "data/Smartseq2dataset_seurat_filtered.RData")


#10x

rm(seurat_merged, seurat_2, seurat_split); gc()


h5ls("data/Chromium_expression_matrix.hdf5")
samples <- h5read("data/Chromium_expression_matrix.hdf5", "/data/samples")
genes   <- h5read("data/Chromium_expression_matrix.hdf5", "/data/gene")
metadata <- read.csv("data/CTX_Hip_anno_10x.csv")
rownames(metadata) <- metadata$exp_component_name
metadata <- subset(metadata, class_label %in% c("Glutamatergic", "GABAergic"))
regions_pass_DESeq2 <- c("VISp","VIS","SSs","SSp","RSPv","RSP","MOp","HIP","ALM","AI","ACA")
metadata <- subset(metadata, region_label %in% regions_pass_DESeq2)

df <- as.data.frame(table(metadata$cluster_label))
df <- df[df$Freq > 10, ]
metadata <- subset(metadata, cluster_label %in% df$Var1)

subset_metadata <- metadata[metadata$subclass_label %in% subclass_pass, ]
nmice10x <- unique(subset_metadata$donor_label)
nclusters10x <- length(unique(subset_metadata$cluster_label))

set.seed(123)
n_target <- min(100000, nrow(subset_metadata))
metadata_subsampled <- subset_metadata[sample(nrow(subset_metadata), n_target), ]

index <- which(samples %in% metadata_subsampled$exp_component_name)
mtx   <- h5read("data/Chromium_expression_matrix.hdf5", name = "data/counts", index = list(index, NULL))

mtx <- t(mtx)
colnames(mtx) <- samples[index]
rownames(mtx) <- genes
mtx = as.data.frame(mtx)

metadata_subsampled <- metadata_subsampled[colnames(mtx), ]
stopifnot(all(rownames(metadata_subsampled) == colnames(mtx)))
metadata_subsampled <- metadata_subsampled[colnames(mtx), ]
stopifnot(nrow(metadata_subsampled) == ncol(mtx))
stopifnot(all(rownames(metadata_subsampled) == colnames(mtx)))
seurat_merged <- CreateSeuratObject(counts = mtx, meta.data = metadata_subsampled)

seurat_merged <- NormalizeData(seurat_merged)

RP <- read.csv("data/RP.csv")
RP = as.character(RP$x)

referende_data <- as.data.frame(seurat_merged@assays$RNA$data)
ribosomal_data <- referende_data[RP, ]

cell_ribosomal_sum <- Matrix::colSums(ribosomal_data)
seurat_merged <- AddMetaData(seurat_merged, metadata = cell_ribosomal_sum, col.name = "RibosomalSum")

datos_modificados <- apply(ribosomal_data, 2, function(x) ifelse(x > 0, 1, 0))
suma_columnas <- colSums(datos_modificados)
seurat_merged <- AddMetaData(seurat_merged, metadata = suma_columnas, col.name = "Detected_RP")

seurat_merged <- SCTransform(seurat_merged)
seurat_merged <- RunPCA(seurat_merged)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:50)

save(seurat_merged, file = "data/seurat_merged_filtered_final.RData")

rm(seurat_merged, mtx, metadata_subsampled, subset_metadata); gc()
