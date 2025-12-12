
library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(RColorBrewer)
library(DESeq2)
load("data/Smartseq2dataset_seurat_filtered.RData")
seurat_split <- SplitObject(seurat_merged, split.by = "region_label")
regions <- names(seurat_split)
remove(seurat_merged)
seurat_split_filtered = list()
for (i in seq(1,11)) {
  counts <- as.matrix(seurat_split[[i]]@assays$RNA@counts)
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  seurat_split_filtered[[i]] <- CreateSeuratObject(counts = counts[keep_genes, ], meta.data = seurat_split[[i]]@meta.data)
}
remove(seurat_split)
seurat_split = seurat_split_filtered
remove(seurat_split_filtered)
gc()


for (z in seq(1,11))   {
  metadata = seurat_split[[z]]@meta.data
  counts_filtered = as.data.frame(seurat_split[[z]]@assays$RNA@layers$counts)
  colnames(counts_filtered) = colnames(seurat_split[[z]])
  rownames(counts_filtered) = rownames(seurat_split[[z]])
  subclass_3 = unique(metadata$subclass_label)
  subclass_nombres = gsub("/","_", subclass_3)
  counts1 = list()
  counts2 = list()
  for(i in seq(1,length(subclass_3)-1)){
    for(h in seq(i+1,length(subclass_3))){
      rats1 = metadata[metadata$subclass_label==subclass_3[i], "donor_label"]
      rats1_stats = as.data.frame(table(as.data.frame(rats1)))
      rats1_stats_mayor_a10 =rats1_stats[rats1_stats$Freq>5,]
      rats2 = metadata[metadata$subclass_label==subclass_3[h], "donor_label"]
      rats2_stats = as.data.frame(table(as.data.frame(rats2)))
      rats2_stats_mayor_a10 =rats2_stats[rats2_stats$Freq>5,]
      if (length(rats1_stats_mayor_a10$rats1)>2 & length(rats2_stats_mayor_a10$rats2)>2) {
        for(p in seq(1,length(rats1_stats_mayor_a10$rats1))){
          counts1[[p]] = counts_filtered[,rownames(metadata[metadata$donor_label==rats1_stats_mayor_a10$rats1[p] & metadata$subclass_label==subclass_3[i],  ])]
        }
        for(p in seq(1,length(rats2_stats_mayor_a10$rats2))){
          counts2[[p]] = counts_filtered[,rownames(metadata[metadata$donor_label==rats2_stats_mayor_a10$rats2[p] & metadata$subclass_label==subclass_3[h],  ])]
        }
        DESEq_table = matrix(nrow = length(rownames(counts_filtered)), ncol = length(rats1_stats_mayor_a10$rats1)+length(rats2_stats_mayor_a10$rats2))
        nombres_columnas = 
          colnames(DESEq_table) = append(gsub(" ","",paste(subclass_3[i],"_",c(1:length(rats1_stats_mayor_a10$rats1)))), gsub(" ","",paste(subclass_3[h],"_",c(1:length(rats2_stats_mayor_a10$rats2)))))
        rownames(DESEq_table) = rownames(counts_filtered)
        for(k in seq(1,length(rats1_stats_mayor_a10$rats1))){
          DESEq_table[,k] = rowSums(counts1[[k]])
        }
        for(m in seq(1,length(rats2_stats_mayor_a10$rats2))){
          DESEq_table[,length(rats1_stats_mayor_a10$rats1)+m] = rowSums(counts2[[m]])
        }
        DESEq_table = round(DESEq_table)
        columnas = data.frame(cell = append(rep(subclass_3[i],length(rats1_stats_mayor_a10$rats1)), rep(subclass_3[h], length(rats2_stats_mayor_a10$rats2))))
        rownames(columnas) = colnames(DESEq_table)
        Deseq_object = DESeqDataSetFromMatrix(countData=DESEq_table, colData=columnas,design=~cell)
        colData(Deseq_object)$cell <- factor(colData(Deseq_object)$cell, levels=c(subclass_3[i],subclass_3[h]))
        Deseq_object = DESeq(Deseq_object)
        Deseq_res = results(Deseq_object)
        resultados = as.data.frame(Deseq_res)
        dir_to_create <- paste("data/Smartseq_DEGs/", as.character(regions[z]), sep = "")
        if (!dir.exists(dir_to_create)) {
          dir.create(dir_to_create, recursive = TRUE)
        }
        write.csv(resultados, paste("data/Smartseq_DEGs/", as.character(regions[z]),"/",as.character(subclass_nombres[i]),"vs",as.character(subclass_nombres[h]),".csv", sep = ""))
      }
    }
  }
}
