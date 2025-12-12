library(rhdf5)
samples = h5read("data/Chromium_expression_matrix.hdf5", "/data/samples")
genes = h5read("data/Chromium_expression_matrix.hdf5", "/data/gene")
metadata = read.csv("data/CTX_Hip_anno_10x.csv")
rownames(metadata) = metadata$exp_component_name
total_regions = unique(as.character(metadata$region_label))
metadata = subset(metadata, subset = class_label %in% c("Glutamatergic", "GABAergic"))
regions_pass_DESeq2 = c("VISp", "VIS", "SSs", "SSp", "RSPv", "RSP", "MOp", "HIP", "ALM", "AI", "ACA")
metadata = subset(metadata, subset = region_label %in% regions_pass_DESeq2)
unique(as.character(metadata$class_label))
df = as.data.frame(table(metadata$cluster_label))
df = df[df$Freq>10, ]
metadata = subset(metadata, subset = cluster_label %in% df$Var1)
subclass = unique(metadata$subclass_label)

metadata_list = split(metadata, f = metadata$region_label)
#Filter >10cells/cluster 
for (i in seq(1,8)) {
  df =as.data.frame(table(metadata_list[[i]]$cluster_label))
  df = df[df$Freq>10, ]
  metadata_list[[i]] = subset(metadata_list[[i]], subset = cluster_label %in% df$Var1)
}
remove(df)
regions = array()
for (i in seq(1,8)) {
  regions = append(regions, unique(metadata_list[[i]]$region_label))
}
regions = regions[-1]


index = list()
for (i in seq(1,8)) {
  index[[i]] = which(samples%in%metadata_list[[i]]$exp_component_name)
  samplenames = samples[index[[i]]]
  mtx = h5read("data_/Chromium_expression_matrix.hdf5", name = "data/counts", index = list(index[[i]], NULL))
  mtx = as.data.frame(mtx)
  mtx = t(mtx)
  colnames(mtx) = samplenames
  rownames(mtx) = genes
  library(Seurat)
  library(ggplot2)
  RP = read.csv("data/RP.csv")
  RP = as.character(RP$x)
  library(DESeq2)
  subclass_3 = unique(metadata_list[[i]]$subclass_label)
  subclass_nombres = gsub("/","_", subclass_3)
  counts1 = list()
  counts2 = list()
  metadata = metadata_list[[i]]
  metadata = as.data.frame(metadata)
  for(j in seq(1,length(subclass_3)-1)){
    for(h in seq(j+1,length(subclass_3))){
      rats1 = metadata[metadata$subclass_label==subclass_3[j], "donor_label"]
      rats1_stats = as.data.frame(table(as.data.frame(rats1)))
      rats1_stats_mayor_a10 =rats1_stats[rats1_stats$Freq>5,]
      rats2 = metadata[metadata$subclass_label==subclass_3[h], "donor_label"]
      rats2_stats = as.data.frame(table(as.data.frame(rats2)))
      rats2_stats_mayor_a10 =rats2_stats[rats2_stats$Freq>5,]
      if (length(rats1_stats_mayor_a10$rats1)>2 & length(rats2_stats_mayor_a10$rats2)>2) {
        for(p in seq(1,length(rats1_stats_mayor_a10$rats1))){
          counts1[[p]] = mtx[,rownames(metadata[metadata$donor_label==rats1_stats_mayor_a10$rats1[p] & metadata$subclass_label==subclass_3[j],  ])]
        }
        for(p in seq(1,length(rats2_stats_mayor_a10$rats2))){
          counts2[[p]] = mtx[,rownames(metadata[metadata$donor_label==rats2_stats_mayor_a10$rats2[p] & metadata$subclass_label==subclass_3[h],  ])]
        }
        DESEq_table = matrix(nrow = length(rownames(mtx)), ncol = length(rats1_stats_mayor_a10$rats1)+length(rats2_stats_mayor_a10$rats2))
        nombres_columnas = 
          colnames(DESEq_table) = append(gsub(" ","",paste(subclass_3[j],"_",c(1:length(rats1_stats_mayor_a10$rats1)))), gsub(" ","",paste(subclass_3[h],"_",c(1:length(rats2_stats_mayor_a10$rats2)))))
        rownames(DESEq_table) = rownames(mtx)
        for(k in seq(1,length(rats1_stats_mayor_a10$rats1))){
          DESEq_table[,k] = rowSums(counts1[[k]])
        }
        for(m in seq(1,length(rats2_stats_mayor_a10$rats2))){
          DESEq_table[,length(rats1_stats_mayor_a10$rats1)+m] = rowSums(counts2[[m]])
        }
        DESEq_table = round(DESEq_table)
        DESEq_table = as.data.frame(DESEq_table)
        columnas = data.frame(cell = append(rep(subclass_3[j],length(rats1_stats_mayor_a10$rats1)), rep(subclass_3[h], length(rats2_stats_mayor_a10$rats2))))
        rownames(columnas) = colnames(DESEq_table)
        Deseq_object = DESeqDataSetFromMatrix(countData=DESEq_table, colData=columnas, design=~cell)
        colData(Deseq_object)$cell <- factor(colData(Deseq_object)$cell, levels=c(subclass_3[j],subclass_3[h]))
        Deseq_object = DESeq(Deseq_object)
        Deseq_res = results(Deseq_object)
        resultados = as.data.frame(Deseq_res)
        dir_to_create <- paste("data/Chromium_DEGs/", as.character(regions[i]), sep = "")
        if (!dir.exists(dir_to_create)) {
          dir.create(dir_to_create, recursive = TRUE)
        }
        write.csv(resultados, paste("data/Chromium_DEGs/", as.character(regions[i]),"/",as.character(subclass_nombres[j]),"vs",as.character(subclass_nombres[h]),".csv", sep = ""))
      }
    }
  }
}



#Analysis for PL-ILA-ORB (same region as Stress Dataset)

samples <- h5read("data/Chromium_expression_matrix.hdf5", "/data/samples")
genes <- h5read("data/Chromium_expression_matrix.hdf5", "/data/gene")
metadata <- read.csv("data/CTX_Hip_anno_10x.csv")
rownames(metadata) <- metadata$exp_component_name
metadata <- subset(metadata, class_label %in% c("Glutamatergic", "GABAergic"))
target_region <- "PL-ILA-ORB"
metadata_region <- subset(metadata, region_label == target_region)
if(nrow(metadata_region) == 0){
  stop("No hay celdas anotadas como PL-ILA-ORB en metadata.")
}
df <- as.data.frame(table(metadata_region$cluster_label))
df <- df[df$Freq > 10, ]
metadata_region <- subset(metadata_region, cluster_label %in% df$Var1)
subclass <- unique(metadata_region$subclass_label)
idx <- which(samples %in% metadata_region$exp_component_name)
sample_names <- samples[idx]
mtx <- h5read("data/Chromium_expression_matrix.hdf5", name = "data/counts",
              index = list(idx, NULL))
mtx <- t(as.data.frame(mtx))
colnames(mtx) <- sample_names
rownames(mtx) <- genes
subclass_list <- unique(metadata_region$subclass_label)
subclass_names_clean <- gsub("/", "_", subclass_list)
out_dir <- paste0("data/Chromium_DEGs_PL_ILA_ORB/", target_region)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

metadata_r <- metadata_region 

for(j in seq_len(length(subclass_list)-1)){
  for(h in seq(j+1, length(subclass_list))){
    rats1 <- metadata_r[metadata_r$subclass_label == subclass_list[j], "donor_label"]
    rats2 <- metadata_r[metadata_r$subclass_label == subclass_list[h], "donor_label"]
    rats1_tab <- as.data.frame(table(rats1))
    rats2_tab <- as.data.frame(table(rats2))
    rats1_tab <- rats1_tab[rats1_tab$Freq > 5, ]
    rats2_tab <- rats2_tab[rats2_tab$Freq > 5, ]
    if(nrow(rats1_tab) < 2 | nrow(rats2_tab) < 2){
      next
    }
    counts1 <- list()
    counts2 <- list()
    for(p in seq_len(nrow(rats1_tab))){
      donor <- rats1_tab$rats1[p]
      cells <- rownames(metadata_r[metadata_r$donor_label == donor &
                                     metadata_r$subclass_label == subclass_list[j], ])
      counts1[[p]] <- mtx[, cells, drop = FALSE]
    }
    for(p in seq_len(nrow(rats2_tab))){
      donor <- rats2_tab$rats2[p]
      cells <- rownames(metadata_r[metadata_r$donor_label == donor &
                                     metadata_r$subclass_label == subclass_list[h], ])
      counts2[[p]] <- mtx[, cells, drop = FALSE]
    }
    DE_mtx <- matrix(
      nrow = nrow(mtx),
      ncol = nrow(rats1_tab) + nrow(rats2_tab)
    )
    rownames(DE_mtx) <- rownames(mtx)
    col_count <- 1
    for(k in seq_len(nrow(rats1_tab))){
      DE_mtx[, col_count] <- rowSums(counts1[[k]])
      col_count <- col_count + 1
    }
    for(k in seq_len(nrow(rats2_tab))){
      DE_mtx[, col_count] <- rowSums(counts2[[k]])
      col_count <- col_count + 1
    }
    
    DE_mtx <- round(DE_mtx)
    col_df <- data.frame(
      cell = c(
        rep(subclass_list[j], nrow(rats1_tab)),
        rep(subclass_list[h], nrow(rats2_tab))
      )
    )
    rownames(col_df) <- paste0("sample", seq_len(ncol(DE_mtx)))
    dds <- DESeqDataSetFromMatrix(DE_mtx, colData = col_df, design = ~ cell)
    dds$cell <- factor(dds$cell, levels = c(subclass_list[j], subclass_list[h]))
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds))
    out_path <- paste0(
      out_dir, "/",
      subclass_names_clean[j], "_vs_", subclass_names_clean[h], ".csv"
    )
    write.csv(res, out_path)
  }
}
