library(Seurat)
library(DESeq2)
library(dplyr)
library(stringr)
library(data.table)

seurat_adult <- readRDS("data/Normalized_Classified_Seurat_Adult")
seurat_aged <- readRDS("data/Normalized_Classified_Seurat_Aged")
meta_paper <- fread("data/metadata_aging.csv")
run_pairwise_deseq <- function(seurat_object, meta_paper_df, group_name, base_output_path, min_donors = 3) {
  output_dir <- file.path(base_output_path, paste0("DESeq_", group_name, "_Between_Subclasses"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  mtx <- as.data.frame(seurat_object@assays$RNA$counts)
  metadata <- seurat_object@meta.data
  metadata$library_prep <- sapply(rownames(metadata), function(x) strsplit(x, "-")[[1]][2])
  metadata$barcode <- rownames(metadata)
  merged_metadata <- merge(metadata, meta_paper_df, by = "library_prep", all.x = TRUE)
  rownames(merged_metadata) <- merged_metadata$barcode
  excluir <- c("MOL NN_4", "CA2-FC-IG Glut_2", "OPC NN_1", "COP NN_1", "MFOL NN_3", "NFOL NN_2")
  metadata_final <- merged_metadata[!merged_metadata$supertype_label %in% excluir, ]
  mtx_final <- mtx[, rownames(metadata_final)]
  subclasses_to_compare <- unique(as.character(metadata_final$predicted.id))
  for (j in 1:(length(subclasses_to_compare) - 1)) {
    for (h in (j + 1):length(subclasses_to_compare)) {
      subclass1_name <- subclasses_to_compare[j]
      subclass2_name <- subclasses_to_compare[h]
      meta_sub1 <- metadata_final[metadata_final$predicted.id == subclass1_name, ]
      meta_sub2 <- metadata_final[metadata_final$predicted.id == subclass2_name, ]
      donors1 <- as.data.frame(table(meta_sub1$external_donor_name)) %>% filter(Freq > 5)
      donors2 <- as.data.frame(table(meta_sub2$external_donor_name)) %>% filter(Freq > 5)
      if (nrow(donors1) < min_donors || nrow(donors2) < min_donors) {
        next 
      }
      tryCatch({
        counts1 <- sapply(donors1$Var1, function(donor) rowSums(mtx_final[, rownames(meta_sub1[meta_sub1$external_donor_name == donor, ])]))
        counts2 <- sapply(donors2$Var1, function(donor) rowSums(mtx_final[, rownames(meta_sub2[meta_sub2$external_donor_name == donor, ])]))
        deseq_counts <- cbind(counts1, counts2)
        deseq_metadata <- data.frame(
          condition = c(rep(subclass1_name, ncol(counts1)), rep(subclass2_name, ncol(counts2)))
        )
        rownames(deseq_metadata) <- colnames(deseq_counts)
        dds <- DESeqDataSetFromMatrix(
          countData = round(deseq_counts),
          colData = deseq_metadata,
          design = ~ condition
        )
        dds <- DESeq(dds) 
        res <- results(dds, contrast = c("condition", subclass1_name, subclass2_name))
        res_df <- as.data.frame(res)
        subclass1_clean <- gsub("/", "_", subclass1_name)
        subclass2_clean <- gsub("/", "_", subclass2_name)
        filename <- paste0(subclass1_clean, "vs", subclass2_clean, ".csv")
        write.csv(res_df, file.path(output_dir, filename))
      }, error = function(e) {
      }) 
    }
  }
}

base_output_path <- "data/Aging/"

run_pairwise_deseq(
  seurat_object = seurat_adult, 
  meta_paper_df = meta_paper,
  group_name = "Adult", 
  base_output_path = base_output_path,
  min_donors = 3 
)
run_pairwise_deseq(
  seurat_object = seurat_aged, 
  meta_paper_df = meta_paper,
  group_name = "Aged", 
  base_output_path = base_output_path,
  min_donors = 3
)

#DEGs_AdultvsAged


seurat_adult <- readRDS("data/Normalized_Classified_Seurat_Adult")
seurat_aged <- readRDS("data/Normalized_Classified_Seurat_Aged")
meta_paper <- fread("data/metadata_aging.csv")

run_aging_deseq_per_subclass <- function(seurat_adult, seurat_aged, meta_paper_df, base_output_path, min_donors = 3) {
  output_dir <- file.path(base_output_path, "DESeq_Adult_vs_Aged_per_Subclass")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  mtx_adult <- as.data.frame(seurat_adult@assays$RNA$counts)
  metadata_adult <- seurat_adult@meta.data
  metadata_adult$library_prep <- sapply(rownames(metadata_adult), function(x) strsplit(x, "-")[[1]][2])
  metadata_adult$barcode <- rownames(metadata_adult)
  merged_metadata_adult <- merge(metadata_adult, meta_paper_df, by = "library_prep", all.x = TRUE)
  rownames(merged_metadata_adult) <- merged_metadata_adult$barcode
  excluir <- c("MOL NN_4", "CA2-FC-IG Glut_2", "OPC NN_1", "COP NN_1", "MFOL NN_3", "NFOL NN_2")
  metadata_final_adult <- merged_metadata_adult[!merged_metadata_adult$supertype_label %in% excluir, ]
  mtx_final_adult <- mtx_adult[, rownames(metadata_final_adult)]
  mtx_aged <- as.data.frame(seurat_aged@assays$RNA$counts)
  metadata_aged <- seurat_aged@meta.data
  metadata_aged$library_prep <- sapply(rownames(metadata_aged), function(x) strsplit(x, "-")[[1]][2])
  metadata_aged$barcode <- rownames(metadata_aged)
  merged_metadata_aged <- merge(metadata_aged, meta_paper_df, by = "library_prep", all.x = TRUE)
  rownames(merged_metadata_aged) <- merged_metadata_aged$barcode
  metadata_final_aged <- merged_metadata_aged[!merged_metadata_aged$supertype_label %in% excluir, ]
  mtx_final_aged <- mtx_aged[, rownames(metadata_final_aged)]
  remove(seurat_aged,seurat_adult)
  gc()
  subclasses_comunes <- intersect(
    unique(as.character(metadata_final_adult$predicted.id)),
    unique(as.character(metadata_final_aged$predicted.id))
  )
  for (subclase in subclasses_comunes) {
    meta_sub_adult <- metadata_final_adult[metadata_final_adult$predicted.id == subclase, ]
    meta_sub_aged <- metadata_final_aged[metadata_final_aged$predicted.id == subclase, ]
    donors_adult <- as.data.frame(table(meta_sub_adult$external_donor_name)) %>% filter(Freq > 5)
    donors_aged <- as.data.frame(table(meta_sub_aged$external_donor_name)) %>% filter(Freq > 5)
    if (nrow(donors_adult) < min_donors || nrow(donors_aged) < min_donors) {
      cat(paste("    -> Omitido: no hay suficientes donantes (Adult:", nrow(donors_adult), ", Aged:", nrow(donors_aged), ").\n"))
      next
    }
    tryCatch({
      counts_adult <- sapply(donors_adult$Var1, function(d) rowSums(mtx_final_adult[, rownames(meta_sub_adult[meta_sub_adult$external_donor_name == d, ])]))
      counts_aged <- sapply(donors_aged$Var1, function(d) rowSums(mtx_final_aged[, rownames(meta_sub_aged[meta_sub_aged$external_donor_name == d, ])]))
      deseq_counts <- cbind(counts_adult, counts_aged)
      deseq_metadata <- data.frame(
        age_group = c(rep("Adult", ncol(counts_adult)), rep("Aged", ncol(counts_aged)))
      )
      rownames(deseq_metadata) <- colnames(deseq_counts)
      dds <- DESeqDataSetFromMatrix(
        countData = round(deseq_counts),
        colData = deseq_metadata,
        design = ~ age_group
      )
      dds$age_group <- relevel(dds$age_group, ref = "Adult")
      dds <- DESeq(dds)
      res <- results(dds)
      res_df <- as.data.frame(res)
      subclase_clean <- gsub("/", "_", subclase)
      filename <- paste0(subclase_clean, "_Aged_vs_Adult.csv")
      write.csv(res_df, file.path(output_dir, filename))
    }, error = function(e) {
      cat(paste("    -> ERROR al procesar subclase:", subclase, "\n"))
      cat(paste("       Mensaje del error:", conditionMessage(e), "\n"))
    }) 
  }
}

base_output_path <- "data/Aging/"

run_aging_deseq_per_subclass(
  seurat_adult = seurat_adult,
  seurat_aged = seurat_aged,
  meta_paper_df = meta_paper,
  base_output_path = base_output_path,
  min_donors = 3
)