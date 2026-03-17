library(DESeq2)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)
library(stringr)
library(data.table)

seurat_adult <- readRDS("data/Normalized_Classified_Seurat_Adult")
seurat_aged  <- readRDS("data/Normalized_Classified_Seurat_Aged")
meta_paper   <- fread("data/metadata_aging.csv")

metadata_adult <- seurat_adult@meta.data
metadata_adult$library_prep <- sapply(rownames(metadata_adult), function(x) strsplit(x, "-")[[1]][2])
metadata_adult$barcode <- rownames(metadata_adult)
merged_adult <- merge(metadata_adult, meta_paper, by = "library_prep", all.x = TRUE)
rownames(merged_adult) <- merged_adult$barcode
excluir <- c("MOL NN_4", "CA2-FC-IG Glut_2", "OPC NN_1", "COP NN_1", "MFOL NN_3", "NFOL NN_2")
metadata_final_adult <- merged_adult[!merged_adult$supertype_label %in% excluir, ]
mtx_adult <- as.data.frame(seurat_adult@assays$RNA$counts)
mtx_final_adult <- mtx_adult[, rownames(metadata_final_adult)]
remove(seurat_adult); gc()

metadata_aged <- seurat_aged@meta.data
metadata_aged$library_prep <- sapply(rownames(metadata_aged), function(x) strsplit(x, "-")[[1]][2])
metadata_aged$barcode <- rownames(metadata_aged)
merged_aged <- merge(metadata_aged, meta_paper, by = "library_prep", all.x = TRUE)
rownames(merged_aged) <- merged_aged$barcode
metadata_final_aged <- merged_aged[!merged_aged$supertype_label %in% excluir, ]
mtx_aged <- as.data.frame(seurat_aged@assays$RNA$counts)
mtx_final_aged <- mtx_aged[, rownames(metadata_final_aged)]
remove(seurat_aged); gc()

genes_subclasses <- list(
  Uba52_L23  = c(gene = "Uba52", subclass = "L2/3 IT CTX"),
  Uba52_L6CT = c(gene = "Uba52", subclass = "L6 CT CTX"),
  Uba52_L6IT = c(gene = "Uba52", subclass = "L6 IT CTX"),
  Rpl36al    = c(gene = "Rpl36al", subclass = "L4 RSP-ACA")
)

get_deseq_padj_aging <- function(gene, subclass) {
  subclass_clean <- gsub("/", "_", subclass)
  fname <- paste0("data/Aging/DESeq_Adult_vs_Aged_per_Subclass/",
                  subclass_clean, "_Aged_vs_Adult.csv")
  df   <- fread(fname)
  padj <- df[df$V1 == gene, padj]
  return(padj)
}

run_deseq_get_normcounts_aging <- function(subclass, gene,
                                           mtx_final_adult, metadata_final_adult,
                                           mtx_final_aged,  metadata_final_aged) {
  
  meta_sub_adult <- metadata_final_adult[metadata_final_adult$predicted.id == subclass, ]
  meta_sub_aged  <- metadata_final_aged[metadata_final_aged$predicted.id == subclass, ]
  
  donors_adult <- as.data.frame(table(meta_sub_adult$external_donor_name)) %>% filter(Freq > 5)
  donors_aged  <- as.data.frame(table(meta_sub_aged$external_donor_name))  %>% filter(Freq > 5)
  
  n1 <- nrow(donors_adult)
  n2 <- nrow(donors_aged)
  
  counts1 <- list()
  counts2 <- list()
  
  for (p in seq_len(n1)) {
    counts1[[p]] <- mtx_final_adult[, rownames(meta_sub_adult[meta_sub_adult$external_donor_name == donors_adult$Var1[p], ]), drop = FALSE]
  }
  for (p in seq_len(n2)) {
    counts2[[p]] <- mtx_final_aged[, rownames(meta_sub_aged[meta_sub_aged$external_donor_name == donors_aged$Var1[p], ]), drop = FALSE]
  }
  
  DESEq_table <- matrix(nrow = nrow(mtx_final_adult), ncol = n1 + n2)
  colnames(DESEq_table) <- c(paste0(gsub(" ", "", subclass), "_", seq_len(n1), "_adult"),
                             paste0(gsub(" ", "", subclass), "_", seq_len(n2), "_aged"))
  rownames(DESEq_table) <- rownames(mtx_final_adult)
  
  for (k in seq_len(n1)) DESEq_table[, k]      <- rowSums(counts1[[k]])
  for (m in seq_len(n2)) DESEq_table[, n1 + m] <- rowSums(counts2[[m]])
  
  DESEq_table <- round(as.data.frame(DESEq_table))
  
  labels   <- c(rep(paste0(subclass, "_adult"), n1),
                rep(paste0(subclass, "_aged"),  n2))
  columnas <- data.frame(age_group = labels)
  rownames(columnas) <- colnames(DESEq_table)
  
  dds <- DESeqDataSetFromMatrix(countData = DESEq_table,
                                colData   = columnas,
                                design    = ~ age_group)
  colData(dds)$age_group <- factor(colData(dds)$age_group,
                                   levels = c(paste0(subclass, "_adult"),
                                              paste0(subclass, "_aged")))
  dds <- DESeq(dds)
  
  norm_counts <- counts(dds, normalized = TRUE)
  gene_norm   <- norm_counts[gene, ]
  
  result <- data.frame(
    replica   = colnames(DESEq_table),
    age_group = factor(c(rep("Adult", n1), rep("Aged", n2)),
                       levels = c("Adult", "Aged")),
    norm_expr = as.numeric(gene_norm)
  )
  
  return(result)
}

make_dotplot_panel_aging <- function(gene, subclass,
                                     mtx_final_adult, metadata_final_adult,
                                     mtx_final_aged,  metadata_final_aged) {
  
  pseudo <- run_deseq_get_normcounts_aging(subclass, gene,
                                           mtx_final_adult, metadata_final_adult,
                                           mtx_final_aged,  metadata_final_aged)
  padj   <- get_deseq_padj_aging(gene, subclass)
  
  stars <- ifelse(padj < 0.0001, "****",
                  ifelse(padj < 0.001,  "***",
                         ifelse(padj < 0.01,   "**",
                                ifelse(padj < 0.05,   "*", "ns"))))
  
  p_label <- if (padj < 0.001) {
    "p-adj < 0.001"
  } else {
    paste0("p-adj = ", round(padj, 3))
  }
  
  group_means <- pseudo %>%
    group_by(age_group) %>%
    summarise(norm_expr = mean(norm_expr), .groups = "drop")
  
  max_y     <- max(pseudo$norm_expr, na.rm = TRUE)
  pos_barra <- max_y * 1.07
  pos_stars <- max_y * 1.13
  pos_padj  <- max_y * 1.22
  
  p <- ggplot(pseudo, aes(x = age_group, y = norm_expr)) +
    geom_crossbar(data = group_means,
                  aes(y = norm_expr, ymin = norm_expr, ymax = norm_expr),
                  width = 0.32, linewidth = 0.45, fatten = 0, color = "black") +
    geom_jitter(width = 0.06, size = 1.2, alpha = 0.9, color = "black") +
    geom_segment(aes(x = 1, xend = 2, y = pos_barra, yend = pos_barra),
                 color = "black", linewidth = 0.25, inherit.aes = FALSE) +
    geom_segment(aes(x = 1, xend = 1,
                     y = pos_barra - max_y * 0.025, yend = pos_barra),
                 color = "black", linewidth = 0.25, inherit.aes = FALSE) +
    geom_segment(aes(x = 2, xend = 2,
                     y = pos_barra - max_y * 0.025, yend = pos_barra),
                 color = "black", linewidth = 0.25, inherit.aes = FALSE) +
    annotate("text", x = 1.5, y = pos_stars,
             label = stars, size = 2.5, hjust = 0.5) +
    annotate("text", x = 1.5, y = pos_padj,
             label = p_label, size = 1.8, hjust = 0.5, color = "grey30") +
    scale_x_discrete(labels = c("Adult" = "Adult", "Aged" = "Aged")) +
    coord_cartesian(ylim = c(NA, max_y * 1.35)) +
    theme_classic() +
    theme(
      legend.position  = "none",
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(size = 6),
      axis.text        = element_text(size = 6),
      axis.text.x      = element_text(angle = 25, hjust = 1),
      plot.title       = element_text(size = 7, face = "italic", hjust = 0.5),
      plot.subtitle    = element_text(size = 6, hjust = 0.5, color = "grey40"),
      axis.line        = element_line(linewidth = 0.25),
      axis.ticks       = element_line(linewidth = 0.25),
      plot.margin      = margin(4, 4, 2, 4)
    ) +
    labs(
      title    = gene,
      subtitle = subclass,
      y        = "DESeq2 norm. counts"
    )
  
  return(p)
}

panels <- lapply(genes_subclasses, function(x) {
  make_dotplot_panel_aging(
    gene              = x["gene"],
    subclass          = x["subclass"],
    mtx_final_adult   = mtx_final_adult,
    metadata_final_adult = metadata_final_adult,
    mtx_final_aged    = mtx_final_aged,
    metadata_final_aged  = metadata_final_aged
  )
})

fig_col <- wrap_plots(panels, ncol = 1)

print(fig_col)
ggsave(
  filename = "Figuras/Figure7D_New_RP_genes_aging_dotplot_column.svg",
  plot     = fig_col,
  width    = 1.6,
  height   = 5.6,
  units    = "in"
)

ggsave(
  filename = "Figuras/Figure7D_New_RP_genes_aging_dotplot_column.png",
  plot     = fig_col,
  width    = 1.6,
  height   = 5.6,
  units    = "in",
  dpi      = 300
)
