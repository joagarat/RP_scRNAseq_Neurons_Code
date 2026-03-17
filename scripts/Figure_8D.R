library(DESeq2)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)
library(stringr)
library(data.table)

RP <- read.csv("data/RP.csv")
RP <- as.character(RP$x)
RP <- gsub("Rack1", "Gnb2l1", x = RP)

combined_subclass_annotated_0_7_confidence_filtered <- LoadSeuratRds("data/Hing_2024_PFC_stress/SeuratProject.Rds")
mtx      <- as.data.frame(combined_subclass_annotated_0_7_confidence_filtered@assays$RNA$counts)
metadata <- combined_subclass_annotated_0_7_confidence_filtered@meta.data
metadata$replica_group <- str_extract(metadata$replica, "^[^_]+")
metadata <- as.data.frame(metadata)

metadata_high <- metadata[metadata$replica_group == "high", ]
metadata_low  <- metadata[metadata$replica_group == "low",  ]

genes_subclasses <- list(
  Rpl15 = "Pvalb",
  Uba52 = "L5 NP CTX",
  Rps2  = "Vip",
  Rps27 = "Vip"
)

stress_colors <- c("low" = "black", "high" = "black")

get_deseq_padj <- function(gene, subclass) {
  fname <- paste0("data/Hing_2024_PFC_stress_DEGs_High_vs_Low/",
                  subclass, "high_vs",
                  subclass, "low.csv")
  df   <- fread(fname)
  padj <- df[df$V1 == gene, padj]
  return(padj)
}

run_deseq_get_normcounts <- function(subclass, gene, mtx, metadata_high, metadata_low, metadata) {
  
  counts1 <- list()
  counts2 <- list()
  
  rats1       <- metadata_high[metadata_high$predicted.id == subclass, "replica"]
  rats1_stats <- as.data.frame(table(as.data.frame(rats1)))
  rats1_stats_mayor_a10 <- rats1_stats[rats1_stats$Freq > 5, ]
  
  rats2       <- metadata_low[metadata_low$predicted.id == subclass, "replica"]
  rats2_stats <- as.data.frame(table(as.data.frame(rats2)))
  rats2_stats_mayor_a10 <- rats2_stats[rats2_stats$Freq > 5, ]
  
  stopifnot(length(rats1_stats_mayor_a10$rats1) > 2,
            length(rats2_stats_mayor_a10$rats2) > 2)
  
  for (p in seq_along(rats1_stats_mayor_a10$rats1)) {
    counts1[[p]] <- mtx[, rownames(metadata[metadata_high$replica == rats1_stats_mayor_a10$rats1[p] &
                                              metadata_high$predicted.id == subclass, ])]
  }
  for (p in seq_along(rats2_stats_mayor_a10$rats2)) {
    counts2[[p]] <- mtx[, rownames(metadata_low[metadata_low$replica == rats2_stats_mayor_a10$rats2[p] &
                                                  metadata_low$predicted.id == subclass, ])]
  }
  
  n1 <- length(rats1_stats_mayor_a10$rats1)
  n2 <- length(rats2_stats_mayor_a10$rats2)
  
  DESEq_table <- matrix(nrow = nrow(mtx), ncol = n1 + n2)
  colnames(DESEq_table) <- append(
    gsub(" ", "", paste0(subclass, "_", 1:n1, "_high")),
    gsub(" ", "", paste0(subclass, "_", 1:n2, "_low"))
  )
  rownames(DESEq_table) <- rownames(mtx)
  
  for (k in seq_len(n1)) DESEq_table[, k]      <- rowSums(counts1[[k]])
  for (m in seq_len(n2)) DESEq_table[, n1 + m] <- rowSums(counts2[[m]])
  
  DESEq_table <- round(as.data.frame(DESEq_table))
  
  labels   <- append(rep(paste0(subclass, "_high"), n1),
                     rep(paste0(subclass, "_low"),  n2))
  columnas <- data.frame(cell = labels)
  rownames(columnas) <- colnames(DESEq_table)
  
  Deseq_object <- DESeqDataSetFromMatrix(countData = DESEq_table,
                                         colData   = columnas,
                                         design    = ~cell)
  colData(Deseq_object)$cell <- factor(colData(Deseq_object)$cell,
                                       levels = c(paste0(subclass, "_high"),
                                                  paste0(subclass, "_low")))
  Deseq_object <- DESeq(Deseq_object)
  
  norm_counts <- counts(Deseq_object, normalized = TRUE)
  gene_norm   <- norm_counts[gene, ]
  
  result <- data.frame(
    replica      = colnames(DESEq_table),
    stress_group = factor(ifelse(grepl("_high", colnames(DESEq_table)), "high", "low"),
                          levels = c("low", "high")),
    norm_expr    = as.numeric(gene_norm)
  )
  
  return(result)
}

make_dotplot_panel <- function(gene, subclass, mtx, metadata_high, metadata_low, metadata, colors) {
  
  pseudo <- run_deseq_get_normcounts(subclass, gene, mtx, metadata_high, metadata_low, metadata)
  padj   <- get_deseq_padj(gene, subclass)
  
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
    group_by(stress_group) %>%
    summarise(norm_expr = mean(norm_expr), .groups = "drop")
  
  max_y     <- max(pseudo$norm_expr, na.rm = TRUE)
  pos_barra <- max_y * 1.07
  pos_stars <- max_y * 1.13
  pos_padj  <- max_y * 1.22
  
  p <- ggplot(pseudo, aes(x = stress_group, y = norm_expr, color = stress_group)) +
    geom_crossbar(data = group_means,
                  aes(y = norm_expr, ymin = norm_expr, ymax = norm_expr),
                  width = 0.32, linewidth = 0.45, fatten = 0) +
    geom_jitter(width = 0.06, size = 1.2, alpha = 0.9) +

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
    scale_color_manual(values = colors) +
    scale_x_discrete(labels = c("low" = "Low", "high" = "High")) +
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

panels <- mapply(
  FUN      = make_dotplot_panel,
  gene     = names(genes_subclasses),
  subclass = unlist(genes_subclasses),
  MoreArgs = list(mtx           = mtx,
                  metadata_high = metadata_high,
                  metadata_low  = metadata_low,
                  metadata      = metadata,
                  colors        = stress_colors),
  SIMPLIFY = FALSE
)

fig_col <- wrap_plots(panels, ncol = 1)

ggsave(
  filename = "Figuras/Figura8D_New.svg",
  plot     = fig_col,
  width    = 1.6,
  height   = 5.6,
  units    = "in"
)

ggsave(
  filename = "Figuras/Figura8D_New.png",
  plot     = fig_col,
  width    = 1.6,
  height   = 5.6,
  units    = "in",
  dpi      = 300
)
