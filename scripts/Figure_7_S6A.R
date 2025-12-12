library(ggplot2)
library(patchwork)
library(purrr)
library(data.table)
library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dendextend) 
library(sigclust2)
library(tidyr)
library(RColorBrewer)
library(svglite)
library(DESeq2)
library(stringr)
library(ggrepel)

DEG_dir_adult = "data/Aging/DESeq_Adult_Between_Subclasses/"
filelist = list.files(path = DEG_dir_adult, pattern="*.csv$")
df_input_list_Allen <- lapply(filelist, fread)
names(df_input_list_Allen) <- gsub(filelist, pattern="\\..*", replacement="")
df_Allen <- bind_rows(df_input_list_Allen, .id = "id")
df_sin_na_Allen <- na.omit(df_Allen)
DEG_dir_aged = "data/Aging/DESeq_Aged_Between_Subclasses/"
filelist = list.files(path = DEG_dir_aged, pattern="*.csv$")
df_input_list_Aged <- lapply(filelist, fread)
names(df_input_list_Aged) <- gsub(filelist, pattern="\\..*", replacement="")
df_Aged <- bind_rows(df_input_list_Aged, .id = "id")
df_sin_na_Aged <- na.omit(df_Aged)
ordenar_id <- function(id) {
  subtipos <- unlist(strsplit(id, "vs")) 
  id_ordenado <- paste(sort(subtipos), collapse = "vs")  
  invertido <- ifelse(id == id_ordenado, FALSE, TRUE) 
  return(c(id_ordenado, invertido))  
}

df_sin_na_Allen <- df_sin_na_Allen %>%
  mutate(id_info = map(id, ordenar_id)) %>%
  mutate(id_estandarizado = sapply(id_info, `[`, 1),
         invertido = as.logical(sapply(id_info, `[`, 2))) %>%
  select(-id_info)

df_sin_na_Aged <- df_sin_na_Aged %>%
  mutate(id_info = map(id, ordenar_id)) %>%
  mutate(id_estandarizado = sapply(id_info, `[`, 1),
         invertido = as.logical(sapply(id_info, `[`, 2))) %>%
  select(-id_info)   

tabla_comun_sin_na <- inner_join(df_sin_na_Allen, df_sin_na_Aged, 
                                 by = c("id_estandarizado", "V1"), 
                                 suffix = c("_Adult", "_Aged"))


tabla_comun_sin_na <- tabla_comun_sin_na %>%
  mutate(log2FoldChange_Adult = ifelse(invertido_Adult, -log2FoldChange_Adult, log2FoldChange_Adult),
         log2FoldChange_Aged = ifelse(invertido_Aged, -log2FoldChange_Aged, log2FoldChange_Aged))
tabla_comun_sin_na <- tabla_comun_sin_na %>%
  mutate(FoldChange_Adult = 2^log2FoldChange_Adult,
         FoldChange_Aged = 2^log2FoldChange_Aged)
RP <- read.csv("data/RP.csv")
RP <- as.character(RP$x)
RP <- gsub("Rack1", "Gnb2l1", x = RP)

correlacion_global <- cor(tabla_comun_sin_na$log2FoldChange_Adult, tabla_comun_sin_na$log2FoldChange_Aged, method = "spearman")

tabla_comun_sin_na <- tabla_comun_sin_na %>%
  mutate(gene_type = ifelse(V1 %in% RP, "RP DEGs", "Other DEGs")) %>%
  arrange(desc(gene_type))
max_abs_logfc <- max(abs(c(tabla_comun_sin_na$log2FoldChange_Adult, tabla_comun_sin_na$log2FoldChange_Aged)), na.rm = TRUE) * 1.05
print(table(tabla_comun_sin_na$gene_type))
tabla_comun_sin_na <- tabla_comun_sin_na %>%
  mutate(gene_type = ifelse(toupper(V1) %in% toupper(RP), "RP DEGs", "Other DEGs")) %>%
  arrange(gene_type)
max_abs_logfc <- max(abs(c(tabla_comun_sin_na$log2FoldChange_Adult, tabla_comun_sin_na$log2FoldChange_Aged)), na.rm = TRUE) * 1.05

panel_a_global_corr <- ggplot(tabla_comun_sin_na, aes(x = log2FoldChange_Adult, y = log2FoldChange_Aged)) +
  geom_point(aes(color = gene_type), alpha = 0.7, size = 2) +
  scale_color_manual(name = "Gene Type", values = c("RP DEGs" = "#e41a1c", "Other DEGs" = "grey80")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-max_abs_logfc, max_abs_logfc), ylim = c(-max_abs_logfc, max_abs_logfc)) +
  annotate("text", x = -max_abs_logfc * 0.95, y = max_abs_logfc * 0.95, 
           label = paste("Spearman's ρ =", round(correlacion_global, 2)), hjust = 0, vjust = 1, size = 5) +
  theme_minimal(base_size = 12) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  theme(legend.position = c(0.85, 0.2), plot.title = element_text(face="bold")) +
  labs(title = "A) Global DGE Correlation", x = "Subclass-Level log2FC (Adult)", y = "Subclass-Level log2FC (Aged)")

print(panel_a_global_corr)

tabla_comun_sin_na_final <- tabla_comun_sin_na %>%
  mutate(
    gene_type = ifelse(toupper(V1) %in% toupper(RP), "RP", "Non-RP"),
    gene_type = factor(gene_type, levels = c("Non-RP", "RP"))
  ) %>%
  arrange(gene_type)

max_abs_logfc <- max(abs(c(tabla_comun_sin_na_final$log2FoldChange_Adult, tabla_comun_sin_na_final$log2FoldChange_Aged)), na.rm = TRUE) * 1.05

correlacion_global <- cor(
  tabla_comun_sin_na_final$log2FoldChange_Adult,
  tabla_comun_sin_na_final$log2FoldChange_Aged,
  method = "spearman", 
  use = "pairwise.complete.obs"
)

max_abs_logfc_global <- max(abs(c(tabla_comun_sin_na_final$log2FoldChange_Adult, tabla_comun_sin_na_final$log2FoldChange_Aged)), na.rm = TRUE) * 1.05

plot_a_global <- ggplot(tabla_comun_sin_na_final, aes(x = log2FoldChange_Adult, y = log2FoldChange_Aged)) +
  geom_point(aes(color = gene_type), alpha = 0.6, size = 0.8) + # <-- Ajustado a size = 1.2
  scale_color_manual(
    name = "Gene Set", 
    values = c("RP" = "#e41a1c", "Non-RP" = "grey80")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-max_abs_logfc_global, max_abs_logfc_global), ylim = c(-max_abs_logfc_global, max_abs_logfc_global)) +
  annotate(
    "text", 
    x = -max_abs_logfc_global * 0.95, 
    y = max_abs_logfc_global * 0.95, 
    label = paste("Spearman's ρ =", round(correlacion_global, 2)), 
    hjust = 0, vjust = 1, size = 2.5 
  ) +
  theme_classic() + 
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,-5,-5,-5),
    aspect.ratio = 1,
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7)
  ) +

  guides(color = guide_legend(override.aes = list(size = 3))) + 
  labs(title = NULL, x = "Subclass-Level Log2FC (Adult)", y = "Subclass-Level Log2FC (Aged)") 

logfc_threshold <- log2(1.5)
padj_threshold <- 0.05

scatter_data_rps <- tabla_comun_sin_na %>%
  filter(toupper(V1) %in% toupper(RP)) %>%
  mutate(
    color_category = case_when(
      padj_Adult < padj_threshold & abs(log2FoldChange_Adult) > logfc_threshold &
        padj_Aged < padj_threshold & abs(log2FoldChange_Aged) > logfc_threshold &
        sign(log2FoldChange_Adult) == sign(log2FoldChange_Aged) ~ "Concordant RP",
      padj_Adult < padj_threshold & abs(log2FoldChange_Adult) > logfc_threshold &
        padj_Aged < padj_threshold & abs(log2FoldChange_Aged) > logfc_threshold &
        sign(log2FoldChange_Adult) != sign(log2FoldChange_Aged) ~ "Discordant RP", 
    ),
    color_category = factor(color_category, levels = c("Non-Replicated RP", "Discordant RP", "Concordant RP"))
  ) %>%
  arrange(color_category)

max_abs_logfc_rp <- max(abs(c(scatter_data_rps$log2FoldChange_Adult, scatter_data_rps$log2FoldChange_Aged)), na.rm = TRUE) * 1.05

correlacion_rp <- cor(scatter_data_rps$log2FoldChange_Adult, scatter_data_rps$log2FoldChange_Aged, method = "spearman")

plot_b_rp <- ggplot(scatter_data_rps, aes(x = log2FoldChange_Adult, y = log2FoldChange_Aged)) +
  geom_point(aes(color = color_category), size = 0.8, alpha = 0.8) + # <-- Ajustado a size = 1.2
  
  scale_color_manual(
    name = "RP Concordance", 
    values = c("Concordant RP" = "#e41a1c", "Discordant RP" = "#377eb8", "Non-Replicated RP" = "grey80")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-max_abs_logfc_rp, max_abs_logfc_rp), ylim = c(-max_abs_logfc_rp, max_abs_logfc_rp)) +
  annotate(
    "text", 
    x = -max_abs_logfc_rp * 0.95, 
    y = max_abs_logfc_rp * 0.95, 
    label = paste("Spearman's ρ =", round(correlacion_rp, 2)), 
    hjust = 0, vjust = 1, size = 2.5
  ) +
  theme_classic() + 
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,-5,-5,-5),
    aspect.ratio = 1,
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7)
  ) +
  guides(color = guide_legend(
    nrow = 2, 
    override.aes = list(size = 3) 
  )) +
  labs(title = NULL, x = "Subclass-Level Log2FC (Adult)", y = NULL) 

scatters_side_by_side <- plot_a_global + plot_b_rp 

scatters_side_by_side

ggsave(
  "Figuras/Figura_7/Figura_7AB-2.png", 
  plot = scatters_side_by_side, 
  width = 4.0,  
  height = 2.5, 
  units = "in",
  dpi = 300
)

seurat_adult <- readRDS("data/Normalized_Classified_Seurat_Adult") 
seurat_aged <- readRDS("data/Normalized_Classified_Seurat_Aged")
RP <- read.csv("data/RP.csv")
RP <- as.character(RP$x)

calculate_avg_expression_matrix <- function(seurat_object, rp_genes) {
  avg_expr <- AverageExpression(
    object = seurat_object,
    features = rp_genes,
    group.by = "predicted.id",
    assay = "RNA"
  )
  return(avg_expr$RNA)
}

m_adult <- calculate_avg_expression_matrix(seurat_adult, RP)
m_aged <- calculate_avg_expression_matrix(seurat_aged, RP)
common_genes <- intersect(rownames(m_adult), rownames(m_aged))
common_subtypes <- intersect(colnames(m_adult), colnames(m_aged))
m_adult <- m_adult[common_genes, common_subtypes]
m_aged <- m_aged[common_genes, common_subtypes]
colnames(m_adult) <- paste0(colnames(m_adult), "_adult")
colnames(m_aged) <- paste0(colnames(m_aged), "_aged")
m_combined <- cbind(m_adult, m_aged)

m_combined_zscore_by_col <- scale(m_combined)
m_combined_zscore_by_col[is.na(m_combined_zscore_by_col)] <- 0

hc_cols <- hclust(dist(t(m_combined_zscore_by_col)), method = "ward.D2")

dend <- as.dendrogram(hc_cols)

k <- 9 
clusters <- cutree(hc_cols, k = k)
print(table(clusters))
cluster_composition <- data.frame(
  sample = names(clusters),
  cluster_id = clusters
) %>%
  mutate(
    subclass = gsub("_adult|_aged", "", sample),
    age = ifelse(grepl("_adult", sample), "Adult", "Aged")
  )

print(as.data.frame(
  cluster_composition %>%
    group_by(cluster_id) %>%
    summarise(
      n_samples = n(),
      subclases_presentes = paste(sort(unique(subclass)), collapse = ", ")
    )
))

desired_cluster_order <- c(4,1,6,3,5,7,8,9,2) 

ordered_labels <- labels(dend)
clusters_in_dend_order <- clusters[ordered_labels] 
weights_map <- setNames(1:length(desired_cluster_order), desired_cluster_order)
reordering_weights <- weights_map[as.character(clusters_in_dend_order)]

dend_reordered <- reorder(dend, wts = reordering_weights, agglo.fun = "mean")

final_column_order <- labels(dend_reordered)
m_reordered_for_shc <- m_combined_zscore_by_col[, final_column_order]
shc_result_reordered <- shc(as.matrix(t(m_reordered_for_shc)), metric = "euclidean", linkage = "ward.D2")
column_names <- colnames(m_combined_zscore_by_col)
column_names
inhibitory_subclasses <- c("Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip")
excitatory_subclasses <- c("Car3", "CR", "L2  IT ENTl", "L2/3 IT CTX", "L2/3 IT ENTl", "L2/3 IT PPP", "L4 RSP-ACA", "L4/5 IT CTX", "L5 IT CTX", "L5 IT TPE-ENT", "L5 NP CTX", "L5 PT CTX", "L6 CT CTX", "L6 IT CTX", "L6b CTX")


column_names <- colnames(m_combined_zscore_by_col)

annotation_df <- data.frame(sample = column_names, row.names = column_names) %>%
  separate(sample, into = c("Subclass", "Age"), sep = "_", remove = FALSE) %>% 
  mutate(
    Class = case_when(
      Subclass %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip") ~ "GABAergic",
      TRUE ~ "Glutamatergic"
    ),
    SigClust = as.numeric(factor(Subclass)) %% 8 + 1 
  )

if (!exists("base_sigclust_colors")) {
  base_sigclust_colors <- RColorBrewer::brewer.pal(8, "Set2")
  names(base_sigclust_colors) <- 1:8
}

subclass_to_cluster_map <- annotation_df %>%
  count(Subclass, SigClust) %>%
  group_by(Subclass) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup()

final_subclass_colors <- list()
for (cluster_id in unique(subclass_to_cluster_map$SigClust)) {
  base_color <- tryCatch(base_sigclust_colors[as.character(cluster_id)], error = function(e) "#CCCCCC")
  if(is.na(base_color)) base_color <- "#CCCCCC"
  
  subclasses_in_cluster <- subclass_to_cluster_map %>%
    filter(SigClust == cluster_id) %>%
    pull(Subclass)
  n_shades <- length(subclasses_in_cluster)
  color_generator <- colorRampPalette(c(adjustcolor(base_color, red.f = 1.2, green.f = 1.2, blue.f = 1.2), adjustcolor(base_color, red.f = 0.7, green.f = 0.7, blue.f = 0.7)))
  shades <- color_generator(n_shades)
  names(shades) <- subclasses_in_cluster
  final_subclass_colors[[length(final_subclass_colors) + 1]] <- shades
}
subclass_colors <- unlist(final_subclass_colors)

class_colors <- c("GABAergic" = "#882255", "Glutamatergic" = "#009988", "Unknown" = "grey")

age_colors <- c("adult" = "#E69F00", "aged" = "#56B4E9") 

annotation_colors <- list(
  Class = class_colors,
  Subclass = subclass_colors,
  Age = age_colors
)

top_annotation <- HeatmapAnnotation(
  df = annotation_df %>% select(Class, Subclass), 
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  show_legend = TRUE,
  annotation_legend_param = list(
    Subclass = list(ncol = 3)
  )
)

bottom_annotation <- HeatmapAnnotation(
  df = annotation_df %>% select(Age), 
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  show_legend = TRUE
)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

svglite("Figuras/Figura_7/Figura_7C_Heatmap.svg", 
        width = 10.5,  # Ancho completo
        height = 9.5)   # Alto

ht <- Heatmap(m_combined_zscore_by_col,
              name = "Z-score",
              col = col_fun,
              top_annotation = top_annotation,
              bottom_annotation = bottom_annotation,
              cluster_columns = rev(dend_reordered), 
              cluster_rows = TRUE,
              show_column_names = FALSE, 
              row_names_gp = gpar(fontsize = 8),
              show_row_dend = TRUE,
              show_column_dend = TRUE
)

draw(ht, 
     heatmap_legend_side = "left",     
     annotation_legend_side = "bottom",
     merge_legend = TRUE 
)

dev.off()

seurat_adult <- readRDS("Normalized_Classified_Seurat_Adult")
seurat_aged <- readRDS("Normalized_Classified_Seurat_Aged")
meta_paper <- fread("data/metadata_aging.csv")

RP <- read.csv("data/RP.csv")
RP <- as.character(RP$x)
RP <- gsub("Rack1", "Gnb2l1", x = RP)
path_aging_effect <- "data/Aging/DESeq_Adult_vs_Aged_per_Subclass/"
filelist_aging_effect <- list.files(path = path_aging_effect, pattern="*.csv$")

if (length(filelist_aging_effect) == 0) {
  stop("No se encontraron archivos .csv en la ruta especificada: ", path_aging_effect)
}

df_list_aging_effect <- lapply(filelist_aging_effect, function(f) fread(file.path(path_aging_effect, f)))
names(df_list_aging_effect) <- gsub(filelist_aging_effect, pattern="\\..*", replacement="")
df_aging_effect <- bind_rows(df_list_aging_effect, .id = "subclass_comparison")
df_aging_effect_clean <- na.omit(df_aging_effect)

volcano_data_aging <- df_aging_effect_clean %>%
  mutate(
    is_RP = toupper(V1) %in% toupper(RP),
    color_category = case_when(
      padj < 0.05 & abs(log2FoldChange) > 0.58 & is_RP ~ "RP DEGs",
      padj < 0.05 & abs(log2FoldChange) > 0.58 & !is_RP ~ "Non-RP DEGs",
      TRUE ~ "Not Significant"
    ),
    color_category = factor(color_category, 
                            levels = c("Not Significant", "Non-RP DEGs", "RP DEGs"))
  ) %>%
  arrange(color_category) 
genes_a_etiquetar <- volcano_data_aging %>% 
  filter(color_category == "RP DEGs") %>%
  arrange(padj) %>%
  slice_head(n = 10)

volcano_plot_aging <- ggplot(volcano_data_aging, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color_category), size = 0.8, alpha = 0.7) +
  scale_color_manual(
    name = NULL, 
    values = c(
      "RP DEGs" = "#e41a1c", 
      "Non-RP DEGs" = "#FCCDE5", 
      "Not Significant" = "grey80"
    )
  ) +
  geom_text_repel(
    data = genes_a_etiquetar, 
    aes(label = V1), 
    size = 2,               
    box.padding = 0.3,      
    point.padding = 0.1,    
    min.segment.length = 0, 
    max.overlaps = Inf,
    segment.size = 0.2
  ) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.3) +
  theme_classic() + 
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,-5,-5,-5), 
    axis.title = element_text(size = 9), 
    axis.text = element_text(size = 7),
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = NULL,
    x = "Log2FC (Aged/Adult)", 
    y = "-log10(p-adj)"
  )

print(volcano_plot_aging)

ggsave(
  "Figuras/Figura_7/Figura_7_Volcano.png", 
  plot = volcano_plot_aging, 
  width = 3, 
  height = 2.5, 
  units = "in",
  dpi = 300 
)

svglite(
  "Figuras/Figura_S6/Figura_S6A.svg", 
  width = 15, 
  height = 10
)

plot(shc_result_reordered, hang = -1, main = "Hierarchical Clustering with Significance (Reordered)")
dev.off()
