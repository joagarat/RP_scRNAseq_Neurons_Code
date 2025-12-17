library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(RColorBrewer)
library(svglite)
library(data.table)
library(sigclust2)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(tibble)
library(patchwork)
library(ggrepel)
library(cowplot)

load("data/Smartseq2dataset_seurat_filtered.RData")
seurat_split <- SplitObject(seurat_merged, split.by = "region_label")
seurat_split <- lapply(seurat_split, function(x) {
  DefaultAssay(x) <- "RNA"
  x
})
subclass_pass <- list(
  c("L5 PT CTX","L5 IT CTX","L4/5 IT CTX","L6 IT CTX","L6 CT CTX","L5 NP CTX","Pvalb","Vip","L2/3 IT CTX","Lamp5","Sst","Sst Chodl","Sncg","Car3","L6b CTX"),
  c("L4/5 IT CTX","L2/3 IT CTX","Lamp5","L5 IT CTX","L6 IT CTX","Vip","L6 CT CTX","L5 PT CTX","Car3","L5 NP CTX"),
  c("Car3","L2/3 IT CTX"),
  c("Vip","L5 NP CTX","L4/5 IT CTX","Sst","L6 IT CTX","L5 IT CTX","Lamp5","L2/3 IT CTX","L6 CT CTX","L6b CTX"),
  c("L5 PT CTX","L4 RSP-ACA"),
  c("L2/3 IT CTX","L6 IT CTX"),
  c("Sst","L5 NP CTX","L5 IT CTX","Vip","Lamp5","L4/5 IT CTX","Pvalb","L2/3 IT CTX","L6 IT CTX","L6 CT CTX","L6b CTX","Sncg"),
  c("DG","CA1-ProS","CA3","Lamp5","Vip","Sncg","Pvalb","Sst"),
  c("Sst","Sst Chodl","Sncg","Vip","L6 CT CTX","L4/5 IT CTX","L5 IT CTX","L2/3 IT CTX","Lamp5","L5 NP CTX","Pvalb","L6 IT CTX","L5 PT CTX"),
  c("L2/3 IT CTX","L4/5 IT CTX"),
  c("L4/5 IT CTX","L2/3 IT CTX","L5 IT CTX","L2/3 IT PPP","L6 CT CTX","L5 NP CTX","L6 IT CTX","L6b CTX","Sst","Lamp5","Vip","Pvalb")
)
RP = read.csv("data/RP.csv")
RP = as.character(RP$x)
RP =gsub("Rack1", "Gnb2l1", x = RP)
regions <- c("VISp","VIS","SSs","SSp","RSPv","RSP","MOp","HIP","ALM","AI","ACA")
lengths_of_vectors <- sapply(subclass_pass, length)
total_elements <- sum(lengths_of_vectors)
m <- matrix(nrow = 84,ncol = 88)
modified_vectors <- lapply(seq_along(subclass_pass), function(i) {
  paste(subclass_pass[[i]], regions[i], sep = "_")
})
final_vector <- unlist(modified_vectors)
colnames(m)<-final_vector
rownames(m)<-RP
z=0
for (i in seq(1,11)) {
  metadata_i = as.data.frame(seurat_split[[i]]@meta.data)
  data <- GetAssayData(object = seurat_split[[i]], layer = "data")
  data = as.data.frame(data)
  for (p in seq_len(length(subclass_pass[[i]]))) {
    z=z+1
    cells <- metadata_i[metadata_i$subclass_label%in%subclass_pass[[i]][p], ]
    datacellRP <- data[RP  ,colnames(data)%in%rownames(cells) ]
    m[ ,z] <- rowMeans(datacellRP)
  }
}
m_scaled = scale(m)
shc_result <- shc(as.matrix(t(m_scaled)), metric="euclidean", linkage="ward.D2")
hclust_tree <- shc_result$hc_dat 
k_clusters <- 18 
cluster_assignments <- cutree(hclust_tree, k = k_clusters)
inhibitory_subclasses <- c("Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip")
excitatory_subclasses <- c("L5 PT CTX", "L5 IT CTX", "L4/5 IT CTX", "L6 IT CTX", "L6 CT CTX", "L5 NP CTX", "L2/3 IT CTX", "Car3", "L6b CTX", "L4 RSP-ACA", "DG", "CA1-ProS", "CA3", "L2/3 IT PPP")
column_names <- colnames(m_scaled)
annotation_df <- data.frame(sample = column_names, row.names = column_names) %>%
  tidyr::separate(sample, into = c("subclass", "Region"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(
    Class = case_when(
      subclass %in% inhibitory_subclasses ~ "GABAergic",
      subclass %in% excitatory_subclasses ~ "Glutamatergic",
      TRUE ~ "Unknown"
    ),
    SigClust = as.factor(cluster_assignments)
  ) %>%
  dplyr::select(SigClust, subclass, Class, Region) %>%
  dplyr::rename(Subclass = subclass)
unique_clusters_for_base <- sort(unique(annotation_df$SigClust))
n_clusters_for_base <- length(unique_clusters_for_base)
sigclust_color_pool <- c(RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Accent"), RColorBrewer::brewer.pal(9, "Pastel1"), RColorBrewer::brewer.pal(9, "Set1"))
base_sigclust_colors <- unique(sigclust_color_pool)[1:n_clusters_for_base]
names(base_sigclust_colors) <- unique_clusters_for_base
subclass_to_cluster_map <- annotation_df %>%
  dplyr::count(Subclass, SigClust) %>%
  group_by(Subclass) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup()
final_subclass_colors <- list()
for (cluster_id in unique(subclass_to_cluster_map$SigClust)) {
  base_color <- base_sigclust_colors[as.character(cluster_id)]
  subclasses_in_cluster <- subclass_to_cluster_map %>%
    filter(SigClust == cluster_id) %>%
    pull(Subclass) 
  n_shades <- length(subclasses_in_cluster)
  color_generator <- colorRampPalette(c(adjustcolor(base_color, red.f = 1.2, green.f = 1.2, blue.f = 1.2), adjustcolor(base_color, red.f = 0.7, green.f = 0.7, blue.f = 0.7)))
  shades <- color_generator(n_shades)
  names(shades) <- subclasses_in_cluster
  final_subclass_colors[[length(final_subclass_colors) + 1]] <- shades
}
final_subclass_colors <- unlist(final_subclass_colors)
class_colors <- c("GABAergic" = "#882255", "Glutamatergic" = "#009988", "Unknown" = "grey")
unique_regions <- sort(unique(annotation_df$Region))
n_regions <- length(unique_regions)
region_colors <- RColorBrewer::brewer.pal(n_regions, "Set3")
names(region_colors) <- unique_regions
subclass_colors <- final_subclass_colors
annotation_colors <- list(
  Class = class_colors,
  Region = region_colors,
  Subclass = subclass_colors
)
top_annotation <- HeatmapAnnotation(
  df = annotation_df %>% select(Subclass, Class), 
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)
bottom_annotation <- HeatmapAnnotation(
  df = annotation_df %>% select(Region),
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
library(svglite)
dir.create("Figuras/Figura_3", recursive = TRUE, showWarnings = FALSE)
svglite("Figuras/Figura_3/Fig_3A.svg", 
        width = 6.5,  # Ancho final
        height = 12)  # Alto
ht <- Heatmap(m_scaled,
              name = "Z-score",
              col = col_fun,
              top_annotation = top_annotation,
              bottom_annotation = bottom_annotation,
              cluster_columns = hclust_tree,
              cluster_rows = TRUE,
              row_names_gp = gpar(fontsize = 8), 
              show_column_names = FALSE,
              show_row_dend = FALSE,
              show_column_dend = TRUE
)
draw(ht, 
     heatmap_legend_side = "left",    
     annotation_legend_side = "bottom" 
)

dev.off()
shc_result <- shc(as.matrix(t(m_scaled)), metric="euclidean", linkage="ward.D2")
plot(shc_result, hang=.1)
dir.create("Figuras/Figura_S3", recursive = TRUE, showWarnings = FALSE)
svglite("Figuras/Figura_S3/Fig_S3a.svg", 
        width = 15, 
        height = 7.5,
        pointsize = 8
)
plot(shc_result, hang=.1)
dev.off()
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
df = as.data.frame(m)
df$mean_all = rowMeans(df[, 1:ncol(m)], na.rm = TRUE)
column_names = colnames(df)
subtypes = unique(gsub("_.*", "", colnames(m)))
results_list <- list()
for (subtype in subtypes) {
  subtype_columns <- grep(paste0("^", subtype, "_"), column_names, value = TRUE)
  if(length(subtype_columns) > 0) {
    mean_subtype <- rowMeans(df[, subtype_columns, drop = FALSE], na.rm = TRUE)
    model_data <- data.frame(Average_All = df$mean_all, Mean_Subtype = mean_subtype)
    model_data_clean <- na.omit(model_data)
    if(nrow(model_data_clean) > 1) {
      model <- lm(Mean_Subtype ~ Average_All, data = model_data_clean)
      temp_residuals <- rep(NA, nrow(df))
      names(temp_residuals) <- rownames(df)
      temp_residuals[rownames(model_data_clean)] <- residuals(model)
      results_list[[subtype]] <- temp_residuals
    }
  }
}
residuals_df <- do.call(rbind, results_list)
colnames(residuals_df) <- rownames(df)
global_residual_mean <- mean(as.matrix(residuals_df), na.rm = TRUE)
global_residual_sd <- sd(as.matrix(residuals_df), na.rm = TRUE)
specificity_scores_df <- (residuals_df - global_residual_mean) / global_residual_sd
library(ComplexHeatmap)
library(svglite)
ht_drawn <- draw(Heatmap(as.matrix(m_scaled), cluster_rows = TRUE, show_heatmap_legend = FALSE))
row_order <- rev(row_order(ht_drawn))
gene_order_levels <- rownames(m_scaled)[row_order]
scores_long <- as.data.frame(specificity_scores_df) %>%
  rownames_to_column(var = "Subtype") %>%
  pivot_longer(
    cols = -Subtype,
    names_to = "Gene",
    values_to = "Specificity_Score"
  )
scores_long$Gene <- factor(scores_long$Gene, levels = gene_order_levels)
p_final_boxplot_final_style <- ggplot(scores_long, aes(x = Gene, y = Specificity_Score)) +
  geom_boxplot(coef = 0, fill = "#D6EAF8") + 
  geom_hline(yintercept = c(-2.5, 2.5), linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = NULL, 
    x = NULL,     
    y = "Specificity Score"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(
      angle = 90, 
      vjust = 0.5, 
      hjust = 1,
      size = 8 
    ) 
  )

output_dir <- "Figuras/Figura_3/"

svg(
  paste0(output_dir, "Fig_3B.svg"), 
  width = 2.5,  
  height = 11.2  
)
print(p_final_boxplot_final_style)
dev.off()

#S3B

lengths_of_vectors <- sapply(subclass_pass, length)
total_elements <- sum(lengths_of_vectors)
m <- matrix(nrow = 84,ncol = 88)
modified_vectors <- lapply(seq_along(subclass_pass), function(i) {
  paste(subclass_pass[[i]], regions[i], sep = "_")
})
final_vector <- unlist(modified_vectors)
colnames(m)<-final_vector
rownames(m)<-RP
z=0
for (i in seq(1,11)) {
  metadata_i = as.data.frame(seurat_split[[i]]@meta.data)
  data <- GetAssayData(object = seurat_split[[i]], layer = "data")
  data = as.data.frame(data)
  for (p in seq_len(length(subclass_pass[[i]]))) {
    z=z+1
    cells <- metadata_i[metadata_i$subclass_label%in%subclass_pass[[i]][p], ]
    datacellRP <- data[RP  ,colnames(data)%in%rownames(cells) ]
    m[ ,z] <- rowMeans(datacellRP)
  }
}

df = as.data.frame(m)
df$mean_all = rowMeans(df[, 1:ncol(m)], na.rm = TRUE)
column_names = colnames(df)
subtypes = unique(gsub("_.*", "", colnames(m)))
results_list <- list()
for (subtype in subtypes) {
  subtype_columns <- grep(paste0("^", subtype, "_"), column_names, value = TRUE)
  if(length(subtype_columns) > 0) {
    mean_subtype <- rowMeans(df[, subtype_columns, drop = FALSE], na.rm = TRUE)
    model_data <- data.frame(Average_All = df$mean_all, Mean_Subtype = mean_subtype)
    model_data_clean <- na.omit(model_data)
    if(nrow(model_data_clean) > 1) {
      model <- lm(Mean_Subtype ~ Average_All, data = model_data_clean)
      temp_residuals <- rep(NA, nrow(df))
      names(temp_residuals) <- rownames(df)
      temp_residuals[rownames(model_data_clean)] <- residuals(model)
      results_list[[subtype]] <- temp_residuals
    }
  }
}

residuals_df <- do.call(rbind, results_list)
colnames(residuals_df) <- rownames(df)

global_residual_mean <- mean(as.matrix(residuals_df), na.rm = TRUE)
global_residual_sd <- sd(as.matrix(residuals_df), na.rm = TRUE)

specificity_scores_df <- (residuals_df - global_residual_mean) / global_residual_sd

original_subtype_names <- rownames(residuals_df)
clean_subtype_names <- gsub("[ /]", "_", original_subtype_names)
name_map <- setNames(original_subtype_names, clean_subtype_names)
residuals_df_clean <- residuals_df
rownames(residuals_df_clean) <- clean_subtype_names
residuals_df_clean = as.data.frame(residuals_df_clean)
umbral_sd <- 2.5 * sd(as.matrix(residuals_df), na.rm = TRUE)
significant_points <- residuals_df_clean %>%
  rownames_to_column(var = "Subtype") %>% 
  pivot_longer(cols = -Subtype, names_to = "Gene", values_to = "Residual") %>%
  filter(abs(Residual) > umbral_sd)

significant_subtypes_clean <- unique(significant_points$Subtype)

if (length(significant_subtypes_clean) == 0) {
  cat(paste("No se encontraron genes con residuales >", round(umbral_sd, 2), "unidades.\n"))
} else {
  cat(paste("Generando scatterplots para", length(significant_subtypes_clean), "subtipos con genes específicos.\n"))
  plots_list <- list()
  plotting_df_base <- as.data.frame(m)
  plotting_df_base$mean_all <- rowMeans(plotting_df_base, na.rm = TRUE)
  for (i in seq_along(original_subtype_names)) {
    original_name <- original_subtype_names[i]
    clean_name <- clean_subtype_names[i]
    cols_to_avg <- grep(paste0("^", original_name, "_"), colnames(m), value = TRUE)
    if (length(cols_to_avg) > 0) {
      plotting_df_base[[paste0("mean_", clean_name)]] <- rowMeans(m[, cols_to_avg, drop=FALSE], na.rm = TRUE)
    }
  }
  plotting_df_base <- plotting_df_base %>% rownames_to_column(var = "Gene")
  for (current_clean_subtype in significant_subtypes_clean) {
    significant_genes_for_this_subtype <- significant_points %>%
      filter(Subtype == current_clean_subtype) %>%
      pull(Gene)
    data_for_plot <- plotting_df_base %>%
      select(Gene, mean_all, mean_specific = paste0("mean_", current_clean_subtype)) %>%
      mutate(is_significant = Gene %in% significant_genes_for_this_subtype)
    original_title_name <- name_map[current_clean_subtype]
    p <- ggplot(data_for_plot, aes(x = mean_all, y = mean_specific)) +
      geom_point(data = ~filter(., !is_significant), color = "grey70", alpha = 0.7) +
      geom_point(data = ~filter(., is_significant), color = "red", size = 2.5) +
      geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
      geom_text_repel(
        data = ~filter(., is_significant),
        aes(label = Gene),
        size = 2.8,                # ≈ 8 pt mínimo
        max.overlaps = Inf,
        box.padding = 0.5,
        segment.color = "grey50"
      ) +
      labs(
        x = "Average RP Expression",
        y = paste(original_title_name, "RP Expression")
      ) +
      theme_minimal_grid(font_size = 8)   # nunca menor a 8
    
    plots_list[[current_clean_subtype]] <- p
  }
  output_dir <- "Figuras/Figura_S3/"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  n_plots <- length(plots_list)
  cols <- min(5, n_plots)

  max_width <- 7.5
  fig_width <- min(max_width, cols * 3)  
  fig_height <- fig_width
  
  combined_plots <- wrap_plots(plots_list, ncol = cols)
  
  svglite(paste0(output_dir, "FigS3B.svg"),
          width = fig_width, height = fig_height)
  print(combined_plots)
  dev.off()
}


