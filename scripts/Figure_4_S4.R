library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(RColorBrewer)
library(svglite)
library(tidyr)
library(ggpubr)   
library(patchwork)
library(scales) 
library(cowplot)
load("data/Smartseq2dataset_seurat_filtered.RData")
seurat_split <- SplitObject(seurat_merged, split.by = "region_label")
RP = read.csv("data/RP.csv")
RP = as.character(RP$x)
RP =gsub("Rack1", "Gnb2l1", x = RP)
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
regions_pass_DESeq2 <- c("VISp","VIS","SSs","SSp","RSPv","RSP","MOp","HIP","ALM","AI","ACA")
geneswithparalogs <- c("Rpl7", "Rpl7l1", "Rpl22", "Rpl22l1", "Rps27", "Rps27l")
procesar_seurat <- function(seurat_obj, genes, subclasses, region) {
  datos <- FetchData(seurat_obj, c(genes, "subclass_label", "donor_label"))
  datos_long <- datos %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
    mutate(presence = expression > 0)
  datos_filtered <- datos_long %>%
    filter(subclass_label %in% subclasses)
  resultado <- datos_filtered %>%
    group_by(subclass_label, donor_label, gene) %>%
    summarise(presencia_pct = mean(presence) * 100, .groups = 'drop') %>%
    mutate(region = region)
  return(resultado)
}

resultados <- lapply(1:length(seurat_split), function(i) {
  procesar_seurat(seurat_split[[i]], geneswithparalogs, subclass_pass[[i]], names(seurat_split)[i])
})

resultados_finales <- bind_rows(resultados)

results_Rpl7_Rpl7l1 <- resultados_finales %>%
  filter(gene %in% c("Rpl7", "Rpl7l1")) %>%
  compare_means(presencia_pct ~ gene, data = ., method = "wilcox.test", group.by = "subclass_label")

results_Rpl22_Rpl22l1 <- resultados_finales %>%
  filter(gene %in% c("Rpl22", "Rpl22l1")) %>%
  compare_means(presencia_pct ~ gene, data = ., method = "wilcox.test", group.by = "subclass_label")

results_Rps27_Rps27l <- resultados_finales %>%
  filter(gene %in% c("Rps27", "Rps27l")) %>%
  compare_means(presencia_pct ~ gene, data = ., method = "wilcox.test", group.by = "subclass_label")

stat_results_combined <- bind_rows(
  results_Rpl7_Rpl7l1 %>% mutate(paralog_pair = "Rpl7 / Rpl7l1"),
  results_Rpl22_Rpl22l1 %>% mutate(paralog_pair = "Rpl22 / Rpl22l1"),
  results_Rps27_Rps27l %>% mutate(paralog_pair = "Rps27 / Rps27l")
)

data_with_pairs <- resultados_finales %>%
  mutate(paralog_pair = case_when(
    gene %in% c("Rpl7", "Rpl7l1")   ~ "Rpl7 / Rpl7l1",
    gene %in% c("Rpl22", "Rpl22l1") ~ "Rpl22 / Rpl22l1",
    gene %in% c("Rps27", "Rps27l")  ~ "Rps27 / Rps27l"
  ))

y_positions <- data_with_pairs %>%
  group_by(paralog_pair, subclass_label) %>%
  summarise(max_y = max(presencia_pct, na.rm = TRUE), .groups = 'drop')

stat_results_formatted_manual <- stat_results_combined %>%
  left_join(y_positions, by = c("paralog_pair", "subclass_label")) %>%
  mutate(y.position = max_y + 5) # Ajusta el offset (5) si es necesario

paralog_colors <- c("Rpl7"="#1f78b4", "Rpl7l1"="#a6cee3", "Rpl22"="#33a02c", "Rpl22l1"="#b2df8a", "Rps27"="#e31a1c", "Rps27l"="#fb9a99")

mean_values <- data_with_pairs %>%
  group_by(paralog_pair, subclass_label, gene) %>%
  summarise(mean_pct = mean(presencia_pct, na.rm = TRUE), .groups = 'drop')

mean_diffs <- mean_values %>%
  group_by(paralog_pair, subclass_label) %>%
  summarise(
    mean_difference = abs(mean_pct[1] - mean_pct[2]),
    .groups = 'drop'
  )

stat_results_with_diff <- stat_results_formatted_manual %>%
  left_join(mean_diffs, by = c("paralog_pair", "subclass_label"))

stat_results_simple_asterisk <- stat_results_with_diff %>%
  filter(
    p.adj <= 0.05 &    
      mean_difference >= 5  
  ) %>%
  mutate(label = "*")

p_boxplot <- ggplot(
  data_with_pairs,
  aes(x = subclass_label, y = presencia_pct, fill = gene, color = gene)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.8, linewidth = 0.7) +
  
  stat_pvalue_manual(
    stat_results_simple_asterisk,
    x = "subclass_label",
    label = "label",
    tip.length = 0.01
  ) +
  facet_wrap(~ paralog_pair, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = paralog_colors) +
  scale_color_manual(values = paralog_colors) +
  scale_y_continuous(limits = c(NA, 105), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "% of Expressing Cells",
    fill = "Gene",
    color = "Gene" 
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic() + 
  theme(
    legend.position = "top",
    panel.spacing.x = unit(0.5, "lines"), 
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7), 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 7),
    strip.text = element_text(face = "bold.italic", size = 10) 
  )

print(p_boxplot)

output_dir <- "Figuras/Figura_4/" 
output_filename <- "Fig_4A.svg"

ggsave(
  filename = file.path(output_dir, output_filename), 
  plot = p_boxplot, 
  width = 7.5,
  height = 3,   
  units = "in"
)

features_ordered <- c("Rpl22", "Rpl22l1", "Rpl7", "Rpl7l1", "Rps27", "Rps27l")
DefaultAssay(seurat_merged) <- "RNA"   
dot_plot_agregado <- DotPlot(
  seurat_merged,
  features = features_ordered, 
  group.by = "subclass_label",
  scale = FALSE, dot.min = 0.25, scale.min = 65, col.max = 1.25, 
  dot.scale = 2
) + 
  labs(
    y = "Neuronal Subclass", 
    x = NULL,
    color = "Average Expression",     
    size = "Percentage Expression" 
  ) + 
  scale_colour_gradient2(
    low = "blue", 
    mid = "lightgray", 
    high = "red", 
    midpoint = 0.75,
    limits = c(NA, 1.0), 
    oob = scales::squish 
  ) +
  guides(
    colour = guide_colorbar(
      barheight = unit(2.5, "lines"), 
      barwidth = unit(0.5, "lines")  
    ),
    size = guide_legend(
      keyheight = unit(0.5, "lines") 
    )
  ) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size = 7), 
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 7), 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 7)  
  )

print(dot_plot_agregado)

output_dir <- "Figuras/Figura_4/" 
output_filename <- "Fig_4B.svg"

ggsave(
  filename = file.path(output_dir, output_filename), 
  plot = dot_plot_agregado, 
  width = 7.5, 
  height = 2.5,  
  units = "in"
)


output_dir <- "Figuras/Figura_S4/DotPlots_Individuales_SinEscalar_SVG"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

geneswithparalogs <- c("Rpl7", "Rpl7l1", "Rpl22", "Rpl22l1", "Rps27", "Rps27l")

dotplots_sin_leyenda <- list()
for (i in seq_along(seurat_split)) {
  dot_plot <- DotPlot(
    seurat_split[[i]],
    features = geneswithparalogs,
    group.by = "subclass_label",
    scale = FALSE, dot.min = 0.25, scale.min = 65, col.max = 1.25,
    dot.scale = 4
  ) +
    labs(title = regions_pass_DESeq2[i]) +
    scale_colour_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 0.75) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
  dotplots_sin_leyenda[[i]] <- dot_plot
}

legend_plot <- DotPlot(
  seurat_split[[1]],
  features = geneswithparalogs,
  group.by = "subclass_label",
  scale = FALSE, dot.min = 0.25, scale.min = 65, col.max = 1.25,
  dot.scale = 4
) +
  scale_colour_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 0.75) +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )

legend <- cowplot::get_legend(legend_plot)
ncol <- 3
n_plots <- length(dotplots_sin_leyenda)
total_slots <- ceiling((n_plots + 1) / ncol) * ncol
n_fillers <- total_slots - (n_plots + 1)
dotplots_con_leyenda <- c(dotplots_sin_leyenda, rep(list(NULL), n_fillers), list(legend))
final_plot <- cowplot::plot_grid(plotlist = dotplots_con_leyenda, ncol = ncol)

ggsave(
  filename = file.path(output_dir, "DotPlots_Combined_LegendRight.svg"),
  plot = final_plot,
  width = 7.5,
  height = ceiling(total_slots / ncol) * 3.5,
  units = "in"
)

