library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

load("data/Smartseq2dataset_seurat_filtered.RData")
p <- VlnPlot(seurat_merged, features = "Global RP expression", group.by = "subclass_label", pt.size = 0) + 
  labs(title = NULL, y = "Aggregated RP Expression") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10), 
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )
print(p)
dir.create("Figuras/Figura_1", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = "Figuras/Figura_1/Fig1_G.svg",
  plot = p,
  width = 3.4,
  height = 2.5,
  units = "in"
)
p <- VlnPlot(seurat_merged, features = "nCount_RNA", group.by = "subclass_label", pt.size = 0, y.max = 3E6) + 
  labs(title = NULL, y = "nCounts") + 
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10), 
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) 
  )
print(p)
dir.create("Figuras/Figura_S1", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = "Figuras/Figura_S1/Fig_S1A.svg", 
  plot = p,
  width = 3.75,  
  height = 2.5,   
  units = "in"
)
p <- VlnPlot(seurat_merged, features = "nFeature_RNA", group.by = "subclass_label", pt.size = 0) + 
  labs(title = NULL, y = "Detected Genes") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10), 
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) 
  )
print(p)
ggsave(
  filename = "Figuras/Figura_S1/Fig_S1B.svg", 
  plot = p,
  width = 3.75, 
  height = 2.5,
  units = "in"
)

data <- FetchData(seurat_merged, vars = c("Detected_RP", "subclass_label"))
data_counts <- data %>%
  group_by(subclass_label, Detected_RP) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(subclass_label) %>%
  mutate(total = sum(count),  
         freq = count / total)  

color_palette_func <- colorRampPalette(c("#440154FF", "#21908CFF", "#FDE725FF"))
discrete_colors <- color_palette_func(85)
names(discrete_colors) <- 0:84
main_plot <- ggplot(data_counts, aes(x = subclass_label, y = freq, fill = as.factor(Detected_RP))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = discrete_colors) +
  labs(
    x = "Identity",
    y = "Proportion of Cells"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,  
      size = 8
    )
  )
legend_data <- data.frame(y = 0:84)
legend_plot <- ggplot(legend_data, aes(x = 1, y = y, fill = as.factor(y))) +
  geom_raster() + 
  scale_fill_manual(values = discrete_colors) +
  scale_y_continuous(breaks = c(0, 84), expand = c(0, 0)) +
  labs(y = "Detected RP") +
  theme_void() + 
  theme(
    axis.text.y = element_text(size = 8, hjust = 0),
    axis.title.y = element_text(
      angle = 90, 
      hjust = 0.5, 
      size = 10 
    ),
    legend.position = "none" 
  )
final_plot <- main_plot + legend_plot + plot_layout(widths = c(20, 1))
print(final_plot)
ggsave(
  filename = "Figuras/Figura_1/Fig1_E.svg",
  plot = final_plot,
  width = 3.75, 
  height = 2.5,
  units = "in"
)

p1 <- DimPlot(seurat_merged, group.by = "subclass_label", raster = TRUE) +
  labs(
    title = NULL,
    color = "Identity"
  ) +
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 2)
  )) +
  theme_classic() +
  theme(
    aspect.ratio = 1, 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"), 
    legend.spacing.y = unit(0.05, "cm"), 
    legend.margin = margin(t = 0, b = 0) 
  )
print(p1)
ggsave(
  filename = "Figuras/Figura_1/Fig1_C.svg",
  plot = p1,
  width = 3.75, 
  height = 2.5, 
  units = "in"
)
p2 <- DimPlot(seurat_merged, group.by = "region_label", raster = TRUE) +
  labs(
    title = NULL,         
    color = "Brain Region" 
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 2) 
  )) +
  theme_classic() +
  theme(
    aspect.ratio = 1, 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"), 
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(t = 0, b = 0) 
  )
print(p2)
ggsave(
  filename = "Figuras/Figura_1/Fig1_A.svg", 
  plot = p2,
  width = 3.75, 
  height = 2, 
  units = "in"
)

#10x
load("data/seurat_merged_filtered_final.RData")

p <- VlnPlot(seurat_merged, features = "RibosomalSum", group.by = "subclass_label", pt.size = 0) + 
  labs(title = NULL, y = "Aggregated RP Expression") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10), 
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

print(p)

ggsave(
  filename = "Figuras/Figura_1/Fig1_H.svg",
  plot = p,
  width = 3.4,
  height = 2.5,
  units = "in"
)


p <- VlnPlot(seurat_merged, features = "nCount_RNA", group.by = "subclass_label", pt.size = 0, y.max = 60000)+
 
  labs(title = NULL, y = "nUMI") + 
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10), 
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) 
  )

print(p)

ggsave(
  filename = "Figuras/Figura_S1/Fig_S1C.svg", 
  plot = p,
  width = 3.75,  
  height = 2.5,   
  units = "in"
)

p <- VlnPlot(seurat_merged, features = "nFeature_RNA", group.by = "subclass_label", pt.size = 0) + 
  labs(title = NULL, y = "Detected Genes") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10), 
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) 
  )
print(p)


ggsave(
  filename = "Figuras/Figura_S1/Fig_S1D.svg", 
  plot = p,
  width = 3.75,  # <-- Coincide con las anteriores
  height = 2.5,   # <-- Coincide con las anteriores
  units = "in"
)

data <- FetchData(seurat_merged, vars = c("Detected_RP", "subclass_label"))
data_counts <- data %>%
  group_by(subclass_label, Detected_RP) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(subclass_label) %>%
  mutate(total = sum(count), 
         freq = count / total)  

color_palette_func <- colorRampPalette(c("#440154FF", "#21908CFF", "#FDE725FF"))
discrete_colors <- color_palette_func(85)
names(discrete_colors) <- 0:84

main_plot <- ggplot(data_counts, aes(x = subclass_label, y = freq, fill = as.factor(Detected_RP))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = discrete_colors) +
  labs(
    x = "Identity",
    y = "Proportion of Cells"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,  
      size = 8
    )
  )
legend_data <- data.frame(y = 0:84)

legend_plot <- ggplot(legend_data, aes(x = 1, y = y, fill = as.factor(y))) +
  geom_raster() + 
  scale_fill_manual(values = discrete_colors) +
  scale_y_continuous(breaks = c(0, 84), expand = c(0, 0)) +
  labs(y = "Detected RP") +
  theme_void() + 
  theme(
    axis.text.y = element_text(size = 8, hjust = 0), 
    axis.title.y = element_text(
      angle = 90, 
      hjust = 0.5, 
      size = 10 
    ),
    legend.position = "none" 
  )
final_plot <- main_plot + legend_plot + plot_layout(widths = c(20, 1))
print(final_plot)

ggsave(
  filename = "Figuras/Figura_1/Fig1_F.svg",
  plot = final_plot,
  width = 3.75,  
  height = 2.5,   
  units = "in"
)

p1 <- DimPlot(seurat_merged, group.by = "subclass_label", raster = TRUE) +
  labs(
    title = NULL,
    color = "Identity"
  ) +
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 2)
  )) +
  
  theme_classic() +
  theme(
    aspect.ratio = 1, 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"), 
    legend.spacing.y = unit(0.05, "cm"), 
    legend.margin = margin(t = 0, b = 0) 
  )

print(p1)
ggsave(
  filename = "Figuras/Figura_1/Fig1_D.svg",
  plot = p1,
  width = 3.75, 
  height = 2.5, 
  units = "in"
)
p2 <- DimPlot(seurat_merged, group.by = "region_label", raster = TRUE) +
  labs(
    title = NULL,         
    color = "Brain Region" 
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 2) 
  )) +
  theme_classic() +
  theme(
    aspect.ratio = 1, 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"), 
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(t = 0, b = 0) 
  )

print(p2)
ggsave(
  filename = "Figuras/Figura_1/Fig1_B.svg", 
  plot = p2,
  width = 3.75, 
  height = 2, 
  units = "in"
)
