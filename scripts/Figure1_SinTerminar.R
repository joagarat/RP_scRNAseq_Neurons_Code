#Smartseq

library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(RColorBrewer)
library(svglite)
setwd("/mnt/Disk01/jgarat/Articulo_SingleCell/")
load("Smartseq/seurat_split_ClasificacionActual.RData")
regions = array()
for (i in seq(1,21)) {
  regions = append(regions, unique(seurat_split[[i]]$region_label))
}
regions = regions[-1]

#Analisis con seuratsplit
#Filter non-neuronal cells
for (i in seq(1,21)) {
  seurat_split[[i]] = subset(seurat_split[[i]], subset = class_label %in% c("Glutamatergic", "GABAergic"))
}
#Filter >10cells/cluster 
for (i in seq(1,21)) {
  df =as.data.frame(table(seurat_split[[i]]$cluster_label))
  df = df[df$Freq>10, ]
  seurat_split[[i]] = subset(seurat_split[[i]], subset = cluster_label %in% df$Var1)
}
remove(df)

RP = read.csv("/mnt/Disk01/jgarat/Maestria/SingleCell/wholecortexallen/RP.csv")
RP = as.character(RP$x)
RP =gsub("Rack1", "Gnb2l1", x = RP)
seurat_split[[1]]@meta.data
#igualando criterios de filtrado, me quedo con las regiones que pudieron ser analizadas por DESeq2

regions_pass_DESeq2 = c("VISp", "VIS", "SSs", "SSp", "RSPv", "RSP", "MOp", "HIP", "ALM", "AI", "ACA")
for (i in seq_along(seurat_split)) {
  nombre_nuevo <- seurat_split[[i]]$region_label[1] 
  names(seurat_split)[i] <- nombre_nuevo
}

seurat_split = seurat_split[c("VISp", "VIS", "SSs", "SSp", "RSPv", "RSP", "MOp", "HIP", "ALM", "AI", "ACA")]
regions = names(seurat_split)



#Numero de celulas y de ratones HACERLO DESPUES DE FILTRAR SUBCLASS_PASS
cells = 0
for (i in seq_along(seurat_split)) {
  cells = cells + length(colnames(seurat_split[[i]]))
}
nmice = 0
mice = array()
for (i in seq_along(seurat_split)) {
  mice = append(mice, unique(seurat_split[[i]]$donor_label))
}
mice = unique(mice)

subclass_pass =list()
for (z in seq(1,11))   { metadata = seurat_split[[z]]@meta.data
subclass_3 = unique(metadata$subclass_label)
subclass_nombres = gsub("/","_", subclass_3)
subclass_pass[[z]] = array()
for(i in seq(1,length(subclass_3))){
  rats1 = metadata[metadata$subclass_label==subclass_3[i], "donor_label"]
  rats1_stats = as.data.frame(table(as.data.frame(rats1)))
  rats1_stats_mayor_a10 =rats1_stats[rats1_stats$Freq>5,]
  if (length(rats1_stats_mayor_a10$rats1)>2) {
    subclass_pass[[z]] = append(subclass_pass[[z]], subclass_3[i])
  }
}
subclass_pass[[z]] = subclass_pass[[z]][-1]
}

#Filter subclasses with less than 2 mice
for (i in seq(1,11)) {
  seurat_split[[i]] = subset(seurat_split[[i]], subset = class_label %in% c("Glutamatergic", "GABAergic"))
}
for (i in seq(1,11)) {
  seurat_split[[i]] = subset(seurat_split[[i]], subset = subclass_label %in% subclass_pass[[i]])
}

subclass_pass = unique(unlist(subclass_pass))
length(subclass_pass)

clusters = array()
for (i in seq_along(seurat_split)) {
  clusters = append(clusters, unique(seurat_split[[i]]$cluster_label))
}
length(clusters)


#Filter >10cells/cluster 
for (i in seq(1,11)) {
  df =as.data.frame(table(seurat_split[[i]]$cluster_label))
  df = df[df$Freq>10, ]
  seurat_split[[i]] = subset(seurat_split[[i]], subset = cluster_label %in% df$Var1)
}

#Composicion de subclases por región
library(dittoSeq)
library(scCustomize)
for (i in names(seurat_split)) {
  seurat_split[[i]] <- RenameCells(seurat_split[[i]],
                                   add.cell.id = i)
}


seurat_merged <- merge(
  x = seurat_split[[1]], 
  y = seurat_split[2:length(seurat_split)], 
  add.cell.ids = regions
)

#QC
#grafica expresion de genes RP
seurat_merged = NormalizeData(seurat_merged)
referende_data <- as.data.frame(seurat_merged@assays$RNA$data)
ribosomal_data <- referende_data[RP, ]
cell_ribosomal_sum <- Matrix::colSums(ribosomal_data)
seurat_merged <- AddMetaData(object = seurat_merged, metadata = cell_ribosomal_sum, col.name = "Global RP expression")

# El VlnPlot base
p <- VlnPlot(seurat_merged, features = "Global RP expression", group.by = "subclass_label", pt.size = 0) + 
  
  # Le quitamos el título generado por defecto y ponemos el del eje Y
  labs(title = NULL, y = "Aggregated RP Expression") +
  
  # Aplicamos un tema limpio
  theme_classic() +
  
  # --- AQUÍ LA MAGIA ---
  theme(
    legend.position = "none",
    
    # Ajusta el título del eje Y (ej: "Aggregated RP Expression")
    axis.title.y = element_text(size = 10), 
    
    axis.title.x = element_text(size = 10),
    
    # Ajusta el texto de los números/marcas del eje Y
    axis.text.y = element_text(size = 8),
    
    # Ajusta el texto de las etiquetas del eje X (las rotadas)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) # <-- TAMAÑO AÑADIDO
  )

# Para visualizar el plot
print(p)

# Tu ggsave se mantiene igual.
# El tamaño del texto (en puntos) es relativo al tamaño de guardado (en pulgadas).
ggsave(
  filename = "Figuras/Figura_1/Fig1_G.svg",
  plot = p,
  width = 3.4,
  height = 2.5,
  units = "in"
)


#QC
datos_modificados <- apply(ribosomal_data, 2, function(x) ifelse(x > 0, 1, 0))
suma_columnas <- colSums(datos_modificados)
seurat_merged <- AddMetaData(object = seurat_merged, metadata = suma_columnas, col.name = "Detected_RP")

# El VlnPlot base
p <- VlnPlot(seurat_merged, features = "nCount_RNA", group.by = "subclass_label", pt.size = 0, y.max = 3E6) + 
  
  # Le quitamos el título generado por defecto y ponemos el del eje Y
  labs(title = NULL, y = "nCounts") + # <-- CAMBIADO COMO PEDISTE
  
  # Aplicamos un tema limpio
  theme_classic() +
  
  # --- COPIAMOS LA ESTÉTICA DE FIG1 G ---
  theme(
    legend.position = "none",
    
    # Ajusta el título del eje Y
    axis.title.y = element_text(size = 10), 
    
    # Ajusta el título del eje X
    axis.title.x = element_text(size = 10),
    
    # Ajusta el texto de los números/marcas del eje Y
    axis.text.y = element_text(size = 8),
    
    # Ajusta el texto de las etiquetas del eje X (las rotadas)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) 
  )

# Para visualizar el plot
print(p)

# --- AJUSTAMOS EL TAMAÑO DE GUARDADO PARA QUE COINCIDA CON FIG1 G ---
ggsave(
  filename = "Figuras/Figura_S1/Fig_S1A.svg", 
  plot = p,
  width = 3.75,  # <-- Coincide con Fig1 G
  height = 2.5,   # <-- Coincide con Fig1 G
  units = "in"
)

# El VlnPlot base
p <- VlnPlot(seurat_merged, features = "nFeature_RNA", group.by = "subclass_label", pt.size = 0) + 
  
  # Le quitamos el título generado por defecto
  labs(title = NULL, y = "Detected Genes") +
  
  # Aplicamos un tema limpio
  theme_classic() +
  
  # --- APLICAMOS LA MISMA ESTÉTICA QUE LAS ANTERIORES ---
  theme(
    legend.position = "none",
    
    # Ajusta el título del eje Y
    axis.title.y = element_text(size = 10), 
    
    # Ajusta el título del eje X
    axis.title.x = element_text(size = 10),
    
    # Ajusta el texto de los números/marcas del eje Y
    axis.text.y = element_text(size = 8),
    
    # Ajusta el texto de las etiquetas del eje X (las rotadas)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) 
  )

# Para visualizar el plot
print(p)

# --- AJUSTAMOS EL TAMAÑO DE GUARDADO PARA QUE COINCIDA ---
ggsave(
  filename = "Figuras/Figura_S1/Fig_S1B.svg", 
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
  mutate(total = sum(count),  # Suma total de células en el cluster
         freq = count / total)  # Frecuencia relativa de cada valor de mi_variable

# Cargar las librerías necesarias
library(ggplot2)
library(patchwork) # Para unir los gráficos


# --- PASO 1: Crear una paleta de 85 colores ---
# (Esta parte se mantiene igual)
color_palette_func <- colorRampPalette(c("#440154FF", "#21908CFF", "#FDE725FF"))
discrete_colors <- color_palette_func(85)
names(discrete_colors) <- 0:84

# --- PASO 2: Crear el gráfico principal (con la nueva estética) ---

main_plot <- ggplot(data_counts, aes(x = subclass_label, y = freq, fill = as.factor(Detected_RP))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = discrete_colors) +
  labs(
    x = "Identity",
    y = "Proportion of Cells"
  ) +
  theme_classic() +
  theme(
    # --- APLICAMOS LA ESTÉTICA ESTÁNDAR ---
    legend.position = "none",
    
    # Títulos de Ejes
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    
    # Etiquetas de Ejes
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,  
      size = 8  # <-- TAMAÑO AUMENTADO
    )
  )

# --- PASO 3: Crear la leyenda personalizada (con la nueva estética) ---

legend_data <- data.frame(y = 0:84)

legend_plot <- ggplot(legend_data, aes(x = 1, y = y, fill = as.factor(y))) +
  geom_raster() + 
  scale_fill_manual(values = discrete_colors) +
  scale_y_continuous(breaks = c(0, 84), expand = c(0, 0)) +
  labs(y = "Detected RP") +
  theme_void() + 
  theme(
    # --- APLICAMOS LA ESTÉTICA ESTÁNDAR A LA LEYENDA ---
    
    # Etiquetas (0 y 84)
    axis.text.y = element_text(size = 8, hjust = 0), # <-- TAMAÑO AUMENTADO
    
    # Título ("Detected RP")
    axis.title.y = element_text(
      angle = 90, 
      hjust = 0.5, 
      size = 10 # <-- TAMAÑO AUMENTADO
    ),
    legend.position = "none" 
  )

# --- PASO 4: Unir los gráficos y guardarlos (con el nuevo tamaño) ---

final_plot <- main_plot + legend_plot + plot_layout(widths = c(20, 1))
print(final_plot)

# Guardar el gráfico final
ggsave(
  filename = "Figuras/Figura_1/Fig1_E.svg",
  plot = final_plot,
  width = 3.75,  # <-- COINCIDE CON EL RESTO
  height = 2.5,   # <-- COINCIDE CON EL RESTO
  units = "in"
)


seurat_merged = SCTransform(seurat_merged)
seurat_merged = RunPCA(seurat_merged)
seurat_merged = RunUMAP(seurat_merged, dims = 1:50)



p1 <- DimPlot(seurat_merged, group.by = "subclass_label", raster = TRUE) +
  labs(
    title = NULL,
    color = "Identity"
  ) +
  
  # --- AQUÍ LA MODIFICACIÓN ---
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 2) # <-- FUERZA EL TAMAÑO
  )) +
  
  theme_classic() +
  theme(
    # (El resto de tu theme se mantiene igual)
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
# ... (tu ggsave se mantiene igual)
ggsave(
  filename = "Figuras/Figura_1/Fig1_C.svg",
  plot = p1,
  width = 3.75, 
  height = 2.5, 
  units = "in"
)


# --- Gráfico 2: UMAP por Región Cerebral ---

p2 <- DimPlot(seurat_merged, group.by = "region_label", raster = TRUE) +
  labs(
    title = NULL,         
    color = "Brain Region" 
  ) +
  
  # --- AQUÍ LA MODIFICACIÓN ---
  guides(color = guide_legend(
    override.aes = list(size = 2) # <-- FUERZA EL MISMO TAMAÑO
  )) +
  
  theme_classic() +
  theme(
    # (El resto de tu theme se mantiene igual)
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
# ... (tu ggsave se mantiene igual)
ggsave(
  filename = "Figuras/Figura_1/Fig1_A.svg", 
  plot = p2,
  width = 3.75, 
  height = 2, 
  units = "in"
)


# --- INICIO: Bloque de Análisis Estadístico LMM ---

# 1. Cargar las librerías necesarias para el modelo
# (Asegúrate de tenerlas instaladas: install.packages(c("lme4", "lmerTest", "emmeans")))
library(lme4)
library(emmeans)  # Para los tests post-hoc (comparaciones por pares)
library(dplyr)    # Para la manipulación de datos (ya lo usas)

# 2. Preparar los datos para el Pseudobulking
# Es CRUCIAL usar los conteos CRUDOS (raw counts), no los normalizados

# Obtenemos la matriz de CONTEOS CRUDOS del assay "RNA"
# (Tu VlnPlot usaba datos normalizados, para esto necesitamos los crudos)
counts_matrix <- GetAssayData(seurat_merged, assay = "RNA", slot = "counts")

# Verificamos qué genes de tu lista 'RP' están en la matriz
# (Esto es por seguridad, por si algún gen RP no pasó los filtros)
rp_genes_in_matrix <- intersect(RP, rownames(counts_matrix))

# Filtramos la matriz de conteos crudos solo para esos genes RP
rp_counts_matrix <- counts_matrix[rp_genes_in_matrix, ]

# Sumamos los conteos de RP crudos por CÉLULA
cell_rp_sum_raw <- Matrix::colSums(rp_counts_matrix)

# Obtenemos el metadata relevante de Seurat
metadata <- seurat_merged@meta.data

# Creamos un data.frame único para el pseudobulking
# 'nCount_RNA' ya contiene los conteos totales crudos por célula
data_for_pb <- data.frame(
  cell_id = rownames(metadata),
  subclass = metadata$subclass_label, # El factor que queremos testear
  donor = metadata$donor_label,     # El factor aleatorio (ratón)
  total_counts_raw = metadata$nCount_RNA,
  rp_counts_raw = cell_rp_sum_raw
)

# 3. Realizar el Pseudobulking y Normalización
cat("--- Realizando Pseudobulking ---\n")

data_pseudobulk <- data_for_pb %>%
  group_by(subclass, donor) %>%
  summarise(
    # Sumamos todos los conteos de RP de ESE ratón en ESE clúster
    Suma_RP_Raw = sum(as.numeric(rp_counts_raw), na.rm = TRUE), 
    # Sumamos todos los conteos totales de ESE ratón en ESE clúster
    Suma_Total_Raw = sum(as.numeric(total_counts_raw), na.rm = TRUE),
    .groups = 'drop' # Desagrupamos la tabla
  )

# Calculamos log2(CPM+1) sobre los datos agregados
data_pseudobulk <- data_pseudobulk %>%
  mutate(
    RP_CPM = (Suma_RP_Raw / Suma_Total_Raw) * 1e6,
    log2_RP_CPM = log2(RP_CPM + 1)
  )

cat("Datos agregados listos para el modelo:\n")
print(head(data_pseudobulk))

# 4. Ejecutar el Modelo Lineal Mixto (LMM)
cat("\n--- Ejecutando Modelo Lineal Mixto (LMM) ---\n")

# Aseguramos que R entienda que son variables categóricas
data_pseudobulk$subclass <- as.factor(data_pseudobulk$subclass)
data_pseudobulk$donor <- as.factor(data_pseudobulk$donor)

# Corremos el modelo
# Queremos modelar la expresión de RP (log2_RP_CPM)
# en función de la subclase (subclass)
# PERO controlando por la variación inicial de cada ratón (1 | donor)
modelo_lmm <- lmer(log2_RP_CPM ~ subclass + (1 | donor), data = data_pseudobulk)

# Vemos un resumen del modelo.
# Es útil para ver cuánta varianza es explicada por el ratón (donor)
cat("Resumen del Modelo LMM:\n")
print(summary(modelo_lmm))

# Un test ANOVA sobre el modelo nos dice si 'subclass' es significativo EN GENERAL
cat("\nTest ANOVA para el efecto global de 'subclass':\n")
print(anova(modelo_lmm))

# 5. Obtener comparaciones por pares (Fold Change y p-adj)
cat("\n--- Calculando Comparaciones por Pares (Post-Hoc) ---\n")

# Calculamos las medias marginales estimadas (los valores promedio por subclase)
medias_marginales <- emmeans(modelo_lmm, ~ subclass)

# Realizamos TODAS las comparaciones por pares, con corrección Tukey para múltiples tests
comparaciones_pares <- pairs(medias_marginales, adjust = "tukey")

# Convertimos a data.frame para ver/guardar fácilmente
resultados_finales <- as.data.frame(comparaciones_pares)

cat("Resultados Finales (Fold Change y p-adj):\n")
print(resultados_finales)


ggplot(data_pseudobulk, aes(x = subclass, y = log2_RP_CPM)) +
  # El boxplot muestra la distribución (mediana, cuartiles)
  geom_boxplot(outlier.shape = NA) + 
  
  # El jitter muestra los puntos reales (cada punto es un ratón)
  # Le quitamos el 'color = donor' del 'aes()' para que no genere leyenda
  geom_jitter(width = 0.2, height = 0, alpha = 0.7, color = "black") + 
  
  theme_classic() +
  labs(
    title = "Métrica del Pseudobulk (Lo que el LMM está testeando)",
    x = "Subclase",
    y = "log2(Agregated RP CPM + 1)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none" # <-- Aquí está la solución
  )

# Puedes guardar esta tabla
 write.csv(resultados_finales, "Figuras/Supplementary_Tables_Extra/Supplementary_Table_LMM_RP_expression_comparisons_Smartseq.csv")

# --- FIN: Bloque de Análisis Estadístico LMM ---

# Un cambio del 50% es un Fold Change (FC) de 1.5
# El 'estimate' está en log2, así que calculamos el umbral:
log2fc_threshold <- log2(1.5) 

cat(paste0("\n--- Resultados Filtrados (p-adj < 0.05 Y |log2FC| > ", round(log2fc_threshold, 3), ") ---\n"))

resultados_filtrados <- resultados_finales %>%
  filter(
    p.value < 0.05,            # Filtro de p-valor ajustado (la columna 'p.value' ya está ajustada por Tukey)
    abs(estimate) > log2fc_threshold  # Filtro para >50% de cambio (abs() cubre subida o bajada)
  )

#10x

setwd("/mnt/Disk01/jgarat/Articulo_SingleCell/10x/")
library(rhdf5)

load("../Smartseq/seurat_split_ClasificacionActual.RData")

#Filtrar con RscriptSmartseq y hacer subclass_pass del smartseq
regions = array()
for (i in seq(1,21)) {
  regions = append(regions, unique(seurat_split[[i]]$region_label))
}
regions = regions[-1]

#Analisis con seuratsplit
#Filter non-neuronal cells
for (i in seq(1,21)) {
  seurat_split[[i]] = subset(seurat_split[[i]], subset = class_label %in% c("Glutamatergic", "GABAergic"))
}
#Filter >10cells/cluster 
for (i in seq(1,21)) {
  df =as.data.frame(table(seurat_split[[i]]$cluster_label))
  df = df[df$Freq>10, ]
  seurat_split[[i]] = subset(seurat_split[[i]], subset = cluster_label %in% df$Var1)
}
remove(df)

RP = read.csv("/mnt/Disk01/jgarat/Maestria/SingleCell/wholecortexallen/RP.csv")
RP = as.character(RP$x)
RP =gsub("Rack1", "Gnb2l1", x = RP)
seurat_split[[1]]@meta.data
#igualando criterios de filtrado, me quedo con las regiones que pudieron ser analizadas por DESeq2

regions_pass_DESeq2 = c("VISp", "VIS", "SSs", "SSp", "RSPv", "RSP", "MOp", "HIP", "ALM", "AI", "ACA")
for (i in seq_along(seurat_split)) {
  nombre_nuevo <- seurat_split[[i]]$region_label[1] 
  names(seurat_split)[i] <- nombre_nuevo
}

seurat_split = seurat_split[c("VISp", "VIS", "SSs", "SSp", "RSPv", "RSP", "MOp", "HIP", "ALM", "AI", "ACA")]
regions = names(seurat_split)



#Numero de celulas y de ratones HACERLO DESPUES DE FILTRAR SUBCLASS_PASS
cells = 0
for (i in seq_along(seurat_split)) {
  cells = cells + length(colnames(seurat_split[[i]]))
}
nmice = 0
mice = array()
for (i in seq_along(seurat_split)) {
  mice = append(mice, unique(seurat_split[[i]]$donor_label))
}
mice = unique(mice)

subclass_pass =list()
for (z in seq(1,11))   { metadata = seurat_split[[z]]@meta.data
subclass_3 = unique(metadata$subclass_label)
subclass_nombres = gsub("/","_", subclass_3)
subclass_pass[[z]] = array()
for(i in seq(1,length(subclass_3))){
  rats1 = metadata[metadata$subclass_label==subclass_3[i], "donor_label"]
  rats1_stats = as.data.frame(table(as.data.frame(rats1)))
  rats1_stats_mayor_a10 =rats1_stats[rats1_stats$Freq>5,]
  if (length(rats1_stats_mayor_a10$rats1)>2) {
    subclass_pass[[z]] = append(subclass_pass[[z]], subclass_3[i])
  }
}
subclass_pass[[z]] = subclass_pass[[z]][-1]
}

#Filter subclasses with less than 2 mice
for (i in seq(1,11)) {
  seurat_split[[i]] = subset(seurat_split[[i]], subset = class_label %in% c("Glutamatergic", "GABAergic"))
}
for (i in seq(1,11)) {
  seurat_split[[i]] = subset(seurat_split[[i]], subset = subclass_label %in% subclass_pass[[i]])
}

subclass_pass = unique(unlist(subclass_pass))
length(subclass_pass)

remove(seurat_split)

gc()


h5ls("expression_matrix.hdf5")
samples = h5read("expression_matrix.hdf5", "/data/samples")
genes = h5read("expression_matrix.hdf5", "/data/gene")
metadata = read.csv("CTX_Hip_anno_10x/CTX_Hip_anno_10x.csv")
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


subset_metadata <- metadata[metadata$subclass_label %in% subclass_pass, ]
nmice10x = unique(subset_metadata$donor_label)
length(nmice10x)
dim(subset_metadata)
nclusters10x = length(unique(subset_metadata$cluster_label))
#Subsample para poder hacer analisis global
metadata_subsampled = subset_metadata[sample(nrow(subset_metadata), 100000), ]
index = which(samples%in%metadata_subsampled$exp_component_name)
mtx = h5read("expression_matrix.hdf5", name = "data/counts", index = list(index, NULL))
mtx = as.data.frame(mtx)
mtx = t(mtx)
samplenames = samples[index]
colnames(mtx) = samplenames
rownames(mtx) = genes
mtx = as.data.frame(mtx)
library(Seurat)
library(ggplot2)
seurat = CreateSeuratObject(counts = mtx, meta.data = metadata_subsampled)
seurat_merged = seurat
seurat_merged = NormalizeData(seurat_merged)
#QC
RP = read.csv("/mnt/Disk01/jgarat/Maestria/SingleCell/wholecortexallen/RP.csv")
RP = as.character(RP$x)
#grafica expresion de genes RP
library(dittoSeq)
referende_data <- as.data.frame(seurat_merged@assays$RNA$data)
ribosomal_data <- referende_data[RP, ]
cell_ribosomal_sum <- Matrix::colSums(ribosomal_data)
seurat_merged <- AddMetaData(object = seurat_merged, metadata = cell_ribosomal_sum, col.name = "RibosomalSum")

# El VlnPlot base
p <- VlnPlot(seurat_merged, features = "RibosomalSum", group.by = "subclass_label", pt.size = 0) + 
  
  # Le quitamos el título generado por defecto y ponemos el del eje Y
  labs(title = NULL, y = "Aggregated RP Expression") +
  
  # Aplicamos un tema limpio
  theme_classic() +
  
  # --- AQUÍ LA MAGIA ---
  theme(
    legend.position = "none",
    
    # Ajusta el título del eje Y (ej: "Aggregated RP Expression")
    axis.title.y = element_text(size = 10), 
    
    axis.title.x = element_text(size = 10),
    
    # Ajusta el texto de los números/marcas del eje Y
    axis.text.y = element_text(size = 8),
    
    # Ajusta el texto de las etiquetas del eje X (las rotadas)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) # <-- TAMAÑO AÑADIDO
  )

# Para visualizar el plot
print(p)
setwd("/mnt/Disk01/jgarat/Articulo_SingleCell/")
# Tu ggsave se mantiene igual.
# El tamaño del texto (en puntos) es relativo al tamaño de guardado (en pulgadas).
ggsave(
  filename = "Figuras/Figura_1/Fig1_H.svg",
  plot = p,
  width = 3.4,
  height = 2.5,
  units = "in"
)


datos_modificados <- apply(ribosomal_data, 2, function(x) ifelse(x > 0, 1, 0))
suma_columnas <- colSums(datos_modificados)
seurat_merged <- AddMetaData(object = seurat_merged, metadata = suma_columnas, col.name = "Detected_RP")

# El VlnPlot base
p <- VlnPlot(seurat_merged, features = "nCount_RNA", group.by = "subclass_label", pt.size = 0, y.max = 60000)+
  # Le quitamos el título generado por defecto y ponemos el del eje Y
  labs(title = NULL, y = "nUMI") + # <-- CAMBIADO COMO PEDISTE
  
  # Aplicamos un tema limpio
  theme_classic() +
  
  # --- COPIAMOS LA ESTÉTICA DE FIG1 G ---
  theme(
    legend.position = "none",
    
    # Ajusta el título del eje Y
    axis.title.y = element_text(size = 10), 
    
    # Ajusta el título del eje X
    axis.title.x = element_text(size = 10),
    
    # Ajusta el texto de los números/marcas del eje Y
    axis.text.y = element_text(size = 8),
    
    # Ajusta el texto de las etiquetas del eje X (las rotadas)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) 
  )

# Para visualizar el plot
print(p)

# --- AJUSTAMOS EL TAMAÑO DE GUARDADO PARA QUE COINCIDA CON FIG1 G ---
ggsave(
  filename = "Figuras/Figura_S1/Fig_S1C.svg", 
  plot = p,
  width = 3.75,  # <-- Coincide con Fig1 G
  height = 2.5,   # <-- Coincide con Fig1 G
  units = "in"
)


# El VlnPlot base
p <- VlnPlot(seurat_merged, features = "nFeature_RNA", group.by = "subclass_label", pt.size = 0) + 
  
  # Le quitamos el título generado por defecto
  labs(title = NULL, y = "Detected Genes") +
  
  # Aplicamos un tema limpio
  theme_classic() +
  
  # --- APLICAMOS LA MISMA ESTÉTICA QUE LAS ANTERIORES ---
  theme(
    legend.position = "none",
    
    # Ajusta el título del eje Y
    axis.title.y = element_text(size = 10), 
    
    # Ajusta el título del eje X
    axis.title.x = element_text(size = 10),
    
    # Ajusta el texto de los números/marcas del eje Y
    axis.text.y = element_text(size = 8),
    
    # Ajusta el texto de las etiquetas del eje X (las rotadas)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8) 
  )

# Para visualizar el plot
print(p)

# --- AJUSTAMOS EL TAMAÑO DE GUARDADO PARA QUE COINCIDA ---
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
  mutate(total = sum(count),  # Suma total de células en el cluster
         freq = count / total)  # Frecuencia relativa de cada valor de mi_variable


# Cargar las librerías necesarias
library(ggplot2)
library(patchwork) # Para unir los gráficos


# --- PASO 1: Crear una paleta de 85 colores ---
# (Esta parte se mantiene igual)
color_palette_func <- colorRampPalette(c("#440154FF", "#21908CFF", "#FDE725FF"))
discrete_colors <- color_palette_func(85)
names(discrete_colors) <- 0:84

# --- PASO 2: Crear el gráfico principal (con la nueva estética) ---

main_plot <- ggplot(data_counts, aes(x = subclass_label, y = freq, fill = as.factor(Detected_RP))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = discrete_colors) +
  labs(
    x = "Identity",
    y = "Proportion of Cells"
  ) +
  theme_classic() +
  theme(
    # --- APLICAMOS LA ESTÉTICA ESTÁNDAR ---
    legend.position = "none",
    
    # Títulos de Ejes
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    
    # Etiquetas de Ejes
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,  
      size = 8  # <-- TAMAÑO AUMENTADO
    )
  )

# --- PASO 3: Crear la leyenda personalizada (con la nueva estética) ---

legend_data <- data.frame(y = 0:84)

legend_plot <- ggplot(legend_data, aes(x = 1, y = y, fill = as.factor(y))) +
  geom_raster() + 
  scale_fill_manual(values = discrete_colors) +
  scale_y_continuous(breaks = c(0, 84), expand = c(0, 0)) +
  labs(y = "Detected RP") +
  theme_void() + 
  theme(
    # --- APLICAMOS LA ESTÉTICA ESTÁNDAR A LA LEYENDA ---
    
    # Etiquetas (0 y 84)
    axis.text.y = element_text(size = 8, hjust = 0), # <-- TAMAÑO AUMENTADO
    
    # Título ("Detected RP")
    axis.title.y = element_text(
      angle = 90, 
      hjust = 0.5, 
      size = 10 # <-- TAMAÑO AUMENTADO
    ),
    legend.position = "none" 
  )

# --- PASO 4: Unir los gráficos y guardarlos (con el nuevo tamaño) ---

final_plot <- main_plot + legend_plot + plot_layout(widths = c(20, 1))
print(final_plot)

# Guardar el gráfico final
ggsave(
  filename = "Figuras/Figura_1/Fig1_F.svg",
  plot = final_plot,
  width = 3.75,  # <-- COINCIDE CON EL RESTO
  height = 2.5,   # <-- COINCIDE CON EL RESTO
  units = "in"
)

seurat_merged = SCTransform(seurat_merged)
seurat_merged = RunPCA(seurat_merged)
seurat_merged = RunUMAP(seurat_merged, dims = 1:50)



p1 <- DimPlot(seurat_merged, group.by = "subclass_label", raster = TRUE) +
  labs(
    title = NULL,
    color = "Identity"
  ) +
  
  # --- AQUÍ LA MODIFICACIÓN ---
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 2) # <-- FUERZA EL TAMAÑO
  )) +
  
  theme_classic() +
  theme(
    # (El resto de tu theme se mantiene igual)
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
# ... (tu ggsave se mantiene igual)
ggsave(
  filename = "Figuras/Figura_1/Fig1_D.svg",
  plot = p1,
  width = 3.75, 
  height = 2.5, 
  units = "in"
)


# --- Gráfico 2: UMAP por Región Cerebral ---

p2 <- DimPlot(seurat_merged, group.by = "region_label", raster = TRUE) +
  labs(
    title = NULL,         
    color = "Brain Region" 
  ) +
  
  # --- AQUÍ LA MODIFICACIÓN ---
  guides(color = guide_legend(
    override.aes = list(size = 2) # <-- FUERZA EL MISMO TAMAÑO
  )) +
  
  theme_classic() +
  theme(
    # (El resto de tu theme se mantiene igual)
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
# ... (tu ggsave se mantiene igual)
ggsave(
  filename = "Figuras/Figura_1/Fig1_B.svg", 
  plot = p2,
  width = 3.75, 
  height = 2, 
  units = "in"
)

# --- INICIO: Bloque de Análisis Estadístico LMM ---

# 1. Cargar las librerías necesarias para el modelo
# (Asegúrate de tenerlas instaladas: install.packages(c("lme4", "lmerTest", "emmeans")))
library(lme4)
library(emmeans)  # Para los tests post-hoc (comparaciones por pares)
library(dplyr)    # Para la manipulación de datos (ya lo usas)

# 2. Preparar los datos para el Pseudobulking
# Es CRUCIAL usar los conteos CRUDOS (raw counts), no los normalizados

# Obtenemos la matriz de CONTEOS CRUDOS del assay "RNA"
# (Tu VlnPlot usaba datos normalizados, para esto necesitamos los crudos)
counts_matrix <- GetAssayData(seurat_merged, assay = "RNA", slot = "counts")

# Verificamos qué genes de tu lista 'RP' están en la matriz
# (Esto es por seguridad, por si algún gen RP no pasó los filtros)
rp_genes_in_matrix <- intersect(RP, rownames(counts_matrix))

# Filtramos la matriz de conteos crudos solo para esos genes RP
rp_counts_matrix <- counts_matrix[rp_genes_in_matrix, ]

# Sumamos los conteos de RP crudos por CÉLULA
cell_rp_sum_raw <- Matrix::colSums(rp_counts_matrix)

# Obtenemos el metadata relevante de Seurat
metadata <- seurat_merged@meta.data

# Creamos un data.frame único para el pseudobulking
# 'nCount_RNA' ya contiene los conteos totales crudos por célula
data_for_pb <- data.frame(
  cell_id = rownames(metadata),
  subclass = metadata$subclass_label, # El factor que queremos testear
  donor = metadata$donor_label,     # El factor aleatorio (ratón)
  total_counts_raw = metadata$nCount_RNA,
  rp_counts_raw = cell_rp_sum_raw
)

# 3. Realizar el Pseudobulking y Normalización
cat("--- Realizando Pseudobulking ---\n")

data_pseudobulk <- data_for_pb %>%
  group_by(subclass, donor) %>%
  summarise(
    # Sumamos todos los conteos de RP de ESE ratón en ESE clúster
    Suma_RP_Raw = sum(as.numeric(rp_counts_raw), na.rm = TRUE), 
    # Sumamos todos los conteos totales de ESE ratón en ESE clúster
    Suma_Total_Raw = sum(as.numeric(total_counts_raw), na.rm = TRUE),
    .groups = 'drop' # Desagrupamos la tabla
  )

# Calculamos log2(CPM+1) sobre los datos agregados
data_pseudobulk <- data_pseudobulk %>%
  mutate(
    RP_CPM = (Suma_RP_Raw / Suma_Total_Raw) * 1e6,
    log2_RP_CPM = log2(RP_CPM + 1)
  )

cat("Datos agregados listos para el modelo:\n")
print(head(data_pseudobulk))

# 4. Ejecutar el Modelo Lineal Mixto (LMM)
cat("\n--- Ejecutando Modelo Lineal Mixto (LMM) ---\n")

# Aseguramos que R entienda que son variables categóricas
data_pseudobulk$subclass <- as.factor(data_pseudobulk$subclass)
data_pseudobulk$donor <- as.factor(data_pseudobulk$donor)

# Corremos el modelo
# Queremos modelar la expresión de RP (log2_RP_CPM)
# en función de la subclase (subclass)
# PERO controlando por la variación inicial de cada ratón (1 | donor)
modelo_lmm <- lmer(log2_RP_CPM ~ subclass + (1 | donor), data = data_pseudobulk)

# Vemos un resumen del modelo.
# Es útil para ver cuánta varianza es explicada por el ratón (donor)
cat("Resumen del Modelo LMM:\n")
print(summary(modelo_lmm))

# Un test ANOVA sobre el modelo nos dice si 'subclass' es significativo EN GENERAL
cat("\nTest ANOVA para el efecto global de 'subclass':\n")
print(anova(modelo_lmm))

# 5. Obtener comparaciones por pares (Fold Change y p-adj)
cat("\n--- Calculando Comparaciones por Pares (Post-Hoc) ---\n")

# Calculamos las medias marginales estimadas (los valores promedio por subclase)
medias_marginales <- emmeans(modelo_lmm, ~ subclass)

# Realizamos TODAS las comparaciones por pares, con corrección Tukey para múltiples tests
comparaciones_pares <- pairs(medias_marginales, adjust = "tukey")

# Convertimos a data.frame para ver/guardar fácilmente
resultados_finales <- as.data.frame(comparaciones_pares)

cat("Resultados Finales (Fold Change y p-adj):\n")
print(resultados_finales)


ggplot(data_pseudobulk, aes(x = subclass, y = log2_RP_CPM)) +
  # El boxplot muestra la distribución (mediana, cuartiles)
  geom_boxplot(outlier.shape = NA) + 
  
  # El jitter muestra los puntos reales (cada punto es un ratón)
  # Le quitamos el 'color = donor' del 'aes()' para que no genere leyenda
  geom_jitter(width = 0.2, height = 0, alpha = 0.7, color = "black") + 
  
  theme_classic() +
  labs(
    title = "Métrica del Pseudobulk (Lo que el LMM está testeando)",
    x = "Subclase",
    y = "log2(Agregated RP CPM + 1)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none" # <-- Aquí está la solución
  )

# Puedes guardar esta tabla
write.csv(resultados_finales, "Figuras/Supplementary_Tables_Extra/Supplementary_Table_LMM_RP_expression_comparisons_10XGenomics.csv")

# --- FIN: Bloque de Análisis Estadístico LMM ---

# Un cambio del 50% es un Fold Change (FC) de 1.5
# El 'estimate' está en log2, así que calculamos el umbral:
log2fc_threshold <- log2(1.5) 

cat(paste0("\n--- Resultados Filtrados (p-adj < 0.05 Y |log2FC| > ", round(log2fc_threshold, 3), ") ---\n"))

resultados_filtrados <- resultados_finales %>%
  filter(
    p.value < 0.05,            # Filtro de p-valor ajustado (la columna 'p.value' ya está ajustada por Tukey)
    abs(estimate) > log2fc_threshold  # Filtro para >50% de cambio (abs() cubre subida o bajada)
  )
