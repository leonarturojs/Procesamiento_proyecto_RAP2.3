# ═══════════════════════════════════════════════════════════════════════════
# HEATMAP ERF-VII - Arabidopsis PEP1 vs Mock
# Dataset: ERFVII_PEP1_TPM_Zscore.xlsx
# Orden: Control 3H→24H, luego PEP1 3H→24H (R1,R2,R3)
# Etiquetas: multilinea verticales (angle_col = 0)
# Escala: Z-score clampeado a [-2, +2]
# Sin barras de anotacion superiores
# ═══════════════════════════════════════════════════════════════════════════

library(pheatmap)
library(RColorBrewer)
library(openxlsx)

# =============================================================================
# PASO 1: CARGAR DATOS
# =============================================================================

ruta <- "ERFVII_PEP1_TPM_Zscore.xlsx"   # <-- CAMBIA ESTO

zscore_df <- read.xlsx(ruta, sheet = "Zscore")
metadata  <- read.xlsx(ruta, sheet = "Metadata")

rownames(zscore_df) <- zscore_df$Gen
zscore_df           <- zscore_df[, -1]

cat("Genes:", nrow(zscore_df), "\n")
cat("Muestras:", ncol(zscore_df), "\n\n")

# =============================================================================
# PASO 2: ORDEN DE COLUMNAS — Control primero (3→24H), luego PEP1 (3→24H, R1R2R3)
# =============================================================================

tiempo_order <- c("3H", "6H", "12H", "24H")

mock_samples <- metadata[metadata$condicion == "Mock", ]
pep1_samples <- metadata[metadata$condicion == "PEP1", ]

mock_samples <- mock_samples[order(match(mock_samples$tiempo, tiempo_order)), ]
pep1_samples <- pep1_samples[order(match(pep1_samples$tiempo, tiempo_order)), ]

col_order <- c(mock_samples$label, pep1_samples$label)

zscore_matrix           <- as.matrix(zscore_df[, col_order])
rownames(zscore_matrix) <- rownames(zscore_df)

# =============================================================================
# PASO 3: CLAMP Z-score a [-2, +2]
# =============================================================================

zscore_matrix <- pmax(pmin(zscore_matrix, 2), -2)

# =============================================================================
# PASO 4: ETIQUETAS MULTILINEA VERTICALES
# Control\n3H        PEP1\n3H\nR1
# =============================================================================

conds_ordered  <- metadata$condicion[match(col_order, metadata$label)]
tiempos_ordered <- metadata$tiempo[match(col_order, metadata$label)]

col_labels <- mapply(function(cond, tiempo, label) {
  if (cond == "Mock") {
    paste0("Control\n", tiempo)
  } else {
    rep_n <- sub(".*_R(\\d+)$", "R\\1", label)
    paste0("PEP1\n", tiempo, "\n", rep_n)
  }
}, cond = conds_ordered, tiempo = tiempos_ordered, label = col_order)

# =============================================================================
# PASO 5: PALETA DE COLORES
# ↓↓↓ MODIFICA AQUÍ TUS COLORES ↓↓↓
# =============================================================================

color_negativo <- "#1a62a8"   # valores bajos  (azul)
color_medio    <- "#FFFFFF"   # valor cero     (blanco)
color_positivo <- "#ae1024"   # valores altos  (rojo)

heat_pal <- colorRampPalette(c(color_negativo, color_medio, color_positivo))(101)

# =============================================================================
# PASO 6: GAPS
# Control(n_mock cols) | PEP1·3H(3) | PEP1·6H(3) | PEP1·12H(3) | PEP1·24H(3)
# =============================================================================

n_mock <- nrow(mock_samples)   # normalmente 4
gaps_vals <- c(
  n_mock,
  n_mock + 3,
  n_mock + 6,
  n_mock + 9
)

# =============================================================================
# PASO 7: GUARDAR HEATMAP
# =============================================================================

output_dir <- file.path(getwd(), "ERF-VII_PEP1_Heatmaps")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_genes   <- nrow(zscore_matrix)
n_samples <- ncol(zscore_matrix)

filename <- file.path(output_dir, "ERF-VII_PEP1_vs_Control.png")
png(filename, width = max(1400, 800 + n_samples * 120),
    height = max(900, 600 + n_genes * 60),
    res = 150)

pheatmap(
  zscore_matrix,
  color          = heat_pal,
  breaks         = seq(-2, 2, length.out = 102),
  scale          = "none",
  cluster_cols   = FALSE,
  cluster_rows   = FALSE,
  annotation_col = NA,              # sin barras de color superiores
  labels_col     = col_labels,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize       = 13,
  fontsize_row   = 13,
  fontsize_col   = max(8, min(12, 180 / n_samples)),
  border_color   = "grey90",
  main           = paste0(
    "ERF-VII — Arabidopsis thaliana PEP1 vs Control\n",
    "Z-score de log2(TPM+1) | ", n_genes, " genes | ", n_samples, " muestras\n",
    "Control (3H->24H) -> PEP1 (3H->24H)"
  ),
  cellwidth      = 65,
  cellheight     = 40,
  angle_col      = 0,               # etiquetas verticales sin inclinar
  legend_breaks  = c(-2, -1, 0, 1, 2),
  legend_labels  = c("-2", "-1", "0", "+1", "+2"),
  gaps_col       = gaps_vals,
  gaps_row       = 1
)

dev.off()
cat("Heatmap guardado:", filename, "\n")