# =============================================================================
# Diagrama de Venn — Hipoxia vs Patógenos
# Posicionamiento geométrico estricto en cada zona
# =============================================================================

library(ggplot2)
library(ggforce)

# =============================================================================
# DATOS
# =============================================================================

hipoxia   <- c("PCO3", "ACBP6", "ADH1", "ACBP2", "PDC1", "HRE1", "HRE2")

patogenos <- c("RCD1", "WRKY33", "RAP2.3", "PR2", "RGA", "DOF3.4", "BAP1",
               "WRKY48", "WRKY7", "RPM1", "WRKY15", "PDF1.4", "WRKY57",
               "RAP2.1", "Vestitone\nreductase")

ambas     <- c("NPR1", "EDS1", "PAD4", "ACBP4", "MPK5", "TGA4", "JAR1",
               "ETR2", "ETR1", "EIN2", "EBF2", "EBF1", "EIN4", "WRKY28",
               "CBP60G", "RAP2.12", "RAP2.2", "BAD1", "NAC92")

col_hip  <- "#D6EAF8"
col_pat  <- "#FADBD8"
col_both <- "#EAD5E8"

# =============================================================================
# CENTROS Y RADIO
# cx_L = centro círculo hipoxia, cx_R = centro círculo patógenos
# Los círculos se solapan: distancia entre centros = 2 * cx_R
# La intersección vertical ocurre en x = 0
# =============================================================================

r    <- 2.8    # radio
cx_L <- -1.6   # centro izquierdo
cx_R <-  1.6   # centro derecho
# La intersección ocurre entre x_int_left y x_int_right:
# x donde círculo L intersecta eje Y: x_int = (r² - cx_R² + cx_L²)/(2*(cx_L-cx_R))
# Pero más simple: zona exclusiva L = x < 0, zona exclusiva R = x > 0
# La frontera vertical de cada zona la calculamos:
# Círculo L: (x - cx_L)² + y² = r²  → en y=0: x = cx_L ± r
# Intersección empieza donde los dos círculos se cruzan
# x_cross = (cx_L + cx_R)/2 = 0  (por simetría cuando cx_L = -cx_R)
# → zona exclusiva L: x < x_cross = 0 Y dentro del círculo L
# → zona exclusiva R: x > x_cross = 0 Y dentro del círculo R

# Función: genera una cuadrícula de puntos DENTRO de una zona
# zona "L_only": dentro círculo L Y fuera círculo R
# zona "R_only": dentro círculo R Y fuera círculo L
# zona "both":   dentro ambos círculos

in_circle <- function(x, y, cx, cy, r) {
  (x - cx)^2 + (y - cy)^2 <= r^2
}

# Genera candidatos y filtra los que caen en la zona correcta
sample_zone <- function(n, zone, cx_L, cx_R, r,
                        x_range, y_range, seed = 42) {
  set.seed(seed)
  pts <- data.frame(x = numeric(0), y = numeric(0))
  # Generar muchos candidatos en cuadrícula densa
  xs <- seq(x_range[1], x_range[2], length.out = 200)
  ys <- seq(y_range[1], y_range[2], length.out = 200)
  grid <- expand.grid(x = xs, y = ys)
  inL  <- in_circle(grid$x, grid$y, cx_L, 0, r)
  inR  <- in_circle(grid$x, grid$y, cx_R, 0, r)
  if (zone == "L_only") valid <- grid[inL & !inR, ]
  if (zone == "R_only") valid <- grid[inR & !inL, ]
  if (zone == "both")   valid <- grid[inL &  inR, ]
  # Seleccionar n puntos distribuidos (dividir en subregiones)
  if (nrow(valid) < n) stop("Zona muy pequeña")
  # Distribuir equidistantes mediante k-means
  set.seed(seed)
  km <- kmeans(valid, centers = n, nstart = 20, iter.max = 100)
  centers <- as.data.frame(km$centers)
  names(centers) <- c("x","y")
  return(centers)
}

# --- Zona exclusiva Hipoxia (izquierda, fuera del círculo R) ---
df_hip_pts <- sample_zone(length(hipoxia), "L_only",
                          cx_L, cx_R, r,
                          x_range = c(cx_L - r, cx_R - r - 0.05),
                          y_range = c(-r * 0.95, r * 0.95))
df_hip <- data.frame(label = hipoxia,
                     x     = df_hip_pts$x,
                     y     = df_hip_pts$y)

# --- Zona exclusiva Patógenos (derecha, fuera del círculo L) ---
df_pat_pts <- sample_zone(length(patogenos), "R_only",
                          cx_L, cx_R, r,
                          x_range = c(cx_L + r + 0.05, cx_R + r),
                          y_range = c(-r * 0.95, r * 0.95))
df_pat <- data.frame(label = patogenos,
                     x     = df_pat_pts$x,
                     y     = df_pat_pts$y)

# --- Zona de intersección ---
df_both_pts <- sample_zone(length(ambas), "both",
                           cx_L, cx_R, r,
                           x_range = c(-1.5, 1.5),
                           y_range = c(-r * 0.92, r * 0.92))
df_both <- data.frame(label = ambas,
                      x     = df_both_pts$x,
                      y     = df_both_pts$y)

# =============================================================================
# PLOT
# =============================================================================

p <- ggplot() +
  
  geom_circle(aes(x0 = cx_L, y0 = 0, r = r),
              fill = col_hip, color = NA, alpha = 0.7) +
  geom_circle(aes(x0 = cx_R, y0 = 0, r = r),
              fill = col_pat, color = NA, alpha = 0.7) +
  # Zona de intersección más saturada (mezcla visual)
  geom_circle(aes(x0 = cx_L, y0 = 0, r = r),
              fill = col_both, color = NA, alpha = 0.3) +
  geom_circle(aes(x0 = cx_R, y0 = 0, r = r),
              fill = col_both, color = NA, alpha = 0.3) +
  
  # Genes hipoxia
  geom_text(data = df_hip,
            aes(x = x, y = y, label = label),
            size = 3.5, fontface = "bold",
            color = "#1A5276", hjust = 0.5, lineheight = 0.85) +
  
  # Genes patógenos
  geom_text(data = df_pat,
            aes(x = x, y = y, label = label),
            size = 3.5, fontface = "bold",
            color = "#922B21", hjust = 0.5, lineheight = 0.85) +
  
  # Genes compartidos
  geom_text(data = df_both,
            aes(x = x, y = y, label = label),
            size = 3.2, fontface = "bold.italic",
            color = "#4A235A", hjust = 0.5, lineheight = 0.85) +
  
  # Títulos de zonas
  annotate("text", x = cx_L - 0.6, y = r + 0.35,
           label = paste0("Hipoxia\n(n = ", length(hipoxia), ")"),
           size = 5.5, fontface = "bold", color = "#1A5276", hjust = 0.5) +
  
  annotate("text", x = cx_R + 0.6, y = r + 0.35,
           label = paste0("Patógenos\n(n = ", length(patogenos), ")"),
           size = 5.5, fontface = "bold", color = "#922B21", hjust = 0.5) +
  
  annotate("text", x = 0, y = r + 0.35,
           label = paste0("Ambas\n(n = ", length(ambas), ")"),
           size = 5, fontface = "bold", color = "#4A235A", hjust = 0.5) +
  
  labs(title = "Genes asociados a Respuesta Hipóxica y/o Respuesta a Patógenos") +
  
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-3.5, 3.8)) +
  theme_void(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 15, hjust = 0.5,
                                   color = "#1B2631",
                                   margin = margin(b = 10)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(20, 20, 20, 20)
  )

print(p)

ggsave("Venn_Hipoxia_Patogenos.pdf", plot = p,
       width = 16, height = 9, device = cairo_pdf)

ggsave("Venn_Hipoxia_Patogenos.png", plot = p,
       width = 16, height = 9, dpi = 300)

message("Listo: Venn_Hipoxia_Patogenos.pdf / .png")