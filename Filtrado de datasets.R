# =============================================================================
# PASO 1: Leer y organizar los 16 archivos de expresión de PEP1 (raíz At)
# =============================================================================
library(readxl)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Define la ruta donde están los archivos xlsx
# -----------------------------------------------------------------------------
ruta_archivos <- "PEP1 raíz At/"   # <-- MODIFICA ESTO

# -----------------------------------------------------------------------------
# 2. Tabla de metadatos
# -----------------------------------------------------------------------------
metadatos <- data.frame(
  archivo = c(
    "GSM8001046_Mock1_24H.xlsx",
    "GSM8001047_PEP1_24H_1.xlsx",
    "GSM8001048_PEP1_24H_2.xlsx",
    "GSM8001049_PEP1_24H_3.xlsx",
    "GSM8001050_Mock2_12H.xlsx",
    "GSM8001051_PEP1_12H_1.xlsx",
    "GSM8001052_PEP1_12H_2.xlsx",
    "GSM8001053_PEP1_12H_3.xlsx",
    "GSM8001054_Mock3_3H.xlsx",
    "GSM8001055_PEP1_3H_1.xlsx",
    "GSM8001056_PEP1_3H_2.xlsx",
    "GSM8001057_PEP1_3H_3.xlsx",
    "GSM8001058_Mock4_6H.xlsx",
    "GSM8001059_PEP1_6H_1.xlsx",
    "GSM8001060_PEP1_6H_2.xlsx",
    "GSM8001061_PEP1_6H_3.xlsx"
  ),
  muestra = c(
    "Mock_24H",
    "PEP1_24H_rep1", "PEP1_24H_rep2", "PEP1_24H_rep3",
    "Mock_12H",
    "PEP1_12H_rep1", "PEP1_12H_rep2", "PEP1_12H_rep3",
    "Mock_3H",
    "PEP1_3H_rep1",  "PEP1_3H_rep2",  "PEP1_3H_rep3",
    "Mock_6H",
    "PEP1_6H_rep1",  "PEP1_6H_rep2",  "PEP1_6H_rep3"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 3. Función para leer cada archivo
# -----------------------------------------------------------------------------
leer_muestra <- function(archivo, nombre_muestra, ruta) {
  ruta_completa <- file.path(ruta, archivo)
  
  df <- read_excel(
    ruta_completa,
    col_types = c(
      "text",    # Name
      "numeric", # Chromosome
      "text",    # Region
      "numeric", # Expression value
      "numeric", # TPM
      "numeric", # RPKM
      "text",    # Gene ID
      "numeric", # Gene length
      "text",    # Biotype
      "numeric", # Unique gene reads
      "numeric"  # Total gene reads
    )
  )
  
  df_limpio <- df %>%
    select(Gene_ID = `Gene ID`, !!nombre_muestra := TPM)
  
  return(df_limpio)
}

# -----------------------------------------------------------------------------
# 4. Leer todos los archivos y unirlos por Gene ID
# -----------------------------------------------------------------------------
cat("Leyendo archivos...\n")

lista_dfs <- mapply(
  FUN      = leer_muestra,
  archivo  = metadatos$archivo,
  nombre_muestra = metadatos$muestra,
  MoreArgs = list(ruta = ruta_archivos),
  SIMPLIFY = FALSE
)

datos_completos <- Reduce(function(df1, df2) full_join(df1, df2, by = "Gene_ID"),
                          lista_dfs)

cat("✓ Datos unidos correctamente.\n")
cat("  Genes totales:  ", nrow(datos_completos), "\n")
cat("  Columnas:       ", ncol(datos_completos), "\n")

# -----------------------------------------------------------------------------
# 5. Reordenar columnas por tiempo (3H → 6H → 12H → 24H)
# -----------------------------------------------------------------------------
orden_columnas <- c(
  "Gene_ID",
  "Mock_3H",  "PEP1_3H_rep1",  "PEP1_3H_rep2",  "PEP1_3H_rep3",
  "Mock_6H",  "PEP1_6H_rep1",  "PEP1_6H_rep2",  "PEP1_6H_rep3",
  "Mock_12H", "PEP1_12H_rep1", "PEP1_12H_rep2", "PEP1_12H_rep3",
  "Mock_24H", "PEP1_24H_rep1", "PEP1_24H_rep2", "PEP1_24H_rep3"
)

datos_completos <- datos_completos[, orden_columnas]

# -----------------------------------------------------------------------------
# 6. Vista previa
# -----------------------------------------------------------------------------
cat("\nPrimeras filas del dataset organizado:\n")
print(head(datos_completos))

cat("\nEstructura:\n")
glimpse(datos_completos)

# -----------------------------------------------------------------------------
# 7. Guardar CSV completo
# -----------------------------------------------------------------------------
write.csv(datos_completos, "PEP1_raiz_TPM_organizado.csv", row.names = FALSE)
cat("\n✓ Archivo guardado: PEP1_raiz_TPM_organizado.csv\n")


# =============================================================================
# PASO 1B: Filtrar genes de interés
# =============================================================================

# -----------------------------------------------------------------------------
# 8. Lista de nombres simbólicos de los genes de interés
#    NOTA: Gene_ID en tus archivos puede estar como símbolo (HRE1, ADH1...)
#    o como locus AGI (AT1G72360, etc.). Ajusta genes_interes según corresponda.
# -----------------------------------------------------------------------------
genes_interes <- c(
  "HRE1", "HRE2",
  "RAP2.12", "RAP2.2", "RAP2.3",
  "PCO2",
  "WRKY33",
  "ADH1",
  "PRT1", "PRT6",
  "PDF1.2",
  "NPR1", "PR1",
  "FRK1",
  "RBOHD"
)

# -----------------------------------------------------------------------------
# 9. Filtrar el dataset completo por los genes de interés
# -----------------------------------------------------------------------------
datos_genes_interes <- datos_completos %>%
  filter(Gene_ID %in% genes_interes)

# Reporte de cuántos se encontraron vs. cuántos se buscaron
genes_encontrados <- datos_genes_interes$Gene_ID
genes_no_encontrados <- setdiff(genes_interes, genes_encontrados)

cat("\n========================================\n")
cat("  FILTRO DE GENES DE INTERÉS\n")
cat("========================================\n")
cat("  Genes buscados:    ", length(genes_interes), "\n")
cat("  Genes encontrados: ", nrow(datos_genes_interes), "\n")

if (length(genes_no_encontrados) > 0) {
  cat("\n  ⚠ Genes NO encontrados en el dataset:\n")
  cat("   ", paste(genes_no_encontrados, collapse = ", "), "\n")
  cat("  → Verifica si en tus archivos usan locus AGI (ej. AT1G72360)\n")
  cat("    en lugar del nombre simbólico, y ajusta genes_interes.\n")
} else {
  cat("\n  ✓ Todos los genes fueron encontrados.\n")
}

cat("\nExpresión TPM de los genes de interés:\n")
print(datos_genes_interes)

# -----------------------------------------------------------------------------
# 10. Guardar CSV filtrado
# -----------------------------------------------------------------------------
write.csv(datos_genes_interes,
          "PEP1_raiz_TPM_genes_interes.csv",
          row.names = FALSE)

cat("\n✓ Archivo guardado: PEP1_raiz_TPM_genes_interes.csv\n")
cat("  → Usa este archivo para graficar o analizar solo los genes de interés.\n")
# =============================================================================
# PASO 2: Calcular log2(TPM+1), promediar réplicas PEP1 y exportar a Excel
# Input:  PEP1_raiz_TPM_genes_interes.csv
# Output: PEP1_raiz_genes_interes_promedio.xlsx
# =============================================================================

library(dplyr)
library(readr)
library(writexl)

# -----------------------------------------------------------------------------
# 1. Cargar el archivo filtrado
# -----------------------------------------------------------------------------
datos <- read_csv("PEP1_raiz_TPM_genes_interes.csv")

# -----------------------------------------------------------------------------
# 2. Transformar a log2(TPM + 1) en todas las columnas de expresión
# -----------------------------------------------------------------------------
cols_expresion <- colnames(datos)[!(colnames(datos) %in% c("Gene_name", "Gene_ID"))]

datos_log2 <- datos %>%
  mutate(across(all_of(cols_expresion), ~ log2(. + 1)))

cat("✓ Transformación log2(TPM+1) aplicada.\n")

# -----------------------------------------------------------------------------
# 3. Calcular promedio de las 3 réplicas PEP1 por tiempo
#    y mantener los valores Mock individuales
# -----------------------------------------------------------------------------
resultado <- datos_log2 %>%
  mutate(
    # Promedios PEP1
    PEP1_3H_mean  = rowMeans(select(., PEP1_3H_rep1,  PEP1_3H_rep2,  PEP1_3H_rep3)),
    PEP1_6H_mean  = rowMeans(select(., PEP1_6H_rep1,  PEP1_6H_rep2,  PEP1_6H_rep3)),
    PEP1_12H_mean = rowMeans(select(., PEP1_12H_rep1, PEP1_12H_rep2, PEP1_12H_rep3)),
    PEP1_24H_mean = rowMeans(select(., PEP1_24H_rep1, PEP1_24H_rep2, PEP1_24H_rep3))
  ) %>%
  # Seleccionar columnas finales: Gene_name, Gene_ID, Mocks y promedios PEP1
  select(
    Gene_name, Gene_ID,
    Mock_3H,  PEP1_3H_mean,
    Mock_6H,  PEP1_6H_mean,
    Mock_12H, PEP1_12H_mean,
    Mock_24H, PEP1_24H_mean
  )

cat("✓ Promedios calculados.\n")

# -----------------------------------------------------------------------------
# 4. Vista previa
# -----------------------------------------------------------------------------
cat("\nResultado final:\n")
print(resultado)

# -----------------------------------------------------------------------------
# 5. Exportar a Excel
# -----------------------------------------------------------------------------
write_xlsx(resultado, "PEP1_raiz_genes_interes_promedio.xlsx")

cat("\n✓ Archivo guardado: PEP1_raiz_genes_interes_promedio.xlsx\n")
cat("  Columnas: Gene_name | Gene_ID | Mock_3H | PEP1_3H_mean | Mock_6H | PEP1_6H_mean | Mock_12H | PEP1_12H_mean | Mock_24H | PEP1_24H_mean\n")
