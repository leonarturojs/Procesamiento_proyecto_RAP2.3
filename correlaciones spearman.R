# ============================================================
#  Correlaciones de Spearman: RAP2.3 (Glyma.09G041500.1)
#  vs TODOS los genes — Dataset Gmax waterlogging (TPM)
#  24 muestras: Control (3h, 6h, 12h, 24h x R1-R3) +
#               Hipoxia (3h, 6h, 12h, 24h x R1-R3)
# ============================================================

# --- Paquetes ---
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

tsv_file <- "Gmax_waterlogging_mRNA_expression_TPM (1).tsv"
out_file <- "RAP2.3_Gmax_spearman_all_genes.csv"

cat("Directorio de trabajo:", getwd(), "\n")
cat("Archivo de entrada existe:", file.exists(tsv_file), "\n\n")

# ============================================================
# 1. Leer el archivo TSV
#    El header tiene saltos de línea dentro de las celdas
#    (formato especial con comillas), por eso se parsea manualmente
# ============================================================

cat("Leyendo archivo...\n")
raw_lines <- readLines(tsv_file)

# El header ocupa las primeras 25 líneas:
# Muestras 1-12  = Control  (3h/R1 ... 24h/R3)
# Muestras 13-24 = Hipoxia  (3h/R1 ... 24h/R3)
header_lines <- raw_lines[1:25]

col_suffixes <- sapply(header_lines[-1], function(x) {
  gsub('"', '', strsplit(x, "\t")[[1]][1])
})
prefixes  <- c(rep("Ctrl", 12), rep("Hipo", 12))
col_names <- c("Transcript_ID", paste(prefixes, col_suffixes, sep = "_"))

cat("Columnas detectadas:", length(col_names) - 1, "muestras\n")

# Leer datos desde línea 26
cat("Cargando matriz de expresión...\n")
df <- fread(
  text      = paste(raw_lines[26:length(raw_lines)], collapse = "\n"),
  sep       = "\t",
  header    = FALSE,
  col.names = col_names,
  na.strings = c("NA", "")
)
cat(sprintf("Dataset cargado: %d genes x %d muestras\n", nrow(df), ncol(df) - 1))

# ============================================================
# 2. Extraer vector de RAP2.3
# ============================================================

RAP23_ID <- "Glyma.09G041500.1"
rap_row  <- df[Transcript_ID == RAP23_ID]

if (nrow(rap_row) == 0) stop(paste("No se encontro el gen:", RAP23_ID))

expr_cols <- names(df)[2:ncol(df)]
rap_vals  <- as.numeric(rap_row[1, ..expr_cols])
cat(sprintf("RAP2.3 encontrado. TPM mediana = %.2f\n\n", median(rap_vals)))

# ============================================================
# 3. Filtrar genes con expresion nula
# ============================================================

cat("Filtrando genes sin expresion (todos TPM = 0)...\n")
row_sums    <- rowSums(df[, ..expr_cols])
df_filtered <- df[row_sums > 0]
cat(sprintf("Genes con expresion detectada: %d\n\n", nrow(df_filtered)))

# ============================================================
# 4. Calcular correlacion de Spearman vs RAP2.3
# ============================================================

cat("Calculando correlaciones de Spearman (aprox. 1-2 min)...\n")

gene_matrix <- as.matrix(df_filtered[, ..expr_cols])
rownames(gene_matrix) <- df_filtered$Transcript_ID

n         <- nrow(gene_matrix)
rho_vec   <- numeric(n)
pval_vec  <- numeric(n)
rap_rank  <- rank(rap_vals)   # rangos de RAP2.3, calculados una sola vez

for (i in seq_len(n)) {
  gene_rank    <- rank(gene_matrix[i, ])
  r            <- cor(rap_rank, gene_rank)
  rho_vec[i]   <- r
  t_stat       <- r * sqrt(22 / (1 - r^2))
  pval_vec[i]  <- 2 * pt(-abs(t_stat), df = 22)
}

# Ajuste por multiples pruebas (Benjamini-Hochberg)
padj_vec <- p.adjust(pval_vec, method = "BH")

# ============================================================
# 5. Ensamblar y guardar tabla
# ============================================================

results <- data.frame(
  Transcript_ID = df_filtered$Transcript_ID,
  rho           = round(rho_vec,  6),
  pval          = round(pval_vec, 8),
  padj_BH       = round(padj_vec, 8),
  signif        = ifelse(pval_vec < 0.001, "***",
                         ifelse(pval_vec < 0.01,  "**",
                                ifelse(pval_vec < 0.05,  "*", "ns"))),
  stringsAsFactors = FALSE
)

results <- results[order(-results$rho), ]
rownames(results) <- NULL

# Guardar CSV
write.csv(results, file = out_file, row.names = FALSE)

# Confirmar que se guardó correctamente
if (file.exists(out_file)) {
  cat(sprintf("\n✔ Archivo guardado exitosamente:\n  %s/%s\n", getwd(), out_file))
  cat(sprintf("  Tamaño: %.1f MB\n", file.info(out_file)$size / 1e6))
} else {
  cat("\n✘ ERROR: el archivo no se pudo guardar. Verifica los permisos del directorio.\n")
}

# ============================================================
# 6. Resumen en consola
# ============================================================

cat("\n=== RESUMEN ===\n")
cat(sprintf("Total genes analizados:             %d\n", nrow(results)))
cat(sprintf("Positivos significativos (p<0.05):  %d\n", sum(results$pval < 0.05 & results$rho > 0)))
cat(sprintf("Negativos significativos (p<0.05):  %d\n", sum(results$pval < 0.05 & results$rho < 0)))
cat(sprintf("No significativos:                  %d\n", sum(results$pval >= 0.05)))

cat("\nTop 10 correlaciones POSITIVAS:\n")
print(head(results, 10))

cat("\nTop 10 correlaciones NEGATIVAS:\n")
print(tail(results, 10))