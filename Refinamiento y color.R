# =============================================================================
# ANÁLISIS COMPLETO DE PROTEÍNA
# InterPro (JSON) + AlphaFold3 (PDB) → Validación → PyMOL
# =============================================================================
# Paquetes necesarios:
#   install.packages(c("jsonlite", "tidyverse", "bio3d"))
# =============================================================================
 
library(jsonlite)
library(tidyverse)
library(bio3d)
 
# =============================================================================
# 0. CONFIGURACIÓN — CAMBIA SOLO ESTOS VALORES
# =============================================================================
 
ARCHIVO_JSON <- "JSONPvRAP23.json"
ARCHIVO_PDB  <- "fold_pv2_3_model_0.cif"
UMBRAL_PLDDT <- 70        # residuos con pLDDT < este valor se marcan en rojo
DIR_SALIDA   <- "./"      # carpeta donde se guardarán los archivos de salida
 
# =============================================================================
# 1. CARGAR Y PROCESAR JSON DE INTERPRO
# =============================================================================
 
cat("── [1/5] Cargando JSON de InterPro...\n")
 
data <- fromJSON(ARCHIVO_JSON, simplifyVector = FALSE)
 
extract_matches <- function(result_item) {
  map_df(result_item$matches, function(match) {
    nombre      <- match$signature$accession
    descripcion <- match$signature$name
    posiciones  <- map_df(match$locations, function(loc) {
      data.frame(start = loc$start, end = loc$end)
    })
    posiciones$dominio <- nombre
    posiciones$desc    <- descripcion
    return(posiciones)
  })
}
 
tabla_regiones <- extract_matches(data$results[[1]])
 
# Agrupar rangos del mismo dominio (evita que selecciones se pisen en PyMOL)
#colores_pymol <- c("salmon", "brown", "orange", "purple", "magenta",
                  # "red", "yellow", "cyan", "green", "blue", "pink", "lightgreen")
 
tabla_agrupada <- tabla_regiones %>%
  mutate(tag = gsub("[^[:alnum:]]", "_", dominio)) %>%
  group_by(dominio, tag, desc) %>%
  summarise(
    resi_str = paste0("resi ", start, "-", end, collapse = " + "),
    start_min = min(start),
    end_max   = max(end),
    .groups   = "drop"
  ) %>%
  mutate(color = colores_pymol[((row_number() - 1) %% length(colores_pymol)) + 1])
 
cat("   Dominios encontrados:", nrow(tabla_agrupada), "\n")
print(tabla_agrupada %>% select(dominio, desc, start_min, end_max, color))
 
# =============================================================================
# 2. CARGAR PDB DE ALPHAFOLD3 Y EXTRAER pLDDT
# =============================================================================
 
cat("\n── [2/5] Cargando PDB de AlphaFold3...\n")
 
pdb <- read.pdb(ARCHIVO_PDB)
 
# En AlphaFold3 el B-factor contiene el pLDDT (0-100)
ca_idx  <- atom.select(pdb, "calpha")
plddt   <- pdb$atom$b[ca_idx$atom]
residuos <- pdb$atom$resno[ca_idx$atom]
 
cat("   Residuos totales:", length(residuos), "\n")
cat("   pLDDT promedio:  ", round(mean(plddt), 2), "\n")
cat("   Residuos < ", UMBRAL_PLDDT, ": ",
    sum(plddt < UMBRAL_PLDDT), " (",
    round(100 * mean(plddt < UMBRAL_PLDDT), 1), "%)\n", sep = "")
 
# =============================================================================
# 3. VALIDACIÓN ESTRUCTURAL
# =============================================================================
 
cat("\n── [3/5] Generando gráficos de validación...\n")
 
pdf(file.path(DIR_SALIDA, "validacion_modelo.pdf"), width = 12, height = 8)
 
# ── Panel 1: pLDDT por residuo ──────────────────────────────────────────────
par(mfrow = c(1, 2))
 
# Color por umbral de confianza
col_plddt <- case_when(
  plddt >= 90 ~ "#1a73e8",   # muy alta — azul
  plddt >= 70 ~ "#34a853",   # alta     — verde
  plddt >= 50 ~ "#fbbc04",   # media    — amarillo
  TRUE        ~ "#ea4335"    # baja     — rojo
)
# Extraer pLDDT por átomo y residuos por token
atom_plddt   <- as.numeric(unlist(data_af3$atom_plddts))
token_resids <- as.numeric(unlist(data_af3$token_res_ids))

# atom_plddts tiene uno por átomo (1845 valores)
# token_res_ids tiene uno por residuo (241 valores)
# Necesitamos emparejarlos — usar token_res_ids como índice de residuos
# y tomar el pLDDT promedio por residuo

atom_chains  <- unlist(data_af3$atom_chain_ids)
token_chains <- unlist(data_af3$token_chain_ids)

# Crear tabla de pLDDT por residuo promediando átomos del mismo residuo
# Primero necesitamos saber a qué residuo pertenece cada átomo
# Usamos los primeros 241 tokens (residuos de la proteína, cadena A)
plddt    <- as.numeric(unlist(data_af3$atom_plddts))
residuos <- seq_along(token_resids)  # índice secuencial 1:241

# Verificar
cat("Residuos:", length(residuos), "\n")
cat("Rango pLDDT:", min(plddt[1:length(residuos)]), "-", 
    max(plddt[1:length(residuos)]), "\n")

# Tomar solo un pLDDT por residuo (el primero de cada grupo = Cα)
# Los átomos están ordenados por residuo en AF3
plddt_por_residuo <- plddt[1:length(token_resids)]
residuos_num      <- token_resids

# Limpiar NAs
validos  <- !is.na(residuos_num) & !is.na(plddt_por_residuo)
residuos <- residuos_num[validos]
plddt    <- plddt_por_residuo[validos]

cat("Residuos válidos:", length(residuos), "\n")
cat("Rango pLDDT:", min(plddt), "-", max(plddt), "\n")
plot(residuos, plddt,
     type  = "l", col = "blue", lwd = 1,
     xlab  = "Número de residuo",
     ylab  = "pLDDT",
     main  = "Confianza por residuo (AlphaFold3)",
     ylim  = c(0, 100))
points(residuos, plddt, pch = 20, cex = 0.6, col = col_plddt)
abline(h = 90, col = "#1a73e8", lty = 2, lwd = 1.5)
abline(h = 70, col = "#34a853", lty = 2, lwd = 1.5)
abline(h = 50, col = "#fbbc04", lty = 2, lwd = 1.5)
legend("bottomleft",
       legend = c("Muy alta (≥90)", "Alta (70-90)", "Media (50-70)", "Baja (<50)"),
       col    = c("#1a73e8", "#34a853", "#fbbc04", "#ea4335"),
       pch    = 20, cex = 0.8, bty = "n")
 
# Sombrear regiones de dominios
for (i in seq_len(nrow(tabla_agrupada))) {
  rect(tabla_agrupada$start_min[i], 0,
       tabla_agrupada$end_max[i],  100,
       col = adjustcolor("purple", alpha.f = 0.08),
       border = NA)
  text(x      = (tabla_agrupada$start_min[i] + tabla_agrupada$end_max[i]) / 2,
       y      = 102,
       labels = tabla_agrupada$tag[i],
       cex    = 0.55, srt = 45, xpd = TRUE)
}
 
# ── Panel 2: Ramachandran ────────────────────────────────────────────────────
tor <- torsion.pdb(pdb)
 
# Clasificar residuos en regiones del Ramachandran
en_alpha <- (tor$phi > -160 & tor$phi < -40) & (tor$psi > -70  & tor$psi < 50)
en_beta  <- (tor$phi > -170 & tor$phi < -50) & (tor$psi > 90   | tor$psi < -170)
 
col_rama <- ifelse(en_alpha, "#1a73e8",
            ifelse(en_beta,  "#34a853", "#ea4335"))
 
plot(tor$phi, tor$psi,
     pch  = 20, cex = 0.7,
     col  = col_rama,
     xlim = c(-180, 180), ylim = c(-180, 180),
     xlab = "Phi (°)", ylab = "Psi (°)",
     main = "Diagrama de Ramachandran")
abline(h = 0, v = 0, col = "pink")
 
# Regiones canónicas aproximadas
rect(-160, -70, -40,  50,  border = "#1a73e8", lty = 2, lwd = 1.5)  # alfa
rect(-170,  90, -50,  180, border = "#34a853", lty = 2, lwd = 1.5)  # beta
 
legend("topright",
       legend = c(sprintf("α-hélice (%d)", sum(en_alpha, na.rm = TRUE)),
                  sprintf("β-hoja   (%d)", sum(en_beta,  na.rm = TRUE)),
                  sprintf("Otro     (%d)", sum(!en_alpha & !en_beta, na.rm = TRUE))),
       col    = c("#1a73e8", "#34a853", "#ea4335"),
       pch    = 20, cex = 0.9, bty = "n")
 
dev.off()
cat("   Guardado: validacion_modelo.pdf\n")
 
# =============================================================================
# 4. MODIFICAR B-FACTORS PARA COLOREAR DOMINIOS EN PYMOL
# =============================================================================
 
cat("\n── [4/5] Codificando dominios en B-factors del PDB...\n")
 
pdb_out <- pdb
 
# Primero, asignar pLDDT como base (dominio = 0)
pdb_out$atom$b <- 0
 
# Asignar índice numérico a cada dominio
for (i in seq_len(nrow(tabla_agrupada))) {
  # Expandir todos los rangos del dominio (puede tener varios segmentos)
  rangos_raw <- tabla_regiones %>%
    filter(dominio == tabla_agrupada$dominio[i])
 
  for (j in seq_len(nrow(rangos_raw))) {
    residuos_dom <- rangos_raw$start[j]:rangos_raw$end[j]
    idx <- which(pdb_out$atom$resno %in% residuos_dom)
    pdb_out$atom$b[idx] <- i
  }
}
 
write.pdb(pdb_out, file.path(DIR_SALIDA, "modelo_dominios_bfactor.pdb"))
cat("   Guardado: modelo_dominios_bfactor.pdb\n")
cat("   → En PyMOL: load modelo_dominios_bfactor.pdb\n")
cat("               spectrum b, rainbow\n")
 
# =============================================================================
# 5. GENERAR SCRIPT .PML PARA PYMOL
# =============================================================================
 
cat("\n── [5/5] Generando script PyMOL (.pml)...\n")
 
nombre_objeto <- tools::file_path_sans_ext(basename(ARCHIVO_PDB))
 
comandos_pymol <- tabla_agrupada %>%
  mutate(
    sel = paste0("select ", tag, ", ", resi_str),
    col = paste0("color ",  color, ", ", tag)
  )
 
# Residuos de baja confianza
resno_bajos <- residuos[plddt < UMBRAL_PLDDT]
sel_baja_confianza <- if (length(resno_bajos) > 0) {
  rangos <- paste0("resi ", resno_bajos, collapse = " + ")
  c(paste0("select baja_confianza, ", rangos),
    "color gray50, baja_confianza",
    "show sticks, baja_confianza")
} else {
  "# No hay residuos con pLDDT bajo el umbral"
}
 
lineas_pml <- c(
  paste0("# Script generado automáticamente — ", Sys.time()),
  paste0("# Proteína: ", nombre_objeto),
  paste0("# Dominios: ", nrow(tabla_agrupada)),
  paste0("# Umbral pLDDT: ", UMBRAL_PLDDT),
  "",
  paste0("load ", ARCHIVO_PDB),
  "hide everything",
  "show cartoon",
  "color gray80",
  "",
  "# ── Dominios de InterPro ──",
  comandos_pymol$sel,
  comandos_pymol$col,
  "",
  paste0("# ── Residuos con pLDDT < ", UMBRAL_PLDDT, " ──"),
  sel_baja_confianza,
  "",
  "# ── Estética final ──",
  "set cartoon_fancy_helices, 1",
  "set cartoon_flat_sheets, 1",
  "set ray_shadows, 0",
  "deselect",
  "zoom"
)
 
writeLines(lineas_pml, file.path(DIR_SALIDA, "mis_regiones.pml"))
cat("   Guardado: mis_regiones.pml\n")
 
# =============================================================================
# RESUMEN FINAL
# =============================================================================
 
cat("\n", strrep("=", 60), "\n")
cat("  ANÁLISIS COMPLETADO\n")
cat(strrep("=", 60), "\n")
cat("  Archivos generados en:", normalizePath(DIR_SALIDA), "\n\n")
cat("  validacion_modelo.pdf         → gráficos de calidad\n")
cat("  modelo_dominios_bfactor.pdb   → PDB con dominios en B-factor\n")
cat("  mis_regiones.pml              → script listo para PyMOL\n\n")
cat("  Resumen de dominios:\n")
tabla_agrupada %>%
  select(dominio, desc, start_min, end_max, color) %>%
  rename(inicio = start_min, fin = end_max) %>%
  { print(as.data.frame(.)); . } %>%
  invisible()
cat(strrep("=", 60), "\n")
# =============================================================================
# 5. GENERAR SCRIPT .PML PARA PYMOL
# =============================================================================

cat("\n── [5/5] Generando script PyMOL (.pml)...\n")

comandos_pymol <- tabla_agrupada %>%
  mutate(
    sel = paste0("select ", tag, ", ", resi_str),
    col = paste0("color ",  color, ", ", tag)
  )

resno_bajos <- residuos[plddt < UMBRAL_PLDDT]
sel_baja_confianza <- if (length(resno_bajos) > 0) {
  rangos <- paste0("resi ", resno_bajos, collapse = " + ")
  c(paste0("select baja_confianza, ", rangos),
    "color gray50, baja_confianza",
    "show sticks, baja_confianza")
} else {
  "# No hay residuos con pLDDT bajo el umbral"
}

lineas_pml <- c(
  paste0("# Script generado — ", Sys.time()),
  "",
  paste0("load ", ARCHIVO_PDB),          # carga el CIF original completo
  "hide everything",
  "show cartoon, polymer.protein",       # ← solo proteína, oculta nanotubo
  "color gray80, polymer.protein",
  "",
  "# ── Dominios de InterPro ──",
  comandos_pymol$sel,
  comandos_pymol$col,
  "",
  paste0("# ── Residuos con pLDDT < ", UMBRAL_PLDDT, " ──"),
  sel_baja_confianza,
  "",
  "# ── Estética final ──",
  "set cartoon_fancy_helices, 1",
  "set cartoon_flat_sheets, 1",
  "set ray_shadows, 0",
  "deselect",
  "zoom polymer.protein"                 # ← zoom solo a la proteína
)

writeLines(lineas_pml, file.path(DIR_SALIDA, "mis_regiones.pml"))
cat("   Guardado: mis_regiones.pml\n")

