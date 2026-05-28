library(jsonlite)
library(tidyverse)

data <- fromJSON("iprscan5-R20260323-214108-0128-69166529-p2m.json", simplifyVector = FALSE)

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

# ── CORRECCIONES ────────────────────────────────────────────────────────────────

colores <- c("salmon", "brown", "orange", "purple", "magenta", "red",
             "yellow", "cyan", "green", "blue", "wheat", "slate")

# 1. Agrupar rangos del MISMO dominio con "+" para una sola selección
tabla_agrupada <- tabla_regiones %>%
  mutate(tag = gsub("[^[:alnum:]]", "_", dominio)) %>%
  group_by(dominio, tag, desc) %>%
  summarise(
    resi_str = paste0("resi ", start, "-", end, collapse = " + "),
    .groups  = "drop"
  ) %>%
  # 2. Asignar color por índice de fila (no con %%, que empieza en 0 y puede repetir)
  mutate(color = colores[((row_number() - 1) %% length(colores)) + 1])

# 3. Generar comandos PyMOL
comandos <- tabla_agrupada %>%
  mutate(
    sel = paste0("select ", tag, ", ", resi_str),
    col = paste0("color ",  color, ", ", tag)
  )

# 4. Escribir el .pml con un "show sticks" opcional al final
lineas <- c(
  "# Regiones de dominio generadas automáticamente",
  comandos$sel,
  comandos$col,
  "deselect"          # limpia la selección activa al final
)

writeLines(lineas, "mis_regiones.pml")
cat("Archivo 'mis_regiones.pml' generado con éxito.\n")
