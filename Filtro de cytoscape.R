library(tidyverse)

# Cargar el archivo TSV
#df <- read_tsv("AtHRE1.tsv"),("AtHRE2.tsv")
archivos <- list.files(pattern = "\\.tsv$", full.names = TRUE)
df <- do.call(rbind, lapply(archivos, read_tsv))

# Renombrar columnas si tienen caracteres especiales (el # en la primera)
colnames(df)[1] <- "node"

# --- OPCIÓN 1: Buscar identifiers con una palabra específica ---
palabra_buscar <- "defense","pathogen","immune"  # <-- cambia esta palabra

df_filtrado <- df %>%
  filter(str_detect(`term description`, regex(palabra_buscar, ignore_case = TRUE)))

# Ver los identifiers únicos que tienen esa palabra
unique(df_filtrado$identifier)


# --- OPCIÓN 2: Ver qué palabras clave aparecen por identifier ---
# (útil para explorar qué términos tiene cada identifier)

resumen_por_identifier <- df %>%
  group_by(identifier) %>%
  summarise(
    n_terminos = n(),
    terminos_unicos = n_distinct(`term description`),
    todas_las_descripciones = paste(unique(`term description`), collapse = " | ")
  ) %>%
  arrange(desc(terminos_unicos))


# --- OPCIÓN 3: Buscar múltiples palabras y ver qué identifiers las tienen ---
palabras <- c("Pathogen", "immune", "resistance")  

resultado <- map_dfr(palabras, function(p) {
  df %>%
    filter(str_detect(`term description`, regex(p, ignore_case = TRUE))) %>%
    distinct(identifier) %>%
    mutate(palabra_encontrada = p)
})

# Ver tabla: qué identifier tiene qué palabras
resultado_wide <- resultado %>%
  mutate(presente = TRUE) %>%
  pivot_wider(names_from = palabra_encontrada, values_from = presente, values_fill = FALSE)



# --- OPCIÓN 4: Identifiers que tienen TODAS las palabras del vector ---
identifiers_con_todas <- resultado %>%
  group_by(identifier) %>%
  summarise(n_palabras = n_distinct(palabra_encontrada)) %>%
  filter(n_palabras == length(palabras)) %>%
  pull(identifier)

install.packages("writexl") # Solo una vez
library(writexl)
write_xlsx(resultado_wide, "filtros_cytoscape_final.xlsx")

