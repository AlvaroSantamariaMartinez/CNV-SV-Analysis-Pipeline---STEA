# =============================================================================
# Install_Dependencies.R
# Instalación de todas las dependencias del repositorio CNV/SV Analysis Pipeline
# Scripts:
#   - STEA CNVSV Analysis Pipeline.R
#   - FiltroB_STEA.R
#
# Uso: ejecutar este script una sola vez antes de lanzar la aplicación.
# Requiere R >= 4.3 y conexión a internet.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. PAQUETES DE CRAN
# -----------------------------------------------------------------------------
cran_packages <- c(
  # Interfaz Shiny y UI
  "shiny",        # Framework principal de la aplicación web
  "bslib",        # Temas Bootstrap modernos y layout cards/navbars
  "bsicons",      # Iconos Bootstrap para bslib
  "DT",           # Tablas interactivas (DataTables)
  "plotly",       # Gráficos interactivos
  "shinyFiles",   # Selector de archivos/carpetas nativo en Shiny
  "shinyjs",      # JavaScript helpers para Shiny (enable/disable, runjs…)

  # Ejecución de procesos externos
  "processx",     # Lanzar y monitorear procesos R externos (pipeline backend)

  # Lectura y escritura de archivos
  "openxlsx",     # Crear y escribir archivos Excel (.xlsx) con estilos
  "readxl",       # Leer archivos Excel (.xlsx / .xls)
  "data.table",   # Lectura rápida de TSVs (fread) y manipulación eficiente

  # Manipulación de datos
  "dplyr",        # Verbos de transformación de datos (filter, mutate, bind_rows…)
  "tidyr",        # Reshaping de datos (pivot, unnest…)
  "stringr",      # Manipulación de cadenas de texto (str_trim, str_detect…)

  # Visualización
  "ggplot2",      # Gráficos estáticos base
  "gridExtra",    # Composición de múltiples gráficos ggplot2
  "scales",       # Escalas y formateo de ejes en ggplot2

  # Comunicación HTTP y datos JSON
  "httr",         # Peticiones HTTP (consultas a APIs externas: ClinVar, OMIM…)
  "jsonlite",     # Parseo y serialización de JSON

  # IDE / entorno de desarrollo
  "rstudioapi"    # Detectar ruta del script activo en RStudio
)

cat("=============================================================================\n")
cat("Instalando paquetes CRAN...\n")
cat("=============================================================================\n")

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) == 0) {
    cat("✓ Todos los paquetes CRAN ya están instalados.\n")
    return(invisible(NULL))
  }
  cat(paste0("  Instalando: ", paste(missing_pkgs, collapse = ", "), "\n"))
  install.packages(
    missing_pkgs,
    dependencies = TRUE,
    repos        = "https://cloud.r-project.org"
  )
}

install_if_missing(cran_packages)


# -----------------------------------------------------------------------------
# 2. PAQUETES DE BIOCONDUCTOR
# -----------------------------------------------------------------------------
bioc_packages <- c(
  "GenomicRanges",  # Representación y operaciones sobre rangos genómicos
                    # (findOverlaps, pintersect, IRanges, GRanges…)
  "IRanges",        # Dependencia de GenomicRanges (rangos de enteros)
  "S4Vectors",      # Infraestructura S4 compartida por Bioconductor
  "BiocGenerics"    # Generics compartidos por todos los paquetes Bioconductor
)

cat("\n=============================================================================\n")
cat("Instalando paquetes Bioconductor...\n")
cat("=============================================================================\n")

# Instalar BiocManager si no está disponible
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("  Instalando BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[!bioc_packages %in% rownames(installed.packages())]
if (length(missing_bioc) == 0) {
  cat("✓ Todos los paquetes Bioconductor ya están instalados.\n")
} else {
  cat(paste0("  Instalando: ", paste(missing_bioc, collapse = ", "), "\n"))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}


# -----------------------------------------------------------------------------
# 3. VERIFICACIÓN FINAL
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("Verificación de dependencias...\n")
cat("=============================================================================\n")

all_packages <- c(cran_packages, bioc_packages)
installed_ok  <- all_packages[all_packages %in% rownames(installed.packages())]
not_installed <- setdiff(all_packages, installed_ok)

cat(paste0("✓ Instalados correctamente: ", length(installed_ok), "/", length(all_packages), "\n"))

if (length(not_installed) > 0) {
  cat("\n⚠  Los siguientes paquetes NO se pudieron instalar:\n")
  cat(paste0("   - ", not_installed, collapse = "\n"), "\n")
  cat("\nIntenta instalarlos manualmente:\n")
  cat('  install.packages(c("', paste(intersect(not_installed, cran_packages), collapse = '", "'), '"))\n')
  if (length(intersect(not_installed, bioc_packages)) > 0) {
    cat('  BiocManager::install(c("', paste(intersect(not_installed, bioc_packages), collapse = '", "'), '"))\n')
  }
} else {
  cat("\n✓ Todas las dependencias están disponibles.\n")
  cat("  El pipeline está listo para ejecutarse.\n")
}

cat("=============================================================================\n")
