# =============================================================================
# app.R — Shiny interface for the Genomic CNV/SV Analysis Pipeline
# =============================================================================

library(shiny)
library(bslib)
library(DT)
library(plotly)
library(shinyFiles)
library(shinyjs)
library(processx)
library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(httr)
library(jsonlite)
library(bsicons)

# ── Operador null-coalescing (debe definirse ANTES de cualquier uso) ──────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Ruta al script del pipeline (mismo directorio que app.R) ─────────────────
PIPELINE_SCRIPT <- tryCatch({
  doc_path <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(doc_path))
    file.path(dirname(doc_path), "FiltroB_STEA.R")
  else
    file.path(getwd(), "FiltroB_STEA.R")
}, error = function(e) {
  file.path(getwd(), "FiltroB_STEA.R")
})
if (!file.exists(PIPELINE_SCRIPT))
  PIPELINE_SCRIPT <- file.path(getwd(), "FiltroB_STEA.R")

# =============================================================================
# DEFAULT CONFIGURATION VALUES
# =============================================================================
CONFIG_DEFAULTS <- list(
  ruta_entrada        = "",
  ruta_salida         = file.path(getwd(), "results"),
  ruta_rangos         = file.path(getwd(), "reference", "reference_ranges.xlsx"),
  ruta_panel          = file.path(getwd(), "reference", "gene_panel.txt"),
  umbral_strict       = 0.70,
  umbral_herencia     = 0.70,
  umbral_herencia_lax = 0.60,
  umbral_jaccard      = 0.50,
  margen_lateral      = 3e6,
  deduplicar          = TRUE,
  n_cores             = max(1L, parallel::detectCores(logical = FALSE) - 1L)
)

#=============================================================================
# PALETA DE COLORES HERENCIA / RANGO
# =============================================================================
HER_COLORS <- c(
  "De novo"           = "#CC0000",
  "Paternal"           = "#0066CC",
  "Maternal"           = "#CC6600",
  "Combined"          = "#6600CC",
  "Paternal (Probable)"= "#4499DD",
  "Maternal (Probable)"= "#DD9944",
  "Combined (Probable)"="#9944DD",
  "Desconocido"       = "#D5D8DC",
  "Sin datos"         = "#AAAAAA"
)
RANGO_COLORS <- c(
  "Strict" = "#155724",
  "Wide"   = "#856404",
  "Outside"    = "#9C0006"
)

# Longitudes cromosomicas hg38 (pb)
CHR_LENGTHS <- c(
  "1"=248956422, "2"=242193529, "3"=198295559, "4"=190214555,
  "5"=181538259, "6"=170805979, "7"=159345973, "8"=145138636,
  "9"=138394717, "10"=133797422,"11"=135086622,"12"=133275309,
  "13"=114364328,"14"=107043718,"15"=101991189,"16"=90338345,
  "17"=83257441, "18"=80373285, "19"=58617616, "20"=64444167,
  "21"=46709983, "22"=50818468, "X"=156040895, "Y"=57227415
)
CHR_ORDER <- c(as.character(1:22), "X", "Y")

# =============================================================================
# HELPERS
# =============================================================================

abrir_archivo_nativo <- function(ruta) {
  if (!file.exists(ruta)) return(FALSE)
  os <- .Platform$OS.type
  if (os == "windows") {
    shell.exec(normalizePath(ruta))
  } else if (Sys.info()["sysname"] == "Darwin") {
    system2("open", shQuote(ruta))
  } else {
    system2("xdg-open", shQuote(ruta))
  }
  TRUE
}

cargar_todos_resultados <- function(ruta_salida) {
  message("\n--- CARGANDO RESULTADOS ---")
  ruta_norm <- normalizePath(trimws(ruta_salida), winslash = "/", mustWork = FALSE)
  
  archivos <- list.files(ruta_norm,
                         pattern    = "_Analisis_Completo\\.xlsx$",
                         recursive  = TRUE,
                         full.names = TRUE)
  
  if (length(archivos) == 0) {
    g1 <- Sys.glob(file.path(ruta_norm, "*_Analisis_Completo.xlsx"))
    g2 <- Sys.glob(file.path(ruta_norm, "*",  "*_Analisis_Completo.xlsx"))
    g3 <- Sys.glob(file.path(ruta_norm, "*", "*", "*_Analisis_Completo.xlsx"))
    archivos <- unique(c(g1, g2, g3))
  }
  
  if (length(archivos) == 0) {
    message("⚠️  No Complete_Analysis.xlsx files found")
    return(list(CNVs = NULL, SVs = NULL))
  }
  
  leer_hoja <- function(archivo, hoja) {
    tryCatch({
      df <- readxl::read_excel(archivo, sheet = hoja) |> as.data.frame()
      if (nrow(df) > 0) df[] <- lapply(df, as.character)
      df$._archivo_ <- archivo
      df
    }, error = function(e) {
      NULL
    })
  }
  
  cnvs_list <- Filter(Negate(is.null), lapply(archivos, leer_hoja, hoja = "CNVs"))
  svs_list  <- Filter(Negate(is.null), lapply(archivos, leer_hoja, hoja = "SVs"))
  
  cnvs <- if (length(cnvs_list) > 0) bind_rows(cnvs_list) else data.frame()
  svs  <- if (length(svs_list)  > 0) bind_rows(svs_list)  else data.frame()
  
  list(
    CNVs = if (nrow(cnvs) > 0) cnvs else NULL,
    SVs  = if (nrow(svs)  > 0) svs  else NULL
  )
}

her_base <- function(x) {
  x <- sub("\\s*\\[.*", "", as.character(x))
  x <- trimws(sub("\\s*\\(Probable\\)", "", x))
  x
}

# Detecta si un vector de valores de En_Panel_Genes es TRUE/positivo,
# cubriendo todas las variantes (booleano, "TRUE", "Yes", "Yes", "Yes", "Yes",
# "1", "yes") de forma robusta e insensible a tildes y mayusculas.
es_sfari_col <- function(x) {
  x_norm <- toupper(trimws(iconv(as.character(x), to = "ASCII//TRANSLIT")))
  x_norm %in% c("TRUE", "Yes", "1", "YES") | (x == TRUE)
}

# =============================================================================
# SISTEMA DE TagS — HELPERS
# =============================================================================
COLS_CLAVE <- c("ID_Familia", "SV_chrom", "SV_start", "SV_end", "SV_type")

hacer_clave_variante <- function(df) {
  n <- if (is.null(df)) 0L else nrow(df)
  if (n == 0L) return(character(0))
  if (!all(COLS_CLAVE %in% names(df)))
    return(rep(NA_character_, n))
  # Incluir Annotation_mode para que filas split y full sean claves independientes
  mode_col <- if ("Annotation_mode" %in% names(df))
    toupper(trimws(as.character(df$Annotation_mode)))
  else
    rep("FULL", n)
  # Normalizar start/end a entero para evitar notacion cientifica entre recargas
  start_norm <- as.character(suppressWarnings(as.integer(as.numeric(df$SV_start))))
  end_norm   <- as.character(suppressWarnings(as.integer(as.numeric(df$SV_end))))
  paste(
    as.character(df$ID_Familia),
    toupper(trimws(gsub("CHR", "", as.character(df$SV_chrom)))),
    start_norm,
    end_norm,
    toupper(as.character(df$SV_type)),
    mode_col,
    sep = "|"
  )
}

detectar_similares_interindividuales <- function(df, umbral_sim = 0.70, verbose = FALSE) {
  if (is.null(df) || nrow(df) == 0) {
    if (verbose) return(list(claves = character(0), n_pares = 0L, n_variantes = 0L))
    return(character(0))
  }
  if (!all(COLS_CLAVE %in% names(df))) {
    if (verbose) return(list(claves = character(0), n_pares = 0L, n_variantes = 0L))
    return(character(0))
  }
  
  if ("Annotation_mode" %in% names(df))
    df <- df[!is.na(df$Annotation_mode) & df$Annotation_mode == "full", , drop = FALSE]
  
  if (nrow(df) <= 1) {
    if (verbose) return(list(claves = character(0), n_pares = 0L, n_variantes = 0L))
    return(character(0))
  }
  
  keys      <- hacer_clave_variante(df)
  chr_vec   <- toupper(trimws(gsub("CHR", "", as.character(df$SV_chrom))))
  sv_up     <- toupper(as.character(df$SV_type))
  tipo_vec  <- ifelse(grepl("DEL|LOSS", sv_up), "DEL",
                      ifelse(grepl("DUP|GAIN", sv_up), "DUP", sv_up))
  start_vec <- suppressWarnings(as.integer(df$SV_start))
  end_vec   <- suppressWarnings(as.integer(df$SV_end))
  len_vec   <- pmax(1L, end_vec - start_vec + 1L, na.rm = TRUE)
  fam_vec   <- as.character(df$ID_Familia)
  
  # Contador de similares por variante (solo activa auto-tag si > 2 similares)
  sim_count <- list()
  n_pares   <- 0L
  grupos    <- unique(paste(chr_vec, tipo_vec, sep = "::"))
  
  for (g in grupos) {
    tipo_g <- strsplit(g, "::", fixed = TRUE)[[1]][2]
    if (tipo_g == "OTHER") next
    idx <- which(paste(chr_vec, tipo_vec, sep = "::") == g)
    if (length(idx) < 2) next
    n_g <- length(idx)
    ij  <- which(upper.tri(matrix(TRUE, n_g, n_g)), arr.ind = TRUE)
    ii  <- idx[ij[, 1]]
    jj  <- idx[ij[, 2]]
    cross <- fam_vec[ii] != fam_vec[jj]
    ii <- ii[cross]; jj <- jj[cross]
    if (length(ii) == 0) next
    ovl <- pmax(0L, pmin(end_vec[ii], end_vec[jj]) - pmax(start_vec[ii], start_vec[jj]) + 1L)
    sim <- ovl / pmax(len_vec[ii], len_vec[jj])
    ok <- !is.na(sim) & sim >= umbral_sim
    if (any(ok)) {
      n_pares <- n_pares + sum(ok)
      for (p in which(ok)) {
        ki <- keys[ii[p]]; kj <- keys[jj[p]]
        sim_count[[ki]] <- union(sim_count[[ki]], kj)
        sim_count[[kj]] <- union(sim_count[[kj]], ki)
      }
    }
  }
  
  # Solo auto-etiqueta variantes con MAS de 2 similares interindividuales
  flagged <- names(Filter(function(v) length(v) > 2, sim_count))
  
  if (verbose)
    list(claves = flagged, n_pares = n_pares, n_variantes = length(flagged))
  else
    flagged
}

# Detecta variantes del dataset completo similares a una clave dada (para tageo manual)
detectar_similares_a_clave <- function(clave, df_full, umbral_sim = 0.70) {
  if (is.null(df_full) || nrow(df_full) == 0) return(character(0))
  if (!all(COLS_CLAVE %in% names(df_full))) return(character(0))
  
  partes <- strsplit(clave, "|", fixed = TRUE)[[1]]
  if (length(partes) < 5L) return(character(0))
  
  fam_ref   <- partes[1]
  chr_ref   <- partes[2]
  start_ref <- suppressWarnings(as.integer(partes[3]))
  end_ref   <- suppressWarnings(as.integer(partes[4]))
  tipo_ref  <- partes[5]
  if (is.na(start_ref) || is.na(end_ref)) return(character(0))
  len_ref   <- max(1L, end_ref - start_ref + 1L)
  
  if ("Annotation_mode" %in% names(df_full))
    df_full <- df_full[!is.na(df_full$Annotation_mode) & df_full$Annotation_mode == "full", , drop = FALSE]
  if (nrow(df_full) == 0) return(character(0))
  
  keys_full  <- hacer_clave_variante(df_full)
  chr_vec    <- toupper(trimws(gsub("CHR", "", as.character(df_full$SV_chrom))))
  sv_up      <- toupper(as.character(df_full$SV_type))
  tipo_vec   <- ifelse(grepl("DEL|LOSS", sv_up), "DEL",
                       ifelse(grepl("DUP|GAIN", sv_up), "DUP", sv_up))
  tipo_ref_n <- ifelse(grepl("DEL|LOSS", tipo_ref), "DEL",
                       ifelse(grepl("DUP|GAIN", tipo_ref), "DUP", tipo_ref))
  start_vec  <- suppressWarnings(as.integer(df_full$SV_start))
  end_vec    <- suppressWarnings(as.integer(df_full$SV_end))
  len_vec    <- pmax(1L, end_vec - start_vec + 1L, na.rm = TRUE)
  fam_vec    <- as.character(df_full$ID_Familia)
  
  # Solo candidatos: mismo cromosoma, mismo tipo, distinta familia
  candidatos <- which(chr_vec == chr_ref & tipo_vec == tipo_ref_n & fam_vec != fam_ref)
  if (length(candidatos) == 0) return(character(0))
  
  ovl <- pmax(0L, pmin(end_ref, end_vec[candidatos]) - pmax(start_ref, start_vec[candidatos]) + 1L)
  sim <- ovl / pmax(len_ref, len_vec[candidatos])
  
  ok <- !is.na(sim) & sim >= umbral_sim
  keys_full[candidatos[ok]]
}

encontrar_archivo_hijo <- function(ruta_entrada, archivo_excel, id_familia, modalidad) {
  nombre_sub <- basename(dirname(normalizePath(archivo_excel, mustWork = FALSE)))
  dir_hi <- file.path(ruta_entrada, nombre_sub, paste0(id_familia, "HI"))
  
  if (!dir.exists(dir_hi)) {
    candidatos <- list.dirs(file.path(ruta_entrada, nombre_sub), full.names = TRUE, recursive = FALSE)
    dir_hi <- candidatos[grepl(paste0(id_familia, "HI$"), candidatos, ignore.case = TRUE)][1]
  }
  if (is.na(dir_hi) || !dir.exists(dir_hi)) return(NA_character_)
  
  tsvs <- list.files(dir_hi, pattern = "annotated\\.tsv$", full.names = TRUE, recursive = TRUE)
  if (length(tsvs) == 0) return(NA_character_)
  
  patron <- if (toupper(modalidad) == "CNV") "HI[.]CNVs" else "HI[.]SVs"
  match  <- tsvs[grepl(patron, basename(tsvs), ignore.case = TRUE)]
  if (length(match) > 0) return(match[1])
  tsvs[1]
}

# =============================================================================
# UI
# =============================================================================
ui <- page_navbar(
  id = "main_navbar",
  title = tags$span(
    tags$img(src = "https://www.r-project.org/logo/Rlogo.svg",
             height = "26px", style = "margin-right:8px; vertical-align:middle;"),
    "CNV/SV Analysis Pipeline"
  ),
  theme = bs_theme(
    bootswatch   = "flatly",
    primary      = "#2C6FAC",
    base_font    = font_google("Inter"),
    heading_font = font_google("Inter"),
    font_scale   = 0.92
  ),
  header = tagList(
    useShinyjs(),
    # ── Toggle modo oscuro ──────────────────────────────────────────────────
    div(class = "d-flex align-items-center me-2",
        bslib::input_dark_mode(id = "dm_toggle", mode = "light")
    ),
    # ── Semantic CSS adaptable to dark-mode ─────────────────────────────────
    tags$style(HTML("
      /* ── Variables de tema ───────────────────────────────────────────── */
      :root {
        --dm-panel-bg:          #f8f9fa;
        --dm-panel-border:      #dee2e6;
        --dm-info-bg:           #f0f7ff;
        --dm-info-border:       #cce0ff;
        --dm-info-left:         #f0f7ff;
        --dm-success-bg:        #f0fff4;
        --dm-success-left:      #27AE60;
        --dm-warn-bg:           #fff8f0;
        --dm-warn-border:       #ffe082;
        --dm-warn-left:         #E67E22;
        --dm-purple-bg:         #f8f0ff;
        --dm-purple-left:       #8E44AD;
        --dm-muted-left:        #6c757d;
        --dm-fenotipo-bg:       #F5EEF8;
        --dm-fenotipo-border:   #D2B4DE;
        --dm-fenotipo-input-b:  #C39BD3;
        --dm-preview-ft-bg:     #EDE0F5;
        --dm-preview-ft-color:  #4A235A;
        --dm-preview-ft-border: #C39BD3;
        --dm-preview-ft-e-bg:   #F5EEF8;
        --dm-preview-ft-e-col:  #9B59B6;
        --dm-preview-nt-bg:     #D6EAF8;
        --dm-preview-nt-color:  #1A5276;
        --dm-preview-nt-border: #85C1E9;
        --dm-preview-nt-e-bg:   #EBF5FB;
        --dm-preview-nt-e-col:  #2E86C1;
        --dm-tag-manual-bg:     #d4edda;
        --dm-tag-manual-left:   #28a745;
        --dm-tag-auto-bg:       #e8f5e9;
        --dm-tag-auto-left:     #66bb6a;
        --dm-loc-badge-bg:      #cce9f5;
        --dm-loc-badge-color:   #000000;
        --dm-loc-badge-border:  #99d3ec;
        --dm-stat-blue-bg:      #EBF5FB;
        --dm-stat-gray-bg:      #F9F9F9;
        --dm-muted-text:        #555;
        --dm-dark-text:         #444;
        --dm-h0-bg:             #e8f4fd;
        --dm-uso-bg:            #f0fff0;
      }
      [data-bs-theme='dark'] {
        --dm-panel-bg:          #1e2530;
        --dm-panel-border:      #3a4555;
        --dm-info-bg:           #1a2840;
        --dm-info-border:       #2a4060;
        --dm-info-left:         #1a2840;
        --dm-success-bg:        #1a2e22;
        --dm-success-left:      #2ecc71;
        --dm-warn-bg:           #2e2010;
        --dm-warn-border:       #7a5800;
        --dm-warn-left:         #f39c12;
        --dm-purple-bg:         #251535;
        --dm-purple-left:       #a569bd;
        --dm-muted-left:        #6c757d;
        --dm-fenotipo-bg:       #2a1f38;
        --dm-fenotipo-border:   #6c3483;
        --dm-fenotipo-input-b:  #8e44ad;
        --dm-preview-ft-bg:     #3d2260;
        --dm-preview-ft-color:  #d7b4f5;
        --dm-preview-ft-border: #8e44ad;
        --dm-preview-ft-e-bg:   #2a1f38;
        --dm-preview-ft-e-col:  #a569bd;
        --dm-preview-nt-bg:     #1a3050;
        --dm-preview-nt-color:  #aed6f1;
        --dm-preview-nt-border: #2e86c1;
        --dm-preview-nt-e-bg:   #1a2535;
        --dm-preview-nt-e-col:  #5dade2;
        --dm-tag-manual-bg:     #1a3325;
        --dm-tag-manual-left:   #27ae60;
        --dm-tag-auto-bg:       #1a2e20;
        --dm-tag-auto-left:     #52be80;
        --dm-loc-badge-bg:      #1a3050;
        --dm-loc-badge-color:   #aed6f1;
        --dm-loc-badge-border:  #2e86c1;
        --dm-stat-blue-bg:      #1a2840;
        --dm-stat-gray-bg:      #1e2530;
        --dm-muted-text:        #adb5bd;
        --dm-dark-text:         #ced4da;
        --dm-h0-bg:             #1a2840;
        --dm-uso-bg:            #1a2e22;
      }

      /* ── Semantic classes ────────────────────────────────────────────── */
      .dm-panel           { background:var(--dm-panel-bg)   !important; border-color:var(--dm-panel-border)   !important; }
      .dm-info-panel      { background:var(--dm-info-bg)    !important; border-color:var(--dm-info-border)    !important; }
      .dm-info-left       { background:var(--dm-info-left)  !important; border-left-color:#2C6FAC             !important; }
      .dm-success-left    { background:var(--dm-success-bg) !important; border-left-color:var(--dm-success-left) !important; }
      .dm-warn-left       { background:var(--dm-warn-bg)    !important; border-left-color:var(--dm-warn-left)  !important; }
      .dm-purple-left     { background:var(--dm-purple-bg)  !important; border-left-color:var(--dm-purple-left) !important; }
      .dm-muted-left      { background:var(--dm-panel-bg)   !important; border-left-color:var(--dm-muted-left) !important; }
      .dm-fenotipo-panel  { background:var(--dm-fenotipo-bg)!important; border-color:var(--dm-fenotipo-border) !important; }
      .dm-tag-manual      { background:var(--dm-tag-manual-bg)!important; border-left-color:var(--dm-tag-manual-left)!important; }
      .dm-tag-auto        { background:var(--dm-tag-auto-bg) !important; border-left-color:var(--dm-tag-auto-left) !important; }
      .dm-loc-badge       { background:var(--dm-loc-badge-bg)!important; color:var(--dm-loc-badge-color)!important; border-color:var(--dm-loc-badge-border)!important; }
      .dm-stat-blue       { background:var(--dm-stat-blue-bg)!important; }
      .dm-stat-gray       { background:var(--dm-stat-gray-bg)!important; }
      .dm-h0-panel        { background:var(--dm-h0-bg)      !important; border-left-color:#2C6FAC !important; }
      .dm-uso-panel       { background:var(--dm-uso-bg)     !important; border-left-color:var(--dm-success-left) !important; }
      .dm-warn-panel      { background:var(--dm-warn-bg)    !important; border-color:var(--dm-warn-border) !important; }
      .dm-muted-txt       { color:var(--dm-muted-text)      !important; }
      .dm-dark-txt        { color:var(--dm-dark-text)       !important; }

      /* Previews fenotipo */
      .dm-preview-ft      { background:var(--dm-preview-ft-bg)   !important;
                            color:var(--dm-preview-ft-color)     !important;
                            border-color:var(--dm-preview-ft-border) !important; }
      .dm-preview-ft-e    { background:var(--dm-preview-ft-e-bg) !important;
                            color:var(--dm-preview-ft-e-col)     !important;
                            border-color:var(--dm-fenotipo-input-b) !important; }
      /* Previews notas */
      .dm-preview-nt      { background:var(--dm-preview-nt-bg)   !important;
                            color:var(--dm-preview-nt-color)     !important;
                            border-color:var(--dm-preview-nt-border) !important; }
      .dm-preview-nt-e    { background:var(--dm-preview-nt-e-bg) !important;
                            color:var(--dm-preview-nt-e-col)     !important;
                            border-color:var(--dm-preview-nt-border) !important; }

      /* Input fenotipo */
      .dm-fenotipo-input  { border-color:var(--dm-fenotipo-input-b) !important; }

      /* Test result panels */
      .dm-result-panel    { background:var(--dm-panel-bg) !important;
                            border-color:var(--dm-panel-border) !important; }

      /* About page cards */
      .dm-about-guide     { background:var(--dm-info-left)  !important;
                            border-left-color:#2C6FAC !important; }
      .dm-about-success   { background:var(--dm-success-bg) !important;
                            border-left-color:var(--dm-success-left) !important; }
      .dm-about-warn      { background:var(--dm-warn-bg)    !important;
                            border-left-color:var(--dm-warn-left) !important; }
      .dm-about-purple    { background:var(--dm-purple-bg)  !important;
                            border-left-color:var(--dm-purple-left) !important; }
      .dm-about-muted     { background:var(--dm-panel-bg)   !important;
                            border-left-color:var(--dm-muted-left) !important; }
      .dm-about-code      { background:var(--dm-panel-bg)   !important;
                            border-color:var(--dm-panel-border) !important; }

      /* Dark mode: plot backgrounds */
      [data-bs-theme='dark'] .plot-container { background:transparent !important; }

      /* Dark mode: DTtable fixes */
      [data-bs-theme='dark'] table.dataTable thead th,
      [data-bs-theme='dark'] table.dataTable tbody td { color: var(--bs-body-color); }
      [data-bs-theme='dark'] .dataTables_wrapper { color: var(--bs-body-color); }

      /* Wizard: inactive cards in dark */
      [data-bs-theme='dark'] .wiz-card-inactive {
        border-color: #3a4555 !important;
        background: #1e2530 !important;
      }
      [data-bs-theme='dark'] .wiz-card-inactive .fw-bold { color: #ced4da !important; }
      [data-bs-theme='dark'] .wiz-card-inactive .text-muted { color: #6c757d !important; }

      /* Clasif panel */
      [data-bs-theme='dark'] .dm-clasif-panel {
        background: #1e2530 !important;
        border-color: #3a4555 !important;
      }
      [data-bs-theme='dark'] .dm-clasif-panel .dm-dark-txt { color: #ced4da !important; }
    "))
  ),
  
  nav_panel(
    title = tagList(bsicons::bs_icon("gear-fill"), " Configuration"),
    value = "tab_config",
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header(class = "bg-primary text-white", bsicons::bs_icon("folder2-open"), " Rutas de archivos"),
        card_body(
          textInput("cfg_entrada",  "📂 Input folder (families HI/PA/MA)", value = CONFIG_DEFAULTS$ruta_entrada, width = "100%"),
          textInput("cfg_salida",   "📁 Output folder (results)", value = CONFIG_DEFAULTS$ruta_salida,  width = "100%"),
          textInput("cfg_rangos",   "📋 Reference ranges file (.xlsx)", value = CONFIG_DEFAULTS$ruta_rangos,  width = "100%"),
          textInput("cfg_panel",    "🧬 Panel de genes (.txt)", value = CONFIG_DEFAULTS$ruta_panel,   width = "100%"),
          actionButton("btn_validar_rutas", "✔ Validar rutas", class = "btn-outline-primary btn-sm mt-1")
        )
      ),
      card(
        card_header(class = "bg-primary text-white", bsicons::bs_icon("sliders"), " Thresholds and parameters"),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            numericInput("cfg_umbral_strict", "Umbral estricto", value = CONFIG_DEFAULTS$umbral_strict, min = 0, max = 1, step = 0.05),
            numericInput("cfg_umbral_her", "Sim herencia (Fuerte)", value = CONFIG_DEFAULTS$umbral_herencia, min = 0, max = 1, step = 0.05)
          ),
          layout_columns(
            col_widths = c(6, 6),
            numericInput("cfg_umbral_her_lax", "Sim herencia (Probable)", value = CONFIG_DEFAULTS$umbral_herencia_lax, min = 0, max = 1, step = 0.05),
            numericInput("cfg_umbral_jaccard", "Umbral Jaccard", value = CONFIG_DEFAULTS$umbral_jaccard, min = 0, max = 1, step = 0.05)
          ),
          layout_columns(
            col_widths = c(6, 6),
            numericInput("cfg_margen", "Margen lateral (Mb)", value = CONFIG_DEFAULTS$margen_lateral / 1e6, min = 0, max = 20, step = 0.5),
            numericInput("cfg_cores", "Processing cores", value = CONFIG_DEFAULTS$n_cores, min = 1, max = parallel::detectCores(), step = 1)
          ),
          checkboxInput("cfg_dedup", "Deduplicar regiones similares", value = CONFIG_DEFAULTS$deduplicar),
          hr(),
          div(class = "d-flex align-items-center gap-2",
              actionButton("btn_guardar_config", "💾 Save configuration", class = "btn-success"),
              actionButton("btn_reset_config",   "↺ Restore defaults", class = "btn-outline-secondary btn-sm"),
              
              # --- NEW DARK MODE BUTTON ---
              tags$span(class = "ms-auto"), # Empuja el interruptor a la derecha
              tags$span(class = "text-muted small fw-bold mb-0", "Tema visual:"),
              input_dark_mode(id = "modo_noche", mode = "light")
              # ----------------------------------
          ),
          br(),
          uiOutput("ui_validacion_rutas")
        )
      )
    )
  ),
  
  nav_panel(
    title = tagList(bsicons::bs_icon("play-circle-fill"), " Execution"),
    value = "tab_run",
    layout_columns(
      col_widths = c(4, 8),
      card(
        card_header(class = "bg-success text-white", bsicons::bs_icon("rocket-takeoff"), " Launch analysis"),
        card_body(
          p(class = "text-muted small", "Ejecuta el pipeline completo sobre todas las familias detectadas en la carpeta de entrada."),
          hr(),
          div(
            class = "d-grid gap-2",
            actionButton("btn_run_todo",    "▶ Ejecutar cohorte completa", class = "btn-success btn-lg"),
            actionButton("btn_run_carpeta", "▶ Ejecutar carpeta concreta", class = "btn-outline-success"),
            textInput("inp_carpeta_ruta", NULL, placeholder = "Ruta de carpeta batch", width = "100%"),
            actionButton("btn_run_familia", "▶ Ejecutar familia concreta", class = "btn-outline-success"),
            textInput("inp_familia_id", NULL, placeholder = "ID familia (p.ej. G04PMP069)", width = "100%")
          ),
          hr(),
          actionButton("btn_stop", "⏹ Detener", class = "btn-danger btn-sm", disabled = TRUE),
          hr(),
          div(class = "d-flex align-items-center gap-2", uiOutput("ui_run_status"))
        )
      ),
      card(
        card_header(class = "bg-dark text-white", bsicons::bs_icon("terminal"), " Log en tiempo real"),
        card_body(
          class = "p-0", style = "background:#1e1e1e;",
          tags$div(
            id    = "log_container",
            style = paste("height:480px; overflow-y:auto; padding:14px;", "font-family:'Fira Mono','Courier New',monospace;", "font-size:12.5px; color:#d4d4d4; white-space:pre-wrap;"),
            textOutput("log_output", inline = FALSE)
          )
        )
      )
    )
  ),
  
  nav_panel(
    title = tagList(bsicons::bs_icon("table"), " Resultados"),
    value = "tab_results",
    div(
      style = "display:flex; gap:8px; align-items:flex-start;",
      
      # ── Caja filtros (redimensionable) ──────────────────────────────────────
      div(
        style = "resize:both; overflow:auto; min-width:160px; max-width:40vw; min-height:300px; max-height:90vh; width:220px; flex-shrink:0;",
        card(
          style = "height:100%; width:100%; min-height:0;",
          card_header(class = "bg-info text-white", bsicons::bs_icon("funnel-fill"), " Filtros"),
          card_body(
            style = "overflow-y:auto;",
            actionButton("btn_cargar_resultados", "🔄 Cargar / Actualizar", class = "btn-info btn-sm w-100"),
            hr(),
            radioButtons("filtro_modalidad", "Modalidad", choices  = c("CNVs", "SVs", "Ambas"), selected = "Ambas", inline = TRUE),
            hr(),
            selectizeInput("filtro_familia", "Familia(s)", choices  = NULL, multiple = TRUE, options  = list(placeholder = "Todas")),
            selectizeInput("filtro_herencia", "Tipo herencia", choices  = NULL, multiple = TRUE, options  = list(placeholder = "Todas")),
            selectizeInput("filtro_sv_type", "Tipo de variante", choices  = NULL, multiple = TRUE, options  = list(placeholder = "Todos")),
            selectizeInput("filtro_rango", "Tipo rango", choices  = c("Strict","Wide","Outside"), multiple = TRUE, options  = list(placeholder = "Todos")),
            selectInput("filtro_flagged", "🚩 Tags",
                        choices  = c("Todas" = "todas", "Con Tag" = "con", "Sin Tag" = "sin"),
                        selected = "todas"),
            selectInput("filtro_clasif", "🏷 Classification",
                        choices  = c("Todas" = "todas", "Sin clasificar" = "sin",
                                     "Benign" = "Benign", "VUS" = "VUS",
                                     "Pathogenic" = "Pathogenic", "In Progress" = "In Progress"),
                        selected = "todas"),
            selectInput("filtro_sexo", "♂/♀ Sexo",
                        choices  = c("All" = "all", "Male ♂" = "M", "Female ♀" = "F", "Unknown" = "ND"),
                        selected = "todos"),
            checkboxInput("filtro_sfari", "Gene panel only", value = FALSE),
            textInput("filtro_gen", "🔍 Search gene", placeholder = "E.g.: SHANK3, NRXN1...", width = "100%"),
            hr(),
            sliderInput("filtro_score", "AnnotSV Score (min.)", min = -1, max = 1, value = -1, step = 0.05),
            checkboxInput("filtro_solo_con_score", "🔢 SCORE — ocultar sin valor", value = FALSE),
            hr(),
            tags$details(
              tags$summary(style="cursor:pointer; font-weight:600; font-size:0.88em; color:#2C6FAC; user-select:none;",
                           "⚙ Visible columns"),
              div(class="mt-1",
                  checkboxGroupInput("sel_columnas", NULL,
                                     choices  = c("ID_Familia","SV_chrom","SV_start","SV_end","SV_length","SV_type",
                                                  "Tipo_Herencia","Genotipo_Hijo","Tipo_Rango","AnnotSV_ranking_score",
                                                  "En_Panel_Genes","Gene_name","Illumina_DRAGEN.exact.counts",
                                                  "Illumina_DRAGEN.similar.counts","Referencia_Match",
                                                  "Sim_MaxLen","Jaccard","Confianza_Region","Annotation_mode","Modalidad"),
                                     selected = c("ID_Familia","SV_chrom","SV_start","SV_end","SV_length","SV_type",
                                                  "Tipo_Herencia","Genotipo_Hijo","Tipo_Rango","AnnotSV_ranking_score",
                                                  "En_Panel_Genes","Gene_name","Illumina_DRAGEN.exact.counts",
                                                  "Illumina_DRAGEN.similar.counts","Referencia_Match",
                                                  "Sim_MaxLen","Jaccard","Confianza_Region")
                  ),
                  div(class="d-flex gap-1 mt-1",
                      actionButton("btn_aplicar_columnas", "✔ Aplicar", class="btn-primary btn-sm"),
                      actionButton("btn_reset_columnas",   "↺ Reset",   class="btn-outline-secondary btn-sm")
                  )
              )
            ),
            hr(),
            div(class = "d-grid gap-2",
                downloadButton("btn_descargar_xlsx", "⬇ Descargar Excel filtrado", class = "btn-primary btn-sm"),
                downloadButton("btn_descargar_csv",  "⬇ Descargar CSV filtrado", class = "btn-outline-secondary btn-sm"),
                hr(),
                tags$p(class="text-muted small mb-1", "📄 Clinical report per family"),
                selectInput("pdf_familia_sel", NULL, choices = c("— Selecciona familia —" = ""), width = "100%"),
                downloadButton("btn_pdf_familia", "📄 Generate PDF report", class = "btn-outline-primary btn-sm w-100")
            )
          )
        )
      ),
      
      # ── Columna derecha: tabla + panel detalle ──────────────────────────────
      div(
        style = "display:flex; flex-direction:column; gap:8px; flex:1; min-width:300px;",
        
        # Bloque tabla (redimensionable)
        div(
          style = "resize:both; overflow:auto; min-width:300px; max-width:90vw; min-height:200px; max-height:90vh;",
          card(
            style = "height:100%; width:100%; min-height:0;",
            card_header(
              class = "d-flex justify-content-between align-items-center",
              div(bsicons::bs_icon("grid-3x3"), " Tabla de variantes"),
              uiOutput("ui_tabla_info")
            ),
            card_body(
              class = "p-2",
              DTOutput("tabla_variantes", height = "820px")
            )
          )
        ),
        
        # Bloque detalle variante (solo visible al seleccionar fila, redimensionable)
        uiOutput("ui_panel_detalle")
      )
    )
  ),
  
  # ============================================================
  # TAB: GRAPHS (formerly Statistics)
  # ============================================================
  nav_panel(
    title = tagList(bsicons::bs_icon("bar-chart-line-fill"), " Graphs"),
    value = "tab_graficos",
    layout_columns(
      col_widths = c(2, 10),
      card(
        card_header(class = "bg-primary text-white",
                    bsicons::bs_icon("sliders"), " Opciones"),
        card_body(
          div(
            class = "p-2 rounded mb-2",
            style = "background:var(--dm-info-bg); border:1px solid var(--dm-info-border);",
            tags$b(class = "small d-block mb-1", "\U0001f4ca Fuente de datos"),
            tags$small(class = "text-muted",
                       "Los gr\u00e1ficos reflejan la selecci\u00f3n activa de ",
                       tags$b("Resultados"), ": modalidad, familias, herencia, rango, SFARI, gen, score\u2026")
          ),
          hr(),
          uiOutput("ui_est_resumen_filtro"),
          hr(),
          actionButton("btn_est_refresh", "\U0001f504 Recalcular", class = "btn-info btn-sm w-100"),
          div(class = "mt-2",
              actionButton("btn_ir_resultados", "\u2190 Ir a Resultados",
                           class = "btn-outline-secondary btn-sm w-100"))
        )
      ),
      div(
        style = "display:flex; flex-direction:column; gap:14px;",
        uiOutput("ui_est_kpis"),
        layout_columns(
          col_widths = c(4, 8),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("pie-chart-fill"), " Inheritance distribution"),
            card_body(class = "p-1", plotlyOutput("est_plot_herencia_pie", height = "300px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("bar-chart-steps"), " Variantes por familia (apilado por herencia)"),
            card_body(class = "p-1", plotlyOutput("est_plot_carga_familia", height = "300px"))
          )
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("bar-chart"), " AnnotSV score distribution"),
            card_body(class = "p-1", plotlyOutput("est_plot_scores", height = "300px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("speedometer2"), " Score por tipo de herencia (boxplot)"),
            card_body(class = "p-1", plotlyOutput("est_plot_score_her", height = "300px"))
          )
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("diagram-3"), " Variantes por cromosoma"),
            card_body(class = "p-1", plotlyOutput("est_plot_por_chr", height = "300px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("grid-3x2"), " Mapa cromosoma x familia"),
            card_body(class = "p-1", plotlyOutput("est_plot_heatmap_chr", height = "300px"))
          )
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("rulers"), " Size distribution (bp)"),
            card_body(class = "p-1", plotlyOutput("est_plot_tamano_hist", height = "300px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("distribute-vertical"), " Size by inheritance (boxplot)"),
            card_body(class = "p-1", plotlyOutput("est_plot_tamano_her", height = "300px"))
          )
        ),
        card(
          full_screen = TRUE,
          card_header(bsicons::bs_icon("graph-up-arrow"), " Score vs size (scatter)"),
          card_body(class = "p-1", plotlyOutput("est_plot_score_tamano", height = "300px"))
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("trophy"), " Top 25 most affected genes"),
            card_body(class = "p-1", plotlyOutput("est_plot_top_genes", height = "360px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("table"), " Resumen por gen (top 30)"),
            card_body(class = "p-1", DTOutput("est_tabla_genes", height = "340px"))
          )
        ),
        layout_columns(
          col_widths = c(4, 4, 4),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("bullseye"), " Tipo de rango"),
            card_body(class = "p-1", plotlyOutput("est_plot_rango", height = "260px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("copy"), " Tipo de variante (DEL / DUP)"),
            card_body(class = "p-1", plotlyOutput("est_plot_tipo_sv", height = "260px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("person-lines-fill"), " Score medio por familia"),
            card_body(class = "p-1", plotlyOutput("est_plot_score_familia", height = "260px"))
          )
        ),
        layout_columns(
          col_widths = c(4, 4, 4),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("gender-ambiguous"), " Distribuci\u00f3n por sexo"),
            card_body(class = "p-1", plotlyOutput("est_plot_sexo_pie", height = "260px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("bar-chart-fill"), " Carga de variantes por sexo"),
            card_body(class = "p-1", plotlyOutput("est_plot_carga_sexo", height = "260px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("speedometer"), " Score por sexo (boxplot)"),
            card_body(class = "p-1", plotlyOutput("est_plot_score_sexo", height = "260px"))
          )
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("bar-chart-steps"), " Herencia por sexo (apilado)"),
            card_body(class = "p-1", plotlyOutput("est_plot_herencia_sexo", height = "300px"))
          ),
          card(
            full_screen = TRUE,
            card_header(bsicons::bs_icon("rulers"), " Tama\u00f1o de variante por sexo (boxplot)"),
            card_body(class = "p-1", plotlyOutput("est_plot_tamano_sexo", height = "300px"))
          )
        )
      )
    )
  ),
  
  nav_menu(
    title = tagList(bsicons::bs_icon("layers"), " Visualisation"),
    align = "left",
    
    nav_panel(
      title = tagList(bsicons::bs_icon("map"), " Ideograma"),
      value = "tab_ideograma",
      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          div(bsicons::bs_icon("map"), " Interactive chromosomal ideogram — hg38"),
          div(
            class = "d-flex align-items-center gap-2",
            tags$span(class = "text-muted small", "Color by inheritance · Hover for detail · Scroll to zoom")
          )
        ),
        card_body(
          class = "p-2",
          plotlyOutput("plot_ideograma", height = "720px")
        )
      )
    ),
    
    
    nav_panel(
      title = tagList(bsicons::bs_icon("columns-gap"), " Comparador"),
      value = "tab_comparador",
      layout_columns(
        col_widths = c(3, 9),
        card(
          card_header(class="bg-primary text-white", bsicons::bs_icon("sliders"), " Selection"),
          card_body(
            div(class="d-flex align-items-center gap-1 mb-1",
                div(class="rounded-circle flex-shrink-0",
                    style="width:12px;height:12px;background:#2C6FAC;"),
                tags$b("Familia A", style="color:#2C6FAC;")),
            selectizeInput("comp_fam_a", NULL, choices=NULL, options=list(placeholder="Selecciona...")),
            div(class="d-flex align-items-center gap-1 mb-1 mt-2",
                div(class="rounded-circle flex-shrink-0",
                    style="width:12px;height:12px;background:#E67E22;"),
                tags$b("Familia B", style="color:#E67E22;")),
            selectizeInput("comp_fam_b", NULL, choices=NULL, options=list(placeholder="Selecciona...")),
            hr(),
            radioButtons("comp_modalidad", "Modalidad", choices=c("CNVs","SVs","Ambas"), selected="Ambas", inline=TRUE),
            numericInput("comp_umbral", "Umbral solapamiento (%)", value=50, min=10, max=100, step=5),
            hr(),
            actionButton("btn_comparar", "⚖ Compare families", class="btn-primary w-100"),
            hr(),
            uiOutput("ui_comp_resumen")
          )
        ),
        div(
          style="display:flex; flex-direction:column; gap:8px;",
          card(
            card_header(bsicons::bs_icon("bar-chart-steps"), " Comparative genomic view"),
            card_body(class="p-2",
                      uiOutput("ui_comp_placeholder"),
                      plotlyOutput("plot_comparador_genomico", height="520px")
            )
          ),
          card(
            card_header(
              class="d-flex justify-content-between align-items-center",
              div(bsicons::bs_icon("table"), " Detalle de variantes"),
              uiOutput("ui_comp_leyenda")
            ),
            card_body(class="p-2", DTOutput("tabla_comparador", height="340px"))
          )
        )
      )
    ),
    
    nav_panel(
      title = tagList(bsicons::bs_icon("diagram-3"), " Genes"),
      value = "tab_genes",
      layout_columns(
        col_widths = c(8, 4),
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center",
            div(bsicons::bs_icon("map"), " Ideograma de genes afectados"),
            div(class = "d-flex align-items-center gap-2",
                tags$span(class = "badge bg-secondary", uiOutput("ui_genes_badge", inline=TRUE)),
                checkboxInput("genes_solo_sfari", "Panel only", value = FALSE, width = "auto")
            )
          ),
          card_body(class = "p-2",
                    uiOutput("ui_genes_placeholder"),
                    plotlyOutput("plot_genes_ideograma", height = "700px")
          )
        ),
        card(
          card_header(bsicons::bs_icon("info-circle"), " Detalle del gen seleccionado"),
          card_body(class = "p-2",
                    uiOutput("ui_gen_detalle")
          )
        )
      )
    ),
    
    # ── HPOs ─────────────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(bsicons::bs_icon("clipboard2-pulse"), " HPOs"),
      value = "tab_hpos",
      tags$head(tags$script(HTML("
        function hpoSelGen(gen) {
          Shiny.setInputValue('hpo_gen_seleccionado', {gen: gen, t: Date.now()}, {priority: 'event'});
        }
      "))),
      layout_columns(
        col_widths = c(3, 6, 3),
        gap = "10px",
        card(
          height = "90vh",
          card_header(class = "bg-primary text-white d-flex align-items-center gap-2",
                      style = "font-size:0.9em;",
                      bsicons::bs_icon("people-fill"), " Familias y genes"),
          card_body(class = "p-2 d-flex flex-column gap-2", style = "overflow:auto;",
                    div(class = "mb-1",
                        tags$label(class = "form-label fw-bold small", "Seleccionar familias"),
                        selectizeInput("hpo_familias_sel", label = NULL, choices = NULL, multiple = TRUE,
                                       options = list(placeholder = "Elige una o varias familias…",
                                                      plugins = list("remove_button")), width = "100%")),
                    div(class = "d-flex justify-content-between align-items-center mb-1",
                        tags$span(class = "fw-bold small", "\U0001f9ec Genes afectados"),
                        uiOutput("ui_hpo_n_genes", inline = TRUE)),
                    div(class = "mb-1",
                        checkboxInput("hpo_solo_sfari", "Gene panel only", value = FALSE, width = "100%")),
                    div(class = "mb-1",
                        tags$input(id = "hpo_buscar_gen", class = "form-control form-control-sm",
                                   type = "text", placeholder = "Filtrar gen\u2026")),
                    uiOutput("ui_hpo_genes_lista")
          )
        ),
        card(
          height = "90vh",
          card_header(class = "d-flex align-items-center justify-content-between gap-2",
                      style = "font-size:0.9em;",
                      div(class = "d-flex align-items-center gap-2",
                          bsicons::bs_icon("clipboard2-pulse"), " Visor HPO"),
                      uiOutput("ui_hpo_visor_badge")),
          card_body(class = "p-2", style = "overflow:auto; height:calc(90vh - 54px);",
                    uiOutput("ui_hpo_visor"))
        ),
        card(
          height = "90vh",
          card_header(class = "d-flex align-items-center gap-2", style = "font-size:0.9em;",
                      bsicons::bs_icon("journal-text"), " Descripci\u00f3n HPO del paciente"),
          card_body(class = "p-2 d-flex flex-column gap-2", style = "overflow:auto;",
                    div(class = "mb-1",
                        tags$label(class = "form-label fw-bold small", "Familia activa para notas"),
                        selectizeInput("hpo_fam_nota_sel", label = NULL, choices = NULL,
                                       options = list(placeholder = "Selecciona una familia\u2026"), width = "100%")),
                    uiOutput("ui_hpo_nota_panel")
          )
        )
      )
    )
    
  ),  # end nav_menu Visualisation
  
  
  # ============================================================
  # TAB: STATISTICAL TESTS
  # ============================================================
  nav_panel(
    title = tagList(bsicons::bs_icon("calculator-fill"), " Statistical Tests"),
    value = "tab_tests",
    
    layout_columns(
      col_widths = c(3, 9),
      
      # Left panel: test configuration
      card(
        card_header(class = "bg-primary text-white",
                    bsicons::bs_icon("sliders2"), " Test configuration"),
        card_body(
          style = "overflow-y:auto; max-height:88vh;",
          
          # ── Selector de modo ──────────────────────────────────────────────
          div(class = "d-flex justify-content-center mb-2",
              radioButtons("test_modo_ui", NULL,
                           choices  = c("\U0001f9ed Guiado" = "guiado",
                                        "\u2699\ufe0f Experto"  = "experto"),
                           selected = "guiado", inline = TRUE)
          ),
          
          # ── Technical selector (always in DOM for reactivity, visible only in expert mode) ──
          conditionalPanel(
            condition = "input.test_modo_ui === \'experto\'",
            selectInput("test_tipo", "\U0001f52c Tipo de test",
                        choices = c(
                          "Wilcoxon / Mann-Whitney (2 grupos)"  = "wilcoxon",
                          "Kruskal-Wallis (>=3 grupos)"          = "kruskal",
                          "Chi-cuadrado de independencia"        = "chisq",
                          "Test exacto de Fisher"                = "fisher",
                          "Correlacion de Spearman"              = "spearman",
                          "Correlacion de Pearson"               = "pearson",
                          "Normalidad - Shapiro-Wilk"            = "shapiro",
                          "Kolmogorov-Smirnov (2 muestras)"      = "ks",
                          "Test binomial de proporciones"        = "prop_binom",
                          "Comparacion de proporciones 2 grupos" = "prop2",
                          "Carga por familia (Poisson)"          = "burden"
                        ),
                        selected = "wilcoxon", width = "100%"
            ),
            hr()
          ),
          
          # ── Wizard guiado ─────────────────────────────────────────────────
          uiOutput("ui_test_wizard_cards"),
          
          hr(),
          uiOutput("ui_test_params"),
          hr(),
          actionButton("btn_ejecutar_test", "\u25b6 Ejecutar test",
                       class = "btn-success w-100 fw-bold"),
          hr(),
          uiOutput("ui_test_ayuda")
        )
      ),
      
      # Panel derecho: resultados
      div(
        style = "display:flex; flex-direction:column; gap:12px;",
        uiOutput("ui_test_resultado"),
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center",
            div(bsicons::bs_icon("graph-up"), " Visualizacion del test"),
            uiOutput("ui_test_grafico_titulo")
          ),
          card_body(class = "p-1",
                    plotlyOutput("plot_test_resultado", height = "400px")
          )
        ),
        card(
          card_header(bsicons::bs_icon("table"), " Datos utilizados en el test"),
          card_body(class = "p-1",
                    DTOutput("tabla_test_datos", height = "220px")
          )
        )
      )
    )
  ),  # end nav_panel statistical tests
  
  nav_panel(
    title = tagList(bsicons::bs_icon("clock-history"), " Historial"),
    value = "tab_historial",
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center",
        div(bsicons::bs_icon("clock-history"), " Session action history"),
        actionButton("btn_limpiar_historial", "🗑 Clear history", class = "btn-outline-danger btn-sm")
      ),
      card_body(
        class = "p-2",
        DTOutput("tabla_historial", height = "600px")
      )
    )
  ),
  
  nav_panel(
    title = tagList(bsicons::bs_icon("info-circle"), " About"),
    value = "tab_about",
    
    div(
      class = "container-fluid py-3",
      style = "max-width:1100px; margin:0 auto;",
      
      # ── Cabecera ────────────────────────────────────────────────────────────
      div(
        class = "p-4 mb-4 rounded-3",
        style = "background: linear-gradient(135deg, #1A2980 0%, #2C6FAC 100%); color:white;",
        div(class = "d-flex align-items-center gap-3",
            div(style = "font-size:3em;", "\U0001f9ec"),
            div(
              h2(class = "mb-1", style = "font-weight:800; letter-spacing:-0.02em;",
                 "CNV/SV Analysis Pipeline 2014 FiltroB"),
              p(class = "mb-0", style = "opacity:0.85; font-size:1.05em;",
                "Interfaz de an\u00e1lisis automatizado de variantes de n\u00famero de copia (CNV) ",
                "y variantes estructurales (SV).")
            )
        )
      ),
      
      # ── Dos columnas principales ─────────────────────────────────────────────
      layout_columns(
        col_widths = c(7, 5),
        
        # ── Left column: usage guide ──────────────────────────────────
        div(
          
          # USAGE GUIDE
          card(
            card_header(
              class = "fw-bold dm-about-guide",
              bsicons::bs_icon("map"), " Gu\u00eda de uso paso a paso"
            ),
            card_body(
              class = "p-3",
              style = "font-size:0.92em;",
              
              # Paso 1
              div(class = "d-flex gap-3 mb-3",
                  div(class = "rounded-circle text-white d-flex align-items-center justify-content-center flex-shrink-0",
                      style = "width:32px;height:32px;background:#2C6FAC;font-weight:700;", "1"),
                  div(
                    tags$b("\u2699\ufe0f Configuraci\u00f3n del pipeline"),
                    tags$p(class="mb-0 text-muted",
                           "Ve a la pesta\u00f1a ", tags$b("Configuraci\u00f3n"), " e introduce las rutas de la",
                           " carpeta de entrada (con subcarpetas HI/PA/MA por familia), la carpeta de salida,",
                           " el archivo de rangos de referencia (.xlsx) y el panel de genes (.txt).",
                           " Ajusta los umbrales de similitud y el n\u00famero de n\u00facleos. Pulsa ",
                           tags$b("\U0001f4be Guardar configuraci\u00f3n"), " antes de continuar.")
                  )
              ),
              
              # Paso 2
              div(class = "d-flex gap-3 mb-3",
                  div(class = "rounded-circle text-white d-flex align-items-center justify-content-center flex-shrink-0",
                      style = "width:32px;height:32px;background:#2C6FAC;font-weight:700;", "2"),
                  div(
                    tags$b("\u25b6\ufe0f Ejecuci\u00f3n del an\u00e1lisis"),
                    tags$p(class="mb-0 text-muted",
                           "Ve a ", tags$b("Ejecuci\u00f3n"), ". Puedes lanzar la cohorte completa,",
                           " una carpeta batch concreta o una familia individual (introduce su ID, p.\u00a0ej. ",
                           tags$code("G06PlS056"), ").",
                           " El log de ejecuci\u00f3n se actualiza en tiempo real.",
                           " Usa ", tags$b("\u23f9 Detener"), " si necesitas cancelar el proceso.")
                  )
              ),
              
              # Paso 3
              div(class = "d-flex gap-3 mb-3",
                  div(class = "rounded-circle text-white d-flex align-items-center justify-content-center flex-shrink-0",
                      style = "width:32px;height:32px;background:#2C6FAC;font-weight:700;", "3"),
                  div(
                    tags$b("\U0001f504 Carga de resultados"),
                    tags$p(class="mb-0 text-muted",
                           "En la pesta\u00f1a ", tags$b("Resultados"), ", pulsa ",
                           tags$b("\U0001f504 Cargar / Actualizar"), " para importar todos los archivos",
                           " ", tags$code("_Complete_Analysis.xlsx"), " encontrados en la carpeta de salida.",
                           " La app detectar\u00e1 autom\u00e1ticamente variantes similares entre individuos (umbral configurable).")
                  )
              ),
              
              # Paso 4
              div(class = "d-flex gap-3 mb-3",
                  div(class = "rounded-circle text-white d-flex align-items-center justify-content-center flex-shrink-0",
                      style = "width:32px;height:32px;background:#2C6FAC;font-weight:700;", "4"),
                  div(
                    tags$b("\U0001f50d Filtrado y revisi\u00f3n"),
                    tags$p(class="mb-0 text-muted",
                           "Filtra por modalidad (CNV/SV), herencia, tipo de rango, score AnnotSV, presencia en panel de genes,",
                           " gen, clasificaci\u00f3n cl\u00ednica o Tag. Haz clic en una fila para ver el",
                           " panel de detalle con acceso directo a FRANKLIN, ClinVar, gnomAD,",
                           " DECIPHER, VarSome, UCSC, Ensembl y OMIM.")
                  )
              ),
              
              # Paso 5
              div(class = "d-flex gap-3 mb-3",
                  div(class = "rounded-circle text-white d-flex align-items-center justify-content-center flex-shrink-0",
                      style = "width:32px;height:32px;background:#2C6FAC;font-weight:700;", "5"),
                  div(
                    tags$b("\U0001f3f7\ufe0f Anotaci\u00f3n cl\u00ednica"),
                    tags$p(class="mb-0 text-muted",
                           "En el panel de detalle puedes: asignar una ",
                           tags$b("clasificaci\u00f3n cl\u00ednica"),
                           " (Ben\u00edgna / VUS / Patog\u00e9nica / En Progreso), a\u00f1adir ",
                           tags$b("notas libres"), " (con opci\u00f3n de replicar a variantes similares),",
                           " marcar con ", tags$b("Tag \U0001f6a9"), " para seguimiento,",
                           " y registrar el ", tags$b("sexo del probando"), " (\u2642/\u2640).",
                           " Todo se guarda autom\u00e1ticamente en disco.")
                  )
              ),
              
              # Paso 6
              div(class = "d-flex gap-3 mb-3",
                  div(class = "rounded-circle text-white d-flex align-items-center justify-content-center flex-shrink-0",
                      style = "width:32px;height:32px;background:#2C6FAC;font-weight:700;", "6"),
                  div(
                    tags$b("\U0001f4ca Visualizaci\u00f3n y estad\u00edsticos"),
                    tags$p(class="mb-0 text-muted",
                           "Explora la cohorte en la pesta\u00f1a ", tags$b("Gr\u00e1ficos"),
                           " (KPIs, distribuciones, heatmap cromosoma\u00d7familia, top genes, etc.).",
                           " Usa ", tags$b("Tests Estad\u00edsticos"), " para realizar an\u00e1lisis formales",
                           " (Wilcoxon, Kruskal-Wallis, Chi-cuadrado, Spearman, Shapiro-Wilk, burden\u2026).",
                           " La pesta\u00f1a ", tags$b("Visualizaci\u00f3n"),
                           " incluye ideograma hg38, comparador de familias e ideograma de genes.")
                  )
              ),
              
              # Paso 7
              div(class = "d-flex gap-3",
                  div(class = "rounded-circle text-white d-flex align-items-center justify-content-center flex-shrink-0",
                      style = "width:32px;height:32px;background:#2C6FAC;font-weight:700;", "7"),
                  div(
                    tags$b("\u2b07\ufe0f Exportaci\u00f3n"),
                    tags$p(class="mb-0 text-muted",
                           "Desde la pesta\u00f1a Resultados puedes descargar las variantes filtradas en formato ",
                           tags$b("Excel (.xlsx)"), " o ", tags$b("CSV"),
                           ". Tambi\u00e9n puedes generar un ", tags$b("informe PDF cl\u00ednico"),
                           " por familia, que incluye variantes, anotaciones, notas y clasificaciones.")
                  )
              )
            )
          ),
          
          # ESTRUCTURA DE CARPETAS
          card(
            class = "mt-3",
            card_header(
              class = "fw-bold dm-about-success",
              bsicons::bs_icon("folder2-open"), " Estructura de carpetas esperada"
            ),
            card_body(
              class = "p-3",
              tags$pre(
                class = "p-3 rounded small",
                class = "dm-about-code p-2 rounded", style="font-size:0.82em; line-height:1.6;",
                "Carpeta de entrada/\n",
                "\u251c\u2500\u2500 Batch_01/\n",
                "\u2502   \u251c\u2500\u2500 G06PlS056HI/      \u2190 probando (hijo)\n",
                "\u2502   \u251c\u2500\u2500 G06PlS056PA/      \u2190 padre\n",
                "\u2502   \u2514\u2500\u2500 G06PlS056MA/      \u2190 madre\n",
                "\u251c\u2500\u2500 Batch_02/\n",
                "\u2502   \u251c\u2500\u2500 G07XY001HI/\n",
                "\u2502   \u251c\u2500\u2500 069-070PA/        \u2190 progenitor compartido\n",
                "\u2502   \u2514\u2500\u2500 069-070MA/\n",
                "...\n\n",
                "Carpeta de salida/\n",
                "\u251c\u2500\u2500 G06PlS056_Analisis_Completo.xlsx\n",
                "\u251c\u2500\u2500 G07XY001_Analisis_Completo.xlsx\n",
                "\u251c\u2500\u2500 .flags_variantes.rds\n",
                "\u251c\u2500\u2500 .notas_variantes.rds\n",
                "\u251c\u2500\u2500 .clasificaciones_variantes.rds\n",
                "\u2514\u2500\u2500 .sexos_familias.rds"
              ),
              tags$small(class="text-muted",
                         "Los archivos ocultos (", tags$code(".rds"), ") contienen las anotaciones cl\u00ednicas ",
                         "y se cargan autom\u00e1ticamente junto con los resultados.")
            )
          )
        ),
        
        # ── Right column: legends and technical reference ──────────────────
        div(
          
          # LEYENDA HERENCIA
          card(
            card_header(
              class = "fw-bold dm-about-warn",
              bsicons::bs_icon("diagram-2"), " Leyenda: Tipo de herencia"
            ),
            card_body(class = "p-3",
                      lapply(list(
                        list("#CC0000", "De novo",             "Variante ausente en ambos progenitores. Mayor relevancia cl\u00ednica."),
                        list("#0066CC", "Paternal",             "Heredada del padre (confirmada)."),
                        list("#CC6600", "Maternal",             "Heredada de la madre (confirmada)."),
                        list("#6600CC", "Combined",            "Presente en ambos progenitores."),
                        list("#4499DD", "Paternal (Probable)",  "Probablemente paterna; similitud alta pero no exacta."),
                        list("#DD9944", "Maternal (Probable)",  "Probablemente materna; similitud alta pero no exacta."),
                        list("#9944DD", "Combined (Probable)", "Probablemente en ambos; similitudes detectadas."),
                        list("#6c757d", "Desconocido",         "Sin datos de progenitores o similitud insuficiente.")
                      ), function(x) {
                        div(class = "d-flex align-items-start gap-2 mb-2",
                            div(style = paste0("width:14px;height:14px;border-radius:3px;background:",x[[1]],
                                               ";flex-shrink:0;margin-top:3px;")),
                            div(tags$b(x[[2]]), tags$br(),
                                tags$span(class="text-muted", style="font-size:0.82em;", x[[3]]))
                        )
                      })
            )
          ),
          
          # RANGE + CLASSIFICATION LEGEND
          card(
            class = "mt-3",
            card_header(
              class = "fw-bold dm-about-success",
              bsicons::bs_icon("bullseye"), " Leyenda: Tipo de rango y clasificaci\u00f3n"
            ),
            card_body(class = "p-3",
                      tags$b("Tipo de rango"),
                      tags$hr(class="my-1"),
                      lapply(list(
                        list("#155724", "Strict", "La variante solapa completamente con un rango de referencia."),
                        list("#856404", "Wide",   "La variante solapa parcialmente (supera el umbral Jaccard)."),
                        list("#9C0006", "Outside",    "Sin solapamiento significativo con ning\u00fan rango de referencia.")
                      ), function(x) {
                        div(class="d-flex align-items-start gap-2 mb-2",
                            div(style=paste0("width:14px;height:14px;border-radius:3px;background:",x[[1]],";flex-shrink:0;margin-top:3px;")),
                            div(tags$b(x[[2]]),tags$br(),
                                tags$span(class="text-muted",style="font-size:0.82em;",x[[3]]))
                        )
                      }),
                      tags$b(class="mt-2 d-block","Clasificaci\u00f3n cl\u00ednica"),
                      tags$hr(class="my-1"),
                      lapply(list(
                        list("#27AE60", "\U0001f7e2 Ben\u00edgna",      "Sin relevancia cl\u00ednica aparente."),
                        list("#95A5A6", "\u2b1c VUS",            "Variante de significado incierto."),
                        list("#E74C3C", "\U0001f534 Patog\u00e9nica", "Relevancia cl\u00ednica confirmada o altamente probable."),
                        list("#8E44AD", "\U0001f7e3 En Progreso",  "Bajo revisi\u00f3n o pendiente de datos adicionales.")
                      ), function(x) {
                        div(class="d-flex align-items-start gap-2 mb-2",
                            div(style=paste0("width:14px;height:14px;border-radius:3px;background:",x[[1]],";flex-shrink:0;margin-top:3px;")),
                            div(tags$b(x[[2]]),tags$br(),
                                tags$span(class="text-muted",style="font-size:0.82em;",x[[3]]))
                        )
                      })
            )
          ),
          
          # COLUMNAS CLAVE
          card(
            class = "mt-3",
            card_header(
              class = "fw-bold dm-about-purple",
              bsicons::bs_icon("table"), " Columnas clave de resultados"
            ),
            card_body(class="p-3",
                      style="font-size:0.83em;",
                      tags$table(
                        class="table table-sm table-hover",
                        style="margin-bottom:0;",
                        tags$thead(tags$tr(tags$th("Columna"),tags$th("Descripci\u00f3n"))),
                        tags$tbody(
                          tags$tr(tags$td(tags$code("SV_chrom/start/end")),tags$td("Coordenadas gen\u00f3micas (hg38)")),
                          tags$tr(tags$td(tags$code("SV_length")),tags$td("Tama\u00f1o de la variante en pares de bases")),
                          tags$tr(tags$td(tags$code("SV_type")),tags$td("Tipo: DEL (deleci\u00f3n), DUP (duplicaci\u00f3n), INV\u2026")),
                          tags$tr(tags$td(tags$code("Tipo_Herencia")),tags$td("Patr\u00f3n de herencia calculado por el pipeline")),
                          tags$tr(tags$td(tags$code("Tipo_Rango")),tags$td("Solapamiento con rangos de referencia (Estricto/Amplio/Fuera)")),
                          tags$tr(tags$td(tags$code("AnnotSV_ranking_score")),tags$td("Score de patogenicidad AnnotSV (-2\u20132; >0.5 = probable patog\u00e9nica)")),
                          tags$tr(tags$td(tags$code("En_Panel_Genes")),tags$td("Si alg\u00fan gen afectado pertenece al panel SFARI de genes de autismo")),
                          tags$tr(tags$td(tags$code("Gene_name")),tags$td("Genes de OMIM anotados en la regi\u00f3n")),
                          tags$tr(tags$td(tags$code("Referencia_Match")),tags$td("Solapamiento con la base de referencia interna")),
                          tags$tr(tags$td(tags$code("Jaccard")),tags$td("\u00cdndice de similitud Jaccard con el rango de referencia")),
                          tags$tr(tags$td(tags$code("Sim_MaxLen")),tags$td("M\u00e1xima similitud por longitud frente a otros individuos")),
                          tags$tr(tags$td(tags$code("Confianza_Region")),tags$td("Nivel de confianza de la asignaci\u00f3n de rango")),
                          tags$tr(tags$td(tags$code("Annotation_mode")),tags$td(tags$span(tags$code("full"), ": variante completa; ",
                                                                                          tags$code("split"), ": fragmento de variante anotado por gen")))
                        )
                      )
            )
          ),
          
          # TECHNICAL STACK
          card(
            class = "mt-3",
            card_header(
              class = "fw-bold dm-about-muted",
              bsicons::bs_icon("cpu"), " Stack t\u00e9cnico"
            ),
            card_body(class="p-3",
                      div(class="d-flex flex-wrap gap-2",
                          lapply(c(
                            "R \u2265 4.3","Shiny","bslib","DT","plotly","dplyr",
                            "openxlsx","processx","GenomicRanges","AnnotSV","scales",
                            "grid","bsicons"
                          ), function(pkg) {
                            tags$span(class="badge",
                                      style="background:#2C6FAC; font-size:0.85em; font-weight:500;",
                                      pkg)
                          })
                      ),
                      hr(),
                      tags$small(class="text-muted",
                                 "\U0001f4cc El pipeline de an\u00e1lisis (", tags$code("FiltroB_STEA.R"),
                                 ") se ejecuta en un proceso R separado v\u00eda ",
                                 tags$code("processx"), ", lo que permite monitorizar el progreso en",
                                 " tiempo real sin bloquear la interfaz."
                      )
            )
          )
        )
      )
    )
  ),
  
  nav_spacer(),
  nav_item(uiOutput("ui_badge_global"))
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {
  COLOR_A <- "#2C6FAC"
  COLOR_B <- "#E67E22"
  COLOR_AB <- "#27AE60"
  
  # ── Carpeta centralizada de logs ─────────────────────────────────────────
  APP_LOGS_DIR <- file.path(getwd(), "logs")
  if (!dir.exists(APP_LOGS_DIR)) dir.create(APP_LOGS_DIR, recursive = TRUE)
  # ─────────────────────────────────────────────────────────────────────────
  rv <- reactiveValues(
    proceso           = NULL,
    log_lines         = character(0),
    pipeline_activo   = FALSE,
    resultados        = list(CNVs = NULL, SVs = NULL),
    config_guardada   = CONFIG_DEFAULTS,
    fila_seleccionada = NULL,
    fila_data         = NULL,
    flags             = character(0),
    flags_auto        = character(0),
    clasificaciones   = list(),      # lista nombrada: clave -> "Benign"|"VUS"|"Pathogenic"|"In Progress"
    notas             = list(),      # lista nombrada: clave -> texto
    historial         = list(),      # lista de eventos: timestamp, tipo, clave, detalle
    sexos               = list(),      # lista nombrada: ID_Familia -> "M" | "F" | NULL
    fenotipos           = list(),      # lista nombrada: ID_Familia -> texto libre de fenotipo
    hpo_notas           = list(),      # named list: ID_Family -> patient HPO description
    pendiente_seleccion = NULL         # list(fam, sta, end) — seleccion pendiente en tabla
  )
  # ── GUARDAR EL TEXTO "FANTASMA" MIENTRAS SE TECLEA ──
  observeEvent(input$input_nota_variante, {
    rv$nota_fantasma <- input$input_nota_variante
  }, ignoreInit = TRUE)
  
  observeEvent(input$input_fenotipo_familia, {
    rv$feno_fantasma <- input$input_fenotipo_familia
  }, ignoreInit = TRUE)
  
  # When switching to a DIFFERENT variant, we do want to clear stale values 
  # para que no se pegue el texto a la nueva variante
  observeEvent(rv$fila_data, {
    rv$nota_fantasma <- NULL
    rv$feno_fantasma <- NULL
  }, priority = 10)
  # ── Initial load at session start (runs once) ──────────────────
  local({
    cargado <- FALSE
    observe({
      if (cargado) return()
      cargado <<- TRUE
      
      # Cargar flags
      ruta_fl <- file.path(APP_LOGS_DIR, ".flags_variantes.rds")
      if (file.exists(ruta_fl)) {
        saved <- tryCatch(readRDS(ruta_fl), error = function(e) NULL)
        if (!is.null(saved) && is.character(saved)) {
          n_pipes  <- lengths(regmatches(saved, gregexpr("|", saved, fixed = TRUE)))
          sin_modo <- saved[n_pipes == 4L]
          con_modo <- if (length(sin_modo) > 0) paste0(sin_modo, "|FULL") else character(0)
          rv$flags <- union(rv$flags, c(saved[n_pipes != 4L], con_modo))
        }
      }
      
      # Cargar notas
      ruta_n <- file.path(APP_LOGS_DIR, ".notas_variantes.rds")
      if (file.exists(ruta_n)) {
        saved_n <- tryCatch(readRDS(ruta_n), error = function(e) NULL)
        if (!is.null(saved_n) && is.list(saved_n)) rv$notas <- saved_n
      }
      ruta_sx <- file.path(APP_LOGS_DIR, ".sexos_familias.rds")
      if (file.exists(ruta_sx)) {
        saved_sx <- tryCatch(readRDS(ruta_sx), error = function(e) NULL)
        if (!is.null(saved_sx) && is.list(saved_sx)) rv$sexos <- saved_sx
      }
      ruta_cl <- file.path(APP_LOGS_DIR, ".clasificaciones_variantes.rds")
      if (file.exists(ruta_cl)) {
        saved_cl <- tryCatch(readRDS(ruta_cl), error = function(e) NULL)
        if (!is.null(saved_cl) && is.list(saved_cl)) rv$clasificaciones <- saved_cl
      }
      ruta_ft <- file.path(APP_LOGS_DIR, ".fenotipos_familias.rds")
      if (file.exists(ruta_ft)) {
        saved_ft <- tryCatch(readRDS(ruta_ft), error = function(e) NULL)
        if (!is.null(saved_ft) && is.list(saved_ft)) rv$fenotipos <- saved_ft
      }
      ruta_hn <- file.path(APP_LOGS_DIR, ".hpo_notas_familias.rds")
      if (file.exists(ruta_hn)) {
        saved_hn <- tryCatch(readRDS(ruta_hn), error = function(e) NULL)
        if (!is.null(saved_hn) && is.list(saved_hn)) rv$hpo_notas <- saved_hn
      }
    })
  })
  
  # ── RUTAS CENTRALIZADAS A LA CARPETA LOGS ──
  ruta_flags_file     <- reactive({ file.path(APP_LOGS_DIR, ".flags_variantes.rds") })
  ruta_notas_file     <- reactive({ file.path(APP_LOGS_DIR, ".notas_variantes.rds") })
  ruta_sexos_file     <- reactive({ file.path(APP_LOGS_DIR, ".sexos_familias.rds") })
  ruta_clasif_file    <- reactive({ file.path(APP_LOGS_DIR, ".clasificaciones_variantes.rds") })
  ruta_fenotipos_file <- reactive({ file.path(APP_LOGS_DIR, ".fenotipos_familias.rds") })
  ruta_hpo_notas_file <- reactive({ file.path(APP_LOGS_DIR, ".hpo_notas_familias.rds") })
  
  # Auto-guardado reactivo: cada vez que rv$flags cambia, persiste en disco
  observe({
    flags_actuales <- rv$flags
    ruta_fl <- tryCatch(ruta_flags_file(), error = function(e) NULL)
    if (!is.null(ruta_fl) && nchar(trimws(ruta_fl)) > 0) {
      tryCatch({
        if (!dir.exists(dirname(ruta_fl))) dir.create(dirname(ruta_fl), recursive = TRUE)
        saveRDS(flags_actuales, ruta_fl)
      }, error = function(e) {
        showNotification(paste0("⚠️ Could not save Tags: ", e$message), type = "warning", duration = 5)
      })
    }
  })
  
  # Auto-guardado reactivo de notas
  observe({
    notas_actuales <- rv$notas
    ruta_n <- tryCatch(ruta_notas_file(), error = function(e) NULL)
    if (!is.null(ruta_n) && nchar(trimws(ruta_n)) > 0) {
      tryCatch({
        if (!dir.exists(dirname(ruta_n))) dir.create(dirname(ruta_n), recursive = TRUE)
        saveRDS(notas_actuales, ruta_n)
      }, error = function(e) {})
    }
  })
  observe({
    sexos_actuales <- rv$sexos
    ruta_sx <- tryCatch(ruta_sexos_file(), error = function(e) NULL)
    if (!is.null(ruta_sx) && nchar(trimws(ruta_sx)) > 0) {
      tryCatch({
        if (!dir.exists(dirname(ruta_sx))) dir.create(dirname(ruta_sx), recursive = TRUE)
        saveRDS(sexos_actuales, ruta_sx)
      }, error = function(e) {})
    }
  })
  # Auto-guardado reactivo de clasificaciones
  observe({
    cl <- rv$clasificaciones
    ruta_cl <- tryCatch(ruta_clasif_file(), error = function(e) NULL)
    if (!is.null(ruta_cl) && nchar(trimws(ruta_cl)) > 0) {
      tryCatch({
        if (!dir.exists(dirname(ruta_cl))) dir.create(dirname(ruta_cl), recursive = TRUE)
        saveRDS(cl, ruta_cl)
      }, error = function(e) {})
    }
  })
  # Auto-guardado reactivo de fenotipos
  observe({
    ft <- rv$fenotipos
    ruta_ft <- tryCatch(ruta_fenotipos_file(), error = function(e) NULL)
    if (!is.null(ruta_ft) && nchar(trimws(ruta_ft)) > 0) {
      tryCatch({
        if (!dir.exists(dirname(ruta_ft))) dir.create(dirname(ruta_ft), recursive = TRUE)
        saveRDS(ft, ruta_ft)
      }, error = function(e) {})
    }
  })
  # Auto-guardado reactivo de notas HPO
  observe({
    hn <- rv$hpo_notas
    ruta_hn <- tryCatch(ruta_hpo_notas_file(), error = function(e) NULL)
    if (!is.null(ruta_hn) && nchar(trimws(ruta_hn)) > 0) {
      tryCatch({
        if (!dir.exists(dirname(ruta_hn))) dir.create(dirname(ruta_hn), recursive = TRUE)
        saveRDS(hn, ruta_hn)
      }, error = function(e) {})
    }
  })
  
  observeEvent(input$btn_reset_config, {
    updateTextInput(session,   "cfg_entrada",       value = CONFIG_DEFAULTS$ruta_entrada)
    updateTextInput(session,   "cfg_salida",        value = CONFIG_DEFAULTS$ruta_salida)
    updateTextInput(session,   "cfg_rangos",        value = CONFIG_DEFAULTS$ruta_rangos)
    updateTextInput(session,   "cfg_panel",         value = CONFIG_DEFAULTS$ruta_panel)
    updateNumericInput(session,"cfg_umbral_strict", value = CONFIG_DEFAULTS$umbral_strict)
    updateNumericInput(session,"cfg_umbral_her",    value = CONFIG_DEFAULTS$umbral_herencia)
    updateNumericInput(session,"cfg_umbral_her_lax",value = CONFIG_DEFAULTS$umbral_herencia_lax)
    updateNumericInput(session,"cfg_umbral_jaccard",value = CONFIG_DEFAULTS$umbral_jaccard)
    updateNumericInput(session,"cfg_margen",        value = CONFIG_DEFAULTS$margen_lateral / 1e6)
    updateNumericInput(session,"cfg_cores",         value = CONFIG_DEFAULTS$n_cores)
    updateCheckboxInput(session,"cfg_dedup",        value = CONFIG_DEFAULTS$deduplicar)
    showNotification("Configuration restored to default values.", type = "message")
  })
  
  observeEvent(input$btn_guardar_config, {
    rv$config_guardada <- list(
      ruta_entrada        = input$cfg_entrada,
      ruta_salida         = input$cfg_salida,
      ruta_rangos         = input$cfg_rangos,
      ruta_panel          = input$cfg_panel,
      umbral_strict       = input$cfg_umbral_strict,
      umbral_herencia     = input$cfg_umbral_her,
      umbral_herencia_lax = input$cfg_umbral_her_lax,
      umbral_jaccard      = input$cfg_umbral_jaccard,
      margen_lateral      = input$cfg_margen * 1e6,
      deduplicar          = input$cfg_dedup,
      n_cores             = input$cfg_cores
    )
    showNotification("✅ Configuration saved.", type = "message", duration = 3)
  })
  
  output$ui_validacion_rutas <- renderUI({
    req(input$btn_validar_rutas)
    isolate({
      checks <- list(
        list(label = "Carpeta entrada",  ok = dir.exists(input$cfg_entrada)),
        list(label = "Carpeta salida",   ok = dir.exists(input$cfg_salida) ||
               tryCatch({ dir.create(input$cfg_salida, recursive=TRUE); TRUE }, error=function(e) FALSE)),
        list(label = "Archivo rangos",   ok = file.exists(input$cfg_rangos)),
        list(label = "Gene panel",      ok = file.exists(input$cfg_panel))
      )
      tags$ul(class = "list-unstyled mt-2", lapply(checks, function(ch) {
        tags$li(if (ch$ok) tags$span(class="text-success", "✔ ", ch$label) else tags$span(class="text-danger",  "✖ ", ch$label, " — no encontrado"))
      }))
    })
  })
  
  construir_script_run <- function(id_familia = NULL, ruta_carpeta = NULL) {
    cfg <- rv$config_guardada
    # Normalizar todas las rutas a barras forward para evitar secuencias de escape en el script generado
    cfg$ruta_entrada <- gsub("\\", "/", cfg$ruta_entrada, fixed = TRUE)
    cfg$ruta_salida  <- gsub("\\", "/", cfg$ruta_salida,  fixed = TRUE)
    cfg$ruta_rangos  <- gsub("\\", "/", cfg$ruta_rangos,  fixed = TRUE)
    cfg$ruta_panel   <- gsub("\\", "/", cfg$ruta_panel,   fixed = TRUE)
    lineas_config <- paste0(
      'CONFIG <- list(\n',
      '  ruta_entrada        = ', shQuote(cfg$ruta_entrada),        ',\n',
      '  ruta_salida         = ', shQuote(cfg$ruta_salida),         ',\n',
      '  ruta_rangos         = ', shQuote(cfg$ruta_rangos),         ',\n',
      '  ruta_panel          = ', shQuote(cfg$ruta_panel),          ',\n',
      '  umbral_strict       = ', cfg$umbral_strict,                ',\n',
      '  umbral_herencia     = ', cfg$umbral_herencia,              ',\n',
      '  umbral_herencia_lax = ', cfg$umbral_herencia_lax,          ',\n',
      '  umbral_jaccard      = ', cfg$umbral_jaccard,               ',\n',
      '  margen_lateral      = ', cfg$margen_lateral,               ',\n',
      '  deduplicar_regiones = ', toupper(cfg$deduplicar),          ',\n',
      '  n_cores             = ', cfg$n_cores,                      '\n',
      ')\n'
    )
    
    script_lines <- readLines(PIPELINE_SCRIPT, warn = FALSE)
    idx_insercion <- grep("procesar_tarea <- function", script_lines)[1]
    
    if (!is.null(id_familia) && nchar(trimws(id_familia)) > 0) {
      # Filtrar por familia concreta
      id_buscado <- trimws(id_familia)
      if (!is.na(idx_insercion)) {
        filtro_lines <- c(
          '',
          '# --- BLOQUE INYECTADO POR LA APP SHINY ---',
          paste0('id_objetivo <- "', id_buscado, '"'),
          'cat("\\n[APP] 🔎 Filtrando para procesar solo la familia:", id_objetivo, "\\n")',
          'if (exists("tareas")) {',
          '  tareas <- tareas[id_familia == id_objetivo]',
          '  cat("[APP] 📋 Tasks found after filtering:", nrow(tareas), "\\n")',
          '  if(nrow(tareas) == 0) {',
          '    cat("[APP] ❌ ERROR: El ID proporcionado no se encuentra en la ruta de entrada.\\n")',
          '    stop("ID de familia no reconocido por la estructura del script.")',
          '  }',
          '}',
          '# ------------------------------------------',
          ''
        )
        script_lines <- append(script_lines, filtro_lines, after = idx_insercion - 1)
      }
    } else if (!is.null(ruta_carpeta) && nchar(trimws(ruta_carpeta)) > 0) {
      # Filtrar por carpeta: usar nombre_sub igual que familia usa id_familia
      carpeta_norm   <- gsub("\\", "/", trimws(ruta_carpeta), fixed = TRUE)
      carpeta_nombre <- basename(gsub("/$", "", carpeta_norm))
      if (!is.na(idx_insercion)) {
        filtro_lines <- c(
          '',
          '# --- BLOQUE INYECTADO POR LA APP SHINY ---',
          paste0('carpeta_objetivo <- "', carpeta_nombre, '"'),
          'cat("\\n[APP] 📂 Filtering to process only subfolder:", carpeta_objetivo, "\\n")',
          'if (exists("tareas")) {',
          '  tareas <- tareas[nombre_sub == carpeta_objetivo]',
          '  cat("[APP] 📋 Tasks found after filtering:", nrow(tareas), "\\n")',
          '  if(nrow(tareas) == 0) {',
          '    cat("[APP] ❌ ERROR: No se encontraron familias con nombre_sub ==", carpeta_objetivo, "\\n")',
          '    stop("Subcarpeta no reconocida por la estructura del script.")',
          '  }',
          '}',
          '# ------------------------------------------',
          ''
        )
        script_lines <- append(script_lines, filtro_lines, after = idx_insercion - 1)
      }
    }
    
    tmp <- tempfile(fileext = ".R")
    writeLines(c(lineas_config, script_lines), tmp)
    tmp
  }
  
  lanzar_pipeline <- function(id_familia = NULL, ruta_carpeta = NULL) {
    if (rv$pipeline_activo) {
      showNotification("An analysis is already running.", type = "warning")
      return()
    }
    rv$log_lines <- character(0)
    rv$pipeline_activo <- TRUE
    shinyjs::disable("btn_run_todo")
    shinyjs::disable("btn_run_carpeta")
    shinyjs::disable("btn_run_familia")
    shinyjs::enable("btn_stop")
    
    script_tmp <- construir_script_run(id_familia, ruta_carpeta)
    rv$proceso <- processx::process$new(command = "Rscript", args = script_tmp, stdout = "|", stderr = "|")
    showNotification("🚀 Pipeline started.", type = "message", duration = 3)
  }
  
  observeEvent(input$btn_run_todo,    { lanzar_pipeline(NULL) })
  observeEvent(input$btn_run_carpeta, {
    ruta <- trimws(input$inp_carpeta_ruta)
    if (!nchar(ruta) || !dir.exists(ruta)) {
      showNotification("❌ Please enter a valid folder path.", type = "error", duration = 5)
      return()
    }
    lanzar_pipeline(id_familia = NULL, ruta_carpeta = ruta)
  })
  observeEvent(input$btn_run_familia, { lanzar_pipeline(input$inp_familia_id) })
  
  observeEvent(input$btn_stop, {
    if (!is.null(rv$proceso) && rv$proceso$is_alive()) {
      rv$proceso$kill()
      rv$log_lines <- c(rv$log_lines, "\n⏹ Proceso detenido por el usuario.")
    }
    rv$pipeline_activo <- FALSE
    shinyjs::enable("btn_run_todo")
    shinyjs::enable("btn_run_carpeta")
    shinyjs::enable("btn_run_familia")
    shinyjs::disable("btn_stop")
  })
  
  observe({
    req(rv$pipeline_activo, !is.null(rv$proceso))
    invalidateLater(600, session)
    
    proc <- rv$proceso
    nuevas <- tryCatch({ c(proc$read_output_lines(), proc$read_error_lines()) }, error = function(e) character(0))
    if (length(nuevas) > 0) rv$log_lines <- c(rv$log_lines, nuevas)
    
    if (!proc$is_alive()) {
      rv$pipeline_activo <- FALSE
      shinyjs::enable("btn_run_todo")
      shinyjs::enable("btn_run_carpeta")
      shinyjs::enable("btn_run_familia")
      shinyjs::disable("btn_stop")
      rv$log_lines <- c(rv$log_lines, paste0("\n✅ Process finished — exit code: ", proc$get_exit_status()))
      showNotification("✅ Pipeline completado.", type = "message", duration = 5)
    }
  })
  
  output$log_output <- renderText({ paste(tail(rv$log_lines, 400), collapse = "\n") })
  
  observeEvent(rv$log_lines, {
    shinyjs::runjs("var el = document.getElementById('log_container'); if(el) el.scrollTop = el.scrollHeight;")
  })
  
  output$ui_run_status <- renderUI({
    if (rv$pipeline_activo) {
      tagList(tags$span(class = "spinner-border spinner-border-sm text-success"), tags$span(" Analysis running...", class = "text-success fw-bold ms-2"))
    } else if (length(rv$log_lines) > 0) {
      last <- tail(rv$log_lines, 1)
      color <- if (grepl("✅|FINALIZADO", last)) "text-success" else "text-secondary"
      tags$span(class = color, last)
    } else {
      tags$span(class = "text-muted", "En espera")
    }
  })
  
  observeEvent(input$btn_cargar_resultados, {
    withProgress(message = "Cargando resultados...", value = 0.5, {
      cfg             <- rv$config_guardada
      cfg$ruta_salida <- trimws(input$cfg_salida)
      
      if (!dir.exists(cfg$ruta_salida)) {
        showNotification(paste0("❌ Carpeta de salida no encontrada:\n", cfg$ruta_salida), type = "error", duration = 8)
        return()
      }
      
      datos <- tryCatch(
        cargar_todos_resultados(cfg$ruta_salida),
        error = function(e) {
          showNotification(paste0("❌ Error cargando resultados: ", e$message), type = "error", duration = 8)
          list(CNVs = NULL, SVs = NULL)
        }
      )
      rv$resultados <- datos
      
      todas <- bind_rows(datos$CNVs, datos$SVs)
      if (nrow(todas) > 0) {
        fams   <- sort(unique(as.character(todas$ID_Familia)))
        hers   <- sort(unique(her_base(todas$Tipo_Herencia)))
        hers   <- hers[!is.na(hers) & hers != "NA"]
        tipos  <- sort(unique(as.character(todas$SV_type)))
        
        updateSelectizeInput(session, "filtro_sv_type", choices = tipos, server = TRUE)
        updateSelectizeInput(session, "filtro_familia",  choices = fams, server = TRUE)
        updateSelectizeInput(session, "comp_fam_a", choices = fams, server = TRUE)
        updateSelectizeInput(session, "comp_fam_b", choices = fams, server = TRUE)
        updateSelectInput(session, "pdf_familia_sel", choices = c("— Selecciona familia —" = "", setNames(fams, fams)))
        updateSelectizeInput(session, "filtro_herencia", choices = hers, server = TRUE)
        
        scores <- suppressWarnings(as.numeric(todas$AnnotSV_ranking_score))
        scores <- scores[is.finite(scores)]
        if (length(scores) > 0)
          updateSliderInput(session, "filtro_score", min = floor(min(scores) * 10) / 10, max = ceiling(max(scores) * 10) / 10, value = floor(min(scores) * 10) / 10)
      }
      
      incProgress(0.3, message = "Detectando variantes similares...")
      umbral <- isolate(input$cfg_umbral_strict) %||% 0.70
      res_auto <- tryCatch(
        detectar_similares_interindividuales(todas, umbral_sim = umbral, verbose = TRUE),
        error = function(e) list(claves = character(0), n_pares = 0L, n_variantes = 0L)
      )
      auto          <- res_auto$claves
      rv$flags_auto <- auto
      
      ruta_fl <- file.path(APP_LOGS_DIR, ".flags_variantes.rds")
      if (file.exists(ruta_fl)) {
        saved  <- tryCatch(readRDS(ruta_fl), error = function(e) NULL)
        if (!is.null(saved) && is.character(saved)) {
          n_pipes  <- lengths(regmatches(saved, gregexpr("|", saved, fixed = TRUE)))
          sin_modo <- saved[n_pipes == 4L]
          con_modo <- if (length(sin_modo) > 0) paste0(sin_modo, "|FULL") else character(0)
          migradas <- c(saved[n_pipes != 4L], con_modo)
          rv$flags <- union(rv$flags, migradas)
        }
      }
      ruta_n <- file.path(APP_LOGS_DIR, ".notas_variantes.rds")
      if (file.exists(ruta_n)) {
        saved_n <- tryCatch(readRDS(ruta_n), error = function(e) NULL)
        if (!is.null(saved_n) && is.list(saved_n)) rv$notas <- saved_n
      }
      
      if (length(auto) > 0)
        showNotification(
          paste0("🟢 Auto-detection: ", res_auto$n_variantes, " variant(s) involved in ",
                 res_auto$n_pares, " par(es) de similitud entre individuos."),
          type = "message", duration = 8
        )
    })
    showNotification(paste0("Cargados: ", nrow(rv$resultados$CNVs %||% data.frame()), " CNVs, ", nrow(rv$resultados$SVs  %||% data.frame()), " SVs"), type = "message", duration = 4)
  })
  
  datos_filtrados <- reactive({
    datos <- rv$resultados
    if (is.null(datos$CNVs) && is.null(datos$SVs)) return(data.frame())
    
    modalidad <- input$filtro_modalidad
    df <- switch(modalidad,
                 "CNVs"  = datos$CNVs %||% data.frame(),
                 "SVs"   = datos$SVs  %||% data.frame(),
                 "Ambas" = {
                   cnv <- if (!is.null(datos$CNVs)) mutate(datos$CNVs, Modalidad = "CNV") else data.frame()
                   sv  <- if (!is.null(datos$SVs))  mutate(datos$SVs,  Modalidad = "SV")  else data.frame()
                   bind_rows(cnv, sv)
                 }
    )
    if (nrow(df) == 0) return(df)
    
    if (length(input$filtro_familia) > 0) df <- df[df$ID_Familia %in% input$filtro_familia, ]
    sexo_sel <- input$filtro_sexo
    if (!is.null(sexo_sel) && sexo_sel != "todos") {
      todas_fams_df <- unique(as.character(df$ID_Familia))
      fams_conocidas <- names(rv$sexos)[vapply(rv$sexos, function(s) nzchar(s %||% ""), logical(1))]
      if (sexo_sel == "ND") {
        df <- df[df$ID_Familia %in% setdiff(todas_fams_df, fams_conocidas), , drop = FALSE]
      } else {
        fams_sexo <- names(rv$sexos)[vapply(rv$sexos, function(s) identical(s, sexo_sel), logical(1))]
        df <- df[df$ID_Familia %in% fams_sexo, , drop = FALSE]
      }
    }
    if (length(input$filtro_herencia) > 0) df <- df[her_base(df$Tipo_Herencia) %in% input$filtro_herencia, ]
    if (length(input$filtro_sv_type) > 0) df <- df[df$SV_type %in% input$filtro_sv_type, ]
    if (length(input$filtro_rango) > 0) df <- df[df$Tipo_Rango %in% input$filtro_rango, ]
    if (input$filtro_sfari) df <- df[!is.na(df$En_Panel_Genes) & es_sfari_col(df$En_Panel_Genes), ]
    gen_busq <- trimws(input$filtro_gen)
    if (nzchar(gen_busq) && "Gene_name" %in% names(df)) {
      patron <- paste(strsplit(gen_busq, "[,; ]+")[[1]], collapse = "|")
      df <- df[grepl(patron, df$Gene_name, ignore.case = TRUE, perl = TRUE), , drop = FALSE]
    }
    flag_sel <- input$filtro_flagged
    if (!is.null(flag_sel) && flag_sel != "todas") {
      claves_df <- hacer_clave_variante(df)
      if (length(claves_df) == nrow(df)) {
        tiene_Tag <- !is.na(claves_df) & (claves_df %in% rv$flags | claves_df %in% rv$flags_auto)
        if (flag_sel == "con")
          df <- df[tiene_Tag,  , drop = FALSE]
        else if (flag_sel == "sin")
          df <- df[!tiene_Tag, , drop = FALSE]
      }
    }
    clasif_sel <- input$filtro_clasif
    if (!is.null(clasif_sel) && clasif_sel != "todas") {
      claves_df <- hacer_clave_variante(df)
      if (length(claves_df) == nrow(df)) {
        clasif_df <- vapply(claves_df, function(k) {
          if (is.na(k)) "" else rv$clasificaciones[[k]] %||% ""
        }, character(1))
        if (clasif_sel == "sin")
          df <- df[!nzchar(clasif_df), , drop = FALSE]
        else
          df <- df[clasif_df == clasif_sel, , drop = FALSE]
      }
    }
    score_col <- suppressWarnings(as.numeric(df$AnnotSV_ranking_score))
    solo_con_score <- isTRUE(input$filtro_solo_con_score)
    if (solo_con_score) {
      df <- df[!is.na(score_col) & score_col >= input$filtro_score, ]
    } else {
      df <- df[is.na(score_col) | score_col >= input$filtro_score, ]
    }
    
    df
  })
  
  output$ui_tabla_info <- renderUI({
    n <- nrow(datos_filtrados())
    badge_color <- if (n == 0) "bg-secondary" else "bg-primary"
    span(class = paste("badge", badge_color, "fs-6"), n, "variantes")
  })
  
  COLS_DISPONIBLES <- c(
    "ID_Familia","SV_chrom","SV_start","SV_end","SV_length","SV_type",
    "Tipo_Herencia","Genotipo_Hijo","Tipo_Rango","AnnotSV_ranking_score",
    "En_Panel_Genes","Gene_name","Illumina_DRAGEN.exact.counts","Illumina_DRAGEN.similar.counts",
    "Referencia_Match","Sim_MaxLen","Jaccard","Confianza_Region","Annotation_mode","Modalidad"
  )
  COLS_DEFAULT <- c(
    "ID_Familia","SV_chrom","SV_start","SV_end","SV_length","SV_type",
    "Tipo_Herencia","Genotipo_Hijo","Tipo_Rango","AnnotSV_ranking_score",
    "En_Panel_Genes","Gene_name","Illumina_DRAGEN.exact.counts","Illumina_DRAGEN.similar.counts",
    "Referencia_Match","Sim_MaxLen","Jaccard","Confianza_Region"
  )
  
  cols_usuario <- reactiveVal(COLS_DEFAULT)
  
  observeEvent(input$btn_aplicar_columnas, {
    sel <- input$sel_columnas
    if (length(sel) > 0) {
      cols_usuario(sel)
      showNotification("✅ Columnas actualizadas.", type="message", duration=2)
    }
  })
  
  observeEvent(input$btn_reset_columnas, {
    cols_usuario(COLS_DEFAULT)
    updateCheckboxGroupInput(session, "sel_columnas", selected=COLS_DEFAULT)
    showNotification("Columnas restauradas a valores por defecto.", type="message", duration=2)
  })
  
  output$tabla_variantes <- renderDT({
    df <- datos_filtrados()
    if (nrow(df) == 0) return(datatable(data.frame(Mensaje = "Sin variantes con los filtros actuales"), options = list(dom = "t"), rownames = FALSE))
    
    archivo_col <- if ("._archivo_" %in% names(df)) df$`._archivo_` else rep(NA, nrow(df))
    cols_vis <- intersect(cols_usuario(), names(df))
    df_show  <- df[, cols_vis, drop = FALSE]
    df_show$`._archivo_` <- archivo_col
    
    claves_df <- tryCatch(hacer_clave_variante(df), error = function(e) rep(NA_character_, nrow(df)))
    if (length(claves_df) != nrow(df)) claves_df <- rep(NA_character_, nrow(df))
    es_auto   <- !is.na(claves_df) & claves_df %in% rv$flags_auto
    es_manual <- !is.na(claves_df) & claves_df %in% rv$flags
    
    hacer_btn <- function(clave, auto, manual) {
      if (is.na(clave)) return("")
      if (manual && auto) {
        icono <- "🚩"; titulo <- "Manually flagged + similar across individuals"; bg <- "#5cb85c"
      } else if (manual) {
        icono <- "🚩"; titulo <- "Manually flagged — click to remove tag"; bg <- "#5cb85c"
      } else if (auto) {
        icono <- "🟢"; titulo <- "Similar entre individuos (auto) — clic para marcar manualmente"; bg <- "#d4edda"
      } else {
        icono <- "🏳"; titulo <- "Unflagged — click to add tag"; bg <- "transparent"
      }
      sprintf(
        '<button class="flag-btn" data-clave="%s" title="%s" style="background:%s;border:1px solid #ccc;border-radius:4px;cursor:pointer;font-size:1.1em;padding:1px 5px;line-height:1.4;">%s</button>',
        htmltools::htmlEscape(clave, attribute = TRUE), titulo, bg, icono
      )
    }
    df_show$Tag <- mapply(hacer_btn, claves_df, es_auto, es_manual, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    
    hacer_clasif_badge <- function(clave) {
      if (is.na(clave)) return("")
      val <- rv$clasificaciones[[clave]] %||% ""
      cfg <- switch(val,
                    "Benign"     = list(bg="#27AE60", title="Benign"),
                    "VUS"         = list(bg="#95A5A6", title="VUS"),
                    "Pathogenic"  = list(bg="#E74C3C", title="Patog\u00e9nica"),
                    "In Progress" = list(bg="#8E44AD", title="In Progress"),
                    list(bg="transparent", title="Sin clasificar")
      )
      border <- if (nzchar(val)) paste0("2px solid ", cfg$bg) else "2px dashed #CCC"
      sprintf(
        '<span title="%s" style="display:inline-block;width:18px;height:18px;border-radius:5px;background:%s;border:%s;vertical-align:middle;"></span>',
        cfg$title, cfg$bg, border
      )
    }
    
    df_show$Clasif <- vapply(claves_df, hacer_clasif_badge, character(1))
    fam_ids_tbl <- as.character(df$ID_Familia)
    df_show$Sexo <- vapply(fam_ids_tbl, function(fid) {
      s <- rv$sexos[[fid]]
      if (is.null(s) || !nzchar(s)) "—"
      else if (s == "M") "♂" else "♀"
    }, character(1))
    df_show <- df_show[, c("Tag", "Clasif", "Sexo", setdiff(names(df_show), c("Tag", "Clasif", "Sexo"))), drop = FALSE]
    df_show$Herencia_Limpia <- trimws(sub("\\s*\\[.*", "", as.character(df_show$Tipo_Herencia)))
    
    idx_archivo <- which(names(df_show) == "._archivo_") - 1
    idx_limpia  <- which(names(df_show) == "Herencia_Limpia") - 1
    idx_Tag <- which(names(df_show) == "Tag") - 0
    idx_clasif  <- which(names(df_show) == "Clasif") - 0
    idx_sexo    <- which(names(df_show) == "Sexo") - 0
    
    js_cb <- "function(settings) { $('#tabla_variantes').off('click.flagbtn').on('click.flagbtn', '.flag-btn', function(e) { e.stopPropagation(); var clave = $(this).data('clave'); Shiny.setInputValue('flag_btn_click', {clave: clave, t: Date.now()}, {priority: 'event'}); }); }"
    
    # ── LECTURA DEL MODO OSCURO ──
    modo_oscuro <- isTRUE(input$modo_noche == "dark")
    txt_color   <- if (modo_oscuro) "white" else NULL
    
    bg_her <- if (modo_oscuro) c("#CC0000","#0066CC","#CC6600","#6600CC","#4499DD","#DD9944","#9944DD","#6c757d") 
    else c("#FFE6E6","#E6F3FF","#FFF0E6","#F0E6FF","#CCE5FF","#FFE5CC","#EBD7FF","#E8E8E8")
    bg_ran <- if (modo_oscuro) c("#155724","#D35400","#9C0006") else c("#D4EDDA","#FFF3CD","#FFC7CE")
    bg_sco <- if (modo_oscuro) c("#27AE60","#E67E22","#C0392B") else c("#D4EDDA","#FFF3CD","#FFC7CE")
    bg_sfa <- if (modo_oscuro) c("#27AE60","#27AE60","#27AE60","#27AE60","#C0392B","#C0392B","#C0392B","#C0392B","#27AE60","#C0392B") 
    else c("#D4EDDA","#D4EDDA","#D4EDDA","#D4EDDA","#FFC7CE","#FFC7CE","#FFC7CE","#FFC7CE","#D4EDDA","#FFC7CE")
    
    datatable(
      df_show,
      rownames   = FALSE, escape = FALSE, selection = "single", filter = "top", extensions = c("Buttons", "Scroller"),
      options    = list(
        stateSave= TRUE,
        dom = "Bfrtip", buttons = list("colvis"), scrollX = TRUE, scrollY = "530px", scroller = TRUE, deferRender = TRUE,
        columnDefs = list(list(visible = FALSE, targets = c(idx_archivo, idx_limpia)), list(targets = idx_Tag, title = "🏳", orderable = FALSE, width = "40px"), list(targets = idx_clasif, title = "🏷", orderable = FALSE, width = "34px", className = "dt-center"), list(targets = idx_sexo, title = "♂/♀", width = "36px", className = "dt-center")),
        initComplete = JS(js_cb)
      ), class = "table table-hover table-sm"
    ) |>
      formatStyle("Tipo_Herencia", valueColumns = "Herencia_Limpia", 
                  backgroundColor = styleEqual(c("De novo","Paternal","Maternal","Combined","Paternal (Probable)","Maternal (Probable)","Combined (Probable)","Desconocido"), bg_her), 
                  color = txt_color, fontWeight = styleEqual("De novo", "bold")) |>
      formatStyle("Tipo_Rango", 
                  backgroundColor = styleEqual(c("Strict","Wide","Outside"), bg_ran), 
                  color = txt_color, fontWeight = styleEqual("Strict", "bold")) |>
      formatStyle("AnnotSV_ranking_score", 
                  backgroundColor = styleInterval(c(0, 0.5), bg_sco), 
                  color = txt_color, fontWeight = if(modo_oscuro) "bold" else "normal") |>
      formatStyle("En_Panel_Genes", 
                  backgroundColor = styleEqual(c("TRUE","Yes","Yes","Yes","FALSE","No","NO","False",TRUE,FALSE), bg_sfa), 
                  color = txt_color, fontWeight = styleEqual(c("TRUE","Yes","Yes","Yes",TRUE), c("bold","bold","bold","bold","bold"))) |>
      (\(tbl) {
        if ("Genotipo_Hijo" %in% names(df_show))
          formatStyle(tbl, "Genotipo_Hijo",
                      backgroundColor = styleEqual(
                        c("0/1", "1/1", "./1", "0/0"),
                        if (modo_oscuro) c("#1A5276","#922B21","#6C3483","#1E8449")
                        else            c("#D6EAF8","#FADBD8","#E8DAEF","#D5F5E3")
                      ),
                      fontWeight = styleEqual(c("1/1","./1"), c("bold","bold")),
                      color = txt_color)
        else tbl
      })()
  }, server = TRUE)
  
  observeEvent(input$tabla_variantes_rows_selected, {
    idx <- input$tabla_variantes_rows_selected
    df  <- datos_filtrados()
    if (length(idx) == 1) {
      if ("._archivo_" %in% names(df)) rv$fila_seleccionada <- df$`._archivo_`[idx]
      rv$fila_data <- df[idx, , drop = FALSE]
    } else {
      rv$fila_seleccionada <- NULL
      rv$fila_data         <- NULL
    }
  })
  
  output$ui_panel_detalle <- renderUI({
    if (is.null(rv$fila_seleccionada)) return(NULL)
    div(
      style = "resize:both; overflow:auto; min-width:300px; max-width:90vw; min-height:120px; max-height:80vh;",
      card(
        style = "height:100%; width:100%; min-height:0;",
        card_header(
          class = "bg-light d-flex align-items-center gap-2",
          style = "font-size:0.9em;",
          bsicons::bs_icon("info-circle"), " Detalle de variante seleccionada"
        ),
        card_body(
          class = "p-2",
          style = "overflow:auto;",
          uiOutput("ui_boton_abrir_excel"),
          uiOutput("ui_locus_compartido")
        )
      )
    )
  })
  
  output$ui_boton_abrir_excel <- renderUI({
    req(rv$fila_seleccionada)
    archivo <- rv$fila_seleccionada
    franklin_btn <- decipher_btn <- clinvar_btn <- gnomad_btn <- varsome_btn <- hijo_btn <- omim_btn <- ucsc_btn <- ensembl_btn <- hpo_btn <- NULL
    
    if (!is.null(rv$fila_data) && nrow(rv$fila_data) > 0) {
      fd <- rv$fila_data
      chr <- trimws(gsub("(?i)^chr", "", as.character(fd$SV_chrom[1]), perl = TRUE))
      chr_raw <- as.character(fd$SV_chrom[1])
      sta <- as.character(fd$SV_start[1])
      end <- as.character(fd$SV_end[1])
      sty <- toupper(trimws(as.character(fd$SV_type[1])))
      
      modalidad_fila <- if ("Modalidad" %in% names(fd)) as.character(fd$Modalidad[1]) else "SV"
      ruta_tsv_hijo  <- tryCatch(encontrar_archivo_hijo(isolate(input$cfg_entrada), archivo, as.character(fd$ID_Familia[1]), modalidad_fila), error = function(e) NA_character_)
      if (!is.na(ruta_tsv_hijo) && file.exists(ruta_tsv_hijo)) {
        hijo_btn <- actionButton("btn_abrir_hijo", "👫 Open proband TSV", class = "btn-outline-primary btn-sm", style = "font-weight:600;")
      }
      
      franklin_btn <- tags$button(class="btn btn-info btn-sm", style="font-weight:600; letter-spacing:0.03em;", type="button", onclick=paste0("Shiny.setInputValue('open_url', {url:'https://franklin.genoox.com/clinical-db/variant/sv/chr", chr, "-", sta, "-", end, "-", sty, "-hg38', t:Date.now()}, {priority:'event'});"), "🧬 FRANKLIN")
      clinvar_btn  <- tags$button(class="btn btn-danger btn-sm", style="font-weight:600; letter-spacing:0.03em;", type="button", onclick=paste0("Shiny.setInputValue('open_url', {url:'https://www.ncbi.nlm.nih.gov/clinvar/?term=Chr", chr_raw, ":", sta, "-", end, "', t:Date.now()}, {priority:'event'});"), "🧪 ClinVar")
      gnomad_btn   <- tags$button(class="btn btn-warning btn-sm", style="font-weight:600; letter-spacing:0.03em;", type="button", onclick=paste0("Shiny.setInputValue('open_url', {url:'https://gnomad.broadinstitute.org/region/", chr_raw, ":", sta, "-", end, "?dataset=gnomad_sv_r4', t:Date.now()}, {priority:'event'});"), "📊 gnomAD")
      varsome_btn  <- tags$button(class="btn btn-secondary btn-sm", style="font-weight:600; letter-spacing:0.03em;", type="button", onclick=paste0("Shiny.setInputValue('open_url', {url:'https://varsome.com/variant/hg38/Chr", chr_raw, ":", sta, "-", end, ":", sty, "?annotation-mode=germline', t:Date.now()}, {priority:'event'});"), "✅ VarSome")
      loc <- paste0(chr_raw, ":", sta, "-", end)
      decipher_btn <- tags$button(class="btn btn-success btn-sm", style="font-weight:600; letter-spacing:0.03em;", type="button", onclick=paste0("Shiny.setInputValue('open_url', {url:'https://www.deciphergenomics.org/browser#q/", loc, "/location/", loc, "', t:Date.now()}, {priority:'event'});"), "🔍 DECIPHER")
      ucsc_btn     <- tags$button(class="btn btn-sm", style="font-weight:600; letter-spacing:0.03em; background:#336699; color:#fff; border-color:#336699;", type="button", onclick=paste0("Shiny.setInputValue('open_url', {url:'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr", chr_raw, ":", sta, "-", end, "', t:Date.now()}, {priority:'event'});"), "🧭 UCSC")
      ensembl_btn  <- tags$button(class="btn btn-sm", style="font-weight:600; letter-spacing:0.03em; background:#8B3A8B; color:#fff; border-color:#8B3A8B;", type="button", onclick=paste0("Shiny.setInputValue('open_url', {url:'https://www.ensembl.org/Homo_sapiens/Location/View?r=", chr_raw, ":", sta, "-", end, "', t:Date.now()}, {priority:'event'});"), "🔬 Ensembl")
      
      # ── OMIM button(s) ──────────────────────────────────────────────────────
      omim_btn <- NULL
      omim_col <- names(fd)[grep("^OMIM_ID$|omim.*id|mim.*number|omim_morbid|phenotype.*omim|omim",
                                 names(fd), ignore.case = TRUE, perl = TRUE)][1]
      if (!is.na(omim_col)) {
        omim_raw <- as.character(fd[[omim_col]][1])
        omim_ids <- unique(regmatches(omim_raw, gregexpr("[0-9]{6,7}", omim_raw))[[1]])
        if (length(omim_ids) == 1) {
          # Single ID → direct button
          omim_btn <- tags$button(
            class   = "btn btn-sm",
            style   = "font-weight:600; letter-spacing:0.03em; background:#7B2D8B; color:#fff; border-color:#7B2D8B;",
            type    = "button",
            onclick = paste0("Shiny.setInputValue('open_url', {url:'https://omim.org/entry/", omim_ids[1], "', t:Date.now()}, {priority:'event'});"),
            "🔬 OMIM"
          )
        } else if (length(omim_ids) > 1) {
          # Multiple IDs → Bootstrap dropdown
          omim_btn <- tags$div(
            class = "btn-group",
            tags$button(
              class          = "btn btn-sm dropdown-toggle",
              style          = "font-weight:600; letter-spacing:0.03em; background:#7B2D8B; color:#fff; border-color:#7B2D8B;",
              type           = "button",
              onclick        = "event.stopPropagation();",
              `data-bs-toggle` = "dropdown",
              `aria-expanded`  = "false",
              paste0("🔬 OMIM (", length(omim_ids), ")")
            ),
            tags$ul(
              class = "dropdown-menu dropdown-menu-end shadow-sm",
              style = "min-width:160px;",
              lapply(omim_ids, function(id) {
                tags$li(
                  tags$button(
                    class   = "dropdown-item d-flex align-items-center gap-2",
                    style   = "font-size:0.9em; background:none; border:none; width:100%; text-align:left; cursor:pointer;",
                    type    = "button",
                    onclick = paste0("Shiny.setInputValue('open_url', {url:'https://omim.org/entry/", id, "', t:Date.now()}, {priority:'event'});"),
                    tags$span(style = "color:#7B2D8B; font-weight:700;", paste0("# ", id)),
                    tags$span(class = "text-muted", "→ omim.org")
                  )
                )
              })
            )
          )
        }
      }
      
      # ── HPO button(s) ────────────────────────────────────────────────────────
      hpo_btn <- NULL
      # Robust detection: looks for any column whose name contains "gene"
      # para cubrir variantes como Gene_name, gene_name, Gene_Name, genes, etc.
      gene_col_hpo <- names(fd)[grep("gene", names(fd), ignore.case = TRUE, perl = TRUE)][1]
      
      if (!is.na(gene_col_hpo) && !is.null(gene_col_hpo)) {
        gene_raw_hpo <- trimws(as.character(fd[[gene_col_hpo]][1]))
        # Separadores usados por AnnotSV: ; / , | y espacios
        genes_hpo <- unique(trimws(unlist(strsplit(gene_raw_hpo, "[;,/|[:space:]]+"))))
        genes_hpo <- genes_hpo[nzchar(genes_hpo) & genes_hpo != "." & 
                                 toupper(genes_hpo) != "NA" & genes_hpo != "-"]
        
        hpo_color <- "#007B7B"   # teal oscuro
        
        if (length(genes_hpo) == 1) {
          # Single gene -> direct button
          hpo_btn <- tags$button(
            class   = "btn btn-sm",
            style   = paste0("font-weight:600; letter-spacing:0.03em; background:", hpo_color,
                             "; color:#fff; border-color:", hpo_color, ";"),
            type    = "button",
            onclick = paste0("Shiny.setInputValue('open_url', {url:'https://hpo.jax.org/browse/search?q=",
                             URLencode(genes_hpo[1], reserved = TRUE),
                             "&navFilter=all', t:Date.now()}, {priority:'event'});"),
            "📍 HPO"
          )
        } else if (length(genes_hpo) > 1) {
          # Varios genes -> dropdown igual que OMIM
          hpo_btn <- tags$div(
            class = "btn-group",
            tags$button(
              class            = "btn btn-sm dropdown-toggle",
              style            = paste0("font-weight:600; letter-spacing:0.03em; background:", hpo_color,
                                        "; color:#fff; border-color:", hpo_color, ";"),
              type             = "button",
              onclick          = "event.stopPropagation();",
              `data-bs-toggle` = "dropdown",
              `aria-expanded`  = "false",
              paste0("🧬 HPO (", length(genes_hpo), ")")
            ),
            tags$ul(
              class = "dropdown-menu dropdown-menu-end shadow-sm",
              style = "min-width:180px;",
              lapply(genes_hpo, function(g) {
                tags$li(
                  tags$button(
                    class   = "dropdown-item d-flex align-items-center gap-2",
                    style   = "font-size:0.9em; background:none; border:none; width:100%; text-align:left; cursor:pointer;",
                    type    = "button",
                    onclick = paste0("Shiny.setInputValue('open_url', {url:'https://hpo.jax.org/browse/search?q=",
                                     URLencode(g, reserved = TRUE),
                                     "&navFilter=all', t:Date.now()}, {priority:'event'});"),
                    tags$span(style = paste0("color:", hpo_color, "; font-weight:700;"), g),
                    tags$span(class = "text-muted", "→ hpo.jax.org")
                  )
                )
              })
            )
          )
        }
      }
    }  # fin if (!is.null(rv$fila_data))
    
    # --- Tag information block for the selected row ---
    fd_flag <- rv$fila_data
    info_Tag <- NULL
    if (!is.null(fd_flag) && nrow(fd_flag) > 0) {
      clave_sel  <- tryCatch(hacer_clave_variante(fd_flag)[1], error = function(e) NA_character_)
      es_manual  <- !is.na(clave_sel) && clave_sel %in% rv$flags
      es_auto_f  <- !is.na(clave_sel) && clave_sel %in% rv$flags_auto
      if (es_manual || es_auto_f) {
        todas_df <- tryCatch({
          bind_rows(
            if (!is.null(rv$resultados$CNVs)) mutate(rv$resultados$CNVs, Modalidad = "CNV") else data.frame(),
            if (!is.null(rv$resultados$SVs))  mutate(rv$resultados$SVs,  Modalidad = "SV")  else data.frame()
          )
        }, error = function(e) data.frame())
        umbral_f <- isolate(input$cfg_umbral_strict) %||% 0.70
        sim_claves <- if (!is.na(clave_sel) && nrow(todas_df) > 0)
          tryCatch(detectar_similares_a_clave(clave_sel, todas_df, umbral_sim = umbral_f),
                   error = function(e) character(0))
        else character(0)
        fams_sim <- if (length(sim_claves) > 0)
          unique(sapply(strsplit(sim_claves, "|", fixed = TRUE), `[`, 1))
        else character(0)
        tipo_tag <- if (es_manual && es_auto_f) "🚩🟢 Manually flagged and similar across individuals"
        else if (es_manual) "🚩 Manually flagged"
        else "🟢 Similar entre individuos (auto)"
        sim_texto <- if (length(fams_sim) > 0)
          tags$span(class = "ms-2", tags$b(paste0(length(sim_claves), " variante(s) similar(es) en: ")), paste(fams_sim, collapse = ", "))
        else tags$span(class = "ms-2 text-muted fst-italic", "sin coincidencias con otros individuos")
        tag_extra_class <- if (es_manual) "dm-tag-manual" else "dm-tag-auto"
        info_Tag <- div(
          class = paste("mt-2 px-2 py-1 rounded small d-flex align-items-start flex-wrap gap-1", tag_extra_class),
          tags$span(class = "fw-bold", tipo_tag), sim_texto)
      } else {
        info_Tag <- div(class = "mt-2 small text-muted",
                        "🏳 unflagged  |  🟢 similar across individuals (auto)  |  🚩 manually flagged")
      }
    } else {
      info_Tag <- div(class = "mt-2 small text-muted",
                      "🏳 unflagged  |  🟢 similar across individuals (auto)  |  🚩 manually flagged")
    }
    
    # --- Clinical classification of selected variant (dropdown) ---
    fd_clasif <- rv$fila_data
    info_clasif <- NULL
    if (!is.null(fd_clasif) && nrow(fd_clasif) > 0) {
      clave_clasif <- tryCatch(hacer_clave_variante(fd_clasif)[1], error = function(e) NA_character_)
      val_clasif   <- if (!is.na(clave_clasif)) rv$clasificaciones[[clave_clasif]] %||% "" else ""
      CLASIF_BG <- c("Benign"="#27AE60","VUS"="#95A5A6","Pathogenic"="#E74C3C","In Progress"="#8E44AD")
      badge_style <- if (nzchar(val_clasif))
        paste0("display:inline-block;width:14px;height:14px;border-radius:4px;",
               "background:", CLASIF_BG[val_clasif], ";vertical-align:middle;margin-right:6px;")
      else
        "display:inline-block;width:14px;height:14px;border-radius:4px;border:2px dashed #CCC;vertical-align:middle;margin-right:6px;"
      info_clasif <- div(
        class = "mt-2 px-2 py-2 rounded",
        class = "dm-clasif-panel rounded px-2 py-1",
        div(class = "d-flex align-items-center gap-2 mb-1",
            tags$span(style = badge_style),
            tags$b(class = "dm-dark-txt", style = "font-size:0.85em;", "Clinical classification"),
            if (nzchar(val_clasif))
              actionButton("btn_quitar_clasif", "× Remove",
                           class = "btn-outline-danger btn-sm ms-auto",
                           style = "font-size:0.78em; padding:1px 7px;")
        ),
        selectInput(
          inputId  = "clasif_detalle_select",
          label    = NULL,
          choices  = c(
            "— Sin clasificar" = "",
            "🟢 Benigna"    = "Benign",
            "⬜ VUS"            = "VUS",
            "🔴 Pathogenic" = "Pathogenic",
            "🟣 En Progreso" = "In Progress"
          ),
          selected = val_clasif,
          width    = "100%"
        )
      )
    }
    # --- Nota de la variante seleccionada ---
    clave_nota <- if (!is.null(rv$fila_data) && nrow(rv$fila_data) > 0)
      tryCatch(hacer_clave_variante(rv$fila_data)[1], error = function(e) NA_character_)
    else NA_character_
    
    # --- Nota de la variante seleccionada ---
    clave_nota <- if (!is.null(rv$fila_data) && nrow(rv$fila_data) > 0)
      tryCatch(hacer_clave_variante(rv$fila_data)[1], error = function(e) NA_character_)
    else NA_character_
    
    # MODIFICADO: Recuperar el texto guardado o el fantasma
    nota_guardada <- if (!is.na(clave_nota) && !is.null(rv$notas[[clave_nota]])) rv$notas[[clave_nota]] else ""
    nota_actual <- isolate(rv$nota_fantasma) %||% nota_guardada
    
    panel_notas <- if (!is.na(clave_nota)) {
      tiene_nota    <- nzchar(trimws(nota_actual))
      preview_nota  <- if (tiene_nota) {
        txt <- trimws(nota_actual)
        if (nchar(txt) > 55) paste0(substr(txt, 1, 55), "\u2026") else txt
      } else NULL
      div(
        class = "mt-2",
        tags$details(
          # Siempre cerrado al renderizar
          tags$summary(
            style = "cursor:pointer; font-weight:600; font-size:0.9em; color:#2C6FAC; user-select:none; list-style:none; display:flex; align-items:center; gap:6px;",
            tags$span("\U0001f4dd Anotaciones de esta variante"),
            # Recuadrito azul con preview
            if (tiene_nota)
              tags$span(
                class = "dm-preview-nt",
                style = "display:inline-block; max-width:220px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; font-size:0.75em; font-weight:400; border-radius:4px; padding:1px 7px; vertical-align:middle; border:1px solid;",
                preview_nota
              )
            else
              tags$span(
                class = "dm-preview-nt-e",
                style = "font-size:0.75em; font-weight:400; font-style:italic; border-radius:4px; padding:1px 7px; vertical-align:middle; border:1px dashed;",
                "sin anotar"
              )
          ),
          div(
            class = "mt-2",
            tags$textarea(
              id    = "input_nota_variante",
              class = "form-control form-control-sm",
              style = "min-height:80px; font-size:0.85em; resize:vertical;",
              placeholder = "Type your notes about this variant here...",
              nota_actual
            ),
            div(
              class = "d-flex gap-2 mt-1 align-items-center",
              actionButton("btn_guardar_nota", "\U0001f4be Guardar nota", class = "btn-primary btn-sm"),
              actionButton("btn_nota_replicar", "\U0001f501 Guardar y replicar a similares", class = "btn-outline-secondary btn-sm"),
              if (tiene_nota)
                actionButton("btn_borrar_nota", "\U0001f5d1 Borrar nota", class = "btn-outline-danger btn-sm")
            )
          )
        )
      )
    } else NULL
    
    fam_sel_det <- if (!is.null(rv$fila_data) && nrow(rv$fila_data) > 0)
      as.character(rv$fila_data$ID_Familia[1]) else NULL
    sexo_actual_det <- if (!is.null(fam_sel_det)) rv$sexos[[fam_sel_det]] %||% "" else ""
    btn_varon <- if (!is.null(fam_sel_det)) {
      active_v <- identical(sexo_actual_det, "M")
      actionButton("btn_sexo_varon",
                   label = if (active_v) "♂ Male ✓" else "♂ Male",
                   class = if (active_v) "btn-primary btn-sm" else "btn-outline-primary btn-sm",
                   title = if (active_v) "Marked — click to unmark" else "Mark as Male")
    } else NULL
    btn_mujer <- if (!is.null(fam_sel_det)) {
      active_f <- identical(sexo_actual_det, "F")
      actionButton("btn_sexo_mujer",
                   label = if (active_f) "♀ Female ✓" else "♀ Female",
                   class = if (active_f) "btn-danger btn-sm" else "btn-outline-danger btn-sm",
                   title = if (active_f) "Marked — click to unmark" else "Marcar como Mujer")
    } else NULL
    tagList(
      div(class = "d-flex align-items-center flex-wrap gap-2 mb-1",
          div(class = "d-flex align-items-center gap-2",
              tags$span(class = "text-muted small", "Variante en: ", tags$b(if (!is.na(archivo)) basename(archivo) else "—")),
              if (!is.na(archivo)) actionButton("btn_abrir_excel_sel", "📂 Open Excel", class = "btn-warning btn-sm",
                                                title = if (file.exists(archivo)) archivo else paste0("⚠ File not found: ", basename(archivo))),
              hijo_btn),
          tags$span(class = "ms-auto"),
          div(class = "d-flex align-items-center gap-2 flex-wrap", franklin_btn, decipher_btn, clinvar_btn, gnomad_btn, varsome_btn, ucsc_btn, ensembl_btn, omim_btn, hpo_btn)
      ),
      div(class = "d-flex align-items-center gap-2 px-2 py-1 rounded mb-1",
          class = "dm-clasif-panel rounded px-2 py-1",
          tags$span(class = "text-muted small",
                    tags$b("Sexo:"),
                    if (!is.null(fam_sel_det)) tags$span(class = "ms-1 text-secondary", paste0("(", fam_sel_det, ")"))
          ),
          btn_varon, btn_mujer,
          if (!is.null(fam_sel_det) && nzchar(sexo_actual_det))
            actionButton("btn_sexo_borrar", "× Unmark",
                         class = "btn-outline-secondary btn-sm",
                         title = "Remove sex information for this family")
      ),
      # ── Family clinical phenotype ─────────────────────────────────────
      local({
        fenotipo_guardado <- if (!is.null(fam_sel_det)) rv$fenotipos[[fam_sel_det]] %||% "" else ""
        fenotipo_actual <- isolate(rv$feno_fantasma) %||% fenotipo_guardado
        tiene_fenotipo  <- nzchar(trimws(fenotipo_actual))
        
        preview_ft <- if (tiene_fenotipo) {
          txt <- trimws(fenotipo_actual)
          if (nchar(txt) > 55) paste0(substr(txt, 1, 55), "…") else txt
        } else NULL
        
        div(
          class = "mt-1 mb-1",
          tags$details(
            tags$summary(
              style = "cursor:pointer; font-weight:600; font-size:0.9em; color:#6C3483; user-select:none; list-style:none; display:flex; align-items:center; gap:6px;",
              tags$span("🧬 Clinical phenotype"),
              if (!is.null(fam_sel_det)) tags$span(class = "text-muted fw-normal", style = "font-size:0.82em;", paste0("(", fam_sel_det, ")")),
              if (tiene_fenotipo)
                tags$span(class = "dm-preview-ft", style = "display:inline-block; max-width:220px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; font-size:0.75em; font-weight:400; border-radius:4px; padding:1px 7px; vertical-align:middle; border:1px solid;", preview_ft)
              else
                tags$span(class = "dm-preview-ft-e", style = "font-size:0.75em; font-weight:400; font-style:italic; border-radius:4px; padding:1px 7px; vertical-align:middle; border:1px dashed;", "sin registrar")
            ),
            div(
              class = "dm-fenotipo-panel rounded mt-2 px-2 py-2",
              if (!is.null(fam_sel_det)) {
                tagList(
                  tags$textarea(
                    id          = "input_fenotipo_familia",
                    class       = "form-control form-control-sm dm-fenotipo-input", 
                    style       = "min-height:80px; font-size:0.85em; resize:vertical;",
                    placeholder = paste0("Describe el fenotipo del probando de la familia ", fam_sel_det, "...\nEj: TEA nivel 2, discapacidad intelectual leve, epilepsia, dismorfias faciales"),
                    fenotipo_actual
                  ),
                  div(
                    class = "d-flex gap-2 mt-1 align-items-center",
                    actionButton("btn_guardar_fenotipo", "💾 Guardar fenotipo", class = "btn-sm", style = "background:#6C3483; color:white; border-color:#6C3483;"),
                    if (tiene_fenotipo) actionButton("btn_borrar_fenotipo", "🗑 Delete", class = "btn-outline-danger btn-sm")
                  )
                )
              } else {
                tags$span(class = "text-muted small fst-italic", "Selecciona una variante para registrar el fenotipo de su familia.")
              }
            )
          )
        )
      }),
      
      # ── Info de Tags (auto / manual) ──────────────────────────────────────
      info_Tag,
      
      # ── Clinical classification (dropdown) ───────────────────────────────
      info_clasif,
      
      # ── Panel de Notas (Aseguramos que se muestre) ──
      panel_notas
      
    ) # <-- Cierre de la etiqueta tagList principal
  }) # <-- Cierre de output$ui_boton_abrir_excel    
  
  # ── Guardar nota ──────────────────────────────────────────────────────────
  observeEvent(input$btn_guardar_nota, {
    fd <- rv$fila_data
    req(!is.null(fd), nrow(fd) > 0)
    clave <- tryCatch(hacer_clave_variante(fd)[1], error = function(e) NA_character_)
    req(!is.na(clave))
    texto <- isolate(input$input_nota_variante)
    rv$notas[[clave]] <- texto
    rv$historial <- c(rv$historial, list(list(
      ts = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tipo = "📝 Note saved", clave = clave,
      detalle = paste0(nchar(trimws(texto)), " caracteres")
    )))
    showNotification("💾 Nota guardada.", type = "message", duration = 2)
  })
  
  # ── Guardar y replicar nota a variantes similares ─────────────────────────
  observeEvent(input$btn_nota_replicar, {
    fd <- rv$fila_data
    req(!is.null(fd), nrow(fd) > 0)
    clave <- tryCatch(hacer_clave_variante(fd)[1], error = function(e) NA_character_)
    req(!is.na(clave))
    texto <- isolate(input$input_nota_variante)
    rv$notas[[clave]] <- texto
    
    todas_df <- tryCatch({
      bind_rows(
        if (!is.null(rv$resultados$CNVs)) mutate(rv$resultados$CNVs, Modalidad = "CNV") else data.frame(),
        if (!is.null(rv$resultados$SVs))  mutate(rv$resultados$SVs,  Modalidad = "SV")  else data.frame()
      )
    }, error = function(e) data.frame())
    umbral <- isolate(input$cfg_umbral_strict) %||% 0.70
    similares <- tryCatch(detectar_similares_a_clave(clave, todas_df, umbral_sim = umbral),
                          error = function(e) character(0))
    n_rep <- 0L
    for (k in similares) {
      rv$notas[[k]] <- texto
      n_rep <- n_rep + 1L
    }
    rv$historial <- c(rv$historial, list(list(
      ts = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tipo = "📝 Note replicated", clave = clave,
      detalle = paste0("replicada a ", n_rep, " variante(s)")
    )))
    msg <- if (n_rep > 0)
      paste0("💾 Nota guardada y replicada a ", n_rep, " variante(s) similar(es).")
    else
      "💾 Nota guardada. No se encontraron variantes similares a las que replicar."
    showNotification(msg, type = "message", duration = 4)
  })
  
  # ── Borrar nota ───────────────────────────────────────────────────────────
  observeEvent(input$btn_borrar_nota, {
    fd <- rv$fila_data
    req(!is.null(fd), nrow(fd) > 0)
    clave <- tryCatch(hacer_clave_variante(fd)[1], error = function(e) NA_character_)
    req(!is.na(clave))
    rv$notas[[clave]] <- NULL
    rv$historial <- c(rv$historial, list(list(
      ts = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tipo = "🗑 Note deleted", clave = clave, detalle = ""
    )))
    showNotification("🗑 Note deleted.", type = "warning", duration = 2)
  })
  
  # ── Fenotipo familiar ─────────────────────────────────────────────────────
  observeEvent(input$btn_guardar_fenotipo, {
    fd <- rv$fila_data
    req(!is.null(fd), nrow(fd) > 0)
    fam   <- as.character(fd$ID_Familia[1])
    texto <- trimws(input$input_fenotipo_familia %||% "")
    rv$fenotipos[[fam]] <- texto
    rv$historial <- c(rv$historial, list(list(
      ts      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      tipo    = "🧬 Fenotipo guardado",
      clave   = fam,
      detalle = if (nchar(texto) > 60) paste0(substr(texto, 1, 60), "…") else texto
    )))
    if (nzchar(texto))
      showNotification(paste0("🧬 Fenotipo guardado para ", fam, "."),
                       type = "message", duration = 2)
    else
      showNotification(paste0("Fenotipo de ", fam, " borrado."),
                       type = "warning", duration = 2)
  })
  
  observeEvent(input$btn_borrar_fenotipo, {
    fd <- rv$fila_data
    req(!is.null(fd), nrow(fd) > 0)
    fam <- as.character(fd$ID_Familia[1])
    rv$fenotipos[[fam]] <- NULL
    rv$historial <- c(rv$historial, list(list(
      ts = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      tipo = "🗑 Phenotype deleted", clave = fam, detalle = ""
    )))
    showNotification(paste0("🗑 Phenotype of ", fam, " deleted."),
                     type = "warning", duration = 2)
  })
  
  # ── Resetear textareas al cambiar de variante / familia ─────────────────────
  # Shiny cachea los valores de inputs por ID; sin este observer, al seleccionar
  # a different row the fields would inherit values from the previous one.
  observeEvent(rv$fila_data, {
    # 1. Fenotipo familiar — se indexa por ID_Familia
    fam <- tryCatch(
      if (!is.null(rv$fila_data) && nrow(rv$fila_data) > 0)
        as.character(rv$fila_data$ID_Familia[1]) else NULL,
      error = function(e) NULL
    )
    val_ft <- if (!is.null(fam)) rv$fenotipos[[fam]] %||% "" else ""
    updateTextAreaInput(session, "input_fenotipo_familia", value = val_ft)
    
    # 2. Variant note — indexed by unique variant key
    clave <- tryCatch(
      if (!is.null(rv$fila_data) && nrow(rv$fila_data) > 0)
        hacer_clave_variante(rv$fila_data)[1] else NA_character_,
      error = function(e) NA_character_
    )
    val_nota <- if (!is.na(clave) && !is.null(rv$notas[[clave]])) rv$notas[[clave]] else ""
    updateTextAreaInput(session, "input_nota_variante", value = val_nota)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  
  # ── Abrir URLs externas con browseURL (evita el visor de RStudio) ────────────
  observeEvent(input$open_url, {
    url <- input$open_url$url
    if (!is.null(url) && nzchar(url)) browseURL(url)
  })
  
  observeEvent(input$btn_abrir_excel_sel, {
    req(rv$fila_seleccionada)
    if (!abrir_archivo_nativo(rv$fila_seleccionada)) showNotification("No se pudo abrir el archivo.", type = "error")
  })
  
  # ---- Sexo del probando ------------------------------------------------
  observeEvent(input$btn_sexo_varon, {
    fd <- rv$fila_data; req(!is.null(fd), nrow(fd) > 0)
    fam <- as.character(fd$ID_Familia[1])
    if (identical(rv$sexos[[fam]], "M")) { rv$sexos[[fam]] <- NULL
    showNotification(paste0("Sexo de ", fam, " desmarcado."), type="message", duration=2)
    } else { rv$sexos[[fam]] <- "M"
    showNotification(paste0("♂ ", fam, " marked as Male."), type="message", duration=2) }
  })
  observeEvent(input$btn_sexo_mujer, {
    fd <- rv$fila_data; req(!is.null(fd), nrow(fd) > 0)
    fam <- as.character(fd$ID_Familia[1])
    if (identical(rv$sexos[[fam]], "F")) { rv$sexos[[fam]] <- NULL
    showNotification(paste0("Sexo de ", fam, " desmarcado."), type="message", duration=2)
    } else { rv$sexos[[fam]] <- "F"
    showNotification(paste0("♀ ", fam, " marcado como Mujer."), type="message", duration=2) }
  })
  observeEvent(input$btn_sexo_borrar, {
    fd <- rv$fila_data; req(!is.null(fd), nrow(fd) > 0)
    fam <- as.character(fd$ID_Familia[1])
    rv$sexos[[fam]] <- NULL
    showNotification(paste0("Sexo de ", fam, " deleted."), type="warning", duration=2)
  })
  
  observeEvent(input$btn_abrir_hijo, {
    req(rv$fila_data)
    fd <- rv$fila_data
    modalidad_fila <- if ("Modalidad" %in% names(fd)) as.character(fd$Modalidad[1]) else "SV"
    ruta_tsv <- tryCatch(encontrar_archivo_hijo(input$cfg_entrada, rv$fila_seleccionada, as.character(fd$ID_Familia[1]), modalidad_fila), error = function(e) NA_character_)
    if (is.na(ruta_tsv) || !file.exists(ruta_tsv)) {
      showNotification(paste0("❌ Proband TSV not found for family ", fd$ID_Familia[1], "."), type = "error", duration = 6)
    } else {
      if (!abrir_archivo_nativo(ruta_tsv)) showNotification("No se pudo abrir el archivo TSV.", type = "error")
    }
  })
  
  # ── Sistema de Tags ────────────────────────────────────────────────────────────
  observeEvent(input$flag_btn_click, {
    payload <- input$flag_btn_click
    req(!is.null(payload), nzchar(payload$clave))
    clave <- payload$clave
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    if (clave %in% rv$flags) {
      rv$flags <- setdiff(rv$flags, clave)
      rv$historial <- c(rv$historial, list(list(ts=ts, tipo="🏳 Tag removed", clave=clave)))
      showNotification("🏳 Tag removed.", type = "warning", duration = 2)
    } else {
      rv$flags <- union(rv$flags, clave)
      todas_df <- tryCatch({
        bind_rows(
          if (!is.null(rv$resultados$CNVs)) mutate(rv$resultados$CNVs, Modalidad = "CNV") else data.frame(),
          if (!is.null(rv$resultados$SVs))  mutate(rv$resultados$SVs,  Modalidad = "SV")  else data.frame()
        )
      }, error = function(e) data.frame())
      umbral <- isolate(input$cfg_umbral_strict) %||% 0.70
      similares <- tryCatch(detectar_similares_a_clave(clave, todas_df, umbral_sim = umbral),
                            error = function(e) character(0))
      # Auto-tag solo si la variante tiene MAS de 2 similares interindividuales
      n_sim <- length(similares)
      if (n_sim > 2) {
        nuevas_auto <- setdiff(similares, c(clave, rv$flags_auto))
        if (length(nuevas_auto) > 0) rv$flags_auto <- union(rv$flags_auto, nuevas_auto)
      }
      rv$historial <- c(rv$historial, list(list(
        ts=ts, tipo="🚩 Tag added", clave=clave,
        detalle=if (n_sim > 0) paste0(n_sim, " similar(es)", if (n_sim > 2) " — auto-marcadas" else " (umbral no alcanzado)") else "sin similares")))
      if (n_sim > 2)
        showNotification(paste0("🚩 Variant flagged. Found ", n_sim,
                                " similar variant(s) in other individuals, automatically flagged."),
                         type = "message", duration = 5)
      else if (n_sim > 0)
        showNotification(paste0("🚩 Variant flagged. Found ", n_sim,
                                " similar(es), pero se necesitan m\u00e1s de 2 para el auto-etiquetado."),
                         type = "message", duration = 4)
      else
        showNotification("🚩 Variant flagged. No similar variants found in other individuals.",
                         type = "message", duration = 3)
    }
  })
  
  
  # ── Clinical classification: saved only when a real value is selected ──
  observeEvent(input$clasif_detalle_select, {
    valor <- input$clasif_detalle_select
    req(nzchar(valor))
    fd <- rv$fila_data
    req(!is.null(fd), nrow(fd) > 0)
    clave <- tryCatch(hacer_clave_variante(fd)[1], error = function(e) NA_character_)
    req(!is.na(clave))
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    CLASIF_LABEL <- c(
      "Benign"="\U0001f7e2 Benigna", "VUS"="\u2b1c VUS",
      "Pathogenic"="\U0001f534 Patog\u00e9nica", "In Progress"="\U0001f7e3 En Progreso"
    )
    rv$clasificaciones[[clave]] <- valor
    label <- CLASIF_LABEL[valor] %||% valor
    rv$historial <- c(rv$historial, list(list(
      ts=ts, tipo=paste0(label, " \u2014 Clasificada"), clave=clave, detalle=valor)))
    showNotification(paste0(label, " guardada."), type="message", duration=2)
  }, ignoreInit = TRUE)
  
  # Explicit button to remove classification ──────────────────────────
  observeEvent(input$btn_quitar_clasif, {
    fd <- rv$fila_data
    req(!is.null(fd), nrow(fd) > 0)
    clave <- tryCatch(hacer_clave_variante(fd)[1], error = function(e) NA_character_)
    req(!is.na(clave))
    rv$clasificaciones[[clave]] <- NULL
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    rv$historial <- c(rv$historial, list(list(
      ts=ts, tipo="\u2715 Clasificaci\u00f3n eliminada", clave=clave, detalle="")))
    showNotification("Clasificaci\u00f3n eliminada.", type="warning", duration=2)
  })
  
  output$btn_descargar_xlsx <- downloadHandler(
    filename = function() paste0("variantes_filtradas_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      df <- datos_filtrados()
      df <- df[, !names(df) %in% "._archivo_", drop = FALSE]
      wb <- createWorkbook()
      addWorksheet(wb, "Variantes")
      writeData(wb, "Variantes", df)
      col_her <- which(names(df) == "Tipo_Herencia")
      col_ran <- which(names(df) == "Tipo_Rango")
      col_sfa <- which(names(df) == "En_Panel_Genes")
      
      if (length(col_sfa) > 0) {
        estilos_sfa <- list("Yes" = createStyle(fgFill="#D4EDDA", fontColour="#155724", textDecoration="bold"), "No" = createStyle(fgFill="#FFC7CE", fontColour="#9C0006", textDecoration="bold"))
        her_vals <- her_base(df$En_Panel_Genes)
        for (nm in names(estilos_sfa)) {
          rows <- which(her_vals == nm)
          if (length(rows) > 0) addStyle(wb, "Variantes", estilos_sfa[[nm]], rows = rows + 1, cols = col_sfa)
        }
      }
      
      if (length(col_her) > 0) {
        estilos_her <- list("De novo" = createStyle(fgFill="#FFE6E6", fontColour="#CC0000", textDecoration="bold"), "Paternal" = createStyle(fgFill="#E6F3FF", fontColour="#0066CC"), "Maternal" = createStyle(fgFill="#FFF0E6", fontColour="#CC6600"), "Combined" = createStyle(fgFill="#F0E6FF", fontColour="#6600CC"), "Paternal (Probable)" = createStyle(fgFill="#CCE5FF", fontColour="#0066CC"), "Maternal (Probable)" = createStyle(fgFill="#FFE5CC", fontColour="#CC6600"), "Combined (Probable)"= createStyle(fgFill="#EBD7FF", fontColour="#6600CC"))
        her_vals <- her_base(df$Tipo_Herencia)
        for (nm in names(estilos_her)) {
          rows <- which(her_vals == nm)
          if (length(rows) > 0) addStyle(wb, "Variantes", estilos_her[[nm]], rows = rows + 1, cols = col_her)
        }
      }
      if (length(col_ran) > 0) {
        estilos_ran <- list("Strict" = createStyle(fgFill="#D4EDDA", fontColour="#155724", textDecoration="bold"), "Wide" = createStyle(fgFill="#FFF3CD", fontColour="#856404"), "Outside" = createStyle(fgFill="#FFC7CE", fontColour="#9C0006", textDecoration="bold"))
        for (nm in names(estilos_ran)) {
          rows <- which(df$Tipo_Rango == nm)
          if (length(rows) > 0) addStyle(wb, "Variantes", estilos_ran[[nm]], rows = rows + 1, cols = col_ran)
        }
      }
      setColWidths(wb, "Variantes", cols = seq_len(ncol(df)), widths = "auto")
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$btn_descargar_csv <- downloadHandler(
    filename = function() paste0("variantes_filtradas_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- datos_filtrados()
      df <- df[, !names(df) %in% "._archivo_", drop = FALSE]
      write.csv(df, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )
  
  # ── Informe PDF por familia ────────────────────────────────────────────────
  output$btn_pdf_familia <- downloadHandler(
    filename = function() {
      fam <- input$pdf_familia_sel
      if (!nzchar(fam)) fam <- "familia"
      paste0("informe_", gsub("[^A-Za-z0-9_-]", "_", fam), "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      fam <- input$pdf_familia_sel
      
      if (!nzchar(fam)) {
        showNotification("Selecciona una familia antes de generar el informe.", type = "warning", duration = 4)
        grDevices::pdf(file, width = 11.69, height = 8.27)
        grid::grid.text("Sin familia seleccionada.", x = 0.5, y = 0.5, gp = grid::gpar(fontsize = 14))
        grDevices::dev.off()
        return(invisible(NULL))
      }
      
      todas_df <- tryCatch({
        bind_rows(
          if (!is.null(rv$resultados$CNVs)) mutate(rv$resultados$CNVs, Modalidad = "CNV") else data.frame(),
          if (!is.null(rv$resultados$SVs))  mutate(rv$resultados$SVs,  Modalidad = "SV")  else data.frame()
        )
      }, error = function(e) data.frame())
      
      if (nrow(todas_df) == 0) {
        showNotification("No hay datos cargados. Usa 'Cargar / Actualizar' primero.", type = "error", duration = 5)
        grDevices::pdf(file, width = 11.69, height = 8.27)
        grid::grid.text("Sin datos cargados.", x = 0.5, y = 0.5, gp = grid::gpar(fontsize = 14))
        grDevices::dev.off()
        return(invisible(NULL))
      }
      
      n_fichas_est <- nrow(tryCatch({
        df_tmp <- todas_df[as.character(todas_df$ID_Familia) == fam, , drop = FALSE]
        if ("Annotation_mode" %in% names(df_tmp))
          df_tmp <- df_tmp[!is.na(df_tmp$Annotation_mode) & df_tmp$Annotation_mode == "full", , drop = FALSE]
        sc_tmp <- suppressWarnings(as.numeric(df_tmp$AnnotSV_ranking_score))
        her_tmp <- her_base(df_tmp$Tipo_Herencia)
        sfari_tmp <- if ("En_Panel_Genes" %in% names(df_tmp)) es_sfari_col(df_tmp$En_Panel_Genes) else rep(FALSE, nrow(df_tmp))
        df_tmp <- df_tmp[(!is.na(sc_tmp) & sc_tmp >= 0) | (!is.na(her_tmp) & her_tmp == "De novo") | sfari_tmp, , drop = FALSE]
        head(df_tmp, 12)
      }, error = function(e) data.frame()))
      n_pasos_total <- 4L + max(n_fichas_est, 0L)  # 4 paginas fijas + 1 paso por variante
      
      grDevices::pdf(file, width = 11.69, height = 8.27, paper = "a4r")
      on.exit(grDevices::dev.off(), add = TRUE)
      
      withProgress(message = "Generando informe PDF…", value = 0, {
        tryCatch({
          
          # ── Datos de la familia ─────────────────────────────────────────────
          df_fam <- todas_df[as.character(todas_df$ID_Familia) == fam, , drop = FALSE]
          if ("Annotation_mode" %in% names(df_fam))
            df_full <- df_fam[!is.na(df_fam$Annotation_mode) & df_fam$Annotation_mode == "full", , drop = FALSE]
          else
            df_full <- df_fam
          if (nrow(df_full) == 0) df_full <- df_fam
          
          # Numeric and derived fields
          df_full$score_n  <- suppressWarnings(as.numeric(df_full$AnnotSV_ranking_score))
          df_full$her_b    <- her_base(df_full$Tipo_Herencia)
          df_full$chr_c    <- toupper(trimws(gsub("(?i)^chr","",as.character(df_full$SV_chrom),perl=TRUE)))
          df_full$sta_n    <- suppressWarnings(as.numeric(df_full$SV_start))
          df_full$end_n    <- suppressWarnings(as.numeric(df_full$SV_end))
          df_full$len_n    <- pmax(0, df_full$end_n - df_full$sta_n)
          df_full$len_lbl  <- sapply(df_full$len_n, function(l) {
            if (is.na(l) || l == 0) "--"
            else if (l >= 1e6) paste0(round(l/1e6, 2), " Mb")
            else paste0(round(l/1e3, 1), " kb")
          })
          sfari_col <- if ("En_Panel_Genes" %in% names(df_full))
            es_sfari_col(df_full$En_Panel_Genes)
          else rep(FALSE, nrow(df_full))
          df_full$en_sfari <- sfari_col
          
          # Unique key per variant (to access classifications and notes)
          claves_fila <- tryCatch(hacer_clave_variante(df_full), error = function(e) rep(NA_character_, nrow(df_full)))
          
          # Recuperar clasificaciones y notas
          clasif_vec <- sapply(claves_fila, function(k) {
            if (is.na(k)) return("--")
            cl <- rv$clasificaciones[[k]]
            if (is.null(cl) || !nzchar(cl)) "--" else cl
          })
          notas_vec <- sapply(claves_fila, function(k) {
            if (is.na(k)) return("")
            n <- rv$notas[[k]]
            if (is.null(n) || !nzchar(trimws(n))) "" else trimws(n)
          })
          df_full$clasificacion <- clasif_vec
          df_full$nota          <- notas_vec
          
          # Sexo del probando
          sexo_fam <- rv$sexos[[fam]]
          sexo_lbl <- if (is.null(sexo_fam) || !nzchar(sexo_fam)) "No registrado"
          else if (sexo_fam == "M") "Varon (masculino)"
          else if (sexo_fam == "F") "Mujer (femenino)"
          else sexo_fam
          fenotipo_fam <- rv$fenotipos[[fam]] %||% ""
          fenotipo_lbl <- if (nzchar(trimws(fenotipo_fam))) fenotipo_fam else "No registrado"
          
          # URLs externas
          df_full$url_franklin <- paste0("https://franklin.genoox.com/clinical-db/variant/sv/chr",
                                         df_full$chr_c,"-",df_full$sta_n,"-",df_full$end_n)
          df_full$url_decipher <- paste0("https://www.deciphergenomics.org/browser#q/",
                                         df_full$chr_c,":",df_full$sta_n,"-",df_full$end_n)
          df_full$url_ucsc     <- paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr",
                                         df_full$chr_c,":",df_full$sta_n,"-",df_full$end_n)
          df_full$url_gnomad   <- paste0("https://gnomad.broadinstitute.org/region/",
                                         df_full$chr_c,":",df_full$sta_n,"-",df_full$end_n,"?dataset=gnomad_sv_r4")
          df_full$url_omim     <- if ("Gene_name" %in% names(df_full))
            paste0("https://www.omim.org/search?search=", gsub(" ","+",(sapply(as.character(df_full$Gene_name),
                                                                               function(g) strsplit(g,"[,/|]")[[1]][1]))))
          else rep("https://www.omim.org", nrow(df_full))
          
          # Paleta herencia
          her_pal <- c(
            "De novo"            = "#C0392B",
            "Paternal"            = "#2980B9",
            "Maternal"            = "#E67E22",
            "Combined"           = "#8E44AD",
            "Paternal (Probable)" = "#5DADE2",
            "Maternal (Probable)" = "#F0A500",
            "Combined (Probable)"= "#A569BD",
            "Desconocido"        = "#7F8C8D",
            "Unknown"        = "#7F8C8D",
            "Sin datos"          = "#95A5A6"
          )
          
          # Colores clasificacion
          clasif_col <- c(
            "Pathogenic" = "#C0392B", "Pathogenic" = "#C0392B",
            "VUS"        = "#E67E22",
            "Benign"    = "#27AE60",
            "In Progress"= "#2980B9",
            "--"         = "#95A5A6"
          )
          
          # Tema ggplot limpio
          tema_pdf <- theme_minimal(base_size = 9) +
            theme(
              plot.title       = element_text(size = 10, face = "bold", color = "#2C3E50", margin = margin(b = 4)),
              plot.subtitle    = element_text(size = 8, color = "#7F8C8D", margin = margin(b = 3)),
              plot.background  = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "#FAFAFA", color = NA),
              panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.4),
              panel.grid.minor = element_blank(),
              axis.text        = element_text(size = 8, color = "#555"),
              axis.title       = element_text(size = 8.5, color = "#444"),
              legend.position  = "bottom",
              legend.key.size  = unit(0.35, "cm"),
              legend.text      = element_text(size = 7.5),
              legend.title     = element_text(size = 8, face = "bold"),
              plot.margin      = margin(6, 8, 6, 8)
            )
          
          # ── Summary statistics ──────────────────────────────────────────────
          n_total    <- nrow(df_full)
          n_denovo   <- sum(df_full$her_b == "De novo", na.rm = TRUE)
          n_sfari    <- sum(df_full$en_sfari, na.rm = TRUE)
          n_patog    <- sum(!is.na(df_full$score_n) & df_full$score_n >= 0.5, na.rm = TRUE)
          n_pospatog <- sum(!is.na(df_full$score_n) & df_full$score_n >= 0 & df_full$score_n < 0.5, na.rm = TRUE)
          n_benign   <- sum(!is.na(df_full$score_n) & df_full$score_n < 0, na.rm = TRUE)
          n_clasif   <- sum(df_full$clasificacion != "--")
          score_max  <- if (any(is.finite(df_full$score_n))) max(df_full$score_n, na.rm = TRUE) else NA
          score_med  <- if (any(is.finite(df_full$score_n))) round(mean(df_full$score_n, na.rm = TRUE), 3) else NA
          n_del      <- sum(grepl("DEL|LOSS", toupper(as.character(df_full$SV_type))), na.rm = TRUE)
          n_dup      <- sum(grepl("DUP|GAIN", toupper(as.character(df_full$SV_type))), na.rm = TRUE)
          tam_med_kb <- if (any(!is.na(df_full$len_n) & df_full$len_n > 0))
            round(median(df_full$len_n[df_full$len_n > 0], na.rm = TRUE) / 1e3, 1)
          else NA
          
          # Helper: dibuja un KPI box
          kpi_box <- function(x, y, w, h, titulo, valor, subtitulo = NULL,
                              bg = "#2C3E50", fg = "white", valor_col = "white") {
            grid::grid.rect(x = x, y = y, width = w, height = h, just = c("left","top"),
                            gp = grid::gpar(fill = bg, col = NA))
            grid::grid.text(toupper(titulo), x = x + 0.01, y = y - 0.005,
                            just = c("left","top"),
                            gp = grid::gpar(fontsize = 6.5, fontface = "bold",
                                            col = if (bg == "white") "#888" else "#CCCCCC",
                                            letterSpacing = 1))
            grid::grid.text(as.character(valor), x = x + w/2, y = y - h*0.42,
                            just = c("center","center"),
                            gp = grid::gpar(fontsize = 18, fontface = "bold", col = valor_col))
            if (!is.null(subtitulo))
              grid::grid.text(subtitulo, x = x + w/2, y = y - h*0.82,
                              just = c("center","center"),
                              gp = grid::gpar(fontsize = 7, col = if (bg == "white") "#888" else "#BBBBBB"))
          }
          
          # Helper: draws section header band
          seccion_header <- function(texto, y_pos, height = 0.055,
                                     bg = "#2C3E50", fg = "white", subtexto = NULL) {
            grid::grid.rect(x = 0, y = y_pos, width = 1, height = height, just = c("left","top"),
                            gp = grid::gpar(fill = bg, col = NA))
            grid::grid.rect(x = 0, y = y_pos, width = 0.006, height = height, just = c("left","top"),
                            gp = grid::gpar(fill = "#F39C12", col = NA))
            grid::grid.text(texto, x = 0.015, y = y_pos - height/2,
                            just = c("left","center"),
                            gp = grid::gpar(fontsize = 10, fontface = "bold", col = fg))
            if (!is.null(subtexto))
              grid::grid.text(subtexto, x = 0.98, y = y_pos - height/2,
                              just = c("right","center"),
                              gp = grid::gpar(fontsize = 8, col = "#BDC3C7"))
          }
          
          # Helper: footer
          pie_pagina <- function(pagina, total_pags = "?") {
            grid::grid.rect(x = 0, y = 0.025, width = 1, height = 0.025, just = c("left","top"),
                            gp = grid::gpar(fill = "#F8F9FA", col = NA))
            grid::grid.text(
              paste0("CNV/SV Pipeline  |  Familia: ", fam,
                     "  |  Generado: ", format(Sys.time(), "%d/%m/%Y %H:%M"),
                     "  |  Genoma: hg38  |  Pipeline: AnnotSV"),
              x = 0.015, y = 0.013, just = c("left","center"),
              gp = grid::gpar(fontsize = 6.5, col = "#888888"))
            grid::grid.text(paste0("Pag. ", pagina), x = 0.985, y = 0.013,
                            just = c("right","center"),
                            gp = grid::gpar(fontsize = 6.5, col = "#888888"))
          }
          
          # ─────────────────────────────────────────────────────────────────────
          # PAGINA 1 — Portada + Resumen ejecutivo
          # ─────────────────────────────────────────────────────────────────────
          incProgress(1 / (n_pasos_total + 1), detail = "Portada y resumen…")
          grid::grid.newpage()
          
          # Banda superior de portada
          grid::grid.rect(x = 0, y = 1, width = 1, height = 0.17,
                          just = c("left","top"), gp = grid::gpar(fill = "#2C3E50", col = NA))
          grid::grid.rect(x = 0, y = 1, width = 1, height = 0.003,
                          just = c("left","top"), gp = grid::gpar(fill = "#F39C12", col = NA))
          
          # Main title
          grid::grid.text("INFORME CLINICO  CNV / SV", x = 0.05, y = 0.965,
                          just = c("left","top"),
                          gp = grid::gpar(fontsize = 22, fontface = "bold", col = "white"))
          grid::grid.text("CNV/SV Pipeline  |  Analisis de variantes de numero de copia y variantes estructurales",
                          x = 0.05, y = 0.925, just = c("left","top"),
                          gp = grid::gpar(fontsize = 10, col = "#BDC3C7"))
          
          # Datos identificativos (lado derecho de cabecera)
          grid::grid.text(paste0("ID Familia:  ", fam), x = 0.96, y = 0.968,
                          just = c("right","top"),
                          gp = grid::gpar(fontsize = 11, fontface = "bold", col = "white"))
          grid::grid.text(paste0("Sexo del probando:  ", sexo_lbl), x = 0.96, y = 0.946,
                          just = c("right","top"),
                          gp = grid::gpar(fontsize = 9, col = "#BDC3C7"))
          grid::grid.text(paste0("Fecha informe:  ", format(Sys.Date(), "%d/%m/%Y")), x = 0.96, y = 0.926,
                          just = c("right","top"),
                          gp = grid::gpar(fontsize = 9, col = "#BDC3C7"))
          # Fenotipo en cabecera (truncado si es largo)
          fenotipo_hdr <- if (nchar(fenotipo_lbl) > 95)
            paste0(substr(fenotipo_lbl, 1, 93), "\u2026") else fenotipo_lbl
          grid::grid.text(paste0("Fenotipo:  ", fenotipo_hdr), x = 0.05, y = 0.908,
                          just = c("left","top"),
                          gp = grid::gpar(fontsize = 8.5, col = "#BDC3C7",
                                          fontface = "italic"))
          
          # ── Fila de KPIs ──────────────────────────────────────────────────────
          kpi_y     <- 0.815; kpi_h <- 0.105; kpi_w <- 0.148; kpi_gap <- 0.012
          kpi_left  <- 0.025
          kpi_data  <- list(
            list(v = n_total,    t = "Total variantes",      s = "familia completa",   bg = "#2980B9"),
            list(v = n_denovo,   t = "De novo",              s = "herencia de novo",   bg = "#C0392B"),
            list(v = n_sfari,    t = "Panel genes",          s = "in gene panel",     bg = "#E67E22"),
            list(v = n_patog,    t = "Patogenicas",          s = "score >= 0.5",       bg = "#922B21"),
            list(v = if(!is.na(score_max)) round(score_max,3) else "--",
                 t = "Score maximo", s = "AnnotSV ranking",  bg = "#1A5276"),
            list(v = paste0(n_del,"/",n_dup), t = "DEL / DUP",
                 s = "tipo de variante",                     bg = "#117A65")
          )
          for (ki in seq_along(kpi_data)) {
            kpi_x <- kpi_left + (ki - 1) * (kpi_w + kpi_gap)
            kd    <- kpi_data[[ki]]
            kpi_box(kpi_x, kpi_y, kpi_w, kpi_h,
                    titulo    = kd$t,
                    valor     = kd$v,
                    subtitulo = kd$s,
                    bg        = kd$bg)
          }
          
          # ── Segunda fila de KPIs ──────────────────────────────────────────────
          kpi_y2 <- kpi_y - kpi_h - 0.01
          kpi_data2 <- list(
            list(v = n_pospatog, t = "Potenc. patogenicas", s = "score 0 - 0.5",       bg = "#D35400"),
            list(v = n_benign,   t = "Benignas",            s = "score < 0",           bg = "#27AE60"),
            list(v = n_clasif,   t = "Clasificadas",        s = "con clasificacion",   bg = "#5D6D7E"),
            list(v = if(!is.na(score_med)) score_med else "--",
                 t = "Score medio", s = "AnnotSV mean",     bg = "#2471A3"),
            list(v = if(!is.na(tam_med_kb)) paste0(tam_med_kb," kb") else "--",
                 t = "Tamano mediano", s = "variante",      bg = "#148F77"),
            list(v = n_pospatog + n_patog, t = "Score >= 0",
                 s = "revision prioritaria",                bg = "#6C3483")
          )
          for (ki in seq_along(kpi_data2)) {
            kpi_x <- kpi_left + (ki - 1) * (kpi_w + kpi_gap)
            kd    <- kpi_data2[[ki]]
            kpi_box(kpi_x, kpi_y2, kpi_w, kpi_h,
                    titulo    = kd$t, valor = kd$v, subtitulo = kd$s, bg = kd$bg)
          }
          
          # ── Tabla de variantes prioritarias (de novo + score alto) ───────────
          df_prio <- df_full[
            (!is.na(df_full$her_b) & df_full$her_b == "De novo") |
              (!is.na(df_full$score_n) & df_full$score_n >= 0.3) |
              df_full$en_sfari,
            , drop = FALSE
          ]
          df_prio <- df_prio[order(-ifelse(is.na(df_prio$score_n), -99, df_prio$score_n)), , drop = FALSE]
          n_prio  <- min(nrow(df_prio), 12)
          
          seccion_header(
            texto    = "Variantes prioritarias (De novo | Score >= 0.3 | SFARI)",
            y_pos    = kpi_y2 - kpi_h - 0.022,
            subtexto = paste0("Mostrando top ", n_prio, " of ", nrow(df_prio))
          )
          
          if (n_prio > 0) {
            df_p_show <- df_prio[seq_len(n_prio), , drop = FALSE]
            df_p_tbl  <- data.frame(
              Region     = paste0("chr", df_p_show$chr_c, ":", format(round(df_p_show$sta_n), big.mark=".", scientific=FALSE),
                                  "-", format(round(df_p_show$end_n), big.mark=".", scientific=FALSE)),
              Tamano     = df_p_show$len_lbl,
              Tipo       = if ("SV_type" %in% names(df_p_show)) as.character(df_p_show$SV_type) else "--",
              Herencia   = df_p_show$her_b,
              Rango      = if ("Tipo_Rango" %in% names(df_p_show)) as.character(df_p_show$Tipo_Rango) else "--",
              Score      = ifelse(is.na(df_p_show$score_n), "--", sprintf("%.3f", df_p_show$score_n)),
              Genes      = if ("Gene_name" %in% names(df_p_show)) {
                sapply(as.character(df_p_show$Gene_name), function(g) {
                  g <- trimws(g); if (is.na(g) || g == "" || g == "NA" || g == ".") return("--")
                  genes_split <- strsplit(g, "[,/|]")[[1]]
                  if (length(genes_split) > 3) paste0(paste(genes_split[1:3], collapse=", "), "...") else g
                })
              } else "--",
              SFARI      = ifelse(df_p_show$en_sfari, "Yes", ""),
              Clasif     = df_p_show$clasificacion,
              stringsAsFactors = FALSE
            )
            y_tbl_start <- kpi_y2 - kpi_h - 0.09
            vp_tbl <- grid::viewport(x = 0.025, y = y_tbl_start,
                                     width = 0.95, height = 0.38,
                                     just = c("left","top"))
            grid::pushViewport(vp_tbl)
            n_rows_p <- nrow(df_p_tbl)
            row_bg <- rep(c("white","#F2F3F4"), length.out = n_rows_p)
            tt_prio <- gridExtra::tableGrob(
              df_p_tbl, rows = NULL,
              theme = gridExtra::ttheme_minimal(
                core    = list(fg_params = list(fontsize = 6.8, hjust = 0, x = 0.03),
                               bg_params = list(fill = row_bg)),
                colhead = list(fg_params = list(fontsize = 7.5, fontface = "bold", col = "#2C3E50"),
                               bg_params = list(fill = "#D6EAF8"))
              )
            )
            grid::grid.draw(tt_prio)
            grid::popViewport()
          } else {
            grid::grid.text("No hay variantes prioritarias para esta familia.",
                            x = 0.5, y = kpi_y2 - kpi_h - 0.13, just = c("center","center"),
                            gp = grid::gpar(fontsize = 9, col = "#888"))
          }
          
          pie_pagina(1)
          
          # ─────────────────────────────────────────────────────────────────────
          # PAGINA 2 — Tabla completa de variantes
          # ─────────────────────────────────────────────────────────────────────
          incProgress(1 / (n_pasos_total + 1), detail = "Tabla completa de variantes…")
          grid::grid.newpage()
          seccion_header("Tabla completa de variantes  —  Familia: ", y_pos = 1,
                         subtexto = paste0(n_total, " variantes | Ordenadas por score desc."))
          
          df_full_ord <- df_full[order(-ifelse(is.na(df_full$score_n), -99, df_full$score_n)), , drop = FALSE]
          n_show_full <- min(nrow(df_full_ord), 30)
          df_t2 <- df_full_ord[seq_len(n_show_full), , drop = FALSE]
          
          df_tabla2 <- data.frame(
            "#"        = seq_len(n_show_full),
            Region     = paste0("chr", df_t2$chr_c, ":",
                                format(round(df_t2$sta_n), big.mark = ".", scientific = FALSE), "-",
                                format(round(df_t2$end_n), big.mark = ".", scientific = FALSE)),
            Tamano     = df_t2$len_lbl,
            Tipo       = if ("SV_type" %in% names(df_t2)) as.character(df_t2$SV_type) else "--",
            Herencia   = df_t2$her_b,
            Rango      = if ("Tipo_Rango" %in% names(df_t2)) as.character(df_t2$Tipo_Rango) else "--",
            Score      = ifelse(is.na(df_t2$score_n), "--", sprintf("%.3f", df_t2$score_n)),
            SFARI      = ifelse(df_t2$en_sfari, "Yes", ""),
            DRAGEN.Ex  = if ("Illumina_DRAGEN.exact.counts"   %in% names(df_t2)) as.character(df_t2[["Illumina_DRAGEN.exact.counts"]])   else "--",
            DRAGEN.Sim = if ("Illumina_DRAGEN.similar.counts" %in% names(df_t2)) as.character(df_t2[["Illumina_DRAGEN.similar.counts"]]) else "--",
            Confianza  = if ("Confianza_Region" %in% names(df_t2)) as.character(df_t2$Confianza_Region) else "--",
            Clasif     = df_t2$clasificacion,
            Genes      = if ("Gene_name" %in% names(df_t2)) {
              sapply(as.character(df_t2$Gene_name), function(g) {
                g <- trimws(g); if (is.na(g) || g == "" || g == "NA" || g == ".") "--"
                else { gs <- strsplit(g, "[,/|]")[[1]]; if (length(gs) > 4) paste0(paste(gs[1:4], collapse=", "), "...") else g }
              })
            } else "--",
            stringsAsFactors = FALSE
          )
          names(df_tabla2)[1] <- "#"
          
          # Row colours by score
          score_num_ord <- df_t2$score_n
          row_bg2 <- sapply(seq_len(n_show_full), function(i) {
            s <- score_num_ord[i]
            if (is.na(s))     return(if (i%%2==0) "#F8F9FA" else "white")
            if (s >= 0.5)     return("#FDEDEC")
            if (s >= 0)       return("#FEF9E7")
            return(if (i%%2==0) "#F8F9FA" else "white")
          })
          # Resaltar de novo
          her_ord <- df_t2$her_b
          row_bg2[!is.na(her_ord) & her_ord == "De novo"] <- "#FDEDEC"
          
          # Reserve space for footer (0.03) and optional note (0.04)
          vp2_y      <- 0.935
          vp2_height <- 0.895
          vp2 <- grid::viewport(x = 0.015, y = vp2_y, width = 0.97, height = vp2_height,
                                just = c("left","top"), clip = "on")
          grid::pushViewport(vp2)
          tt2 <- gridExtra::tableGrob(
            df_tabla2, rows = NULL,
            theme = gridExtra::ttheme_minimal(
              core    = list(fg_params = list(fontsize = 6.5, hjust = 0, x = 0.03),
                             bg_params = list(fill = row_bg2)),
              colhead = list(fg_params = list(fontsize = 7.5, fontface = "bold", col = "#2C3E50"),
                             bg_params = list(fill = "#EAF2FF"))
            )
          )
          # Normalizar anchos y altos al viewport (evita desbordamiento y
          # garantiza que la tabla llene el espacio disponible)
          tt2$widths  <- unit(as.numeric(tt2$widths)  / sum(as.numeric(tt2$widths)),  "npc")
          tt2$heights <- unit(as.numeric(tt2$heights) / sum(as.numeric(tt2$heights)), "npc")
          grid::grid.draw(tt2)
          grid::popViewport()
          
          if (n_show_full < nrow(df_full_ord))
            grid::grid.text(paste0("* Mostrando ", n_show_full, " of ", nrow(df_full_ord),
                                   " variantes. Descarga el Excel para ver todas."),
                            x = 0.5, y = 0.04, just = c("center","center"),
                            gp = grid::gpar(fontsize = 7, col = "#888", fontface = "italic"))
          
          pie_pagina(2)
          
          # ─────────────────────────────────────────────────────────────────────
          # PAGINA 3 — Genes afectados + Herencia + Frecuencias
          # ─────────────────────────────────────────────────────────────────────
          incProgress(1 / (n_pasos_total + 1), detail = "Genes afectados y herencia…")
          grid::grid.newpage()
          seccion_header("Genes afectados y distribucion de herencia", y_pos = 1)
          
          # Tabla de genes
          if ("Gene_name" %in% names(df_full)) {
            genes_raw  <- as.character(df_full$Gene_name)
            genes_raw  <- genes_raw[!is.na(genes_raw) & nzchar(genes_raw) & genes_raw != "." & genes_raw != "NA"]
            genes_vec  <- trimws(unlist(strsplit(genes_raw, "[/|,]")))
            genes_vec  <- genes_vec[nzchar(genes_vec) & genes_vec != "."]
            
            if (length(genes_vec) > 0) {
              gene_tab  <- sort(table(genes_vec), decreasing = TRUE)
              top_genes <- head(gene_tab, 20)
              
              # Abreviaturas compactas para herencia
              abrev_her <- function(hers) {
                lut <- c("De novo"="DN","Paternal (Probable)"="PAT(P)","Maternal (Probable)"="MAT(P)",
                         "Combined (Probable)"="CONJ(P)","Paternal"="PAT","Maternal"="MAT",
                         "Combined"="CONJ","No determinada"="ND","Unknown"="ND")
                parts <- sapply(hers, function(h) { v <- lut[h]; if (is.na(v)) substr(h,1,6) else v })
                txt <- paste(parts, collapse="/")
                if (nchar(txt) > 20) paste0(substr(txt,1,18),"\u2026") else txt
              }
              
              gene_her <- sapply(names(top_genes), function(g) {
                mask <- grepl(g, as.character(df_full$Gene_name), fixed = TRUE)
                hers <- unique(df_full$her_b[mask & !is.na(df_full$her_b)])
                if (length(hers) == 0) "--" else abrev_her(hers)
              })
              gene_sfari <- sapply(names(top_genes), function(g) {
                mask <- grepl(g, as.character(df_full$Gene_name), fixed = TRUE)
                if (any(df_full$en_sfari[mask], na.rm = TRUE)) "\u2605" else ""
              })
              gene_score_max <- sapply(names(top_genes), function(g) {
                mask <- grepl(g, as.character(df_full$Gene_name), fixed = TRUE)
                sc <- df_full$score_n[mask]; sc <- sc[is.finite(sc)]
                if (length(sc) == 0) "--" else sprintf("%.3f", max(sc))
              })
              gene_nfam <- sapply(names(top_genes), function(g) {
                mask <- grepl(g, as.character(df_full$Gene_name), fixed = TRUE)
                length(unique(df_full$ID_Familia[mask]))
              })
              
              df_genes_tbl <- data.frame(
                Gen      = names(top_genes),
                N        = as.integer(top_genes),
                Fam      = gene_nfam,
                SFARI    = gene_sfari,
                Her      = gene_her,
                ScoreMax = gene_score_max,
                stringsAsFactors = FALSE
              )
              
              seccion_header("Top genes mas frecuentes", y_pos = 0.94,
                             bg = "#1A5276",
                             subtexto = paste0(length(genes_vec), " genes anotados  /  ",
                                               length(gene_tab), " genes unicos"))
              
              row_bg_g <- ifelse(df_genes_tbl$SFARI == "\u2605", "#FFFBEA",
                                 ifelse(seq_len(nrow(df_genes_tbl)) %% 2 == 0, "#F4F6F7", "white"))
              score_num_g <- suppressWarnings(as.numeric(df_genes_tbl$ScoreMax))
              score_col_g <- ifelse(!is.na(score_num_g) & score_num_g >= 0.5, "#922B21",
                                    ifelse(!is.na(score_num_g) & score_num_g >= 0, "#D35400", "#444444"))
              
              vp_gen <- grid::viewport(x = 0.015, y = 0.885, width = 0.44, height = 0.815,
                                       just = c("left","top"))
              grid::pushViewport(vp_gen)
              n_cols_g <- ncol(df_genes_tbl)
              score_idx <- which(names(df_genes_tbl) == "ScoreMax")
              tt_gen <- gridExtra::tableGrob(
                df_genes_tbl, rows = NULL,
                theme = gridExtra::ttheme_minimal(
                  core    = list(fg_params = list(fontsize = 7, hjust = 0, x = 0.05),
                                 bg_params = list(fill = row_bg_g)),
                  colhead = list(fg_params = list(fontsize = 7.5, fontface = "bold", col = "#1A5276"),
                                 bg_params = list(fill = "#D6EAF8"))
                )
              )
              # Colorear columna ScoreMax
              for (ri in seq_len(nrow(df_genes_tbl))) {
                idx <- (ri + 1) * n_cols_g - (n_cols_g - score_idx)
                tryCatch(
                  tt_gen$grobs[[idx]]$gp <- grid::gpar(fontsize=7, col=score_col_g[ri], fontface="bold"),
                  error = function(e) NULL
                )
              }
              grid::grid.draw(tt_gen)
              grid::popViewport()
            }
          }
          
          # Tabla resumen de herencia
          her_tab <- as.data.frame(table(Herencia = df_full$her_b), stringsAsFactors = FALSE)
          her_tab <- her_tab[her_tab$Freq > 0, ]
          her_tab <- her_tab[order(-her_tab$Freq), ]
          if (nrow(her_tab) > 0) {
            her_tab$Pct_total <- paste0(round(her_tab$Freq / sum(her_tab$Freq) * 100, 1), "%")
            n_sc_her <- sapply(her_tab$Herencia, function(h) {
              sc <- df_full$score_n[!is.na(df_full$her_b) & df_full$her_b == h]
              if (length(sc[is.finite(sc)]) == 0) "--"
              else sprintf("%.3f", median(sc, na.rm = TRUE))
            })
            her_tab$Score_med <- n_sc_her
            n_sfari_her <- sapply(her_tab$Herencia, function(h) {
              sum(df_full$en_sfari[!is.na(df_full$her_b) & df_full$her_b == h], na.rm = TRUE)
            })
            her_tab$N_SFARI <- n_sfari_her
            names(her_tab) <- c("Herencia","N variantes","% total","Score mediano","N SFARI")
            
            seccion_header("Distribucion por herencia", y_pos = 0.94,
                           bg = "#117A65")
            
            vp_her <- grid::viewport(x = 0.475, y = 0.885, width = 0.51, height = 0.45,
                                     just = c("left","top"))
            grid::pushViewport(vp_her)
            row_bg_h <- sapply(her_tab$Herencia, function(h) {
              col <- her_pal[h]; if (is.na(col)) "#F8F9FA"
              else paste0(col, "22")
            })
            tt_her <- gridExtra::tableGrob(
              her_tab, rows = NULL,
              theme = gridExtra::ttheme_minimal(
                core    = list(fg_params = list(fontsize = 7.5, hjust = 0, x = 0.03),
                               bg_params = list(fill = row_bg_h)),
                colhead = list(fg_params = list(fontsize = 8, fontface = "bold", col = "#117A65"),
                               bg_params = list(fill = "#D5F5E3"))
              )
            )
            grid::grid.draw(tt_her)
            grid::popViewport()
          }
          
          # Tabla de frecuencias DRAGEN
          has_dragen <- "Illumina_DRAGEN.exact.counts" %in% names(df_full)
          if (has_dragen) {
            df_dragen <- data.frame(
              Region     = paste0("chr", df_full$chr_c, ":",
                                  format(round(df_full$sta_n), big.mark=".", scientific=FALSE)),
              Tamano     = df_full$len_lbl,
              Score      = ifelse(is.na(df_full$score_n),"--", sprintf("%.3f",df_full$score_n)),
              DRAGEN.Ex  = as.character(df_full[["Illumina_DRAGEN.exact.counts"]]),
              DRAGEN.Sim = if ("Illumina_DRAGEN.similar.counts" %in% names(df_full))
                as.character(df_full[["Illumina_DRAGEN.similar.counts"]]) else "--",
              Ref.Match  = if ("Referencia_Match" %in% names(df_full)) as.character(df_full$Referencia_Match) else "--",
              Sim.MaxLen = if ("Sim_MaxLen" %in% names(df_full)) as.character(df_full$Sim_MaxLen) else "--",
              Confianza  = if ("Confianza_Region" %in% names(df_full)) as.character(df_full$Confianza_Region) else "--",
              stringsAsFactors = FALSE
            )
            df_dragen <- df_dragen[order(-ifelse(is.na(df_full$score_n),-99,df_full$score_n)),]
            df_dragen <- head(df_dragen, 20)
            
            seccion_header("Frecuencias poblacionales DRAGEN y confianza de region",
                           y_pos = 0.39, bg = "#6C3483",
                           subtexto = "DRAGEN.Ex = n exactas | DRAGEN.Sim = n similares | Ref.Match = solapamiento referencia")
            
            vp_dr <- grid::viewport(x = 0.475, y = 0.335, width = 0.51, height = 0.32,
                                    just = c("left","top"))
            grid::pushViewport(vp_dr)
            row_bg_dr <- rep(c("white","#F5EEF8"), length.out = nrow(df_dragen))
            tt_dr <- gridExtra::tableGrob(
              df_dragen, rows = NULL,
              theme = gridExtra::ttheme_minimal(
                core    = list(fg_params = list(fontsize = 6.5, hjust = 0, x = 0.02),
                               bg_params = list(fill = row_bg_dr)),
                colhead = list(fg_params = list(fontsize = 7, fontface = "bold", col = "#6C3483"),
                               bg_params = list(fill = "#E8DAEF"))
              )
            )
            grid::grid.draw(tt_dr)
            grid::popViewport()
          }
          
          pie_pagina(3)
          
          # ─────────────────────────────────────────────────────────────────────
          # PAGINA 4 — Graficas estadisticas
          # ─────────────────────────────────────────────────────────────────────
          incProgress(1 / (n_pasos_total + 1), detail = "Statistical plots…")
          grid::grid.newpage()
          seccion_header(paste0("Graficas estadisticas  —  Familia: ", fam), y_pos = 1,
                         subtexto = paste0(n_total, " variantes | ",
                                           sum(!is.na(df_full$score_n)), " con score AnnotSV"))
          
          # Plot 1: Score distribution
          g_scores <- tryCatch({
            sc <- df_full$score_n[is.finite(df_full$score_n)]
            if (length(sc) > 0) {
              df_sc <- data.frame(
                score = sc,
                categoria = cut(sc, breaks = c(-Inf, 0, 0.5, Inf),
                                labels = c("Benigna (< 0)", "Incierta (0-0.5)", "Patogenica (>= 0.5)"),
                                include.lowest = TRUE)
              )
              ggplot(df_sc, aes(x = score, fill = categoria)) +
                geom_histogram(binwidth = 0.05, color = "white", linewidth = 0.3) +
                geom_vline(xintercept = 0,   color = "#E74C3C", linetype = "dashed", linewidth = 0.8) +
                geom_vline(xintercept = 0.5, color = "#C0392B", linetype = "dashed", linewidth = 0.8) +
                scale_fill_manual(values = c("#3498DB","#F39C12","#C0392B"), name = "") +
                labs(title = "Distribucion de scores AnnotSV",
                     subtitle = paste0("n=", length(sc), " | mediana=", round(median(sc),3)),
                     x = "Score AnnotSV", y = "N variantes") +
                tema_pdf + theme(legend.position = "bottom")
            } else {
              ggplot() + annotate("text",x=.5,y=.5,label="Sin scores") + tema_pdf +
                labs(title="Scores AnnotSV")
            }
          }, error = function(e) ggplot()+annotate("text",x=.5,y=.5,label="Error")+tema_pdf)
          
          # Grafico 2: Herencia
          g_herencia <- tryCatch({
            tbl_h <- as.data.frame(table(her = df_full$her_b))
            tbl_h <- tbl_h[!is.na(tbl_h$her) & tbl_h$her != "" & tbl_h$Freq > 0, ]
            if (nrow(tbl_h) > 0) {
              tbl_h$pct  <- round(tbl_h$Freq / sum(tbl_h$Freq) * 100, 1)
              cols_h <- her_pal[as.character(tbl_h$her)]
              cols_h[is.na(cols_h)] <- "#95A5A6"
              names(cols_h) <- as.character(tbl_h$her)
              ggplot(tbl_h, aes(x = reorder(her, -Freq), y = Freq, fill = her)) +
                geom_col(color = "white", linewidth = 0.3, width = 0.65) +
                geom_text(aes(label = paste0(Freq, "\n(", pct, "%)")), vjust = -0.3, size = 2.5) +
                scale_fill_manual(values = cols_h, guide = "none") +
                labs(title = "Tipo de herencia",
                     subtitle = paste0(n_total, " variantes en total"),
                     x = "", y = "N variantes") +
                tema_pdf + theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 7))
            } else ggplot()+annotate("text",x=.5,y=.5,label="Sin datos")+tema_pdf+labs(title="Herencia")
          }, error = function(e) ggplot()+annotate("text",x=.5,y=.5,label="Error")+tema_pdf)
          
          # Grafico 3: Tipo SV (DEL/DUP)
          g_tipo <- tryCatch({
            df_tp <- df_full
            sv_up <- toupper(as.character(df_tp$SV_type))
            df_tp$tipo_agrup <- ifelse(grepl("DEL|LOSS",sv_up),"DEL",
                                       ifelse(grepl("DUP|GAIN",sv_up),"DUP","OTHER"))
            tbl_t <- as.data.frame(table(tipo = df_tp$tipo_agrup))
            tbl_t <- tbl_t[tbl_t$Freq > 0, ]
            if (nrow(tbl_t) > 0) {
              cols_t <- c(DEL="#C0392B",DUP="#2980B9",OTRO="#7F8C8D")
              c_fill <- cols_t[as.character(tbl_t$tipo)]; c_fill[is.na(c_fill)] <- "#888"
              names(c_fill) <- as.character(tbl_t$tipo)
              tbl_t$pct <- round(tbl_t$Freq/sum(tbl_t$Freq)*100,1)
              ggplot(tbl_t, aes(x = tipo, y = Freq, fill = tipo)) +
                geom_col(color = "white", width = 0.55) +
                geom_text(aes(label = paste0(Freq,"\n(",pct,"%)")), vjust = -0.3, size = 2.8) +
                scale_fill_manual(values = c_fill, guide = "none") +
                labs(title = "Tipo de variante", subtitle = "DEL = delecion | DUP = duplicacion",
                     x = "", y = "N variantes") +
                tema_pdf
            } else ggplot()+annotate("text",x=.5,y=.5,label="Sin datos")+tema_pdf+labs(title="Tipo SV")
          }, error = function(e) ggplot()+annotate("text",x=.5,y=.5,label="Error")+tema_pdf)
          
          # Grafico 4: Score por herencia boxplot
          g_boxplot <- tryCatch({
            df_bx <- df_full[is.finite(df_full$score_n) & !is.na(df_full$her_b), ]
            if (nrow(df_bx) >= 2) {
              cols_bx <- her_pal[unique(as.character(df_bx$her_b))]
              cols_bx[is.na(cols_bx)] <- "#95A5A6"
              ggplot(df_bx, aes(x = her_b, y = score_n, fill = her_b)) +
                geom_boxplot(outlier.shape = 21, outlier.size = 1.5, width = 0.55, color = "#444") +
                geom_jitter(width = 0.1, alpha = 0.55, size = 1.5, aes(color = her_b)) +
                scale_fill_manual(values = cols_bx, guide = "none") +
                scale_color_manual(values = cols_bx, guide = "none") +
                geom_hline(yintercept = 0.5, color="#C0392B", linetype="dashed", linewidth=0.6) +
                geom_hline(yintercept = 0,   color="#E67E22", linetype="dashed", linewidth=0.5) +
                labs(title = "Score AnnotSV por herencia",
                     subtitle = "Linea roja=0.5 | Linea naranja=0",
                     x = "", y = "Score") +
                tema_pdf + theme(axis.text.x = element_text(angle = 22, hjust = 1, size = 7))
            } else {
              ggplot(df_full[is.finite(df_full$score_n),], aes(x = "Todas", y = score_n)) +
                geom_jitter(width=0.1, color="#2980B9", size=2, alpha=0.7) +
                geom_hline(yintercept=0.5, color="#C0392B", linetype="dashed", linewidth=0.6) +
                labs(title = "Score AnnotSV", x = "", y = "Score") + tema_pdf
            }
          }, error = function(e) ggplot()+annotate("text",x=.5,y=.5,label="Error")+tema_pdf)
          
          # Grafico 5: Rango
          g_rango <- tryCatch({
            tbl_r <- if ("Tipo_Rango" %in% names(df_full)) {
              as.data.frame(table(rango = df_full$Tipo_Rango))
            } else data.frame()
            if (nrow(tbl_r) > 0) {
              tbl_r <- tbl_r[tbl_r$Freq > 0, ]
              cols_r <- c("Strict"="#155724","Wide"="#856404","Outside"="#922B21")
              c_f <- cols_r[as.character(tbl_r$rango)]; c_f[is.na(c_f)] <- "#888"
              names(c_f) <- as.character(tbl_r$rango)
              tbl_r$pct <- round(tbl_r$Freq/sum(tbl_r$Freq)*100,1)
              ggplot(tbl_r, aes(x = reorder(rango,-Freq), y = Freq, fill = rango)) +
                geom_col(color = "white", width = 0.5) +
                geom_text(aes(label = paste0(Freq,"\n(",pct,"%)")), vjust = -0.3, size = 2.8) +
                scale_fill_manual(values = c_f, guide = "none") +
                labs(title = "Rango de la variante", subtitle = "Estricto > Amplio > Fuera de rango",
                     x = "", y = "N variantes") +
                tema_pdf
            } else ggplot()+annotate("text",x=.5,y=.5,label="Sin datos")+tema_pdf+labs(title="Rango")
          }, error = function(e) ggplot()+annotate("text",x=.5,y=.5,label="Error")+tema_pdf)
          
          # Grafico 6: Ideograma
          g_ideograma <- tryCatch({
            df_id <- df_full[df_full$chr_c %in% CHR_ORDER & !is.na(df_full$sta_n) & !is.na(df_full$end_n), ]
            df_id$chr_f  <- factor(df_id$chr_c, levels = rev(CHR_ORDER))
            chr_df <- data.frame(
              chr = factor(CHR_ORDER, levels = rev(CHR_ORDER)),
              len = as.numeric(CHR_LENGTHS[CHR_ORDER])
            )
            df_id$min_w <- pmax(df_id$len_n,
                                chr_df$len[match(df_id$chr_c, CHR_ORDER)] * 0.02, na.rm = TRUE)
            cols_id <- her_pal[as.character(df_id$her_b)]; cols_id[is.na(cols_id)] <- "#95A5A6"
            ggplot() +
              geom_segment(data = chr_df, aes(x = 0, xend = len, y = chr, yend = chr),
                           color = "#E0E0E0", linewidth = 3.5, lineend = "round") +
              geom_segment(data = df_id,
                           aes(x = sta_n, xend = sta_n + min_w, y = chr_f, yend = chr_f,
                               color = her_b),
                           linewidth = 4.5, lineend = "round", alpha = 0.88) +
              scale_color_manual(values = her_pal, na.value = "#95A5A6", name = "Herencia") +
              scale_x_continuous(labels = function(x) paste0(round(x/1e6), "Mb")) +
              labs(title = paste0("Ideograma cromosomico  —  ", nrow(df_id), " variantes"),
                   subtitle = "Referencia: hg38",
                   x = "Posicion genomica", y = "") +
              tema_pdf + theme(legend.position = "right", axis.text.y = element_text(size = 6.5),
                               legend.text = element_text(size = 6.5), legend.key.size = unit(0.3,"cm"))
          }, error = function(e) {
            ggplot() + annotate("text",x=.5,y=.5,label="Error generando ideograma") + tema_pdf
          })
          
          # Layout: ideograma arriba (ancho completo), 5 graficas abajo (grid 3+2)
          vp_upper <- grid::viewport(x = 0.01, y = 0.935, width = 0.98, height = 0.39, just = c("left","top"))
          grid::pushViewport(vp_upper)
          print(g_ideograma, newpage = FALSE)
          grid::popViewport()
          
          grob_scores  <- ggplot2::ggplotGrob(g_scores)
          grob_herencia<- ggplot2::ggplotGrob(g_herencia)
          grob_tipo    <- ggplot2::ggplotGrob(g_tipo)
          grob_rango   <- ggplot2::ggplotGrob(g_rango)
          grob_box     <- ggplot2::ggplotGrob(g_boxplot)
          
          gridExtra::grid.arrange(
            grob_scores, grob_herencia, grob_tipo,
            grob_rango,  grob_box,
            ncol = 5, nrow = 1,
            widths = c(1,1,1,1,1),
            vp = grid::viewport(x = 0.01, y = 0.515, width = 0.98, height = 0.465,
                                just = c("left","top"))
          )
          
          pie_pagina(4)
          
          # ─────────────────────────────────────────────────────────────────────
          # PAGINAS 5+ — Ficha detallada por variante prioritaria
          # ─────────────────────────────────────────────────────────────────────
          df_fichas <- df_full[
            (!is.na(df_full$score_n) & df_full$score_n >= 0) |
              (!is.na(df_full$her_b)   & df_full$her_b == "De novo") |
              df_full$en_sfari,
            , drop = FALSE
          ]
          df_fichas <- df_fichas[order(-ifelse(is.na(df_fichas$score_n),-99,df_fichas$score_n)), , drop=FALSE]
          df_fichas <- head(df_fichas, 12)
          
          for (i in seq_len(nrow(df_fichas))) {
            row <- df_fichas[i, ]
            incProgress(1 / (n_pasos_total + 1),
                        detail = paste0("Variante ", i, " of ", nrow(df_fichas), "…"))
            grid::grid.newpage()
            
            # Classification and header colour
            clasif_val <- if (!is.null(row$clasificacion)) as.character(row$clasificacion) else "--"
            cat_score  <- if (!is.na(row$score_n) && row$score_n >= 0.5)  "PATOGENICA"
            else if (!is.na(row$score_n) && row$score_n >= 0) "POTENCIALMENTE PATOGENICA"
            else if (!is.na(row$her_b) && row$her_b == "De novo") "DE NOVO"
            else if (!is.na(row$en_sfari) && row$en_sfari)        "GEN SFARI"
            else "REVISION"
            col_cat    <- switch(cat_score,
                                 "PATOGENICA"                 = "#922B21",
                                 "POTENCIALMENTE PATOGENICA"  = "#D35400",
                                 "DE NOVO"                    = "#C0392B",
                                 "GEN SFARI"                  = "#E67E22",
                                 "#2C3E50")
            her_color  <- her_pal[as.character(row$her_b)]
            if (is.na(her_color)) her_color <- "#7F8C8D"
            
            region_str <- paste0("chr", row$chr_c, ":",
                                 format(round(row$sta_n), big.mark=".", scientific=FALSE), " - ",
                                 format(round(row$end_n), big.mark=".", scientific=FALSE))
            
            genes_str <- if ("Gene_name" %in% names(row) && nzchar(as.character(row$Gene_name))
                             && as.character(row$Gene_name) != "NA") as.character(row$Gene_name) else "--"
            genes_short <- if (nchar(genes_str) > 120) paste0(substr(genes_str, 1, 120), "...") else genes_str
            
            nota_str   <- if (nchar(row$nota) > 0) row$nota else "(sin notas)"
            nota_short <- if (nchar(nota_str) > 200) paste0(substr(nota_str, 1, 200), "...") else nota_str
            
            # ── Cabecera de ficha ──────────────────────────────────────────────
            grid::grid.rect(x=0, y=1, width=1, height=0.13, just=c("left","top"),
                            gp=grid::gpar(fill=col_cat, col=NA))
            grid::grid.rect(x=0, y=1, width=1, height=0.003, just=c("left","top"),
                            gp=grid::gpar(fill="#F39C12", col=NA))
            grid::grid.rect(x=0.985, y=1, width=0.015, height=0.13, just=c("left","top"),
                            gp=grid::gpar(fill=her_color, col=NA))
            
            grid::grid.text(paste0("Variante #", i, "  |  ", cat_score),
                            x=0.025, y=0.975, just=c("left","top"),
                            gp=grid::gpar(fontsize=15, fontface="bold", col="white"))
            grid::grid.text(paste0(region_str, "  |  ", row$len_lbl,
                                   "  |  Familia: ", fam,
                                   "  |  Sexo: ", sexo_lbl),
                            x=0.025, y=0.940, just=c("left","top"),
                            gp=grid::gpar(fontsize=9, col="#FDFEFE"))
            grid::grid.text(
              paste0("Score AnnotSV: ", if(!is.na(row$score_n)) sprintf("%.4f", row$score_n) else "--",
                     "   |   Clasificacion: ", clasif_val,
                     "   |   SFARI: ", if (row$en_sfari) "Yes" else "No"),
              x=0.025, y=0.907, just=c("left","top"),
              gp=grid::gpar(fontsize=9, col="white", fontface="bold"))
            
            # ── Cuerpo en dos columnas ──────────────────────────────────────────
            # COLUMNA IZQUIERDA: campos clinicos
            y0 <- 0.845; dy <- 0.063; col1_x <- 0.025; col2_x <- 0.51
            
            field <- function(label, value, x, y, label_w = 0.22) {
              grid::grid.rect(x=x, y=y, width=0.45, height=dy*0.85, just=c("left","top"),
                              gp=grid::gpar(fill="#F8F9FA", col="#E5E7E9", lwd=0.5))
              grid::grid.text(label, x=x+0.008, y=y-0.008, just=c("left","top"),
                              gp=grid::gpar(fontsize=7.5, fontface="bold", col="#2C3E50"))
              val_display <- if (is.na(value) || value == "" || value == "NA") "--"
              else if (nchar(as.character(value)) > 75)
                paste0(substr(as.character(value),1,75),"...")
              else as.character(value)
              grid::grid.text(val_display, x=x+label_w, y=y-0.008, just=c("left","top"),
                              gp=grid::gpar(fontsize=8, col="#444444"))
            }
            
            # Campos columna izquierda
            campos_izq <- list(
              list("Region (hg38)",         region_str),
              list("Tamano",                row$len_lbl),
              list("Tipo de variante",      if("SV_type"%in%names(row)) as.character(row$SV_type) else "--"),
              list("Herencia",              if(!is.na(row$her_b)) row$her_b else "--"),
              list("Tipo de Herencia (raw)",if("Tipo_Herencia"%in%names(row)) as.character(row$Tipo_Herencia) else "--"),
              list("Rango anotacion",       if("Tipo_Rango"%in%names(row)) as.character(row$Tipo_Rango) else "--"),
              list("Score AnnotSV",         if(!is.na(row$score_n)) sprintf("%.4f", row$score_n) else "--"),
              list("Clasificacion clinica", clasif_val),
              list("Gen SFARI",             if(row$en_sfari) "SI - en panel SFARI" else "No"),
              list("Modalidad",             if("Modalidad"%in%names(row)) as.character(row$Modalidad) else "--"),
              list("Fenotipo familiar",     {
                ft <- rv$fenotipos[[fam]] %||% ""
                if (nzchar(trimws(ft)))
                  if (nchar(ft) > 78) paste0(substr(ft,1,76),"\u2026") else ft
                else "No registrado"
              })
            )
            
            for (j in seq_along(campos_izq)) {
              field(campos_izq[[j]][[1]], campos_izq[[j]][[2]],
                    x = col1_x, y = y0 - (j-1)*dy)
            }
            
            # Campos columna derecha
            campos_der <- list(
              list("Genes afectados",       genes_short),
              list("DRAGEN exact",          if("Illumina_DRAGEN.exact.counts"%in%names(row))
                as.character(row[["Illumina_DRAGEN.exact.counts"]]) else "--"),
              list("DRAGEN similares",      if("Illumina_DRAGEN.similar.counts"%in%names(row))
                as.character(row[["Illumina_DRAGEN.similar.counts"]]) else "--"),
              list("Referencia match",      if("Referencia_Match"%in%names(row))
                as.character(row$Referencia_Match) else "--"),
              list("Sim MaxLen",            if("Sim_MaxLen"%in%names(row)) as.character(row$Sim_MaxLen) else "--"),
              list("Jaccard",               if("Jaccard"%in%names(row)) as.character(row$Jaccard) else "--"),
              list("Confianza region",      if("Confianza_Region"%in%names(row)) as.character(row$Confianza_Region) else "--"),
              list("DECIPHER",              paste0("deciphergenomics.org/browser#q/",row$chr_c,":",row$sta_n,"-",row$end_n)),
              list("Franklin",              paste0("franklin.genoox.com/  chr",row$chr_c,"-",row$sta_n)),
              list("gnomAD SV (r4)",        paste0("gnomad.broadinstitute.org/region/",row$chr_c,":",row$sta_n,"-",row$end_n)),
              list("UCSC Browser",          paste0("genome.ucsc.edu/hgTracks?db=hg38&position=chr",row$chr_c,":",row$sta_n)),
              list("OMIM",                  paste0("omim.org/search?search=", strsplit(genes_str,"[,/|]")[[1]][1]))
            )
            
            for (j in seq_along(campos_der)) {
              field(campos_der[[j]][[1]], campos_der[[j]][[2]],
                    x = col2_x, y = y0 - (j-1)*dy)
            }
            
            # ── Bloque de notas ────────────────────────────────────────────────
            nota_y <- y0 - length(campos_izq) * dy - 0.015
            grid::grid.rect(x=0.025, y=nota_y, width=0.95, height=0.1,
                            just=c("left","top"),
                            gp=grid::gpar(fill="#EBF5FB", col="#2980B9", lwd=0.7))
            grid::grid.text("Notas clinicas:", x=0.035, y=nota_y-0.01,
                            just=c("left","top"),
                            gp=grid::gpar(fontsize=8, fontface="bold", col="#2980B9"))
            grid::grid.text(nota_short, x=0.035, y=nota_y-0.035,
                            just=c("left","top"), default.units="npc",
                            gp=grid::gpar(fontsize=7.5, col="#333333"))
            
            # ── Web screenshot pages ────────────────────────────────────────
            # Each variant generates 2 screenshot pages (B and C):
            #   Page B: UCSC (full width) + gnomAD | DECIPHER | Franklin
            #   Page C: ClinVar | VarSome | Ensembl  +  OMIM | URL table
            pagina_ficha <- (i - 1) * 3 + 5   # each variant: 1 clinical + 2 screenshots = 3 pages
            pie_pagina(pagina_ficha)
            
            # ── Construir todas las URLs ─────────────────────────────────────
            sty_up  <- toupper(trimws(if("SV_type"%in%names(row)) as.character(row$SV_type) else "CNV"))
            gene1   <- trimws(strsplit(genes_str, "[,/|]")[[1]][1])
            if (is.na(gene1) || gene1 == "--") gene1 <- ""
            
            url_ucsc     <- paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38",
                                   "&hideTracks=1&knownGene=pack&rmsk=dense",
                                   "&position=chr", row$chr_c, ":", row$sta_n, "-", row$end_n)
            url_gnomad   <- paste0("https://gnomad.broadinstitute.org/region/",
                                   row$chr_c, ":", row$sta_n, "-", row$end_n, "?dataset=gnomad_sv_r4")
            url_decipher <- paste0("https://www.deciphergenomics.org/browser#q/",
                                   row$chr_c, ":", row$sta_n, "-", row$end_n)
            url_franklin <- paste0("https://franklin.genoox.com/clinical-db/variant/sv/chr",
                                   row$chr_c, "-", row$sta_n, "-", row$end_n, "-", sty_up, "-hg38")
            url_clinvar  <- paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=",
                                   row$chr_c, "%5Bchr%5D+AND+",
                                   row$sta_n, "%3A", row$end_n, "%5Bchrpos37%5D")
            url_varsome  <- paste0("https://varsome.com/variant/hg38/Chr",
                                   row$chr_c, ":", row$sta_n, "-", row$end_n,
                                   ":", sty_up, "?annotation-mode=germline")
            url_ensembl  <- paste0("https://www.ensembl.org/Homo_sapiens/Location/View?r=",
                                   row$chr_c, ":", row$sta_n, "-", row$end_n)
            url_omim     <- if (nzchar(gene1))
              paste0("https://www.omim.org/search?search=", URLencode(gene1, reserved=TRUE))
            else "https://www.omim.org"
            
            # ── Helpers ───────────────────────────────────────────────────────
            # ── Localiza Edge, Chrome o Chromium (multiplataforma) ────────────────
            encontrar_edge <- function() {
              # 1. Variable de entorno explicita (maxima prioridad)
              env_path <- Sys.getenv("EDGE_PATH", unset = "")
              if (nzchar(env_path) && file.exists(env_path))
                return(normalizePath(env_path, winslash = "/"))
              
              # 2. Rutas estandar de instalacion en Windows
              local_appdata <- Sys.getenv("LOCALAPPDATA", unset = "")
              prog_files    <- Sys.getenv("PROGRAMFILES",  unset = "C:/Program Files")
              prog_files86  <- Sys.getenv("PROGRAMFILES(X86)", unset = "C:/Program Files (x86)")
              
              rutas_std <- c(
                file.path(prog_files,    "Microsoft/Edge/Application/msedge.exe"),
                file.path(prog_files86,  "Microsoft/Edge/Application/msedge.exe"),
                file.path(local_appdata, "Microsoft/Edge/Application/msedge.exe")
              )
              for (p in rutas_std) {
                if (nzchar(p) && file.exists(p))
                  return(normalizePath(p, winslash = "/"))
              }
              
              # 3. Registro de Windows (cubre instalaciones no estandar)
              if (.Platform$OS.type == "windows") {
                reg_keys <- c(
                  "HKLM\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\App Paths\\msedge.exe",
                  "HKCU\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\App Paths\\msedge.exe"
                )
                for (key in reg_keys) {
                  ruta_reg <- tryCatch(
                    readRegistry(key)[["(Default)"]],
                    error = function(e) NULL
                  )
                  if (!is.null(ruta_reg) && nzchar(ruta_reg) && file.exists(ruta_reg))
                    return(normalizePath(ruta_reg, winslash = "/"))
                }
                
                # 4. where.exe como ultimo recurso
                where_result <- tryCatch(
                  system2("where", "msedge", stdout = TRUE, stderr = FALSE),
                  error = function(e) character(0)
                )
                for (p in where_result) {
                  if (nzchar(p) && file.exists(p))
                    return(normalizePath(p, winslash = "/"))
                }
              }
              
              NULL  # Edge no encontrado
            }
            
            ultimo_error_webshot <- ""   # se sobreescribe si hay error de captura
            
            capturar_web <- function(url, w=1440, h=750, delay=4, zoom=1.5, clip=NULL) {
              if (!requireNamespace("webshot2", quietly=TRUE) ||
                  !requireNamespace("png",      quietly=TRUE)) return(NULL)
              
              # If we find an explicit browser, we pass it to chromote.
              # If not, we let chromote attempt its own auto-detection.
              edge_path  <- encontrar_edge()
              old_chrome <- Sys.getenv("CHROMOTE_CHROME", unset = NA_character_)
              if (!is.null(edge_path)) {
                Sys.setenv(CHROMOTE_CHROME = edge_path)
                on.exit({
                  if (is.na(old_chrome)) Sys.unsetenv("CHROMOTE_CHROME")
                  else                   Sys.setenv(CHROMOTE_CHROME = old_chrome)
                }, add = TRUE)
              }
              
              tmp <- tempfile(fileext = ".png")
              tryCatch({
                webshot2::webshot(url, file = tmp, vwidth = w, vheight = h,
                                  delay = delay, zoom = zoom, cliprect = clip)
                if (file.exists(tmp) && file.info(tmp)$size > 8000) {
                  img <- png::readPNG(tmp)
                  try(file.remove(tmp), silent = TRUE)
                  img
                } else {
                  try(file.remove(tmp), silent = TRUE)
                  NULL
                }
              }, error = function(e) {
                try(file.remove(tmp), silent = TRUE)
                ultimo_error_webshot <<- conditionMessage(e)
                NULL
              })
            }
            
            draw_shot <- function(img, x, y, w, h, titulo, dominio,
                                  col_tit="#2C3E50", tit_size=8, err_msg="") {
              bar_h <- 0.036
              grid::grid.rect(x=x, y=y, width=w, height=h, just=c("left","top"),
                              gp=grid::gpar(fill="#F4F6F7", col="#BBBBBB", lwd=0.6))
              grid::grid.rect(x=x, y=y, width=w, height=bar_h, just=c("left","top"),
                              gp=grid::gpar(fill=col_tit, col=NA))
              grid::grid.text(titulo, x=x+0.008, y=y-bar_h/2, just=c("left","center"),
                              gp=grid::gpar(fontsize=tit_size, fontface="bold", col="white"))
              grid::grid.text(dominio, x=x+w-0.006, y=y-bar_h/2, just=c("right","center"),
                              gp=grid::gpar(fontsize=6, col="#DDDDDD", fontface="italic"))
              img_h <- h - bar_h
              if (!is.null(img)) {
                vp <- grid::viewport(x=x, y=y-bar_h, width=w, height=img_h, just=c("left","top"))
                grid::pushViewport(vp)
                grid::grid.raster(img, width=1, height=1, interpolate=TRUE)
                grid::popViewport()
              } else {
                mx <- x + w/2; my <- y - bar_h - img_h/2
                grid::grid.rect(x=x+0.004, y=y-bar_h, width=w-0.008, height=img_h,
                                just=c("left","top"),
                                gp=grid::gpar(fill="white", col="#E8E8E8", lwd=0.4))
                grid::grid.text("\u26a0  Captura no disponible",
                                x=mx, y=my+0.04, just="center",
                                gp=grid::gpar(fontsize=8.5, col="#999999", fontface="italic"))
                motivo <- if (nzchar(trimws(err_msg))) {
                  if (nchar(err_msg) > 90) paste0(substr(err_msg, 1, 88), "\u2026") else err_msg
                } else {
                  "Abre la URL manualmente con el navegador"
                }
                grid::grid.text(motivo, x=mx, y=my+0.01, just="center",
                                gp=grid::gpar(fontsize=6.5, col="#BBBBBB"))
                grid::grid.text(paste0("\u2192 ", dominio), x=mx, y=my-0.02, just="center",
                                gp=grid::gpar(fontsize=6.5, col="#AAAAAA", fontface="italic"))
              }
            }
            
            cab_capturas <- function(subtitulo) {
              grid::grid.rect(x=0, y=1, width=1, height=0.055, just=c("left","top"),
                              gp=grid::gpar(fill=col_cat, col=NA))
              grid::grid.rect(x=0, y=1, width=1, height=0.003, just=c("left","top"),
                              gp=grid::gpar(fill="#F39C12", col=NA))
              grid::grid.text(
                paste0("Capturas web  |  Variante #", i, "  |  ", region_str, "  |  ", row$len_lbl),
                x=0.025, y=0.967, just=c("left","top"),
                gp=grid::gpar(fontsize=10, fontface="bold", col="white"))
              grid::grid.text(
                paste0(subtitulo, "  |  Score: ",
                       if(!is.na(row$score_n)) sprintf("%.4f",row$score_n) else "--",
                       "  |  Familia: ", fam, "  |  ", format(Sys.time(),"%d/%m/%Y %H:%M")),
                x=0.025, y=0.934, just=c("left","top"),
                gp=grid::gpar(fontsize=7.5, col="white"))
            }
            
            has_webshot <- requireNamespace("webshot2", quietly=TRUE) &&
              requireNamespace("png", quietly=TRUE)
            
            aviso_no_ws <- function() {
              if (!has_webshot) {
                grid::grid.rect(x=0.015, y=0.882, width=0.97, height=0.038,
                                just=c("left","top"),
                                gp=grid::gpar(fill="#FFF3CD", col="#F0AD4E", lwd=0.7))
                grid::grid.text(
                  "Para capturas reales: install.packages('webshot2'); install.packages('png')  |  Edge se detecta automaticamente (rutas estandar + registro Windows). Ruta manual: variable EDGE_PATH",
                  x=0.5, y=0.863, just="center",
                  gp=grid::gpar(fontsize=7, col="#856404"))
              }
            }
            
            pw3  <- 0.314; gap3 <- 0.013   # 1/3 panel width and gap
            pw2l <- 0.47;  pw2r <- 0.497   # anchos panel 50/50
            
            # ════════════════════════════════════════════════════════════════
            # PAGE B: UCSC (full width) + gnomAD | DECIPHER | Franklin
            # ════════════════════════════════════════════════════════════════
            grid::grid.newpage()
            cab_capturas("1/2 — UCSC · gnomAD · DECIPHER · Franklin")
            aviso_no_ws()
            
            img_ucsc     <- if(has_webshot) capturar_web(url_ucsc,
                                                         w=1440, h=620, delay=4, zoom=1.8, clip=c(0,55,1440,620)) else NULL
            img_gnomad   <- if(has_webshot) capturar_web(url_gnomad,
                                                         w=1440, h=700, delay=5, zoom=1.5, clip=c(0,70,1440,700)) else NULL
            img_decipher <- if(has_webshot) capturar_web(url_decipher,
                                                         w=1440, h=700, delay=5, zoom=1.4, clip=c(0,55,1440,700)) else NULL
            img_franklin <- if(has_webshot) capturar_web(url_franklin,
                                                         w=1440, h=700, delay=6, zoom=1.4, clip=c(0,55,1440,700)) else NULL
            
            y_top_b  <- if(has_webshot) 0.885 else 0.84
            y_row2_b <- y_top_b - 0.385
            
            draw_shot(img_ucsc,
                      x=0.015, y=y_top_b, w=0.97, h=0.375,
                      titulo   = paste0("UCSC Genome Browser — hg38 | chr", row$chr_c, ":",
                                        format(round(row$sta_n),big.mark=".",scientific=FALSE),
                                        "\u2013", format(round(row$end_n),big.mark=".",scientific=FALSE),
                                        "  (", row$len_lbl, ")"),
                      dominio  = "genome.ucsc.edu",
                      col_tit  = "#1A4A7A", tit_size = 8.5, err_msg = ultimo_error_webshot)
            
            draw_shot(img_gnomad,
                      x=0.015, y=y_row2_b, w=pw3, h=0.35,
                      titulo   = paste0("gnomAD SV r4 — chr", row$chr_c, ":", row$sta_n, "\u2013", row$end_n),
                      dominio  = "gnomad.broadinstitute.org", col_tit = "#1A6B3A", err_msg = ultimo_error_webshot)
            draw_shot(img_decipher,
                      x=0.015+pw3+gap3, y=y_row2_b, w=pw3, h=0.35,
                      titulo   = paste0("DECIPHER — chr", row$chr_c, ":", row$sta_n, "\u2013", row$end_n),
                      dominio  = "deciphergenomics.org", col_tit = "#2D6A27", err_msg = ultimo_error_webshot)
            draw_shot(img_franklin,
                      x=0.015+2*(pw3+gap3), y=y_row2_b, w=pw3, h=0.35,
                      titulo   = paste0("Franklin Genoox — chr", row$chr_c, "\u2013", row$sta_n, " (", sty_up, ")"),
                      dominio  = "franklin.genoox.com", col_tit = "#6B2D6B", err_msg = ultimo_error_webshot)
            
            pie_pagina(pagina_ficha + 1)
            
            # ════════════════════════════════════════════════════════════════
            # PAGE C: ClinVar | VarSome | Ensembl  +  OMIM | URL table
            # ════════════════════════════════════════════════════════════════
            grid::grid.newpage()
            cab_capturas("2/2 — ClinVar · VarSome · Ensembl · OMIM")
            aviso_no_ws()
            
            img_clinvar  <- if(has_webshot) capturar_web(url_clinvar,
                                                         w=1440, h=700, delay=5, zoom=1.4, clip=c(0,60,1440,700)) else NULL
            img_varsome  <- if(has_webshot) capturar_web(url_varsome,
                                                         w=1440, h=700, delay=6, zoom=1.4, clip=c(0,60,1440,700)) else NULL
            img_ensembl  <- if(has_webshot) capturar_web(url_ensembl,
                                                         w=1440, h=700, delay=5, zoom=1.5, clip=c(0,60,1440,700)) else NULL
            img_omim     <- if(has_webshot && nzchar(gene1))
              capturar_web(url_omim, w=1440, h=700, delay=5, zoom=1.4,
                           clip=c(0,55,1440,700)) else NULL
            
            y_top_c  <- if(has_webshot) 0.885 else 0.84
            y_row2_c <- y_top_c - 0.385
            
            draw_shot(img_clinvar,
                      x=0.015, y=y_top_c, w=pw3, h=0.375,
                      titulo   = paste0("ClinVar — chr", row$chr_c, ":", row$sta_n, "\u2013", row$end_n),
                      dominio  = "ncbi.nlm.nih.gov/clinvar", col_tit = "#993333", err_msg = ultimo_error_webshot)
            draw_shot(img_varsome,
                      x=0.015+pw3+gap3, y=y_top_c, w=pw3, h=0.375,
                      titulo   = paste0("VarSome — chr", row$chr_c, ":", row$sta_n, "\u2013", row$end_n, " (", sty_up, ")"),
                      dominio  = "varsome.com", col_tit = "#4A6FA5", err_msg = ultimo_error_webshot)
            draw_shot(img_ensembl,
                      x=0.015+2*(pw3+gap3), y=y_top_c, w=pw3, h=0.375,
                      titulo   = paste0("Ensembl — chr", row$chr_c, ":", row$sta_n, "\u2013", row$end_n),
                      dominio  = "ensembl.org", col_tit = "#8B4500", err_msg = ultimo_error_webshot)
            
            # OMIM (izquierda)
            draw_shot(img_omim,
                      x=0.015, y=y_row2_c, w=pw2l, h=0.355,
                      titulo   = paste0("OMIM — Gen: ", if(nzchar(gene1)) gene1 else "(sin gen anotado)"),
                      dominio  = "omim.org", col_tit = "#7B2D8B", err_msg = ultimo_error_webshot)
            
            # Panel tabla de URLs (derecha)
            tx <- 0.015 + pw2l + gap3; ty <- y_row2_c; tw <- pw2r; th <- 0.355
            grid::grid.rect(x=tx, y=ty, width=tw, height=th, just=c("left","top"),
                            gp=grid::gpar(fill="#F8F9FA", col="#BBBBBB", lwd=0.6))
            grid::grid.rect(x=tx, y=ty, width=tw, height=0.036, just=c("left","top"),
                            gp=grid::gpar(fill="#2C3E50", col=NA))
            grid::grid.text("URLs completas para consulta manual",
                            x=tx+0.008, y=ty-0.018, just=c("left","center"),
                            gp=grid::gpar(fontsize=8, fontface="bold", col="white"))
            
            urls_ref <- list(
              list("UCSC",          "#1A4A7A", url_ucsc),
              list("gnomAD SV r4",  "#1A6B3A", url_gnomad),
              list("DECIPHER",      "#2D6A27", url_decipher),
              list("Franklin",      "#6B2D6B", url_franklin),
              list("ClinVar",       "#993333", url_clinvar),
              list("VarSome",       "#4A6FA5", url_varsome),
              list("Ensembl",       "#8B4500", url_ensembl),
              list("OMIM",          "#7B2D8B", url_omim)
            )
            url_dy <- (th - 0.036 - 0.006) / length(urls_ref)
            for (ui in seq_along(urls_ref)) {
              ur  <- urls_ref[[ui]]
              uy  <- ty - 0.036 - 0.003 - (ui-1) * url_dy
              bg_u <- if(ui %% 2 == 0) "#EAECEE" else "white"
              grid::grid.rect(x=tx+0.003, y=uy, width=tw-0.006, height=url_dy*0.93,
                              just=c("left","top"), gp=grid::gpar(fill=bg_u, col=NA))
              grid::grid.text(ur[[1]], x=tx+0.01, y=uy-url_dy*0.46,
                              just=c("left","center"),
                              gp=grid::gpar(fontsize=7, fontface="bold", col=ur[[2]]))
              url_s <- if(nchar(ur[[3]])>70) paste0(substr(ur[[3]],1,68),"\u2026") else ur[[3]]
              grid::grid.text(url_s, x=tx+0.075, y=uy-url_dy*0.46,
                              just=c("left","center"),
                              gp=grid::gpar(fontsize=5.8, col="#555555"))
            }
            
            pie_pagina(pagina_ficha + 2)
          }
          
          incProgress(1, detail = "Report complete!")
          invisible(NULL)
          
        }, error = function(e) {
          tryCatch({
            grid::grid.newpage()
            grid::grid.text(paste0("Error al generar el informe:\n", conditionMessage(e)),
                            x = 0.5, y = 0.5, gp = grid::gpar(fontsize = 11, col = "#C0392B"))
          }, error = function(e2) NULL)
          showNotification(paste0("Error en informe PDF: ", conditionMessage(e)), type = "error", duration = 8)
        })
      })  # withProgress
    }
  )
  
  # =============================================================================
  # GRAPHS (formerly Statistics) — server logic
  # =============================================================================
  
  est_datos <- reactive({
    input$btn_est_refresh
    df <- datos_filtrados()
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    if ("Annotation_mode" %in% names(df))
      df <- df[!is.na(df$Annotation_mode) & df$Annotation_mode == "full", , drop = FALSE]
    if (nrow(df) == 0) return(data.frame())
    df$her_limpia <- her_base(df$Tipo_Herencia)
    df$score_num  <- suppressWarnings(as.numeric(df$AnnotSV_ranking_score))
    df$start_num  <- suppressWarnings(as.numeric(df$SV_start))
    df$end_num    <- suppressWarnings(as.numeric(df$SV_end))
    df$tamano_kb  <- pmax(0, df$end_num - df$start_num + 1) / 1000
    df$chr_clean  <- toupper(trimws(gsub("(?i)^chr", "", as.character(df$SV_chrom), perl = TRUE)))
    sv_up         <- toupper(as.character(df$SV_type))
    df$tipo_agrup <- ifelse(grepl("DEL|LOSS", sv_up), "DEL",
                            ifelse(grepl("DUP|GAIN", sv_up), "DUP", "OTHER"))
    if ("ID_Familia" %in% names(df)) {
      df$sexo <- sapply(as.character(df$ID_Familia), function(fam) {
        s <- rv$sexos[[fam]]
        if (is.null(s) || !nzchar(s)) "ND" else s
      })
    } else {
      df$sexo <- "ND"
    }
    df
  })
  
  
  output$ui_est_resumen_filtro <- renderUI({
    df    <- datos_filtrados()
    n_tot <- if (!is.null(df)) nrow(df) else 0
    n_full <- if (!is.null(df) && "Annotation_mode" %in% names(df))
      sum(!is.na(df$Annotation_mode) & df$Annotation_mode == "full")
    else n_tot
    n_fam <- if (!is.null(df) && "ID_Familia" %in% names(df))
      length(unique(df$ID_Familia)) else 0
    div(class = "small",
        div(class = "d-flex justify-content-between py-1",
            tags$span(class = "text-muted", "Variantes (full):"), tags$b(n_full)),
        div(class = "d-flex justify-content-between py-1",
            tags$span(class = "text-muted", "Total filas:"),      tags$b(n_tot)),
        div(class = "d-flex justify-content-between py-1",
            tags$span(class = "text-muted", "Familias:"),         tags$b(n_fam))
    )
  })
  
  observeEvent(input$btn_ir_resultados, {
    updateNavbarPage(session, "main_navbar", selected = "tab_results")
  })
  
  output$ui_est_kpis <- renderUI({
    df <- est_datos()
    if (nrow(df) == 0)
      return(div(class="alert alert-warning", "\u26a0 Sin datos. Carga resultados primero."))
    n_var    <- nrow(df)
    n_fam    <- length(unique(as.character(df$ID_Familia)))
    n_denovo <- sum(!is.na(df$her_limpia) & df$her_limpia == "De novo")
    pct_dn   <- if (n_var > 0) round(n_denovo/n_var*100,1) else 0
    n_sfari  <- if ("En_Panel_Genes" %in% names(df))
      sum(es_sfari_col(df$En_Panel_Genes)) else NA
    pct_sfar <- if (!is.na(n_sfari) && n_var>0) round(n_sfari/n_var*100,1) else NA
    sc_med   <- round(mean(df$score_num,na.rm=TRUE),3)
    sc_med_s <- if (is.nan(sc_med)) "\u2014" else as.character(sc_med)
    tam_med  <- round(median(df$tamano_kb,na.rm=TRUE),1)
    # Sexo
    n_varones <- if ("sexo" %in% names(df)) length(unique(df$ID_Familia[df$sexo=="M"])) else 0
    n_mujeres <- if ("sexo" %in% names(df)) length(unique(df$ID_Familia[df$sexo=="F"])) else 0
    n_sex_nd  <- n_fam - n_varones - n_mujeres
    kpi <- function(icono, titulo, valor, color="#2C6FAC", subtitulo=NULL) {
      div(class="card h-100", style=paste0("border-left:5px solid ",color,"; border-radius:8px;"),
          div(class="card-body p-3",
              div(class="d-flex justify-content-between align-items-start",
                  div(div(class="text-muted small fw-bold text-uppercase mb-1",
                          style="font-size:0.72em; letter-spacing:0.05em;", titulo),
                      div(style=paste0("font-size:1.9em; font-weight:700; color:",color,";"), valor)),
                  div(style=paste0("font-size:2em; color:",color,"; opacity:0.2;"), icono)),
              if (!is.null(subtitulo)) div(class="text-muted mt-1",style="font-size:0.78em;",subtitulo)))
    }
    tagList(
      layout_columns(col_widths=rep(2,6),
                     kpi("\U0001f9ec","Variantes",     format(n_var,big.mark="."),"#2C6FAC"),
                     kpi("\U0001f46a","Familias",      format(n_fam,big.mark="."),"#27AE60"),
                     kpi("\u26a1","De novo",           paste0(n_denovo," (",pct_dn,"%)"), "#CC0000","del total"),
                     kpi("\u2b50","Con gen SFARI",
                         if(is.na(pct_sfar)) "\u2014" else paste0(n_sfari," (",pct_sfar,"%)"),
                         "#E67E22","in gene panel"),
                     kpi("\U0001f4ca","Score medio",   sc_med_s,"#8E44AD","AnnotSV score"),
                     kpi("\U0001f4cf","Tam. mediano",  paste0(tam_med," Kb"),"#16A085","variante")
      ),
      layout_columns(col_widths=c(4,4,4),
                     kpi("\u2642","Varones (familias)", format(n_varones,big.mark="."), "#2980B9",
                         if(n_fam>0) paste0(round(n_varones/n_fam*100,1),"% de familias") else NULL),
                     kpi("\u2640","Mujeres (familias)", format(n_mujeres,big.mark="."), "#C0392B",
                         if(n_fam>0) paste0(round(n_mujeres/n_fam*100,1),"% de familias") else NULL),
                     kpi("\u2753","Sexo desconocido",   format(n_sex_nd,big.mark="."),  "#7F8C8D",
                         "familias sin sexo registrado")
      )
    )
  })
  
  output$est_plot_scores <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0)
    scores <- df$score_num[is.finite(df$score_num)]; req(length(scores)>0)
    plot_ly(x=~scores, type="histogram", nbinsx=40,
            marker=list(color="#2C6FAC",line=list(color="white",width=0.5)),
            hovertemplate="Score:%{x:.3f}<br>N:%{y}<extra></extra>") |>
      layout(xaxis=list(title="AnnotSV ranking score"), yaxis=list(title="N variantes"),
             shapes=list(
               list(type="line",x0=0,x1=0,y0=0,y1=1,yref="paper",line=list(color="orange",dash="dash",width=1.5)),
               list(type="line",x0=0.5,x1=0.5,y0=0,y1=1,yref="paper",line=list(color="red",dash="dash",width=1.5))),
             margin=list(t=20,b=50,l=50,r=10), plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_herencia_pie <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0,"her_limpia" %in% names(df))
    tab <- sort(table(df$her_limpia[!is.na(df$her_limpia)]),decreasing=TRUE)
    cols <- HER_COLORS[names(tab)]; cols[is.na(cols)] <- "#AAAAAA"
    plot_ly(labels=names(tab),values=as.numeric(tab),type="pie",hole=0.4,
            marker=list(colors=unname(cols),line=list(color="#FFFFFF",width=1.5)),
            textinfo="label+percent",
            hovertemplate="<b>%{label}</b><br>%{value} (%{percent})<extra></extra>") |>
      layout(showlegend=FALSE,margin=list(t=10,b=10,l=10,r=10),paper_bgcolor="#FAFAFA")
  })
  
  output$est_plot_carga_familia <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0)
    her_levs <- c("De novo","Paternal","Maternal","Combined",
                  "Paternal (Probable)","Maternal (Probable)","Combined (Probable)","Unknown")
    df$her_limpia[is.na(df$her_limpia)] <- "Unknown"
    df$her_f <- factor(df$her_limpia,levels=c(her_levs,setdiff(unique(df$her_limpia),her_levs)))
    conteo <- df |> group_by(ID_Familia,her_f) |> summarise(n=n(),.groups="drop")
    orden_fam <- conteo |> group_by(ID_Familia) |> summarise(total=sum(n),.groups="drop") |>
      arrange(desc(total)) |> pull(ID_Familia)
    # Add ♂/♀ symbol to X-axis label
    fam_sex <- if ("sexo" %in% names(df))
      setNames(df$sexo[!duplicated(df$ID_Familia)], df$ID_Familia[!duplicated(df$ID_Familia)])
    else setNames(rep("ND", length(orden_fam)), orden_fam)
    x_labels <- sapply(orden_fam, function(f) {
      s <- fam_sex[f]
      if (is.na(s)||s=="ND") f
      else paste0(f, if(s=="M") " \u2642" else " \u2640")
    })
    fig <- plot_ly()
    for (h in levels(conteo$her_f)) {
      sub <- conteo[conteo$her_f==h,]; if(nrow(sub)==0) next
      col <- HER_COLORS[h]; if(is.na(col)) col <- "#AAAAAA"
      sub_x <- x_labels[match(sub$ID_Familia, orden_fam)]
      fig <- add_trace(fig,type="bar",x=sub_x,y=sub$n,name=h,
                       marker=list(color=col),
                       hovertemplate=paste0("<b>",h,"</b><br>%{x}<br>%{y}<extra></extra>"))
    }
    layout(fig,barmode="stack",
           xaxis=list(title="",categoryorder="array",categoryarray=unname(x_labels),
                      tickangle=-45,tickfont=list(size=9)),
           yaxis=list(title="N variantes"),
           legend=list(orientation="h",y=-0.35,x=0,font=list(size=10)),
           margin=list(b=110,l=50,r=10,t=10),plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_score_her <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0)
    df2 <- df[!is.na(df$score_num) & !is.na(df$her_limpia),]; req(nrow(df2)>0)
    her_ord <- names(sort(tapply(df2$score_num,df2$her_limpia,median,na.rm=TRUE),decreasing=TRUE))
    fig <- plot_ly()
    for (h in her_ord) {
      sub <- df2[df2$her_limpia==h,]; if(nrow(sub)==0) next
      col <- HER_COLORS[h]; if(is.na(col)) col <- "#AAAAAA"
      fig <- add_trace(fig,type="box",y=sub$score_num,name=h,
                       marker=list(color=col,size=3),line=list(color=col),
                       fillcolor=paste0(col,"55"),
                       hovertemplate=paste0("<b>",h,"</b><br>%{y:.3f}<extra></extra>"))
    }
    layout(fig,yaxis=list(title="AnnotSV score"),xaxis=list(title="",tickangle=-30,tickfont=list(size=9)),
           showlegend=FALSE,margin=list(b=80,l=60,r=10,t=10),plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_por_chr <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0,"chr_clean" %in% names(df))
    df2 <- df[df$chr_clean %in% CHR_ORDER,]; req(nrow(df2)>0)
    df2$chr_f <- factor(df2$chr_clean,levels=CHR_ORDER)
    conteo <- df2 |> group_by(chr_f,tipo_agrup) |> summarise(n=n(),.groups="drop")
    fig <- plot_ly()
    for (tp in c("DEL","DUP","OTHER")) {
      sub <- conteo[conteo$tipo_agrup==tp,]; if(nrow(sub)==0) next
      col <- switch(tp,DEL="#CC0000",DUP="#0066CC","#888888")
      fig <- add_trace(fig,type="bar",x=sub$chr_f,y=sub$n,name=tp,marker=list(color=col),
                       hovertemplate=paste0("<b>",tp,"</b><br>chr%{x}<br>%{y}<extra></extra>"))
    }
    layout(fig,barmode="stack",
           xaxis=list(title="Cromosoma",categoryorder="array",categoryarray=CHR_ORDER,tickfont=list(size=10)),
           yaxis=list(title="N variantes"),legend=list(orientation="h",y=-0.25,x=0),
           margin=list(b=60,l=50,r=10,t=10),plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_heatmap_chr <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0,"chr_clean" %in% names(df))
    df2 <- df[df$chr_clean %in% CHR_ORDER,]; req(nrow(df2)>0)
    conteo <- df2 |> group_by(ID_Familia,chr_clean) |> summarise(n=n(),.groups="drop")
    fams <- sort(unique(conteo$ID_Familia))
    chrs <- CHR_ORDER[CHR_ORDER %in% unique(conteo$chr_clean)]
    mat  <- matrix(0,nrow=length(fams),ncol=length(chrs),dimnames=list(fams,chrs))
    for (k in seq_len(nrow(conteo))) mat[conteo$ID_Familia[k],conteo$chr_clean[k]] <- conteo$n[k]
    plot_ly(z=mat,x=paste0("chr",chrs),y=fams,type="heatmap",
            colorscale=list(c(0,"#FFFFFF"),c(0.01,"#D6EAF8"),c(0.5,"#2C6FAC"),c(1,"#1A2980")),
            hovertemplate="Familia:<b>%{y}</b><br>%{x}<br>N:<b>%{z}</b><extra></extra>") |>
      layout(xaxis=list(title="",tickfont=list(size=9),tickangle=-45),
             yaxis=list(title="",tickfont=list(size=9),autorange="reversed"),
             margin=list(b=80,l=90,r=20,t=10),paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_tamano_hist <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0)
    df2 <- df[!is.na(df$tamano_kb) & df$tamano_kb>0,]; req(nrow(df2)>0)
    fig <- plot_ly()
    for (tp in c("DEL","DUP")) {
      sub <- df2[df2$tipo_agrup==tp,]; if(nrow(sub)==0) next
      col <- if(tp=="DEL") "#CC0000" else "#0066CC"
      fig <- add_trace(fig,type="histogram",x=log10(sub$tamano_kb),name=tp,
                       opacity=0.75,marker=list(color=col,line=list(color="#FFFFFF",width=0.5)),nbinsx=30,
                       hovertemplate=paste0("<b>",tp,"</b><br>log10(Kb):%{x:.2f}<br>N:%{y}<extra></extra>"))
    }
    layout(fig,barmode="overlay",
           xaxis=list(title="Size (log10 Kb)",tickvals=c(-2,-1,0,1,2,3),
                      ticktext=c("10pb","100pb","1Kb","10Kb","100Kb","1Mb")),
           yaxis=list(title="N variantes"),legend=list(orientation="h",y=-0.2,x=0),
           margin=list(b=60,l=50,r=10,t=10),plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_tamano_her <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0)
    df2 <- df[!is.na(df$tamano_kb) & df$tamano_kb>0 & !is.na(df$her_limpia),]; req(nrow(df2)>0)
    her_ord <- names(sort(tapply(df2$tamano_kb,df2$her_limpia,median,na.rm=TRUE),decreasing=TRUE))
    fig <- plot_ly()
    for (h in her_ord) {
      sub <- df2[df2$her_limpia==h,]; if(nrow(sub)==0) next
      col <- HER_COLORS[h]; if(is.na(col)) col <- "#AAAAAA"
      fig <- add_trace(fig,type="box",y=sub$tamano_kb,name=h,
                       marker=list(color=col,size=3),line=list(color=col),fillcolor=paste0(col,"55"),
                       hovertemplate=paste0("<b>",h,"</b><br>Kb:%{y:.1f}<extra></extra>"))
    }
    layout(fig,yaxis=list(title="Size (Kb)",type="log"),
           xaxis=list(title="",tickangle=-30,tickfont=list(size=9)),
           showlegend=FALSE,margin=list(b=80,l=60,r=10,t=10),plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_score_tamano <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0)
    df2 <- df[!is.na(df$score_num) & !is.na(df$tamano_kb) & df$tamano_kb>0 & !is.na(df$her_limpia),]
    req(nrow(df2)>0)
    sexo_symbol <- c("M"="triangle-up","F"="circle","ND"="square")
    df2$sym <- sexo_symbol[df2$sexo]; df2$sym[is.na(df2$sym)] <- "square"
    fig <- plot_ly()
    for (h in unique(df2$her_limpia)) {
      sub <- df2[df2$her_limpia==h,]
      col <- HER_COLORS[h]; if(is.na(col)) col <- "#AAAAAA"
      fig <- add_trace(fig,type="scatter",mode="markers",x=log10(sub$tamano_kb),y=sub$score_num,name=h,
                       marker=list(color=paste0(col,"99"),size=7,symbol=sub$sym,
                                   line=list(color=col,width=0.5)),
                       text=paste0("Familia:",sub$ID_Familia,"<br>",h,"<br>",
                                   round(sub$tamano_kb,1)," Kb \u00b7 Score:",round(sub$score_num,3),
                                   "<br>Sexo:",ifelse(sub$sexo=="M","\u2642 Var\u00f3n",
                                                      ifelse(sub$sexo=="F","\u2640 Mujer","ND"))),
                       hovertemplate="%{text}<extra></extra>")
    }
    layout(fig,
           xaxis=list(title="Tama\u00f1o (log10 Kb)",tickvals=c(-2,-1,0,1,2,3),
                      ticktext=c("10pb","100pb","1Kb","10Kb","100Kb","1Mb")),
           yaxis=list(title="AnnotSV score"),
           legend=list(orientation="h",y=-0.25,x=0,font=list(size=10)),
           annotations=list(list(x=1,y=1.05,xref="paper",yref="paper",showarrow=FALSE,align="right",
                                 text="\u25b2 Var\u00f3n &nbsp; \u25cf Mujer &nbsp; \u25a0 ND",
                                 font=list(size=10,color="#555"))),
           margin=list(b=80,l=60,r=10,t=30),plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_top_genes <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0,"Gene_name" %in% names(df))
    gen_col <- as.character(df$Gene_name)
    gen_col <- gen_col[!is.na(gen_col) & nzchar(gen_col) & gen_col!="."]
    req(length(gen_col)>0)
    genes_vec <- trimws(unlist(strsplit(gen_col,"[/|,]")))
    genes_vec <- genes_vec[nzchar(genes_vec) & genes_vec!="."]
    tab <- sort(table(genes_vec),decreasing=TRUE); top <- head(tab,25)
    col_vec <- if ("En_Panel_Genes" %in% names(df)) {
      sfari_g <- toupper(unique(trimws(unlist(strsplit(
        as.character(df$Gene_name[es_sfari_col(df$En_Panel_Genes)]),"[/|,]")))))
      ifelse(toupper(names(top)) %in% sfari_g,"#E74C3C","#2C6FAC")
    } else rep("#2C6FAC",length(top))
    plot_ly(x=as.numeric(top),y=factor(names(top),levels=rev(names(top))),type="bar",orientation="h",
            marker=list(color=rev(col_vec),line=list(color="#FFFFFF",width=0.5)),
            hovertemplate="<b>%{y}</b><br>%{x} variantes<extra></extra>") |>
      layout(xaxis=list(title="N variantes"),yaxis=list(title="",tickfont=list(size=10)),
             margin=list(l=100,b=50,r=20,t=35),
             annotations=list(list(x=1,y=1.06,xref="paper",yref="paper",showarrow=FALSE,align="right",
                                   text="<span style='color:#E74C3C'>\u25a0</span> SFARI &nbsp; <span style='color:#2C6FAC'>\u25a0</span> No SFARI",
                                   font=list(size=10))),
             plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  output$est_tabla_genes <- renderDT({
    dg <- genes_agregados()
    if (is.null(dg) || nrow(dg)==0)
      return(datatable(data.frame(Info="Sin datos"),options=list(dom="t"),rownames=FALSE))
    # filtro SFARI ya aplicado por datos_filtrados()
    req(nrow(dg)>0)
    df_show <- head(dg[,c("gen_solo","chr_clean","N_familias","N_variantes","Herencias","Score_medio","En_SFARI")],30)
    names(df_show) <- c("Gen","Chr","Familias","Variantes","Herencias","Score medio","SFARI")
    datatable(df_show,rownames=FALSE,
              options=list(dom="ftp",pageLength=15,scrollX=TRUE,order=list(list(3,"desc"))),
              class="table table-sm table-hover") |>
      formatStyle("SFARI",backgroundColor=styleEqual(c(TRUE,FALSE),c("#FDECEA","#FFFFFF")),
                  fontWeight=styleEqual(TRUE,"bold"))
  },server=TRUE)
  
  output$est_plot_rango <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0,"Tipo_Rango" %in% names(df))
    tab <- table(df$Tipo_Rango[!is.na(df$Tipo_Rango)])
    cols <- RANGO_COLORS[names(tab)]; cols[is.na(cols)] <- "#AAAAAA"
    plot_ly(labels=names(tab),values=as.numeric(tab),type="pie",
            marker=list(colors=unname(cols),line=list(color="#FFFFFF",width=1.5)),
            textinfo="label+percent",hovertemplate="<b>%{label}</b><br>%{value}<extra></extra>") |>
      layout(showlegend=FALSE,margin=list(t=10,b=10,l=10,r=10),paper_bgcolor="#FAFAFA")
  })
  
  output$est_plot_tipo_sv <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0)
    tab <- sort(table(df$tipo_agrup),decreasing=TRUE)
    cols <- c(DEL="#CC0000",DUP="#0066CC",OTRO="#888888")
    col_vec <- cols[names(tab)]; col_vec[is.na(col_vec)] <- "#AAAAAA"
    plot_ly(labels=names(tab),values=as.numeric(tab),type="pie",
            marker=list(colors=unname(col_vec),line=list(color="#FFFFFF",width=1.5)),
            textinfo="label+percent",hovertemplate="<b>%{label}</b><br>%{value}<extra></extra>") |>
      layout(showlegend=FALSE,margin=list(t=10,b=10,l=10,r=10),paper_bgcolor="#FAFAFA")
  })
  
  output$est_plot_score_familia <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0)
    df2 <- df[!is.na(df$score_num),]; req(nrow(df2)>0)
    sc_fam <- df2 |> group_by(ID_Familia) |>
      summarise(score_med=mean(score_num,na.rm=TRUE),n_var=n(),.groups="drop") |>
      arrange(desc(score_med))
    plot_ly(x=sc_fam$score_med,
            y=factor(sc_fam$ID_Familia,levels=rev(sc_fam$ID_Familia)),
            type="bar",orientation="h",
            marker=list(
              color=scales::col_numeric(c("#D4EDDA","#FFF3CD","#FFC7CE"),
                                        domain=range(sc_fam$score_med))(sc_fam$score_med),
              line=list(color="#AAAAAA",width=0.4)),
            text=paste0(round(sc_fam$score_med,3)," (n=",sc_fam$n_var,")"),
            textposition="outside",
            hovertemplate="<b>%{y}</b><br>Score medio:%{x:.3f}<extra></extra>") |>
      layout(xaxis=list(title="Score medio AnnotSV"),yaxis=list(title="",tickfont=list(size=9)),
             margin=list(l=90,b=50,r=70,t=10),plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
  })
  
  # ── Plots stratified by sex ─────────────────────────────────────
  SEX_COLORS <- c("M" = "#2980B9", "F" = "#C0392B", "ND" = "#95A5A6")
  SEX_LABELS <- c("M" = "Var\u00f3n \u2642", "F" = "Mujer \u2640", "ND" = "Desconocido")
  
  output$est_plot_sexo_pie <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0, "sexo" %in% names(df))
    fam_sex <- df[!duplicated(df$ID_Familia), c("ID_Familia","sexo")]
    tab <- table(fam_sex$sexo)
    etiq <- SEX_LABELS[names(tab)]; etiq[is.na(etiq)] <- names(tab)[is.na(etiq)]
    cols <- SEX_COLORS[names(tab)]; cols[is.na(cols)] <- "#AAAAAA"
    plot_ly(labels=unname(etiq), values=as.numeric(tab), type="pie", hole=0.4,
            marker=list(colors=unname(cols), line=list(color="#FFFFFF",width=1.5)),
            textinfo="label+percent",
            hovertemplate="<b>%{label}</b><br>%{value} familias (%{percent})<extra></extra>") |>
      layout(showlegend=FALSE, margin=list(t=10,b=10,l=10,r=10), paper_bgcolor="#FAFAFA")
  })
  
  output$est_plot_carga_sexo <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0, "sexo" %in% names(df))
    her_levs <- c("De novo","Paternal","Maternal","Combined",
                  "Paternal (Probable)","Maternal (Probable)","Combined (Probable)","Unknown")
    df$her_f <- factor(her_base(df$Tipo_Herencia),
                       levels=c(her_levs, setdiff(unique(her_base(df$Tipo_Herencia)), her_levs)))
    df$sexo_lbl <- ifelse(df$sexo=="M","Var\u00f3n \u2642",
                          ifelse(df$sexo=="F","Mujer \u2640","Desconocido"))
    conteo <- df |> group_by(sexo_lbl, her_f) |> summarise(n=n(), .groups="drop")
    fig <- plot_ly()
    for (h in her_levs) {
      sub <- conteo[as.character(conteo$her_f)==h,]; if(nrow(sub)==0) next
      col <- HER_COLORS[h]; if(is.na(col)) col <- "#AAAAAA"
      fig <- add_trace(fig, type="bar", x=sub$sexo_lbl, y=sub$n, name=h,
                       marker=list(color=col),
                       hovertemplate=paste0("<b>",h,"</b><br>%{x}<br>%{y}<extra></extra>"))
    }
    layout(fig, barmode="stack",
           xaxis=list(title=""),
           yaxis=list(title="N variantes"),
           legend=list(orientation="h",y=-0.35,x=0,font=list(size=10)),
           margin=list(b=110,l=50,r=10,t=10), plot_bgcolor="#FAFAFA", paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_score_sexo <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0, "sexo" %in% names(df))
    df2 <- df[!is.na(df$score_num) & df$sexo != "ND", ]
    req(nrow(df2)>0)
    fig <- plot_ly()
    for (s in c("M","F")) {
      sub <- df2[df2$sexo==s,]; if(nrow(sub)==0) next
      col <- SEX_COLORS[s]; lbl <- SEX_LABELS[s]
      fig <- add_trace(fig, type="box", y=sub$score_num, name=lbl,
                       marker=list(color=col,size=4), line=list(color=col),
                       fillcolor=paste0(col,"55"),
                       hovertemplate=paste0("<b>",lbl,"</b><br>Score:%{y:.3f}<extra></extra>"))
    }
    p_val <- tryCatch({
      m <- df2$score_num[df2$sexo=="M"]; f <- df2$score_num[df2$sexo=="F"]
      if(length(m)>1 && length(f)>1) wilcox.test(m,f)$p.value else NA
    }, error=function(e) NA)
    titulo <- if(!is.na(p_val)) paste0("Wilcoxon p=",formatC(p_val,digits=3,format="g"),
                                       if(p_val<0.05)" *" else " ns") else "AnnotSV score por sexo"
    layout(fig, yaxis=list(title="AnnotSV score"), xaxis=list(title=""),
           showlegend=TRUE,
           title=list(text=titulo, font=list(size=12,
                                             color=if(!is.na(p_val)&&p_val<0.05)"#CC0000" else "#555")),
           margin=list(b=60,l=60,r=10,t=50), plot_bgcolor="#FAFAFA", paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_herencia_sexo <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0, "sexo" %in% names(df))
    df$her_limpia2 <- her_base(df$Tipo_Herencia)
    df$sexo_lbl <- ifelse(df$sexo=="M","Var\u00f3n \u2642",
                          ifelse(df$sexo=="F","Mujer \u2640","Desconocido"))
    conteo <- df |>
      filter(!is.na(her_limpia2)) |>
      group_by(sexo_lbl, her_limpia2) |>
      summarise(n=n(), .groups="drop")
    fig <- plot_ly()
    for (h in names(HER_COLORS)) {
      sub <- conteo[conteo$her_limpia2==h,]; if(nrow(sub)==0) next
      col <- HER_COLORS[h]
      fig <- add_trace(fig, type="bar", x=sub$sexo_lbl, y=sub$n, name=h,
                       marker=list(color=col),
                       hovertemplate=paste0("<b>",h,"</b><br>%{x}<br>%{y}<extra></extra>"))
    }
    layout(fig, barmode="stack",
           xaxis=list(title="", categoryorder="array",
                      categoryarray=c("Var\u00f3n \u2642","Mujer \u2640","Desconocido")),
           yaxis=list(title="N variantes"),
           legend=list(orientation="h",y=-0.3,x=0,font=list(size=10)),
           margin=list(b=100,l=50,r=10,t=10), plot_bgcolor="#FAFAFA", paper_bgcolor="#FFFFFF")
  })
  
  output$est_plot_tamano_sexo <- renderPlotly({
    df <- est_datos(); req(nrow(df)>0, "sexo" %in% names(df))
    df2 <- df[!is.na(df$tamano_kb) & df$tamano_kb>0 & df$sexo!="ND",]
    req(nrow(df2)>0)
    fig <- plot_ly()
    for (s in c("M","F")) {
      sub <- df2[df2$sexo==s,]; if(nrow(sub)==0) next
      col <- SEX_COLORS[s]; lbl <- SEX_LABELS[s]
      fig <- add_trace(fig, type="box", y=sub$tamano_kb, name=lbl,
                       marker=list(color=col,size=3), line=list(color=col),
                       fillcolor=paste0(col,"55"),
                       hovertemplate=paste0("<b>",lbl,"</b><br>Kb:%{y:.1f}<extra></extra>"))
    }
    p_val <- tryCatch({
      m <- df2$tamano_kb[df2$sexo=="M"]; f <- df2$tamano_kb[df2$sexo=="F"]
      if(length(m)>1 && length(f)>1) wilcox.test(m,f)$p.value else NA
    }, error=function(e) NA)
    titulo <- if(!is.na(p_val)) paste0("Wilcoxon p=",formatC(p_val,digits=3,format="g"),
                                       if(p_val<0.05)" *" else " ns") else "Tama\u00f1o por sexo"
    layout(fig, yaxis=list(title="Tama\u00f1o (Kb)",type="log"), xaxis=list(title=""),
           showlegend=TRUE,
           title=list(text=titulo, font=list(size=12,
                                             color=if(!is.na(p_val)&&p_val<0.05)"#CC0000" else "#555")),
           margin=list(b=60,l=60,r=10,t=50), plot_bgcolor="#FAFAFA", paper_bgcolor="#FFFFFF")
  })
  
  output$ui_badge_global <- renderUI({
    n_cnv <- nrow(rv$resultados$CNVs %||% data.frame())
    n_sv  <- nrow(rv$resultados$SVs  %||% data.frame())
    if (n_cnv + n_sv == 0) return(tags$span(class="badge bg-secondary ms-2", "Sin datos cargados"))
    tags$span(class = "badge bg-success ms-2", paste0(n_cnv, " CNVs · ", n_sv, " SVs cargados"))
  })
  
  # ── Chromosomal ideogram ──────────────────────────────────────────────────
  # ── Comparador de familias ────────────────────────────────────────────────
  comp_data <- eventReactive(input$btn_comparar, {
    req(nzchar(input$comp_fam_a), nzchar(input$comp_fam_b))
    fam_a  <- input$comp_fam_a
    fam_b  <- input$comp_fam_b
    mod    <- input$comp_modalidad
    umbral <- (input$comp_umbral %||% 50) / 100
    datos  <- rv$resultados
    
    todas <- switch(mod,
                    "CNVs"  = datos$CNVs %||% data.frame(),
                    "SVs"   = datos$SVs  %||% data.frame(),
                    "Ambas" = bind_rows(
                      if (!is.null(datos$CNVs)) mutate(datos$CNVs, Modalidad="CNV") else data.frame(),
                      if (!is.null(datos$SVs))  mutate(datos$SVs,  Modalidad="SV")  else data.frame()
                    )
    )
    if (nrow(todas) == 0) return(NULL)
    
    df_ab <- todas[todas$ID_Familia %in% c(fam_a, fam_b), , drop=FALSE]
    if ("Annotation_mode" %in% names(df_ab))
      df_ab <- df_ab[!is.na(df_ab$Annotation_mode) & df_ab$Annotation_mode == "full", , drop=FALSE]
    if (nrow(df_ab) == 0) return(NULL)
    
    chr_v <- toupper(trimws(gsub("(?i)^chr","",as.character(df_ab$SV_chrom),perl=TRUE)))
    sta_v <- suppressWarnings(as.numeric(df_ab$SV_start))
    end_v <- suppressWarnings(as.numeric(df_ab$SV_end))
    fam_v <- as.character(df_ab$ID_Familia)
    len_v <- pmax(1L, end_v - sta_v + 1L, na.rm=TRUE)
    
    compartida <- logical(nrow(df_ab))
    for (i in seq_len(nrow(df_ab))) {
      if (is.na(sta_v[i]) || is.na(end_v[i])) next
      otra <- fam_v != fam_v[i] & chr_v == chr_v[i] & !is.na(sta_v) & !is.na(end_v)
      if (!any(otra)) next
      ovl <- pmax(0, pmin(end_v[i], end_v[otra]) - pmax(sta_v[i], sta_v[otra]) + 1L)
      sim <- ovl / pmax(len_v[i], len_v[otra])
      compartida[i] <- any(sim >= umbral, na.rm=TRUE)
    }
    
    df_ab$chr_clean   <- chr_v
    df_ab$sta_num     <- sta_v
    df_ab$end_num     <- end_v
    df_ab$Compartida  <- compartida
    df_ab$Estado      <- ifelse(compartida, "Compartida", paste0("Solo ", fam_v))
    df_ab$fam_v       <- fam_v
    df_ab$her_b       <- her_base(df_ab$Tipo_Herencia)
    
    list(df=df_ab, fam_a=fam_a, fam_b=fam_b)
  })
  
  # Tarjetas resumen
  output$ui_comp_resumen <- renderUI({
    cd <- comp_data(); req(!is.null(cd))
    df <- cd$df; fa <- cd$fam_a; fb <- cd$fam_b
    n_a    <- sum(df$fam_v == fa)
    n_b    <- sum(df$fam_v == fb)
    n_comp <- sum(df$Compartida)
    n_solo_a <- sum(df$fam_v == fa & !df$Compartida)
    n_solo_b <- sum(df$fam_v == fb & !df$Compartida)
    mk <- function(val, label, bg, col) {
      div(class="text-center rounded p-2 mb-1",
          style=paste0("background:", bg, "; color:", col, ";"),
          tags$div(style="font-size:1.5em; font-weight:700;", val),
          tags$div(style="font-size:0.75em;", label))
    }
    tagList(
      mk(n_a,      paste0("Variantes ", fa), "#EBF5FB", "#2C6FAC"),
      mk(n_b,      paste0("Variantes ", fb), "#FEF5E7", "#E67E22"),
      mk(n_comp,   "Shared (≥threshold)", "#D4EDDA", "#155724"),
      mk(n_solo_a, paste0("Solo en ", fa),   "#EBF5FB", "#2C6FAC"),
      mk(n_solo_b, paste0("Solo en ", fb),   "#FEF5E7", "#E67E22")
    )
  })
  
  output$ui_comp_placeholder <- renderUI({
    if (is.null(comp_data()))
      div(class="text-center text-muted py-5",
          bsicons::bs_icon("bar-chart-steps", size="2em"), br(),
          tags$span("Selecciona dos familias y pulsa ", tags$b("Comparar familias")))
    else NULL
  })
  
  output$ui_comp_leyenda <- renderUI({
    cd <- comp_data(); if (is.null(cd)) return(NULL)
    div(class="d-flex gap-3 align-items-center flex-wrap",
        div(class="d-flex align-items-center gap-1",
            div(style="width:12px;height:12px;background:#2C6FAC;border-radius:2px;"),
            tags$small(cd$fam_a)),
        div(class="d-flex align-items-center gap-1",
            div(style="width:12px;height:12px;background:#E67E22;border-radius:2px;"),
            tags$small(cd$fam_b)),
        div(class="d-flex align-items-center gap-1",
            div(style="width:12px;height:12px;background:#27AE60;border-radius:2px;"),
            tags$small("Compartida"))
    )
  })
  
  # Comparative genomic plot: A/B tracks per chromosome
  output$plot_comparador_genomico <- renderPlotly({
    cd <- comp_data(); req(!is.null(cd))
    df <- cd$df; fa <- cd$fam_a; fb <- cd$fam_b
    
    chrs_presentes <- intersect(CHR_ORDER, unique(df$chr_clean))
    req(length(chrs_presentes) > 0)
    
    chr_idx  <- setNames(seq_along(chrs_presentes), chrs_presentes)
    n_chrs   <- length(chrs_presentes)
    TRACK_H  <- 0.35   # mitad de la altura de cada pista (A arriba, B abajo)
    COLOR_A  <- "#2C6FAC"
    COLOR_B  <- "#E67E22"
    COLOR_AB <- "#27AE60"
    
    fig <- plot_ly()
    
    # Fondo de cromosomas
    for (chr in chrs_presentes) {
      len <- CHR_LENGTHS[chr] %||% 1e8
      y_c <- chr_idx[chr]
      fig <- add_trace(fig, type="bar", orientation="h",
                       x=len, y=y_c, base=0, width=TRACK_H*2+0.05,
                       marker=list(color="#EEEEEE", line=list(color="#CCCCCC", width=0.5)),
                       hoverinfo="none", showlegend=FALSE)
    }
    
    # Function to add a track for a family
    add_family_track <- function(fig, df_fam, fam_id, y_offset, color, name, show_leg) {
      if (nrow(df_fam) == 0) return(fig)
      for (chr in chrs_presentes) {
        dfc <- df_fam[df_fam$chr_clean == chr & !is.na(df_fam$sta_num) & !is.na(df_fam$end_num), ]
        if (nrow(dfc) == 0) next
        y_c    <- chr_idx[chr] + y_offset
        widths <- pmax(dfc$end_num - dfc$sta_num, CHR_LENGTHS[chr] * 0.012, na.rm=TRUE)
        cols   <- ifelse(dfc$Compartida, COLOR_AB, color)
        for (j in seq_len(nrow(dfc))) {
          fig <- add_trace(fig, type="bar", orientation="h",
                           x=widths[j], y=y_c, base=dfc$sta_num[j], width=TRACK_H,
                           marker=list(color=paste0(cols[j],"CC"), line=list(color=cols[j], width=1.2)),
                           name=if (j==1) name else name, legendgroup=name,
                           showlegend=(j==1 && show_leg),
                           text=paste0(
                             "<b>", fam_id, "</b>  ", dfc$SV_type[j], "<br>",
                             "chr", chr, ":", format(round(dfc$sta_num[j]), big.mark=".", scientific=FALSE),
                             "–", format(round(dfc$end_num[j]), big.mark=".", scientific=FALSE), "<br>",
                             "Herencia: ", dfc$her_b[j], "<br>",
                             if (dfc$Compartida[j]) "<b style='color:#27AE60'>✅ Compartida con la otra familia</b>"
                             else "Solo en esta familia"
                           ),
                           hovertemplate="%{text}<extra></extra>"
          )
        }
      }
      fig
    }
    
    df_a <- df[df$fam_v == fa, ]
    df_b <- df[df$fam_v == fb, ]
    fig <- add_family_track(fig, df_a, fa, +TRACK_H/2, COLOR_A, fa, TRUE)
    fig <- add_family_track(fig, df_b, fb, -TRACK_H/2, COLOR_B, fb, TRUE)
    
    # Central separator line per chromosome
    for (chr in chrs_presentes) {
      len <- CHR_LENGTHS[chr] %||% 1e8
      fig <- add_trace(fig, type="scatter", mode="lines",
                       x=c(0, len), y=c(chr_idx[chr], chr_idx[chr]),
                       line=list(color="#CCCCCC", width=0.5, dash="dot"),
                       hoverinfo="none", showlegend=FALSE)
    }
    
    layout(fig,
           barmode  = "overlay",
           xaxis    = list(title="Genomic position (bp)", tickformat=",.0f",
                           showgrid=TRUE, gridcolor="#F5F5F5"),
           yaxis    = list(title="", tickvals=seq_along(chrs_presentes),
                           ticktext=paste0("chr", chrs_presentes),
                           autorange="reversed", tickfont=list(size=11)),
           legend   = list(orientation="h", y=-0.1, x=0, font=list(size=11),
                           title=list(text="<b>Familia</b>")),
           annotations = list(
             list(x=0.02, y=1.01, xref="paper", yref="paper", showarrow=FALSE,
                  text=paste0("<b style='color:", COLOR_A, "'>▲ ", fa, "</b>  |  ",
                              "<b style='color:", COLOR_B, "'>▼ ", fb, "</b>  |  ",
                              "<b style='color:", COLOR_AB, "'>■ Compartida</b>"),
                  font=list(size=11), align="left")
           ),
           margin = list(l=70, r=20, t=35, b=60),
           plot_bgcolor="#FAFAFA", paper_bgcolor="#FFFFFF", height=520
    )
  })
  
  output$tabla_comparador <- renderDT({
    cd <- comp_data()
    if (is.null(cd))
      return(datatable(data.frame(Mensaje="Selecciona dos familias y pulsa Comparar."),
                       options=list(dom="t"), rownames=FALSE))
    df  <- cd$df; fa <- cd$fam_a; fb <- cd$fam_b
    
    df$Estado <- ifelse(df$Compartida, "✅ Compartida",
                        ifelse(df$fam_v == fa,
                               paste0("🔵 Solo ", fa),
                               paste0("🟠 Solo ", fb)))
    
    cols_vis <- intersect(c("ID_Familia","Estado","SV_chrom","SV_start","SV_end",
                            "SV_length","SV_type","Tipo_Herencia","Tipo_Rango",
                            "AnnotSV_ranking_score","Gene_name"), names(df))
    df_show <- df[order(df$SV_chrom, df$sta_num), cols_vis, drop=FALSE]
    
    # ── LECTURA DEL MODO OSCURO ──
    modo_oscuro <- isTRUE(input$modo_noche == "dark")
    txt_color   <- if (modo_oscuro) "white" else NULL
    
    bg_est <- if (modo_oscuro) c("#27AE60", "#2C6FAC", "#E67E22") else c("#D4EDDA", "#EBF5FB", "#FEF5E7")
    bg_her <- if (modo_oscuro) c("#CC0000","#0066CC","#CC6600","#6600CC","#4499DD","#DD9944","#9944DD","#6c757d","#6c757d")
    else c("#FFE6E6","#E6F3FF","#FFF0E6","#F0E6FF","#CCE5FF","#FFE5CC","#EBD7FF","#E8E8E8","#E8E8E8")
    
    datatable(df_show, rownames=FALSE, filter="top",
              options=list(dom="frtip", pageLength=20, scrollX=TRUE,
                           scrollY="280px", scroller=TRUE),
              extensions="Scroller",
              class="table table-hover table-sm") |>
      formatStyle("Estado",
                  backgroundColor = styleEqual(c("✅ Compartida", paste0("🔵 Solo ", fa), paste0("🟠 Solo ", fb)), bg_est),
                  color = txt_color, fontWeight = styleEqual("✅ Compartida", "bold")) |>
      formatStyle("ID_Familia",
                  color = styleEqual(c(fa, fb), c(COLOR_A, COLOR_B)),
                  fontWeight = "bold") |>
      formatStyle("Tipo_Herencia",
                  backgroundColor = styleEqual(
                    c("De novo","Paternal","Maternal","Combined",
                      "Paternal (Probable)","Maternal (Probable)","Combined (Probable)",
                      "Desconocido","Unknown"), bg_her),
                  color = txt_color, fontWeight = styleEqual("De novo", "bold"))
  }, server=TRUE)
  
  # ── Vista de locus compartido ──────────────────────────────────────────────
  output$ui_locus_compartido <- renderUI({
    fd <- rv$fila_data
    if (is.null(fd) || nrow(fd) == 0) return(NULL)
    
    chr_ref <- toupper(trimws(gsub("(?i)^chr", "", as.character(fd$SV_chrom[1]), perl=TRUE)))
    sta_ref <- suppressWarnings(as.numeric(fd$SV_start[1]))
    end_ref <- suppressWarnings(as.numeric(fd$SV_end[1]))
    if (is.na(sta_ref) || is.na(end_ref)) return(NULL)
    
    todas_df <- tryCatch({
      bind_rows(
        if (!is.null(rv$resultados$CNVs)) mutate(rv$resultados$CNVs, Modalidad="CNV") else data.frame(),
        if (!is.null(rv$resultados$SVs))  mutate(rv$resultados$SVs,  Modalidad="SV")  else data.frame()
      )
    }, error = function(e) data.frame())
    if (nrow(todas_df) == 0) return(NULL)
    if ("Annotation_mode" %in% names(todas_df))
      todas_df <- todas_df[!is.na(todas_df$Annotation_mode) & todas_df$Annotation_mode == "full", , drop=FALSE]
    
    chr_vec <- toupper(trimws(gsub("(?i)^chr", "", as.character(todas_df$SV_chrom), perl=TRUE)))
    sta_vec <- suppressWarnings(as.numeric(todas_df$SV_start))
    end_vec <- suppressWarnings(as.numeric(todas_df$SV_end))
    margen  <- max((end_ref - sta_ref) * 0.25, 75000)
    mask    <- chr_vec == chr_ref & !is.na(sta_vec) & !is.na(end_vec) &
      sta_vec < (end_ref + margen) & end_vec > (sta_ref - margen)
    df_loc  <- todas_df[mask, , drop=FALSE]
    if (nrow(df_loc) == 0) return(NULL)
    n_fams  <- length(unique(df_loc$ID_Familia))
    if (n_fams < 2) return(NULL)
    
    # dynamic height based on number of families
    plot_h <- max(180, n_fams * 45 + 80)
    
    # Genes affected in the region (if available)
    genes_reg <- if ("Gene_name" %in% names(fd)) {
      gn <- trimws(as.character(fd$Gene_name[1]))
      gn <- gsub("[./]", ", ", gn)
      if (nzchar(gn) && gn != "NA") gn else NULL
    } else NULL
    
    len_ref <- end_ref - sta_ref
    len_fmt <- if (len_ref >= 1e6) paste0(round(len_ref/1e6, 2), " Mb")
    else paste0(round(len_ref/1e3, 1), " kb")
    
    tagList(
      hr(),
      div(class="d-flex align-items-start justify-content-between flex-wrap gap-1 mb-1",
          div(
            tags$small(class="fw-bold text-primary",
                       paste0("🔍 Shared locus — ", n_fams, " family/families")),
            tags$small(class="text-muted d-block",
                       paste0("chr", chr_ref, ":",
                              format(round(sta_ref - margen), big.mark=".", scientific=FALSE), "–",
                              format(round(end_ref + margen), big.mark=".", scientific=FALSE),
                              " · variante: ", len_fmt))
          ),
          if (!is.null(genes_reg))
            div(tags$small(
              class = "dm-loc-badge rounded", style="border-radius:4px; padding:2px 7px; font-size:0.82em; font-weight:600; border:1px solid;",
              paste0("Genes: ", genes_reg)
            ))
      ),
      div(class="text-muted", style="font-size:0.78em; margin-top:2px;",
          "🖱 Click a bar to load that variant in the detail panel"),
      plotlyOutput("plot_locus", height = paste0(plot_h, "px"))
    )
  })
  output$plot_locus <- renderPlotly({
    fd <- rv$fila_data
    req(!is.null(fd), nrow(fd) > 0)
    
    chr_ref <- toupper(trimws(gsub("(?i)^chr", "", as.character(fd$SV_chrom[1]), perl=TRUE)))
    sta_ref <- suppressWarnings(as.numeric(fd$SV_start[1]))
    end_ref <- suppressWarnings(as.numeric(fd$SV_end[1]))
    req(!is.na(sta_ref), !is.na(end_ref))
    fam_ref <- as.character(fd$ID_Familia[1])
    
    todas_df <- tryCatch({
      bind_rows(
        if (!is.null(rv$resultados$CNVs)) mutate(rv$resultados$CNVs, Modalidad="CNV") else data.frame(),
        if (!is.null(rv$resultados$SVs))  mutate(rv$resultados$SVs,  Modalidad="SV")  else data.frame()
      )
    }, error = function(e) data.frame())
    req(nrow(todas_df) > 0)
    if ("Annotation_mode" %in% names(todas_df))
      todas_df <- todas_df[!is.na(todas_df$Annotation_mode) & todas_df$Annotation_mode == "full", , drop=FALSE]
    
    chr_vec <- toupper(trimws(gsub("(?i)^chr", "", as.character(todas_df$SV_chrom), perl=TRUE)))
    sta_vec <- suppressWarnings(as.numeric(todas_df$SV_start))
    end_vec <- suppressWarnings(as.numeric(todas_df$SV_end))
    margen  <- max((end_ref - sta_ref) * 0.25, 75000)
    
    mask    <- chr_vec == chr_ref & !is.na(sta_vec) & !is.na(end_vec) &
      sta_vec < (end_ref + margen) & end_vec > (sta_ref - margen)
    df_loc  <- todas_df[mask, , drop=FALSE]
    req(nrow(df_loc) > 0)
    
    df_loc$her_b    <- her_base(df_loc$Tipo_Herencia)
    df_loc$sta_num  <- suppressWarnings(as.numeric(df_loc$SV_start))
    df_loc$end_num  <- suppressWarnings(as.numeric(df_loc$SV_end))
    df_loc$score_n  <- suppressWarnings(as.numeric(df_loc$AnnotSV_ranking_score))
    df_loc$fam_id   <- as.character(df_loc$ID_Familia)
    df_loc$es_ref   <- df_loc$fam_id == fam_ref & df_loc$sta_num == sta_ref & df_loc$end_num == end_ref
    df_loc$sv_type  <- if ("SV_type" %in% names(df_loc)) as.character(df_loc$SV_type) else rep("", nrow(df_loc))
    df_loc$gen_name <- if ("Gene_name" %in% names(df_loc)) as.character(df_loc$Gene_name) else rep("", nrow(df_loc))
    
    fams_ord <- c(fam_ref, sort(setdiff(unique(df_loc$fam_id), fam_ref)))
    df_loc$fam_y <- match(df_loc$fam_id, fams_ord)
    n_fams <- length(fams_ord)
    
    # Build shapes: alternating rows + reference region shadow
    shapes_list <- c(
      # Filas alternas de fondo
      lapply(seq_along(fams_ord), function(i) {
        if (i %% 2 == 0)
          list(type="rect", xref="x", yref="y",
               x0=sta_ref - margen, x1=end_ref + margen,
               y0=i - 0.45, y1=i + 0.45,
               fillcolor="rgba(245,245,245,0.6)",
               line=list(width=0), layer="below")
        else NULL
      }),
      # Reference region shadow
      list(list(type="rect", xref="x", yref="y",
                x0=sta_ref, x1=end_ref,
                y0=0.4, y1=n_fams + 0.6,
                fillcolor="rgba(44,111,172,0.07)",
                line=list(color="rgba(44,111,172,0.4)", width=1, dash="dot"),
                layer="below"))
    )
    shapes_list <- Filter(Negate(is.null), shapes_list)
    
    fig <- plot_ly(source = "plot_locus")
    
    
    for (her in unique(df_loc$her_b)) {
      dh  <- df_loc[!is.na(df_loc$her_b) & df_loc$her_b == her, , drop=FALSE]
      col <- HER_COLORS[her]; if (is.na(col)) col <- "#888888"
      # Reference variant thicker and more opaque
      widths <- pmax(dh$end_num - dh$sta_num, margen * 0.03)
      bar_w  <- ifelse(dh$es_ref, 0.65, 0.45)
      opac   <- ifelse(dh$es_ref, "EE", "99")
      score_lbl <- ifelse(is.na(dh$score_n), "—", as.character(round(dh$score_n, 2)))
      len_lbl <- sapply(dh$end_num - dh$sta_num, function(l) {
        if (is.na(l)) "—"
        else if (l >= 1e6) paste0(round(l/1e6,2)," Mb")
        else paste0(round(l/1e3,1)," kb")
      })
      # customdata: identificador de variante para navegacion al hacer click
      cdata <- paste0(dh$fam_id, "|||", dh$sta_num, "|||", dh$end_num)
      fig <- add_trace(fig,
                       type="bar", orientation="h",
                       x     = widths,
                       y     = dh$fam_y,
                       base  = dh$sta_num,
                       width = bar_w,
                       customdata = cdata,
                       marker = list(color=paste0(col, opac), line=list(color=col, width=ifelse(dh$es_ref,2,1))),
                       name=her, legendgroup=her, showlegend=TRUE,
                       text=paste0(
                         "<b>", dh$fam_id, "</b>", ifelse(dh$es_ref, " ⬅ variante seleccionada", ""), "<br>",
                         "Tipo: ", dh$sv_type, " · Herencia: ", her, "<br>",
                         "chr", chr_ref, ":", format(round(dh$sta_num), big.mark=".", scientific=FALSE),
                         "–", format(round(dh$end_num), big.mark=".", scientific=FALSE), " (", len_lbl, ")<br>",
                         "Score AnnotSV: <b>", score_lbl, "</b><br>",
                         ifelse(nzchar(trimws(dh$gen_name)), paste0("Genes: ", dh$gen_name), "")
                       ),
                       hovertemplate="%{text}<extra></extra>"
      )
    }
    
    # Etiquetas de tipo SV centradas en cada barra
    annots <- lapply(seq_len(nrow(df_loc)), function(i) {
      row <- df_loc[i, ]
      mid <- (row$sta_num + row$end_num) / 2
      lbl <- if (nzchar(row$sv_type)) row$sv_type else ""
      if (!nzchar(lbl)) return(NULL)
      list(x=mid, y=row$fam_y, text=lbl, showarrow=FALSE,
           font=list(size=9, color="#333"), xanchor="center", yanchor="middle")
    })
    annots <- Filter(Negate(is.null), annots)
    
    fig <- event_register(fig, "plotly_click")
    layout(fig,
           barmode="overlay",
           xaxis=list(
             title="Genomic position (bp)", tickformat=",.0f",
             showgrid=TRUE, gridcolor="#EEEEEE",
             range=list(sta_ref - margen, end_ref + margen)
           ),
           yaxis=list(
             title="", tickvals=seq_along(fams_ord), ticktext=fams_ord,
             autorange="reversed", tickfont=list(size=11),
             range=list(0.4, n_fams + 0.6)
           ),
           legend=list(orientation="h", y=-0.22, x=0, font=list(size=10),
                       title=list(text="<b>Herencia</b>")),
           annotations=annots,
           shapes=shapes_list,
           margin=list(l=90, r=15, t=10, b=65),
           plot_bgcolor="#FAFAFA", paper_bgcolor="#FFFFFF"
    )
  })
  
  
  
  # ---- Click en locus compartido: navegar a la variante clicada --------
  observeEvent(event_data("plotly_click", source = "plot_locus"), {
    ed <- event_data("plotly_click", source = "plot_locus")
    req(!is.null(ed), !is.null(ed$customdata))
    partes <- strsplit(ed$customdata, "\\|\\|\\|")[[1]]
    req(length(partes) == 3)
    click_fam <- partes[1]
    click_sta <- suppressWarnings(as.numeric(partes[2]))
    click_end <- suppressWarnings(as.numeric(partes[3]))
    req(!is.na(click_sta), !is.na(click_end))
    
    # Buscar la variante en todos los datos (sin filtros)
    todas_df <- tryCatch({
      bind_rows(
        if (!is.null(rv$resultados$CNVs)) mutate(rv$resultados$CNVs, Modalidad="CNV") else data.frame(),
        if (!is.null(rv$resultados$SVs))  mutate(rv$resultados$SVs,  Modalidad="SV")  else data.frame()
      )
    }, error = function(e) data.frame())
    req(nrow(todas_df) > 0)
    
    sta_v  <- suppressWarnings(as.numeric(todas_df$SV_start))
    end_v  <- suppressWarnings(as.numeric(todas_df$SV_end))
    fam_v  <- as.character(todas_df$ID_Familia)
    mode_v <- if ("Annotation_mode" %in% names(todas_df)) todas_df$Annotation_mode else rep("full", nrow(todas_df))
    
    idx <- which(fam_v == click_fam &
                   !is.na(sta_v) & abs(sta_v - click_sta) < 2 &
                   !is.na(end_v) & abs(end_v - click_end) < 2 &
                   (is.na(mode_v) | mode_v == "full"))
    
    if (length(idx) == 0) {
      showNotification(paste0("Variante de ", click_fam, " no encontrada (puede estar oculta por filtros activos)."),
                       type = "warning", duration = 4)
      return()
    }
    
    # Cargar datos en el panel de detalle
    fila <- todas_df[idx[1], , drop = FALSE]
    if ("._archivo_" %in% names(fila)) rv$fila_seleccionada <- fila$`._archivo_`[1]
    rv$fila_data <- fila
    
    # Guardar la seleccion pendiente: el observer reactivo de abajo la resolvera
    # una vez que datos_filtrados() haya recomputado con los filtros correctos
    rv$pendiente_seleccion <- list(fam = click_fam, sta = click_sta, end = click_end)
    
    # Limpiar filtro de familia si la variante quedaria excluida
    fams_actuales <- isolate(input$filtro_familia)
    if (!is.null(fams_actuales) && length(fams_actuales) > 0 && !(click_fam %in% fams_actuales)) {
      updateSelectizeInput(session, "filtro_familia", selected = character(0))
    }
    
    # Navegar al tab Resultados — la tabla se re-renderiza y el observer de abajo
    # seleccionara la fila cuando datos_filtrados() haya recomputado
    updateNavbarPage(session, "main_navbar", selected = "tab_results")
    
    showNotification(
      paste0("🔍 Variant from ", click_fam, " loaded — navigating to table..."),
      type = "message", duration = 3
    )
  })
  
  # ---- Observer reactivo: selecciona la fila en la tabla cuando datos_filtrados()
  #      ya incluye la variante pendiente (se dispara tras cualquier cambio de filtros) ----
  observe({
    ps <- rv$pendiente_seleccion
    req(!is.null(ps))
    
    df_filtrado <- datos_filtrados()   # reactivo: se reevalua tras los cambios de filtro
    req(nrow(df_filtrado) > 0)
    
    sta_f  <- suppressWarnings(as.numeric(df_filtrado$SV_start))
    end_f  <- suppressWarnings(as.numeric(df_filtrado$SV_end))
    fam_f  <- as.character(df_filtrado$ID_Familia)
    mode_f <- if ("Annotation_mode" %in% names(df_filtrado)) df_filtrado$Annotation_mode else rep("full", nrow(df_filtrado))
    
    fila_idx <- which(fam_f == ps$fam &
                        !is.na(sta_f) & abs(sta_f - ps$sta) < 2 &
                        !is.na(end_f) & abs(end_f - ps$end) < 2 &
                        (is.na(mode_f) | mode_f == "full"))
    
    req(length(fila_idx) > 0)   # si la variante aun no esta visible, esperar
    
    # Limpiar la seleccion pendiente antes de actuar para no repetir
    rv$pendiente_seleccion <- NULL
    
    dt_proxy <- dataTableProxy("tabla_variantes")
    shinyjs::delay(150, {
      selectRows(dt_proxy, fila_idx[1])
    })
  })
  
  # Reactive: tabla de genes agregada (usada en ideograma y detalle)
  genes_agregados <- reactive({
    df <- datos_filtrados()
    if (nrow(df) == 0 || !"Gene_name" %in% names(df)) return(NULL)
    if (!all(c("SV_chrom","SV_start","SV_end") %in% names(df))) return(NULL)
    
    gen_col <- as.character(df$Gene_name)
    n_genes <- lengths(strsplit(gen_col, "[/|,]"))
    idx_rep <- rep(seq_len(nrow(df)), n_genes)
    
    # Construir data.frame limpio desde cero para evitar problemas de dplyr con rownames
    df_exp <- data.frame(
      gene_solo  = trimws(unlist(strsplit(gen_col, "[/|,]"))),
      chr_clean  = toupper(trimws(gsub("(?i)^chr", "", as.character(df$SV_chrom[idx_rep]), perl=TRUE))),
      sta_num    = suppressWarnings(as.numeric(df$SV_start[idx_rep])),
      end_num    = suppressWarnings(as.numeric(df$SV_end[idx_rep])),
      ID_Familia = as.character(df$ID_Familia[idx_rep]),
      her_b      = her_base(df$Tipo_Herencia[idx_rep]),
      es_sfari   = if ("En_Panel_Genes" %in% names(df))
        es_sfari_col(df$En_Panel_Genes[idx_rep])
      else rep(FALSE, length(idx_rep)),
      score_num  = suppressWarnings(as.numeric(df$AnnotSV_ranking_score[idx_rep])),
      stringsAsFactors = FALSE
    )
    df_exp$pos_mid <- (df_exp$sta_num + df_exp$end_num) / 2
    df_exp <- df_exp[nzchar(df_exp$gene_solo) & df_exp$gene_solo != "." &
                       df_exp$chr_clean %in% CHR_ORDER &
                       !is.na(df_exp$pos_mid), , drop=FALSE]
    if (nrow(df_exp) == 0) return(NULL)
    
    # Agregar por gen+cromosoma con base R
    claves <- paste(df_exp$gene_solo, df_exp$chr_clean, sep="__")
    claves_unicas <- unique(claves)
    
    result <- do.call(rbind, lapply(claves_unicas, function(k) {
      rows <- df_exp[claves == k, , drop=FALSE]
      data.frame(
        gen_solo    = rows$gene_solo[1],
        chr_clean   = rows$chr_clean[1],
        pos_media   = mean(rows$pos_mid, na.rm=TRUE),
        N_familias  = length(unique(rows$ID_Familia)),
        N_variantes = nrow(rows),
        Familias    = paste(sort(unique(rows$ID_Familia)), collapse=", "),
        Herencias   = paste(sort(unique(na.omit(rows$her_b))), collapse=" / "),
        En_SFARI    = any(rows$es_sfari, na.rm=TRUE),
        Score_medio = round(mean(rows$score_num, na.rm=TRUE), 3),
        stringsAsFactors = FALSE
      )
    }))
    
    result[order(-result$N_familias, -result$N_variantes), , drop=FALSE]
  })
  
  output$ui_genes_badge <- renderUI({
    dg <- genes_agregados()
    n  <- if (is.null(dg)) 0 else n_distinct(dg$gen_solo)
    tags$span(paste0(n, " genes"))
  })
  
  output$ui_genes_placeholder <- renderUI({
    if (is.null(genes_agregados()))
      div(class="text-center text-muted py-5",
          bsicons::bs_icon("map", size="2em"), br(),
          tags$span("Carga resultados para ver el ideograma de genes"))
    else NULL
  })
  
  output$plot_genes_ideograma <- renderPlotly({
    dg <- genes_agregados(); req(!is.null(dg))
    if (isTRUE(input$genes_solo_sfari)) dg <- dg[dg$En_SFARI, , drop=FALSE]
    req(nrow(dg) > 0)
    
    # Size and colour scale by number of families
    max_fam  <- max(dg$N_familias, 1)
    dg$sz    <- scales::rescale(dg$N_familias, to=c(6, 18), from=c(1, max_fam))
    dg$color <- ifelse(dg$En_SFARI, "#E74C3C", "#2C6FAC")   # rojo=SFARI, azul=no SFARI
    # More families → more opaque
    dg$alpha <- scales::rescale(dg$N_familias, to=c(0.55, 1), from=c(1, max_fam))
    dg$color_hex <- mapply(function(c, a) {
      rgb_v <- col2rgb(c)/255
      sprintf("rgba(%d,%d,%d,%.2f)", round(rgb_v[1]*255), round(rgb_v[2]*255), round(rgb_v[3]*255), a)
    }, dg$color, dg$alpha)
    
    chr_idx <- setNames(seq_along(CHR_ORDER), CHR_ORDER)
    
    fig <- plot_ly(source = "plot_genes_ideograma")
    
    # Barras de cromosomas
    for (chr in CHR_ORDER) {
      if (!chr %in% names(CHR_LENGTHS)) next
      fig <- add_trace(fig, type="bar", orientation="h",
                       x=CHR_LENGTHS[chr], y=chr_idx[chr], base=0, width=0.6,
                       marker=list(color="#E8E8E8", line=list(color="#CCCCCC", width=0.7)),
                       hoverinfo="none", showlegend=FALSE)
    }
    
    # Genes no SFARI
    dg_n <- dg[!dg$En_SFARI, ]
    if (nrow(dg_n) > 0) {
      fig <- add_trace(fig,
                       type="scatter", mode="markers",
                       x    = dg_n$pos_media,
                       y    = chr_idx[dg_n$chr_clean],
                       customdata = dg_n$gen_solo,
                       marker = list(
                         size   = dg_n$sz,
                         color  = dg_n$color_hex,
                         symbol = "line-ns",
                         line   = list(color=dg_n$color_hex, width=dg_n$sz)
                       ),
                       name          = "Gen afectado",
                       legendgroup   = "normal",
                       showlegend    = TRUE,
                       text = paste0(
                         "<b>", dg_n$gen_solo, "</b><br>",
                         "Familias afectadas: <b>", dg_n$N_familias, "</b><br>",
                         dg_n$Familias, "<br>",
                         "Herencias: ", dg_n$Herencias, "<br>",
                         "Score medio: ", dg_n$Score_medio
                       ),
                       hovertemplate = "%{text}<extra></extra>"
      )
    }
    
    # Panel genes (above, different symbol)
    dg_s <- dg[dg$En_SFARI, ]
    if (nrow(dg_s) > 0) {
      fig <- add_trace(fig,
                       type="scatter", mode="markers",
                       x    = dg_s$pos_media,
                       y    = chr_idx[dg_s$chr_clean],
                       customdata = dg_s$gen_solo,
                       marker = list(
                         size   = dg_s$sz + 2,
                         color  = dg_s$color_hex,
                         symbol = "line-ns",
                         line   = list(color=dg_s$color_hex, width=dg_s$sz + 2)
                       ),
                       name          = "Gen SFARI",
                       legendgroup   = "sfari",
                       showlegend    = TRUE,
                       text = paste0(
                         "<b>⭐ ", dg_s$gen_solo, " (panel)</b><br>",
                         "Familias afectadas: <b>", dg_s$N_familias, "</b><br>",
                         dg_s$Familias, "<br>",
                         "Herencias: ", dg_s$Herencias, "<br>",
                         "Score medio: ", dg_s$Score_medio
                       ),
                       hovertemplate = "%{text}<extra></extra>"
      )
    }
    
    layout(fig,
           barmode = "overlay",
           xaxis   = list(title="Genomic position (bp)", tickformat=",.0f",
                          showgrid=TRUE, gridcolor="#F0F0F0"),
           yaxis   = list(title="", tickvals=seq_along(CHR_ORDER),
                          ticktext=paste0("chr", CHR_ORDER),
                          autorange="reversed", tickfont=list(size=11)),
           legend  = list(orientation="h", y=-0.08, x=0, font=list(size=11)),
           annotations = list(list(
             x=0, y=1.02, xref="paper", yref="paper", showarrow=FALSE, align="left",
             text="<span style='color:#2C6FAC'>● Affected gene</span>   <span style='color:#E74C3C'>● Panel gene</span>   Size = no. affected families · Click for detail",
             font=list(size=11)
           )),
           margin = list(l=65, r=15, t=35, b=60),
           plot_bgcolor="#FAFAFA", paper_bgcolor="#FFFFFF", height=700
    ) |> event_register("plotly_click")
  })
  
  # Clic en gen → panel de detalle (via event_data con customdata = nombre del gen)
  # Usamos un reactiveVal para compartir el gen seleccionado entre ui y tabla
  gen_seleccionado <- reactiveVal(NULL)
  
  observeEvent(event_data("plotly_click", source = "plot_genes_ideograma"), {
    ed <- event_data("plotly_click", source = "plot_genes_ideograma")
    req(!is.null(ed))
    # customdata contains the gene name directly (assigned in add_trace)
    gen_name <- ed$customdata
    if (!is.null(gen_name) && length(gen_name) > 0 && nzchar(gen_name[1])) {
      gen_seleccionado(gen_name[1])
    }
  })
  
  output$ui_gen_detalle <- renderUI({
    gen_name <- gen_seleccionado()
    dg       <- genes_agregados()
    
    if (is.null(gen_name) || is.null(dg))
      return(div(class="text-center text-muted py-4 small",
                 bsicons::bs_icon("cursor"), br(), "Haz clic en un gen para ver el detalle"))
    
    gen_info <- dg[dg$gen_solo == gen_name, , drop=FALSE]
    if (nrow(gen_info) == 0)
      return(div(class="text-center text-muted py-4 small", paste0("Gen no encontrado: ", gen_name)))
    gen_info <- gen_info[1, ]
    
    sfari_tag <- if (isTRUE(gen_info$En_SFARI))
      tags$span(class="badge ms-1", style="background:#E74C3C;", "\u2b50 SFARI")
    else
      tags$span(class="badge bg-secondary ms-1", "No SFARI")
    
    tagList(
      div(class="mb-3",
          h5(class="mb-1 d-flex align-items-center",
             tags$span(style="font-size:1.2em; font-weight:700;", gen_info$gen_solo),
             sfari_tag),
          tags$small(class="text-muted",
                     paste0("chr", gen_info$chr_clean, " \u00b7 ",
                            format(round(gen_info$pos_media), big.mark=".", scientific=FALSE), " pb"))
      ),
      layout_columns(
        col_widths = c(6, 6),
        div(class="dm-stat-blue rounded p-2 text-center",
            tags$div(style="font-size:1.8em; font-weight:700; color:#2C6FAC;", gen_info$N_familias),
            tags$small("familia(s)")),
        div(class="dm-stat-gray rounded p-2 text-center",
            tags$div(class="dm-muted-txt fw-bold", style="font-size:1.8em;", gen_info$N_variantes),
            tags$small("variante(s)"))
      ),
      div(class="mt-2 mb-1",
          tags$b("Herencias: "),
          tags$span(class="text-muted small", gen_info$Herencias)),
      div(class="mb-2",
          tags$b("Score medio: "),
          tags$span(class="text-muted small",
                    if (is.nan(gen_info$Score_medio)) "\u2014" else as.character(gen_info$Score_medio))),
      hr(class="my-2"),
      tags$b("Familias afectadas:", class="small"),
      div(class="mt-1",
          lapply(strsplit(gen_info$Familias, ", ")[[1]], function(fam) {
            tags$span(class="badge me-1 mb-1",
                      style="background:#2C6FAC; font-size:0.82em;", fam)
          })
      ),
      hr(class="my-2"),
      tags$b("Variantes:", class="small"),
      div(class="mt-1",
          DTOutput("dt_gen_detalle_filas", height="280px"))
    )
  })
  
  output$dt_gen_detalle_filas <- renderDT({
    gen_name <- gen_seleccionado()
    if (is.null(gen_name)) return(NULL)
    df_full <- datos_filtrados()
    if (!"Gene_name" %in% names(df_full)) return(NULL)
    gen_col  <- as.character(df_full$Gene_name)
    n_g      <- lengths(strsplit(gen_col, "[/|,]"))
    df_exp   <- df_full[rep(seq_len(nrow(df_full)), n_g), , drop=FALSE]
    df_exp$gene_solo <- trimws(unlist(strsplit(gen_col, "[/|,]")))
    df_rows  <- df_exp[df_exp$gene_solo == gen_name, , drop=FALSE]
    cols_det <- intersect(c("ID_Familia","SV_chrom","SV_start","SV_end","SV_type",
                            "Tipo_Herencia","Tipo_Rango","AnnotSV_ranking_score"), names(df_rows))
    
    # ── LECTURA DEL MODO OSCURO ──
    modo_oscuro <- isTRUE(input$modo_noche == "dark")
    txt_color   <- if (modo_oscuro) "white" else NULL
    bg_her      <- if (modo_oscuro) c("#CC0000","#0066CC","#CC6600","#6600CC","#6c757d") else c("#FFE6E6","#E6F3FF","#FFF0E6","#F0E6FF","#E8E8E8")
    
    datatable(df_rows[, cols_det, drop=FALSE],
              rownames=FALSE,
              options=list(dom="t", pageLength=20, scrollX=TRUE),
              class="table table-sm table-hover") |>
      formatStyle("Tipo_Herencia",
                  backgroundColor = styleEqual(c("De novo","Paternal","Maternal","Combined","Desconocido"), bg_her),
                  color = txt_color, fontWeight = styleEqual("De novo", "bold"))
  }, server=TRUE)
  
  # ── Historial de acciones ──────────────────────────────────────────────────
  output$tabla_historial <- renderDT({
    hist <- rv$historial
    if (length(hist) == 0)
      return(datatable(data.frame(Mensaje = "No actions recorded in this session."),
                       options = list(dom = "t"), rownames = FALSE))
    df_h <- do.call(rbind, lapply(rev(hist), function(x) {
      data.frame(
        Timestamp = x$ts,
        Accion    = x$tipo,
        Variante  = x$clave,
        Detalle   = if (!is.null(x$detalle)) x$detalle else "",
        stringsAsFactors = FALSE
      )
    }))
    datatable(df_h, rownames = FALSE, options = list(
      dom = "frtip", pageLength = 25, order = list(list(0, "desc")),
      columnDefs = list(list(width = "160px", targets = 0), list(width = "160px", targets = 1))
    ), class = "table table-hover table-sm table-striped")
  })
  
  observeEvent(input$btn_limpiar_historial, {
    rv$historial <- list()
    showNotification("Historial limpiado.", type = "warning", duration = 2)
  })
  
  output$plot_ideograma <- renderPlotly({
    df <- datos_filtrados()
    req(nrow(df) > 0)
    
    if ("Annotation_mode" %in% names(df))
      df <- df[!is.na(df$Annotation_mode) & df$Annotation_mode == "full", , drop = FALSE]
    req(nrow(df) > 0)
    
    df$chr_clean <- toupper(trimws(gsub("(?i)^chr", "", as.character(df$SV_chrom), perl = TRUE)))
    df$start_num <- suppressWarnings(as.numeric(df$SV_start))
    df$end_num   <- suppressWarnings(as.numeric(df$SV_end))
    df$her_base  <- her_base(df$Tipo_Herencia)
    df$score_lbl <- ifelse(is.na(df$AnnotSV_ranking_score), "—", as.character(df$AnnotSV_ranking_score))
    df <- df[df$chr_clean %in% CHR_ORDER & !is.na(df$start_num) & !is.na(df$end_num), , drop = FALSE]
    req(nrow(df) > 0)
    
    # Vertical index per chromosome (Y axis)
    chr_idx <- setNames(seq_along(CHR_ORDER), CHR_ORDER)
    
    fig <- plot_ly()
    
    # Barras de cromosomas (fondo gris)
    for (chr in CHR_ORDER) {
      if (!chr %in% names(CHR_LENGTHS)) next
      fig <- add_trace(fig,
                       type = "bar", orientation = "h",
                       x = CHR_LENGTHS[chr], y = chr_idx[chr],
                       base = 0,
                       width = 0.55,
                       marker = list(color = "#E0E0E0", line = list(color = "#BDBDBD", width = 0.8)),
                       hoverinfo = "none", showlegend = FALSE
      )
    }
    
    # Variantes coloreadas por herencia
    her_levels <- setdiff(unique(df$her_base), NA)
    for (her in her_levels) {
      df_h <- df[!is.na(df$her_base) & df$her_base == her, , drop = FALSE]
      if (nrow(df_h) == 0) next
      col <- HER_COLORS[her]; if (is.na(col)) col <- "#AAAAAA"
      
      # Minimum visible width: 2% of chromosome to not lose small variants
      widths <- pmax(df_h$end_num - df_h$start_num,
                     CHR_LENGTHS[df_h$chr_clean] * 0.015, na.rm = TRUE)
      
      fig <- add_trace(fig,
                       type = "bar", orientation = "h",
                       x    = widths,
                       y    = chr_idx[df_h$chr_clean],
                       base = df_h$start_num,
                       width = 0.55,
                       marker = list(
                         color = paste0(col, "BB"),
                         line  = list(color = col, width = 1.2)
                       ),
                       name          = her,
                       legendgroup   = her,
                       showlegend    = TRUE,
                       text = paste0(
                         "<b>Familia:</b> ",   df_h$ID_Familia,   "<br>",
                         "<b>Region:</b> chr", df_h$chr_clean, ":", format(df_h$start_num, big.mark=".", scientific=FALSE),
                         "–", format(df_h$end_num, big.mark=".", scientific=FALSE), "<br>",
                         "<b>Tipo:</b> ",      df_h$SV_type,       "<br>",
                         "<b>Herencia:</b> ",  her,                 "<br>",
                         "<b>Score:</b> ",     df_h$score_lbl
                       ),
                       hovertemplate = "%{text}<extra></extra>"
      )
    }
    
    layout(fig,
           barmode = "overlay",
           xaxis = list(
             title      = "Genomic position (bp)",
             tickformat = ",.0f",
             showgrid   = TRUE,
             gridcolor  = "#F0F0F0"
           ),
           yaxis = list(
             title      = "",
             tickvals   = seq_along(CHR_ORDER),
             ticktext   = paste0("chr", CHR_ORDER),
             autorange  = "reversed",
             tickfont   = list(size = 11)
           ),
           legend = list(
             title       = list(text = "<b>Herencia</b>"),
             orientation = "h",
             y           = -0.12,
             x           = 0,
             font        = list(size = 11)
           ),
           margin  = list(l = 65, r = 20, t = 15, b = 70),
           plot_bgcolor  = "#FAFAFA",
           paper_bgcolor = "#FFFFFF",
           height  = 720
    )
  })
  
  
  # =============================================================================
  # STATISTICAL TESTS — server logic
  # =============================================================================
  
  # Reactive: clean dataset for tests (without Results tab filters)
  test_datos_base <- reactive({
    datos <- rv$resultados
    if (is.null(datos$CNVs) && is.null(datos$SVs)) return(data.frame())
    cnv <- if (!is.null(datos$CNVs)) mutate(datos$CNVs, Modalidad="CNV") else data.frame()
    sv  <- if (!is.null(datos$SVs))  mutate(datos$SVs,  Modalidad="SV")  else data.frame()
    df  <- bind_rows(cnv, sv)
    if (nrow(df) == 0) return(df)
    if ("Annotation_mode" %in% names(df))
      df <- df[!is.na(df$Annotation_mode) & df$Annotation_mode == "full", , drop=FALSE]
    df$her_limpia  <- her_base(df$Tipo_Herencia)
    df$score_num   <- suppressWarnings(as.numeric(df$AnnotSV_ranking_score))
    df$start_num   <- suppressWarnings(as.numeric(df$SV_start))
    df$end_num     <- suppressWarnings(as.numeric(df$SV_end))
    df$tamano_kb   <- pmax(0, df$end_num - df$start_num + 1) / 1000
    df$chr_clean   <- toupper(trimws(gsub("(?i)^chr","",as.character(df$SV_chrom),perl=TRUE)))
    sv_up          <- toupper(as.character(df$SV_type))
    df$tipo_agrup  <- ifelse(grepl("DEL|LOSS",sv_up),"DEL",
                             ifelse(grepl("DUP|GAIN",sv_up),"DUP","OTHER"))
    df$en_sfari    <- if ("En_Panel_Genes" %in% names(df))
      es_sfari_col(df$En_Panel_Genes) else FALSE
    # Numeric encodings of categorical variables (for numeric tests)
    df$her_limpia_num  <- as.integer(factor(df$her_limpia))
    df$tipo_agrup_num  <- as.integer(factor(df$tipo_agrup, levels=c("DEL","DUP","OTHER")))
    df$chr_clean_num   <- match(df$chr_clean, CHR_ORDER)
    df$Modalidad_num   <- as.integer(factor(df$Modalidad))
    df$en_sfari_num    <- as.integer(df$en_sfari)
    df
  })
  
  # Variables continuas disponibles
  VARS_CONTINUAS <- c(
    "Score AnnotSV"   = "score_num",
    "Size (Kb)"     = "tamano_kb"
  )
  # Available categorical variables
  VARS_CATEG <- c(
    "Herencia"        = "her_limpia",
    "Tipo variante"   = "tipo_agrup",
    "Cromosoma"       = "chr_clean",
    "Modalidad (CNV/SV)" = "Modalidad",
    "Gen SFARI"       = "en_sfari"
  )
  # Categorical variables encoded as integers (for numeric tests)
  VARS_CATEG_AS_NUM <- c(
    "Herencia (\u2192 n\u00ba)"          = "her_limpia_num",
    "Tipo variante (\u2192 n\u00ba)"     = "tipo_agrup_num",
    "Cromosoma (\u2192 n\u00ba)"         = "chr_clean_num",
    "Modalidad (\u2192 n\u00ba)"         = "Modalidad_num",
    "Gen SFARI (\u2192 n\u00ba)"         = "en_sfari_num"
  )
  
  # =============================================================================
  # WIZARD GUIADO — tarjetas de pregunta en lenguaje natural
  # =============================================================================
  
  # Wizard question definitions (ordered by frequency of use)
  WIZARD_PREGUNTAS <- list(
    list(
      id     = "wiz_wilcoxon",  tipo  = "wilcoxon",
      icono  = "\U0001f4ca",   color = "#2C6FAC",
      titulo = "\u00bfDifieren dos grupos?",
      desc   = "Compara Score o Tama\u00f1o entre dos tipos de herencia, DEL vs DUP, SFARI vs no\u2026",
      tip    = "El m\u00e1s usado en estudios de CNV"
    ),
    list(
      id     = "wiz_kruskal",   tipo  = "kruskal",
      icono  = "\U0001f4ca",   color = "#5B6CBD",
      titulo = "\u00bfDifieren varios grupos?",
      desc   = "Compara Score o Tama\u00f1o entre todos los tipos de herencia a la vez",
      tip    = "Generaliza Wilcoxon a 3+ grupos"
    ),
    list(
      id     = "wiz_chisq",     tipo  = "chisq",
      icono  = "\U0001f517",   color = "#27AE60",
      titulo = "\u00bfSe asocian dos variables categ\u00f3ricas?",
      desc   = "Ej: \u00bfla herencia y el tipo de variante (DEL/DUP) van juntos? \u00bfSFARI y herencia?",
      tip    = "Chi\u00b2 si n grande \u00b7 Fisher si hay celdas peque\u00f1as"
    ),
    list(
      id     = "wiz_spearman",  tipo  = "spearman",
      icono  = "\U0001f4c8",   color = "#E67E22",
      titulo = "\u00bfCorrelacionan dos variables num\u00e9ricas?",
      desc   = "Ej: \u00bflas variantes m\u00e1s grandes tienen mayor Score AnnotSV?",
      tip    = "No asume normalidad \u00b7 robusto a valores extremos"
    ),
    list(
      id     = "wiz_burden",    tipo  = "burden",
      icono  = "\U0001f9ec",   color = "#8E44AD",
      titulo = "\u00bfTienen mayor carga gen\u00f3mica ciertos grupos?",
      desc   = "Compara el n\u00ba de variantes por familia: SFARI vs no, De novo vs heredada\u2026",
      tip    = "Est\u00e1ndar en estudios de CNV en TEA"
    ),
    list(
      id     = "wiz_prop_binom", tipo = "prop_binom",
      icono  = "\U0001f3af",   color = "#CC0000",
      titulo = "\u00bfLa proporci\u00f3n observada es la esperada?",
      desc   = "Ej: \u00bfla tasa de De novo (~15-20\u0025 en TEA) concuerda con esta cohorte?",
      tip    = "Contrasta contra un valor de referencia"
    ),
    list(
      id     = "wiz_shapiro",   tipo  = "shapiro",
      icono  = "\U0001f4d0",   color = "#7F8C8D",
      titulo = "\u00bfSiguen los datos una distribuci\u00f3n normal?",
      desc   = "Paso previo para saber si usar tests param\u00e9tricos (t-test) o no param\u00e9tricos (Wilcoxon)",
      tip    = "Usar siempre antes de Pearson o ANOVA"
    )
  )
  
  output$ui_test_wizard_cards <- renderUI({
    if (!isTRUE(input$test_modo_ui == "guiado")) return(NULL)
    df <- test_datos_base()
    
    if (nrow(df) == 0)
      return(div(class = "alert alert-warning small p-2 mt-1",
                 bsicons::bs_icon("exclamation-triangle"), " Carga resultados primero."))
    
    tipo_actual <- isolate(input$test_tipo) %||% "wilcoxon"
    n_vars      <- nrow(df)
    
    tagList(
      # ── Cabecera ──────────────────────────────────────────────────────────
      div(class = "mb-2",
          tags$p(class = "fw-bold mb-1", style = "font-size:0.88em; color:#2C6FAC;",
                 "\U0001f9ed \u00bfQu\u00e9 quieres investigar?"),
          tags$p(class = "text-muted mb-0", style = "font-size:0.78em;",
                 paste0("Datos disponibles: ", n_vars, " variantes"))
      ),
      
      # ── Tarjetas de pregunta ──────────────────────────────────────────────
      lapply(WIZARD_PREGUNTAS, function(q) {
        activo <- identical(tipo_actual, q$tipo)
        div(
          class = "mb-1",
          tags$button(
            id        = paste0("btn_", q$id),
            class     = paste0(
              "btn w-100 text-start p-2 ",
              if (activo) "btn-primary shadow-sm" else "btn-outline-secondary wiz-card-inactive"
            ),
            style     = paste0(
              "border-radius:8px; line-height:1.3; font-size:0.82em;",
              if (activo) paste0(" border-color:", q$color, "; background-color:", q$color, ";")
              else        paste0(" border-color:#dee2e6;")
            ),
            onclick   = paste0(
              "Shiny.setInputValue(\'wiz_sel\', \'", q$tipo, "\',",
              "{priority:\'event\'});"
            ),
            div(class = "d-flex align-items-start gap-2",
                div(style = "font-size:1.25em; line-height:1; min-width:1.5em; text-align:center;",
                    q$icono),
                div(
                  div(class = if (activo) "fw-bold text-white" else "fw-bold text-dark",
                      style = "font-size:0.95em;", q$titulo),
                  div(class = if (activo) "text-white opacity-75" else "text-muted",
                      style = "font-size:0.88em; white-space:normal;", q$desc),
                  if (activo)
                    div(class = "mt-1",
                        tags$span(class = "badge",
                                  style = "background:rgba(255,255,255,0.25); font-size:0.75em;",
                                  paste0("\u2139\ufe0f ", q$tip)))
                  else
                    div(class = "text-muted", style = "font-size:0.75em; margin-top:2px;",
                        paste0("\u2139\ufe0f ", q$tip))
                )
            )
          )
        )
      }),
      
      # ── Enlace a modo experto ─────────────────────────────────────────────
      div(class = "mt-2 text-center",
          tags$small(
            class = "text-muted",
            "Tests avanzados (Fisher, KS, Pearson\u2026): ",
            tags$a(
              href    = "#",
              style   = "color:#2C6FAC; font-size:0.85em;",
              onclick = paste0(
                "Shiny.setInputValue(\'test_modo_ui_js\',\'experto\',",
                "{priority:\'event\'});return false;"
              ),
              "cambiar a modo experto"
            )
          )
      )
    )
  })
  
  # ── Observer: click en tarjeta del wizard ─────────────────────────────────
  observeEvent(input$wiz_sel, {
    req(input$wiz_sel)
    updateSelectInput(session, "test_tipo", selected = input$wiz_sel)
  }, ignoreNULL = TRUE)
  
  # ── Observer: enlace "modo experto" desde el wizard ───────────────────────
  observeEvent(input$test_modo_ui_js, {
    req(input$test_modo_ui_js == "experto")
    updateRadioButtons(session, "test_modo_ui", selected = "experto")
  }, ignoreNULL = TRUE)
  
  # ── Dynamic parameter UI based on selected test ────────────────────────
  output$ui_test_params <- renderUI({
    tipo <- input$test_tipo
    df   <- test_datos_base()
    if (nrow(df) == 0)
      return(div(class="alert alert-warning small", "Carga resultados primero."))
    
    usar_categ_num <- isTRUE(input$test_categ_como_num)
    vars_num <- if (usar_categ_num) c(VARS_CONTINUAS, VARS_CATEG_AS_NUM) else VARS_CONTINUAS
    
    her_vals  <- sort(unique(na.omit(df$her_limpia)))
    tipo_vals <- sort(unique(na.omit(df$tipo_agrup)))
    chr_vals  <- CHR_ORDER[CHR_ORDER %in% unique(df$chr_clean)]
    mod_vals  <- sort(unique(na.omit(as.character(df$Modalidad))))
    
    opciones_grp <- function(col) {
      switch(col,
             her_limpia  = her_vals,
             tipo_agrup  = tipo_vals,
             chr_clean   = chr_vals,
             Modalidad   = mod_vals,
             en_sfari    = c("TRUE","FALSE"),
             character(0)
      )
    }
    
    if (tipo %in% c("wilcoxon","ks","prop2")) {
      tagList(
        div(class="mb-2",
            checkboxInput("test_categ_como_num",
                          tags$span(bsicons::bs_icon("123"), " Include factors as numerics"),
                          value=isTRUE(input$test_categ_como_num))),
        selectInput("test_var_y", "Numeric variable",
                    choices=vars_num, selected=isolate(input$test_var_y) %||% "score_num", width="100%"),
        selectInput("test_var_grupo", "Variable de grupo",
                    choices=VARS_CATEG, selected="her_limpia", width="100%"),
        uiOutput("ui_test_grupo_ab")
      )
    } else if (tipo == "kruskal") {
      tagList(
        div(class="mb-2",
            checkboxInput("test_categ_como_num",
                          tags$span(bsicons::bs_icon("123"), " Include factors as numerics"),
                          value=isTRUE(input$test_categ_como_num))),
        selectInput("test_var_y", "Numeric variable",
                    choices=vars_num, selected=isolate(input$test_var_y) %||% "score_num", width="100%"),
        selectInput("test_var_grupo", "Variable de grupo",
                    choices=VARS_CATEG, selected="her_limpia", width="100%"),
        checkboxGroupInput("test_grupos_multi", "Grupos a incluir",
                           choices  = her_vals,
                           selected = her_vals,
                           inline   = FALSE
        )
      )
    } else if (tipo %in% c("chisq","fisher")) {
      tagList(
        selectInput("test_var_fila", "Variable filas",
                    choices=VARS_CATEG, selected="her_limpia", width="100%"),
        selectInput("test_var_col", "Variable columnas",
                    choices=VARS_CATEG[VARS_CATEG != "her_limpia"], selected="tipo_agrup", width="100%"),
        numericInput("test_chisq_minexp", "Minimum expected frequency (cell filter)",
                     value=5, min=1, step=1, width="100%")
      )
    } else if (tipo %in% c("spearman","pearson")) {
      tagList(
        div(class="mb-2",
            checkboxInput("test_categ_como_num",
                          tags$span(bsicons::bs_icon("123"), " Include factors as numerics"),
                          value=isTRUE(input$test_categ_como_num))),
        selectInput("test_var_x", "Variable X",
                    choices=vars_num, selected=isolate(input$test_var_x) %||% "tamano_kb", width="100%"),
        selectInput("test_var_y2", "Variable Y",
                    choices=vars_num, selected=isolate(input$test_var_y2) %||% "score_num", width="100%"),
        checkboxInput("test_corr_log_x", "Log10 on X (size)", value=TRUE),
        hr(),
        tags$small(class="text-muted",
                   "Pairs with NA or infinite values are automatically excluded.")
      )
    } else if (tipo == "shapiro") {
      tagList(
        div(class="mb-2",
            checkboxInput("test_categ_como_num",
                          tags$span(bsicons::bs_icon("123"), " Include factors as numerics"),
                          value=isTRUE(input$test_categ_como_num))),
        selectInput("test_var_shap", "Variable a testear",
                    choices=vars_num, selected=isolate(input$test_var_shap) %||% "score_num", width="100%"),
        selectInput("test_shap_grupo_var", "Separar por grupo (opcional)",
                    choices=c("Ninguno"="ninguno", VARS_CATEG), selected="ninguno", width="100%"),
        tags$small(class="text-muted",
                   "Shapiro-Wilk requiere n entre 3 y 5000. Grupos con n fuera de rango se omiten.")
      )
    } else if (tipo == "prop_binom") {
      tagList(
        selectInput("test_prop_var", "Success variable",
                    choices=VARS_CATEG, selected="en_sfari", width="100%"),
        uiOutput("ui_test_prop_valor"),
        numericInput("test_prop_p0", "Expected proportion (H0)",
                     value=0.10, min=0.001, max=0.999, step=0.01, width="100%"),
        tags$small(class="text-muted",
                   "Tests whether the observed proportion differs significantly from p0.")
      )
    } else if (tipo == "burden") {
      tagList(
        selectInput("test_burden_grupo", "Comparison group",
                    choices=VARS_CATEG, selected="en_sfari", width="100%"),
        uiOutput("ui_test_burden_niveles"),
        tags$small(class="text-muted",
                   "Compares the mean number of variants per family between groups (Wilcoxon on counts).")
      )
    }
  })
  
  # Sub-UI: selection of groups A and B for 2-group tests
  output$ui_test_grupo_ab <- renderUI({
    grp_var <- input$test_var_grupo
    req(!is.null(grp_var))
    df <- test_datos_base(); req(nrow(df) > 0)
    vals <- sort(unique(na.omit(as.character(df[[grp_var]]))))
    tagList(
      selectInput("test_grupo_a", "Grupo A",
                  choices=vals, selected=vals[1], width="100%"),
      selectInput("test_grupo_b", "Grupo B",
                  choices=vals, selected=if(length(vals)>=2) vals[2] else vals[1], width="100%")
    )
  })
  
  # Sub-UI: success value for prop_binom
  output$ui_test_prop_valor <- renderUI({
    var <- input$test_prop_var; req(!is.null(var))
    df  <- test_datos_base(); req(nrow(df) > 0)
    vals <- sort(unique(na.omit(as.character(df[[var]]))))
    selectInput("test_prop_exito", "Value considered 'success'",
                choices=vals, selected=vals[1], width="100%")
  })
  
  # Sub-UI: niveles para burden
  output$ui_test_burden_niveles <- renderUI({
    var <- input$test_burden_grupo; req(!is.null(var))
    df  <- test_datos_base(); req(nrow(df) > 0)
    vals <- sort(unique(na.omit(as.character(df[[var]]))))
    tagList(
      selectInput("test_burden_a", "Grupo A", choices=vals, selected=vals[1], width="100%"),
      selectInput("test_burden_b", "Grupo B",
                  choices=vals, selected=if(length(vals)>=2) vals[2] else vals[1], width="100%")
    )
  })
  
  # ── Ayuda contextual por test ─────────────────────────────────────────────────
  AYUDAS <- list(
    wilcoxon  = list(
      nombre = "Test de Wilcoxon-Mann-Whitney",
      cuando = "Compara medianas de una variable continua entre dos grupos independientes sin asumir normalidad.",
      h0     = "H0: The distributions of both groups are identical (equal medians).",
      uso    = "Ideal for comparing scores or sizes across inheritances (De novo vs Paternal, DEL vs DUP...).",
      ref    = "Mann & Whitney (1947) · No paramétrico · Alternativa al t-test de Student."
    ),
    kruskal   = list(
      nombre = "Test de Kruskal-Wallis",
      cuando = "Generalisation of Wilcoxon for 3 or more independent groups.",
      h0     = "H0: The distributions of all groups are identical.",
      uso    = "Compare scores across all inheritance types simultaneously.",
      ref    = "Kruskal & Wallis (1952) · No paramétrico · Alternativa al ANOVA de 1 factor."
    ),
    chisq     = list(
      nombre = "Chi-cuadrado de independencia",
      cuando = "Detects association between two categorical variables with large samples.",
      h0     = "H0: The two variables are independent (no association).",
      uso    = "Asociar tipo de herencia con tipo de variante (DEL/DUP), cromosoma, etc.",
      ref    = "Pearson (1900) · Requiere frecuencias esperadas >= 5 por celda."
    ),
    fisher    = list(
      nombre = "Test exacto de Fisher",
      cuando = "Association between two binary variables or tables with small cells (n < 20).",
      h0     = "H0: Las dos variables son independientes.",
      uso    = "More precise than chi-squared when there are cells with few observations.",
      ref    = "Fisher (1922) · Exacto (no asintótico) · Recomendado para tamaños muestrales pequeños."
    ),
    spearman  = list(
      nombre = "Correlacion de Spearman",
      cuando = "Mide la asociacion monotona entre dos variables continuas sin asumir normalidad.",
      h0     = "H0: rho = 0 (no hay correlacion monotona).",
      uso    = "Correlacion entre tamanio de variante y score AnnotSV.",
      ref    = "Spearman (1904) · No parametrico · Robusto a valores extremos."
    ),
    pearson   = list(
      nombre = "Correlacion de Pearson",
      cuando = "Mide la asociacion lineal entre dos variables continuas normalmente distribuidas.",
      h0     = "H0: r = 0 (no hay correlacion lineal).",
      uso    = "Usar solo si ambas variables son aproximadamente normales (verificar con Shapiro-Wilk).",
      ref    = "Pearson (1895) · Parametrico · Sensible a valores extremos."
    ),
    shapiro   = list(
      nombre = "Test de normalidad de Shapiro-Wilk",
      cuando = "Contrasta si una muestra procede de una distribucion normal.",
      h0     = "H0: Los datos siguen una distribucion normal.",
      uso    = "Paso previo para decidir entre tests parametricos (t-test, ANOVA) y no parametricos (Wilcoxon, Kruskal).",
      ref    = "Shapiro & Wilk (1965) · Potente para n < 50 · Para n grande, usar QQ-plot."
    ),
    ks        = list(
      nombre = "Test de Kolmogorov-Smirnov (2 muestras)",
      cuando = "Compara la distribucion completa de dos muestras (no solo la mediana).",
      h0     = "H0: Ambas muestras proceden de la misma distribucion continua.",
      uso    = "Detectar si DEL y DUP tienen distribuciones de tamano o score distintas en su forma global.",
      ref    = "Smirnov (1948) · Sensible a diferencias en forma, media y varianza."
    ),
    prop_binom = list(
      nombre = "Test binomial de proporciones",
      cuando = "Tests whether the observed proportion of a category differs from a reference value.",
      h0     = "H0: La proporcion observada es igual a p0 (valor esperado teorico).",
      uso    = "Ejemplo: ¿la tasa de variantes de novo (~15-20% esperada en TEA) es distinta en esta cohorte?",
      ref    = "Test binomial exacto · No asume normalidad · Valido para cualquier n."
    ),
    prop2     = list(
      nombre = "Comparacion de proporciones (2 grupos)",
      cuando = "Contrasta si dos grupos tienen la misma proporcion de una caracteristica.",
      h0     = "H0: La proporcion del evento es igual en ambos grupos.",
      uso    = "Ejemplo: ¿los varones tienen mayor tasa de variantes de novo que las mujeres?",
      ref    = "prop.test de R (correccion de Yates) · Equivalente al chi-cuadrado 2x2."
    ),
    burden    = list(
      nombre = "Analisis de carga genomica (Burden)",
      cuando = "Compara el numero de variantes por familia entre dos subgrupos.",
      h0     = "H0: La carga media de variantes es igual en ambos grupos.",
      uso    = "Ejemplo: ¿Las familias con variante de novo tienen mas variantes en total? ¿Las SFARI vs no SFARI?",
      ref    = "Wilcoxon sobre conteos por individuo · Robusto · Estandar en estudios de CNV en TEA."
    )
  )
  
  output$ui_test_ayuda <- renderUI({
    tipo <- input$test_tipo; req(!is.null(tipo))
    ay   <- AYUDAS[[tipo]]; if (is.null(ay)) return(NULL)
    
    # p-value interpretation guide (always visible)
    guia_p <- div(
      class = "mt-2 p-2 rounded",
      class = "dm-warn-panel rounded mt-2 p-2", style="font-size:0.80em; border:1px solid;",
      tags$b(class="dm-dark-txt", "\U0001f4d6 C\u00f3mo interpretar el p-valor:"),
      div(class="mt-1 d-flex flex-column gap-1",
          div(class="d-flex align-items-center gap-2",
              tags$span(class="badge bg-danger", "p < 0.001"),
              tags$span("Diferencia muy significativa (***)")),
          div(class="d-flex align-items-center gap-2",
              tags$span(class="badge bg-warning text-dark", "p < 0.01"),
              tags$span("Diferencia significativa (**)")),
          div(class="d-flex align-items-center gap-2",
              tags$span(class="badge bg-primary", "p < 0.05"),
              tags$span("Diferencia marginalmente significativa (*)")),
          div(class="d-flex align-items-center gap-2",
              tags$span(class="badge bg-secondary", "p \u2265 0.05"),
              tags$span("Sin evidencia suficiente de diferencia (ns)"))
      ),
      div(class="mt-1 text-muted", style="font-size:0.9em;",
          "\u26a0\ufe0f Un p-valor significativo indica que la diferencia observada no es atribuible al azar,",
          " but does not imply causality or clinical relevance.")
    )
    
    tagList(
      div(
        class = "p-2 rounded mt-1",
        class = "dm-info-panel rounded p-2 mt-1", style="font-size:0.83em;",
        div(class="d-flex align-items-center gap-1 mb-1",
            bsicons::bs_icon("info-circle-fill", class="text-primary"),
            tags$b(style="color:#2C6FAC;", ay$nombre)),
        tags$span(class="text-muted", ay$cuando), br(), br(),
        div(class="p-1 rounded mb-1",
            class="dm-h0-panel rounded p-1 mb-1", style="border-left:3px solid;",
            tags$b("H\u2080: "), tags$span(class="dm-muted-txt", ay$h0)),
        div(class="p-1 rounded",
            class="dm-uso-panel rounded p-1", style="border-left:3px solid;",
            tags$b("\U0001f4a1 T\u00edpico: "), tags$span(class="dm-muted-txt", ay$uso)),
        br(),
        tags$i(class="text-muted", style="font-size:0.9em;", ay$ref)
      ),
      guia_p
    )
  })
  
  # ── Reactive: runs test when button is pressed ──────────────────────────────
  test_resultado <- eventReactive(input$btn_ejecutar_test, {
    tipo <- input$test_tipo
    df   <- test_datos_base()
    if (nrow(df) == 0) return(list(error="Sin datos cargados."))
    
    tryCatch({
      
      if (tipo == "wilcoxon") {
        var_y <- input$test_var_y; ga <- input$test_grupo_a; gb <- input$test_grupo_b
        grp   <- input$test_var_grupo
        y_a   <- df[[var_y]][as.character(df[[grp]]) == ga & is.finite(df[[var_y]])]
        y_b   <- df[[var_y]][as.character(df[[grp]]) == gb & is.finite(df[[var_y]])]
        if (length(y_a) < 3 || length(y_b) < 3) return(list(error="Groups too small (n<3)."))
        res   <- wilcox.test(y_a, y_b, conf.int=TRUE)
        list(tipo=tipo, test=res,
             grupos=list(A=list(nombre=ga,vals=y_a),B=list(nombre=gb,vals=y_b)),
             var_y=var_y, grp=grp)
        
      } else if (tipo == "kruskal") {
        var_y  <- input$test_var_y; grp <- input$test_var_grupo
        grupos <- input$test_grupos_multi
        df2    <- df[as.character(df[[grp]]) %in% grupos & is.finite(df[[var_y]]),]
        if (nrow(df2) < 6) return(list(error="Insuficientes datos tras el filtro."))
        res <- kruskal.test(df2[[var_y]] ~ factor(df2[[grp]]))
        list(tipo=tipo, test=res, df2=df2, var_y=var_y, grp=grp, grupos=grupos)
        
      } else if (tipo %in% c("chisq","fisher")) {
        v1 <- input$test_var_fila; v2 <- input$test_var_col
        df2 <- df[!is.na(df[[v1]]) & !is.na(df[[v2]]),]
        tbl <- table(df2[[v1]], df2[[v2]])
        if (any(dim(tbl) < 2)) return(list(error="Alguna variable tiene menos de 2 niveles."))
        if (tipo == "chisq") {
          res <- chisq.test(tbl)
        } else {
          if (prod(dim(tbl)) > 100) return(list(error="Tabla demasiado grande para Fisher exacto. Usa Chi-cuadrado."))
          res <- fisher.test(tbl, simulate.p.value=(prod(dim(tbl))>4))
        }
        list(tipo=tipo, test=res, tabla=tbl, v1=v1, v2=v2)
        
      } else if (tipo %in% c("spearman","pearson")) {
        vx <- input$test_var_x; vy <- input$test_var_y2
        log_x <- isTRUE(input$test_corr_log_x) && vx == "tamano_kb"
        x <- df[[vx]]; y <- df[[vy]]
        if (log_x) x <- log10(x)
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < 5) return(list(error="Menos de 5 pares validos."))
        metodo <- if (tipo=="spearman") "spearman" else "pearson"
        res <- cor.test(x[ok], y[ok], method=metodo)
        list(tipo=tipo, test=res, x=x[ok], y=y[ok], vx=vx, vy=vy, log_x=log_x)
        
      } else if (tipo == "shapiro") {
        var  <- input$test_var_shap
        grp_var <- input$test_shap_grupo_var
        if (grp_var == "ninguno") {
          vals <- df[[var]][is.finite(df[[var]])]
          if (length(vals) < 3 || length(vals) > 5000)
            return(list(error=paste0("n=",length(vals)," fuera de rango [3,5000].")))
          res <- shapiro.test(vals)
          list(tipo=tipo, tests=list(list(grp="Todos",res=res,vals=vals)), var=var)
        } else {
          niveles <- sort(unique(na.omit(as.character(df[[grp_var]]))))
          results <- lapply(niveles, function(g) {
            v <- df[[var]][as.character(df[[grp_var]])==g & is.finite(df[[var]])]
            if (length(v) < 3 || length(v) > 5000) return(NULL)
            list(grp=g, res=shapiro.test(v), vals=v)
          })
          results <- Filter(Negate(is.null), results)
          if (length(results)==0) return(list(error="No group has n in [3,5000]."))
          list(tipo=tipo, tests=results, var=var, grp_var=grp_var)
        }
        
      } else if (tipo == "ks") {
        var_y <- input$test_var_y; ga <- input$test_grupo_a; gb <- input$test_grupo_b
        grp   <- input$test_var_grupo
        y_a   <- df[[var_y]][as.character(df[[grp]]) == ga & is.finite(df[[var_y]])]
        y_b   <- df[[var_y]][as.character(df[[grp]]) == gb & is.finite(df[[var_y]])]
        if (length(y_a) < 3 || length(y_b) < 3) return(list(error="Groups too small."))
        res   <- ks.test(y_a, y_b)
        list(tipo=tipo, test=res,
             grupos=list(A=list(nombre=ga,vals=y_a),B=list(nombre=gb,vals=y_b)),
             var_y=var_y, grp=grp)
        
      } else if (tipo == "prop_binom") {
        var <- input$test_prop_var; exito <- input$test_prop_exito; p0 <- input$test_prop_p0
        vals  <- as.character(df[[var]])
        x <- sum(!is.na(vals) & vals == exito)
        n <- sum(!is.na(vals))
        if (n < 5) return(list(error="Menos de 5 observaciones validas."))
        res <- binom.test(x, n, p=p0)
        list(tipo=tipo, test=res, x=x, n=n, p0=p0, exito=exito, var=var)
        
      } else if (tipo == "prop2") {
        grp <- input$test_var_grupo; var_y <- input$test_var_y
        ga  <- input$test_grupo_a;   gb    <- input$test_grupo_b
        # Success = value > global median (derived binary variable)
        med_global <- median(df[[var_y]], na.rm=TRUE)
        df2 <- df[as.character(df[[grp]]) %in% c(ga,gb) & is.finite(df[[var_y]]),]
        if (nrow(df2) < 6) return(list(error="Datos insuficientes."))
        exito <- df2[[var_y]] > med_global
        n_a <- sum(as.character(df2[[grp]])==ga);  x_a <- sum(exito[as.character(df2[[grp]])==ga])
        n_b <- sum(as.character(df2[[grp]])==gb);  x_b <- sum(exito[as.character(df2[[grp]])==gb])
        if (n_a<1 || n_b<1) return(list(error="Un grupo esta vacio."))
        res <- prop.test(c(x_a,x_b), c(n_a,n_b))
        list(tipo=tipo, test=res, ga=ga, gb=gb, n_a=n_a, n_b=n_b, x_a=x_a, x_b=x_b,
             med_global=med_global, var_y=var_y, grp=grp)
        
      } else if (tipo == "burden") {
        grp_var <- input$test_burden_grupo
        ga <- input$test_burden_a; gb <- input$test_burden_b
        conta <- function(g) {
          fams <- unique(df$ID_Familia[as.character(df[[grp_var]])==g])
          sapply(fams, function(f) sum(as.character(df[[grp_var]])==g & df$ID_Familia==f))
        }
        ca <- conta(ga); cb <- conta(gb)
        if (length(ca)<2 || length(cb)<2) return(list(error="Grupos con menos de 2 familias."))
        res <- wilcox.test(ca, cb, conf.int=TRUE)
        list(tipo=tipo, test=res,
             grupos=list(A=list(nombre=ga,vals=ca),B=list(nombre=gb,vals=cb)),
             grp_var=grp_var)
      }
      
    }, error = function(e) list(error=paste("Error al ejecutar el test:", e$message)))
  })
  
  # ── Formateo p-value con estrellas ────────────────────────────────────────────
  fmt_p <- function(p) {
    if (is.null(p) || is.na(p)) return("NA")
    estrellas <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
    paste0(formatC(p, digits=4, format="g"), " ", estrellas)
  }
  semaforo_p <- function(p) {
    if (is.null(p) || is.na(p)) return("bg-secondary")
    if (p < 0.001) "bg-danger" else if (p < 0.01) "bg-warning text-dark" else if (p < 0.05) "bg-primary" else "bg-secondary"
  }
  
  # ── Renderizado del resultado ─────────────────────────────────────────────────
  output$ui_test_resultado <- renderUI({
    res <- test_resultado()
    if (is.null(res)) return(div(class="alert alert-info", "Configura el test y pulsa ▶ Ejecutar."))
    if (!is.null(res$error)) return(div(class="alert alert-danger", res$error))
    
    tipo <- res$tipo
    
    # Construir tarjetas de resultado
    tarjeta <- function(titulo, valor, clase="bg-primary") {
      div(class=paste("card text-white h-100", clase),
          div(class="card-body p-3 text-center",
              div(style="font-size:0.75em; font-weight:600; text-transform:uppercase; opacity:0.85;", titulo),
              div(style="font-size:1.6em; font-weight:700; margin-top:4px;", valor)))
    }
    
    if (tipo %in% c("wilcoxon","ks","burden")) {
      tst <- res$test; p <- tst$p.value
      ga <- res$grupos$A$nombre; gb <- res$grupos$B$nombre
      med_a <- round(median(res$grupos$A$vals, na.rm=TRUE), 3)
      med_b <- round(median(res$grupos$B$vals, na.rm=TRUE), 3)
      efecto <- if (!is.null(tst$estimate)) round(as.numeric(tst$estimate), 3) else "—"
      ci_txt <- if (!is.null(tst$conf.int))
        paste0("[", round(tst$conf.int[1],3), ", ", round(tst$conf.int[2],3), "]") else "—"
      tagList(
        card(card_header(class="bg-success text-white fw-bold",
                         bsicons::bs_icon("check-circle"), " Resultados"),
             card_body(
               layout_columns(col_widths=c(3,3,3,3),
                              tarjeta("p-valor", fmt_p(p), semaforo_p(p)),
                              tarjeta("W statistic", round(tst$statistic,2), "bg-dark"),
                              tarjeta(paste0("Mediana ", ga), med_a, "bg-info"),
                              tarjeta(paste0("Mediana ", gb), med_b, "bg-info")
               ),
               div(class="dm-result-panel rounded border p-3 mt-3",
                   tags$b("Location difference estimator: "), efecto, br(),
                   tags$b("Intervalo de confianza 95%: "), ci_txt, br(), br(),
                   tags$b("Interpretation: "),
                   if (p < 0.05)
                     tags$span(style="color:#CC0000; font-weight:600;",
                               paste0("Se rechaza H0 (p=", formatC(p,digits=3,format="g"), "). ",
                                      "There are statistically significant differences between '",
                                      ga, "' y '", gb, "'."))
                   else
                     tags$span(class="dm-muted-txt",
                               paste0("No se rechaza H0 (p=", formatC(p,digits=3,format="g"), "). ",
                                      "No hay evidencia suficiente de diferencias entre '",
                                      ga, "' y '", gb, "'."))
               )
             )
        )
      )
      
    } else if (tipo == "kruskal") {
      tst <- res$test; p <- tst$p.value
      meds <- sapply(res$grupos, function(g) {
        vals <- res$df2[[res$var_y]][as.character(res$df2[[res$grp]])==g & is.finite(res$df2[[res$var_y]])]
        round(median(vals, na.rm=TRUE), 3)
      })
      tagList(
        card(card_header(class="bg-success text-white fw-bold",
                         bsicons::bs_icon("check-circle"), " Resultados"),
             card_body(
               layout_columns(col_widths=c(4,4,4),
                              tarjeta("p-valor", fmt_p(p), semaforo_p(p)),
                              tarjeta("Chi² (H)", round(tst$statistic,2), "bg-dark"),
                              tarjeta("g.l.", tst$parameter, "bg-secondary")
               ),
               div(class="dm-result-panel rounded border p-3 mt-3",
                   tags$b("Medianas por grupo:"),
                   tags$ul(lapply(names(meds), function(g) tags$li(tags$b(g,":"), meds[[g]]))),
                   tags$b("Interpretation: "),
                   if (p < 0.05)
                     tags$span(style="color:#CC0000; font-weight:600;",
                               paste0("Se rechaza H0. Al menos un grupo difiere significativamente (p=",
                                      formatC(p,digits=3,format="g"), "). Considera post-hoc por pares (Wilcoxon + Bonferroni)."))
                   else
                     tags$span(class="dm-muted-txt",
                               paste0("No se rechaza H0 (p=", formatC(p,digits=3,format="g"),
                                      "). No hay diferencias significativas entre los grupos."))
               )
             )
        )
      )
      
    } else if (tipo %in% c("chisq","fisher")) {
      tst <- res$test; p <- tst$p.value
      tagList(
        card(card_header(class="bg-success text-white fw-bold",
                         bsicons::bs_icon("check-circle"), " Resultados"),
             card_body(
               layout_columns(col_widths=c(4,4,4),
                              tarjeta("p-valor", fmt_p(p), semaforo_p(p)),
                              tarjeta(if(tipo=="chisq") "Chi² stat" else "OR estimado",
                                      if(tipo=="chisq") round(tst$statistic,3)
                                      else if(!is.null(tst$estimate)) round(tst$estimate,3) else "—",
                                      "bg-dark"),
                              tarjeta("g.l.", if(tipo=="chisq") tst$parameter else "—", "bg-secondary")
               ),
               div(class="dm-result-panel rounded mt-3",
                   tags$b("Tabla de contingencia observada:"), br(),
                   renderTable(as.data.frame.matrix(res$tabla)),
                   br(), tags$b("Interpretation: "),
                   if (p < 0.05)
                     tags$span(style="color:#CC0000;font-weight:600;",
                               paste0("Significant association (p=", formatC(p,digits=3,format="g"),
                                      "). Las variables no son independientes."))
                   else
                     tags$span(class="dm-muted-txt",
                               paste0("No significant association (p=", formatC(p,digits=3,format="g"),
                                      "). Las variables parecen independientes."))
               )
             )
        )
      )
      
    } else if (tipo %in% c("spearman","pearson")) {
      tst <- res$test; p <- tst$p.value; rho <- round(tst$estimate, 4)
      ci_txt <- if (!is.null(tst$conf.int))
        paste0("[", round(tst$conf.int[1],3),", ",round(tst$conf.int[2],3),"]") else "—"
      fuerza <- if(abs(rho) >= 0.7) "strong" else if(abs(rho) >= 0.4) "moderate" else "weak"
      direc  <- if(rho > 0) "positiva" else "negativa"
      tagList(
        card(card_header(class="bg-success text-white fw-bold",
                         bsicons::bs_icon("check-circle"), " Resultados"),
             card_body(
               layout_columns(col_widths=c(4,4,4),
                              tarjeta("p-valor", fmt_p(p), semaforo_p(p)),
                              tarjeta(if(tipo=="spearman") "rho de Spearman" else "r de Pearson", rho,
                                      if(abs(rho)>=0.4)"bg-primary" else "bg-secondary"),
                              tarjeta("IC 95%", ci_txt, "bg-dark")
               ),
               div(class="dm-result-panel rounded mt-3",
                   tags$b("Interpretation: "),
                   if (p < 0.05)
                     tags$span(style="color:#2C6FAC; font-weight:600;",
                               paste0("Correlacion ", fuerza, " ", direc, " (rho=", rho, ", p=",
                                      formatC(p,digits=3,format="g"), ")."))
                   else
                     tags$span(class="dm-muted-txt",
                               paste0("No hay correlacion significativa (rho=", rho, ", p=",
                                      formatC(p,digits=3,format="g"), ")."))
               )
             )
        )
      )
      
    } else if (tipo == "shapiro") {
      resultados_sw <- res$tests
      tagList(
        card(card_header(class="bg-success text-white fw-bold",
                         bsicons::bs_icon("check-circle"), " Resultados"),
             card_body(
               DTOutput("tabla_shapiro_detalle"),
               div(class="dm-result-panel rounded mt-2 p-2", style="font-size:0.85em;",
                   tags$b("Interpretation guide: "),
                   "p > 0.05 → no se rechaza normalidad (distribución compatible con Normal). ",
                   "p ≤ 0.05 → se rechaza normalidad → usar tests no paramétricos (Wilcoxon, Kruskal-Wallis).")
             )
        )
      )
      
    } else if (tipo == "prop_binom") {
      tst <- res$test; p <- tst$p.value
      p_obs <- round(res$x / res$n, 4)
      ci_txt <- paste0("[", round(tst$conf.int[1],3), ", ", round(tst$conf.int[2],3), "]")
      tagList(
        card(card_header(class="bg-success text-white fw-bold",
                         bsicons::bs_icon("check-circle"), " Resultados"),
             card_body(
               layout_columns(col_widths=c(3,3,3,3),
                              tarjeta("p-valor", fmt_p(p), semaforo_p(p)),
                              tarjeta("Observed proportion", p_obs, "bg-primary"),
                              tarjeta("H0 proportion", res$p0, "bg-secondary"),
                              tarjeta("IC 95%", ci_txt, "bg-dark")
               ),
               div(class="dm-result-panel rounded mt-3",
                   tags$b(paste0("Successes: ", res$x, " of ", res$n, " (", round(p_obs*100,1), "%)")), br(),
                   tags$b("Interpretation: "),
                   if (p < 0.05)
                     tags$span(style="color:#CC0000;font-weight:600;",
                               paste0("La proporcion observada (", p_obs, ") difiere significativamente de p0=",
                                      res$p0, " (p=", formatC(p,digits=3,format="g"), ")."))
                   else
                     tags$span(class="dm-muted-txt",
                               paste0("No hay diferencia significativa respecto a p0=", res$p0,
                                      " (p=", formatC(p,digits=3,format="g"), ")."))
               )
             )
        )
      )
      
    } else if (tipo == "prop2") {
      tst <- res$test; p <- tst$p.value
      p_a <- round(res$x_a/res$n_a, 3); p_b <- round(res$x_b/res$n_b, 3)
      tagList(
        card(card_header(class="bg-success text-white fw-bold",
                         bsicons::bs_icon("check-circle"), " Resultados"),
             card_body(
               layout_columns(col_widths=c(4,4,4),
                              tarjeta("p-valor", fmt_p(p), semaforo_p(p)),
                              tarjeta(paste0("Prop. ", res$ga), p_a, "bg-info"),
                              tarjeta(paste0("Prop. ", res$gb), p_b, "bg-info")
               ),
               div(class="dm-result-panel rounded mt-3",
                   tags$b("Success defined as: "), paste0(res$var_y, " > mediana global (", round(res$med_global,3), ")"), br(),
                   tags$b(res$ga,": "), paste0(res$x_a,"/",res$n_a," = ",p_a), br(),
                   tags$b(res$gb,": "), paste0(res$x_b,"/",res$n_b," = ",p_b), br(), br(),
                   tags$b("Interpretation: "),
                   if (p < 0.05)
                     tags$span(style="color:#CC0000;font-weight:600;",
                               paste0("Las proporciones son significativamente distintas (p=",
                                      formatC(p,digits=3,format="g"), ")."))
                   else
                     tags$span(class="dm-muted-txt",
                               paste0("No hay diferencia significativa en proporciones (p=",
                                      formatC(p,digits=3,format="g"), ")."))
               )
             )
        )
      )
    }
  })
  
  # Tabla shapiro detalle
  output$tabla_shapiro_detalle <- renderDT({
    res <- test_resultado(); req(!is.null(res), res$tipo=="shapiro")
    df_sw <- do.call(rbind, lapply(res$tests, function(x) {
      data.frame(Grupo=x$grp, n=length(x$vals),
                 W=round(x$res$statistic,4), p_valor=formatC(x$res$p.value,digits=4,format="g"),
                 Normal=if(x$res$p.value>0.05) "✓ Yes" else "✗ No",
                 stringsAsFactors=FALSE)
    }))
    datatable(df_sw, rownames=FALSE, options=list(dom="t", pageLength=20),
              class="table table-sm table-hover") |>
      formatStyle("Normal",
                  color=styleEqual(c("✓ Yes","✗ No"), c("#155724","#721c24")),
                  fontWeight="bold")
  }, server=FALSE)
  
  # ── Plot title ─────────────────────────────────────────────────────────
  output$ui_test_grafico_titulo <- renderUI({
    res <- test_resultado(); if (is.null(res) || !is.null(res$error)) return(NULL)
    tipo <- res$tipo
    lbl  <- switch(tipo,
                   wilcoxon  = "Comparative boxplot with significance annotation",
                   kruskal   = "Boxplot por grupos",
                   chisq     = , fisher = "Mosaico de la tabla de contingencia",
                   spearman  = , pearson = "Scatter plot with trend line",
                   shapiro   = "QQ-plot (cuantil-cuantil) por grupo",
                   ks        = "Cumulative distribution functions (ECDF)",
                   prop_binom = "Proportion bars with exact binomial CI",
                   prop2     = "Comparative proportion bars",
                   burden    = "Boxplot de carga por familia",
                   ""
    )
    tags$small(class="text-muted", lbl)
  })
  
  # ── Test plot ──────────────────────────────────────────────────────────
  output$plot_test_resultado <- renderPlotly({
    res <- test_resultado()
    if (is.null(res) || !is.null(res$error)) {
      return(plot_ly() |> layout(
        annotations=list(list(text="Run the test to see the plot",
                              x=0.5,y=0.5,xref="paper",yref="paper",showarrow=FALSE,
                              font=list(size=16,color="#AAAAAA"))),
        xaxis=list(visible=FALSE), yaxis=list(visible=FALSE),
        paper_bgcolor="#FAFAFA"))
    }
    tipo <- res$tipo
    
    if (tipo %in% c("wilcoxon","burden","ks")) {
      ga <- res$grupos$A$nombre; gb <- res$grupos$B$nombre
      ya <- res$grupos$A$vals;   yb <- res$grupos$B$vals
      p  <- res$test$p.value
      p_lbl <- paste0("p = ", formatC(p, digits=3, format="g"),
                      if(p<0.001)" ***" else if(p<0.01)" **" else if(p<0.05)" *" else " ns")
      
      if (tipo == "ks") {
        # ECDF
        df_ecdf <- rbind(
          data.frame(val=sort(ya), g=ga, stringsAsFactors=FALSE),
          data.frame(val=sort(yb), g=gb, stringsAsFactors=FALSE)
        )
        fig <- plot_ly()
        for (grp in c(ga,gb)) {
          sub <- df_ecdf[df_ecdf$g==grp,]
          n   <- nrow(sub)
          fig <- add_trace(fig, x=sub$val, y=seq_len(n)/n, type="scatter", mode="lines",
                           name=grp, line=list(width=2),
                           hovertemplate=paste0(grp,"<br>x=%{x:.3f}<br>F(x)=%{y:.3f}<extra></extra>"))
        }
        layout(fig, xaxis=list(title=res$var_y %||% "valor"),
               yaxis=list(title="F(x) acumulada", range=c(0,1)),
               legend=list(orientation="h",y=-0.2),
               annotations=list(list(x=0.5,y=1.05,xref="paper",yref="paper",
                                     showarrow=FALSE,text=p_lbl,font=list(size=13,color=if(p<0.05)"#CC0000"else"#555"))),
               plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF",margin=list(t=50))
      } else {
        # Boxplot + puntos
        fig <- plot_ly()
        for (grp in c(ga,gb)) {
          vals <- if(grp==ga) ya else yb
          col  <- if(grp==ga) "#2C6FAC" else "#E67E22"
          fig  <- add_trace(fig, y=vals, name=grp, type="box",
                            boxpoints="all", jitter=0.35, pointpos=0,
                            marker=list(size=4,color=paste0(col,"88"),
                                        line=list(color=col,width=0.5)),
                            line=list(color=col), fillcolor=paste0(col,"33"),
                            hovertemplate=paste0(grp,"<br>%{y:.3f}<extra></extra>"))
        }
        y_max <- max(c(ya,yb), na.rm=TRUE)
        layout(fig,
               yaxis=list(title=res$var_y %||% "valor"),
               xaxis=list(title=""),
               showlegend=TRUE,
               annotations=list(list(
                 x=0.5, y=y_max * 1.05, xref="paper", yref="y",
                 showarrow=FALSE, text=p_lbl,
                 font=list(size=14, color=if(p<0.05)"#CC0000"else"#555"))),
               margin=list(t=60,b=50),
               plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
      }
      
    } else if (tipo == "kruskal") {
      var_y <- res$var_y; grp <- res$grp; grupos <- res$grupos
      fig <- plot_ly()
      for (g in grupos) {
        vals <- res$df2[[var_y]][as.character(res$df2[[grp]])==g & is.finite(res$df2[[var_y]])]
        col  <- HER_COLORS[g] %||% "#888888"
        fig  <- add_trace(fig, y=vals, name=g, type="box",
                          boxpoints="all", jitter=0.35, pointpos=0,
                          marker=list(size=4,color=paste0(col,"88"),line=list(color=col,width=0.5)),
                          line=list(color=col), fillcolor=paste0(col,"33"))
      }
      p <- res$test$p.value
      layout(fig,
             yaxis=list(title=var_y), xaxis=list(title=""),
             showlegend=FALSE,
             title=list(text=paste0("Kruskal-Wallis · p=",formatC(p,digits=3,format="g"),
                                    if(p<0.05)" ***"else" ns"),
                        font=list(size=13,color=if(p<0.05)"#CC0000"else"#555")),
             margin=list(t=60,b=60,l=60),
             plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
      
    } else if (tipo %in% c("chisq","fisher")) {
      tbl <- res$tabla
      # Heatmap de residuos estandarizados
      if (tipo=="chisq" && !is.null(res$test$residuals)) {
        z_mat <- res$test$residuals
      } else {
        exp_mat <- outer(rowSums(tbl), colSums(tbl)) / sum(tbl)
        z_mat   <- (tbl - exp_mat) / sqrt(exp_mat + 1e-9)
      }
      plot_ly(z=z_mat, x=colnames(tbl), y=rownames(tbl), type="heatmap",
              colorscale=list(c(0,"#0066CC"),c(0.5,"#FFFFFF"),c(1,"#CC0000")),
              zmid=0,
              text=matrix(paste0("Obs:",tbl,"\nRes:",round(z_mat,2)), nrow=nrow(tbl)),
              hovertemplate="Fila:%{y}<br>Col:%{x}<br>%{text}<extra></extra>") |>
        layout(xaxis=list(title=res$v2,tickangle=-30),
               yaxis=list(title=res$v1,autorange="reversed"),
               title=list(text="Residuos estandarizados de Pearson (rojo=mayor de lo esperado)",
                          font=list(size=11)),
               margin=list(t=60,b=80,l=100),paper_bgcolor="#FFFFFF")
      
    } else if (tipo %in% c("spearman","pearson")) {
      x <- res$x; y <- res$y
      rho <- round(res$test$estimate,3); p <- res$test$p.value
      # Trend line (regression)
      fit <- lm(y ~ x)
      x_seq <- seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length.out=100)
      y_fit <- coef(fit)[1] + coef(fit)[2] * x_seq
      x_lbl <- if(isTRUE(res$log_x)) paste0("log10(",res$vx,")") else res$vx
      fig <- plot_ly(x=x, y=y, type="scatter", mode="markers", name="Variantes",
                     marker=list(color="#2C6FAC", size=5, opacity=0.6,
                                 line=list(color="#2C6FAC",width=0.4)),
                     hovertemplate=paste0(x_lbl,"=%{x:.3f}<br>",res$vy,"=%{y:.3f}<extra></extra>")) |>
        add_trace(x=x_seq, y=y_fit, type="scatter", mode="lines", name="Tendencia",
                  line=list(color="#CC0000",width=2,dash="dash")) |>
        layout(xaxis=list(title=x_lbl), yaxis=list(title=res$vy),
               legend=list(orientation="h",y=-0.2),
               title=list(text=paste0(if(tipo=="spearman")"rho"else"r","=",rho,
                                      " · p=",formatC(p,digits=3,format="g"),
                                      if(p<0.05)" *"else" ns"),
                          font=list(size=13,color=if(p<0.05)"#2C6FAC"else"#555")),
               margin=list(t=60), plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
      
    } else if (tipo == "shapiro") {
      # QQ-plots por grupo
      fig <- plot_ly()
      colores <- c("#2C6FAC","#CC0000","#27AE60","#E67E22","#8E44AD","#16A085")
      for (i in seq_along(res$tests)) {
        it  <- res$tests[[i]]
        v   <- sort(it$vals)
        n   <- length(v)
        col <- colores[((i-1) %% length(colores)) + 1]
        probs <- (seq_len(n) - 0.375) / (n + 0.25)
        q_teo <- qnorm(probs)
        fig <- add_trace(fig, x=q_teo, y=v, type="scatter", mode="markers", name=it$grp,
                         marker=list(color=col,size=5,opacity=0.7))
      }
      x_range <- c(-3,3)
      fig <- add_trace(fig, x=x_range, y=x_range*sd(res$tests[[1]]$vals)+mean(res$tests[[1]]$vals),
                       type="scatter", mode="lines", name="Theoretical normal",
                       line=list(color="#AAAAAA",dash="dash",width=1.5), showlegend=TRUE)
      layout(fig, xaxis=list(title="Theoretical normal quantiles"),
             yaxis=list(title="Cuantiles observados"),
             legend=list(orientation="h",y=-0.2),
             title=list(text="Desviación de la línea → no normalidad", font=list(size=11)),
             margin=list(t=50,b=80), plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
      
    } else if (tipo == "prop_binom") {
      p_obs <- res$x / res$n; p0 <- res$p0
      ci    <- binom.test(res$x, res$n)$conf.int
      df_bar <- data.frame(
        etiqueta = c("Observada","Esperada (H0)"),
        prop     = c(p_obs, p0),
        color    = c("#2C6FAC","#AAAAAA")
      )
      plot_ly(df_bar, x=~etiqueta, y=~prop, type="bar",
              marker=list(color=~color, line=list(color="#FFFFFF",width=1)),
              text=~paste0(round(prop*100,1),"%"), textposition="outside",
              hovertemplate="%{x}<br>%{y:.4f}<extra></extra>") |>
        add_trace(x=c("Observada","Observada"), y=c(ci[1],ci[2]),
                  type="scatter", mode="lines", line=list(color="#CC0000",width=3),
                  name="IC 95%", hovertemplate="IC 95%: %{y:.3f}<extra></extra>") |>
        layout(yaxis=list(title="Proportion", range=c(0, max(p_obs,p0)*1.3)),
               xaxis=list(title=""), showlegend=TRUE,
               legend=list(orientation="h",y=-0.2),
               margin=list(t=30,b=60), plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
      
    } else if (tipo == "prop2") {
      ga <- res$ga; gb <- res$gb
      p_a <- res$x_a/res$n_a; p_b <- res$x_b/res$n_b
      ci_a <- binom.test(res$x_a,res$n_a)$conf.int
      ci_b <- binom.test(res$x_b,res$n_b)$conf.int
      df_b <- data.frame(g=c(ga,gb), p=c(p_a,p_b),
                         lo=c(ci_a[1],ci_b[1]), hi=c(ci_a[2],ci_b[2]))
      fig <- plot_ly(df_b, x=~g, y=~p, type="bar",
                     marker=list(color=c("#2C6FAC","#E67E22"),
                                 line=list(color="#FFFFFF",width=1)),
                     text=~paste0(round(p*100,1),"%"), textposition="outside",
                     hovertemplate="%{x}<br>%{y:.4f}<extra></extra>") |>
        add_trace(type="scatter", mode="markers+lines",
                  x=df_b$g, y=df_b$lo, name="IC inf",
                  marker=list(color="#CC0000",symbol="line-ew",size=12,
                              line=list(color="#CC0000",width=2)), showlegend=FALSE) |>
        add_trace(type="scatter", mode="markers",
                  x=df_b$g, y=df_b$hi, name="IC sup",
                  marker=list(color="#CC0000",symbol="line-ew",size=12,
                              line=list(color="#CC0000",width=2)), showlegend=FALSE)
      pval_txt <- paste0("p=",formatC(res$test$p.value,digits=3,format="g"),
                         if(res$test$p.value<0.05)" *"else" ns")
      layout(fig, yaxis=list(title="Proportion",range=c(0,max(df_b$hi)*1.3)),
             xaxis=list(title=""), showlegend=FALSE,
             annotations=list(list(x=0.5,y=max(df_b$hi)*1.2,xref="paper",yref="y",
                                   showarrow=FALSE,text=pval_txt,
                                   font=list(size=14,color=if(res$test$p.value<0.05)"#CC0000"else"#555"))),
             margin=list(t=60,b=60), plot_bgcolor="#FAFAFA",paper_bgcolor="#FFFFFF")
    }
  })
  
  # ── Tabla de datos del test ───────────────────────────────────────────────────
  output$tabla_test_datos <- renderDT({
    res <- test_resultado()
    if (is.null(res) || !is.null(res$error)) return(NULL)
    tipo <- res$tipo
    df_show <- tryCatch({
      if (tipo %in% c("wilcoxon","burden")) {
        rbind(
          data.frame(Grupo=res$grupos$A$nombre, Valor=res$grupos$A$vals, stringsAsFactors=FALSE),
          data.frame(Grupo=res$grupos$B$nombre, Valor=res$grupos$B$vals, stringsAsFactors=FALSE)
        )
      } else if (tipo == "kruskal") {
        df2 <- res$df2
        data.frame(Grupo=as.character(df2[[res$grp]]), Valor=df2[[res$var_y]],
                   stringsAsFactors=FALSE)
      } else if (tipo %in% c("chisq","fisher")) {
        as.data.frame.matrix(res$tabla)
      } else if (tipo %in% c("spearman","pearson")) {
        data.frame(X=res$x, Y=res$y, stringsAsFactors=FALSE)
      } else if (tipo == "shapiro") {
        do.call(rbind, lapply(res$tests, function(x)
          data.frame(Grupo=x$grp, W=round(x$res$statistic,4),
                     p=formatC(x$res$p.value,digits=4,format="g"),
                     Normal=x$res$p.value>0.05, stringsAsFactors=FALSE)))
      } else if (tipo == "ks") {
        rbind(
          data.frame(Grupo=res$grupos$A$nombre, Valor=res$grupos$A$vals),
          data.frame(Grupo=res$grupos$B$nombre, Valor=res$grupos$B$vals))
      } else if (tipo == "prop_binom") {
        data.frame(Variable=res$var, Exito=res$exito, n=res$n, Exitos=res$x,
                   Prop_obs=round(res$x/res$n,4), p0=res$p0)
      } else if (tipo == "prop2") {
        data.frame(Grupo=c(res$ga,res$gb), n=c(res$n_a,res$n_b),
                   Exitos=c(res$x_a,res$x_b),
                   Prop=round(c(res$x_a/res$n_a, res$x_b/res$n_b),4))
      }
    }, error=function(e) data.frame(Info="No disponible"))
    datatable(df_show, rownames=FALSE, options=list(dom="ftp",pageLength=10,scrollX=TRUE),
              class="table table-sm table-hover")
  }, server=TRUE)
  
  # =============================================================================
  # HPO TAB — SERVER LOGIC
  # =============================================================================
  
  hpo_familias_disponibles <- reactive({
    res  <- rv$resultados
    fams <- character(0)
    if (!is.null(res$CNVs) && nrow(res$CNVs) > 0 && "ID_Familia" %in% names(res$CNVs))
      fams <- union(fams, as.character(res$CNVs$ID_Familia))
    if (!is.null(res$SVs) && nrow(res$SVs) > 0 && "ID_Familia" %in% names(res$SVs))
      fams <- union(fams, as.character(res$SVs$ID_Familia))
    sort(unique(fams))
  })
  
  observe({
    fams <- hpo_familias_disponibles()
    updateSelectizeInput(session, "hpo_familias_sel",
                         choices  = fams,
                         selected = if (length(fams) > 0) fams[1] else NULL,
                         server   = TRUE)
    updateSelectizeInput(session, "hpo_fam_nota_sel",
                         choices  = fams,
                         selected = if (length(fams) > 0) fams[1] else NULL,
                         server   = TRUE)
  })
  
  hpo_genes_tabla <- reactive({
    fams_sel  <- input$hpo_familias_sel
    if (is.null(fams_sel) || length(fams_sel) == 0) return(data.frame())
    solo_sfari <- isTRUE(input$hpo_solo_sfari)
    res <- rv$resultados
    dfs <- list()
    if (!is.null(res$CNVs) && nrow(res$CNVs) > 0) dfs[["CNVs"]] <- res$CNVs
    if (!is.null(res$SVs)  && nrow(res$SVs)  > 0) dfs[["SVs"]]  <- res$SVs
    if (length(dfs) == 0) return(data.frame())
    filas <- bind_rows(lapply(dfs, function(df) {
      df_f <- df[as.character(df$ID_Familia) %in% fams_sel, , drop = FALSE]
      if ("Annotation_mode" %in% names(df_f))
        df_f <- df_f[!is.na(df_f$Annotation_mode) & df_f$Annotation_mode == "full", , drop = FALSE]
      df_f
    }))
    if (is.null(filas) || nrow(filas) == 0) return(data.frame())
    # Buscar columna de genes de forma robusta
    gene_col <- names(filas)[grep("gene", names(filas), ignore.case = TRUE)][1]
    if (is.na(gene_col)) return(data.frame())
    if (solo_sfari && "En_Panel_Genes" %in% names(filas))
      filas <- filas[es_sfari_col(filas$En_Panel_Genes), , drop = FALSE]
    gene_raw <- toupper(trimws(as.character(filas[[gene_col]])))
    genes_v  <- unique(trimws(unlist(strsplit(
      paste(gene_raw[!is.na(gene_raw) & gene_raw != "" & gene_raw != "."], collapse = ";"),
      "[;,/|[:space:]]+"
    ))))
    genes_v <- sort(genes_v[nzchar(genes_v) & genes_v != "." & toupper(genes_v) != "NA"])
    if (length(genes_v) == 0) return(data.frame())
    data.frame(gen = genes_v, stringsAsFactors = FALSE)
  })
  
  output$ui_hpo_n_genes <- renderUI({
    n <- nrow(hpo_genes_tabla())
    if (n == 0) tags$span(class = "badge bg-secondary", "0")
    else        tags$span(class = "badge bg-primary",   n)
  })
  
  output$ui_hpo_genes_lista <- renderUI({
    df_g <- hpo_genes_tabla()
    if (nrow(df_g) == 0)
      return(div(class = "text-muted small fst-italic p-2",
                 "Sin genes para las familias seleccionadas."))
    busq  <- trimws(input$hpo_buscar_gen %||% "")
    genes <- df_g$gen
    if (nzchar(busq))
      genes <- genes[grepl(busq, genes, ignore.case = TRUE)]
    if (length(genes) == 0)
      return(div(class = "text-muted small fst-italic p-2", "Sin coincidencias."))
    hpo_color <- "#007B7B"
    div(
      class = "d-flex flex-wrap gap-1 p-1",
      style = "max-height:calc(90vh - 230px); overflow-y:auto;",
      lapply(genes, function(g) {
        tags$button(
          class   = "btn btn-sm",
          style   = paste0("font-size:0.78em; font-weight:600; background:", hpo_color,
                           "; color:#fff; border-color:", hpo_color,
                           "; border-radius:14px; padding:2px 10px;"),
          type    = "button",
          title   = paste0("Ver HPO para ", g),
          onclick = paste0("hpoSelGen('", g, "');"),
          g
        )
      })
    )
  })
  
  # ── HPO API: estado reactivo ────────────────────────────────────────────────
  rv$hpo_api_gen    <- NULL
  rv$hpo_api_data   <- NULL
  rv$hpo_api_estado <- "idle"  # idle | cargando | ok | error
  rv$hpo_api_msg    <- ""
  
  # Function: calls HPO API with SSL fallback and detailed error messages
  .hpo_get_safe <- function(url) {
    # Attempt 1: standard SSL verification
    resp <- tryCatch(
      httr::GET(url, httr::timeout(15),
                httr::add_headers(`User-Agent` = "R/Shiny HPO viewer")),
      error = function(e) list(err = conditionMessage(e))
    )
    # If SSL fails, retry without verification
    if (is.list(resp) && !is.null(resp$err)) {
      resp <- tryCatch(
        httr::GET(url, httr::timeout(15),
                  httr::config(ssl_verifypeer = FALSE),
                  httr::add_headers(`User-Agent` = "R/Shiny HPO viewer")),
        error = function(e) list(err = conditionMessage(e))
      )
    }
    resp
  }
  
  # NEW VERSION - Replaces the .hpo_buscar_gen function
  .hpo_buscar_gen <- function(gen) {
    library(httr)
    library(jsonlite)
    
    gen <- toupper(trimws(gen))
    
    # Attempt 1: GraphQL API (more modern)
    resp <- tryCatch({
      url <- "https://hpo.jax.org/api/hpo/graphql"
      query <- sprintf('{genes(query: "%s", limit: 1) {geneId geneName geneSymbol}}', gen)
      
      POST(url,
           body = toJSON(list(query = query)),
           content_type("application/json"),
           timeout(10),
           config(ssl_verifypeer = 0))
    }, error = function(e) {
      message("GraphQL failed, trying alternative...")
      NULL
    })
    
    # Fallback: MyGene.info (muy estable)
    if (is.null(resp) || status_code(resp) != 200) {
      resp <- tryCatch({
        url <- sprintf("https://mygene.info/v3/query?q=symbol:%s&species=human&size=1", gen)
        GET(url, timeout(10))
      }, error = function(e) NULL)
      
      if (is.null(resp) || status_code(resp) != 200) {
        return(list(ok = FALSE, 
                    msg = paste0("No se pudo conectar con HPO ni MyGene.info.\n",
                                 "Check your internet connection."))
      )}
      
      data <- content(resp, "parsed")
      if (length(data$hits) == 0) {
        return(list(ok = FALSE, msg = paste0("Gen «", gen, "» no encontrado")))
      }
      
      hit <- data$hits[[1]]
      return(list(
        ok = TRUE,
        gen = hit$symbol %||% gen,
        gene_id = hit$entrezgene,
        data = hit
      ))
    }
    
    # Si GraphQL funciona
    data <- content(resp, "parsed")
    if (is.null(data$data) || length(data$data$genes) == 0) {
      return(list(ok = FALSE, msg = "Gen no encontrado en HPO"))
    }
    
    g <- data$data$genes[[1]]
    list(ok = TRUE, gen = g$geneSymbol, gene_id = g$geneId, data = g)
  }
  
  observeEvent(input$hpo_gen_seleccionado, {
    gen <- trimws(input$hpo_gen_seleccionado$gen %||% "")
    if (!nzchar(gen)) return()
    rv$hpo_api_gen    <- gen
    rv$hpo_api_estado <- "cargando"
    rv$hpo_api_data   <- NULL
    rv$hpo_api_msg    <- ""
    result <- .hpo_buscar_gen(gen)
    if (isTRUE(result$ok)) {
      rv$hpo_api_gen    <- result$gen
      rv$hpo_api_data   <- result$data
      rv$hpo_api_estado <- "ok"
    } else {
      rv$hpo_api_estado <- "error"
      rv$hpo_api_msg    <- result$msg
    }
  })
  
  output$ui_hpo_visor_badge <- renderUI({
    gen    <- rv$hpo_api_gen
    estado <- rv$hpo_api_estado
    if (is.null(gen) || !nzchar(gen)) return(NULL)
    color <- switch(estado, ok = "#007B7B", cargando = "#6c757d", error = "#dc3545", "#6c757d")
    tags$span(class = "badge",
              style = paste0("background:", color, "; color:#fff; font-size:0.82em;"),
              paste0("\U0001f9ec ", gen))
  })
  
  output$ui_hpo_visor <- renderUI({
    estado    <- rv$hpo_api_estado
    gen       <- rv$hpo_api_gen
    hpo_color <- "#007B7B"
    
    if (is.null(gen) || !nzchar(gen)) {
      return(div(
        class = "d-flex flex-column align-items-center justify-content-center text-center gap-3",
        style = "min-height:500px;",
        tags$span(style = "font-size:3.5em;", "🧬"),
        tags$h5(class = "text-muted fw-normal", "Selecciona un gen"),
        tags$p(class = "text-muted small mb-0",
               "Haz clic en cualquier gen de la columna izquierda"),
        tags$p(class = "text-muted small",
               "to see its HPO terms on the official website.")
      ))
    }
    
    # HPO search URL
    url_hpo <- paste0("https://hpo.jax.org/browse/search?q=",
                      URLencode(gen, reserved = TRUE),
                      "&navFilter=all")
    
    # BLOQUE 1: El Iframe
    visor_iframe <- div(
      div(class = "d-flex align-items-center justify-content-between gap-2 mb-2 p-2 rounded",
          style = paste0("background:", hpo_color, "15; border:1px solid ", hpo_color, "40;"),
          div(
            tags$span(class = "fw-bold", gen),
            tags$small(class = "text-muted d-block", "HPO Search")
          ),
          tags$a(href = url_hpo, target = "_blank", class = "btn btn-sm",
                 style = paste0("background:", hpo_color, "; color:#fff;"),
                 "↗ Open in new tab")
      ),
      tags$iframe(
        src = url_hpo,
        style = "width:100%; height:550px; border:1px solid #ddd; border-radius:6px;",
        frameborder = "0",
        allow = "fullscreen"
      ),
      div(class = "text-muted small mt-2 p-2 rounded",
          style = "background:#f5f5f5;",
          tags$span(class = "fw-bold", "🛠 Info:"), " ",
          "The page loads directly from hpo.jax.org. ",
          "If content is not visible, open in new tab.")
    )
    
    # BLOQUE 2: Renderizado resultado OK (API)
    visor_api <- NULL
    if (estado == "ok" && !is.null(rv$hpo_api_data)) {
      d          <- rv$hpo_api_data
      gene_name  <- d$geneName     %||% gen
      entrez_id  <- d$entrezGeneId %||% ""
      cat_map    <- d$catTermsMap  %||% list()
      enf_list   <- d$dbDiseases   %||% list()
      n_terms    <- sum(vapply(cat_map, length, integer(1)))
      
      enf_ui <- if (length(enf_list) > 0) {
        tagList(
          tags$h6(class = "mt-3 mb-1 fw-bold small text-uppercase",
                  style = "color:#8E44AD; letter-spacing:.05em;",
                  "🧬 Enfermedades asociadas"),
          tags$ul(class = "list-group list-group-flush mb-2",
                  style = "max-height:180px; overflow-y:auto;",
                  lapply(enf_list, function(e) {
                    dis_id   <- e$dbId        %||% ""
                    dis_name <- e$diseaseName %||% ""
                    dis_url  <- if (grepl("^OMIM", dis_id))
                      paste0("https://omim.org/entry/", gsub("\\D", "", dis_id))
                    else
                      paste0("https://hpo.jax.org/browse/disease/", dis_id)
                    tags$li(class = "list-group-item py-1 px-2 small",
                            style = "border-left:3px solid #8E44AD;",
                            tags$a(href = dis_url, target = "_blank",
                                   style = "text-decoration:none; color:#8E44AD; font-weight:600;",
                                   dis_id), " ", dis_name)
                  }))
        )
      } else NULL
      
      cat_uis <- lapply(names(cat_map), function(cat_name) {
        terminos <- cat_map[[cat_name]]
        if (length(terminos) == 0) return(NULL)
        div(class = "mb-3",
            tags$h6(class = "fw-bold mb-1",
                    style = "color:#555; text-transform:uppercase; letter-spacing:.04em; font-size:0.72em;",
                    paste0(cat_name, " (", length(terminos), ")")),
            div(lapply(terminos, function(t) {
              hp_id   <- t$ontologyId %||% ""
              hp_name <- t$name       %||% ""
              tags$a(
                href   = paste0("https://hpo.jax.org/browse/term/", hp_id),
                target = "_blank",
                class  = "badge text-decoration-none me-1 mb-1",
                style  = paste0("background:", hpo_color, "1A; color:", hpo_color,
                                "; border:1px solid ", hpo_color, "55; font-size:0.76em; font-weight:500;"),
                paste0(hp_id, " – ", hp_name)
              )
            }))
        )
      })
      
      visor_api <- div(
        class="mt-4", # iframe separator
        div(class = "p-3 rounded mb-3",
            style = paste0("background:", hpo_color, "10; border:2px solid ", hpo_color, "40;"),
            div(class = "d-flex align-items-start justify-content-between gap-2",
                div(
                  tags$h5(class = "mb-0 fw-bold", style = paste0("color:", hpo_color, ";"), gene_name),
                  if (nzchar(as.character(entrez_id)))
                    tags$small(class = "text-muted", paste0("NCBI Gene ID: ", entrez_id))
                ),
                tags$a(href   = paste0("https://hpo.jax.org/browse/search?q=",
                                       URLencode(gene_name, TRUE), "&navFilter=all"),
                       target = "_blank",
                       class  = "btn btn-sm flex-shrink-0",
                       style  = paste0("background:", hpo_color, "; color:#fff; border-color:", hpo_color, ";"),
                       "↗ Ver en HPO")
            ),
            tags$small(class = "text-muted d-block mt-1",
                       paste0(n_terms, " terms in ", length(cat_map),
                              " categories • ", length(enf_list), " disease(s)"))
        ),
        enf_ui,
        if (length(cat_uis) > 0)
          tagList(
            tags$h6(class = "mt-1 mb-2 fw-bold small text-uppercase",
                    style = paste0("color:", hpo_color, "; letter-spacing:.05em;"),
                    "📋 Términos HPO por categoría"),
            cat_uis
          )
      )
    }
    
    # Devolver ambos elementos en la interfaz
    return(tagList(visor_iframe, visor_api))
    
  }) # ✔️ CIERRE CORRECTO DEL renderUI
  
  output$ui_hpo_nota_panel <- renderUI({
    fam <- input$hpo_fam_nota_sel
    if (is.null(fam) || !nzchar(fam))
      return(div(class = "text-muted small fst-italic p-2",
                 "Selecciona una familia para a\u00f1adir su descripci\u00f3n HPO."))
    nota_actual <- rv$hpo_notas[[fam]] %||% ""
    tiene_nota  <- nzchar(trimws(nota_actual))
    # Genes de referencia de esta familia
    res  <- rv$resultados
    dfs  <- list()
    if (!is.null(res$CNVs)) dfs[["CNVs"]] <- res$CNVs
    if (!is.null(res$SVs))  dfs[["SVs"]]  <- res$SVs
    genes_fam <- tryCatch({
      filas <- bind_rows(lapply(dfs, function(df) {
        df_f <- df[as.character(df$ID_Familia) == fam, , drop = FALSE]
        if ("Annotation_mode" %in% names(df_f))
          df_f <- df_f[!is.na(df_f$Annotation_mode) & df_f$Annotation_mode == "full", , drop = FALSE]
        df_f
      }))
      if (is.null(filas) || nrow(filas) == 0) return(character(0))
      gc <- names(filas)[grep("gene", names(filas), ignore.case = TRUE)][1]
      if (is.na(gc)) return(character(0))
      gr <- toupper(trimws(as.character(filas[[gc]])))
      gv <- unique(trimws(unlist(strsplit(paste(gr[!is.na(gr) & gr != "" & gr != "."], collapse=";"), "[;,/|[:space:]]+"))))
      sort(gv[nzchar(gv) & gv != "." & toupper(gv) != "NA"])
    }, error = function(e) character(0))
    
    div(class = "d-flex flex-column gap-2",
        if (length(genes_fam) > 0)
          div(class = "p-2 rounded",
              style = "background:var(--dm-info-bg); border:1px solid var(--dm-info-border); font-size:0.82em;",
              tags$b(paste0("Genes de ", fam, ":")), tags$br(),
              tags$span(class = "text-muted",
                        if (length(genes_fam) <= 20) paste(genes_fam, collapse = " \u00b7 ")
                        else paste0(paste(genes_fam[1:20], collapse = " \u00b7 "),
                                    " \u2026 (+", length(genes_fam)-20, ")")))
        else NULL,
        div(
          tags$label(class = "form-label fw-bold small",
                     paste0("\U0001f4cb Descripci\u00f3n HPO \u2014 ", fam)),
          tags$textarea(id = "hpo_input_nota", class = "form-control form-control-sm",
                        style = "min-height:200px; font-size:0.85em; resize:vertical;",
                        placeholder = paste0("Escribe aqu\u00ed los HPOs observados en ", fam,
                                             "...\nEjemplo:\nHP:0001250 - Seizures\nHP:0001999 - Abnormal facial shape"),
                        nota_actual)
        ),
        div(class = "d-flex gap-2 flex-wrap",
            actionButton("hpo_btn_guardar_nota", "\U0001f4be Guardar descripci\u00f3n",
                         class = "btn-primary btn-sm"),
            if (tiene_nota)
              actionButton("hpo_btn_borrar_nota", "\U0001f5d1 Borrar",
                           class = "btn-outline-danger btn-sm")
        ),
        if (tiene_nota) {
          div(class = "p-2 rounded small",
              style = "background:var(--dm-success-bg); border-left:3px solid var(--dm-success-left); white-space:pre-wrap; font-size:0.8em;",
              tags$b("Guardado: "),
              if (nchar(nota_actual) > 200) paste0(substr(nota_actual,1,200),"\u2026") else nota_actual
          )
        } else {
          NULL
        }
    )
  })
  
  observeEvent(input$hpo_btn_guardar_nota, {
    fam   <- input$hpo_fam_nota_sel
    texto <- trimws(input$hpo_input_nota %||% "")
    if (is.null(fam) || !nzchar(fam)) return()
    rv$hpo_notas[[fam]] <- texto
    showNotification(paste0("\U0001f4be Descripci\u00f3n HPO guardada para ", fam),
                     type = "message", duration = 3)
  })
  
  observeEvent(input$hpo_btn_borrar_nota, {
    fam <- input$hpo_fam_nota_sel
    if (is.null(fam) || !nzchar(fam)) return()
    rv$hpo_notas[[fam]] <- NULL
    updateTextAreaInput(session, "hpo_input_nota", value = "")
    showNotification(paste0("\U0001f5d1 Descripci\u00f3n HPO eliminada para ", fam),
                     type = "warning", duration = 3)
  })
}


shinyApp(ui, server)