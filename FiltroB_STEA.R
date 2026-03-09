library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(openxlsx)
library(parallel)          # procesamiento paralelo de familias
library(GenomicRanges)     # overlap preciso: width, prop_referencia, Jaccard

# =============================================================================
# 1. CONFIGURACIÓN GENERAL
# =============================================================================
CONFIG <- list(
  ruta_entrada    = "",                                  # Set to your input folder (with HI/PA/MA subfolders per family),
  ruta_salida     = file.path(getwd(), "results"),           # Output folder (created automatically),
  ruta_rangos     = file.path(getwd(), "reference", "reference_ranges.xlsx"),  # Reference ranges file,
  ruta_panel      = file.path(getwd(), "reference", "gene_panel.txt"),           # Gene panel file (.txt, comma-separated),
  umbral_strict        = 0.70,
  umbral_herencia      = 0.70,   # sim mínimo para herencia "Fuerte"   [sim = ovl / max(len_q, len_r)]
  umbral_herencia_lax  = 0.60,   # sim mínimo para herencia "Probable" (rescatada por Jaccard)
  umbral_jaccard       = 0.50,   # Jaccard mínimo requerido en el nivel "Probable"
  margen_lateral       = 3e6,    # 3 Mb de desviación máxima permitida por cada lado
  umbral_longitud_similar = 0.70,  # ratio mínimo min(len_h,len_p)/max(len_h,len_p) para match por gen sin solapamiento
  deduplicar_regiones  = TRUE,   # TRUE = colapsar variantes similares (mismo chr+tipo, sim >= umbral_herencia)
  # FALSE = conservar todas las variantes sin colapso
  n_cores         = max(1L, detectCores(logical = FALSE) - 1L)  # núcleos físicos - 1
)
if(!dir.exists(CONFIG$ruta_salida)) dir.create(CONFIG$ruta_salida, recursive = TRUE)

# =============================================================================
# 2. CARGA DE REFERENCIAS (CON VALIDACIÓN)
# =============================================================================
cat(">>> Cargando bases de datos de referencia...\n")

# Cargar Rangosr
if(!file.exists(CONFIG$ruta_rangos)) stop("❌ No se encuentra: ", CONFIG$ruta_rangos)

df_ref_raw <- read_excel(CONFIG$ruta_rangos)
dt_ref <- data.table(
  ref_name           = as.character(df_ref_raw[[1]]),
  ref_type           = toupper(str_trim(df_ref_raw[[2]])),
  ref_chr            = toupper(str_trim(gsub("CHR", "", df_ref_raw[[3]]))),
  ref_start_amplio   = as.numeric(df_ref_raw[[5]]),
  ref_start_estricto = as.numeric(df_ref_raw[[6]]),
  ref_end_estricto   = as.numeric(df_ref_raw[[7]]),
  ref_end_amplio     = as.numeric(df_ref_raw[[8]])
)
dt_ref <- dt_ref[!is.na(ref_start_amplio) & !is.na(ref_end_amplio) & 
                   !is.na(ref_name) & ref_name != "" & ref_name != "Nombre_CNV"]
cat("✓ Rangos cargados:", nrow(dt_ref), "\n")

# Cargar Genes SFARI
if(!file.exists(CONFIG$ruta_panel)) stop("❌ No se encuentra: ", CONFIG$ruta_panel)

lineas_panel <- readLines(CONFIG$ruta_panel, warn = FALSE)
genes_panel <- unique(toupper(trimws(unlist(strsplit(paste(lineas_panel, collapse = ","), ",")))))
genes_panel <- genes_panel[genes_panel != ""]
# Hash set: convierte el vector en environment para lookup O(1) en lugar de O(n)
genes_panel_set <- new.env(hash = TRUE, parent = emptyenv())
for(.g in genes_panel) assign(.g, TRUE, envir = genes_panel_set)
rm(.g)
cat("✓ Genes del panel cargados:", length(genes_panel), "\n")

# =============================================================================
# 3. FUNCIÓN DE DEDUPLICACIÓN POR REGIONES SIMILARES
# =============================================================================
# Lógica: dentro de cada resultado de familia, si dos o más variantes solapan
# con sim_maxlen >= umbral_herencia en el mismo cromosoma y mismo tipo de SV/CNV,
# se consideran "la misma región". De ese grupo se conservan:
#   - Las filas con Tipo_Herencia distinta de "De novo", si las hay.
#   - Si todas son "De novo" o NA, se conserva la de mayor ranking_numeric.
# Las filas "split" siguen siempre a su fila "full" correspondiente.
#
# Casos especiales de herencia:
#   - Paterna + Materna similares → se fusionan en Conjunta (no se descarta ninguna)
#   - Herencias distintas no fusionables (p.ej. Paterna vs De novo) → NO se colapsan,
#     son eventos biológicamente distintos aunque sus coordenadas solapen
#   - Igual herencia o ambas De novo/NA → colapso normal, se conserva la de mayor ranking
#
# Algoritmo greedy por pares directos (no Union-Find):
#   Solo se elimina una variante si tiene sim DIRECTO >= umbral con su representante.
#   Evita colapsos transitivos: A~B y B~C no implica colapsar A con C si sim(A,C)<umbral.
# ---------------------------------------------------------------------------

# Función auxiliar: extraer la familia base de herencia (sin genotipo ni "(Probable)")
.her_base <- function(x) {
  if (is.na(x)) return(NA_character_)
  # Quitar el sufijo " [GT]" o " [GT | GT]" y "(Probable)"
  x <- trimws(sub("\\s*\\[.*", "", x))
  x <- trimws(sub("\\s*\\(Probable\\)", "", x))
  x
}

deduplicar_regiones_similares <- function(dt, umbral_sim = CONFIG$umbral_herencia) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  
  tiene_modo  <- "Annotation_mode" %in% colnames(dt)
  tiene_annot <- "AnnotSV_ID"      %in% colnames(dt)
  
  # Trabajar solo sobre filas "full"; las "split" se reconstruyen al final por AnnotSV_ID
  dt_full <- if (tiene_modo) dt[Annotation_mode == "full"] else copy(dt)
  if (nrow(dt_full) <= 1) return(dt)
  
  # --- Construir metadatos por fila ---
  n         <- nrow(dt_full)
  fn_row    <- seq_len(n)
  sv_up     <- toupper(as.character(dt_full$SV_type))
  tipo_vec  <- fcase(
    grepl("DEL|LOSS", sv_up), "DEL",
    grepl("DUP|GAIN", sv_up), "DUP",
    default = "OTRO"
  )
  start_vec <- as.integer(dt_full$SV_start)
  end_vec   <- as.integer(dt_full$SV_end)
  len_vec   <- end_vec - start_vec + 1L
  chr_vec   <- toupper(trimws(gsub("CHR", "", as.character(dt_full$SV_chrom))))
  
  her_vec   <- as.character(dt_full$Tipo_Herencia)
  her_base  <- vapply(her_vec, .her_base, character(1L))  # familia base sin decoración
  is_denovo <- !is.na(her_vec) & her_base %in% c("De novo", "Desconocida")
  rank_num  <- suppressWarnings(as.numeric(dt_full$ranking_numeric))
  rank_num[is.na(rank_num)] <- -Inf
  tiebreak  <- if (tiene_annot) as.character(dt_full$AnnotSV_ID) else as.character(fn_row)
  
  # --- Calcular todos los pares con solapamiento ---
  gr <- GRanges(
    seqnames = chr_vec,
    ranges   = IRanges(start = start_vec, end = end_vec)
  )
  
  hits <- findOverlaps(gr, gr, type = "any")
  hits <- hits[queryHits(hits) < subjectHits(hits)]
  if (length(hits) == 0) return(dt)
  
  hi <- queryHits(hits)
  hj <- subjectHits(hits)
  
  # Filtrar: mismo tipo conocido
  tipo_ok <- tipo_vec[hi] == tipo_vec[hj] & tipo_vec[hi] != "OTRO"
  hi <- hi[tipo_ok];  hj <- hj[tipo_ok]
  if (length(hi) == 0) return(dt)
  
  # Calcular sim
  ovl_w <- pmax(0L, pmin(end_vec[hj], end_vec[hi]) - pmax(start_vec[hj], start_vec[hi]) + 1L)
  sim_v <- ovl_w / pmax(len_vec[hi], len_vec[hj])
  
  # Filtrar por umbral
  ok    <- sim_v >= umbral_sim
  hi    <- hi[ok];  hj <- hj[ok];  sim_v <- sim_v[ok]
  if (length(hi) == 0) return(dt)
  
  # Ordenar: sim desc, desempate determinista
  ord   <- order(-sim_v, tiebreak[hi], tiebreak[hj])
  hi    <- hi[ord];  hj <- hj[ord];  sim_v <- sim_v[ord]
  
  # --- Greedy: procesar cada par ---
  eliminado <- logical(n)
  # fusionar[i] = índice j cuya herencia debe fusionarse en i (caso Paterna+Materna→Conjunta)
  fusionar  <- integer(n)
  
  for (k in seq_along(hi)) {
    i <- hi[k];  j <- hj[k]
    if (eliminado[i] || eliminado[j]) next
    
    hi_base <- her_base[i]
    hj_base <- her_base[j]
    
    # ---- Caso 1: Paterna + Materna (o viceversa) → fusionar en Conjunta ----
    # Son la misma variante vista desde dos progenitores distintos.
    # Se conserva i, se marca j para eliminar, y se anota que i debe heredar
    # también la información de j para construir la etiqueta Conjunta.
    es_fusion_conjunta <- (
      (!is.na(hi_base) & !is.na(hj_base)) &&
        ((hi_base == "Paterna" & hj_base == "Materna") ||
           (hi_base == "Materna" & hj_base == "Paterna"))
    )
    if (es_fusion_conjunta) {
      # Asegurar que i sea Paterna y j Materna para unificar la lógica de fusión
      if (hi_base == "Materna") { tmp <- i; i <- j; j <- tmp }
      eliminado[j]  <- TRUE
      fusionar[i]   <- j    # i (Paterna) absorbe j (Materna) → Conjunta
      next
    }
    
    # ---- Caso 2: herencias distintas NO fusionables → NO colapsar ----
    # Ejemplos: Paterna vs De novo, Materna vs Conjunta, etc.
    # Son eventos biológicamente distintos; se conservan ambas.
    her_i_conocida <- !is.na(hi_base) & !hi_base %in% c("De novo", "Desconocida")
    her_j_conocida <- !is.na(hj_base) & !hj_base %in% c("De novo", "Desconocida")
    if (her_i_conocida && her_j_conocida && hi_base != hj_base) {
      cat(paste0("      [Dedup] Par no colapsado: herencias distintas (",
                 hi_base, " vs ", hj_base, ") con sim=", round(sim_v[k], 4), "\n"))
      next
    }
    
    # ---- Caso 3: colapso normal ----
    # Misma herencia, ambas De novo/NA, o una conocida y otra sin información.
    # Prioridad: (1) no-denovo > denovo, (2) mayor rank, (3) tiebreak determinista
    i_gana <- if (!is_denovo[i] && is_denovo[j]) {
      TRUE
    } else if (is_denovo[i] && !is_denovo[j]) {
      FALSE
    } else if (rank_num[i] != rank_num[j]) {
      rank_num[i] > rank_num[j]
    } else {
      tiebreak[i] <= tiebreak[j]
    }
    
    if (i_gana) eliminado[j] <- TRUE else eliminado[i] <- TRUE
  }
  
  n_eliminadas <- sum(eliminado)
  if (n_eliminadas == 0L) return(dt)
  
  cat(paste0("      [Dedup] Colapsadas ", n_eliminadas,
             " variante(s) (mismo chr+tipo, sim >= ", umbral_sim, ")\n"))
  
  # --- Aplicar fusiones Conjunta ---
  # Para cada i con fusionar[i] > 0: i es Paterna, j=fusionar[i] es Materna.
  # Construir la etiqueta Conjunta combinando genotipos de ambas.
  for (i in which(fusionar > 0L)) {
    j       <- fusionar[i]
    her_i   <- her_vec[i]   # "Paterna [GT]" o "Paterna (Probable) [GT]"
    her_j   <- her_vec[j]   # "Materna [GT]" o "Materna (Probable) [GT]"
    
    # Extraer genotipo de cada etiqueta: lo que hay entre [ ]
    gt_pa  <- regmatches(her_i, regexpr("(?<=\\[)[^\\]]+(?=\\])", her_i, perl = TRUE))
    gt_ma  <- regmatches(her_j, regexpr("(?<=\\[)[^\\]]+(?=\\])", her_j, perl = TRUE))
    gt_pa  <- if (length(gt_pa) > 0) gt_pa else "?"
    gt_ma  <- if (length(gt_ma) > 0) gt_ma else "?"
    
    # Nivel de confianza: si alguna de las dos es Probable → Conjunta (Probable)
    es_probable <- grepl("Probable", her_i, fixed = TRUE) ||
      grepl("Probable", her_j, fixed = TRUE)
    etiqueta <- if (es_probable) {
      paste0("Conjunta (Probable) [", gt_pa, " | ", gt_ma, "]")
    } else {
      paste0("Conjunta [", gt_pa, " | ", gt_ma, "]")
    }
    
    dt_full[fn_row == i, Tipo_Herencia := etiqueta]
    
    # Fusionar también IDs de coincidencia si existen las columnas
    if ("IDs_Coincidencia_PA" %in% colnames(dt_full) &&
        "IDs_Coincidencia_MA" %in% colnames(dt_full)) {
      # i era Paterna: tiene IDs_PA pero no IDs_MA → tomar IDs_MA de j (que era Materna)
      ids_ma_j <- dt_full$IDs_Coincidencia_MA[fn_row == j]
      ids_pa_i <- dt_full$IDs_Coincidencia_PA[fn_row == i]
      # Si j (Materna) tenía IDs_PA porque también era Conjunta parcial, conservar ambos
      ids_pa_j <- dt_full$IDs_Coincidencia_PA[fn_row == j]
      pa_final <- paste(na.omit(unique(c(ids_pa_i, ids_pa_j))), collapse = ";")
      pa_final <- if (nchar(pa_final) == 0) NA_character_ else pa_final
      ma_final <- if (!is.na(ids_ma_j) && nchar(ids_ma_j) > 0) ids_ma_j else NA_character_
      dt_full[fn_row == i, `:=`(IDs_Coincidencia_PA = pa_final,
                                IDs_Coincidencia_MA = ma_final)]
    }
  }
  
  # --- Reconstruir resultado ---
  filas_conservar <- fn_row[!eliminado]
  
  if (tiene_modo && tiene_annot) {
    ids_conservar <- dt_full$AnnotSV_ID[filas_conservar]
    # Propagar Tipo_Herencia fusionado a las filas "split" correspondientes
    dt_result <- dt[AnnotSV_ID %in% ids_conservar]
    # Actualizar herencia en dt_result para los IDs que fueron fusionados
    ids_fusionados <- dt_full$AnnotSV_ID[fn_row %in% which(fusionar > 0L)]
    if (length(ids_fusionados) > 0) {
      her_fusion_map <- dt_full[fn_row %in% which(fusionar > 0L),
                                .(AnnotSV_ID, Tipo_Herencia,
                                  IDs_Coincidencia_PA, IDs_Coincidencia_MA)]
      cols_her <- intersect(c("Tipo_Herencia", "IDs_Coincidencia_PA", "IDs_Coincidencia_MA"),
                            colnames(dt_result))
      for (col in cols_her) {
        idx <- match(dt_result$AnnotSV_ID, her_fusion_map$AnnotSV_ID)
        vals <- her_fusion_map[[col]][idx]
        dt_result[!is.na(idx), (col) := vals[!is.na(idx)]]
      }
    }
    return(dt_result)
  } else {
    dt_full[, .fn_row_tmp := fn_row]
    dt_result <- dt_full[.fn_row_tmp %in% filas_conservar]
    dt_result[, .fn_row_tmp := NULL]
    return(dt_result)
  }
}

# =============================================================================
# 4. FUNCIÓN OPTIMIZADA DE PROCESAMIENTO
# =============================================================================
procesar_modalidad <- function(ruta_hi, ruta_pa, ruta_ma, MODO, dt_ref, genes_panel) {
  
  # --- VALIDACIÓN DE ENTRADA ---
  if(is.na(ruta_hi) || !file.exists(ruta_hi)) return(NULL)
  
  cat(paste0("   -> Analizando ", MODO, "...\n"))
  
  # --- CARGA OPTIMIZADA ---
  # fill=TRUE: AnnotSV puede generar campos con contenido complejo (comas, puntos y coma,
  # saltos de línea en fenotipos OMIM, etc.) que desplazan columnas en filas concretas.
  # fill=TRUE tolera filas con distinto número de campos en lugar de abortar o desplazar.
  dt_annotsv <- tryCatch({
    fread(ruta_hi, header = TRUE, quote = "", fill = TRUE, showProgress = FALSE)
  }, error = function(e) {
    cat(paste("      [ERROR leyendo archivo]:", e$message, "\n"))
    return(NULL)
  })
  
  if(is.null(dt_annotsv) || nrow(dt_annotsv) == 0) return(NULL)
  
  # Snapshot completo: índice de fila + todas las columnas originales, ANTES de cualquier filtro.
  # Al final del pipeline se usa para recuperar los valores reales de cualquier columna
  # que se haya perdido en los pasos intermedios.
  dt_annotsv[, .src_row_ := .I]
  dt_annotsv_src  <- copy(dt_annotsv)        # todas las columnas, todos los valores
  cols_originales <- setdiff(names(dt_annotsv), ".src_row_")
  cat(paste0("      Columnas leídas: ", length(cols_originales), "\n"))
  
  # Cargar progenitores de forma segura
  dt_padre <- if(!is.na(ruta_pa) && file.exists(ruta_pa)) {
    tryCatch(fread(ruta_pa, header=TRUE, quote="", fill=TRUE, showProgress=FALSE), 
             error = function(e) data.table())
  } else data.table()
  
  dt_madre <- if(!is.na(ruta_ma) && file.exists(ruta_ma)) {
    tryCatch(fread(ruta_ma, header=TRUE, quote="", fill=TRUE, showProgress=FALSE), 
             error = function(e) data.table())
  } else data.table()
  
  dt_annotsv[, ranking_numeric := suppressWarnings(as.numeric(as.character(AnnotSV_ranking_score)))]
  
  # --- FILTRO DE FRECUENCIA POBLACIONAL (PRIMER FILTRO) ---
  # Condición 1: frecuencia en bases de datos de CNVs benignos según tipo de SV
  #   DEL → B_loss_AFmax | DUP → B_gain_AFmax | INS → B_ins_AFmax | INV → B_inv_AFmax
  #   Pasa si la celda está vacía/NA o el valor ≤ 0.01
  # Condición 2: frecuencia en cohorte Illumina DRAGEN (similar counts)
  #   Valor de Illumina_DRAGEN.similar.counts dividido entre 1061
  #   Pasa si la celda está vacía/NA o el cociente ≤ 0.01
  # Condición 3: conteo exacto en cohorte Illumina DRAGEN (exact counts)
  #   Valor directo de Illumina_DRAGEN.exact.counts (entero, sin parseo ni normalización)
  #   Pasa si la celda está vacía/NA o el valor ≤ 10
  # Se conservan las filas que cumplen LAS TRES condiciones (celdas vacías/NA pasan)
  {
    # Tipo de SV por fila (vectorizado)
    sv_type_upper <- toupper(as.character(dt_annotsv$SV_type))
    tipo_sv <- fcase(
      grepl("DEL", sv_type_upper), "DEL",
      grepl("DUP", sv_type_upper), "DUP",
      grepl("INS", sv_type_upper), "INS",
      grepl("INV", sv_type_upper), "INV",
      default = "OTRO"
    )
    
    # Mapa tipo → columna AF
    col_af_map <- c(DEL = "B_loss_AFmax", DUP = "B_gain_AFmax",
                    INS = "B_ins_AFmax",  INV = "B_inv_AFmax")
    
    # Extraer valor AF según tipo de cada fila
    extraer_af <- function(dt, tipo_vec, col_map) {
      af_vals <- rep(NA_real_, nrow(dt))
      for(tipo in names(col_map)) {
        col <- col_map[tipo]
        rows <- which(tipo_vec == tipo)
        if(length(rows) > 0 && col %in% colnames(dt)) {
          af_vals[rows] <- suppressWarnings(as.numeric(as.character(dt[[col]][rows])))
        }
      }
      af_vals
    }
    
    af_vals     <- extraer_af(dt_annotsv, tipo_sv, col_af_map)
    mask_af     <- is.na(af_vals) | af_vals <= 0.01
    
    # Condición 2: Illumina DRAGEN similar counts / 1061
    # La columna tiene formato "DEL=2343;DUP=9485;INS=748;INV=7374;TRA=738"
    # Se extrae el valor del tipo de SV de cada fila, se suma con sus posibles
    # sub-tipos (p.ej. todos los que contienen "DEL") y se divide entre 1061.
    # Pasa si la celda está vacía/NA o el cociente ≤ 0.01
    if("Illumina_DRAGEN.similar.counts" %in% colnames(dt_annotsv)) {
      dragen_str <- as.character(dt_annotsv$Illumina_DRAGEN.similar.counts)
      
      # Para cada fila: parsear el string key=value y sumar los valores cuya
      # clave contenga el tipo de SV de esa fila (p.ej. tipo "DEL" coincide
      # con claves "DEL", "DEL_ALU", etc.)
      dragen_counts <- mapply(function(sv_str, sv_tipo) {
        # Celda vacía o tipo desconocido → NA (pasa el filtro)
        if(is.na(sv_str) || sv_str == "" || sv_str == "." || sv_tipo == "OTRO")
          return(NA_real_)
        # Separar pares key=value
        pares <- unlist(strsplit(sv_str, ";", fixed = TRUE))
        # Extraer clave y valor de cada par
        claves  <- sub("=.*$", "", pares)
        valores <- suppressWarnings(as.numeric(sub("^[^=]+=", "", pares)))
        # Sumar solo los valores cuya clave contenga el tipo de SV de la fila
        idx_match <- grepl(sv_tipo, claves, fixed = TRUE)
        if(!any(idx_match)) return(NA_real_)
        sum(valores[idx_match], na.rm = TRUE)
      }, dragen_str, tipo_sv, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      
      dragen_freq <- as.numeric(dragen_counts) / 1061
      mask_dragen <- is.na(dragen_freq) | dragen_freq <= 0.01
    } else {
      mask_dragen <- rep(TRUE, nrow(dt_annotsv))   # columna ausente → pasa
    }
    
    # Condición 3: Illumina DRAGEN exact counts
    # Conteo directo (entero) — no requiere parseo ni normalización.
    # Pasa si la celda está vacía/NA o el valor ≤ 10.
    if("Illumina_DRAGEN.exact.counts" %in% colnames(dt_annotsv)) {
      exact_counts <- suppressWarnings(
        as.numeric(as.character(dt_annotsv$Illumina_DRAGEN.exact.counts))
      )
      mask_exact <- is.na(exact_counts) | exact_counts <= 10
    } else {
      mask_exact <- rep(TRUE, nrow(dt_annotsv))   # columna ausente → pasa
    }
    
    mask_freq   <- mask_af & mask_dragen & mask_exact
    n_antes     <- nrow(dt_annotsv)
    dt_annotsv  <- dt_annotsv[mask_freq]
    cat(paste0("      Filtro frecuencia: ", n_antes, " → ", nrow(dt_annotsv), " variantes\n"))
  }
  
  if(nrow(dt_annotsv) == 0) {
    cat("      (Sin variantes tras filtro de frecuencia)\n")
    return(NULL)
  }
  
  # --- FILTRO DE TAMAÑO MÍNIMO (solo SVs) ---
  # Se aplica únicamente al modo SV. Elimina variantes cuyo tamaño absoluto
  # sea menor de 50 pb (DEL, DUP, INS, INV).
  # El tamaño se calcula como abs(SV_end - SV_start). Se usa valor absoluto
  # porque en algunas notaciones SV_end puede ser menor que SV_start (p.ej. INS).
  # Filas con coordenadas inválidas (NA) se conservan: el filtro de calidad
  # posterior las descartará sin generar pérdidas silenciosas.
  if(MODO == "SV") {
    sv_size  <- abs(as.numeric(dt_annotsv$SV_end) - as.numeric(dt_annotsv$SV_start))
    mask_tam <- is.na(sv_size) | sv_size >= 50
    n_antes_tam     <- nrow(dt_annotsv)
    dt_annotsv      <- dt_annotsv[mask_tam]
    n_filtradas_tam <- n_antes_tam - nrow(dt_annotsv)
    cat(paste0("      Filtro tamaño (<50 pb): ", n_antes_tam, " → ",
               nrow(dt_annotsv), " variantes",
               if(n_filtradas_tam > 0) paste0(" (eliminadas: ", n_filtradas_tam, ")") else "",
               "\n"))
    if(nrow(dt_annotsv) == 0) {
      cat("      (Sin variantes tras filtro de tamaño)\n")
      return(NULL)
    }
  }
  # --- MARCADO DE GENES SFARI (TOTALMENTE VECTORIZADO) ---
  {
    gene_strings <- toupper(trimws(as.character(dt_annotsv$Gene_name)))
    gene_strings[is.na(gene_strings) | gene_strings == "" | gene_strings == "."] <- ""
    gene_list <- strsplit(gene_strings, "[,;[:space:]]+", fixed = FALSE)
    dt_annotsv[, En_Panel_Genes := vapply(gene_list, function(gs) {
      gs <- gs[nzchar(gs) & gs != "."]
      length(gs) > 0L && any(vapply(gs, exists, logical(1L),
                                    envir = genes_panel_set, inherits = FALSE))
    }, logical(1L))]
  }
  
  # --- FILTRO DE CALIDAD OPTIMIZADO ---
  if(!"Samples_ID" %in% colnames(dt_annotsv)) {
    cat("      ⚠ Columna Samples_ID no encontrada\n")
    return(NULL)
  }
  
  extraer_mask_calidad <- function(dt, modo) {
    if(nrow(dt) == 0 || !"Samples_ID" %in% colnames(dt)) return(logical(0))
    sample_ids    <- as.character(dt$Samples_ID)
    raw_gt_values <- rep(NA_character_, nrow(dt))
    for(sid in unique(sample_ids)) {
      rows <- which(sample_ids == sid)
      if(sid %in% colnames(dt)) raw_gt_values[rows] <- as.character(dt[[sid]][rows])
    }
    # Filtro QUALITY: pasa si la columna no existe, el valor es NA, o el valor > 30.
    # Se calcula una sola vez y se aplica en ambas ramas (CNV y SV).
    quality_col <- if("QUALITY" %in% colnames(dt)) {
      suppressWarnings(as.numeric(as.character(dt$QUALITY)))
    } else {
      rep(NA_real_, nrow(dt))
    }
    mask_quality_score <- is.na(quality_col) | quality_col > 30
    if(modo == "CNV") {
      genotype_part <- sub(":.*", "", raw_gt_values)
      filter_col    <- if("FILTER" %in% colnames(dt)) as.character(dt$FILTER) else rep("PASS", nrow(dt))
      return((genotype_part %in% c("0/1", "1/1", "./1")) &
               (filter_col %in% c("PASS", "High Quality")) &
               mask_quality_score)
    } else {
      parts      <- strsplit(raw_gt_values, ":", fixed = TRUE)
      gt_part    <- vapply(parts, `[`, character(1L), 1L)
      ft_part    <- vapply(parts, function(x) if(length(x) >= 2L) x[2L] else NA_character_, character(1L))
      filter_col <- if("FILTER" %in% colnames(dt)) as.character(dt$FILTER) else rep("PASS", nrow(dt))
      annot_col  <- if("Annotation_mode" %in% colnames(dt)) as.character(dt$Annotation_mode) else rep("full", nrow(dt))
      # Nota: se filtra por Annotation_mode == "full" aquí para SVs porque las filas
      # "split" no tienen formato GT:FT en su columna de muestra y no deben evaluarse
      # por calidad individualmente — su inclusión se hereda de la fila "full" asociada.
      return((gt_part %in% c("0/1", "1/1", "./1")) &
               (ft_part %in% c("PASS", "High Quality")) &   # consistente con modo CNV
               (filter_col %in% c("PASS", "High Quality")) &
               (annot_col  == "full") &
               mask_quality_score)
    }
  }
  
  mask_quality <- extraer_mask_calidad(dt_annotsv, MODO)
  dt_fase1 <- dt_annotsv[mask_quality]
  
  if(nrow(dt_fase1) == 0) {
    cat("      (Sin variantes tras filtro calidad)\n")
    return(NULL)
  }
  
  # Filtro de progenitores: se trabaja únicamente con filas "full" para no contar
  # la misma SV múltiples veces. Adicionalmente, se aplica filtro de genotipo
  # (heterocigoto / homocigoto) para excluir variantes de baja calidad o sin
  # llamada válida, lo que reduce falsos positivos de herencia entre batches.
  filtrar_progenitor <- function(dt) {
    if(nrow(dt) == 0) return(data.table())
    dt_f <- if("Annotation_mode" %in% colnames(dt)) dt[Annotation_mode == "full"] else dt
    if(nrow(dt_f) == 0) return(data.table())
    # Aplicar filtro de genotipo: conservar solo variantes con GT hetero/homo válido
    if("Samples_ID" %in% colnames(dt_f)) {
      gt_vec_prog <- rep(NA_character_, nrow(dt_f))
      sids <- as.character(dt_f$Samples_ID)
      for(sid in unique(sids)) {
        rows <- which(sids == sid)
        if(sid %in% colnames(dt_f)) {
          raw <- as.character(dt_f[[sid]][rows])
          gt_vec_prog[rows] <- sub(":.*", "", raw)
        }
      }
      mask_gt <- gt_vec_prog %in% c("0/1", "1/1", "./1")
      dt_f <- dt_f[mask_gt]
    }
    dt_f
  }
  
  dt_padre_filtrado <- filtrar_progenitor(copy(dt_padre))
  dt_madre_filtrado <- filtrar_progenitor(copy(dt_madre))
  
  # --- FILTRO DE SOLAPAMIENTO CON GenomicRanges ---
  # Métrica de similitud: sim = overlap / max(len_query, len_ref)
  # Equivalente a: sim = 1 - d, donde d = 1 - (overlap / max(len_q, len_r))
  # según: ov <- pintersect(gr[qh], gr[sh])
  #        max_len <- pmax(width(gr[qh]), width(gr[sh]))
  #        d <- 1 - (width(ov) / max_len)
  dt_fase1[, idx_original := .I]
  
  dt_query <- dt_fase1[, .(
    query_id    = idx_original,
    ref_chr     = toupper(trimws(gsub("CHR", "", as.character(SV_chrom)))),
    query_start = as.integer(SV_start),
    query_end   = as.integer(SV_end),
    query_type  = fcase(
      grepl("DEL", toupper(SV_type)), "DEL",
      grepl("DUP", toupper(SV_type)), "DUP",
      default = "OTRO"
    )
  )]
  dt_query[, query_length := query_end - query_start + 1L]
  
  # Construir GRanges para query y referencia
  gr_query <- GRanges(
    seqnames     = dt_query$ref_chr,
    ranges       = IRanges(start = dt_query$query_start, end = dt_query$query_end),
    query_id     = dt_query$query_id,
    query_type   = dt_query$query_type,
    query_length = dt_query$query_length
  )
  
  gr_ref <- GRanges(
    seqnames           = dt_ref$ref_chr,
    ranges             = IRanges(start = dt_ref$ref_start_amplio, end = dt_ref$ref_end_amplio),
    ref_name           = dt_ref$ref_name,
    ref_type           = dt_ref$ref_type,
    ref_length         = dt_ref$ref_end_amplio - dt_ref$ref_start_amplio + 1L,
    ref_start_estricto = dt_ref$ref_start_estricto,
    ref_end_estricto   = dt_ref$ref_end_estricto
  )
  
  hits <- findOverlaps(gr_query, gr_ref, type = "any")
  
  if(length(hits) > 0) {
    q_idx <- queryHits(hits)
    s_idx <- subjectHits(hits)
    
    q_start <- start(gr_query)[q_idx];  q_end <- end(gr_query)[q_idx]
    r_start <- start(gr_ref)[s_idx];    r_end <- end(gr_ref)[s_idx]
    q_len   <- gr_query$query_length[q_idx]
    r_len   <- gr_ref$ref_length[s_idx]
    
    # --- Lógica de overlap: d = 1 - (width(pintersect) / pmax(width_q, width_r)) ---
    # pintersect equivalente vectorizado: ancho de la intersección entre cada par
    ovl_width <- pmax(0L, pmin(r_end, q_end) - pmax(r_start, q_start) + 1L)
    max_len   <- pmax(q_len, r_len)                        # pmax(width_query, width_ref)
    sim_maxlen <- ovl_width / max_len                      # sim = 1 - d = overlap / max_len
    # Jaccard se mantiene como métrica complementaria
    union_len  <- q_len + r_len - ovl_width
    jaccard    <- ovl_width / union_len
    
    dt_overlap <- data.table(
      query_id           = gr_query$query_id[q_idx],
      query_type         = gr_query$query_type[q_idx],
      query_start        = q_start,
      query_end          = q_end,
      query_length       = q_len,
      ref_name           = gr_ref$ref_name[s_idx],
      ref_type           = gr_ref$ref_type[s_idx],
      ref_start_amplio   = r_start,
      ref_end_amplio     = r_end,
      ref_length         = r_len,
      ref_start_estricto = gr_ref$ref_start_estricto[s_idx],
      ref_end_estricto   = gr_ref$ref_end_estricto[s_idx],
      # --- métricas de overlap ---
      overlap_width      = ovl_width,                      # pb de intersección
      prop_referencia    = ovl_width / r_len,              # fracción cubierta del rango de referencia
      sim_maxlen         = sim_maxlen,                     # overlap / max(len_q, len_r)  [métrica principal]
      jaccard            = jaccard,                        # Jaccard = intersección / unión [métrica secundaria]
      desv_izq           = r_start - q_start,   # positivo: query empieza antes que ref (desborde izq)
      desv_der           = q_end   - r_end        # positivo: query acaba después que ref (desborde der)
    )
    
    # Asignar nivel de confianza usando sim_maxlen como criterio principal
    # (equivalente a 1 - d del código de referencia):
    #   Fuerte:   sim >= umbral_herencia
    #   Probable: sim >= umbral_herencia_lax  AND  Jaccard >= umbral_jaccard
    dt_overlap[, Confianza_Region := fcase(
      sim_maxlen >= CONFIG$umbral_herencia,                                                "Fuerte",
      sim_maxlen >= CONFIG$umbral_herencia_lax & jaccard >= CONFIG$umbral_jaccard,         "Probable",
      default = NA_character_
    )]
    
    # CORRECCIÓN: margen_lateral limita el desplazamiento en CUALQUIER dirección.
    # desv_izq y desv_der pueden ser negativos (query fuera del borde de la ref)
    # por lo que hay que aplicar abs() — sin él, valores negativos siempre pasan
    # el filtro `<= margen` aunque el desplazamiento sea enorme.
    dt_valid <- dt_overlap[
      !is.na(Confianza_Region) &
        (query_type == ref_type | ref_type == "DEL/DUP") &
        abs(desv_izq) <= CONFIG$margen_lateral &
        abs(desv_der) <= CONFIG$margen_lateral
    ]
  } else {
    dt_valid <- data.table()
  }
  
  if(nrow(dt_valid) > 0) {
    dt_valid[, `:=`(
      strict_start = pmax(query_start, ref_start_estricto),
      strict_end   = pmin(query_end,   ref_end_estricto)
    )]
    dt_valid[, strict_len := pmax(0, strict_end - strict_start + 1)]
    dt_valid[, pct_strict := strict_len / query_length]
    dt_valid[, Tipo_Rango := fifelse(
      pct_strict >= CONFIG$umbral_strict,
      "Estricto",
      "Amplio"
    )]
    
    # Seleccionar el mejor match por query: Fuerte primero, desempate por mayor sim_maxlen,
    # luego Jaccard, y finalmente ref_name (alfabético) para garantizar determinismo entre batches
    nivel_ord_region <- c(Fuerte = 1L, Probable = 2L)
    dt_valid[, nivel_num := nivel_ord_region[Confianza_Region]]
    dt_best <- dt_valid[order(query_id, nivel_num, -sim_maxlen, -jaccard, ref_name), .SD[1], by = query_id]
    
    # CORRECCIÓN CRÍTICA: dt_best$query_id contiene valores de idx_original (índices
    # de dt_annotsv), NO posiciones de fila en dt_fase1. Usar indexación directa
    # dt_fase1[dt_best$query_id] tomaría posiciones incorrectas cuando dt_fase1 es
    # un subconjunto filtrado de dt_annotsv. Se debe hacer un join por idx_original.
    pos_en_fase1 <- match(dt_best$query_id, dt_fase1$idx_original)
    dt_fase2_region <- dt_fase1[pos_en_fase1]
    dt_fase2_region[, `:=`(
      Referencia_Match  = dt_best$ref_name,
      Overlap_Width     = dt_best$overlap_width,                          # pb solapados
      Prop_Referencia   = round(dt_best$prop_referencia * 100, 1),        # % del rango de referencia cubierto
      Sim_MaxLen        = round(dt_best$sim_maxlen, 4),                   # overlap / max(len_q, len_r)
      Jaccard           = round(dt_best$jaccard, 4),                      # índice de Jaccard
      Confianza_Region  = dt_best$Confianza_Region,                       # Fuerte / Probable
      Tipo_Rango        = dt_best$Tipo_Rango
    )]
    ids_en_region <- dt_best$query_id
  } else {
    dt_fase2_region <- dt_fase1[0]
    ids_en_region   <- integer(0)
  }
  
  todos_idx  <- dt_fase1$idx_original
  fuera_idx  <- setdiff(todos_idx, ids_en_region)
  
  dt_fase2_fuera <- dt_fase1[idx_original %in% fuera_idx & En_Panel_Genes == TRUE]
  if(nrow(dt_fase2_fuera) > 0) {
    dt_fase2_fuera[, `:=`(
      Referencia_Match = NA_character_,
      Overlap_Width    = NA_integer_,
      Prop_Referencia  = NA_real_,
      Sim_MaxLen       = NA_real_,
      Jaccard          = NA_real_,
      Confianza_Region = NA_character_,
      Tipo_Rango       = "Fuera"
    )]
  }
  
  dt_fase2 <- rbindlist(list(dt_fase2_region, dt_fase2_fuera), fill = TRUE)
  
  if(nrow(dt_fase2) == 0) {
    cat("      (Sin variantes tras filtro de regiones + SFARI)\n")
    return(NULL)
  }
  
  # --- ANÁLISIS DE HERENCIA ---
  dt_fase2[, Tipo_Herencia      := NA_character_]
  dt_fase2[, IDs_Coincidencia_PA := NA_character_]
  dt_fase2[, IDs_Coincidencia_MA := NA_character_]
  
  norm_prog <- function(dt) {
    if(nrow(dt) == 0 || !"SV_chrom" %in% colnames(dt)) return(data.table())
    id_col <- if("AnnotSV_ID" %in% colnames(dt)) as.character(dt$AnnotSV_ID) else as.character(seq_len(nrow(dt)))
    
    # Clasificar tipo: cubre tanto SVs (DEL/DUP) como CNVs (loss/gain)
    sv_type_up <- toupper(as.character(dt$SV_type))
    tipo_vec <- fcase(
      grepl("DEL|LOSS", sv_type_up), "DEL",
      grepl("DUP|GAIN", sv_type_up), "DUP",
      default = "OTRO"
    )
    
    # Extraer genotipo del progenitor desde su columna de muestra (GT:FT:... → GT)
    gt_vec <- rep(NA_character_, nrow(dt))
    if("Samples_ID" %in% colnames(dt)) {
      sids <- as.character(dt$Samples_ID)
      n_encontradas <- 0L
      for(sid in unique(sids)) {
        rows <- which(sids == sid)
        if(sid %in% colnames(dt)) {
          raw <- as.character(dt[[sid]][rows])
          gt_vec[rows] <- sub(":.*", "", raw)
          n_encontradas <- n_encontradas + length(rows)
        }
      }
      # CORRECCIÓN: advertir si ningún Samples_ID coincide con una columna de muestra.
      # En ese caso gt_vec queda todo NA, los genotipos en el output serán "?" y la
      # herencia podría clasificarse con información incompleta.
      if(n_encontradas == 0L) {
        cat("      ⚠ [norm_prog] Ningún Samples_ID del progenitor coincide con una columna de muestra.",
            "Genotipos no disponibles para este progenitor.\n")
      }
    }
    
    # Extraer nombres de gen normalizados (pueden ser múltiples, separados por coma/punto y coma)
    gene_raw <- if ("Gene_name" %in% colnames(dt)) {
      toupper(trimws(as.character(dt$Gene_name)))
    } else {
      rep(NA_character_, nrow(dt))
    }
    gene_raw[is.na(gene_raw) | gene_raw == "" | gene_raw == "."] <- NA_character_
    
    dt_norm <- data.table(
      annot_id = id_col,
      chr      = toupper(str_trim(gsub("CHR", "", as.character(dt$SV_chrom)))),
      start    = as.numeric(dt$SV_start),
      end      = as.numeric(dt$SV_end),
      type     = tipo_vec,
      genotype = gt_vec,
      gene     = gene_raw   # genes del progenitor para match por gen+longitud
    )
    # Eliminar filas con coordenadas inválidas o tipo desconocido
    # (tipo "OTRO" generaría falsos matches si dos variantes desconocidas solapan)
    dt_norm <- dt_norm[!is.na(start) & !is.na(end) & type != "OTRO"]
    setkey(dt_norm, chr, start, end)
    return(dt_norm)
  }
  
  calc_herencia <- function(hijo, prog) {
    # Herencia basada exclusivamente en coordenadas genómicas (chr, start, end) y tipo de SV.
    # Métrica principal: sim = overlap / max(len_hijo, len_prog)
    #   equivalente a 1 - d, donde d = 1 - (width(pintersect) / pmax(width_q, width_r))
    # Dos niveles de confianza:
    #   Fuerte:   sim >= umbral_herencia
    #   Probable: sim >= umbral_herencia_lax  AND  Jaccard >= umbral_jaccard
    n <- nrow(hijo)
    resultado <- list(
      heredado  = rep(FALSE, n),
      confianza = rep(NA_character_, n),
      ids_match = rep(NA_character_, n),
      genotype  = rep(NA_character_, n)
    )
    if(nrow(prog) == 0 || n == 0) return(resultado)
    
    gr_hijo <- GRanges(
      seqnames = hijo$chr,
      ranges   = IRanges(start = hijo$start, end = hijo$end),
      row_idx  = hijo$row_idx,
      type     = hijo$type,
      len      = hijo$end - hijo$start + 1L
    )
    gr_prog <- GRanges(
      seqnames = prog$chr,
      ranges   = IRanges(start = prog$start, end = prog$end),
      annot_id = prog$annot_id,
      type     = prog$type,
      genotype = prog$genotype,
      len      = prog$end - prog$start + 1L
    )
    
    hits <- findOverlaps(gr_hijo, gr_prog, type = "any")
    if(length(hits) == 0) return(resultado)
    
    h_idx <- queryHits(hits)
    p_idx <- subjectHits(hits)
    
    # Conservar solo hits con mismo tipo de SV y excluir tipo "OTRO"
    tipo_ok <- gr_hijo$type[h_idx] == gr_prog$type[p_idx] &
      gr_hijo$type[h_idx] != "OTRO"
    if(!any(tipo_ok)) return(resultado)
    h_idx <- h_idx[tipo_ok]
    p_idx <- p_idx[tipo_ok]
    
    # --- Lógica de overlap: d = 1 - (width(pintersect) / pmax(width_hijo, width_prog)) ---
    # Vectorizado sobre los pares de hits filtrados
    h_start <- start(gr_hijo)[h_idx];  h_end <- end(gr_hijo)[h_idx]
    p_start <- start(gr_prog)[p_idx];  p_end <- end(gr_prog)[p_idx]
    h_len   <- gr_hijo$len[h_idx]
    p_len   <- gr_prog$len[p_idx]
    
    # pintersect vectorizado: ancho de la intersección de cada par
    ovl_width <- pmax(0L, pmin(h_end, p_end) - pmax(h_start, p_start) + 1L)
    max_len   <- pmax(h_len, p_len)                        # pmax(width_hijo, width_prog)
    sim_maxlen <- ovl_width / max_len                      # sim = 1 - d
    # Jaccard complementario
    union_len  <- h_len + p_len - ovl_width
    jaccard    <- ovl_width / union_len
    
    ov <- data.table(
      row_idx    = gr_hijo$row_idx[h_idx],
      annot_id   = gr_prog$annot_id[p_idx],
      genotype   = gr_prog$genotype[p_idx],
      sim_maxlen = sim_maxlen,
      jaccard    = jaccard
    )
    
    # Asignar nivel de confianza
    ov[, confianza := fcase(
      sim_maxlen >= CONFIG$umbral_herencia,                                          "Fuerte",
      sim_maxlen >= CONFIG$umbral_herencia_lax & jaccard >= CONFIG$umbral_jaccard,   "Probable",
      default = NA_character_
    )]
    
    ov_ok <- ov[!is.na(confianza)]
    if(nrow(ov_ok) == 0) return(resultado)
    
    # Por cada hijo: elegir el match de mayor confianza; desempate por mayor sim_maxlen
    # y luego por annot_id (alfabético) para garantizar determinismo entre batches
    nivel_ord <- c(Fuerte = 1L, Probable = 2L)
    ov_ok[, nivel_num := nivel_ord[confianza]]
    match_por_hijo <- ov_ok[order(row_idx, nivel_num, -sim_maxlen, -jaccard, annot_id), .SD[1], by = row_idx]
    match_por_hijo <- match_por_hijo[, .(
      heredado  = TRUE,
      confianza = confianza,
      ids_match = annot_id,
      genotype  = genotype
    ), by = row_idx]
    
    pos    <- match(hijo$row_idx, match_por_hijo$row_idx)
    validos <- !is.na(pos)
    resultado$heredado[validos]  <- match_por_hijo$heredado[pos[validos]]
    resultado$confianza[validos] <- match_por_hijo$confianza[pos[validos]]
    resultado$ids_match[validos] <- match_por_hijo$ids_match[pos[validos]]
    resultado$genotype[validos]  <- match_por_hijo$genotype[pos[validos]]
    
    # --- Paso 2: Match por gen + longitud similar (sin solapamiento requerido) ---
    # Se aplica SOLO a filas que no han sido emparejadas en el paso 1 (overlap).
    # Condiciones: mismo cromosoma, mismo tipo de SV, al menos un gen en común,
    # y ratio de longitudes min/max >= umbral_longitud_similar.
    # Todas las coincidencias encontradas aquí se clasifican como "Probable",
    # ya que la ausencia de solapamiento espacial implica menor certeza biológica.
    filas_sin_match <- which(!resultado$heredado)
    
    if (length(filas_sin_match) > 0 &&
        "gene" %in% names(hijo) && "gene" %in% names(prog) && nrow(prog) > 0) {
      
      # Precalcular longitudes del progenitor una sola vez
      prog_len <- prog$end - prog$start + 1L
      
      for (fi in filas_sin_match) {
        h_chr  <- hijo$chr[fi]
        h_type <- hijo$type[fi]
        if (is.na(h_type) || h_type == "OTRO") next
        
        h_len      <- hijo$end[fi] - hijo$start[fi] + 1L
        h_gene_str <- hijo$gene[fi]
        if (is.na(h_gene_str) || h_gene_str == "") next
        
        # Parsear genes del hijo (separadores: coma, punto y coma, espacio)
        h_genes <- toupper(trimws(unlist(strsplit(h_gene_str, "[,;[:space:]]+"))))
        h_genes <- h_genes[nzchar(h_genes) & h_genes != "."]
        if (length(h_genes) == 0) next
        
        # Filtrar candidatos del progenitor: mismo chr y mismo tipo
        mask_base <- prog$chr == h_chr & prog$type == h_type
        if (!any(mask_base)) next
        cands     <- prog[mask_base, ]
        c_len     <- prog_len[mask_base]
        
        # Filtrar por ratio de longitud
        ratio_len <- pmin(h_len, c_len) / pmax(h_len, c_len)
        mask_long <- ratio_len >= CONFIG$umbral_longitud_similar
        if (!any(mask_long)) next
        cands     <- cands[mask_long, ]
        ratio_ok  <- ratio_len[mask_long]
        
        # Filtrar por gen compartido
        p_genes_list <- strsplit(toupper(trimws(cands$gene)), "[,;[:space:]]+")
        gene_match <- vapply(p_genes_list, function(pg) {
          pg <- pg[nzchar(pg) & pg != "."]
          length(pg) > 0L && any(pg %in% h_genes)
        }, logical(1L))
        
        if (!any(gene_match)) next
        cands    <- cands[gene_match, ]
        ratio_ok <- ratio_ok[gene_match]
        
        # Elegir el mejor candidato: mayor ratio de longitud; desempate por annot_id
        best_idx <- order(-ratio_ok, cands$annot_id)[1L]
        
        resultado$heredado[fi]  <- TRUE
        resultado$confianza[fi] <- "Probable"   # siempre Probable: no hay solapamiento directo
        resultado$ids_match[fi] <- cands$annot_id[best_idx]
        resultado$genotype[fi]  <- cands$genotype[best_idx]
        
        cat(paste0("      [Herencia gen+lon] Fila ", fi,
                   " | gen(es)=", paste(h_genes, collapse = ","),
                   " | ratio_len=", round(ratio_ok[best_idx], 3),
                   " → match ", cands$annot_id[best_idx], "\n"))
      }
    }
    
    return(resultado)
  }
  
  if(nrow(dt_padre_filtrado) > 0 || nrow(dt_madre_filtrado) > 0) {
    
    dt_pa_n <- norm_prog(dt_padre_filtrado)
    dt_ma_n <- norm_prog(dt_madre_filtrado)
    
    tiene_modo <- "Annotation_mode" %in% colnames(dt_fase2)
    filas_her_idx <- if(tiene_modo) which(dt_fase2$Annotation_mode == "full") else seq_len(nrow(dt_fase2))
    
    if(length(filas_her_idx) > 0) {
      dt_hijo <- dt_fase2[filas_her_idx, .(
        row_idx = filas_her_idx,
        chr   = toupper(str_trim(gsub("CHR", "", as.character(SV_chrom)))),
        start = as.numeric(SV_start),
        end   = as.numeric(SV_end),
        type  = fcase(
          grepl("DEL|LOSS", toupper(SV_type)), "DEL",
          grepl("DUP|GAIN", toupper(SV_type)), "DUP",
          default = "OTRO"
        ),
        gene  = {
          g <- toupper(trimws(as.character(Gene_name)))
          g[is.na(g) | g == "" | g == "."] <- NA_character_
          g
        }
      )]
      # No se usa setkey aquí: los joins se hacen vía GRanges (no por clave data.table)
      # y reordenar podría desincronizar row_idx con filas_her_idx en depuraciones futuras
      
      res_pa <- calc_herencia(dt_hijo, dt_pa_n)
      res_ma <- calc_herencia(dt_hijo, dt_ma_n)
      
      her_pa  <- res_pa$heredado
      her_ma  <- res_ma$heredado
      conf_pa <- res_pa$confianza
      conf_ma <- res_ma$confianza
      ids_pa  <- res_pa$ids_match
      ids_ma  <- res_ma$ids_match
      gt_pa   <- res_pa$genotype
      gt_ma   <- res_ma$genotype
      
      # TRUE si el archivo del progenitor tiene datos reales (no vacio)
      tiene_padre <- nrow(dt_padre_filtrado) > 0
      tiene_madre <- nrow(dt_madre_filtrado) > 0
      
      # Etiqueta: origen + confianza (si Probable) + genotipo del progenitor entre corchetes
      fmt_gt  <- function(gt) ifelse(!is.na(gt) & nzchar(gt), paste0(" [", gt, "]"), "")
      fmt_duo <- function(gp, gm) paste0(
        " [", ifelse(!is.na(gp) & nzchar(gp), gp, "?"),
        " | ", ifelse(!is.na(gm) & nzchar(gm), gm, "?"), "]"
      )
      
      # Si ninguno de los dos progenitores tiene la variante, la clasificacion final
      # depende de si ambos estaban disponibles:
      #   - Ambos disponibles y ausentes -> De novo (confirmado)
      #   - Alguno falta               -> Desconocida (no se puede confirmar)
      etiqueta_sin_herencia <- ifelse(tiene_padre & tiene_madre, "De novo", "Desconocida")
      
      tipo_her <- ifelse(
        her_pa & her_ma & !is.na(conf_pa) & conf_pa == "Fuerte" & !is.na(conf_ma) & conf_ma == "Fuerte",
        paste0("Conjunta", fmt_duo(gt_pa, gt_ma)),
        ifelse(
          her_pa & her_ma,
          paste0("Conjunta (Probable)", fmt_duo(gt_pa, gt_ma)),
          ifelse(
            her_pa & !is.na(conf_pa) & conf_pa == "Fuerte",
            paste0("Paterna", fmt_gt(gt_pa)),
            ifelse(
              her_pa,
              paste0("Paterna (Probable)", fmt_gt(gt_pa)),
              ifelse(
                her_ma & !is.na(conf_ma) & conf_ma == "Fuerte",
                paste0("Materna", fmt_gt(gt_ma)),
                ifelse(
                  her_ma,
                  paste0("Materna (Probable)", fmt_gt(gt_ma)),
                  etiqueta_sin_herencia
                )
              )
            )
          )
        )
      )
      
      dt_fase2[filas_her_idx, Tipo_Herencia      := tipo_her]
      dt_fase2[filas_her_idx, IDs_Coincidencia_PA := ids_pa]
      dt_fase2[filas_her_idx, IDs_Coincidencia_MA := ids_ma]
      
      # Propagar de filas full → split por AnnotSV_ID de forma segura
      if(tiene_modo && "AnnotSV_ID" %in% colnames(dt_fase2)) {
        dt_her_map <- dt_fase2[Annotation_mode == "full", .(AnnotSV_ID, Tipo_Herencia, IDs_Coincidencia_PA, IDs_Coincidencia_MA)]
        dt_her_map <- unique(dt_her_map, by = "AnnotSV_ID") 
        
        dt_fase2[Annotation_mode == "split", c("Tipo_Herencia", "IDs_Coincidencia_PA", "IDs_Coincidencia_MA") := {
          idx <- match(AnnotSV_ID, dt_her_map$AnnotSV_ID)
          list(
            dt_her_map$Tipo_Herencia[idx],
            dt_her_map$IDs_Coincidencia_PA[idx],
            dt_her_map$IDs_Coincidencia_MA[idx]
          )
        }]
      }
    }
  } else {
    cat("      ⚠ Sin archivos de progenitores disponibles; herencia no calculada\n")
  }
  
  # --- GENOTIPO DEL HIJO (probando) ---
  # Extraer GT del hijo igual que se hace con los progenitores: usando Samples_ID
  # como nombre de la columna de muestra y tomando el campo antes del primer ":".
  # Se aplica solo a filas "full"; las "split" lo heredan por AnnotSV_ID.
  if ("Samples_ID" %in% colnames(dt_fase2)) {
    gt_hijo_vec <- rep(NA_character_, nrow(dt_fase2))
    sids_hijo   <- as.character(dt_fase2$Samples_ID)
    for (sid in unique(sids_hijo)) {
      rows <- which(sids_hijo == sid)
      if (sid %in% colnames(dt_fase2)) {
        raw <- as.character(dt_fase2[[sid]][rows])
        gt_hijo_vec[rows] <- sub(":.*", "", raw)
      }
    }
    dt_fase2[, Genotipo_Hijo := gt_hijo_vec]
    
    # Propagar de full -> split por AnnotSV_ID (igual que Tipo_Herencia)
    if ("Annotation_mode" %in% colnames(dt_fase2) && "AnnotSV_ID" %in% colnames(dt_fase2)) {
      map_gt <- dt_fase2[Annotation_mode == "full", .(AnnotSV_ID, Genotipo_Hijo)]
      map_gt <- unique(map_gt, by = "AnnotSV_ID")
      dt_fase2[Annotation_mode == "split", Genotipo_Hijo := {
        idx <- match(AnnotSV_ID, map_gt$AnnotSV_ID)
        map_gt$Genotipo_Hijo[idx]
      }]
    }
  } else {
    dt_fase2[, Genotipo_Hijo := NA_character_]
  }
  
  # --- ORDENAMIENTO ---
  if("Annotation_mode" %in% colnames(dt_fase2)) {
    dt_fase2[, is_full := fifelse(Annotation_mode == "full", 0L, 1L)]
  } else {
    dt_fase2[, is_full := fifelse(grepl("[;,]", Gene_name), 0L, 1L)]
  }
  
  dt_final <- dt_fase2[order(AnnotSV_ID, is_full, -ranking_numeric)]
  dt_final[, is_full := NULL]
  
  # --- DEDUPLICACIÓN DE REGIONES SIMILARES ---
  # Colapsa variantes que solapan con sim >= umbral_herencia en el mismo chr y tipo.
  # Si alguna tiene herencia distinta de "De novo", se conserva esa; si todas son
  # "De novo", se queda la de mayor ranking_numeric.
  # Activar/desactivar con CONFIG$deduplicar_regiones.
  if(isTRUE(CONFIG$deduplicar_regiones)) {
    n_antes_dedup <- nrow(dt_final)
    dt_final <- deduplicar_regiones_similares(dt_final)
    n_despues_dedup <- nrow(dt_final)
    if(n_antes_dedup != n_despues_dedup) {
      cat(paste0("      [Dedup] Total: ", n_antes_dedup, " → ", n_despues_dedup, " variantes\n"))
    }
  } else {
    cat("      [Dedup] Colapso desactivado (CONFIG$deduplicar_regiones = FALSE)\n")
  }
  
  # --- RESTAURAR COLUMNAS ORIGINALES CON VALORES REALES ---
  # Para cada columna del TSV original ausente en dt_final, recuperar los valores
  # reales desde el snapshot dt_annotsv_src usando .src_row_ como clave de fila.
  # Esto garantiza que no solo aparezca el nombre sino el valor correcto.
  cols_aux_internas <- c("ranking_numeric", "idx_original", "is_full", ".src_row_")
  cols_esperadas    <- setdiff(cols_originales, cols_aux_internas)
  cols_perdidas     <- setdiff(cols_esperadas, names(dt_final))
  if (length(cols_perdidas) > 0) {
    cat(paste0("      ⚠ Recuperando ", length(cols_perdidas), " columna(s) desde snapshot: ",
               paste(cols_perdidas, collapse = ", "), "\n"))
    pos <- match(dt_final$.src_row_, dt_annotsv_src$.src_row_)
    for (col in cols_perdidas) {
      dt_final[, (col) := dt_annotsv_src[[col]][pos]]
    }
  }
  
  cols_inicio  <- c("Tipo_Herencia", "Genotipo_Hijo", "Tipo_Rango", "En_Panel_Genes", "AnnotSV_ranking_score")
  cols_fin     <- c("IDs_Coincidencia_PA", "IDs_Coincidencia_MA")
  # Eliminar columnas auxiliares internas que no deben aparecer en el Excel
  cols_aux     <- c("idx_original", "ranking_numeric", ".src_row_")
  for (col_aux in cols_aux) {
    if (col_aux %in% names(dt_final)) dt_final[, (col_aux) := NULL]
  }
  cols_inicio_ex <- intersect(cols_inicio, names(dt_final))
  cols_fin_ex    <- intersect(cols_fin,    names(dt_final))
  cols_medio     <- setdiff(names(dt_final), c(cols_inicio, cols_fin))
  
  setcolorder(dt_final, c(cols_inicio_ex, cols_medio, cols_fin_ex))
  
  cat(paste0("      ✓ ", nrow(dt_final), " variantes finales\n"))
  return(dt_final)
}

# =============================================================================
# 5. FUNCIÓN OPTIMIZADA PARA GUARDAR EXCEL
# =============================================================================
sanitizar_df <- function(df) {
  for(j in seq_along(df)) {
    if(is.list(df[[j]])) {
      df[[j]] <- vapply(df[[j]], function(x) {
        if(is.null(x) || length(x) == 0) NA_character_
        else paste(x, collapse = ";")
      }, character(1L))
    }
  }
  for(j in seq_along(df)) {
    if(is.logical(df[[j]])) {
      df[[j]] <- ifelse(is.na(df[[j]]), NA_character_,
                        ifelse(df[[j]], "Sí", "No"))
    }
  }
  for(j in seq_along(df)) {
    if(is.numeric(df[[j]])) {
      df[[j]][!is.finite(df[[j]])] <- NA_real_
    }
  }
  for(j in seq_along(df)) {
    if(is.character(df[[j]])) {
      df[[j]] <- iconv(df[[j]], from = "", to = "UTF-8", sub = "")
    }
  }
  for(j in seq_along(df)) {
    if(is.factor(df[[j]])) df[[j]] <- as.character(df[[j]])
  }
  names(df) <- iconv(names(df), from = "", to = "UTF-8", sub = "")
  names(df) <- gsub("[\\[\\]\\*\\?:/\\\\]", "_", names(df), perl = TRUE)
  # Resolver colisiones de nombres: columnas de AnnotSV como FORMAT[GT], FORMAT[FT], etc.
  # quedan con nombres tipo FORMAT_GT_ tras el gsub. Si dos columnas distintas producen
  # el mismo nombre, writeData() descarta silenciosamente las duplicadas.
  # Se añade sufijo _2, _3... para garantizar unicidad.
  nms <- names(df)
  for (nm in unique(nms[duplicated(nms)])) {
    idx <- which(nms == nm)
    nms[idx] <- paste0(nm, seq_along(idx))
  }
  names(df) <- nms
  df
}

guardar_excel_estilizado <- function(lista_resultados, ruta_salida) {
  wb <- createWorkbook()
  
  st_green  <- createStyle(fgFill = "#C6EFCE", fontColour = "#006100")
  st_orange <- createStyle(fgFill = "#FFEB9C", fontColour = "#9C6500")
  st_red    <- createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006")
  st_denovo <- createStyle(fgFill = "#FFE6E6", fontColour = "#CC0000", textDecoration = "bold")
  st_paterna<- createStyle(fgFill = "#E6F3FF", fontColour = "#0066CC")
  st_materna<- createStyle(fgFill = "#FFF0E6", fontColour = "#CC6600")
  st_conj   <- createStyle(fgFill = "#F0E6FF", fontColour = "#6600CC")
  st_estric <- createStyle(fgFill = "#D4EDDA", fontColour = "#155724", textDecoration = "bold")
  st_amplio <- createStyle(fgFill = "#FFF3CD", fontColour = "#856404")
  st_fuera  <- createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006", textDecoration = "bold")
  
  for(nombre_hoja in names(lista_resultados)) {
    df <- sanitizar_df(as.data.frame(lista_resultados[[nombre_hoja]]))
    cat(paste0("      [Excel] Hoja '", nombre_hoja, "': ", ncol(df), " columnas, ", nrow(df), " filas\n"))
    addWorksheet(wb, nombre_hoja)
    writeData(wb, nombre_hoja, df)
    
    col_her <- which(names(df) == "Tipo_Herencia")
    col_ran <- which(names(df) == "Tipo_Rango")
    col_sfa <- which(names(df) == "En_Panel_Genes")
    col_sco <- which(names(df) == "AnnotSV_ranking_score")
    
    if(length(col_her) > 0) {
      rows_denovo   <- which(!is.na(df$Tipo_Herencia) & df$Tipo_Herencia %in% c("De novo", "Desconocida"))
      rows_paterna  <- which(!is.na(df$Tipo_Herencia) & grepl("^Paterna",  df$Tipo_Herencia))
      rows_materna  <- which(!is.na(df$Tipo_Herencia) & grepl("^Materna",  df$Tipo_Herencia))
      rows_conjunta <- which(!is.na(df$Tipo_Herencia) & grepl("^Conjunta", df$Tipo_Herencia))
      
      if(length(rows_denovo) > 0)   addStyle(wb, nombre_hoja, st_denovo,  rows = rows_denovo + 1,   cols = col_her)
      if(length(rows_paterna) > 0)  addStyle(wb, nombre_hoja, st_paterna, rows = rows_paterna + 1,  cols = col_her)
      if(length(rows_materna) > 0)  addStyle(wb, nombre_hoja, st_materna, rows = rows_materna + 1,  cols = col_her)
      if(length(rows_conjunta) > 0) addStyle(wb, nombre_hoja, st_conj,    rows = rows_conjunta + 1, cols = col_her)
    }
    
    if(length(col_ran) > 0) {
      rows_estricto <- which(!is.na(df$Tipo_Rango) & df$Tipo_Rango == "Estricto")
      rows_amplio   <- which(!is.na(df$Tipo_Rango) & df$Tipo_Rango == "Amplio")
      rows_fuera    <- which(!is.na(df$Tipo_Rango) & df$Tipo_Rango == "Fuera")
      
      if(length(rows_estricto) > 0) addStyle(wb, nombre_hoja, st_estric, rows = rows_estricto + 1, cols = col_ran)
      if(length(rows_amplio) > 0)   addStyle(wb, nombre_hoja, st_amplio, rows = rows_amplio + 1,   cols = col_ran)
      if(length(rows_fuera) > 0)    addStyle(wb, nombre_hoja, st_fuera,  rows = rows_fuera + 1,    cols = col_ran)
    }
    
    if(length(col_sfa) > 0) {
      rows_sfari_yes <- which(!is.na(df[[col_sfa]]) & df[[col_sfa]] == "Sí")
      rows_sfari_no  <- which(!is.na(df[[col_sfa]]) & df[[col_sfa]] == "No")
      
      if(length(rows_sfari_yes) > 0) addStyle(wb, nombre_hoja, st_green, rows = rows_sfari_yes + 1, cols = col_sfa)
      if(length(rows_sfari_no) > 0)  addStyle(wb, nombre_hoja, st_red,   rows = rows_sfari_no + 1,  cols = col_sfa)
    }
    
    if(length(col_sco) > 0) {
      # AnnotSV_ranking_score: negativos = benignos (verde), 0–0.5 = intermedio (naranja), ≥0.5 = patogénico (rojo)
      score_vals <- suppressWarnings(as.numeric(df$AnnotSV_ranking_score))
      rows_verde   <- which(!is.na(score_vals) & score_vals < 0)
      rows_naranja <- which(!is.na(score_vals) & score_vals >= 0   & score_vals < 0.5)
      rows_rojo    <- which(!is.na(score_vals) & score_vals >= 0.5)
      
      if(length(rows_verde) > 0)   addStyle(wb, nombre_hoja, st_green,  rows = rows_verde + 1,   cols = col_sco)
      if(length(rows_naranja) > 0) addStyle(wb, nombre_hoja, st_orange, rows = rows_naranja + 1, cols = col_sco)
      if(length(rows_rojo) > 0)    addStyle(wb, nombre_hoja, st_red,    rows = rows_rojo + 1,    cols = col_sco)
    }
    
    setColWidths(wb, nombre_hoja, cols = 1:ncol(df), widths = "auto")
  }
  
  saveWorkbook(wb, ruta_salida, overwrite = TRUE)
}

# =============================================================================
# 6. BUCLE PRINCIPAL (BATCH PROCESSING)
# =============================================================================
cat("\n=============================================================================\n")
cat("INICIANDO PROCESAMIENTO EN LOTES\n")
cat("=============================================================================\n")

sel <- function(tsv, patron) tsv[grepl(patron, basename(tsv), ignore.case = TRUE)][1]

EXCEPCIONES_PROGENITORES <- list(
  # Formato: id_familia = list(PA = id_progenitor_pa, MA = id_progenitor_ma)
  # El código buscará un directorio cuyo nombre termine en paste0(id_progenitor, sufijo)
  # p.ej. para PA = "G04PMP068" buscará un directorio que acabe en "G04PMP068PA".
  # ⚠ Asegurarse de que el directorio del progenitor exista con ese sufijo exacto.
  "G04PMP069" = list(PA = "G04PMP068", MA = "G04PMP068")
)

buscar_dir_progenitor <- function(id_familia, sufijo, directorios_sub) {
  if(id_familia %in% names(EXCEPCIONES_PROGENITORES)) {
    id_prog <- EXCEPCIONES_PROGENITORES[[id_familia]][[sufijo]]
    if(!is.null(id_prog)) {
      excepcion <- directorios_sub[grepl(paste0(id_prog, sufijo, "$"), basename(directorios_sub), ignore.case = TRUE)][1]
      if(!is.na(excepcion)) {
        cat(paste0("  ℹ Excepción de progenitor aplicada para ", id_familia, sufijo,
                   ": ", basename(excepcion), "
"))
        return(excepcion)
      }
    }
  }
  
  exacto <- directorios_sub[grepl(paste0(id_familia, sufijo, "$"), directorios_sub, ignore.case = TRUE)][1]
  if(!is.na(exacto)) return(exacto)
  
  prefijo       <- sub("\\d+$", "", id_familia)       
  num_hijo      <- sub("^.*?(\\d+)$", "\\1", id_familia)   
  num_hijo_sin0 <- sub("^0+", "", num_hijo)           
  
  patron_compartido <- paste0(
    "^", prefijo, 
    "(\\d+-(", num_hijo, "|", num_hijo_sin0, ")|(", num_hijo, "|", num_hijo_sin0, ")-\\d+)", 
    sufijo, "$"
  )
  
  compartido <- directorios_sub[grepl(patron_compartido, basename(directorios_sub), ignore.case = TRUE)][1]
  
  if(!is.na(compartido)) {
    cat(paste0("  ℹ Progenitor compartido detectado para ", id_familia, sufijo, ": ", basename(compartido), "\n"))
    return(compartido)
  }
  return(NA)
}

procesar_familia <- function(id_familia, directorios_sub, ruta_salida_sub) {
  
  dir_hi <- directorios_sub[grepl(paste0(id_familia, "HI$"), directorios_sub, ignore.case = TRUE)][1]
  dir_pa <- buscar_dir_progenitor(id_familia, "PA", directorios_sub)
  dir_ma <- buscar_dir_progenitor(id_familia, "MA", directorios_sub)
  
  if(is.na(dir_hi) || !dir.exists(dir_hi)) {
    cat("  ⚠ Sin carpeta HI, omitiendo familia\n")
    return(FALSE)
  }
  
  tsv_hi <- get_tsv_cached(dir_hi)
  tsv_pa <- get_tsv_cached(dir_pa)
  tsv_ma <- get_tsv_cached(dir_ma)
  
  if(length(tsv_hi) == 0) return(FALSE)
  
  f_cnv_hi <- sel(tsv_hi, "HI[.]CNVs")
  f_cnv_pa <- sel(tsv_pa, "PA[.]CNVs")
  f_cnv_ma <- sel(tsv_ma, "MA[.]CNVs")
  f_sv_hi  <- sel(tsv_hi, "HI[.]SVs")
  f_sv_pa  <- sel(tsv_pa, "PA[.]SVs")
  f_sv_ma  <- sel(tsv_ma, "MA[.]SVs")
  
  resultados_familia <- list()
  
  if(!is.na(f_cnv_hi)) {
    res_cnv <- tryCatch(
      procesar_modalidad(f_cnv_hi, f_cnv_pa, f_cnv_ma, "CNV", dt_ref, genes_panel),
      error = function(e) { cat(paste("    [ERROR CNV]:", e$message, "\n")); NULL }
    )
    if(!is.null(res_cnv)) {
      res_cnv <- data.table(ID_Familia = id_familia, res_cnv)
      resultados_familia[["CNVs"]] <- res_cnv
    }
  }
  
  if(!is.na(f_sv_hi)) {
    res_sv <- tryCatch(
      procesar_modalidad(f_sv_hi, f_sv_pa, f_sv_ma, "SV", dt_ref, genes_panel),
      error = function(e) { cat(paste("    [ERROR SV]:", e$message, "\n")); NULL }
    )
    if(!is.null(res_sv)) {
      res_sv <- data.table(ID_Familia = id_familia, res_sv)
      resultados_familia[["SVs"]] <- res_sv
    }
  }
  
  if(length(resultados_familia) > 0) {
    if(!dir.exists(ruta_salida_sub)) dir.create(ruta_salida_sub, recursive = TRUE)
    nombre_excel  <- paste0(id_familia, "_Analisis_Completo.xlsx")
    ruta_completa <- file.path(ruta_salida_sub, nombre_excel)
    tryCatch({
      guardar_excel_estilizado(resultados_familia, ruta_completa)
      cat(paste("  ✓ Excel creado:", file.path(basename(ruta_salida_sub), nombre_excel), "\n"))
      return(resultados_familia)
    }, error = function(e) {
      cat(paste("  [ERROR guardando Excel]:", e$message, "\n"))
      return(FALSE)
    })
  } else {
    cat("  ⚠ Sin resultados para esta familia\n")
    return(FALSE)
  }
}

subcarpetas <- list.dirs(CONFIG$ruta_entrada, full.names = TRUE, recursive = FALSE)

# Cache de rutas TSV: se define ANTES de construir tareas y de cualquier llamada
# a procesar_familia, para evitar referencias a objetos no definidos en cualquier
# orden de ejecución o entorno de worker paralelo.
cache_tsv <- new.env(hash = TRUE, parent = emptyenv())
get_tsv_cached <- function(dir) {
  if(is.na(dir) || !dir.exists(dir)) return(character(0))
  key <- normalizePath(dir, mustWork = FALSE)
  if(exists(key, envir = cache_tsv, inherits = FALSE)) return(get(key, envir = cache_tsv))
  f <- list.files(dir, pattern = "annotated\\.tsv$", full.names = TRUE, recursive = TRUE)
  assign(key, f, envir = cache_tsv)
  f
}

tareas <- rbindlist(lapply(subcarpetas, function(subcarpeta) {
  directorios_sub <- list.dirs(subcarpeta, full.names = TRUE, recursive = FALSE)
  if(length(directorios_sub) == 0) return(NULL)
  ids_dir <- basename(directorios_sub)
  
  # CORRECCIÓN: derivar IDs de familia SOLO desde directorios HI.
  # Antes se usaban todos los directorios (HI + PA + MA), lo que provocaba que
  # familias sin carpeta HI (solo PA o MA) se registrasen igualmente como tareas.
  # Ahora: una familia solo existe si tiene un directorio HI; los progenitores
  # se localizan después mediante buscar_dir_progenitor, que ya maneja tanto
  # progenitores únicos como compartidos (patrón nnn-mmm).
  hi_dirs     <- ids_dir[grepl("HI$", ids_dir, ignore.case = TRUE)]
  ids_familia <- unique(sub("HI$", "", hi_dirs, ignore.case = TRUE))
  # Excluir: cadenas vacías (salvaguarda) y cualquier ID residual con patrón
  # de progenitor compartido (dígitos-dígitos), aunque en la práctica no
  # deberían aparecer porque los directorios compartidos terminan en PA/MA, no HI.
  ids_familia <- ids_familia[ids_familia != "" & !grepl("^[^-]+-\\d+(HI|PA|MA)?$", ids_familia, ignore.case = TRUE)]
  if(length(ids_familia) == 0) return(NULL)
  data.table(
    subcarpeta      = subcarpeta,
    nombre_sub      = basename(subcarpeta),
    id_familia      = ids_familia,
    ruta_salida_sub = file.path(CONFIG$ruta_salida, basename(subcarpeta))
  )
}), fill = TRUE)

procesar_tarea <- function(tarea) {
  id_familia      <- tarea$id_familia
  subcarpeta      <- tarea$subcarpeta
  ruta_salida_sub <- tarea$ruta_salida_sub
  directorios_sub <- list.dirs(subcarpeta, full.names = TRUE, recursive = FALSE)
  
  cat(paste0("  ---------------------------------------------------\n"))
  cat(paste0("  [", tarea$nombre_sub, "] Familia: ", id_familia, "\n"))
  
  exito <- procesar_familia(id_familia, directorios_sub, ruta_salida_sub)
  return(exito)
}

if(.Platform$OS.type == "windows" || CONFIG$n_cores == 1L) {
  resultados <- lapply(seq_len(nrow(tareas)), function(i) procesar_tarea(tareas[i]))
} else {
  resultados <- mclapply(seq_len(nrow(tareas)), function(i) procesar_tarea(tareas[i]),
                         mc.cores = CONFIG$n_cores, mc.preschedule = FALSE)
}

n_subcarpetas <- length(unique(tareas$nombre_sub))
n_familias    <- nrow(tareas)
n_exitosas    <- sum(vapply(resultados, is.list, logical(1))) 

cat("=============================================================================\n")
cat("✓ PROCESAMIENTO FINALIZADO\n")
cat(paste0("  Subcarpetas procesadas: ", n_subcarpetas, "\n"))
cat(paste0("  Familias procesadas:    ", n_familias,    "\n"))
cat(paste0("  Archivos generados:     ", n_exitosas,    "\n"))
cat("=============================================================================\n")

# =============================================================================
# 6. RECOPILACIÓN DE VARIANTES DE ALTO IMPACTO (Score >= 0.15)
# =============================================================================
cat("\n=============================================================================\n")
cat("GENERANDO DOCUMENTO RECOPILATORIO (SCORE >= 0.15)\n")
cat("=============================================================================\n")

lista_cnvs <- list()
lista_svs  <- list()

for (res in resultados) {
  if (is.list(res)) {
    if (!is.null(res$CNVs)) lista_cnvs[[length(lista_cnvs) + 1]] <- res$CNVs
    if (!is.null(res$SVs))  lista_svs[[length(lista_svs) + 1]] <- res$SVs
  }
}

dt_todas_cnv <- rbindlist(lista_cnvs, fill = TRUE)
dt_todas_sv  <- rbindlist(lista_svs, fill = TRUE)

filtrar_alto_impacto <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  
  # ranking_numeric es una columna auxiliar interna que se elimina antes del output;
  # aquí se recalcula desde AnnotSV_ranking_score para no depender de ella.
  score_num <- suppressWarnings(as.numeric(as.character(dt$AnnotSV_ranking_score)))
  
  if ("Annotation_mode" %in% colnames(dt)) {
    dt_filtrado <- dt[dt$Annotation_mode == "full" & !is.na(score_num) & score_num >= 0.15]
  } else {
    dt_filtrado <- dt[!is.na(score_num) & score_num >= 0.15]
  }
  
  if(nrow(dt_filtrado) > 0) {
    score_ord <- suppressWarnings(as.numeric(as.character(dt_filtrado$AnnotSV_ranking_score)))
    dt_filtrado <- dt_filtrado[order(dt_filtrado$ID_Familia, -score_ord)]
  }
  
  return(dt_filtrado)
}

dt_alto_cnv <- filtrar_alto_impacto(dt_todas_cnv)
dt_alto_sv  <- filtrar_alto_impacto(dt_todas_sv)

ruta_recopilacion <- file.path(CONFIG$ruta_salida, "Recopilacion_Alto_Impacto")

if (nrow(dt_alto_cnv) > 0 || nrow(dt_alto_sv) > 0) {
  if(!dir.exists(ruta_recopilacion)) dir.create(ruta_recopilacion, recursive = TRUE)
  
  resultados_recopilados <- list()
  if (nrow(dt_alto_cnv) > 0) resultados_recopilados[["CNVs_Prioridad"]] <- dt_alto_cnv
  if (nrow(dt_alto_sv) > 0)  resultados_recopilados[["SVs_Prioridad"]] <- dt_alto_sv
  
  ruta_excel_recop <- file.path(ruta_recopilacion, "Variantes_Score_Mayor_0.15.xlsx")
  tryCatch({
    guardar_excel_estilizado(resultados_recopilados, ruta_excel_recop)
    cat(paste("✓ Documento recopilatorio creado exitosamente en:\n  ", ruta_excel_recop, "\n"))
    cat(paste("  - Filas de CNVs incluidas:", nrow(dt_alto_cnv), "\n"))
    cat(paste("  - Filas de SVs incluidas: ", nrow(dt_alto_sv), "\n"))
  }, error = function(e) {
    cat(paste("❌ Error guardando documento recopilatorio:", e$message, "\n"))
  })
} else {
  cat("⚠ No se encontraron variantes con AnnotSV_ranking_score >= 0.15 en toda la cohorte.\n")
}

cat("=============================================================================\n")  
cat("=============================================================================\n")