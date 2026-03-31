library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(openxlsx)
library(parallel)          # parallel processing of families
library(GenomicRanges)     # precise overlap: width, prop_reference, Jaccard

# =============================================================================
# 1. GENERAL CONFIGURATION
# =============================================================================
# CONFIG is only defined here if it has NOT already been set by the Shiny app.
# When launched from the app, the app prepends its own CONFIG block (with the
# user-supplied paths) before this script, so this block is skipped entirely.
# When running the script standalone, this default block is used instead.
if (!exists("CONFIG")) {
  CONFIG <- list(
    ruta_entrada    = "",                                  # Set to your input folder (with HI/PA/MA subfolders per family),
    ruta_salida     = file.path(getwd(), "results"),       # Output folder (created automatically),
    ruta_rangos     = file.path(getwd(), "reference", "reference_ranges.xlsx"),  # Reference ranges file,
    ruta_panel      = file.path(getwd(), "reference", "gene_panel.txt"),         # Gene panel file (.txt, comma-separated),
    umbral_strict        = 0.70,
    umbral_herencia      = 0.70,   # minimum sim for "Strong" inheritance   [sim = ovl / max(len_q, len_r)]
    umbral_herencia_lax  = 0.60,   # minimum sim for "Probable" inheritance (rescued by Jaccard)
    umbral_jaccard       = 0.50,   # minimum Jaccard required at the "Probable" level
    margen_lateral       = 3e6,    # 3 Mb maximum allowed deviation per side
    umbral_longitud_similar = 0.70,  # minimum ratio min(len_h,len_p)/max(len_h,len_p) for gene-based match without overlap
    umbral_contencion       = 0.80,  # minimum fraction of the PROBAND variant covered by the parent (containment match)
    # Detects: proband CNV fully inside a larger parent CNV — sim_maxlen would be too low but
    # the proband is clearly inherited. e.g. proband 100 kb inside parent 400 kb →
    # sim_maxlen = 0.25 (would fail), but prop_hijo = 1.0 (clearly inherited).
    deduplicar_regiones  = TRUE,   # TRUE = collapse similar variants (same chr+type, sim >= umbral_herencia)
    # FALSE = keep all variants without collapsing
    n_cores         = max(1L, detectCores(logical = FALSE) - 1L)  # physical cores - 1
  )
}
if(!dir.exists(CONFIG$ruta_salida)) dir.create(CONFIG$ruta_salida, recursive = TRUE)

# =============================================================================
# 2. LOADING REFERENCES (WITH VALIDATION)
# =============================================================================
cat(">>> Loading reference databases...\n")

# Load Ranges
if(!file.exists(CONFIG$ruta_rangos)) stop("❌ File not found: ", CONFIG$ruta_rangos)

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
cat("✓ Ranges loaded:", nrow(dt_ref), "\n")

# Load SFARI Genes
if(!file.exists(CONFIG$ruta_panel)) stop("❌ File not found: ", CONFIG$ruta_panel)

lineas_panel <- readLines(CONFIG$ruta_panel, warn = FALSE)
genes_panel <- unique(toupper(trimws(unlist(strsplit(paste(lineas_panel, collapse = ","), ",")))))
genes_panel <- genes_panel[genes_panel != ""]
# Hash set: converts the vector into an environment for O(1) lookup instead of O(n)
genes_panel_set <- new.env(hash = TRUE, parent = emptyenv())
for(.g in genes_panel) assign(.g, TRUE, envir = genes_panel_set)
rm(.g)
cat("✓ Panel genes loaded:", length(genes_panel), "\n")

# =============================================================================
# 3. DEDUPLICATION FUNCTION FOR SIMILAR REGIONS
# =============================================================================
# Logic: within each family result, if two or more variants overlap
# with sim_maxlen >= umbral_herencia on the same chromosome and same SV/CNV type,
# they are considered "the same region". From that group, the following are kept:
#   - Rows with Tipo_Herencia other than "De novo", if any exist.
#   - If all are "De novo" or NA, the one with the highest ranking_numeric is kept.
# "split" rows always follow their corresponding "full" row.
#
# Special inheritance cases:
#   - Paternal + Maternal similar → merged into Combined (neither is discarded)
#   - Non-mergeable different inheritances (e.g. Paternal vs De novo) → NOT collapsed,
#     they are biologically distinct events even if their coordinates overlap
#   - Same inheritance or both De novo/NA → normal collapse, highest-ranking one is kept
#
# Greedy algorithm by direct pairs (not Union-Find):
#   A variant is only removed if it has a DIRECT sim >= threshold with its representative.
#   Avoids transitive collapses: A~B and B~C does not imply collapsing A with C if sim(A,C)<threshold.
# ---------------------------------------------------------------------------

# Helper function: extract the base inheritance family (without genotype or "(Probable)")
.her_base <- function(x) {
  if (is.na(x)) return(NA_character_)
  # Remove the suffix " [GT]" or " [GT | GT]" and "(Probable)"
  x <- trimws(sub("\\s*\\[.*", "", x))
  x <- trimws(sub("\\s*\\(Probable\\)", "", x))
  x
}

deduplicar_regiones_similares <- function(dt, umbral_sim = CONFIG$umbral_herencia) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  
  tiene_modo  <- "Annotation_mode" %in% colnames(dt)
  tiene_annot <- "AnnotSV_ID"      %in% colnames(dt)
  
  # Work only on "full" rows; "split" rows are reconstructed at the end by AnnotSV_ID
  dt_full <- if (tiene_modo) dt[Annotation_mode == "full"] else copy(dt)
  if (nrow(dt_full) <= 1) return(dt)
  
  # --- Build per-row metadata ---
  n         <- nrow(dt_full)
  fn_row    <- seq_len(n)
  sv_up     <- toupper(as.character(dt_full$SV_type))
  tipo_vec  <- fcase(
    grepl("DEL|LOSS", sv_up), "DEL",
    grepl("DUP|GAIN", sv_up), "DUP",
    default = "OTHER"
  )
  start_vec <- as.integer(dt_full$SV_start)
  end_vec   <- as.integer(dt_full$SV_end)
  len_vec   <- end_vec - start_vec + 1L
  chr_vec   <- toupper(trimws(gsub("CHR", "", as.character(dt_full$SV_chrom))))
  
  her_vec   <- as.character(dt_full$Tipo_Herencia)
  her_base  <- vapply(her_vec, .her_base, character(1L))  # base family without decoration
  is_denovo <- !is.na(her_vec) & her_base %in% c("De novo", "Unknown")
  rank_num  <- suppressWarnings(as.numeric(dt_full$ranking_numeric))
  rank_num[is.na(rank_num)] <- -Inf
  tiebreak  <- if (tiene_annot) as.character(dt_full$AnnotSV_ID) else as.character(fn_row)
  
  # --- Compute all pairs with overlap ---
  gr <- GRanges(
    seqnames = chr_vec,
    ranges   = IRanges(start = start_vec, end = end_vec)
  )
  
  hits <- findOverlaps(gr, gr, type = "any")
  hits <- hits[queryHits(hits) < subjectHits(hits)]
  if (length(hits) == 0) return(dt)
  
  hi <- queryHits(hits)
  hj <- subjectHits(hits)
  
  # Filter: same known type
  tipo_ok <- tipo_vec[hi] == tipo_vec[hj] & tipo_vec[hi] != "OTHER"
  hi <- hi[tipo_ok];  hj <- hj[tipo_ok]
  if (length(hi) == 0) return(dt)
  
  # Compute sim
  ovl_w <- pmax(0L, pmin(end_vec[hj], end_vec[hi]) - pmax(start_vec[hj], start_vec[hi]) + 1L)
  sim_v <- ovl_w / pmax(len_vec[hi], len_vec[hj])
  
  # Filter by threshold
  ok    <- sim_v >= umbral_sim
  hi    <- hi[ok];  hj <- hj[ok];  sim_v <- sim_v[ok]
  if (length(hi) == 0) return(dt)
  
  # Sort: sim desc, deterministic tiebreak
  ord   <- order(-sim_v, tiebreak[hi], tiebreak[hj])
  hi    <- hi[ord];  hj <- hj[ord];  sim_v <- sim_v[ord]
  
  # --- Greedy: process each pair ---
  eliminado <- logical(n)
  # fusionar[i] = index j whose inheritance must be merged into i (Paternal+Maternal→Combined case)
  fusionar  <- integer(n)
  
  for (k in seq_along(hi)) {
    i <- hi[k];  j <- hj[k]
    if (eliminado[i] || eliminado[j]) next
    
    hi_base <- her_base[i]
    hj_base <- her_base[j]
    
    # ---- Case 1: Paternal + Maternal (or vice versa) → merge into Combined ----
    # They are the same variant seen from two different parents.
    # i is kept, j is marked for removal, and i is noted to also absorb
    # j's information to build the Combined label.
    es_fusion_conjunta <- (
      (!is.na(hi_base) & !is.na(hj_base)) &&
        ((hi_base == "Paternal" & hj_base == "Maternal") ||
           (hi_base == "Maternal" & hj_base == "Paternal"))
    )
    if (es_fusion_conjunta) {
      # Ensure i is Paternal and j is Maternal to unify merge logic
      if (hi_base == "Maternal") { tmp <- i; i <- j; j <- tmp }
      eliminado[j]  <- TRUE
      fusionar[i]   <- j    # i (Paternal) absorbs j (Maternal) → Combined
      next
    }
    
    # ---- Case 2: non-mergeable different inheritances → DO NOT collapse ----
    # Examples: Paternal vs De novo, Maternal vs Combined, etc.
    # They are biologically distinct events; both are kept.
    her_i_conocida <- !is.na(hi_base) & !hi_base %in% c("De novo", "Unknown")
    her_j_conocida <- !is.na(hj_base) & !hj_base %in% c("De novo", "Unknown")
    if (her_i_conocida && her_j_conocida && hi_base != hj_base) {
      cat(paste0("      [Dedup] Pair not collapsed: different inheritances (",
                 hi_base, " vs ", hj_base, ") with sim=", round(sim_v[k], 4), "\n"))
      next
    }
    
    # ---- Case 3: normal collapse ----
    # Same inheritance, both De novo/NA, or one known and the other without information.
    # Priority: (1) non-denovo > denovo, (2) higher rank, (3) deterministic tiebreak
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
  
  cat(paste0("      [Dedup] Collapsed ", n_eliminadas,
             " variant(s) (same chr+type, sim >= ", umbral_sim, ")\n"))
  
  # --- Apply Combined merges ---
  # For each i with fusionar[i] > 0: i is Paternal, j=fusionar[i] is Maternal.
  # Build the Combined label by combining genotypes from both.
  for (i in which(fusionar > 0L)) {
    j       <- fusionar[i]
    her_i   <- her_vec[i]   # "Paternal [GT]" or "Paternal (Probable) [GT]"
    her_j   <- her_vec[j]   # "Maternal [GT]" or "Maternal (Probable) [GT]"
    
    # Extract genotype from each label: what is between [ ]
    gt_pa  <- regmatches(her_i, regexpr("(?<=\\[)[^\\]]+(?=\\])", her_i, perl = TRUE))
    gt_ma  <- regmatches(her_j, regexpr("(?<=\\[)[^\\]]+(?=\\])", her_j, perl = TRUE))
    gt_pa  <- if (length(gt_pa) > 0) gt_pa else "?"
    gt_ma  <- if (length(gt_ma) > 0) gt_ma else "?"
    
    # Confidence level: if either is Probable → Combined (Probable)
    es_probable <- grepl("Probable", her_i, fixed = TRUE) ||
      grepl("Probable", her_j, fixed = TRUE)
    etiqueta <- if (es_probable) {
      paste0("Combined (Probable) [", gt_pa, " | ", gt_ma, "]")
    } else {
      paste0("Combined [", gt_pa, " | ", gt_ma, "]")
    }
    
    dt_full[fn_row == i, Tipo_Herencia := etiqueta]
    
    # Also merge match IDs if the columns exist
    if ("IDs_Coincidencia_PA" %in% colnames(dt_full) &&
        "IDs_Coincidencia_MA" %in% colnames(dt_full)) {
      # i was Paternal: has IDs_PA but not IDs_MA → take IDs_MA from j (which was Maternal)
      ids_ma_j <- dt_full$IDs_Coincidencia_MA[fn_row == j]
      ids_pa_i <- dt_full$IDs_Coincidencia_PA[fn_row == i]
      # If j (Maternal) had IDs_PA because it was also partially Combined, keep both
      ids_pa_j <- dt_full$IDs_Coincidencia_PA[fn_row == j]
      pa_final <- paste(na.omit(unique(c(ids_pa_i, ids_pa_j))), collapse = ";")
      pa_final <- if (nchar(pa_final) == 0) NA_character_ else pa_final
      ma_final <- if (!is.na(ids_ma_j) && nchar(ids_ma_j) > 0) ids_ma_j else NA_character_
      dt_full[fn_row == i, `:=`(IDs_Coincidencia_PA = pa_final,
                                IDs_Coincidencia_MA = ma_final)]
    }
  }
  
  # --- Reconstruct result ---
  filas_conservar <- fn_row[!eliminado]
  
  if (tiene_modo && tiene_annot) {
    ids_conservar <- dt_full$AnnotSV_ID[filas_conservar]
    # Propagate merged Tipo_Herencia to the corresponding "split" rows
    dt_result <- dt[AnnotSV_ID %in% ids_conservar]
    # Update inheritance in dt_result for IDs that were merged
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
# 4. OPTIMISED PROCESSING FUNCTION
# =============================================================================
procesar_modalidad <- function(ruta_hi, ruta_pa, ruta_ma, MODO, dt_ref, genes_panel) {
  
  # --- INPUT VALIDATION ---
  if(is.na(ruta_hi) || !file.exists(ruta_hi)) return(NULL)
  
  cat(paste0("   -> Analysing ", MODO, "...\n"))
  
  # --- OPTIMISED LOADING ---
  # fill=TRUE: AnnotSV can generate fields with complex content (commas, semicolons,
  # line breaks in OMIM phenotypes, etc.) that shift columns in specific rows.
  # fill=TRUE tolerates rows with different numbers of fields instead of aborting or shifting.
  dt_annotsv <- tryCatch({
    fread(ruta_hi, header = TRUE, quote = "", fill = TRUE, showProgress = FALSE)
  }, error = function(e) {
    cat(paste("      [ERROR reading file]:", e$message, "\n"))
    return(NULL)
  })
  
  if(is.null(dt_annotsv) || nrow(dt_annotsv) == 0) return(NULL)
  
  # Full snapshot: row index + all original columns, BEFORE any filter.
  # At the end of the pipeline it is used to recover the real values of any column
  # that was lost in intermediate steps.
  dt_annotsv[, .src_row_ := .I]
  dt_annotsv_src  <- copy(dt_annotsv)        # all columns, all values
  cols_originales <- setdiff(names(dt_annotsv), ".src_row_")
  cat(paste0("      Columns read: ", length(cols_originales), "\n"))
  
  # Load parents safely
  dt_padre <- if(!is.na(ruta_pa) && file.exists(ruta_pa)) {
    tryCatch(fread(ruta_pa, header=TRUE, quote="", fill=TRUE, showProgress=FALSE),
             error = function(e) data.table())
  } else data.table()
  
  dt_madre <- if(!is.na(ruta_ma) && file.exists(ruta_ma)) {
    tryCatch(fread(ruta_ma, header=TRUE, quote="", fill=TRUE, showProgress=FALSE),
             error = function(e) data.table())
  } else data.table()
  
  dt_annotsv[, ranking_numeric := suppressWarnings(as.numeric(as.character(AnnotSV_ranking_score)))]
  
  # --- POPULATION FREQUENCY FILTER (FIRST FILTER) ---
  # Condition 1: frequency in benign CNV databases by SV type
  #   DEL → B_loss_AFmax | DUP → B_gain_AFmax | INS → B_ins_AFmax | INV → B_inv_AFmax
  #   Passes if the cell is empty/NA or the value is <= 0.01
  # Condition 2: frequency in the Illumina DRAGEN cohort (similar counts)
  #   Value of Illumina_DRAGEN.similar.counts divided by 1061
  #   Passes if the cell is empty/NA or the quotient is <= 0.01
  # Condition 3: exact count in the Illumina DRAGEN cohort (exact counts)
  #   Direct value of Illumina_DRAGEN.exact.counts (integer, no parsing or normalisation)
  #   Passes if the cell is empty/NA or the value is <= 10
  # Rows that meet ALL THREE conditions are kept (empty/NA cells pass)
  {
    # SV type per row (vectorised)
    sv_type_upper <- toupper(as.character(dt_annotsv$SV_type))
    tipo_sv <- fcase(
      grepl("DEL", sv_type_upper), "DEL",
      grepl("DUP", sv_type_upper), "DUP",
      grepl("INS", sv_type_upper), "INS",
      grepl("INV", sv_type_upper), "INV",
      default = "OTHER"
    )
    
    # Type → AF column map
    col_af_map <- c(DEL = "B_loss_AFmax", DUP = "B_gain_AFmax",
                    INS = "B_ins_AFmax",  INV = "B_inv_AFmax")
    
    # Extract AF value by type for each row
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
    
    # Condition 2: Illumina DRAGEN similar counts / 1061
    # The column has format "DEL=2343;DUP=9485;INS=748;INV=7374;TRA=738"
    # The value for the SV type of each row is extracted, summed with possible
    # sub-types (e.g. all containing "DEL") and divided by 1061.
    # Passes if the cell is empty/NA or the quotient is <= 0.01
    if("Illumina_DRAGEN.similar.counts" %in% colnames(dt_annotsv)) {
      dragen_str <- as.character(dt_annotsv$Illumina_DRAGEN.similar.counts)
      
      # For each row: parse the key=value string and sum values whose
      # key contains the SV type of that row (e.g. type "DEL" matches
      # keys "DEL", "DEL_ALU", etc.)
      dragen_counts <- mapply(function(sv_str, sv_tipo) {
        # Empty cell or unknown type → NA (passes the filter)
        if(is.na(sv_str) || sv_str == "" || sv_str == "." || sv_tipo == "OTHER")
          return(NA_real_)
        # Split key=value pairs
        pares <- unlist(strsplit(sv_str, ";", fixed = TRUE))
        # Extract key and value from each pair
        claves  <- sub("=.*$", "", pares)
        valores <- suppressWarnings(as.numeric(sub("^[^=]+=", "", pares)))
        # Sum only values whose key contains the row's SV type
        idx_match <- grepl(sv_tipo, claves, fixed = TRUE)
        if(!any(idx_match)) return(NA_real_)
        sum(valores[idx_match], na.rm = TRUE)
      }, dragen_str, tipo_sv, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      
      dragen_freq <- as.numeric(dragen_counts) / 1061
      mask_dragen <- is.na(dragen_freq) | dragen_freq <= 0.01
    } else {
      mask_dragen <- rep(TRUE, nrow(dt_annotsv))   # column absent → passes
    }
    
    # Condition 3: Illumina DRAGEN exact counts
    # Direct count (integer) — no parsing or normalisation required.
    # Passes if the cell is empty/NA or the value is <= 10.
    if("Illumina_DRAGEN.exact.counts" %in% colnames(dt_annotsv)) {
      exact_counts <- suppressWarnings(
        as.numeric(as.character(dt_annotsv$Illumina_DRAGEN.exact.counts))
      )
      mask_exact <- is.na(exact_counts) | exact_counts <= 10
    } else {
      mask_exact <- rep(TRUE, nrow(dt_annotsv))   # column absent → passes
    }
    
    mask_freq   <- mask_af & mask_dragen & mask_exact
    n_antes     <- nrow(dt_annotsv)
    dt_annotsv  <- dt_annotsv[mask_freq]
    cat(paste0("      Frequency filter: ", n_antes, " → ", nrow(dt_annotsv), " variants\n"))
  }
  
  if(nrow(dt_annotsv) == 0) {
    cat("      (No variants after frequency filter)\n")
    return(NULL)
  }
  
  # --- MINIMUM SIZE FILTER (SVs only) ---
  # Applied only in SV mode. Removes variants whose absolute size
  # is less than 50 bp (DEL, DUP, INS, INV).
  # Size is calculated as abs(SV_end - SV_start). Absolute value is used
  # because in some notations SV_end can be less than SV_start (e.g. INS).
  # Rows with invalid (NA) coordinates are kept: the subsequent quality filter
  # will discard them without generating silent losses.
  if(MODO == "SV") {
    sv_size  <- abs(as.numeric(dt_annotsv$SV_end) - as.numeric(dt_annotsv$SV_start))
    mask_tam <- is.na(sv_size) | sv_size >= 50
    n_antes_tam     <- nrow(dt_annotsv)
    dt_annotsv      <- dt_annotsv[mask_tam]
    n_filtradas_tam <- n_antes_tam - nrow(dt_annotsv)
    cat(paste0("      Size filter (<50 bp): ", n_antes_tam, " → ",
               nrow(dt_annotsv), " variants",
               if(n_filtradas_tam > 0) paste0(" (removed: ", n_filtradas_tam, ")") else "",
               "\n"))
    if(nrow(dt_annotsv) == 0) {
      cat("      (No variants after size filter)\n")
      return(NULL)
    }
  }
  
  # --- SFARI GENE FLAGGING (FULLY VECTORISED) ---
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
  
  # --- OPTIMISED QUALITY FILTER ---
  if(!"Samples_ID" %in% colnames(dt_annotsv)) {
    cat("      \u26a0 Column Samples_ID not found\n")
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
    # QUALITY filter: passes if the column does not exist, the value is NA, or the value > 30.
    # Computed once and applied in both branches (CNV and SV).
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
      # Note: filtered by Annotation_mode == "full" here for SVs because "split" rows
      # do not have GT:FT format in their sample column and should not be evaluated
      # for quality individually — their inclusion is inherited from the associated "full" row.
      return((gt_part %in% c("0/1", "1/1", "./1")) &
               (ft_part %in% c("PASS", "High Quality")) &   # consistent with CNV mode
               (filter_col %in% c("PASS", "High Quality")) &
               (annot_col  == "full") &
               mask_quality_score)
    }
  }
  
  mask_quality <- extraer_mask_calidad(dt_annotsv, MODO)
  dt_fase1 <- dt_annotsv[mask_quality]
  
  if(nrow(dt_fase1) == 0) {
    cat("      (No variants after quality filter)\n")
    return(NULL)
  }
  
  # Parent filter: works only on "full" rows to avoid counting the same
  # SV multiple times. Additionally, a genotype filter is applied
  # (heterozygous / homozygous) to exclude low-quality or missing-call variants,
  # reducing false-positive inheritance across batches.
  filtrar_progenitor <- function(dt) {
    if(nrow(dt) == 0) return(data.table())
    dt_f <- if("Annotation_mode" %in% colnames(dt)) dt[Annotation_mode == "full"] else dt
    if(nrow(dt_f) == 0) return(data.table())
    # Apply genotype filter: keep only variants with valid het/hom GT
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
      # GT filter: keep only variants with valid het/hom calls.
      # Includes 1/0 and 1/. which are equivalent to 0/1 and ./1 in phased/unphased VCFs.
      mask_gt <- gt_vec_prog %in% c("0/1", "1/0", "1/1", "./1", "1/.")
      dt_f <- dt_f[mask_gt]
    }
    dt_f
  }
  
  dt_padre_filtrado <- filtrar_progenitor(copy(dt_padre))
  dt_madre_filtrado <- filtrar_progenitor(copy(dt_madre))
  
  # --- OVERLAP FILTER WITH GenomicRanges ---
  # Similarity metric: sim = overlap / max(len_query, len_ref)
  # Equivalent to: sim = 1 - d, where d = 1 - (overlap / max(len_q, len_r))
  # as per: ov <- pintersect(gr[qh], gr[sh])
  #         max_len <- pmax(width(gr[qh]), width(gr[sh]))
  #         d <- 1 - (width(ov) / max_len)
  dt_fase1[, idx_original := .I]
  
  dt_query <- dt_fase1[, .(
    query_id    = idx_original,
    ref_chr     = toupper(trimws(gsub("CHR", "", as.character(SV_chrom)))),
    query_start = as.integer(SV_start),
    query_end   = as.integer(SV_end),
    query_type  = fcase(
      grepl("DEL", toupper(SV_type)), "DEL",
      grepl("DUP", toupper(SV_type)), "DUP",
      default = "OTHER"
    )
  )]
  dt_query[, query_length := query_end - query_start + 1L]
  
  # Build GRanges for query and reference
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
    
    # --- Overlap logic: d = 1 - (width(pintersect) / pmax(width_q, width_r)) ---
    # Vectorised pintersect equivalent: width of intersection for each pair
    ovl_width <- pmax(0L, pmin(r_end, q_end) - pmax(r_start, q_start) + 1L)
    max_len   <- pmax(q_len, r_len)                        # pmax(width_query, width_ref)
    sim_maxlen <- ovl_width / max_len                      # sim = 1 - d = overlap / max_len
    # Jaccard kept as complementary metric
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
      # --- overlap metrics ---
      overlap_width      = ovl_width,                      # bp of intersection
      prop_referencia    = ovl_width / r_len,              # fraction of the reference range covered
      sim_maxlen         = sim_maxlen,                     # overlap / max(len_q, len_r)  [main metric]
      jaccard            = jaccard,                        # Jaccard = intersection / union [secondary metric]
      desv_izq           = r_start - q_start,   # positive: query starts before ref (left overhang)
      desv_der           = q_end   - r_end       # positive: query ends after ref (right overhang)
    )
    
    # Assign confidence level using sim_maxlen as the primary criterion
    # (equivalent to 1 - d from the reference code):
    #   Strong:   sim >= umbral_herencia
    #   Probable: sim >= umbral_herencia_lax  AND  Jaccard >= umbral_jaccard
    dt_overlap[, Confianza_Region := fcase(
      sim_maxlen >= CONFIG$umbral_herencia,                                                "Strong",
      sim_maxlen >= CONFIG$umbral_herencia_lax & jaccard >= CONFIG$umbral_jaccard,         "Probable",
      default = NA_character_
    )]
    
    # FIX: margen_lateral limits displacement in ANY direction.
    # desv_izq and desv_der can be negative (query outside the ref boundary)
    # so abs() must be applied — without it, negative values always pass
    # the `<= margin` filter even when displacement is enormous.
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
      "Strict",
      "Wide"
    )]
    
    # Select the best match per query: Strong first, tiebreak by highest sim_maxlen,
    # then Jaccard, and finally ref_name (alphabetical) to guarantee determinism across batches
    nivel_ord_region <- c(Strong = 1L, Probable = 2L)
    dt_valid[, nivel_num := nivel_ord_region[Confianza_Region]]
    dt_best <- dt_valid[order(query_id, nivel_num, -sim_maxlen, -jaccard, ref_name), .SD[1], by = query_id]
    
    # CRITICAL FIX: dt_best$query_id contains values of idx_original (indices
    # of dt_annotsv), NOT row positions in dt_fase1. Using direct indexing
    # dt_fase1[dt_best$query_id] would take incorrect positions when dt_fase1 is
    # a filtered subset of dt_annotsv. A join by idx_original must be used.
    pos_en_fase1 <- match(dt_best$query_id, dt_fase1$idx_original)
    dt_fase2_region <- dt_fase1[pos_en_fase1]
    dt_fase2_region[, `:=`(
      Referencia_Match  = dt_best$ref_name,
      Overlap_Width     = dt_best$overlap_width,                          # overlapping bp
      Prop_Referencia   = round(dt_best$prop_referencia * 100, 1),        # % of reference range covered
      Sim_MaxLen        = round(dt_best$sim_maxlen, 4),                   # overlap / max(len_q, len_r)
      Jaccard           = round(dt_best$jaccard, 4),                      # Jaccard index
      Confianza_Region  = dt_best$Confianza_Region,                       # Strong / Probable
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
      Tipo_Rango       = "Outside"
    )]
  }
  
  dt_fase2 <- rbindlist(list(dt_fase2_region, dt_fase2_fuera), fill = TRUE)
  
  if(nrow(dt_fase2) == 0) {
    cat("      (No variants after region + SFARI filter)\n")
    return(NULL)
  }
  
  # --- INHERITANCE ANALYSIS ---
  # ── Pre-extract proband genotype ─────────────────────────────────────────────
  # MUST happen BEFORE the inheritance block so the rescue pass (which checks
  # proband GT to determine rescue priority) has the column available.
  # The later Genotipo_Hijo block (post-sorting) propagates it to "split" rows
  # and is kept for that purpose; it will not re-extract (guarded by column check).
  if ("Samples_ID" %in% colnames(dt_fase2) && !"Genotipo_Hijo" %in% colnames(dt_fase2)) {
    gt_hijo_early <- rep(NA_character_, nrow(dt_fase2))
    sids_early    <- as.character(dt_fase2$Samples_ID)
    for (sid in unique(sids_early)) {
      rows <- which(sids_early == sid)
      if (sid %in% colnames(dt_fase2)) {
        raw <- as.character(dt_fase2[[sid]][rows])
        gt_hijo_early[rows] <- sub(":.*", "", raw)
      }
    }
    dt_fase2[, Genotipo_Hijo := gt_hijo_early]
    rm(gt_hijo_early, sids_early)
  }
  
  dt_fase2[, Tipo_Herencia      := NA_character_]
  dt_fase2[, IDs_Coincidencia_PA := NA_character_]
  dt_fase2[, IDs_Coincidencia_MA := NA_character_]
  
  norm_prog <- function(dt) {
    if(nrow(dt) == 0 || !"SV_chrom" %in% colnames(dt)) return(data.table())
    id_col <- if("AnnotSV_ID" %in% colnames(dt)) as.character(dt$AnnotSV_ID) else as.character(seq_len(nrow(dt)))
    
    # Classify type: covers CNVs (DEL/DUP via loss/gain) and SVs (INS, INV).
    # INS and INV are now recognised so their inheritance can be detected.
    # BND/TRA remain "OTHER" — they lack simple coordinate semantics for overlap.
    sv_type_up <- toupper(as.character(dt$SV_type))
    tipo_vec <- fcase(
      grepl("DEL|LOSS",    sv_type_up), "DEL",
      grepl("DUP|GAIN",    sv_type_up), "DUP",
      grepl("INS|INSERT",  sv_type_up), "INS",
      grepl("INV|INVERS",  sv_type_up), "INV",
      default = "OTHER"
    )
    
    # Extract parent genotype from its sample column (GT:FT:... → GT)
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
      # FIX: warn if no Samples_ID matches a sample column.
      # In that case gt_vec will be all NA, genotypes in the output will be "?" and
      # inheritance could be classified with incomplete information.
      if(n_encontradas == 0L) {
        cat("      \u26a0 [norm_prog] No parent Samples_ID matches a sample column.",
            "Genotypes not available for this parent.\n")
      }
    }
    
    # Extract normalised gene names (can be multiple, separated by comma/semicolon)
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
      gene     = gene_raw   # parent genes for gene+length matching
    )
    # Remove rows with invalid coordinates or unknown type
    # (type "OTHER" would generate false matches if two unknown variants overlap)
    dt_norm <- dt_norm[!is.na(start) & !is.na(end) & type != "OTHER"]
    setkey(dt_norm, chr, start, end)
    return(dt_norm)
  }
  
  # ---------------------------------------------------------------------------
  # calc_herencia_rescue — relaxed-threshold pass for suspicious De novo calls
  # ---------------------------------------------------------------------------
  # Used ONLY as a second pass on proband variants that:
  #   (a) were not matched by calc_herencia, AND
  #   (b) have a homozygous proband genotype (1/1)
  #
  # A truly de novo homozygous CNV/SV is biologically implausible:
  # two independent identical de novo events on both alleles is essentially
  # impossible.  The most likely explanation is that one or both parents carry
  # the same variant with slightly different breakpoints (caller variability),
  # causing the regular sim_maxlen to fall below CONFIG$umbral_herencia_lax.
  #
  # Key difference from calc_herencia:
  #   sim = pmax(overlap/max(len), overlap/min(len))
  # The second term ("contained" metric) rescues the common case where a large
  # parent CNV fully contains a smaller proband CNV:
  #   e.g. proband 100 kb, parent 160 kb, overlap 95 kb
  #   → sim_maxlen = 0.59  (fails Probable threshold of 0.60)
  #   → sim_min    = 0.95  (passes rescue threshold)
  #
  # All matches from this pass are always labelled "Probable", never "Strong".
  # ---------------------------------------------------------------------------
  UMBRAL_RESCUE  <- CONFIG$umbral_herencia_lax * 0.80   # e.g. 0.60 * 0.80 = 0.48
  JACCARD_RESCUE <- CONFIG$umbral_jaccard       * 0.70   # e.g. 0.50 * 0.70 = 0.35
  
  calc_herencia_rescue <- function(hijo, prog,
                                   umbral_sim = UMBRAL_RESCUE,
                                   umbral_jac = JACCARD_RESCUE) {
    n <- nrow(hijo)
    resultado <- list(
      heredado  = rep(FALSE, n),
      confianza = rep(NA_character_, n),
      ids_match = rep(NA_character_, n),
      genotype  = rep(NA_character_, n)
    )
    if (nrow(prog) == 0 || n == 0) return(resultado)
    
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
    if (length(hits) == 0) return(resultado)
    
    h_idx <- queryHits(hits)
    p_idx <- subjectHits(hits)
    
    tipo_ok <- gr_hijo$type[h_idx] == gr_prog$type[p_idx] &
      gr_hijo$type[h_idx] != "OTHER"
    if (!any(tipo_ok)) return(resultado)
    h_idx <- h_idx[tipo_ok]
    p_idx <- p_idx[tipo_ok]
    
    h_start <- start(gr_hijo)[h_idx];  h_end <- end(gr_hijo)[h_idx]
    p_start <- start(gr_prog)[p_idx];  p_end <- end(gr_prog)[p_idx]
    h_len   <- gr_hijo$len[h_idx];     p_len <- gr_prog$len[p_idx]
    
    ovl_w <- pmax(0L, pmin(h_end, p_end) - pmax(h_start, p_start) + 1L)
    
    # Use the more generous of sim_maxlen and sim_contained.
    # sim_contained catches "large parent contains small proband" cases where
    # sim_maxlen would fail the regular threshold despite very high actual overlap.
    sim_maxlen    <- ovl_w / pmax(h_len, p_len)
    sim_contained <- ovl_w / pmin(h_len, p_len)
    sim_v         <- pmax(sim_maxlen, sim_contained)
    jac_v         <- ovl_w / (h_len + p_len - ovl_w)
    
    ov <- data.table(
      row_idx  = gr_hijo$row_idx[h_idx],
      annot_id = gr_prog$annot_id[p_idx],
      genotype = gr_prog$genotype[p_idx],
      sim      = sim_v,
      jac      = jac_v
    )
    ov_ok <- ov[sim >= umbral_sim & jac >= umbral_jac]
    if (nrow(ov_ok) == 0) return(resultado)
    
    best <- ov_ok[order(row_idx, -sim, -jac, annot_id), .SD[1], by = row_idx]
    pos  <- match(hijo$row_idx, best$row_idx)
    val  <- !is.na(pos)
    resultado$heredado[val]  <- TRUE
    resultado$confianza[val] <- "Probable"   # always Probable in rescue pass
    resultado$ids_match[val] <- best$annot_id[pos[val]]
    resultado$genotype[val]  <- best$genotype[pos[val]]
    return(resultado)
  }
  
  calc_herencia <- function(hijo, prog) {
    # Inheritance based exclusively on genomic coordinates (chr, start, end) and SV type.
    # Main metric: sim = overlap / max(len_proband, len_parent)
    #   equivalent to 1 - d, where d = 1 - (width(pintersect) / pmax(width_q, width_r))
    # Two confidence levels:
    #   Strong:   sim >= umbral_herencia
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
    
    # Keep only hits with same SV type and exclude type "OTHER"
    tipo_ok <- gr_hijo$type[h_idx] == gr_prog$type[p_idx] &
      gr_hijo$type[h_idx] != "OTHER"
    if(!any(tipo_ok)) return(resultado)
    h_idx <- h_idx[tipo_ok]
    p_idx <- p_idx[tipo_ok]
    
    # --- Overlap logic: d = 1 - (width(pintersect) / pmax(width_proband, width_parent)) ---
    # Vectorised over filtered hit pairs
    h_start <- start(gr_hijo)[h_idx];  h_end <- end(gr_hijo)[h_idx]
    p_start <- start(gr_prog)[p_idx];  p_end <- end(gr_prog)[p_idx]
    h_len   <- gr_hijo$len[h_idx]
    p_len   <- gr_prog$len[p_idx]
    
    # Vectorised pintersect: width of intersection for each pair
    ovl_width <- pmax(0L, pmin(h_end, p_end) - pmax(h_start, p_start) + 1L)
    max_len   <- pmax(h_len, p_len)                        # pmax(width_proband, width_parent)
    sim_maxlen <- ovl_width / max_len                      # sim = 1 - d
    # prop_hijo: fraction of the PROBAND covered by the parent.
    # Key for containment: proband 100 kb inside parent 400 kb →
    #   sim_maxlen = 0.25 (fails all thresholds), but prop_hijo = 1.0 (clearly inherited).
    # This metric is NOT used for sim_maxlen-style symmetric matching; it is only
    # triggered when prop_hijo >= CONFIG$umbral_contencion (default 0.80).
    prop_hijo  <- ovl_width / h_len
    # Complementary Jaccard
    union_len  <- h_len + p_len - ovl_width
    jaccard    <- ovl_width / union_len
    
    ov <- data.table(
      row_idx    = gr_hijo$row_idx[h_idx],
      annot_id   = gr_prog$annot_id[p_idx],
      genotype   = gr_prog$genotype[p_idx],
      sim_maxlen = sim_maxlen,
      prop_hijo  = prop_hijo,
      jaccard    = jaccard
    )
    
    # Assign confidence level.
    # Three routes to inheritance:
    #   (1) Strong:              sim_maxlen >= umbral_herencia (symmetric, both agree)
    #   (2) Probable (symmetric): sim_maxlen >= umbral_herencia_lax AND Jaccard >= umbral_jaccard
    #   (3) Probable (containment): prop_hijo >= umbral_contencion
    #       → Catches "proband CNV contained within larger parent CNV" where sim_maxlen
    #         is artificially low due to size asymmetry, even though the entire proband
    #         variant is biologically accounted for by the parent variant.
    ov[, confianza := fcase(
      sim_maxlen >= CONFIG$umbral_herencia,                                         "Strong",
      sim_maxlen >= CONFIG$umbral_herencia_lax & jaccard >= CONFIG$umbral_jaccard,  "Probable",
      prop_hijo  >= CONFIG$umbral_contencion,                                        "Probable",
      default = NA_character_
    )]
    
    ov_ok <- ov[!is.na(confianza)]
    if(nrow(ov_ok) == 0) return(resultado)
    
    # For each proband: choose the highest-confidence match; tiebreak by highest sim_maxlen
    # and then annot_id (alphabetical) to guarantee determinism across batches
    nivel_ord <- c(Strong = 1L, Probable = 2L)
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
    
    # --- Step 2: Match by gene + similar length (no overlap required) ---
    # Applied ONLY to rows that were not matched in step 1 (overlap).
    # Conditions: same chromosome, same SV type, at least one gene in common,
    # and length ratio min/max >= umbral_longitud_similar.
    # All matches found here are classified as "Probable",
    # since the absence of spatial overlap implies lower biological certainty.
    filas_sin_match <- which(!resultado$heredado)
    
    if (length(filas_sin_match) > 0 &&
        "gene" %in% names(hijo) && "gene" %in% names(prog) && nrow(prog) > 0) {
      
      # Pre-compute parent lengths once
      prog_len <- prog$end - prog$start + 1L
      
      for (fi in filas_sin_match) {
        h_chr  <- hijo$chr[fi]
        h_type <- hijo$type[fi]
        if (is.na(h_type) || h_type == "OTHER") next
        
        h_len      <- hijo$end[fi] - hijo$start[fi] + 1L
        h_gene_str <- hijo$gene[fi]
        if (is.na(h_gene_str) || h_gene_str == "") next
        
        # Parse proband genes (separators: comma, semicolon, space)
        h_genes <- toupper(trimws(unlist(strsplit(h_gene_str, "[,;[:space:]]+"))))
        h_genes <- h_genes[nzchar(h_genes) & h_genes != "."]
        if (length(h_genes) == 0) next
        
        # Filter parent candidates: same chr and same type
        mask_base <- prog$chr == h_chr & prog$type == h_type
        if (!any(mask_base)) next
        cands     <- prog[mask_base, ]
        c_len     <- prog_len[mask_base]
        
        # Filter by length ratio
        ratio_len <- pmin(h_len, c_len) / pmax(h_len, c_len)
        mask_long <- ratio_len >= CONFIG$umbral_longitud_similar
        if (!any(mask_long)) next
        cands     <- cands[mask_long, ]
        ratio_ok  <- ratio_len[mask_long]
        
        # Filter by shared gene
        p_genes_list <- strsplit(toupper(trimws(cands$gene)), "[,;[:space:]]+")
        gene_match <- vapply(p_genes_list, function(pg) {
          pg <- pg[nzchar(pg) & pg != "."]
          length(pg) > 0L && any(pg %in% h_genes)
        }, logical(1L))
        
        if (!any(gene_match)) next
        cands    <- cands[gene_match, ]
        ratio_ok <- ratio_ok[gene_match]
        
        # Choose the best candidate: highest length ratio; tiebreak by annot_id
        best_idx <- order(-ratio_ok, cands$annot_id)[1L]
        
        resultado$heredado[fi]  <- TRUE
        resultado$confianza[fi] <- "Probable"   # always Probable: no direct overlap
        resultado$ids_match[fi] <- cands$annot_id[best_idx]
        resultado$genotype[fi]  <- cands$genotype[best_idx]
        
        cat(paste0("      [Gene+len inheritance] Row ", fi,
                   " | gene(s)=", paste(h_genes, collapse = ","),
                   " | ratio_len=", round(ratio_ok[best_idx], 3),
                   " \u2192 match ", cands$annot_id[best_idx], "\n"))
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
          grepl("DEL|LOSS",   toupper(SV_type)), "DEL",
          grepl("DUP|GAIN",   toupper(SV_type)), "DUP",
          grepl("INS|INSERT", toupper(SV_type)), "INS",
          grepl("INV|INVERS", toupper(SV_type)), "INV",
          default = "OTHER"
        ),
        gene  = {
          g <- toupper(trimws(as.character(Gene_name)))
          g[is.na(g) | g == "" | g == "."] <- NA_character_
          g
        }
      )]
      # setkey is not used here: joins are done via GRanges (not by data.table key)
      # and reordering could desync row_idx from filas_her_idx in future debugging
      
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
      
      # --- RESCUE PASS for unmatched proband variants ---
      # Genotipo_Hijo was extracted BEFORE this block (early extraction fix), so it
      # is always available here regardless of execution order.
      gt_hijo_rescue <- if ("Genotipo_Hijo" %in% colnames(dt_fase2)) {
        dt_fase2$Genotipo_Hijo[filas_her_idx]
      } else {
        rep(NA_character_, length(filas_her_idx))
      }
      
      # es_hom_denovo: kept solely for the De_novo_Quality diagnostic flag below.
      # It is NO LONGER used to gate the rescue pass (rescue is now extended to all
      # unmatched variants — see es_candidato_rescue).
      es_hom_denovo <- (!her_pa & !her_ma) &
        !is.na(gt_hijo_rescue) & gt_hijo_rescue == "1/1"
      # Previously only triggered for 1/1 proband genotype. Extended to ALL variants
      # not matched by the regular pass, because:
      #   (a) Breakpoint caller variability affects both 0/1 and 1/1 genotypes.
      #   (b) Extreme containment cases (proband much smaller than parent, prop_hijo
      #       at the boundary of umbral_contencion) may still be missed by the main pass.
      #   (c) 1/1 de novo is biologically implausible regardless of size ratio.
      #
      # Rescue uses pmax(sim_maxlen, sim_contained) which specifically recovers the
      # "large parent fully contains small proband" scenario that the main pass may miss
      # when prop_hijo falls just below umbral_contencion.
      es_candidato_rescue <- !her_pa & !her_ma
      
      if (any(es_candidato_rescue)) {
        n_hom <- sum(es_candidato_rescue & !is.na(gt_hijo_rescue) & gt_hijo_rescue == "1/1")
        n_het <- sum(es_candidato_rescue & (!is.na(gt_hijo_rescue) & gt_hijo_rescue != "1/1" |
                                              is.na(gt_hijo_rescue)))
        cat(paste0("      [Rescue] ", sum(es_candidato_rescue),
                   " unmatched variant(s) (1/1: ", n_hom, " | other/unknown: ", n_het,
                   ") — retrying with relaxed thresholds...\n"))
        
        idx_rescue      <- which(es_candidato_rescue)
        dt_hijo_rescue  <- dt_hijo[idx_rescue, ]
        
        resc_pa <- calc_herencia_rescue(dt_hijo_rescue, dt_pa_n)
        resc_ma <- calc_herencia_rescue(dt_hijo_rescue, dt_ma_n)
        
        # Merge rescue results back only where the main pass found nothing
        for (fi in seq_along(idx_rescue)) {
          gi <- idx_rescue[fi]
          if (resc_pa$heredado[fi] && !her_pa[gi]) {
            her_pa[gi]  <- TRUE
            conf_pa[gi] <- resc_pa$confianza[fi]
            ids_pa[gi]  <- resc_pa$ids_match[fi]
            gt_pa[gi]   <- resc_pa$genotype[fi]
          }
          if (resc_ma$heredado[fi] && !her_ma[gi]) {
            her_ma[gi]  <- TRUE
            conf_ma[gi] <- resc_ma$confianza[fi]
            ids_ma[gi]  <- resc_ma$ids_match[fi]
            gt_ma[gi]   <- resc_ma$genotype[fi]
          }
        }
        cat(paste0("      [Rescue] Rescued: PA=", sum(resc_pa$heredado),
                   " MA=", sum(resc_ma$heredado), "\n"))
      }
      
      # TRUE if the parent file has real data (not empty)
      tiene_padre <- nrow(dt_padre_filtrado) > 0
      tiene_madre <- nrow(dt_madre_filtrado) > 0
      
      # Label: origin + confidence (if Probable) + parent genotype in brackets
      fmt_gt  <- function(gt) ifelse(!is.na(gt) & nzchar(gt), paste0(" [", gt, "]"), "")
      fmt_duo <- function(gp, gm) paste0(
        " [", ifelse(!is.na(gp) & nzchar(gp), gp, "?"),
        " | ", ifelse(!is.na(gm) & nzchar(gm), gm, "?"), "]"
      )
      
      # If neither parent carries the variant, the final classification
      # depends on whether both were available:
      #   - Both available and absent -> De novo (confirmed)
      #   - One missing              -> Unknown (cannot be confirmed)
      etiqueta_sin_herencia <- ifelse(tiene_padre & tiene_madre, "De novo", "Unknown")
      
      tipo_her <- ifelse(
        her_pa & her_ma & !is.na(conf_pa) & conf_pa == "Strong" & !is.na(conf_ma) & conf_ma == "Strong",
        paste0("Combined", fmt_duo(gt_pa, gt_ma)),
        ifelse(
          her_pa & her_ma,
          paste0("Combined (Probable)", fmt_duo(gt_pa, gt_ma)),
          ifelse(
            her_pa & !is.na(conf_pa) & conf_pa == "Strong",
            paste0("Paternal", fmt_gt(gt_pa)),
            ifelse(
              her_pa,
              paste0("Paternal (Probable)", fmt_gt(gt_pa)),
              ifelse(
                her_ma & !is.na(conf_ma) & conf_ma == "Strong",
                paste0("Maternal", fmt_gt(gt_ma)),
                ifelse(
                  her_ma,
                  paste0("Maternal (Probable)", fmt_gt(gt_ma)),
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
      
      # --- De_novo_Quality flag ---
      # Adds a diagnostic column to highlight De novo calls that deserve manual review.
      # Does NOT change Tipo_Herencia; it is purely informational.
      #
      #   "CHECK: hom proband"     — proband is 1/1 yet still De novo after rescue pass.
      #                              This is very unusual; likely a caller artefact or a
      #                              genuinely rare homozygous de novo (e.g. UPD + de novo).
      #   "INCOMPLETE: missing parent" — De novo call but at least one parent file is absent,
      #                              so "de novo" cannot be confirmed: it may simply be
      #                              inherited from the missing parent.
      #   NA                       — no concern: either not De novo, or De novo with both
      #                              parents present and diploid proband genotype.
      {
        gt_full_flag <- if ("Genotipo_Hijo" %in% colnames(dt_fase2)) {
          dt_fase2$Genotipo_Hijo[filas_her_idx]
        } else {
          rep(NA_character_, length(tipo_her))
        }
        
        denovo_flag <- rep(NA_character_, length(tipo_her))
        
        # Case A: De novo but proband is homozygous — very suspicious even after rescue
        denovo_flag[
          tipo_her == "De novo" &
            !is.na(gt_full_flag) & gt_full_flag == "1/1"
        ] <- "CHECK: hom proband"
        
        # Case B: De novo but at least one parent file is missing — cannot be confirmed
        denovo_flag[
          tipo_her == "De novo" & (!tiene_padre | !tiene_madre) &
            is.na(denovo_flag)   # don't overwrite Case A
        ] <- "INCOMPLETE: missing parent"
        
        dt_fase2[filas_her_idx, De_novo_Quality := denovo_flag]
      }
      
      # Propagate from full → split rows by AnnotSV_ID safely
      if(tiene_modo && "AnnotSV_ID" %in% colnames(dt_fase2)) {
        dt_her_map <- dt_fase2[Annotation_mode == "full", .(AnnotSV_ID, Tipo_Herencia, IDs_Coincidencia_PA, IDs_Coincidencia_MA, De_novo_Quality)]
        dt_her_map <- unique(dt_her_map, by = "AnnotSV_ID")
        
        dt_fase2[Annotation_mode == "split", c("Tipo_Herencia", "IDs_Coincidencia_PA", "IDs_Coincidencia_MA", "De_novo_Quality") := {
          idx <- match(AnnotSV_ID, dt_her_map$AnnotSV_ID)
          list(
            dt_her_map$Tipo_Herencia[idx],
            dt_her_map$IDs_Coincidencia_PA[idx],
            dt_her_map$IDs_Coincidencia_MA[idx],
            dt_her_map$De_novo_Quality[idx]
          )
        }]
      }
    }
  } else {
    cat("      \u26a0 No parent files available; inheritance not calculated\n")
  }
  
  # --- PROBAND GENOTYPE (split-row propagation) ---
  # GT was already extracted into Genotipo_Hijo before the inheritance block.
  # Here we only need to propagate it from "full" → "split" rows by AnnotSV_ID.
  # If, for any reason, the column is still absent (e.g. no Samples_ID column),
  # we create it here as NA so downstream code never errors.
  if (!"Genotipo_Hijo" %in% colnames(dt_fase2)) {
    dt_fase2[, Genotipo_Hijo := NA_character_]
  }
  # Propagate full → split rows by AnnotSV_ID
  if ("Annotation_mode" %in% colnames(dt_fase2) && "AnnotSV_ID" %in% colnames(dt_fase2)) {
    map_gt <- dt_fase2[Annotation_mode == "full", .(AnnotSV_ID, Genotipo_Hijo)]
    map_gt <- unique(map_gt, by = "AnnotSV_ID")
    dt_fase2[Annotation_mode == "split", Genotipo_Hijo := {
      idx <- match(AnnotSV_ID, map_gt$AnnotSV_ID)
      map_gt$Genotipo_Hijo[idx]
    }]
  }
  
  # --- SORTING ---
  if("Annotation_mode" %in% colnames(dt_fase2)) {
    dt_fase2[, is_full := fifelse(Annotation_mode == "full", 0L, 1L)]
  } else {
    dt_fase2[, is_full := fifelse(grepl("[;,]", Gene_name), 0L, 1L)]
  }
  
  dt_final <- dt_fase2[order(AnnotSV_ID, is_full, -ranking_numeric)]
  dt_final[, is_full := NULL]
  
  # --- DEDUPLICATION OF SIMILAR REGIONS ---
  # Collapses variants overlapping with sim >= umbral_herencia on the same chr and type.
  # If any has inheritance other than "De novo", that one is kept; if all are
  # "De novo", the one with the highest ranking_numeric is kept.
  # Enable/disable with CONFIG$deduplicar_regiones.
  if(isTRUE(CONFIG$deduplicar_regiones)) {
    n_antes_dedup <- nrow(dt_final)
    dt_final <- deduplicar_regiones_similares(dt_final)
    n_despues_dedup <- nrow(dt_final)
    if(n_antes_dedup != n_despues_dedup) {
      cat(paste0("      [Dedup] Total: ", n_antes_dedup, " \u2192 ", n_despues_dedup, " variants\n"))
    }
  } else {
    cat("      [Dedup] Collapsing disabled (CONFIG$deduplicar_regiones = FALSE)\n")
  }
  
  # --- RESTORE ORIGINAL COLUMNS WITH REAL VALUES ---
  # For each column from the original TSV absent in dt_final, recover the real
  # values from the snapshot dt_annotsv_src using .src_row_ as the row key.
  # This guarantees not only that the column name appears but that the correct value is present.
  cols_aux_internas <- c("ranking_numeric", "idx_original", "is_full", ".src_row_")
  cols_esperadas    <- setdiff(cols_originales, cols_aux_internas)
  cols_perdidas     <- setdiff(cols_esperadas, names(dt_final))
  if (length(cols_perdidas) > 0) {
    cat(paste0("      \u26a0 Recovering ", length(cols_perdidas), " column(s) from snapshot: ",
               paste(cols_perdidas, collapse = ", "), "\n"))
    pos <- match(dt_final$.src_row_, dt_annotsv_src$.src_row_)
    for (col in cols_perdidas) {
      dt_final[, (col) := dt_annotsv_src[[col]][pos]]
    }
  }
  
  cols_inicio  <- c("Tipo_Herencia", "De_novo_Quality", "Genotipo_Hijo", "Tipo_Rango", "En_Panel_Genes", "AnnotSV_ranking_score")
  cols_fin     <- c("IDs_Coincidencia_PA", "IDs_Coincidencia_MA")
  # Remove internal auxiliary columns that should not appear in the Excel output
  cols_aux     <- c("idx_original", "ranking_numeric", ".src_row_")
  for (col_aux in cols_aux) {
    if (col_aux %in% names(dt_final)) dt_final[, (col_aux) := NULL]
  }
  cols_inicio_ex <- intersect(cols_inicio, names(dt_final))
  cols_fin_ex    <- intersect(cols_fin,    names(dt_final))
  cols_medio     <- setdiff(names(dt_final), c(cols_inicio, cols_fin))
  
  setcolorder(dt_final, c(cols_inicio_ex, cols_medio, cols_fin_ex))
  
  cat(paste0("      \u2713 ", nrow(dt_final), " final variants\n"))
  return(dt_final)
}

# =============================================================================
# 5. OPTIMISED FUNCTION TO SAVE EXCEL
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
                        ifelse(df[[j]], "Yes", "No"))
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
  # Resolve name collisions: AnnotSV columns such as FORMAT[GT], FORMAT[FT], etc.
  # become names like FORMAT_GT_ after gsub. If two different columns produce
  # the same name, writeData() silently discards duplicates.
  # Suffix _2, _3... is added to guarantee uniqueness.
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
    cat(paste0("      [Excel] Sheet '", nombre_hoja, "': ", ncol(df), " columns, ", nrow(df), " rows\n"))
    addWorksheet(wb, nombre_hoja)
    writeData(wb, nombre_hoja, df)
    
    col_her <- which(names(df) == "Tipo_Herencia")
    col_ran <- which(names(df) == "Tipo_Rango")
    col_sfa <- which(names(df) == "En_Panel_Genes")
    col_sco <- which(names(df) == "AnnotSV_ranking_score")
    col_dnq <- which(names(df) == "De_novo_Quality")
    
    if(length(col_her) > 0) {
      rows_denovo   <- which(!is.na(df$Tipo_Herencia) & df$Tipo_Herencia %in% c("De novo", "Unknown"))
      rows_paterna  <- which(!is.na(df$Tipo_Herencia) & grepl("^Paternal",  df$Tipo_Herencia))
      rows_materna  <- which(!is.na(df$Tipo_Herencia) & grepl("^Maternal",  df$Tipo_Herencia))
      rows_conjunta <- which(!is.na(df$Tipo_Herencia) & grepl("^Combined",  df$Tipo_Herencia))
      
      if(length(rows_denovo) > 0)   addStyle(wb, nombre_hoja, st_denovo,  rows = rows_denovo + 1,   cols = col_her)
      if(length(rows_paterna) > 0)  addStyle(wb, nombre_hoja, st_paterna, rows = rows_paterna + 1,  cols = col_her)
      if(length(rows_materna) > 0)  addStyle(wb, nombre_hoja, st_materna, rows = rows_materna + 1,  cols = col_her)
      if(length(rows_conjunta) > 0) addStyle(wb, nombre_hoja, st_conj,    rows = rows_conjunta + 1, cols = col_her)
    }
    
    if(length(col_ran) > 0) {
      rows_estricto <- which(!is.na(df$Tipo_Rango) & df$Tipo_Rango == "Strict")
      rows_amplio   <- which(!is.na(df$Tipo_Rango) & df$Tipo_Rango == "Wide")
      rows_fuera    <- which(!is.na(df$Tipo_Rango) & df$Tipo_Rango == "Outside")
      
      if(length(rows_estricto) > 0) addStyle(wb, nombre_hoja, st_estric, rows = rows_estricto + 1, cols = col_ran)
      if(length(rows_amplio) > 0)   addStyle(wb, nombre_hoja, st_amplio, rows = rows_amplio + 1,   cols = col_ran)
      if(length(rows_fuera) > 0)    addStyle(wb, nombre_hoja, st_fuera,  rows = rows_fuera + 1,    cols = col_ran)
    }
    
    if(length(col_sfa) > 0) {
      rows_sfari_yes <- which(!is.na(df[[col_sfa]]) & df[[col_sfa]] == "Yes")
      rows_sfari_no  <- which(!is.na(df[[col_sfa]]) & df[[col_sfa]] == "No")
      
      if(length(rows_sfari_yes) > 0) addStyle(wb, nombre_hoja, st_green, rows = rows_sfari_yes + 1, cols = col_sfa)
      if(length(rows_sfari_no) > 0)  addStyle(wb, nombre_hoja, st_red,   rows = rows_sfari_no + 1,  cols = col_sfa)
    }
    
    if(length(col_sco) > 0) {
      # AnnotSV_ranking_score: negative = benign (green), 0–0.5 = intermediate (orange), >=0.5 = pathogenic (red)
      score_vals <- suppressWarnings(as.numeric(df$AnnotSV_ranking_score))
      rows_verde   <- which(!is.na(score_vals) & score_vals < 0)
      rows_naranja <- which(!is.na(score_vals) & score_vals >= 0   & score_vals < 0.5)
      rows_rojo    <- which(!is.na(score_vals) & score_vals >= 0.5)
      
      if(length(rows_verde) > 0)   addStyle(wb, nombre_hoja, st_green,  rows = rows_verde + 1,   cols = col_sco)
      if(length(rows_naranja) > 0) addStyle(wb, nombre_hoja, st_orange, rows = rows_naranja + 1, cols = col_sco)
      if(length(rows_rojo) > 0)    addStyle(wb, nombre_hoja, st_red,    rows = rows_rojo + 1,    cols = col_sco)
    }
    
    if(length(col_dnq) > 0 && "De_novo_Quality" %in% names(df)) {
      # De_novo_Quality: highlight suspicious De novo calls for manual review
      st_dnq_hom  <- createStyle(fgFill = "#FF6B6B", fontColour = "#FFFFFF", textDecoration = "bold")
      st_dnq_miss <- createStyle(fgFill = "#FFD580", fontColour = "#5A3E00")
      rows_hom  <- which(!is.na(df$De_novo_Quality) & df$De_novo_Quality == "CHECK: hom proband")
      rows_miss <- which(!is.na(df$De_novo_Quality) & df$De_novo_Quality == "INCOMPLETE: missing parent")
      if(length(rows_hom)  > 0) addStyle(wb, nombre_hoja, st_dnq_hom,  rows = rows_hom  + 1, cols = col_dnq)
      if(length(rows_miss) > 0) addStyle(wb, nombre_hoja, st_dnq_miss, rows = rows_miss + 1, cols = col_dnq)
    }
    
    setColWidths(wb, nombre_hoja, cols = 1:ncol(df), widths = "auto")
  }
  
  saveWorkbook(wb, ruta_salida, overwrite = TRUE)
}

# =============================================================================
# 6. MAIN LOOP (BATCH PROCESSING)
# =============================================================================
cat("\n=============================================================================\n")
cat("STARTING BATCH PROCESSING\n")
cat("=============================================================================\n")

sel <- function(tsv, patron) tsv[grepl(patron, basename(tsv), ignore.case = TRUE)][1]

EXCEPCIONES_PROGENITORES <- list(
  # Format: family_id = list(PA = pa_parent_id, MA = ma_parent_id)
  # The code will look for a directory whose name ends in paste0(parent_id, suffix)
  # e.g. for PA = "G04PMP068" it will look for a directory ending in "G04PMP068PA".
  # ⚠ Make sure the parent directory exists with that exact suffix.
  "G04PMP069" = list(PA = "G04PMP068", MA = "G04PMP068"), "HVAPM0015"=list(PA="HVAPM0053", MA="HVAPM0053")
)

todas_las_carpetas_cohorte <- list.dirs(CONFIG$ruta_entrada, full.names = TRUE, recursive = TRUE)

# Función actualizada
buscar_dir_progenitor <- function(id_familia, sufijo, directorios_sub, todos_directorios) {
  
  # 1. Búsqueda de Excepciones (Nivel Global)
  if(id_familia %in% names(EXCEPCIONES_PROGENITORES)) {
    id_prog <- EXCEPCIONES_PROGENITORES[[id_familia]][[sufijo]]
    if(!is.null(id_prog)) {
      # Buscamos en TODOS los directorios de la cohorte, no solo en la subcarpeta actual
      excepcion <- todos_directorios[grepl(paste0(id_prog, sufijo, "$"), basename(todos_directorios), ignore.case = TRUE)][1]
      
      if(!is.na(excepcion)) {
        cat(paste0("  ℹ Parent expception applied fo: ", id_familia, sufijo,
                   ": ", basename(excepcion), "\n"))
        return(excepcion)
      } else {
        cat(paste0("  ⚠ Exception definded for folder:", id_prog, sufijo, " Not found in cohort\n"))
      }
    }
  }
  
  # 2. Búsqueda Exacta (Nivel Local)
  exacto <- directorios_sub[grepl(paste0(id_familia, sufijo, "$"), directorios_sub, ignore.case = TRUE)][1]
  if(!is.na(exacto)) return(exacto)
  
  # 3. Búsqueda Compartida (Nivel Local)
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
    cat(paste0("  ℹ Shared parents detected for: ", id_familia, sufijo, ": ", basename(compartido), "\n"))
    return(compartido)
  }
  return(NA)
}

procesar_familia <- function(id_familia, directorios_sub, ruta_salida_sub) {
  
  dir_hi <- directorios_sub[grepl(paste0(id_familia, "HI$"), directorios_sub, ignore.case = TRUE)][1]
  
  # Añadimos 'todas_las_carpetas_cohorte' a las llamadas
  dir_pa <- buscar_dir_progenitor(id_familia, "PA", directorios_sub, todas_las_carpetas_cohorte)
  dir_ma <- buscar_dir_progenitor(id_familia, "MA", directorios_sub, todas_las_carpetas_cohorte)
  
  # (El resto de la función se queda exactamente igual)
  if(is.na(dir_hi) || !dir.exists(dir_hi)) {
    cat("  ⚠ No HI family, moving on\n")
    return(FALSE)
  }
  # ...
  
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
    nombre_excel  <- paste0(id_familia, "_Complete_Analysis.xlsx")
    ruta_completa <- file.path(ruta_salida_sub, nombre_excel)
    tryCatch({
      guardar_excel_estilizado(resultados_familia, ruta_completa)
      cat(paste("  \u2713 Excel created:", file.path(basename(ruta_salida_sub), nombre_excel), "\n"))
      return(resultados_familia)
    }, error = function(e) {
      cat(paste("  [ERROR saving Excel]:", e$message, "\n"))
      return(FALSE)
    })
  } else {
    cat("  \u26a0 No results for this family\n")
    return(FALSE)
  }
}

subcarpetas <- list.dirs(CONFIG$ruta_entrada, full.names = TRUE, recursive = FALSE)

# TSV path cache: defined BEFORE building tasks and before any call
# to procesar_familia, to avoid references to undefined objects in any
# execution order or parallel worker environment.
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
  
  # FIX: derive family IDs ONLY from HI directories.
  # Previously all directories (HI + PA + MA) were used, which caused
  # families without an HI folder (only PA or MA) to be registered as tasks.
  # Now: a family only exists if it has an HI directory; parents
  # are located afterwards via buscar_dir_progenitor, which already handles both
  # unique and shared parents (nnn-mmm pattern).
  hi_dirs     <- ids_dir[grepl("HI$", ids_dir, ignore.case = TRUE)]
  ids_familia <- unique(sub("HI$", "", hi_dirs, ignore.case = TRUE))
  # Exclude: empty strings (safeguard) and any residual ID with a shared-parent
  # pattern (digits-digits), although in practice these should not appear
  # because shared directories end in PA/MA, not HI.
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
  cat(paste0("  [", tarea$nombre_sub, "] Family: ", id_familia, "\n"))
  
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
cat("\u2713 PROCESSING COMPLETE\n")
cat(paste0("  Subfolders processed: ", n_subcarpetas, "\n"))
cat(paste0("  Families processed:   ", n_familias,    "\n"))
cat(paste0("  Files generated:      ", n_exitosas,    "\n"))
cat("=============================================================================\n")

# =============================================================================
# 7. COLLECTION OF HIGH-IMPACT VARIANTS (Score >= 0.15)
# =============================================================================
cat("\n=============================================================================\n")
cat("GENERATING SUMMARY DOCUMENT (SCORE >= 0.15)\n")
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
  
  # ranking_numeric is an internal auxiliary column removed before output;
  # here it is recalculated from AnnotSV_ranking_score to avoid depending on it.
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

ruta_recopilacion <- file.path(CONFIG$ruta_salida, "High_Impact_Summary")

if (nrow(dt_alto_cnv) > 0 || nrow(dt_alto_sv) > 0) {
  if(!dir.exists(ruta_recopilacion)) dir.create(ruta_recopilacion, recursive = TRUE)
  
  resultados_recopilados <- list()
  if (nrow(dt_alto_cnv) > 0) resultados_recopilados[["CNVs_Priority"]] <- dt_alto_cnv
  if (nrow(dt_alto_sv) > 0)  resultados_recopilados[["SVs_Priority"]]  <- dt_alto_sv
  
  ruta_excel_recop <- file.path(ruta_recopilacion, "Variants_Score_Above_0.15.xlsx")
  tryCatch({
    guardar_excel_estilizado(resultados_recopilados, ruta_excel_recop)
    cat(paste("\u2713 Summary document successfully created at:\n  ", ruta_excel_recop, "\n"))
    cat(paste("  - CNV rows included:", nrow(dt_alto_cnv), "\n"))
    cat(paste("  - SV rows included: ", nrow(dt_alto_sv), "\n"))
  }, error = function(e) {
    cat(paste("\u274c Error saving summary document:", e$message, "\n"))
  })
} else {
  cat("\u26a0 No variants with AnnotSV_ranking_score >= 0.15 found across the cohort.\n")
}

cat("=============================================================================\n")
cat("=============================================================================\n")