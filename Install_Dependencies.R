# instalar_dependencias.R
# Run this script ONCE before launching the app
# ─────────────────────────────────────────────────────────────────────────────
# This script installs all R packages required by the CNV/SV Analysis Pipeline.
# Bioconductor packages (GenomicRanges, IRanges) require BiocManager.
# Optional packages are listed at the bottom.
# ─────────────────────────────────────────────────────────────────────────────

# ── CRAN packages ─────────────────────────────────────────────────────────────
pkgs_cran <- c(
  # Core Shiny + UI
  "shiny",          # Web application framework
  "bslib",          # Bootstrap 5 themes
  "bsicons",        # Bootstrap icons
  "shinyjs",        # JavaScript helpers (enable/disable/runjs)
  "shinyFiles",     # Native file/folder picker

  # Tables
  "DT",             # Interactive DataTables

  # Interactive plots
  "plotly",         # Interactive charts

  # Static plots
  "ggplot2",        # Grammar of graphics
  "gridExtra",      # Arrange multiple ggplot panels
  "grid",           # Low-level graphics (base R, listed for clarity)

  # Async pipeline execution
  "processx",       # Subprocesses with real-time stdout/stderr

  # Excel I/O
  "openxlsx",       # Read/write .xlsx with styles
  "readxl",         # Read .xlsx (fast, read-only)

  # Data wrangling
  "dplyr",
  "tidyr",
  "stringr",
  "scales",
  "data.table",     # High-performance data operations in pipeline

  # Web / API
  "httr",           # HTTP requests (HPO API, external databases)
  "jsonlite",       # JSON parsing

  # Parallel processing
  "parallel"        # Built-in, but listed for documentation purposes
)

# ── Bioconductor packages ──────────────────────────────────────────────────────
pkgs_bioc <- c(
  "GenomicRanges",  # Genomic interval arithmetic
  "IRanges"         # Underlying ranges infrastructure
)

# ── Optional packages (uncomment to install) ──────────────────────────────────
pkgs_optional <- c(
  "png",        # Read PNG images (used in some plot exports)
  "webshot2"    # PDF/PNG export of web-based plots (requires Chrome)
)

# ── Installer function ─────────────────────────────────────────────────────────
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  } else {
    message("Already installed: ", pkg)
  }
}

# ── Install ────────────────────────────────────────────────────────────────────
message("\n── Installing CRAN packages ──────────────────────────────────")
invisible(lapply(pkgs_cran, install_if_missing, bioc = FALSE))

message("\n── Installing Bioconductor packages ──────────────────────────")
invisible(lapply(pkgs_bioc, install_if_missing, bioc = TRUE))

# Uncomment the next line to install optional packages:
# invisible(lapply(pkgs_optional, install_if_missing, bioc = FALSE))

message("\n✅ All dependencies installed.")
message("   Launch the app with: shiny::runApp('app.R')")
message("   Or open app.R in RStudio and click 'Run App'.")
