# CNV/SV Analysis Pipeline — FiltroB

A Shiny application for automated analysis, filtering, and visualization of **copy number variants (CNVs)** and **structural variants (SVs)** from trio-based genomic studies.

The pipeline processes [AnnotSV](https://lbgi.fr/AnnotSV/)-annotated output files, applies configurable filtering thresholds, classifies variants by inheritance pattern, and provides an interactive Shiny interface for clinical review.

---

## Features

- **Automated pipeline** (`FiltroB_STEA.R`): filters and annotates CNV/SV calls from trio families (proband + parents), classifies inheritance (de novo, maternal, paternal, joint), and scores variants.
- **Interactive Shiny app** (`app.R`): real-time pipeline execution, variant table with filtering/tagging/classification, genomic ideogram visualization, family comparator, statistical tests, and HPO integration via JAX API.
- **Configurable**: thresholds for similarity, Jaccard overlap, lateral margin, and deduplication are all adjustable from the UI.
- **Persistent annotations**: variant flags, clinical classifications, notes, sex, phenotype, and HPO terms are saved locally between sessions.

---

## Requirements

- **R** ≥ 4.2
- **AnnotSV** output files (`.xlsx`) as input — one per family, with standard AnnotSV columns
- Input folder structure: one subfolder per family, each containing `HI/`, `PA/`, `MA/` subfolders with the proband and parental AnnotSV files

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/cnv-sv-pipeline.git
cd cnv-sv-pipeline
```

### 2. Install R dependencies

Open R or RStudio and run:

```r
source("instalar_dependencias.R")
```

This installs all required CRAN and Bioconductor packages.

### 3. Prepare your reference files

Place your reference files in the `reference/` folder (or configure custom paths in the app):

| File | Description |
|------|-------------|
| `reference/reference_ranges.xlsx` | Reference genomic ranges for overlap classification (see format below) |
| `reference/gene_panel.txt` | Gene panel for variant prioritization (comma-separated gene symbols) |

See [`reference/README_reference_files.md`](reference/README_reference_files.md) for detailed format specifications.

---

## Usage

### Option A — Shiny App (recommended)

```r
shiny::runApp("app.R")
```

Or open `app.R` in RStudio and click **Run App**.

In the **Configuration** tab, set:
- Input folder (with family subfolders)
- Output folder
- Reference ranges file (`.xlsx`)
- Gene panel file (`.txt`)

Then go to **Execution** and run the pipeline.

### Option B — Command line (pipeline only)

Edit the `CONFIG` block at the top of `FiltroB_STEA.R` with your paths, then:

```r
source("FiltroB_STEA.R")
```

---

## Repository Structure

```
cnv-sv-pipeline/
├── app.R                          # Shiny application
├── FiltroB_STEA.R                 # Core analysis pipeline
├── instalar_dependencias.R        # Dependency installer
├── README.md
├── LICENSE
├── .gitignore
└── reference/
    ├── reference_ranges.xlsx      # Example reference ranges
    ├── gene_panel.txt             # Example gene panel
    └── README_reference_files.md  # Format documentation
```

Runtime folders (auto-created, git-ignored):
- `results/` — pipeline output Excel files
- `logs/` — persistent session annotations (flags, notes, classifications)

---

## Reference File Formats

### `reference_ranges.xlsx`

An Excel file where each row defines a named genomic region. Required columns (by position):

| Column | Content |
|--------|---------|
| 1 | Region name (e.g. `22q11.2_deletion`) |
| 2 | Variant type (`DEL` or `DUP`) |
| 3 | Chromosome (e.g. `22` or `chr22`) |
| 4 | *(reserved / ignored)* |
| 5 | Start — wide boundary |
| 6 | Start — strict boundary |
| 7 | End — strict boundary |
| 8 | End — wide boundary |

### `gene_panel.txt`

Plain text file with gene symbols, comma-separated (single or multiple lines):

```
GENE1,GENE2,GENE3,GENE4,...
```

The included example (`reference/gene_panel.txt`) contains SFARI autism-associated genes.

---

## Input Data Format

The pipeline expects **AnnotSV v3** output files in `.xlsx` format, one per individual, placed in:

```
<input_folder>/
  <family_id>/
    HI/   → proband AnnotSV file (one .xlsx)
    PA/   → father AnnotSV file (one .xlsx)
    MA/   → mother AnnotSV file (one .xlsx)
```

---

## Output

For each family, the pipeline produces `<family_id>_Analisis_Completo.xlsx` in the output folder, containing:
- Filtered and ranked variants
- Inheritance classification
- Overlap with reference ranges (Strict / Wide / Outside)
- Gene panel membership flag (`En_Panel_Genes`)
- Similarity scores and Jaccard overlap with parental variants

---

## Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `umbral_strict` | 0.70 | Minimum overlap ratio for strict range classification |
| `umbral_herencia` | 0.70 | Minimum similarity for "Strong" inheritance |
| `umbral_herencia_lax` | 0.60 | Minimum similarity for "Probable" inheritance (rescued by Jaccard) |
| `umbral_jaccard` | 0.50 | Minimum Jaccard index for "Probable" inheritance |
| `margen_lateral` | 3,000,000 | Maximum lateral deviation allowed (bp) |
| `deduplicar` | TRUE | Collapse similar variants within the same family |
| `n_cores` | auto | Number of CPU cores for parallel processing |

---

## License

MIT License — see [LICENSE](LICENSE).

---

## Citation

If you use this pipeline in your work, please cite the repository:

```
STEA CNV/SV Analysis Pipeline. GitHub: https://github.com/your-username/cnv-sv-pipeline
```
