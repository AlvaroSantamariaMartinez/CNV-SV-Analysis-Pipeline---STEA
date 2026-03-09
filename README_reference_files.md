# Reference Files — Format Documentation

This folder contains the reference files required by the CNV/SV Analysis Pipeline.
Both files must be provided by the user before running the pipeline.

---

## `reference_ranges.xlsx` — Reference genomic ranges

An Excel (`.xlsx`) file where each row defines a named genomic region used to
classify variants as overlapping a known region (Strict / Wide) or falling
outside all known regions.

### Column structure (by position, headers are ignored)

| Column | Field | Type | Example |
|--------|-------|------|---------|
| 1 | Region name | string | `22q11.2_deletion` |
| 2 | Variant type | string | `DEL` or `DUP` |
| 3 | Chromosome | string | `22`, `X`, `chr22` |
| 4 | *(reserved, not used)* | — | — |
| 5 | Wide start (bp) | integer | `18648854` |
| 6 | Strict start (bp) | integer | `19017834` |
| 7 | Strict end (bp) | integer | `21800471` |
| 8 | Wide end (bp) | integer | `22113430` |

- **Strict boundaries**: the core region of the syndrome/locus.
- **Wide boundaries**: extended flanking region allowing positional tolerance.
- Rows with missing start/end coordinates are automatically skipped.
- The header row (if present) is filtered out automatically based on non-numeric coordinates.

### Example rows

| Name | Type | Chr | (skip) | Wide_start | Strict_start | Strict_end | Wide_end |
|------|------|-----|--------|-----------|-------------|-----------|---------|
| 22q11.2_DEL | DEL | 22 | | 18648854 | 19017834 | 21800471 | 22113430 |
| 15q11q13_DUP | DUP | 15 | | 22765628 | 23083774 | 28557186 | 29158407 |

---

## `gene_panel.txt` — Gene panel

A plain text file containing gene symbols for variant prioritization.
Variants overlapping any gene in this panel are flagged in the output column `En_Panel_Genes`.

### Format

- Gene symbols separated by commas (`,`)
- Can be on a single line or multiple lines
- Gene symbols are matched case-insensitively
- Blank entries are ignored

### Example

```
SHANK3,NRXN1,CNTNAP2,SYNGAP1,DYRK1A,ADNP,CHD8,
PTEN,TSC1,TSC2,MECP2,FMR1,ANKRD11,ARID1B
```

The included example file contains **SFARI Gene** autism-associated gene symbols
(https://gene.sfari.org). You can replace it with any gene panel relevant to
your study — for example:

- A disease-specific gene panel (e.g., intellectual disability, epilepsy)
- A custom curated list from your laboratory
- An export from a clinical gene panel database (PanelApp, GenCC, etc.)

---

## Providing your own files

1. Place your reference ranges file in this folder as `reference_ranges.xlsx`
   (or point to a custom path in the app's **Configuration** tab).
2. Place your gene panel file in this folder as `gene_panel.txt`
   (or point to a custom path in the app's **Configuration** tab).

Both paths can be changed at any time from the Shiny app UI without restarting.
