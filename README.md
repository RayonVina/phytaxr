# PhyTaxR

<!-- badges: start -->
![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![Version: 0.0.1](https://img.shields.io/badge/version-0.0.1-blue.svg)
<!-- badges: end -->

**PhyTaxR** provides tools for cleaning, standardising, and resolving
phytoplankton taxon names against external taxonomic databases
([WoRMS](https://www.marinespecies.org/) and [GBIF](https://www.gbif.org/)).

> **Status:** alpha — API may change without notice.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("RayonVina/phytaxr")
```

---

## Overview

The package implements a two-stage pipeline:

1. **Cleaning** — a sequence of lightweight, composable functions that
   normalise a raw taxon string into a clean scientific name plus a
   structured `tax_epithet` column holding qualifiers (size info,
   morphotypes, authorship, uncertainty flags, etc.).

2. **Resolution** — step-by-step functions that query WoRMS and GBIF to
   return AphiaIDs, accepted names, synonymy notes, and full taxonomic
   classification from kingdom to forma.

Both stages operate on plain data frames and are designed to be chained
with the `|>` pipe.

---

## Basic usage

```r
library(phytaxr)
library(dplyr)
library(tibble)

# Build a data frame from raw names
df <- tibble(taxon = c(
  "Chaetoceros debilis",
  "Thalassiosira weissflogii",
  "C. diatom sp.",
  "Gymnodinium cf. catenatum"
))

# ── Stage 1: Cleaning ──────────────────────────────────────────────
df_clean <- df |>
  normalize_characters()           |>
  process_taxonomic_prefixes()     |>
  process_incertae_entries()       |>
  process_sp_entries()             |>
  process_bracket_entries()        |>
  move_size_to_epithet()           |>
  process_epithet_entries()        |>
  normalize_infraspecific_ranks()  |>
  remove_dots()                    |>
  move_reproductive_structures()   |>
  move_uncertainty_descriptors()   |>
  move_morphological_descriptors() |>
  move_with_descriptors()          |>
  move_formia_to_epithet()         |>
  move_commas_to_epithet()         |>
  move_authors_to_epithet()        |>
  remove_sp_tokens()               |>
  split_separator_entries()        |>
  process_generic_taxa()           |>
  clean_trailing_hyphens()

# ── Stage 2: Resolution ────────────────────────────────────────────
# 2a. WoRMS exact match
df_res <- search_worms_priority(df_clean)

# 2b. WoRMS taxamatch (for names not resolved above)
df_res <- search_worms_taxamatch(df_res)

# 2c. GBIF strict match (cross-validated against WoRMS)
df_res <- search_gbif_strict(df_res)

# 2d. Fuzzy minor correction (Levenshtein ≤ 3, similarity ≥ 0.85)
df_res <- search_fuzzy_minor(df_res)

# Inspect
dplyr::glimpse(df_res)
```

---

## Cleaning pipeline

Each function takes a data frame with a `taxon_clean` column and returns
it modified. They must be called in order.

| Step | Function | What it does |
|------|----------|--------------|
| 1 | `normalize_characters()` | Encoding, diacritics, invisible spaces |
| 2 | `process_taxonomic_prefixes()` | O./C./F./P. rank prefixes |
| 3 | `process_incertae_entries()` | cf., aff., *incertae sedis*, s.l., etc. |
| 4 | `process_sp_entries()` | sp., spp., sp1, etc. |
| 5 | `process_bracket_entries()` | Bracket-delimited qualifiers `[]` |
| 6 | `move_size_to_epithet()` | Size annotations (µm, mm, ranges) |
| 7 | `process_epithet_entries()` | Parenthesised epithets |
| 8 | `normalize_infraspecific_ranks()` | var., subsp., f., ssp., cv. |
| 9 | `remove_dots()` | Stray dots (protects rank dots) |
| 10a | `move_reproductive_structures()` | cysts, spores, filaments, … |
| 10b | `move_uncertainty_descriptors()` | unidentified, unknown, unclassified, … |
| 10c | `move_morphological_descriptors()` | centric, pennate, fusiform, … |
| 10d | `move_with_descriptors()` | "with …" phrases |
| 10e | `move_formia_to_epithet()` | *forma* designations |
| 10f | `move_commas_to_epithet()` | Post-comma strings |
| 10g | `move_authors_to_epithet()` | Authorships and dates |
| 11 | `remove_sp_tokens()` | Residual sp/spp/ssp in epithet |
| 12 | `split_separator_entries()` | Unfold `/`, `&`, `+`, `or` entries |
| 13 | `process_generic_taxa()` | Vernacular → scientific name |
| 14 | `clean_trailing_hyphens()` | Trailing/isolated hyphens |

---

## Resolution pipeline

| Step | Function | Method |
|------|----------|--------|
| 1 | `search_worms_priority()` | WoRMS exact match via `wm_name2id()` + `wm_record()` |
| 2 | `search_worms_taxamatch()` | WoRMS fuzzy match via `wm_records_taxamatch()` |
| 3 | `search_gbif_strict()` | GBIF name match API, confidence ≥ 99, cross-validated against WoRMS |
| 4 | `search_fuzzy_minor()` | Levenshtein ≤ 3, similarity ≥ 0.85, genus verified in WoRMS |

All four functions share the same interface: they accept a data frame and
return it with resolution columns populated for previously unresolved rows.

---

## Output columns

| Column | Description |
|--------|-------------|
| `taxon` | Raw input string (unchanged) |
| `taxon_clean` | Cleaned scientific name |
| `tax_epithet` | Structured qualifiers extracted during cleaning |
| `uncertain` | Logical flag: name could not be fully normalised |
| `matched_name` | Name returned by the resolver |
| `matched_aphiaid` | AphiaID of the matched name |
| `accepted_name` | Currently accepted name |
| `accepted_aphiaid` | AphiaID of the accepted name |
| `taxonomic_status` | `"accepted"`, `"synonym"`, etc. |
| `resolution_method` | Which step resolved the name |
| `resolution_notes` | Synonymy notes and nomenclatural changes |
| `kingdom` … `forma` | Full taxonomic hierarchy |

---

## Vernacular dictionary

The package ships an internal dictionary of common phytoplankton vernacular
names (e.g. `"diatom"` → `"Bacillariophyceae"`, `"dinoflagellate"` →
`"Dinoflagellata"`). `process_generic_taxa()` uses it automatically.

```r
# List all entries
vernacular_list()

# Look up a name (case-insensitive)
vernacular_lookup("diatom")

# Search by pattern
vernacular_search("diatom")

# Add / update / remove (session only — lost on restart)
vernacular_add("green alga", "Chlorophyta")
vernacular_update("diatom", "Bacillariophyta")
vernacular_remove("green alga")
```

---

## Dependencies

| Package | Role |
|---------|------|
| [`worrms`](https://cran.r-project.org/package=worrms) | WoRMS REST API |
| [`httr`](https://httr.r-lib.org/) | GBIF HTTP requests |
| [`stringr`](https://stringr.tidyverse.org/) / [`stringi`](https://stringi.gagolewski.com/) | String processing |
| [`dplyr`](https://dplyr.tidyverse.org/) / [`tidyr`](https://tidyr.tidyverse.org/) | Data manipulation |
| [`purrr`](https://purrr.tidyverse.org/) / [`furrr`](https://furrr.futureverse.org/) | Iteration and parallelisation |
| [`stringdist`](https://cran.r-project.org/package=stringdist) | Fuzzy string matching |
| [`R.utils`](https://cran.r-project.org/package=R.utils) | Timeout handling |
| [`progress`](https://cran.r-project.org/package=progress) | Progress bars |

---

## Citation

If you use **phytaxr** in published work, please cite it as:

> Rayón Viña, F. (2026). *phytaxr: Phytoplankton Taxonomic Curation Tools*.
> R package version 0.0.1.
> <https://github.com/RayonVina/phytaxr>

---

## License

MIT © Fernando Rayón Viña
