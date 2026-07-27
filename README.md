# PhyTaxR <img src="man/figures/logo.png" align="right" height="139" alt="PhyTaxR logo" />

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![Version: 0.2.9](https://img.shields.io/badge/version-0.2.9-blue.svg)
[![R-CMD-check](https://github.com/RayonVina/phytaxr/actions/workflows/r.yml/badge.svg)](https://github.com/RayonVina/phytaxr/actions/workflows/r.yml)
![Last commit](https://img.shields.io/github/last-commit/RayonVina/phytaxr?style=flat-square&color=b4befe)
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

The package implements a three-stage pipeline:

1. **Cleaning** — a sequence of lightweight, composable functions that
   normalise a raw taxon string into a clean scientific name plus a
   structured `tax_epithet` column holding qualifiers (size info,
   morphotypes, authorship, uncertainty flags, etc.).

2. **Automatic resolution** — step-by-step functions that query WoRMS and GBIF to
   return AphiaIDs, accepted names, synonymy notes, and full taxonomic
   classification from kingdom to forma.

3. **Semi-automatic fuzzy resolution** — an interactive batch review tool
   (`process_fuzzy_batch()`) for names that could not be resolved
   automatically, using fuzzy matching with human-in-the-loop validation.

All stages operate on plain data frames and are designed to be chained
with the `|>` pipe. Resolution functions also accept a bare character
vector directly, and support a custom column name via the `col` argument.

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
  normalize_characters()                 |> # encoding, diacritics, invisible spaces
  process_taxonomic_prefixes()           |> # O./C./F./P. rank prefixes
  process_incertae_entries()             |> # cf., aff., incertae sedis, s.l., etc.
  process_sp_entries()                   |> # sp., spp., sp1, etc.
  process_bracket_entries()              |> # bracket-delimited qualifiers []
  move_size_to_epithet()                 |> # size annotations (µm, mm, ranges)
  process_epithet_entries()              |> # parenthesised epithets
  normalize_infraspecific_ranks()        |> # var., subsp., f., ssp., cv.
  remove_dots()                          |> # stray dots (protects rank dots)
  move_reproductive_structures()         |> # cysts, spores, filaments, …
  move_uncertainty_descriptors()         |> # unidentified, unknown, unclassified, …
  move_morphological_descriptors()       |> # centric, pennate, fusiform, …
  move_with_descriptors()                |> # "with …" phrases
  move_formia_to_epithet()               |> # forma designations
  move_commas_to_epithet()               |> # post-comma strings
  move_authors_to_epithet()              |> # authorships and dates
  remove_sp_tokens()                     |> # residual sp/spp/ssp in epithet
  split_separator_entries()              |> # unfold /, &, +, or entries
  process_generic_taxa()                 |> # vernacular → scientific name
  clean_trailing_hyphens()               |> # trailing/isolated hyphens
  remove_short_interstitial_tokens()     |> # 1–2 char genus abbreviation artefacts
  apply_vernacular_corrections()            # prefix-based name correction from dictionary

# ── Stage 2: Automatic resolution ─────────────────────────────────
# Run all steps in one call (recommended)
df_res <- run_resolution_pipeline(df_clean)

# Or step by step:
df_res <- search_worms_priority(df_clean)     # 1. WoRMS exact
df_res <- search_worms_taxamatch(df_res)       # 2. WoRMS taxamatch
df_res <- search_gbif_strict(df_res)           # 3. GBIF strict
df_res <- resolve_taxonomic_status(df_res)     # 4. Resolve synonyms
df_res <- search_worms_fuzzy_minor(df_res)     # 5. Fuzzy minor corrections
df_res <- get_taxonomy(df_res)                 # 6. Full classification

# ── Stage 3: Semi-automatic fuzzy resolution ───────────────────────
# For names still unresolved after Stage 2, use the interactive batch tool.
# Build vocabularies from already-resolved names (used as candidate pool):
genus_vocab   <- build_genus_vocabulary(df_res)
epithet_vocab <- build_epithet_vocabulary(df_res)

# Launch interactive review — presents fuzzy candidates for each unresolved
# name and prompts the user to accept, skip, or manually assign a match.
# Progress is saved automatically to a checkpoint file so the session can
# be interrupted and resumed without losing work.
df_final <- process_fuzzy_batch(
  df              = df_res,
  genus_vocab     = genus_vocab,
  epithet_vocab   = epithet_vocab,
  taxon_col       = "taxon",
  taxon_clean_col = "taxon_clean",
  batch_size      = 10,
  checkpoint_file = "phytaxr_step3_checkpoint.rds",
  min_similarity  = 0.85,
  max_suggestions = 15,
  edit_max_dist   = 3,
  timeout_sec     = 15
)

# Inspect
dplyr::glimpse(df_final)
```

Resolution functions also work directly on a character vector:

```r
search_worms_priority(c("Chaetoceros debilis", "Thalassiosira weissflogii"))
```

And accept a custom column name via `col`:

```r
# If your taxon column is named "species" instead of "taxon_clean"
search_worms_priority(df, col = "species")
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
| 13 | `process_generic_taxa()` | Vernacular → scientific name (exact match) |
| 14 | `clean_trailing_hyphens()` | Trailing/isolated hyphens |
| 15 | `remove_short_interstitial_tokens()` | 1–2 char genus abbreviation artefacts (e.g. `p`, `p-n`) |
| 16 | `apply_vernacular_corrections()` | Prefix-based name correction from the vernacular dictionary — replaces known patterns at the start of `taxon_clean` with their canonical form, preserving any trailing epithet (e.g. `"pseudo nitzschia delicatissima"` → `"Pseudo-nitzschia delicatissima"`) |

---

## Automatic resolution pipeline

| Step | Function | Method |
|------|----------|--------|
| 1 | `search_worms_priority()` | WoRMS exact match via `wm_name2id()` + `wm_record()` |
| 2 | `search_worms_taxamatch()` | WoRMS fuzzy match via `wm_records_taxamatch()` |
| 3 | `search_gbif_strict()` | GBIF name match API, confidence ≥ 99, cross-validated against WoRMS |
| 4 | `resolve_taxonomic_status()` | Follows `valid_AphiaID` pointer to resolve synonyms to accepted names |
| 5 | `search_worms_fuzzy_minor()` | Levenshtein ≤ 3, similarity ≥ 0.85, genus verified in WoRMS |
| 6 | `get_taxonomy()` | Retrieves full taxonomic hierarchy (kingdom → forma) from WoRMS |

All six functions share the same interface: they accept a data frame (or
character vector) and return it with resolution columns populated for
previously unresolved rows. All steps are also available as a single call
via `run_resolution_pipeline()`.

---

## Semi-automatic fuzzy resolution

`process_fuzzy_batch()` handles names that remain unresolved after the
automatic pipeline — typically rare taxa, misspellings too severe for
automated correction, or names absent from WoRMS/GBIF.

**How it works:**

1. For each unresolved name, candidate matches are generated using
   Levenshtein distance on genus and epithet vocabularies built from
   already-resolved names in the same dataset.
2. Candidates are ranked by similarity score and presented to the user
   in the console, one batch at a time.
3. The user accepts a candidate, skips the entry, or provides a manual
   AphiaID. Accepted matches are written back into the data frame with
   `resolution_method = "fuzzy_manual"`.
4. Progress is checkpointed automatically after each batch. If the
   session is interrupted, the next run resumes from where it stopped —
   no work is lost.

**Key parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `taxon_col` | `"taxon"` | Column with the original verbatim taxon name (display only) |
| `taxon_clean_col` | `"taxon_clean"` | Column with the cleaned taxon name used for WoRMS searches |
| `aphiaid_col` | `"matched_aphiaid"` | Column used to identify unresolved rows (`NA` = unresolved) |
| `flag_col` | `"flag_for_removal"` | Column for removal flags |
| `method_col` | `"resolution_method"` | Column for resolution method |
| `batch_size` | 10 | Names reviewed per interactive round |
| `checkpoint_file` | `"phytaxr_step3_checkpoint.rds"` | Path for auto-save |
| `min_similarity` | 0.85 | Minimum similarity score to show a candidate |
| `max_suggestions` | 15 | Maximum candidates shown per name |
| `longshot_threshold` | 0.30 | Lower threshold for speculative suggestions |
| `edit_max_dist` | 3 | Maximum Levenshtein distance for candidates |
| `timeout_sec` | 15 | Seconds before skipping a name automatically |

**Resuming an interrupted session:**

```r
# On a subsequent run, the checkpoint is detected and loaded automatically.
# The vocabularies should be rebuilt from the current state of the data frame.
genus_vocab   <- build_genus_vocabulary(df_res)
epithet_vocab <- build_epithet_vocabulary(df_res)

df_final <- process_fuzzy_batch(
  df            = df_res,
  genus_vocab   = genus_vocab,
  epithet_vocab = epithet_vocab,
  checkpoint_file = "phytaxr_step3_checkpoint.rds"
)
```

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
| `resolution_method` | Which step resolved the name (`"worms_exact"`, `"worms_taxamatch"`, `"fuzzy_manual"`, etc.) |
| `resolution_notes` | Synonymy notes and nomenclatural changes |
| `kingdom` … `forma` | Full taxonomic hierarchy |

---

## Vernacular dictionary

The package ships an internal dictionary of common phytoplankton vernacular
names (e.g. `"diatom"` → `"Bacillariophyceae"`, `"dinoflagellate"` →
`"Dinoflagellata"`). Two functions use it at different points of the
cleaning pipeline:

- **`process_generic_taxa()`** (step 13) — exact match against the full
  `taxon_clean` string. Handles standalone vernacular names
  (e.g. `"diatoms"` → `"Bacillariophyceae"`).
- **`apply_vernacular_corrections()`** (step 16) — prefix-based search:
  iterates over every `taxon_clean` value and, if a dictionary key is
  found at the *start* of the string (case-insensitive, longest key
  first), replaces it with its canonical form while preserving any
  trailing epithet. Handles genus-level typos and normalisation errors
  introduced during cleaning (e.g. `"pseudo nitzschia delicatissima"` →
  `"Pseudo-nitzschia delicatissima"`).

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

To make dictionary changes permanent, edit `data-raw/vernacular.R` and
re-run `source("data-raw/vernacular.R")`.

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

If you use **PhyTaxR** in published work, please cite it as:

> Rayón Viña, F. (2026). *PhyTaxR: Phytoplankton Taxonomic Curation Tools*.
> R package version 0.2.9.
> <https://github.com/RayonVina/phytaxr>

---

## License

MIT © Fernando Rayón Viña
