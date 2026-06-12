# phytaxr 0.2.8

## Bug fixes

- `process_fuzzy_batch()`: fixed crash (`missing value where TRUE/FALSE needed`)
  when `taxon_clean_col` pointed to a column renamed internally by `resolve_col()`.
  After `resolve_col()` executes, the internal variable `taxon_clean_col` is now
  updated to `"taxon_clean"` so that subsequent row-level lookups find the column.

- `search_worms_fuzzy_suggestions()`: added `trimws(as.character(...))` guard and
  `!is.na()` check on the normalised name to prevent the cascade from crashing on
  `NA` values in `taxon_clean`.

- `prompt_fuzzy_suggestions()`: added `!is.na(original_taxon)` guard before the
  string comparison `original_taxon != taxon_clean`, preventing a crash when the
  original taxon field is `NA`.

## New features

- `process_fuzzy_batch()`: the five key column names are now configurable via
  explicit parameters (`taxon_col`, `taxon_clean_col`, `aphiaid_col`, `flag_col`,
  `method_col`), all with sensible defaults. This allows the function to be run
  on data frames that use different column naming conventions without renaming
  columns beforehand. The old `col` parameter has been removed.

# phytaxr 0.2.7

### New features

* `apply_vernacular_corrections()`: new function that applies the internal
  vernacular dictionary as a pattern-replacement step over `taxon_clean`.
  Unlike the exact-match lookup used in resolution, this function searches
  for dictionary keys as prefixes within the cleaned name (case-insensitive,
  longest-key-first to avoid partial matches) and replaces them with their
  canonical form, preserving any trailing epithet. Intended to be called as
  the final step of the cleaning pipeline (Step 1), so that `taxon_clean`
  is already well-formed before automatic resolution begins.
  Example: `"pseudo nitzschia delicatissima"` →
  `"Pseudo-nitzschia delicatissima"`.

* Vernacular dictionary expanded with *Pseudo-nitzschia* entries: canonical
  genus name (`"pseudo nitzschia"`, `"pseudonitzschia"`) and 10 common
  typographic variants (missing letters, transposed characters).

* Internal housekeeping: `ensure_resolution_schema()`, `lookup_taxon_info()`,
  `lookup_taxonomy_info()`, `fetch_worms_full_record()`,
  `prompt_fuzzy_suggestions()`, `prompt_manual_aphiaid()`,
  `prompt_resolution_notes()`, and `search_worms_fuzzy_suggestions()` are
  now correctly marked as internal (`@keywords internal`) and removed from
  the public namespace.

# phytaxr 0.2.5

## New features

- `remove_short_interstitial_tokens()`: extended to handle infraspecific
  rank markers in interstitial position (Pattern A). Rank tokens such as
  `var.`, `subsp.`, `f.`, `subf.`, `cv.`, `morph.` and all full ICN/ICZN
  variants are now detected, removed from `taxon_clean`, and preserved in
  `tax_epithet`. The existing short genus-initial abbreviation logic
  (Pattern B) is unchanged but now only fires when Pattern A does not
  match, avoiding double-processing.
  Fixes: `thalassiosira var. expecta` → `thalassiosira expecta`,
  `nitzschia subsp. gracilis` → `nitzschia gracilis`,
  `navicula f. angularis` → `navicula angularis`.
- Fuzzy cascade level L4b: replaced the `nchar(token) <= 2` length check
  with a new `is_genus_initial()` helper that requires the token to be a
  prefix of the genus name. This eliminates false positives from 3-character
  tokens such as `var.` and improves detection reliability (77/77 cases
  confirmed in the dataset).

## Bug fixes

- `process_fuzzy_batch()`: entries with `flag_for_removal = TRUE` were
  incorrectly included in `unresolved_idx` because the filter only excluded
  rows with a resolved `matched_aphiaid`. Fixed by adding
  `(is.na(flag_for_removal) | flag_for_removal != TRUE)` to the filter,
  preventing flagged entries from reappearing on every checkpoint resume.
- `globals.R`: added `has_rank_token` and `rank_token` to
  `globalVariables()` to suppress R CMD check NOTEs introduced by the new
  Pattern A in `remove_short_interstitial_tokens()`. The `.data$` pronoun
  was insufficient in R 4.6.0. Also removed the `rlang` import added in a
  previous attempt.

## Infrastructure

- `NEWS.md`: renamed from `changelog.md` to comply with the CRAN standard
  for top-level changelog files.

---

# phytaxr 0.2.0

## New features

- New `remove_short_interstitial_tokens()` cleaning function (last step
  of the pipeline): removes 1–2 character tokens that remain stranded
  between genus and specific epithet after `sp`/`cf` removal (e.g.
  `pseudo-nitzschia p seriata` → `pseudo-nitzschia seriata`). Affected
  entries are automatically marked `uncertain = TRUE`.

## Bug fixes

- `process_fuzzy_batch()`: vocabularies (`genus_vocab`, `epithet_vocab`)
  are now saved to checkpoint in all branches (`flag_for_removal`,
  `accept`, `expert_review`). Previously only `expert_review` saved
  them, losing vocabularies on session resume if the environment had
  been cleared.
- `process_fuzzy_batch()`: added WoRMS fuzzy raw results as numbered
  fallback options in the interactive prompt when the cascade produces
  no suggestions, mirroring the original `master_taxonomy.R` behaviour.
- Fuzzy cascade levels L1b and L1d no longer require `n_words >= 2`,
  enabling genus-only edit-distance correction and WoRMS fuzzy genus
  lookup for single-token taxa (e.g. `niztschia`).
- Fuzzy cascade level L1d: fixed candidate construction to avoid
  trailing whitespace when `n_words == 1`.
- New fuzzy cascade level L4b: detects single/double-letter abbreviated
  epithets (e.g. `dinophysis d rotundata`, `fragilariopsis f. oceanica`)
  and searches genus + suffix directly instead of treating the
  abbreviation as the epithet.
- Fixed `epithet_jw` ranking to use the suffix as reference epithet
  when `words[2]` is an abbreviation, preventing correct matches from
  being buried in the longshot list.
- `ensure_resolution_schema()`: added `aphiaid` and `marginalia`
  columns; function is now exported.
- Closure scoping fix in `fuzzy.R`: `epithet_ref` is now captured
  explicitly in the `map_dbl()` closure to avoid silent `NA` in
  byte-compiled packages.

## Documentation

- `globals.R`: added `has_short_token` to `globalVariables()`.

# phytaxr 0.1.0

## New features

- All taxonomy search functions (`search_worms_priority()`,
  `search_worms_taxamatch()`, `search_gbif_strict()`,
  `resolve_taxonomic_status()`, `search_worms_fuzzy_minor()`,
  `get_taxonomy()`, `run_resolution_pipeline()`) now accept a bare
  character vector in addition to a data frame, and support a custom
  column name via the `col` argument.
- New internal helpers: `names_to_df()`, `resolve_col()`,
  `format_query_result()`.
- `process_fuzzy_batch()` gains a `col` argument with the same
  semantics.

## Bug fixes

- `resolve_col()` now validates the column name and raises an
  informative error; `process_fuzzy_batch()` delegates column
  validation to `resolve_col()` instead of hardcoding `taxon_clean`.

# phytaxr 0.0.1

## Bug fixes

- `lookup_taxonomy_info()`: fixed silent `NA` in `order` field
  (`result$order_name` → `result$order`).
- `calculate_similarity()`: now returns `0` for empty strings in
  addition to `NA`.
- Added `ensure_resolution_schema()` to initialise missing resolution
  columns; applied at the entry point of all six pipeline functions and
  `process_fuzzy_batch()`.
- Removed duplicate definitions of `calculate_similarity()` and
  `verify_genus_exists()` from `strict.R`.

## Documentation

- Regenerated `man/` and `NAMESPACE`; renamed `run_strict_pipeline.Rd`
  to `run_resolution_pipeline.Rd`.
- Improved docs for `lookup_taxonomy_info()`, `calculate_similarity()`,
  `verify_genus_exists()`, and `%||%`.

## Licence

- Corrected DCF format in `LICENSE`; regenerated MIT licence with
  `usethis::use_mit_license()`.