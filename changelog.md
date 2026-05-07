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