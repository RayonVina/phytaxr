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