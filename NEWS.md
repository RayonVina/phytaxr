# phytaxr 0.0.1

## Bug fixes

* `lookup_taxonomy_info()` now correctly accesses `result$order` instead
  of the non-existent `result$order_name` column, which was silently
  returning `NA` for the `order` rank in every taxonomy look-up
  (`R/resolve.R`).

* `calculate_similarity()` now returns `0` when either input is an empty
  string (`""`), not only when it is `NA`. This prevents silent incorrect
  similarity scores for names that are empty after cleaning
  (`R/utils.R`).

## Refactoring

* Removed duplicate definitions of `calculate_similarity()` and
  `verify_genus_exists()` from `R/strict.R`. Both functions are defined
  (and documented) in `R/utils.R`; the duplicates in `strict.R` were
  triggering *R CMD CHECK* warnings about masked/duplicated definitions
  (`R/strict.R`).

* Added `@name op-null-default` to `%||%` in `R/utils.R` so that
  roxygen2 generates a valid `.Rd` filename without pipe characters,
  eliminating the `checkRd` WARNING (`R/utils.R`).

* Converted `LICENSE` to DCF stub format (`YEAR` / `COPYRIGHT HOLDER`)
  as required by R; full MIT text moved to `LICENSE.note`.
