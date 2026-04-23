# test-strict.R
# Tests for the strict resolution pipeline functions (strict.R).
# API calls are avoided; tests cover argument validation and schema only.

# helpers ----------------------------------------------------------------

make_df <- function(col = "taxon_clean", name = "Chaetoceros decipiens") {
  df <- data.frame(x = name, stringsAsFactors = FALSE)
  names(df) <- col
  df
}

# missing col errors -----------------------------------------------------

test_that("search_worms_priority errors on missing col", {
  expect_error(
    search_worms_priority(make_df(), col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

test_that("search_worms_taxamatch errors on missing col", {
  expect_error(
    search_worms_taxamatch(make_df(), col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

test_that("search_gbif_strict errors on missing col", {
  expect_error(
    search_gbif_strict(make_df(), col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

test_that("resolve_taxonomic_status errors on missing col", {
  expect_error(
    resolve_taxonomic_status(make_df(), col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

test_that("search_worms_fuzzy_minor errors on missing col", {
  expect_error(
    search_worms_fuzzy_minor(make_df(), col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

test_that("get_taxonomy errors on missing col", {
  expect_error(
    get_taxonomy(make_df(), col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

test_that("run_resolution_pipeline errors on missing col", {
  expect_error(
    run_resolution_pipeline(make_df(), col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

# custom col accepted without col-related error --------------------------

test_that("run_resolution_pipeline accepts custom col without col error", {
  df <- make_df(col = "clean_name")
  result <- tryCatch(
    run_resolution_pipeline(df, col = "clean_name"),
    error = function(e) e
  )
  # Any error here must NOT be a missing-column error for "clean_name"
  if (inherits(result, "error")) {
    expect_false(
      grepl(
        "Column 'clean_name' not found",
        conditionMessage(result),
        fixed = TRUE
      )
    )
  } else {
    expect_s3_class(result, "data.frame")
  }
})
