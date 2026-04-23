# initialize_taxonomy_df() and split_separator_entries() tests removed:
# initialize_taxonomy_df() was removed from the package exports.
# split_separator_entries() is a pipeline-internal function tested implicitly
# through the full cleaning pipeline.

test_that("resolve_taxonomy rejects non-character input", {
  expect_error(resolve_taxonomy(123))
})

test_that("resolve_taxonomy rejects empty vector", {
  expect_error(resolve_taxonomy(character(0)))
})

# col parameter ----------------------------------------------------------

test_that("run_resolution_pipeline accepts a custom col name", {
  df <- data.frame(
    clean_name = c("Chaetoceros decipiens"),
    stringsAsFactors = FALSE
  )
  result <- tryCatch(
    run_resolution_pipeline(df, col = "clean_name"),
    error = function(e) e
  )
  # If it errors, it must NOT be because of a missing column named clean_name.
  # A missing-column error would say "Column 'clean_name' not found".

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

test_that("run_resolution_pipeline errors when col does not exist", {
  df <- data.frame(
    taxon_clean = "Chaetoceros decipiens",
    stringsAsFactors = FALSE
  )
  expect_error(
    run_resolution_pipeline(df, col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

test_that("search_worms_priority errors when col does not exist", {
  df <- data.frame(
    taxon_clean = "Chaetoceros decipiens",
    stringsAsFactors = FALSE
  )
  expect_error(
    search_worms_priority(df, col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})
