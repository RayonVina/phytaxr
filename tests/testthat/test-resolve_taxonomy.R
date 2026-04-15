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
