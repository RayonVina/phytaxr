# tests/testthat/test-clean.R
# Tests for Stage 1 cleaning functions (clean.R).
# No API calls are made; all tests operate on local data frames or vectors.

# ---------------------------------------------------------------------------
# ensure_cleaning_schema() -- internal helper
# ---------------------------------------------------------------------------

test_that("ensure_cleaning_schema() errors if df is not a data.frame", {
  expect_error(
    phytaxr:::ensure_cleaning_schema(list(taxon = "Foo")),
    regexp = "data.frame"
  )
})

test_that("ensure_cleaning_schema() errors if 'taxon' column is absent", {
  df <- data.frame(name = "Foo", stringsAsFactors = FALSE)
  expect_error(
    phytaxr:::ensure_cleaning_schema(df),
    regexp = "taxon"
  )
})

test_that("ensure_cleaning_schema() initialises taxon_clean from taxon", {
  df <- data.frame(taxon = "Chaetoceros debilis", stringsAsFactors = FALSE)
  result <- phytaxr:::ensure_cleaning_schema(df)
  expect_equal(result$taxon_clean, "Chaetoceros debilis")
})

test_that("ensure_cleaning_schema() does not overwrite existing taxon_clean", {
  df <- data.frame(
    taxon = "Original name",
    taxon_clean = "Already cleaned",
    stringsAsFactors = FALSE
  )
  result <- phytaxr:::ensure_cleaning_schema(df)
  expect_equal(result$taxon_clean, "Already cleaned")
})

test_that("ensure_cleaning_schema() initialises tax_epithet as NA_character_", {
  df <- data.frame(taxon = "Foo bar", stringsAsFactors = FALSE)
  result <- phytaxr:::ensure_cleaning_schema(df)
  expect_true(is.na(result$tax_epithet))
})

test_that("ensure_cleaning_schema() initialises uncertain as FALSE", {
  df <- data.frame(taxon = "Foo bar", stringsAsFactors = FALSE)
  result <- phytaxr:::ensure_cleaning_schema(df)
  expect_false(result$uncertain)
})

test_that("ensure_cleaning_schema() preserves existing uncertain value", {
  df <- data.frame(taxon = "Foo", uncertain = TRUE, stringsAsFactors = FALSE)
  result <- phytaxr:::ensure_cleaning_schema(df)
  expect_true(result$uncertain)
})

# ---------------------------------------------------------------------------
# normalize_characters()
# ---------------------------------------------------------------------------

test_that("normalize_characters() works on a bare taxon-only data frame", {
  df <- data.frame(
    taxon = "Thalassiosira\u00A0weissflogii",
    stringsAsFactors = FALSE
  )
  result <- normalize_characters(df)
  expect_false(grepl("\u00A0", result$taxon_clean))
})

test_that("normalize_characters() collapses multiple spaces", {
  df <- data.frame(taxon = "Chaetoceros  debilis", stringsAsFactors = FALSE)
  result <- normalize_characters(df)
  expect_equal(result$taxon_clean, "Chaetoceros debilis")
})

# ---------------------------------------------------------------------------
# process_sp_entries()
# ---------------------------------------------------------------------------

test_that("process_sp_entries() moves sp. to tax_epithet", {
  df <- data.frame(taxon = "Thalassiosira sp.", stringsAsFactors = FALSE)
  result <- process_sp_entries(df)
  expect_false(grepl("\\bsp\\b", result$taxon_clean))
  expect_false(is.na(result$tax_epithet))
})

test_that("process_sp_entries() works on bare taxon-only data frame", {
  df <- data.frame(taxon = "Chaetoceros sp.", stringsAsFactors = FALSE)
  result <- process_sp_entries(df)
  expect_s3_class(result, "data.frame")
})

# ---------------------------------------------------------------------------
# process_incertae_entries()
# ---------------------------------------------------------------------------

test_that("process_incertae_entries() sets uncertain=TRUE for cf.", {
  df <- data.frame(
    taxon = "Chaetoceros cf. decipiens",
    stringsAsFactors = FALSE
  )
  result <- process_incertae_entries(df)
  expect_true(result$uncertain)
})

# ---------------------------------------------------------------------------
# run_cleaning_pipeline()
# ---------------------------------------------------------------------------

test_that("run_cleaning_pipeline() accepts a character vector", {
  result <- run_cleaning_pipeline(c(
    "Chaetoceros cf. decipiens",
    "Thalassiosira sp."
  ))
  expect_s3_class(result, "data.frame")
})

test_that("run_cleaning_pipeline() accepts a data.frame with default col", {
  df <- data.frame(
    taxon = c("Nitzschia sp.", "Skeletonema costatum"),
    stringsAsFactors = FALSE
  )
  result <- run_cleaning_pipeline(df)
  expect_s3_class(result, "data.frame")
  expect_true("taxon_clean" %in% names(result))
  expect_true("tax_epithet" %in% names(result))
  expect_true("uncertain" %in% names(result))
})

test_that("run_cleaning_pipeline() accepts a custom col name", {
  df <- data.frame(
    raw_name = c("Emiliania huxleyi s.l."),
    stringsAsFactors = FALSE
  )
  result <- run_cleaning_pipeline(df, col = "raw_name")
  expect_s3_class(result, "data.frame")
  expect_true("taxon_clean" %in% names(result))
})

test_that("run_cleaning_pipeline() errors when col does not exist", {
  df <- data.frame(taxon = "Foo", stringsAsFactors = FALSE)
  expect_error(
    run_cleaning_pipeline(df, col = "nonexistent_col"),
    regexp = "nonexistent_col"
  )
})

test_that("run_cleaning_pipeline() errors on non-df non-character input", {
  expect_error(run_cleaning_pipeline(123))
})

test_that("run_cleaning_pipeline() sets uncertain=TRUE for cf. names in data.frame mode", {
  df <- data.frame(taxon = "Chaetoceros cf. debilis", stringsAsFactors = FALSE)
  result <- run_cleaning_pipeline(df)
  expect_true(result$uncertain[1])
})

test_that("run_cleaning_pipeline() removes sp. from taxon_clean in data.frame mode", {
  df <- data.frame(taxon = "Thalassiosira sp.", stringsAsFactors = FALSE)
  result <- run_cleaning_pipeline(df)
  expect_false(grepl("\\bsp\\b", result$taxon_clean[1]))
})
