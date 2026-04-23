# build_genus_vocabulary() -----------------------------------------------

test_that("build_genus_vocabulary() returns sorted unique genera", {
  df <- data.frame(
    genus = c("Thalassiosira", "Chaetoceros", "Thalassiosira", NA),
    stringsAsFactors = FALSE
  )
  result <- build_genus_vocabulary(df)
  expect_type(result, "character")
  expect_equal(result, c("Chaetoceros", "Thalassiosira"))
})

test_that("build_genus_vocabulary() excludes short names (< 4 chars)", {
  df <- data.frame(
    genus = c("Tha", "Chaetoceros", "Ab"),
    stringsAsFactors = FALSE
  )
  result <- build_genus_vocabulary(df)
  expect_equal(result, "Chaetoceros")
})

test_that("build_genus_vocabulary() returns empty vector when all NA", {
  df <- data.frame(genus = NA_character_, stringsAsFactors = FALSE)
  result <- build_genus_vocabulary(df)
  expect_length(result, 0)
})

# build_epithet_vocabulary() ---------------------------------------------

test_that("build_epithet_vocabulary() extracts epithets correctly", {
  df <- data.frame(
    genus = c("Thalassiosira", "Chaetoceros"),
    species = c("Thalassiosira weissflogii", "Chaetoceros debilis"),
    stringsAsFactors = FALSE
  )
  result <- build_epithet_vocabulary(df)
  expect_type(result, "character")
  expect_true("weissflogii" %in% result)
  expect_true("debilis" %in% result)
})

test_that("build_epithet_vocabulary() excludes short epithets (< 4 chars)", {
  df <- data.frame(
    genus = "Genus",
    species = "Genus sp",
    stringsAsFactors = FALSE
  )
  result <- build_epithet_vocabulary(df)
  expect_length(result, 0)
})

test_that("build_epithet_vocabulary() returns empty vector when species is NA", {
  df <- data.frame(
    genus = "Thalassiosira",
    species = NA_character_,
    stringsAsFactors = FALSE
  )
  result <- build_epithet_vocabulary(df)
  expect_length(result, 0)
})

# normalize_taxonomic_name() ---------------------------------------------

test_that("normalize_taxonomic_name() returns NA for NA input", {
  result <- phytaxr:::normalize_taxonomic_name(NA_character_, character(0))
  expect_true(is.na(result))
})

test_that("normalize_taxonomic_name() splits concatenated genus+epithet", {
  vocab <- c("Chaetoceros")
  result <- phytaxr:::normalize_taxonomic_name("Chaetocerosdebilis", vocab)
  expect_equal(result, "Chaetoceros debilis")
})

test_that("normalize_taxonomic_name() returns original if no rule applies", {
  vocab <- c("Thalassiosira")
  result <- phytaxr:::normalize_taxonomic_name("Skeletonema costatum", vocab)
  expect_equal(result, "Skeletonema costatum")
})

# process_fuzzy_batch() col parameter ------------------------------------

test_that("process_fuzzy_batch errors when col does not exist", {
  df <- data.frame(
    taxon_clean = c("Chaetoceros decipiens"),
    matched_aphiaid = NA_integer_,
    stringsAsFactors = FALSE
  )
  expect_error(
    process_fuzzy_batch(
      df,
      genus_vocab = character(0),
      epithet_vocab = character(0),
      col = "nonexistent_col"
    ),
    regexp = "not found|nonexistent_col|Column"
  )
})

test_that("process_fuzzy_batch accepts a custom col on already-resolved rows", {
  df <- data.frame(
    clean_name = c("Chaetoceros decipiens"),
    taxon_clean = c("Chaetoceros decipiens"),
    matched_aphiaid = 1L, # already resolved: loop should not execute
    stringsAsFactors = FALSE
  )
  result <- process_fuzzy_batch(
    df,
    genus_vocab = character(0),
    epithet_vocab = character(0),
    col = "clean_name",
    checkpoint_file = NULL
  )
  expect_s3_class(result, "data.frame")
})
